#include <AMReX_Particles.H>
#include <AMReX.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MacProjector.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_TagBox.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_WriteEBSurface.H>

#include <MyParticleContainer.H>

using namespace amrex;

extern void make_eb_cylinder(const Geometry& geom);
extern void define_velocity(const Real time, const Geometry& geo, Array<MultiFab,AMREX_SPACEDIM>& vel_out, const MultiFab& ccmf);

Real est_time_step(const Geometry& geom, Array<MultiFab,AMREX_SPACEDIM>& vel)
{
    Real cfl = 0.7;

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx      =  geom.CellSize();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        Real est = vel[idim].norm0(0,0,false);
        amrex::Print() << "Max vel in " << idim << "direction is " << est << std::endl;
        dt_est = amrex::min(dt_est, dx[idim]/est);
    }

    ParallelDescriptor::ReduceRealMin(dt_est);

    return dt_est*cfl;
}

void write_plotfile(int step, Real time, const Geometry& geom, MultiFab& plotmf, 
                    MyParticleContainer& pc, int write_ascii)
{
    // Copy processor id into last component of plotfile_mf
    int proc_comp = AMREX_SPACEDIM; 
    for (MFIter mfi(plotmf); mfi.isValid(); ++mfi)
       plotmf[mfi].setVal(ParallelDescriptor::MyProc(),mfi.validbox(),proc_comp,1);

    std::stringstream sstream;
    sstream << "plt" << std::setw(5) << std::setfill('0') << step;
    std::string plotfile_name = sstream.str();
    
#if (AMREX_SPACEDIM == 2)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "xvel", "yvel", "proc" },
                                     geom, time, 0);
#elif (AMREX_SPACEDIM == 3)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "xvel", "yvel", "zvel", "proc" },
                                     geom, time, 0);
#endif

    pc.Checkpoint(plotfile_name, "particles", true); //Write Tracer particles to plotfile

    std::stringstream pstream;
    pstream << "part" << std::setw(5) << std::setfill('0') << step;
    const std::string ascii_filename = pstream.str();

    if (write_ascii)
       pc.WriteAsciiFile(ascii_filename);

}

int main (int argc, char* argv[])
{
    // Turn off amrex-related output
    amrex::SetVerbose(0);

    amrex::Initialize(argc, argv);

    Real strt_time = amrex::second();
    Real eb_strt_time;
    Real eb_stop_time;
    
    {
        int n_cell = 128;
        int max_grid_size = 32;
        std::string particle_file = "";
        Real max_time = 1000.0;
        int max_steps = 100;
        int plot_int  = 1;
        int write_ascii  = 0;
        int use_hypre  = 0;
        Real dt = 0.01;

        Real particle_radius = 0.02;

        amrex::Vector<int> ob_id;

        // read parameters
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("particle_file", particle_file);
            pp.query("max_time", max_time);
            pp.query("max_steps", max_steps);
            pp.query("plot_int", plot_int);
            pp.query("use_hypre", use_hypre);
            pp.query("write_ascii", write_ascii);
            pp.query("time_step", dt);
        }

#ifndef AMREX_USE_HYPRE
        if (use_hypre == 1) 
           amrex::Abort("Cant use hypre if we dont build with USE_HYPRE=TRUE");
#endif

        if (n_cell%8 != 0)
           amrex::Abort("n_cell must be a multiple of 8");

        int n_cell_x = n_cell;
        int n_cell_y = n_cell;
        int n_cell_z = n_cell/8;

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.0,1.0,0.125)});

            Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
            Geometry::Setup(&rb, 0, is_periodic.data());
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)});
            geom.define(domain);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        Array<MultiFab,AMREX_SPACEDIM> vel;
        Array<MultiFab,AMREX_SPACEDIM> beta;
        MultiFab plotfile_mf;

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible

        eb_strt_time = amrex::second();
        make_eb_cylinder(geom);
        eb_stop_time = amrex::second() - eb_strt_time;

        std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > factory =
           makeEBFabFactory(geom, grids, dmap, {4, 4, 2}, EBSupport::full);
        const EBFArrayBoxFactory* ebfact = &(static_cast<amrex::EBFArrayBoxFactory const&>(*factory));

	// Velocities and Beta are face-centered
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].define (amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 1, MFInfo(), *factory);
            beta[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 0, MFInfo(), *factory);
            beta[idim].setVal(1.0);
        }

        // store plotfile variables; velocity, volfrac and processor id
        plotfile_mf.define(grids, dmap, AMREX_SPACEDIM+1, 0, MFInfo(), *factory);

        // Initialize Particles
        MyParticleContainer MyPC(geom, dmap, grids);
        MyPC.InitOnePerCell(0.5, 0.5, 0.5, {});
        // MyPC.InitParticles(particle_file);

        // set initial velocity to u=(1,0,0)
        AMREX_D_TERM(vel[0].setVal(1.0);,
                     vel[1].setVal(0.0);,
                     vel[2].setVal(0.0););

        LPInfo lp_info;

        // If we want to use hypre to solve the full problem we need to not coarsen inside AMReX
        if (use_hypre) 
            lp_info.setMaxCoarseningLevel(0);

        MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // mac velocity
			     MLMG::Location::FaceCenter,       // velocity located on face centers
                             {amrex::GetArrOfConstPtrs(beta)}, // beta
			     MLMG::Location::FaceCenter,       // beta located on face centers
			     MLMG::Location::CellCenter,       // location of velocity divergence
                             {geom},
                             lp_info);                          // structure for passing info to the operator

        // Set bottom-solver to use hypre instead of native BiCGStab 
        if (use_hypre) 
           macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::hypre);

        macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                          LinOpBCType::Neumann,
                                          LinOpBCType::Neumann)},
            {AMREX_D_DECL(LinOpBCType::Neumann,
                          LinOpBCType::Neumann,
                          LinOpBCType::Neumann)});

        Real reltol = 1.e-8;
        Real abstol = 1.e-12;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

        // get max velocity on grid vx_max, vy_max
        Real vx_max = vel[0].max(0, 0);
        Real vy_max = vel[1].max(0, 0);
        // calculate dt_x = dx/vx_max and dt_y
        Real dt_x = geom.CellSize(0)/vx_max;
        Real dt_y = geom.CellSize(1)/vy_max;
        // dt_limit = min(dt_x, dt_y)
        Real dt_limit = std::min(dt_x, dt_y);
        dt = 0.1 * dt_limit;

        if (AMREX_SPACEDIM > 2)
           vel[2].setVal(0.0);

        amrex::Print() << "Writing EB surface" << std::endl;
        WriteEBSurface (grids, dmap, geom, ebfact);

        Real time = 0.0;

        // Write out the initial data
        {
           amrex::Print() << "Creating the initial velocity field " << std::endl;
           define_velocity(time,geom,vel,plotfile_mf);
           macproj.project(reltol, abstol);
           EB_average_face_to_cellcenter(plotfile_mf,0,amrex::GetArrOfConstPtrs(vel));

           amrex::Print() << "Writing the initial data into plt00000\n" << std::endl;
           write_plotfile(0, time, geom, plotfile_mf, MyPC, write_ascii);
        }

        // This computes the first dt
        dt = est_time_step(geom,vel);

        int nstep = 0;

        for (int i = 0; i < max_steps; i++)
        {
            if (time < max_time)
            {
                amrex::Print() << "\nCoarse STEP " << i+1 << " starts ..." << std::endl;

                dt = amrex::min(dt, max_time - time);

                Real t_nph = time + 0.5 * dt;

                define_velocity(t_nph,geom,vel,plotfile_mf);
                macproj.project(reltol, abstol);

                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    vel[idim].FillBoundary(geom.periodicity());
                }

                // Step Particles
                MyPC.AdvectWithUmac(vel.data(), 0, dt);

                MyPC.Redistribute();

                // Increment time
                time += dt;
                nstep++;

                // Write to a plotfile
                if (i%plot_int == 0)
                {
                   average_face_to_cellcenter(plotfile_mf,0,amrex::GetArrOfConstPtrs(vel));
                   write_plotfile(i, time, geom, plotfile_mf, MyPC, write_ascii);
                }

                amrex::Print() << "Coarse STEP " << i+1 << " ends." << " TIME = " << time
//                             << " DT = " << dt << " Sum(Phi) = " << sum_phi << std::endl;
                               << " DT = " << dt << std::endl;

                // Compute lagged dt for next time step based on this half-time velocity
                dt = est_time_step(geom,vel);

            } else {

                // Copy velocity into plotfile
                average_face_to_cellcenter(plotfile_mf,0,amrex::GetArrOfConstPtrs(vel));

                // Write to a plotfile
                write_plotfile(i, time, geom, plotfile_mf, MyPC, write_ascii);
                break;
            }
        }

        // Write plotfile at final tim
        EB_average_face_to_cellcenter(plotfile_mf,0,amrex::GetArrOfConstPtrs(vel));
        write_plotfile(nstep, time, geom, plotfile_mf, MyPC, write_ascii);
    }

    Real stop_time = amrex::second() - strt_time;
    amrex::Print() << "Time to create EB geometry " << eb_stop_time << std::endl;
    amrex::Print() << "Total run time             " << stop_time << std::endl;

    amrex::Finalize();
}
