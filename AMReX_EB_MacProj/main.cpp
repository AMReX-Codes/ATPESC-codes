#include <AMReX_Particles.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MacProjector.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_TracerParticles.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>



using namespace amrex;

void write_plotfile(int step_counter, const auto& geom, const auto& plotmf, const auto& pc)
{
    std::stringstream sstream;
    sstream << "plt" << std::setw(5) << std::setfill('0') << step_counter;
    std::string plotfile_name = sstream.str();

    amrex::Print() << "Writing " << plotfile_name << std::endl;    
    
    EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                {"before-vx", "before-vy",
#if (AMREX_SPACEDIM == 3)
                                        "before-vz",
#endif
                                        "divu-before",
                                        "xvel", "yvel",
#if (AMREX_SPACEDIM == 3)
                                        "after-vz",
#endif
                                        "divu-after"},
                                geom, 0.0, 0);
    pc.Checkpoint(plotfile_name, "Tracer", true); //Write Tracers to plotfile 
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        int verbose = 1;
        int n_cell = 128;
        int max_grid_size = 32;
        std::string initial_tracer_file = "";
        Real max_time = 1.0;
        int max_steps = 100;
        Real time_step = 0.01;

        // read parameters
        {
            ParmParse pp;
            pp.query("verbose", verbose);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("initial_tracer_file", initial_tracer_file);
            pp.query("max_time", max_time);
            pp.query("max_steps", max_steps);
            pp.query("time_step", time_step);
        }
        int n_cell_x = 2*n_cell;

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(2.,1.,1.)});
            Array<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(0,1,1)};
            Geometry::Setup(&rb, 0, isp.data());
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell_x-1,n_cell-1,n_cell-1)});
            geom.define(domain);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible

#if 0
        // build a simple geometry using the "eb2." parameters in the inputs file
        EB2::Build(geom, required_coarsening_level, max_coarsening_level);
#else
        int n_sphere = 2;
        if (n_sphere == 1)
            {
                EB2::SphereIF sphere_0(0.25, {AMREX_D_DECL(0.5,0.5,0.5)}, false);
                auto gshop = EB2::makeShop(sphere_0);
                EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);

            } else {

            EB2::SphereIF sphere_11(0.1, {AMREX_D_DECL(0.3,0.2,0.5)}, false);
            EB2::SphereIF sphere_12(0.1, {AMREX_D_DECL(0.3,0.5,0.5)}, false);
            EB2::SphereIF sphere_13(0.1, {AMREX_D_DECL(0.3,0.8,0.5)}, false);
            auto column_1 = EB2::makeUnion(sphere_11, sphere_12,sphere_13);

            EB2::SphereIF sphere_21(0.1, {AMREX_D_DECL(0.7,0.3,0.5)}, false);
            EB2::SphereIF sphere_22(0.1, {AMREX_D_DECL(0.7,0.6,0.5)}, false);
            EB2::SphereIF sphere_23(0.1, {AMREX_D_DECL(0.7,0.9,0.5)}, false);
            auto column_2 = EB2::makeUnion(sphere_21, sphere_22,sphere_23);

            EB2::SphereIF sphere_31(0.1, {AMREX_D_DECL(1.1,0.2,0.5)}, false);
            EB2::SphereIF sphere_32(0.1, {AMREX_D_DECL(1.1,0.5,0.5)}, false);
            EB2::SphereIF sphere_33(0.1, {AMREX_D_DECL(1.1,0.8,0.5)}, false);
            auto column_3 = EB2::makeUnion(sphere_31, sphere_32,sphere_33);

            auto all = EB2::makeUnion(column_1,column_2,column_3);

            auto gshop  = EB2::makeShop(all);
            EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
        }
#endif
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom);

        // options are basic, volume, or full
        EBSupport ebs = EBSupport::full;

        // number of ghost cells for each of the 3 EBSupport types
        Vector<int> ng_ebs = {2,2,2};

        // This object provides access to the EB database in the format of basic AMReX objects
        // such as BaseFab, FArrayBox, FabArray, and MultiFab
        EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);


        // Initialize Particles
        TracerParticleContainer TracerPC(geom, dmap, grids);
        TracerPC.InitFromAsciiFile(initial_tracer_file, 0);


        // store plotfile variables; velocity-before, div-before, velocity-after, div-after
        MultiFab plotfile_mf;
        plotfile_mf.define(grids, dmap, 2*AMREX_SPACEDIM+2, 0, MFInfo(), factory);

        Array<MultiFab,AMREX_SPACEDIM> vel;
        Array<MultiFab,AMREX_SPACEDIM> beta;

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].define (amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 1, MFInfo(), factory);
            beta[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 0, MFInfo(), factory);
            beta[idim].setVal(1.0);
        }

        // set initial velocity to u=(1,0,0)
        AMREX_D_TERM(vel[0].setVal(1.0);,
                     vel[1].setVal(0.0);,
                     vel[2].setVal(0.0););

        // copy velocity into plotfile
        average_face_to_cellcenter(plotfile_mf,0,amrex::GetArrOfConstPtrs(vel));

        // compute and output divergence, then copy into plofile
        MultiFab divu(grids, dmap, 1, 0, MFInfo(), factory);
        EB_computeDivergence(divu, amrex::GetArrOfConstPtrs(vel), geom);
        amrex::Print() << "\nmax-norm of divu before projection is " << divu.norm0() << "\n" << std::endl;
        plotfile_mf.copy(divu,0,AMREX_SPACEDIM,1);

        MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // mac velocity
                             {amrex::GetArrOfConstPtrs(beta)}, // beta
                             {geom});                          // Geometry

        macproj.setVerbose(verbose);
        macproj.setCGVerbose(verbose);

        macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                          LinOpBCType::Periodic,
                                          LinOpBCType::Periodic)},
            {AMREX_D_DECL(LinOpBCType::Dirichlet,
                          LinOpBCType::Periodic,
                          LinOpBCType::Periodic)});

        Real reltol = 1.e-8;

        // macproj.setBottomSolver(MLMG::BottomSolver::hypre);
        // macproj.setBottomSolver(MLMG::BottomSolver::bicgcg);
        macproj.project(reltol);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].FillBoundary(geom.periodicity());
        }

        // copy velocity into plotfile
        average_face_to_cellcenter(plotfile_mf,AMREX_SPACEDIM+1,amrex::GetArrOfConstPtrs(vel));

        // compute and output divergence, then copy into plofile
        EB_computeDivergence(divu, amrex::GetArrOfConstPtrs(vel), geom);
        amrex::Print() << "\nmax-norm of divu after projection is " << divu.norm0() << "\n" << std::endl;
        plotfile_mf.copy(divu,0,2*AMREX_SPACEDIM+1,1);

        Real time = 0.0;
        for (int i = 0; i < max_steps; i++)
        {
            if (time < max_time) {
                time_step = std::min(time_step, max_time - time);

                amrex::Print() << "\nTimestep " << i << ", Time = " << time << std::endl;
                amrex::Print() << "Advecting particles with Umac for timestep " << time_step << std::endl;
                // Step Particles
                TracerPC.AdvectWithUmac(vel.data(), 0, time_step);

                // Write to a plotfile
                write_plotfile(i, geom, plotfile_mf, TracerPC);

                // Increment time
                time += time_step;
            } else {
                // Write to a plotfile
                write_plotfile(i, geom, plotfile_mf, TracerPC);
                break;
            }
        }
    }

#if 0
    // This is what I'm thinking the particle loop will look like ... 
    {
        MyParticleContainer particles(geom,dm,ba);
        particles.InitParticles();

        for (int step = 0; step < parms.nsteps; ++step)
        {
             amrex::Print() << "    Time step: " <<  step << std::endl;
             amrex::Print() << " Number of particles is " << particles.TotalNumberOfParticles() << std::endl;
             particles.MoveParticles();
             particles.Redistribute();
        }
    }
#endif

    amrex::Finalize();
}
