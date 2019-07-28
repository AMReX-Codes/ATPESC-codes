#include <AMReX_Particles.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_VisMF.H>
#include <AMReX_TagBox.H>
#include <AMReX_ParmParse.H>

#include <MyParticleContainer.H>

using namespace amrex;

void write_plotfile(int step_counter, const auto& geom, const auto& plotmf, const auto& pc)
{
    std::stringstream sstream;
    sstream << "plt" << std::setw(5) << std::setfill('0') << step_counter;
    std::string plotfile_name = sstream.str();
    
    EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                { "proc" },
                                  geom, 0.0, 0);

    pc.Checkpoint(plotfile_name, "Tracer", true); //Write Tracers to plotfile 
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    // Turn off amrex-related output
    amrex::SetVerbose(0);

    {
        int verbose = 0;
        int n_cell = 128;
        int max_grid_size = 32;
        std::string initial_tracer_file = "";
        Real max_time = 1.0;
        int max_steps = 100;
        int plot_int  = 1;
        Real time_step = 0.01;

        Real obstacle_radius = 0.10;
        Real particle_radius = 0.02;

        // read parameters
        {
            ParmParse pp;
            pp.query("verbose", verbose);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("initial_tracer_file", initial_tracer_file);
            pp.query("max_time", max_time);
            pp.query("max_steps", max_steps);
            pp.query("plot_int", plot_int);
            pp.query("time_step", time_step);

            pp.query("obstacle_radius", obstacle_radius);
            pp.query("particle_radius", particle_radius);
        }

        int n_cell_x = n_cell;
        int n_cell_y = n_cell * 8/5;

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.25,2.,1.)});
            Array<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(0,1,1)};
            Geometry::Setup(&rb, 0, isp.data());
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell-1)});
            geom.define(domain);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        MultiFab plotfile_mf;

        amrex::Vector<amrex::RealArray> obstacle_center = {
            {AMREX_D_DECL(0.30,0.3,0.5)},
            {AMREX_D_DECL(0.60,0.3,0.5)},
            {AMREX_D_DECL(0.90,0.3,0.5)},
            {AMREX_D_DECL(0.15,0.7,0.5)},
            {AMREX_D_DECL(0.45,0.7,0.5)},
            {AMREX_D_DECL(0.75,0.7,0.5)},
            {AMREX_D_DECL(1.05,0.7,0.5)},
            {AMREX_D_DECL(0.30,1.1,0.5)}, 
            {AMREX_D_DECL(0.60,1.1,0.5)}, 
            {AMREX_D_DECL(0.90,1.1,0.5)}, 
            {AMREX_D_DECL(0.15,1.5,0.5)},
            {AMREX_D_DECL(0.45,1.5,0.5)},
            {AMREX_D_DECL(0.75,1.5,0.5)},
            {AMREX_D_DECL(1.05,1.5,0.5)}};

        // The "false" below is the boolean that determines if the fluid is inside ("true") or 
        //     outside ("false") the object(s)
        Array<EB2::SphereIF,14> sphere{
            EB2::SphereIF(0.1, obstacle_center[0], false),
            EB2::SphereIF(0.1, obstacle_center[1], false),
            EB2::SphereIF(0.1, obstacle_center[2], false),
            EB2::SphereIF(0.1, obstacle_center[3], false),
            EB2::SphereIF(0.1, obstacle_center[4], false),
            EB2::SphereIF(0.1, obstacle_center[5], false),
            EB2::SphereIF(0.1, obstacle_center[6], false),
            EB2::SphereIF(0.1, obstacle_center[7], false),
            EB2::SphereIF(0.1, obstacle_center[8], false),
            EB2::SphereIF(0.1, obstacle_center[9], false),
            EB2::SphereIF(0.1, obstacle_center[10], false),
            EB2::SphereIF(0.1, obstacle_center[11], false),
            EB2::SphereIF(0.1, obstacle_center[12], false),
            EB2::SphereIF(0.1, obstacle_center[13], false)};

        auto group_1 = EB2::makeUnion(sphere[0],sphere[1],sphere[2]);
        auto group_2 = EB2::makeUnion(sphere[3],sphere[4],sphere[5]);
        auto group_3 = EB2::makeUnion(sphere[6],sphere[7],sphere[8]);
        auto group_4 = EB2::makeUnion(sphere[9],sphere[10],sphere[11]);
        auto group_5 = EB2::makeUnion(sphere[12],sphere[13]);
        auto all     = EB2::makeUnion(group_1,group_2,group_3,group_4,group_5);
        auto gshop9  = EB2::makeShop(all);

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level      = 0;    // typically a huge number so MG coarsens as much as possible
        EB2::Build(gshop9, geom, required_coarsening_level, max_coarsening_level);
   
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
        MyParticleContainer MyPC(geom, dmap, grids);
        MyPC.InitPachinko(initial_tracer_file);

        // Store processor id in the plotfile
        plotfile_mf.define(grids, dmap, 1, 0, MFInfo(), factory);

        amrex::Print() << " \n********************************************************************" << std::endl; 
        amrex::Print() << " Let's advect the particles ... " << std::endl;
        amrex::Print() << "******************************************************************** \n" << std::endl; 

        // copy processor id into plotfile_mf
        int lev = 0;
        for (MFIter mfi = MyPC.MakeMFIter(lev); mfi.isValid(); ++mfi)
            plotfile_mf[mfi].setVal(ParallelDescriptor::MyProc());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

        Real time = 0.0;
        for (int i = 0; i < max_steps; i++)
        {
            if (time < max_time) {
                time_step = std::min(time_step, max_time - time);

                // Step Particles
                MyPC.AdvectPachinko(time_step,obstacle_center,obstacle_radius,particle_radius);

                MyPC.Redistribute();

                // Write to a plotfile
                if (i%plot_int == 0)
                   write_plotfile(i, geom, plotfile_mf, MyPC);

                if (i%100 == 0)
                   amrex::Print() << "Step " << i << ", Time = " << time << std::endl;

                // Increment time
                time += time_step;

            } else 
                break;
        }
    }

    amrex::Finalize();
}
