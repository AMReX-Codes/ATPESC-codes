#include <AMReX.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MacProjector.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <AmrCoreAdv.H>

using namespace amrex;

// Define the EB geometry -- in this case a cylinder 
void
AmrCoreAdv::make_eb_geometry()
{
    {
        int mg_verbose = 0;
        int cg_verbose = 0;
        int use_hypre  = 0;

        Real obstacle_radius = 0.10;

        amrex::Vector<amrex::RealArray> obstacle_center = { {AMREX_D_dECL(0.3,0.2,0.5)} ;

		int direction =  2;
        Real height   = -1.0;  // Putting a negative number for height means it extends beyond the domain

        // The "false" below is the boolean that determines if the fluid is inside ("true") or
        //     outside ("false") the object(s)

        Array<EB2::CylinderIF,9> obstacles{
            EB2::CylinderIF(    obstacle_radius, height, direction, obstacle_center[ 0], false)};

        // An example of how to have multiple obstacles
        // auto group_1 = EB2::makeUnion(obstacles[0],obstacles[1],obstacles[2]);

        auto gshop   = EB2::makeShop(obstacles[0]);

        EB2::Build(gshop, geom[0], required_coarsening_level, max_coarsening_level);

        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[0]);

        // options are basic, volume, or full
        EBSupport ebs = EBSupport::full;

        // number of ghost cells for each of the 3 EBSupport types
        Vector<int> ng_ebs = {2,2,2};

        // This object provides access to the EB database in the format of basic AMReX objects
        // such as BaseFab, FArrayBox, FabArray, and MultiFab
        factory = new EBFArrayBoxFactory(eb_level, geom[0], grids[0], dmap[0], ng_ebs, ebs);
}
