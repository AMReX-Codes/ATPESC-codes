#include <AMReX.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MacProjector.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <AmrCoreAdv.H>

using namespace amrex;

// Project the face-centered velocities to be divergence free, even around an obstacle
void
AmrCoreAdv::mac_project_velocity()
{
    {
        int mg_verbose = 0;
        int cg_verbose = 0;
        int use_hypre  = 0;

        Real obstacle_radius = 0.10;

        // read parameters
        {
            ParmParse pp;
            pp.query("mg_verbose", mg_verbose);
            pp.query("cg_verbose", cg_verbose);
            pp.query("use_hypre", use_hypre);
        }

#ifndef AMREX_USE_HYPRE
        if (use_hypre == 1)
           amrex::Abort("Cant use hypre if we dont build with USE_HYPRE=TRUE");
#endif

        Vector<Array<MultiFab,AMREX_SPACEDIM>> beta(max_level+1);

        // allocate face-centered velocities and face-centered beta coefficient
        for (int lev = 0; lev <= max_level; ++lev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
#ifdef AMREX_USE_EB
                beta[lev][idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)), 
                                       dmap[lev], 1, 0, MFInfo(), factory);
#else
                beta[lev][idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim))); 
#endif
                beta[lev][idim].setVal(1.0);  // set beta to 1
            }
        }

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible

        LPInfo lp_info;

        // If we want to use hypre to solve the full problem we need to not coarsen inside AMReX
        if (use_hypre)
            lp_info.setMaxCoarseningLevel(0);

        MacProjector macproj(GetVecOfArrOfPtrs(facevel),       // normal velocity on faces
                             MLMG::Location::FaceCenter,       // Location of vel
                             GetVecOfArrOfConstPtrs(beta),     // coefficients on faces
                             MLMG::Location::FaceCenter,       // Location of beta
                             MLMG::Location::CellCenter,       // Location of solution variable phi
                             {geom},                           // the geometry object
                             lp_info);                         // structure for passing info to the operator

        // Set bottom-solver to use hypre instead of native BiCGStab
        if (use_hypre)
            macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::hypre);

        // Hard-wire the boundary conditions to be Neumann on the low x-face, Dirichlet
        // on the high x-face, and periodic in the other two directions
        // (the first argument is for the low end, the second is for the high end)
        macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                          LinOpBCType::Dirichlet,
                                          LinOpBCType::Dirichlet)},
                            {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                          LinOpBCType::Dirichlet,
                                          LinOpBCType::Dirichlet)});

        macproj.setVerbose(mg_verbose);
        macproj.getMLMG().setBottomVerbose(cg_verbose);

        // Define the relative tolerance
        Real reltol = 1.e-8;

        // Define the absolute tolerance; note that this argument is optional
        Real abstol = 1.e-15;

        // Solve for phi and subtract from the velocity to make it divergence-free
        // Note that the normal velocities are at face centers (not centroids)
        macproj.project(reltol,abstol);
    }
}
