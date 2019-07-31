#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_Hypre.H>
#endif

namespace ExtTaoBC 
{
   extern amrex::Vector<amrex::Vector<amrex::Real>> ext_dir_bcs;
   const int lower_boundary = 0;
   const int left_boundary = 1;
   const int upper_boundary = 2;
}

class MyTest
{
public:

    MyTest ();

    void solve ();
    void writePlotfile () const;

public: // make these public for cuda
    void initProbPoisson ();
    void initProbABecLaplacian ();
    void get_number_global_bcs(int& num_lower, int& num_left, int& num_upper);
    void get_number_local_bcs(int& num_lower, int& num_left, int& num_upper);    

private:

    void readParameters ();
    void initData ();
    void solvePoisson ();
    void solveABecLaplacian ();

    int max_level = 1;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;
    
    bool composite_solve = true;

    int prob_type = 1;  // 1. Poisson,  2. ABecLaplacian

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int linop_maxorder = 2;
    bool agglomeration = true;
    bool consolidation = true;
    int max_coarsening_level = 30;
    bool use_hypre = false;
    bool use_petsc = false;
    
#ifdef AMREX_USE_HYPRE
    int hypre_interface_i = 1;  // 1. structed, 2. semi-structed, 3. ij
    amrex::Hypre::Interface hypre_interface = amrex::Hypre::Interface::structed;
#endif

    amrex::BCRec bcs; // BCs for 1 component

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    amrex::Vector<amrex::MultiFab> solution;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> adjoint;
    amrex::Vector<amrex::MultiFab> adjoint_rhs;
    amrex::Vector<amrex::MultiFab> exact_solution;
    amrex::Vector<amrex::MultiFab> acoef;
    amrex::Vector<amrex::MultiFab> bcoef;

    amrex::Real ascalar = 1.e-3;
    amrex::Real bscalar = 1.0;

    amrex::Vector<amrex::Real> dfdp; // {dfdp_x_lo, dfdp_x_hi, dfdp_y_lo, dfdp_y_hi, ...}

    amrex::Real (*target_function)(amrex::Real* physical_coordinates);
};

#endif
