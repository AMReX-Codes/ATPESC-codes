#ifndef ADVECTIONDIFFUSION_H
#define ADVECTIONDIFFUSION_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

// user-data structure passed through SUNDIALS to RHS functions
struct ProblemData
{
   amrex::Geometry* geom;
   amrex::Real advCoeffx;
   amrex::Real advCoeffy;
   amrex::Real diffCoeff;
   amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* flux;
};

// user-data structure for problem options
struct ProblemOpt
{
   int n_cell;
   int max_grid_size;
   int plot_int;
   int stepper;
   int cvode_method;
   int arkode_order;
   int nls_method;
   int nls_fp_acc;
   int rhs_adv;
   int rhs_diff;
   amrex::Real rtol;
   amrex::Real atol;
   amrex::Real tfinal;
   amrex::Real dtout;
};

// Run problem
void DoProblem();

// Advection RHS function
int ComputeRhsAdv(amrex::Real t, N_Vector nv_sol, N_Vector nv_rhs,
                  void* data);

// Diffusion RHS function
int ComputeRhsDiff(amrex::Real t, N_Vector nv_sol, N_Vector nv_rhs,
                   void* data);

// Advection-Diffusion RHS function
int ComputeRhsAdvDiff(amrex::Real t, N_Vector nv_sol, N_Vector nv_rhs,
                      void* data);

// Set the ODE initial condition
void FillInitConds2D(amrex::MultiFab& sol, const amrex::Geometry& geom);

// Parse the problem input file
void ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data);

// Decompose the problem in space
void SetUpGeometry(amrex::BoxArray& ba, amrex::Geometry& geom,
                   ProblemOpt& prob_opt, ProblemData& prob_data);

// Advance the solution in time with CVODE
void ComputeSolutionCV(N_Vector nv_sol, ProblemOpt* prob_opt,
                       ProblemData* prob_data);

// Advance the solution in time with ARKode ARKStep
void ComputeSolutionARK(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data);

#endif
