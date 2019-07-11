#ifndef HANDSON3_H
#define HANDSON3_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include "Utilities.h"

// user-data structure for problem options
struct ProblemOpt
{
   int plot_int;
   int arkode_order;
   int nls_method;
   int nls_max_iter;
   int nls_fp_acc;
   int ls_max_iter;
   int rhs_adv;
   amrex::Real rtol;
   amrex::Real atol;
   amrex::Real fixed_dt;
   amrex::Real tfinal;
   amrex::Real dtout;
   int max_steps;
   int write_diag;
   int use_preconditioner;
};

// Run problem
void DoProblem();

// Parse the problem input file
void ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data);

// Advance the solution in time with ARKode ARKStep
void ComputeSolutionARK(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data);

#endif
