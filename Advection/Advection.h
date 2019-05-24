#ifndef ADVECTION_H
#define ADVECTION_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

// user-data structure passed through SUNDIALS to RHS functions
struct AdvectionProblem
{
   amrex::Geometry* geom;
   amrex::Real      advCoeffx;
   amrex::Real      advCoeffy;
};

// Run problem
void DoProblem();

// ODE RHS functions called by SUNDIALS
int ComputeRhsAdv(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                  void* problem);

// Set the ODE initial condition
void FillInitConds2D(amrex::MultiFab& sol,
                     const amrex::Geometry& geom);

// Parse the problem input file
void ParseInputs(int& n_cell, int& max_grid_size, int& stepper,
                 amrex::Real& tfinal, amrex::Real& dtout,
                 int& plot_int, AdvectionProblem& problem);

// Decompose the problem in space
void SetUpGeometry(amrex::BoxArray& ba,
                   amrex::Geometry& geom,
                   AdvectionProblem& problem,
                   int n_cell, int max_grid_size);

// Advance the solution in time with CVODE
void ComputeSolutionCV(N_Vector nv_sol, AdvectionProblem* problem,
                       amrex::Real tfinal, amrex::Real dtout, int plot_int);

// Advance the solution in time with ARKode ARKStep
void ComputeSolutionARK(N_Vector nv_sol, AdvectionProblem* problem,
                        amrex::Real tfinal, amrex::Real dtout, int plot_int);

#endif
