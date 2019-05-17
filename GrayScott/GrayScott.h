#ifndef GRAYSCOTT_H
#define GRAYSCOTT_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

// user-data structure passed through SUNDIALS to RHS functions
struct GrayScottProblem
{
   amrex::Geometry* geom;
   amrex::Real advCoeffU;
   amrex::Real advCoeffV;
   amrex::Real diffCoeffU;
   amrex::Real diffCoeffV;
   amrex::Real A;
   amrex::Real B;
   amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* flux;
};

// Run problem
void DoProblem();

// ODE RHS functions called by SUNDIALS
int ComputeRhsDiff(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                   void* problem);
int ComputeRhsAdv(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                  void* problem);
int ComputeRhsReact(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                    void* problem);

int ComputeRhsAdvDiff(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                      void* problem);
int ComputeRhsAdvReact(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                       void* problem);
int ComputeRhsDiffReact(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                        void* problem);

int ComputeRhsAdvDiffReact(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                           void* problem);

// Set the ODE initial condition
void FillInitConds2D(amrex::MultiFab& sol,
                     const amrex::Geometry& geom);

// Parse the problem input file
void ParseInputs(int& n_cell, int& max_grid_size, int& stepper,
                 amrex::Real& tfinal, amrex::Real& dtout,
                 int& plot_int, GrayScottProblem& problem);

// Decompose the problem in space
void SetUpGeometry(amrex::BoxArray& ba,
                   amrex::Geometry& geom,
                   GrayScottProblem& problem,
                   int n_cell, int max_grid_size);

// Advance the solution in time with CVODE
void ComputeSolutionCV(N_Vector nv_sol, GrayScottProblem* problem,
                       amrex::Real tfinal, amrex::Real dtout, int plot_int);

// Advance the solution in time with ARKode ARKStep
void ComputeSolutionARK(N_Vector nv_sol, GrayScottProblem* problem,
                        amrex::Real tfinal, amrex::Real dtout, int plot_int);

// Advance the solution in time with ARKode MRIStep
void ComputeSolutionMRI(N_Vector nv_sol, GrayScottProblem* problem,
                        amrex::Real tfinal, amrex::Real dtout, int plot_int);

#endif
