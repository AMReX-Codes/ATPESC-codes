#ifndef GRAYSCOTT_H
#define GRAYSCOTT_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

// user-data structure passed through SUNDIALS to RHS functions
struct GrayScottProblem
{
   amrex::Geometry* geom;
   amrex::Real diffCoeffU;
   amrex::Real diffCoeffV;
   amrex::Real A;
   amrex::Real B;
   amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* flux;
};

// Run problem
void DoProblem();

// Compute the diffusion term in the ODE RHS
void ComputeDiffusion(amrex::MultiFab& sol, amrex::MultiFab& diffusion,
                      GrayScottProblem& problem);

// ODE RHS wrapper function called by SUNDIALS
int ComputeDiffusionNV(realtype t, N_Vector nv_sol, N_Vector nv_diffusion,
                       void* problem);

// Compute the reaction term in the ODE RHS
void ComputeReactions2D(amrex::MultiFab& sol,
                        amrex::MultiFab& reactions,
                        GrayScottProblem& problem);

// ODE RHS wrapper function called by SUNDIALS
int ComputeReactionsNV(realtype t, N_Vector nv_sol, N_Vector nv_reactions,
                       void* problem);

// Compute both the diffusion and recation term in the ODE RHS
void ComputeDiffusionReactions2D(amrex::MultiFab& sol, amrex::MultiFab& rhs,
                                 GrayScottProblem& problem);

// ODE RHS wrapper function called by SUNDIALS
int ComputeDiffusionReactionsNV(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                                void* problem);

// Set the problem initial condition
void FillInitConds2D(amrex::MultiFab& sol,
                     const amrex::Geometry& geom);

// Parse the problem input file
void ParseInputs(int& n_cell, int& max_grid_size, int& stepper,
                 amrex::Real& tfinal, amrex::Real& dtout,
                 int& plot_int, GrayScottProblem& problem);

// Set the problem decomposition
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
