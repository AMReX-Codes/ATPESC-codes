#ifndef GRAYSCOTT_H
#define GRAYSCOTT_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

struct GrayScottProblem
{
   amrex::Geometry* geom;
   amrex::Real diffCoeffU;
   amrex::Real diffCoeffV;
   amrex::Real A;
   amrex::Real B;
   amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* flux;
};

void DoProblem();

void ComputeReactions2D(amrex::MultiFab& sol,
                        amrex::MultiFab& reactions,
                        GrayScottProblem& problem);

int ComputeReactionsNV(realtype t,
                       N_Vector sol,
                       N_Vector reactions,
                       void* problem);

void ComputeLaplacian(amrex::MultiFab& sol,
                      amrex::MultiFab& laplacian,
                      GrayScottProblem& problem);

int ComputeLaplacianNV(realtype t,
                       N_Vector sol,
                       N_Vector laplacian,
                       void* problem);

void FillInitConds2D(amrex::MultiFab& sol,
                     const amrex::Geometry& geom);

void ParseInputs(int& n_cell, int& max_grid_size, int& stepper, int& nsteps,
                 int& plot_int, GrayScottProblem& problem);

void SetUpGeometry(amrex::BoxArray& ba,
                   amrex::Geometry& geom,
                   GrayScottProblem& problem,
                   int n_cell, int max_grid_size);

void ComputeSolutionMRI(N_Vector sol_old,
                        N_Vector sol_new,
                        GrayScottProblem* problem,
                        int nsteps, int plot_int);

void ComputeSolutionARK(N_Vector sol_old,
                        N_Vector sol_new,
                        GrayScottProblem* problem,
                        int nsteps, int plot_int);

#endif
