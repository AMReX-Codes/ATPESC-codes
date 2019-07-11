#ifndef REACTIONS_H
#define REACTIONS_H

#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

void ComputeReactionsGS(amrex::MultiFab& sol_mf,
                        amrex::MultiFab& react_mf,
                        amrex::Real A,
                        amrex::Real B);

#endif
