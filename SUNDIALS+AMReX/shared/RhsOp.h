#ifndef RHSOP_H
#define RHSOP_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

/* void ComputeAdvection(amrex::MultiFab& sol,
                      amrex::MultiFab& advection,
                      amrex::Geometry& geom,
                      int comp, amrex::Real advCoeff); */

void ComputeAdvectionUpwind(amrex::MultiFab& sol,
                            amrex::MultiFab& advection,
                            amrex::Geometry& geom,
                            int comp,
                            amrex::Real advCoeffx,
                            amrex::Real advCoeffy);

void ComputeDiffusion(amrex::MultiFab& sol,
                      amrex::MultiFab& diff_mf,
                      amrex::MultiFab& fx_mf,
                      amrex::MultiFab& fy_mf,
                      amrex::Geometry& geom,
                      int comp,
                      amrex::Real diffCoeffx,
                      amrex::Real diffCoeffy);

#endif
