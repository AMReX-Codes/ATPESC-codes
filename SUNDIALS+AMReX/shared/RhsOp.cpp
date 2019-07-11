#include "RhsOp.h"

using namespace amrex;


/* ---------------------------------------------------------------------------
 * Advection RHS functions
 * ---------------------------------------------------------------------------*/

/*
// Assumes ghost cells already filled
// Adds result to adv_mf MultiFab
void ComputeAdvection(MultiFab& sol_mf, MultiFab& adv_mf,
                      Geometry& geom, int comp, Real advCoeff)
{
   const auto dx = geom.CellSize();
   Real twoDxInv = 0.5 / dx[0]; // assume same in all directions
   Real sideCoeff = advCoeff * twoDxInv;

   int c = comp;  // for brevity
   for (MFIter mfi(sol_mf); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& sol_fab = sol_mf.array(mfi);
      Array4<Real> const& adv_fab = adv_mf.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int j = lo.y; j <= hi.y; ++j) {
         for (int i = lo.x; i <= hi.x; ++i) {
            adv_fab(i,j,0,c) -= sideCoeff * (sol_fab(i+1,j,0,c) - sol_fab(i-1,j,0,c));
            adv_fab(i,j,0,c) -= sideCoeff * (sol_fab(i,j+1,0,c) - sol_fab(i,j-1,0,c));
         }
      }
   }
}
*/

// Assumes ghost cells already filled
// Adds result to adv_mf MultiFab
void ComputeAdvectionUpwind(MultiFab& sol_mf, MultiFab& adv_mf, Geometry& geom,
                            int comp, Real advCoeffx, Real advCoeffy)
{
   const auto dx = geom.CellSize();
   Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
   Real dyInv = 1.0 / dx[1]; // assume same over entire mesh
   Real sideCoeffx = advCoeffx * dxInv;
   Real sideCoeffy = advCoeffy * dyInv;

   int c = comp;  // for brevity
   for (MFIter mfi(sol_mf); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& sol_fab = sol_mf.array(mfi);
      Array4<Real> const& adv_fab = adv_mf.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      // x-direction
      if (advCoeffx > 0) {
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               adv_fab(i,j,0,c) -= sideCoeffx * (sol_fab(i,j,0,c) - sol_fab(i-1,j,0,c));
            }
         }
      } else {
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               adv_fab(i,j,0,c) -= sideCoeffx * (sol_fab(i+1,j,0,c) - sol_fab(i,j,0,c));
            }
         }
      }

      // y-direction
      if (advCoeffy > 0) {
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               adv_fab(i,j,0,c) -= sideCoeffy * (sol_fab(i,j,0,c) - sol_fab(i,j-1,0,c));
            }
         }
      } else {
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               adv_fab(i,j,0,c) -= sideCoeffy * (sol_fab(i,j+1,0,c) - sol_fab(i,j,0,c));
            }
         }
      }
   }
}


/* ---------------------------------------------------------------------------
 * Diffusion RHS functions
 * ---------------------------------------------------------------------------*/

// utility functions for computing diffusion
static void ComputeDiffFlux(MultiFab& sol, MultiFab& fx, MultiFab& fy,
                            Geometry& geom, int comp,
                            Real diffCoeffx, Real diffCoeffy);

static void ComputeDivergence(MultiFab& div, MultiFab& fx, MultiFab& fy,
                              Geometry& geom, int comp);

// Assumes ghots cells are already filled
// Adds result to diff_mf
void ComputeDiffusion(MultiFab& sol, MultiFab& diff_mf, MultiFab& fx_mf,
                      MultiFab& fy_mf, Geometry& geom, int comp,
                      Real diffCoeffx, Real diffCoeffy)
{
   ComputeDiffFlux(sol, fx_mf, fy_mf, geom, comp, diffCoeffx, diffCoeffy);
   ComputeDivergence(diff_mf, fx_mf, fy_mf, geom, comp);
}

/* ---------------------------------------------------------------------------
 * Utility Functions
 * ---------------------------------------------------------------------------*/

// Assumes ghost cells already filled
// Overwrites fx_mf and fy_mf MultiFabs
static void ComputeDiffFlux(MultiFab& sol_mf, MultiFab& fx_mf, MultiFab& fy_mf,
                            Geometry& geom, int comp, Real diffCoeffx, Real diffCoeffy)
{
   const auto dx = geom.CellSize();
   Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
   Real dyInv = 1.0 / dx[1]; // assume same over entire mesh
   Real coeffX = diffCoeffx * dxInv;
   Real coeffY = diffCoeffy * dyInv;

   int c = comp;  // for brevity
   for (MFIter mfi(sol_mf); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& sol = sol_mf.array(mfi);
      Array4<Real> const& fx = fx_mf.array(mfi);
      Array4<Real> const& fy = fy_mf.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      // x-flux
      for (int j = lo.y; j <= hi.y; ++j) {
         for (int i = lo.x; i <= hi.x+1; ++i) {
            // always use zero component for flux
            fx(i,j,0,0) = coeffX * (sol(i,j,0,c) - sol(i-1,j,0,c));
         }
      }

      // y-flux
      for (int j = lo.y; j <= hi.y+1; ++j) {
         for (int i = lo.x; i <= hi.x; ++i) {
            // always use zero component for flux
            fy(i,j,0,0) = coeffY * (sol(i,j,0,c) - sol(i,j-1,0,c));
         }
      }
   }
}

// Assumes ghost cells already filled
// Adds result to div_mf MultiFab
static void ComputeDivergence(MultiFab& div_mf, MultiFab& fx_mf,
                              MultiFab& fy_mf, Geometry& geom, int comp)
{
   const auto dx = geom.CellSize();
   Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
   Real dyInv = 1.0 / dx[1]; // assume same over entire mesh

   int c = comp;  // for brevity
   for (MFIter mfi(div_mf); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& div = div_mf.array(mfi);
      Array4<Real> const& fx = fx_mf.array(mfi);
      Array4<Real> const& fy = fy_mf.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int j = lo.y; j <= hi.y; ++j) {
         for (int i = lo.x; i <= hi.x; ++i) {
            // always use zero component for flux
            div(i,j,0,c) += dxInv * (fx(i+1,j,0,0) - fx(i,j,0,0)) +
                            dyInv * (fy(i,j+1,0,0) - fy(i,j,0,0));
         }
      }
   }
}
