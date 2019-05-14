#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

#include "DiffOp.h"

using namespace amrex;


// Assumes ghost cells already filled
void ComputeAdvection(MultiFab& sol, MultiFab& advection,
                      Geometry& geom, int comp, Real advCoeff)
{
   const auto dx = geom.CellSize();
   Real twoDxInv = 0.5 / dx[0]; // assume same in all directions
   Real sideCoeff = advCoeff * twoDxInv;

   int c = comp;  // for brevity
   for (MFIter mfi(sol); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& sol_fab = sol.array(mfi);
      Array4<Real> const& adv_fab = advection.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int j = lo.y; j <= hi.y; ++j) {
         for (int i = lo.x; i <= hi.x; ++i) {
            adv_fab(i,j,0,c) += sideCoeff * (sol_fab(i+1,j,0,c) - sol_fab(i-1,j,0,c));
            adv_fab(i,j,0,c) += sideCoeff * (sol_fab(i,j+1,0,c) - sol_fab(i,j-1,0,c));
         }
      }
   }
}

// Assumes ghost cells already filled
void ComputeDiffFlux(MultiFab& sol_mf, MultiFab& fx_mf, MultiFab& fy_mf,
                     Geometry& geom, int comp, Real diffCoeff)
{
   const auto dx = geom.CellSize();
   Real dxInv = 1.0 / dx[0]; // assume same in all directions
   Real dyInv = 1.0 / dx[1]; // assume same in all directions
   Real coeffX = diffCoeff * dxInv;
   Real coeffY = diffCoeff * dyInv;

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
            //amrex::Print() << "fx(" << i << "," << j << ") = " << fx(i,j,0,0) << std::endl;
         }
      }

      // y-flux
      for (int j = lo.y; j <= hi.y+1; ++j) {
         for (int i = lo.x; i <= hi.x; ++i) {
            // always use zero component for flux
            fy(i,j,0,0) = coeffY * (sol(i,j,0,c) - sol(i,j-1,0,c));
            //amrex::Print() << "fy(" << i << "," << j << ") = " << fy(i,j,0,0) << std::endl;
         }
      }
   }
}

void ComputeDivergence(MultiFab& div_mf, MultiFab& fx_mf, MultiFab& fy_mf,
                       Geometry& geom, int comp)
{
   const auto dx = geom.CellSize();
   Real dxInv = 1.0 / dx[0]; // assume same in all directions
   Real dyInv = 1.0 / dx[1]; // assume same in all directions

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
            div(i,j,0,c) = dxInv * (fx(i+1,j,0,0) - fx(i,j,0,0)) +
                           dyInv * (fy(i,j+1,0,0) - fy(i,j,0,0));
            //amrex::Print() << "div(" << i << "," << j << "," << c << ") = " << div(i,j,0,c) << std::endl;
         }
      }
   }
}
