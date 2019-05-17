#include "Reactions.h"

using namespace amrex;

// GrayScott reactions
void ComputeReactionsGS(MultiFab& sol_mf, MultiFab& react_mf, Real A, Real B)
{
   for (MFIter mfi(sol_mf); mfi.isValid(); ++mfi)
   {
     const Box& bx = mfi.validbox();
     Array4<Real> const& sol_fab = sol_mf.array(mfi);
     Array4<Real> const& react_fab = react_mf.array(mfi);
     const auto lo = lbound(bx);
     const auto hi = ubound(bx);

     for (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) {
         Real temp = sol_fab(i,j,0,0) * sol_fab(i,j,0,1) * sol_fab(i,j,0,1);
         react_fab(i,j,0,0) += A * (1.0 - sol_fab(i,j,0,0)) - temp;
         react_fab(i,j,0,1) += temp - (A + B) * sol_fab(i,j,0,1);
       }
     }
   }

}
