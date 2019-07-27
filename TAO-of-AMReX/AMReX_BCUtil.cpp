#include <AMReX_BCUtil.H>
#include <AMReX_PhysBCFunct.H>
#include "MyTest.H"

namespace amrex
{

namespace
{

void fill_extdir(Box const &bx, Array4<Real> const &dest,
                 const int dcomp, const int numcomp,
                 GeometryData const &geom, const Real time,
                 const BCRec *bcr, const int bcomp,
                 const int orig_comp)
{
    const auto &domain_box = geom.Domain();
    const auto &domain_lo = amrex::lbound(domain_box);
    const auto &domain_hi = amrex::ubound(domain_box);

    const auto &box_lo = amrex::lbound(bx);
    const auto &box_hi = amrex::ubound(bx);

    // fill left edge
    const int k = 0;

    if (box_lo.x < domain_lo.x)
    {
        const int i = box_lo.x;
        for (int j = box_lo.y; j <= box_hi.y; ++j)
        {
            dest(i, j, k, dcomp) = ExtTaoBC::ext_dir_bcs[ExtTaoBC::left_boundary][j];
        }
    }

    if (box_lo.y < domain_lo.y)
    {
        const int j = box_lo.y;
        for (int i = box_lo.x; i <= box_hi.x; ++i)
        {
            dest(i, j, k, dcomp) = ExtTaoBC::ext_dir_bcs[ExtTaoBC::lower_boundary][i];
        }
    }

    if (box_hi.y > domain_hi.y)
    {
        const int j = box_hi.y;
        for (int i = box_lo.x; i <= box_hi.x; ++i)
        {
            dest(i, j, k, dcomp) = ExtTaoBC::ext_dir_bcs[ExtTaoBC::upper_boundary][i];
        }
    }
}
} // namespace

void FillDomainBoundary(MultiFab &phi, const Geometry &geom, const Vector<BCRec> &bc)
{
    if (geom.isAllPeriodic())
        return;
    if (phi.nGrow() == 0)
        return;

    AMREX_ALWAYS_ASSERT(phi.ixType().cellCentered());

    CpuBndryFuncFab cpu_bndry_func(fill_extdir);
    PhysBCFunct<CpuBndryFuncFab> physbcf(geom, bc, cpu_bndry_func);
    physbcf.FillBoundary(phi, 0, phi.nComp(), 0.0, 0);
}

} // namespace amrex
