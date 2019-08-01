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

            const auto box_lo = amrex::lbound(bx);
            const auto box_hi = amrex::ubound(bx);

            // fill left edge
            const int k = 0;

            auto key = make_box_corner_tuple(bx);
            auto &vvr = ExtTaoBC::ext_dir_bcs[key];


            bool aligned_lower = false;
            bool aligned_left  = false;
            bool aligned_upper = false;

            int ivec_lower = 0;
            int ivec_left = 0;
            int ivec_upper = 0;

            auto& vbc_lower = vvr(ExtTaoBC::lower_boundary);
            auto& vbc_left = vvr(ExtTaoBC::left_boundary);
            auto& vbc_upper = vvr(ExtTaoBC::upper_boundary);

            // Lower boundary
            if (box_lo.y < domain_lo.y)
            {
                aligned_lower = true;
                const int j = box_lo.y;
                for (int i = box_lo.x; i <= box_hi.x; ++i)
                {
                    dest(i, j, k, dcomp) = vbc_lower(ivec_lower); ivec_lower++;
                }
            }

            // Left boundary
            if (box_lo.x < domain_lo.x)
            {
                aligned_left = true;
                const int i = box_lo.x;
                for (int j = box_lo.y; j <= box_hi.y; ++j)
                {
                    dest(i, j, k, dcomp) = vbc_left(ivec_left); ivec_left++;
                }
            }

            // Upper boundary
            if (box_hi.y > domain_hi.y)
            {
                aligned_upper = true;
                const int j = box_hi.y;
                for (int i = box_lo.x; i <= box_hi.x; ++i)
                {
                    dest(i, j, k, dcomp) = vbc_upper(ivec_upper); ivec_upper++;
                }
            }

            // Left Lower corner
            if (aligned_left && aligned_lower) {
                dest(box_lo.x, box_lo.y, k) = vbc_left(ivec_left); ivec_left++; 
            }

            // Left Upper corner 
            if (aligned_left && aligned_upper) {
                dest(box_lo.x, box_hi.y, k) = vbc_left(ivec_left); ivec_left++; 
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
