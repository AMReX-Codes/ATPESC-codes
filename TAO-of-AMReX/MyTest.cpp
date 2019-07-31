#include "MyTest.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

namespace ExtTaoBC
{
amrex::Vector<amrex::Vector<amrex::Real>> ext_dir_bcs;
}

MyTest::MyTest()
{
    readParameters();
    initData();
}

void MyTest::solve()
{
    solvePoisson(solution, rhs);
}


void MyTest::get_number_bcs(int& num_lower, int& num_left, int& num_upper, int& num_global_bcs)
{
    get_number_local_bcs(num_lower, num_left, num_upper);

    const DomainBox& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    num_global_bcs = 2 * (domain_hi.x - domain_lo.x + 1); // lower, upper edges
    num_global_bcs += domain_hi.y - domain_lo.y + 1; // left edge
    num_global_bcs += 2; // left/lower and left/upper corners
}

void MyTest::get_number_local_bcs(int& num_lower, int& num_left, int& num_upper)
{
    // Get number of boundary values local to this MPI rank
    const DomainBox& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    num_lower = 0;
    num_left = 0;
    num_upper = 0;
    
    for (MFIter mfi(solution); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            bx_lo = lbound(bx);
            bx_hi = ubound(bx);

            bool aligned_lower = false;
            bool aligned_left  = false;
            bool aligned_upper = false;
      
            // check lower boundary
            if (bx_lo.y == domain_lo.y) {
                aligned_lower = true;
                num_lower += (bx_hi.x - bx_lo.x + 1);
            }

            // check left boundary
            if (bx_lo.x == domain_lo.x) {
                aligned_left = true;
                num_left += (bx_hi.y - bx_lo.y + 1);
            }

            // check upper boundary
            if (bx_hi.y == domain_hi.y) {
                aligned_upper = true;
                num_upper += (bx_hi.x - bx_lo.x + 1);
            }

            // add corners to left edge
            if (aligned_left && aligned_lower)
                num_left++;

            if (aligned_left && aligned_upper)
                num_left++;
        }
}


void MyTest::update_boundary_values(int nb, Real *xb,
                                    int nl, Real *xl,
                                    int nt, Real *xt)
{
    ExtTaoBC::ext_dir_bcs[ExtTaoBC::lower_boundary].resize(nb)
    ExtTaoBC::ext_dir_bcs[ExtTaoBC::left_boundary].resize(nl)
    ExtTaoBC::ext_dir_bcs[ExtTaoBC::upper_boundary].resize(nt)

    std::copy(xb, xb + nb, ExtTaoBC::ext_dir_bcs[ExtTaoBC::lower_boundary].begin());
    std::copy(xl, xl + nl, ExtTaoBC::ext_dir_bcs[ExtTaoBC::left_boundary].begin());
    std::copy(xt, xt + nt, ExtTaoBC::ext_dir_bcs[ExtTaoBC::upper_boundary].begin());

    const int lev = 0;
    solution[lev].FillBoundary(geom[lev].periodicity());
    FillDomainBoundary(solution[lev], geom[lev], {bcs});
}

void MyTest::setup_adjoint_system()
{
    // setup the (dR/du)^T * lambda = - \partial f/\partial u linear system
    // by evaluating \partial f/\partial u = \int_V (u(p) - u_t) dV

    adjoint_rhs = 0.0;

    // for right boundary, adjoint_rhs(cell) = -dfdu = target solution(cell) - poisson solution(cell)

    // AMReX-ify this...
    for cell:
      if cell in right boundary:
        adjoint_rhs(cell) = (target solution(cell) - poisson solution(cell)) * AMREX_D_TERM(geom[0].CellSize[0], *geom[0].CellSize[1], *geom[0].CellSize[2]);
}

void MyTest::solve_adjoint_system()
{
    // solve (dR/du)^T * lambda = - \partial f/\partial u

    solvePoisson(adjoint, adjoint_rhs);
}

void MyTest::set_target_solution(Real (*ftarget)(Real* coords))
{
    target_function = ftarget;
}

Real MyTest::calculate_obj_val()
{
    // Return the objective function

    // loop over right boundary w/o corners

    // sum contributions to f from processors

    // communicate f to processors
}

void MyTest::calculate_opt_gradient(Real* dfdp_tao)
{
    // Update for cellwise BCs TODO
    // Fill dfdp_tao TODO


    // calculate df/dp - \partial f/\partial p = (dR/dp)^T * lambda

    // iterate over the domain and reduce sum the contributions along the inner edges
    const auto prob_lo = geom[0].ProbLoArray();
    const auto prob_hi = geom[0].ProbHiArray();

    Real dfdp_x_lo = ReduceSum(adjoint, 0,
                               [=] AMREX_GPU_HOST_DEVICE(Box const &bx, FArrayBox const &fab) -> Real {
                                   Real redval = 0.0;

                                   const Array4<Real> fabarray = fab.array();
                                   const auto lo = lbound(bx);
                                   const auto hi = ubound(bx);

                                   for (int k = lo.z; k <= hi.z; ++k)
                                   {
                                       for (int j = lo.y; j <= hi.y; ++j)
                                       {
                                           for (int i = lo.x; i <= hi.x; ++i)
                                           {
                                               if (i == prob_lo[0])
                                                   redval += fabarray(i, j, k);
                                           }
                                       }
                                   }

                                   return redval;
                               });

    dfdp[0] = dfdp_x_lo;

    Real dfdp_x_hi = ReduceSum(adjoint, 0,
                               [=] AMREX_GPU_HOST_DEVICE(Box const &bx, FArrayBox const &fab) -> Real {
                                   Real redval = 0.0;

                                   const Array4<Real> fabarray = fab.array();
                                   const auto lo = lbound(bx);
                                   const auto hi = ubound(bx);

                                   for (int k = lo.z; k <= hi.z; ++k)
                                   {
                                       for (int j = lo.y; j <= hi.y; ++j)
                                       {
                                           for (int i = lo.x; i <= hi.x; ++i)
                                           {
                                               if (i == prob_hi[0])
                                                   redval += fabarray(i, j, k);
                                           }
                                       }
                                   }

                                   return redval;
                               });

    dfdp[1] = dfdp_x_hi;

#if (AMREX_SPACEDIM >= 2)
    Real dfdp_y_lo = ReduceSum(adjoint, 0,
                               [=] AMREX_GPU_HOST_DEVICE(Box const &bx, FArrayBox const &fab) -> Real {
                                   Real redval = 0.0;

                                   const Array4<Real> fabarray = fab.array();
                                   const auto lo = lbound(bx);
                                   const auto hi = ubound(bx);

                                   for (int k = lo.z; k <= hi.z; ++k)
                                   {
                                       for (int j = lo.y; j <= hi.y; ++j)
                                       {
                                           for (int i = lo.x; i <= hi.x; ++i)
                                           {
                                               if (i == prob_lo[1])
                                                   redval += fabarray(i, j, k);
                                           }
                                       }
                                   }

                                   return redval;
                               });

    dfdp[2] = dfdp_y_lo;

    Real dfdp_y_hi = ReduceSum(adjoint, 0,
                               [=] AMREX_GPU_HOST_DEVICE(Box const &bx, FArrayBox const &fab) -> Real {
                                   Real redval = 0.0;

                                   const Array4<Real> fabarray = fab.array();
                                   const auto lo = lbound(bx);
                                   const auto hi = ubound(bx);

                                   for (int k = lo.z; k <= hi.z; ++k)
                                   {
                                       for (int j = lo.y; j <= hi.y; ++j)
                                       {
                                           for (int i = lo.x; i <= hi.x; ++i)
                                           {
                                               if (i == prob_hi[1])
                                                   redval += fabarray(i, j, k);
                                           }
                                       }
                                   }

                                   return redval;
                               });

    dfdp[3] = dfdp_y_hi;
#endif

#if (AMREX_SPACEDIM == 3)
    Real dfdp_z_lo = ReduceSum(adjoint, 0,
                               [=] AMREX_GPU_HOST_DEVICE(Box const &bx, FArrayBox const &fab) -> Real {
                                   Real redval = 0.0;

                                   const Array4<Real> fabarray = fab.array();
                                   const auto lo = lbound(bx);
                                   const auto hi = ubound(bx);

                                   for (int k = lo.z; k <= hi.z; ++k)
                                   {
                                       for (int j = lo.y; j <= hi.y; ++j)
                                       {
                                           for (int i = lo.x; i <= hi.x; ++i)
                                           {
                                               if (i == prob_lo[2])
                                                   redval += fabarray(i, j, k);
                                           }
                                       }
                                   }

                                   return redval;
                               });

    dfdp[4] = dfdp_z_lo;

    Real dfdp_z_hi = ReduceSum(adjoint, 0,
                               [=] AMREX_GPU_HOST_DEVICE(Box const &bx, FArrayBox const &fab) -> Real {
                                   Real redval = 0.0;

                                   const Array4<Real> fabarray = fab.array();
                                   const auto lo = lbound(bx);
                                   const auto hi = ubound(bx);

                                   for (int k = lo.z; k <= hi.z; ++k)
                                   {
                                       for (int j = lo.y; j <= hi.y; ++j)
                                       {
                                           for (int i = lo.x; i <= hi.x; ++i)
                                           {
                                               if (i == prob_hi[2])
                                                   redval += fabarray(i, j, k);
                                           }
                                       }
                                   }

                                   return redval;
                               });

    dfdp[5] = dfdp_z_hi;
#endif
}

void MyTest::solvePoisson(amrex::Vector<amrex::MultiFab> &solution,
                          amrex::Vector<amrex::MultiFab> &rhs)
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

    const int nlevels = geom.size();

    if (composite_solve)
    {

        MLPoisson mlpoisson(geom, grids, dmap, info);

        mlpoisson.setMaxOrder(linop_maxorder);

        // This is a 3d problem with Dirichlet BC
        mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)},
                              {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)});

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            mlpoisson.setLevelBC(ilev, &solution[ilev]);
        }

        MLMG mlmg(mlpoisson);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
        if (use_hypre)
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
            mlmg.setHypreInterface(hypre_interface);
        }
#endif
#ifdef AMREX_USE_PETSC
        if (use_petsc)
        {
            mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
        }
#endif

        mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    }
    else
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            MLPoisson mlpoisson({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info);

            mlpoisson.setMaxOrder(linop_maxorder);

            // This is a 3d problem with Dirichlet BC
            mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet)},
                                  {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet)});

            if (ilev > 0)
            {
                mlpoisson.setCoarseFineBC(&solution[ilev - 1], ref_ratio);
            }

            mlpoisson.setLevelBC(0, &solution[ilev]);

            MLMG mlmg(mlpoisson);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
            if (use_hypre)
            {
                mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
                mlmg.setHypreInterface(hypre_interface);
            }
#endif
#ifdef AMREX_USE_PETSC
            if (use_petsc)
            {
                mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
            }
#endif

            mlmg.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
        }
    }
}

void MyTest::readParameters()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("ref_ratio", ref_ratio);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("composite_solve", composite_solve);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("max_coarsening_level", max_coarsening_level);

#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
    pp.query("hypre_interface", hypre_interface_i);
    if (hypre_interface_i == 1)
    {
        hypre_interface = Hypre::Interface::structed;
    }
    else if (hypre_interface_i == 2)
    {
        hypre_interface = Hypre::Interface::semi_structed;
    }
    else
    {
        hypre_interface = Hypre::Interface::ij;
    }
#endif
#ifdef AMREX_USE_PETSC
    pp.query("use_petsc", use_petsc);
#endif
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!(use_hypre && use_petsc),
                                     "use_hypre & use_petsc cannot be both true");
}

void MyTest::initData()
{
    int nlevels = 1;
    geom.resize(nlevels);
    grids.resize(nlevels);
    dmap.resize(nlevels);

    solution.resize(nlevels);
    rhs.resize(nlevels);
    exact_solution.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0., 0., 0.)}, {AMREX_D_DECL(1., 1., 1.)});
    Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0, 0, 0)}, IntVect{AMREX_D_DECL(n_cell - 1, n_cell - 1, n_cell - 1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids[ilev].define(domain);
        grids[ilev].maxSize(max_grid_size);
        domain.grow(-n_cell / 4); // fine level cover the middle of the coarse domain
        domain.refine(ref_ratio);
    }

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        solution[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        exact_solution[ilev].define(grids[ilev], dmap[ilev], 1, 0);
    }

    // set up boundary conditions:
    // - first set all boundary conditions to external Dirichlet
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        bcs.setLo(i, BCType::ext_dir);
        bcs.setHi(i, BCType::ext_dir);
    }

    // - then set outflow to right = first order extrapolation
    bcs.setHi(0, BCType::foextrap);

    ExtTaoBC::ext_dir_bcs.resize(3);
    ExtTaoBC::ext_dir_bcs[0].resize(n_cell + 2);
    ExtTaoBC::ext_dir_bcs[1].resize(n_cell + 2);
    ExtTaoBC::ext_dir_bcs[2].resize(n_cell + 2);

    initProbPoisson();
}
