#include "MyTest.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

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

void MyTest::update_counter()
{
    iteration_counter++;
}

std::string MyTest::get_iteration_filename(std::string filename)
{
    std::stringstream sstream;
    sstream << filename << std::setw(5) << std::setfill('0') << iteration_counter;
    std::string file_with_counter = sstream.str();
    return file_with_counter;
}

void MyTest::write_plotfile()
{
    VisMF::Write(solution, get_iteration_filename("solution"));
    VisMF::Write(adjoint, get_iteration_filename("adjoint"));
    VisMF::Write(adjoint_rhs, get_iteration_filename("adjoint_rhs"));
    VisMF::Write(rhs, get_iteration_filename("rhs"));
    if (iteration_counter == 0)
        VisMF::Write(exact_solution, get_iteration_filename("exact_solution"));
}
void MyTest::get_number_global_bcs(int& num_lower, int& num_left, int& num_upper)
{
    const DomainBox& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    num_lower = (domain_hi.x - domain_lo.x + 1); // lower, upper edges

    num_upper = (domain_hi.x - domain_lo.x + 1); // lower, upper edges

    num_left = domain_hi.y - domain_lo.y + 1; // left edge
    num_left += 2; // left/lower and left/upper corners
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

    solution.FillBoundary(geom.periodicity());
    FillDomainBoundary(solution, geom, {bcs});
}

void MyTest::setup_adjoint_system()
{
    // setup the (dR/du)^T * lambda = - \partial f/\partial u linear system
    // by evaluating \partial f/\partial u = \int_V (u(p) - u_t) dV

    adjoint_rhs = 0.0;

    // for right boundary, adjoint_rhs(cell) = -dfdu = target solution(cell) - poisson solution(cell)
    const DomainBox& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    const int k = 0;

    for (MFIter mfi(solution); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        bx_lo = lbound(bx);
        bx_hi = ubound(bx);

        const auto sol_arr = solution[mfi].array();
        const auto exact_sol_arr = exact_solution[mfi].array();
        auto adj_arr = adjoint_rhs[mfi].array();

        // check if we have part of the right boundary
        if (bx_hi.x == domain_hi.x) {
            const int i = bx_hi.x;
            for (int j = bx_lo.y; j <= bx_hi.y; ++j) {
                adj_arr(i, j, k) = exact_sol_arr(i, j, k) - sol_arr(i, j, k);
                adj_arr(i, j, k) *= AMREX_D_TERM(geom[0].CellSize[0], *geom[0].CellSize[1], *geom[0].CellSize[2]);
            }
        }
    }
}

void MyTest::solve_adjoint_system()
{
    // solve (dR/du)^T * lambda = - \partial f/\partial u

    solvePoisson(adjoint, adjoint_rhs);
}

void MyTest::set_target_solution(Real (*ftarget)(Real* coords))
{
    target_function = ftarget;

    update_target_solution();
}

void MyTest::update_target_solution()
{
    const DomainBox& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    const int k = 0;

    for (MFIter mfi(exact_solution); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        bx_lo = lbound(bx);
        bx_hi = ubound(bx);

        auto exact_sol_arr = exact_solution[mfi].array();

        for (int i = bx_lo.x; i <= bx_hi.x; ++i) {
            for (int j = bx_lo.y; j <= bx_hi.y; ++j) {
                IntVect cell_indices;
                AMREX_D_TERM(cell_indices[0] = i;, cell_indices[1] = j;, cell_indices[2] = k;)
                Vector<Real> cell_location(AMREX_SPACEDIM);
                geom[0].CellCenter(cell_indices, cell_location);
                exact_sol_arr(i, j, k) = target_function(cell_location.dataPtr());
            }
        }
    }
}

Real MyTest::calculate_obj_val()
{
    const DomainBox& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    Real fobj_local = 0.0;

    const int k = 0;

    for (MFIter mfi(solution); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        bx_lo = lbound(bx);
        bx_hi = ubound(bx);

        const auto sol_arr = solution[mfi].array();
        const auto exact_sol_arr = exact_solution[mfi].array();

        // check if we have part of the right boundary
        if (bx_hi.x == domain_hi.x) {
            const int i = bx_hi.x;
            // loop over right edge, excluding corners
            for (int j = std::max(bx_lo.y, domain_lo.y+1); j <= std::min(bx_hi.y, domain_hi.y-1); ++j) {
                fobj_local += 0.5 * std::pow(exact_sol_arr(i, j, k) - sol_arr(i, j, k), 2.0);
            }
        }
    }

    fobj_local *= AMREX_D_TERM(geom[0].CellSize[0], *geom[0].CellSize[1], *geom[0].CellSize[2]);

    Real fobj_global = fobj_local;
    ParallelDescriptor::ReduceRealSum(&fobj_global);

    return fobj_global;
}

void MyTest::calculate_opt_gradient(int nlower, Real* glower,
                                    int nleft, Real* gleft,
                                    int nupper, Real* gupper)
{
    // calculate df/dp - \partial f/\partial p = (dR/dp)^T * lambda
    const DomainBox& domain_bx = geom.Domain();
    const auto domain_lo = lbound(domain_bx);
    const auto domain_hi = ubound(domain_bx);

    for (int i = 0; i < nlower; ++i) glower[i] = 0.0;
    for (int i = 0; i < nleft; ++i) gleft[i] = 0.0;
    for (int i = 0; i < nupper; ++i) gupper[i] = 0.0;

    // indices into dfdp arrays from TAO
    int ilower = 0;
    int ileft = 0;
    int iupper = 0;

    const int k = 0;

    for (MFIter mfi(adjoint); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            bx_lo = lbound(bx);
            bx_hi = ubound(bx);

            bool aligned_lower = false;
            bool aligned_left  = false;
            bool aligned_upper = false;

            const auto adjoint_arr = adjoint[mfi].array();

            // For each iteration of the MFIter, we fill TAO arrays in the same
            // order as we set number of entries in get_number_local_bcs().

            // check lower
            if (bx_lo.y == domain_lo.y) {
                aligned_lower = true;
                const int j = bx_lo.y;
                for (int i = bx_lo.x; i <= bx_hi.x; ++i) {
                    glower[ilower] += adjoint_arr(i, j, k);
                    ilower++;
                }
            }

            // check left (not including corners)
            if (bx_lo.x == domain_lo.x) {
                aligned_left = true;
                const int i = bx_lo.x;
                for (int j = bx_lo.y; j <= bx_hi.y; ++j) {
                    gleft[ileft] += adjoint_arr(i, j, k);
                    ileft++;
                }
            }

            // check upper
            if (bx_hi.y == domain_hi.y) {
                aligned_upper = true;
                const int j = bx_hi.y;
                for (int i = bx_lo.x; i <= bx_hi.x; ++i) {
                    gupper[iupper] += adjoint_arr(i, j, k);
                    iupper++;
                }
            }

            // lower left corner
            if (aligned_left && aligned_lower) {
                gleft[ileft] += adjoint_arr(bx_lo.x, bx_lo.y, k);
                ileft++;
            }

            // upper left corner
            if (aligned_left && aligned_upper) {
                gleft[ileft] += adjoint_arr(bx_lo.x, bx_hi.y, k);
                ileft++;
            }
        }
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

    const int nlevels = 1;

    if (composite_solve)
    {

        MLPoisson mlpoisson({geom}, {grids}, {dmap}, info);

        mlpoisson.setMaxOrder(linop_maxorder);

        // This is a 3d problem with Dirichlet BC
        mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)},
                              {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)});

        mlpoisson.setLevelBC(0, &solution);

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

        mlmg.solve(GetVecOfPtrs({solution}), GetVecOfConstPtrs({rhs}), tol_rel, tol_abs);
    }
    else
    {
        MLPoisson mlpoisson({geom}, {grids}, {dmap}, info);

        mlpoisson.setMaxOrder(linop_maxorder);

        // This is a 3d problem with Dirichlet BC
        mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)},
            {AMREX_D_DECL(LinOpBCType::Dirichlet,
                          LinOpBCType::Dirichlet,
                          LinOpBCType::Dirichlet)});

        mlpoisson.setLevelBC(0, &solution);

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

        mlmg.solve({&solution}, {&rhs}, tol_rel, tol_abs);
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
