/* -----------------------------------------------------------------------------
 * Time Integration and Nonlinear Solvers Hands-on Lessons with SUNDIALS + AMReX
 * 2019-2022 Argonne Training Program in Extreme-Scale Computing
 *
 * Authors (alphabetical):
 * David Gardner (gardner48@llnl.gov)
 * John Loffeld (loffeld1@llnl.gov)
 * Daniel Reynolds (reynolds@smu.edu)
 * Donald Willcox (dewillcox@lbl.gov)
 * -----------------------------------------------------------------------------
 * Implementation file for SUNDIALS + AMReX 2D advection-diffusion example.
 * ---------------------------------------------------------------------------*/

#include "HandsOn.hpp"

using namespace amrex;

// -----------------------------------------------------------------------------
// Advection-Diffusion main
// -----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  // What time is it now?  We'll use this to compute total run time.
  Real strt_time = amrex::second();

  // Set problem data and options
  ProblemData prob_data;
  ProblemOpt  prob_opt;
  if (ParseInputs(prob_opt, prob_data))
  {
    amrex::Finalize();
    return 1;
  }
  PrintSetup(prob_opt, prob_data);

  // Make BoxArray and Geometry
  BoxArray ba;
  Geometry geom;
  SetUpGeometry(ba, geom, prob_data);

  // How Boxes are distrubuted among MPI processes
  DistributionMapping dm(ba);
  prob_data.dmap = &dm;

  // Allocate the solution MultiFab
  int nGhost = 1;  // number of ghost cells for each array
  int nComp  = 1;  // number of components for each array
  MultiFab sol(ba, dm, nComp, nGhost);

  // Allocate the linear solver coefficient MultiFabs
  MultiFab acoef(ba, dm, nComp, nGhost);
  MultiFab bcoef(ba, dm, nComp, nGhost);
  acoef = 1.0;
  bcoef = 1.0;
  prob_data.acoef = &acoef;
  prob_data.bcoef = &bcoef;

  // Build the flux MultiFabs
  Array<MultiFab, AMREX_SPACEDIM> flux;
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
  {
    // flux(dir) has one component, zero ghost cells, and is nodal in
    // direction dir
    BoxArray edge_ba = ba;
    edge_ba.surroundingNodes(dir);
    flux[dir].define(edge_ba, dm, 1, 0);
  }
  prob_data.flux = &flux;

  // Create an N_Vector wrapper for the solution MultiFab
  sunindextype length = nComp * prob_data.n_cell * prob_data.n_cell;
  N_Vector nv_sol     = amrex::sundials::N_VMake_MultiFab(length, &sol);

  // Set the initial condition
  FillInitConds2D(sol, geom);

  // Integrate in time
  ComputeSolution(nv_sol, &prob_opt, &prob_data);

  // Call the timer again and compute the maximum difference between the start
  // time and stop time over all processors
  Real stop_time = amrex::second() - strt_time;
  const int IOProc = ParallelDescriptor::IOProcessorNumber();
  ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

  // Tell the I/O Processor to write out the "run time"
  amrex::Print() << "Run time = " << stop_time << std::endl;

  amrex::Finalize();

  return 0;
}


// -----------------------------------------------------------------------------
// Parse inputs
// -----------------------------------------------------------------------------


int ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data)
{
  // ParmParse is way of reading inputs from the inputs file
  ParmParse pp;

  pp.query("help", prob_opt.help);
  if (prob_opt.help)
  {
    PrintHelp();
    return 1;
  }

  // Problem options
  pp.query("plot_int",           prob_opt.plot_int);
  pp.query("arkode_order",       prob_opt.arkode_order);
  pp.query("rtol",               prob_opt.rtol);
  pp.query("atol",               prob_opt.atol);
  pp.query("fixed_dt",           prob_opt.fixed_dt);
  pp.query("tfinal",             prob_opt.tfinal);
  pp.query("dtout",              prob_opt.dtout);
  pp.query("max_steps",          prob_opt.max_steps);
  pp.query("write_diag",         prob_opt.write_diag);
  pp.query("nls_max_iter",       prob_opt.nls_max_iter);
  pp.query("ls_max_iter",        prob_opt.ls_max_iter);
  pp.query("rhs_adv",            prob_opt.rhs_adv);
  pp.query("use_preconditioner", prob_opt.use_preconditioner);

  // Grid options
  pp.query("n_cell",        prob_data.n_cell);
  pp.query("max_grid_size", prob_data.max_grid_size);

  // Advection and diffusion coefficient values
  pp.query("advCoeffx",     prob_data.advCoeffx);
  pp.query("advCoeffy",     prob_data.advCoeffy);
  pp.query("diffCoeffx",    prob_data.diffCoeffx);
  pp.query("diffCoeffy",    prob_data.diffCoeffy);

  // ParmParse options prefixed with mlmg.
  ParmParse ppmg("mlmg");

  // MLMG Preconditioner options
  ppmg.query("agglomeration",        prob_data.mg_agglomeration);
  ppmg.query("consolidation",        prob_data.mg_consolidation);
  ppmg.query("max_coarsening_level", prob_data.mg_max_coarsening_level);
  ppmg.query("linop_maxorder",       prob_data.mg_linop_maxorder);
  ppmg.query("max_iter",             prob_data.mg_max_iter);
  ppmg.query("max_fmg_iter",         prob_data.mg_max_fmg_iter);
  ppmg.query("verbose",              prob_data.mg_verbose);
  ppmg.query("bottom_verbose",       prob_data.mg_bottom_verbose);
  ppmg.query("use_hypre",            prob_data.mg_use_hypre);
  ppmg.query("hypre_interface",      prob_data.mg_hypre_interface);
  ppmg.query("use_petsc",            prob_data.mg_use_petsc);
  ppmg.query("tol_rel",              prob_data.mg_tol_rel);

  return 0;
}


// -----------------------------------------------------------------------------
// Set initial state
// -----------------------------------------------------------------------------


void FillInitConds2D(MultiFab& sol, const Geometry& geom)
{
  const auto dx = geom.CellSize();
  const auto prob_lo = geom.ProbLo();
  const auto prob_hi = geom.ProbHi();

  Real sigma = 0.1;
  Real a = 1.0/(sigma*sqrt(2*M_PI));
  Real b = -0.5/(sigma*sigma);

  for (MFIter mfi(sol,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    Array4<Real> const& fab = sol[mfi].array();

    amrex::ParallelFor
      (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
       {
         Real y = prob_lo[1] + (((Real) j) + 0.5) * dx[1];
         Real x = prob_lo[0] + (((Real) i) + 0.5) * dx[0];
         Real r = x * x + y * y;
         fab(i,j,k,n) = a * exp(b * r);
       });
  }
}


// -----------------------------------------------------------------------------
// Setup domain
// -----------------------------------------------------------------------------


void SetUpGeometry(BoxArray& ba, Geometry& geom, ProblemData& prob_data)
{
  // Extract problem options
  int n_cell = prob_data.n_cell;
  int max_grid_size = prob_data.max_grid_size;

  IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
  IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
  Box domain(dom_lo, dom_hi); // cell-centered

  // Initialize the boxarray "ba" from the single box "domain"
  ba.define(domain);

  // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a
  // direction
  ba.maxSize(max_grid_size);

  // This defines the physical box, [-1,1] in each direction.
  RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                   {AMREX_D_DECL(1.0, 1.0, 1.0)});

  // This defines a Geometry object
  Vector<int> is_periodic(AMREX_SPACEDIM, 1);  // periodic in all direction
  geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());

  prob_data.geom = &geom;
  prob_data.grid = &ba;
}


// -----------------------------------------------------------------------------
// User-supplied ODE RHS functions for SUNDIALS
// -----------------------------------------------------------------------------


int ComputeRhsAdv(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
  // extract MultiFabs
  MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
  MultiFab* rhs = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_rhs);

  // extract problem data
  ProblemData* prob_data = (ProblemData*) data;
  Geometry* geom = prob_data->geom;
  Real advCoeffx = prob_data->advCoeffx;
  Real advCoeffy = prob_data->advCoeffy;

  // clear the RHS
  *rhs = 0.0;

  // fill ghost cells
  sol->FillBoundary(geom->periodicity());

  // compute advection
  ComputeAdvectionUpwind(*sol, *rhs, *geom, advCoeffx, advCoeffy);

  return 0;
}


int ComputeRhsDiff(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
  // extract MultiFabs
  MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
  MultiFab* rhs = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_rhs);

  // extract problem data
  ProblemData *prob_data = (ProblemData*) data;
  Geometry* geom = prob_data->geom;
  Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
  Real diffCoeffx = prob_data->diffCoeffx;
  Real diffCoeffy = prob_data->diffCoeffy;

  // fill ghost cells
  sol->FillBoundary(geom->periodicity());

  // clear the RHS
  *rhs = 0.0;

  // compute diffusion
  ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom,
                   diffCoeffx, diffCoeffy);

  return 0;
}


int ComputeRhsAdvDiff(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
  // extract MultiFabs
  MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
  MultiFab* rhs = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_rhs);

  // extract problem data
  ProblemData* prob_data = (ProblemData*) data;
  Geometry* geom = prob_data->geom;
  Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
  Real advCoeffx = prob_data->advCoeffx;
  Real advCoeffy = prob_data->advCoeffy;
  Real diffCoeffx = prob_data->diffCoeffx;
  Real diffCoeffy = prob_data->diffCoeffy;

  // clear the RHS
  *rhs = 0.0;

  // fill ghost cells
  sol->FillBoundary(geom->periodicity());

  // compute advection
  ComputeAdvectionUpwind(*sol, *rhs, *geom, advCoeffx, advCoeffy);

  // compute diffusion
  ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom,
                   diffCoeffx, diffCoeffy);

  return 0;
}


// -----------------------------------------------------------------------------
// Utility functions to compute ODE RHS functions
// -----------------------------------------------------------------------------


// Assumes ghost cells already filled, adds result to adv_mf MultiFab
void ComputeAdvectionUpwind(MultiFab& sol_mf, MultiFab& adv_mf, Geometry& geom,
                            Real advCoeffx, Real advCoeffy)
{
  const auto dx = geom.CellSize();
  Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
  Real dyInv = 1.0 / dx[1]; // assume same over entire mesh
  Real sideCoeffx = advCoeffx * dxInv;
  Real sideCoeffy = advCoeffy * dyInv;

  for (MFIter mfi(sol_mf,TilingIfNotGPU); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    Array4<Real> const& sol_fab = sol_mf[mfi].array();
    Array4<Real> const& adv_fab = adv_mf[mfi].array();

    // x-direction
    if (advCoeffx > 0)
    {
      amrex::ParallelFor
        (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
         {
           adv_fab(i,j,k,n) -= sideCoeffx *
             (sol_fab(i,j,k,n) - sol_fab(i-1,j,k,n));
         });
    }
    else
    {
      amrex::ParallelFor
        (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
         {
           adv_fab(i,j,k,n) -= sideCoeffx *
             (sol_fab(i+1,j,k,n) - sol_fab(i,j,k,n));
         });
    }

    // y-direction
    if (advCoeffy > 0)
    {
      amrex::ParallelFor
        (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
         {
           adv_fab(i,j,k,n) -= sideCoeffy *
             (sol_fab(i,j,k,n) - sol_fab(i,j-1,k,n));
         });
    }
    else
    {
      amrex::ParallelFor
        (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
         {
           adv_fab(i,j,k,n) -= sideCoeffy *
             (sol_fab(i,j+1,k,n) - sol_fab(i,j,k,n));
         });
    }
  }
}


// Assumes ghots cells are already filled, adds result to diff_mf
void ComputeDiffusion(MultiFab& sol, MultiFab& diff_mf, MultiFab& fx_mf,
                      MultiFab& fy_mf, Geometry& geom,
                      Real diffCoeffx, Real diffCoeffy)
{
  ComputeDiffFlux(sol, fx_mf, fy_mf, geom, diffCoeffx, diffCoeffy);
  ComputeDivergence(diff_mf, fx_mf, fy_mf, geom);
}


// Assumes ghost cells already filled, overwrites fx_mf and fy_mf MultiFabs
void ComputeDiffFlux(MultiFab& sol_mf, MultiFab& fx_mf, MultiFab& fy_mf,
                     Geometry& geom, Real diffCoeffx, Real diffCoeffy)
{
  const auto dx = geom.CellSize();
  Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
  Real dyInv = 1.0 / dx[1]; // assume same over entire mesh
  Real coeffX = diffCoeffx * dxInv;
  Real coeffY = diffCoeffy * dyInv;

  for (MFIter mfi(sol_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    Array4<Real> const& sol = sol_mf[mfi].array();
    Array4<Real> const& fx = fx_mf[mfi].array();
    Array4<Real> const& fy = fy_mf[mfi].array();

    // x-flux
    amrex::ParallelFor
      (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
       {
         // always use zero component for flux
         fx(i,j,k,0) = coeffX * (sol(i,j,k,n) - sol(i-1,j,k,n));
       });

    // y-flux
    amrex::ParallelFor
      (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
       {
         // always use zero component for flux
         fy(i,j,k,0) = coeffY * (sol(i,j,k,n) - sol(i,j-1,k,n));
       });
  }
}


// Assumes ghost cells already filled, adds result to div_mf MultiFab
void ComputeDivergence(MultiFab& div_mf, MultiFab& fx_mf,
                       MultiFab& fy_mf, Geometry& geom)
{
  const auto dx = geom.CellSize();
  Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
  Real dyInv = 1.0 / dx[1]; // assume same over entire mesh

  for (MFIter mfi(div_mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    Array4<Real> const& div = div_mf[mfi].array();
    Array4<Real> const& fx = fx_mf[mfi].array();
    Array4<Real> const& fy = fy_mf[mfi].array();

    amrex::ParallelFor
      (bx, 1, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
       {
         // always use zero component for flux
         div(i,j,k,n) += (dxInv * (fx(i+1,j,k,0) - fx(i,j,k,0)) +
                          dyInv * (fy(i,j+1,k,0) - fy(i,j,k,0)));
       });
  }
}


// -----------------------------------------------------------------------------
// Preconditioner
// -----------------------------------------------------------------------------


int precondition_solve(amrex::Real tn, N_Vector u, N_Vector fu, N_Vector r,
                       N_Vector z, amrex::Real gamma, amrex::Real delta, int lr,
                       void *user_data)
{
  ProblemData *prob_data = (ProblemData*) user_data;

  auto geom = *(prob_data->geom);
  auto grid = *(prob_data->grid);
  auto dmap = *(prob_data->dmap);
  auto& acoef = *(prob_data->acoef);
  auto& bcoef = *(prob_data->acoef);

  MultiFab* solution = amrex::sundials::N_VGetVectorPointer_MultiFab(z);
  MultiFab* rhs = amrex::sundials::N_VGetVectorPointer_MultiFab(r);

  LPInfo info;
  info.setAgglomeration(prob_data->mg_agglomeration);
  info.setConsolidation(prob_data->mg_consolidation);
  info.setMaxCoarseningLevel(prob_data->mg_max_coarsening_level);

  const Real tol_rel = 1.e-10;
  const Real tol_abs = 0.0;

  const int nlevels = 1;

  const Real ascalar = 1.0;
  const Real bscalar = gamma;

  MLABecLaplacian mlabec({geom}, {grid}, {dmap}, info);

  mlabec.setMaxOrder(prob_data->mg_linop_maxorder);

  // Set periodic BC
  mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                   LinOpBCType::Periodic,
                                   LinOpBCType::Periodic)},
    {AMREX_D_DECL(LinOpBCType::Periodic,
                  LinOpBCType::Periodic,
                  LinOpBCType::Periodic)});

  mlabec.setLevelBC(0, nullptr);

  mlabec.setScalars(ascalar, bscalar);

  mlabec.setACoeffs(0, acoef);

  Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
  {
    const BoxArray& ba = amrex::convert(bcoef.boxArray(),
                                        IntVect::TheDimensionVector(idim));
    face_bcoef[idim].define(ba, bcoef.DistributionMap(), 1, 0);

    switch (idim)
    {
    case 0:
      face_bcoef[idim] = prob_data->diffCoeffx;
    case 1:
      face_bcoef[idim] = prob_data->diffCoeffy;
    }
  }

  mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

  MLMG mlmg(mlabec);
  mlmg.setMaxIter(prob_data->mg_max_iter);
  mlmg.setMaxFmgIter(prob_data->mg_max_fmg_iter);
  mlmg.setVerbose(prob_data->mg_verbose);
  mlmg.setBottomVerbose(prob_data->mg_bottom_verbose);
#ifdef AMREX_USE_HYPRE
  if (prob_data->mg_use_hypre)
  {
    mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    if (prob_data->mg_hypre_interface == 1)
      mlmg.setHypreInterface(amrex::Hypre::Interface::structed);
    else if (prob_data->mg_hypre_interface == 2)
      mlmg.setHypreInterface(amrex::Hypre::Interface::semi_structed);
    else
      mlmg.setHypreInterface(amrex::Hypre::Interface::ij);
  }
#endif
#ifdef AMREX_USE_PETSC
  if (prob_data->mg_use_petsc)
  {
    mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
  }
#endif

  mlmg.solve({solution}, {rhs}, prob_data->mg_tol_rel, tol_abs);

  return 0;
}
