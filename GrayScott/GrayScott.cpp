#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cvode/cvode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "GrayScott.h"
#include "DiffOp.h"
#include "Reactions.h"

#include "NVector_Multifab.h"

using namespace amrex;

int main(int argc, char* argv[])
{
   amrex::Initialize(argc,argv);

   DoProblem();

   amrex::Finalize();
   return 0;
}

void DoProblem()
{
   // What time is it now?  We'll use this to compute total run time.
   Real strt_time = amrex::second();

   // Set problem data and options
   ProblemData prob_data;
   ProblemOpt  prob_opt;
   ParseInputs(prob_opt, prob_data);

   // Make BoxArray and Geometry
   BoxArray ba;
   Geometry geom;
   SetUpGeometry(ba, geom, prob_opt, prob_data);

   // How Boxes are distrubuted among MPI processes
   DistributionMapping dm(ba);

   // Allocate the solution MultiFab
   int nGhost = 1;  // number of ghost cells for each array
   int nComp  = 2;  // number of components for each array
   MultiFab sol(ba, dm, nComp, nGhost);

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
   sunindextype length = nComp * prob_opt.n_cell * prob_opt.n_cell;
   N_Vector nv_sol     = N_VMake_Multifab(length, &sol);

   // Set the initial condition
   FillInitConds2D(sol, geom);

   // Integrate in time
   switch (prob_opt.stepper)
   {
   case 0:
      ComputeSolutionCV(nv_sol, &prob_opt, &prob_data);
      break;
   case 1:
      ComputeSolutionARK(nv_sol, &prob_opt, &prob_data);
      break;
   case 2:
      ComputeSolutionMRI(nv_sol, &prob_opt, &prob_data);
      break;
   default:
      amrex::Print() << "Invalid stepper option" << std::endl;
      return;
   }

   // Call the timer again and compute the maximum difference between the start
   // time and stop time over all processors
   Real stop_time = amrex::second() - strt_time;
   const int IOProc = ParallelDescriptor::IOProcessorNumber();
   ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

   // Tell the I/O Processor to write out the "run time"
   amrex::Print() << "Run time = " << stop_time << std::endl;
}

int ComputeRhsAdv(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Geometry* geom = prob_data->geom;
   Real advCoeffU = prob_data->advCoeffU;
   Real advCoeffV = prob_data->advCoeffV;

   // clear the RHS
   *rhs = 0.0;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // compute advection of u and v
   ComputeAdvection(*sol, *rhs, *geom, 0, advCoeffU);
   ComputeAdvection(*sol, *rhs, *geom, 1, advCoeffV);

   return 0;
}

int ComputeRhsDiff(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Geometry* geom = prob_data->geom;
   Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
   Real diffCoeffU = prob_data->diffCoeffU;
   Real diffCoeffV = prob_data->diffCoeffV;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // clear the RHS
   *rhs = 0.0;

   // compute diffusion of u and v
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 0,
                    diffCoeffU, diffCoeffU);
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 1,
                    diffCoeffV, diffCoeffV);

   return 0;
}

int ComputeRhsReact(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Real A = prob_data->A;
   Real B = prob_data->B;

   // clear the RHS
   *rhs = 0.0;

   // compute reaction term
   ComputeReactionsGS(*sol, *rhs, A, B);

   return 0;
}

int ComputeRhsAdvDiff(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Geometry* geom = prob_data->geom;
   Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
   Real advCoeffU  = prob_data->advCoeffU;
   Real advCoeffV  = prob_data->advCoeffV;
   Real diffCoeffU = prob_data->diffCoeffU;
   Real diffCoeffV = prob_data->diffCoeffV;

   // clear the RHS
   *rhs = 0.0;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // compute advection of u and v
   ComputeAdvection(*sol, *rhs, *geom, 0, advCoeffU);
   ComputeAdvection(*sol, *rhs, *geom, 1, advCoeffV);

   // compute diffusion of u and v
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 0,
                    diffCoeffU, diffCoeffU);
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 1,
                    diffCoeffV, diffCoeffV);

   return 0;
}

int ComputeRhsAdvReact(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Geometry* geom = prob_data->geom;
   Real advCoeffU = prob_data->advCoeffU;
   Real advCoeffV = prob_data->advCoeffV;
   Real A         = prob_data->A;
   Real B         = prob_data->B;

   // clear the RHS
   *rhs = 0.0;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // compute advection of u and v
   ComputeAdvection(*sol, *rhs, *geom, 0, advCoeffU);
   ComputeAdvection(*sol, *rhs, *geom, 1, advCoeffV);

   // compute reaction term
   ComputeReactionsGS(*sol, *rhs, A, B);

   return 0;
}

int ComputeRhsDiffReact(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Geometry* geom = prob_data->geom;
   Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
   Real diffCoeffU = prob_data->diffCoeffU;
   Real diffCoeffV = prob_data->diffCoeffV;
   Real A          = prob_data->A;
   Real B          = prob_data->B;

   // clear the RHS
   *rhs = 0.0;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // compute diffusion of u and v
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 0,
                    diffCoeffU, diffCoeffU);
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 1,
                    diffCoeffV, diffCoeffV);

   // compute reaction term
   ComputeReactionsGS(*sol, *rhs, A, B);

   return 0;
}

int ComputeRhsAdvDiffReact(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Geometry* geom = prob_data->geom;
   Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
   Real advCoeffU  = prob_data->advCoeffU;
   Real advCoeffV  = prob_data->advCoeffV;
   Real diffCoeffU = prob_data->diffCoeffU;
   Real diffCoeffV = prob_data->diffCoeffV;
   Real A          = prob_data->A;
   Real B          = prob_data->B;

   // clear the RHS
   *rhs = 0.0;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // compute advection of u and v
   ComputeAdvection(*sol, *rhs, *geom, 0, advCoeffU);
   ComputeAdvection(*sol, *rhs, *geom, 1, advCoeffV);

   // compute diffusion of u and v
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 0,
                    diffCoeffU, diffCoeffU);
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 1,
                    diffCoeffV, diffCoeffV);

   // compute reaction term
   ComputeReactionsGS(*sol, *rhs, A, B);

   return 0;
}

void FillInitConds2D(MultiFab& sol, const Geometry& geom)
{
   const auto dx = geom.CellSize();
   const auto prob_lo = geom.ProbLo();
   const auto prob_hi = geom.ProbHi();
   for (MFIter mfi(sol); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& fab = sol.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);
      for (int j = lo.y; j <= hi.y; ++j) {
         Real y = prob_lo[1] + (((Real) j) + 0.5) * dx[1];

         for (int i = lo.x; i <= hi.x; ++i) {
            fab(i,j,0,0) = 0.0;
            fab(i,j,0,1) = 0.0;

            Real x = prob_lo[0] + (((Real) i) + 0.5) * dx[0];
            if (x>=1.14 && x<=1.33 && y>=1.14 && y<=1.33)
            {
               fab(i,j,0,0) = 0.5;
               fab(i,j,0,1) = 0.25;
            }
         }
      }
   }
}

void ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data)
{
   // ParmParse is way of reading inputs from the inputs file
   ParmParse pp;

   // The number of cells on each side of a square domain.
   int n_cell = 256;
   pp.query("n_cell", n_cell);
   prob_opt.n_cell = n_cell;

   // The domain is broken into boxes of size max_grid_size
   int max_grid_size = 64;
   pp.query("max_grid_size", max_grid_size);
   prob_opt.max_grid_size = max_grid_size;

   // Enable (>0) or disable (<0) writing output files
   int plot_int = -1; // plots off
   pp.query("plot_int", plot_int);
   prob_opt.plot_int = plot_int;

   // Specify which integration method to use
   // 0 = CVODE
   // 1 = ARKStep
   // 2 = MRIStep
   int stepper = 0;
   pp.query("stepper", stepper);
   prob_opt.stepper = stepper;

   // Specify which CVODE method to use
   int cvode_method = 0; // BDF
   pp.query("cvode_method", cvode_method);
   prob_opt.cvode_method = cvode_method;

   // Specify the ARKode method order
   int arkode_order = 4; // 4th order
   pp.query("arkode_order", arkode_order);
   prob_opt.arkode_order = arkode_order;

   // Specify the nonlinear solver
   int nls_method = 0; // Newton
   pp.query("nls_method", nls_method);
   prob_opt.nls_method = nls_method;

   // Specify the max number of nonlinear iterations
   int nls_max_iter = 3;
   pp.query("nls_max_iter", nls_max_iter);
   prob_opt.nls_max_iter = nls_max_iter;

   // Specify the number of fixed point acceleration vectors
   int nls_fp_acc = 0; // no acceleration
   pp.query("nls_fp_acc", nls_fp_acc);
   prob_opt.nls_fp_acc = nls_fp_acc;

   // Specify the max number of linear iterations
   int ls_max_iter = 5;
   pp.query("ls_max_iter", ls_max_iter);
   prob_opt.ls_max_iter = ls_max_iter;

   // Specify relative and absolute tolerances
   Real rtol = 1.0e-4;
   Real atol = 1.0e-9;
   pp.query("rtol", rtol);
   pp.query("atol", atol);
   prob_opt.rtol = rtol;
   prob_opt.atol = atol;

   // Specify final time for integration
   Real tfinal = 1.0e4;
   pp.query("tfinal", tfinal);
   prob_opt.tfinal = tfinal;

   // Specify output frequency
   Real dtout = tfinal;
   pp.query("dtout", dtout);
   prob_opt.dtout = dtout;

   // Advection coefficients
   Real advCoeffU = 5.0e-4;
   Real advCoeffV = 5.0e-4;
   pp.query("advCoeffU", advCoeffU);
   pp.query("advCoeffV", advCoeffV);
   prob_data.advCoeffU = advCoeffU;
   prob_data.advCoeffV = advCoeffV;

   // Diffusion coefficients
   Real diffCoeffU = 2.0e-5;
   Real diffCoeffV = 1.0e-5;
   pp.query("diffCoeffU", diffCoeffU);
   pp.query("diffCoeffV", diffCoeffV);
   prob_data.diffCoeffU = diffCoeffU;
   prob_data.diffCoeffV = diffCoeffV;

   // Gray-Scott reaction parameters
   Real A = 0.04;
   Real B = 0.06;
   pp.query("A", A);
   pp.query("B", B);
   prob_data.A = A;
   prob_data.B = B;

   // Output problem options and parameters
   amrex::Print()
      << "n_cell        = " << n_cell        << std::endl
      << "max_grid_size = " << max_grid_size << std::endl
      << "plot_int      = " << plot_int      << std::endl
      << "stepper       = " << stepper       << std::endl;

   if (stepper == 0)
      amrex::Print()
         << "cvode_method  = " << cvode_method  << std::endl;
   else
      amrex::Print()
         << "arkode_order  = " << arkode_order << std::endl;

   amrex::Print()
      << "rtol          = " << rtol          << std::endl
      << "atol          = " << atol          << std::endl
      << "tfinal        = " << tfinal        << std::endl
      << "dtout         = " << dtout         << std::endl;

   amrex::Print()
      << "advCoeffU     = " << advCoeffU << std::endl
      << "advCoeffV     = " << advCoeffV << std::endl;

   amrex::Print()
      << "diffCoeffU    = " << diffCoeffU << std::endl
      << "diffCoeffV    = " << diffCoeffV << std::endl;

   amrex::Print()
      << "A             = " << A << std::endl
      << "B             = " << B << std::endl;

}

void SetUpGeometry(BoxArray& ba, Geometry& geom,
                   ProblemOpt& prob_opt, ProblemData& prob_data)
{
   // Extract problem options
   int n_cell = prob_opt.n_cell;
   int max_grid_size = prob_opt.max_grid_size;

   IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
   IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
   Box domain(dom_lo, dom_hi); // cell-centered

   // Initialize the boxarray "ba" from the single box "domain"
   ba.define(domain);

   // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a
   // direction
   ba.maxSize(max_grid_size);

   // This defines the physical box, [-1,1] in each direction.
   RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                    {AMREX_D_DECL(2.5, 2.5, 2.5)});

   // This defines a Geometry object
   Vector<int> is_periodic(AMREX_SPACEDIM, 1);  // periodic in all direction
   geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());

   prob_data.geom = &geom;
}

void ComputeSolutionCV(N_Vector nv_sol, ProblemOpt* prob_opt,
                       ProblemData* prob_data)
{
   // Extract problem data and options
   Geometry* geom         = prob_data->geom;
   int       plot_int     = prob_opt->plot_int;
   int       cvode_method = prob_opt->cvode_method;
   int       nls_method   = prob_opt->nls_method;
   int       nls_max_iter = prob_opt->nls_max_iter;
   int       nls_fp_acc   = prob_opt->nls_fp_acc;
   int       ls_max_iter  = prob_opt->ls_max_iter;
   Real      rtol         = prob_opt->rtol;
   Real      atol         = prob_opt->atol;
   Real      tfinal       = prob_opt->tfinal;
   Real      dtout        = prob_opt->dtout;

   // initial time, number of outputs, and error flag
   Real time = 0.0;
   int  nout = ceil(tfinal/dtout);
   int  ier  = 0;

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"}, *geom, time, 0);
   }

   // Create CVODE memory
   void* cvode_mem = NULL;
   if (cvode_method == 0)
      cvode_mem = CVodeCreate(CV_BDF);
   else
      cvode_mem = CVodeCreate(CV_ADAMS);

   // Initialize CVODE implicit adv-diff-react
   CVodeInit(cvode_mem, ComputeRhsAdvDiffReact, time, nv_sol);

   // Set CVODE options
   CVodeSetUserData(cvode_mem, prob_data);
   CVodeSStolerances(cvode_mem, atol, rtol);

   // Attach nonlinear/linear solvers as needed
   if (nls_method == 0)
   {
      // Create and attach GMRES linear solver for Newton
      SUNLinearSolver LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, ls_max_iter);
      ier = CVodeSetLinearSolver(cvode_mem, LS, NULL);
      if (ier != CVLS_SUCCESS)
      {
         amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
         return;
      }
   }
   else
   {
      // Create and attach fixed point solver
      SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(nv_sol, nls_fp_acc);
      ier = CVodeSetNonlinearSolver(cvode_mem, NLS);
      if (ier != CV_SUCCESS)
      {
         amrex::Print() << "Creation of nonlinear solver unsuccessful" << std::endl;
         return;
      }
   }

   // Set max number of nonlinear iterations
   ier = CVodeSetMaxNonlinIters(cvode_mem, nls_max_iter);
   if (ier != CV_SUCCESS)
   {
      amrex::Print() << "Error setting max number of nonlinear iterations" << std::endl;
      return;
   }

   // Advance the solution in time
   Real tout = time + dtout; // first output time
   Real tret;                // return time
   for (int iout=0; iout < nout; iout++)
   {
      ier = CVode(cvode_mem, tout, nv_sol, &tret, CV_NORMAL);
      if (ier < 0)
      {
         amrex::Print() << "Error in CVODE" << std::endl;
         return;
      }

      // Get integration stats
      long nf_evals;
      CVodeGetNumRhsEvals(cvode_mem, &nf_evals);
      amrex::Print() << "t = " << std::setw(5) << tret
                     << "  rhs evals = " << std::setw(7) << nf_evals
                     << std::endl;

      // Write output
      if (plot_int > 0)
      {
         const std::string& pltfile = amrex::Concatenate("plt", iout+1, 5);
         MultiFab* sol = NV_MFAB(nv_sol);
         WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"}, *geom, tret,
                                  iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }
}

void ComputeSolutionARK(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data)
{
   // Extract problem data and options
   Geometry* geom         = prob_data->geom;
   int       plot_int     = prob_opt->plot_int;
   int       arkode_order = prob_opt->arkode_order;
   int       nls_method   = prob_opt->nls_method;
   int       nls_max_iter = prob_opt->nls_max_iter;
   int       nls_fp_acc   = prob_opt->nls_fp_acc;
   int       ls_max_iter  = prob_opt->ls_max_iter;
   Real      rtol         = prob_opt->rtol;
   Real      atol         = prob_opt->atol;
   Real      tfinal       = prob_opt->tfinal;
   Real      dtout        = prob_opt->dtout;

   // initial time, number of outputs, and error flag
   Real time = 0.0;
   int  nout = ceil(tfinal/dtout);
   int  ier  = 0;

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"}, *geom, time, 0);
   }

   // Create the ARK stepper
   void* arkode_mem = NULL;

   // explicit diffusion, implicit reactions
   // void* arkode_mem = ARKStepCreate(ComputeRhsDiff, ComputeRhsReact,
   //                                  time, nv_sol);

   // explicit reactions, implicit diffusion
   arkode_mem = ARKStepCreate(ComputeRhsReact, ComputeRhsDiff,
                              time, nv_sol);

   // Set ARKStep options
   ARKStepSetUserData(arkode_mem, prob_data);
   ARKStepSStolerances(arkode_mem, atol, rtol);
   ARKStepSetOrder(arkode_mem, arkode_order);

   // Attach nonlinear/linear solvers as needed
   if (nls_method == 0)
   {
      // Create and attach GMRES linear solver for Newton
      SUNLinearSolver LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, ls_max_iter);
      ier = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
      if (ier != ARKLS_SUCCESS)
      {
         amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
         return;
      }
   }
   else
   {
      // Create and attach GMRES linear solver (if implicit and using Newton)
      SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(nv_sol, nls_fp_acc);
      ier = ARKStepSetNonlinearSolver(arkode_mem, NLS);
      if (ier != ARK_SUCCESS)
      {
         amrex::Print() << "Creation of nonlinear solver unsuccessful" << std::endl;
         return;
      }
   }

   // Set max number of nonlinear iterations
   ier = ARKStepSetMaxNonlinIters(arkode_mem, nls_max_iter);
   if (ier != ARK_SUCCESS)
   {
      amrex::Print() << "Error setting max number of nonlinear iterations" << std::endl;
      return;
   }

   // Advance the solution in time
   Real tout = time + dtout; // first output time
   Real tret;                // return time
   for (int iout=0; iout < nout; iout++)
   {
      ier = ARKStepEvolve(arkode_mem, tout, nv_sol, &tret, ARK_NORMAL);
      if (ier < 0)
      {
         amrex::Print() << "Error in ARKStepEvolve" << std::endl;
         return;
      }

      // Get integration stats
      long nfe_evals, nfi_evals;
      ARKStepGetNumRhsEvals(arkode_mem, &nfe_evals, &nfi_evals);
      amrex::Print() << "t = " << std::setw(5) << tret
                     << "  explicit evals = " << std::setw(7) << nfe_evals
                     << "  implicit evals = " << std::setw(7) << nfi_evals
                     << std::endl;

      // Write output
      if (plot_int > 0)
      {
         const std::string& pltfile = amrex::Concatenate("plt", iout+1, 5);
         MultiFab* sol = NV_MFAB(nv_sol);
         WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"}, *geom, tret,
                                  iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }
}

void ComputeSolutionMRI(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data)
{
   // Extract problem data and options
   Geometry* geom         = prob_data->geom;
   int       plot_int     = prob_opt->plot_int;
   int       nls_method   = prob_opt->nls_method;
   int       nls_max_iter = prob_opt->nls_max_iter;
   int       nls_fp_acc   = prob_opt->nls_fp_acc;
   int       ls_max_iter  = prob_opt->ls_max_iter;
   Real      rtol         = prob_opt->rtol;
   Real      atol         = prob_opt->atol;
   Real      tfinal       = prob_opt->tfinal;
   Real      dtout        = prob_opt->dtout;

   // initial time, number of outputs, and error flag
   Real time = 0.0;
   int  nout = ceil(tfinal/dtout);
   int  ier  = 0;

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"}, *geom, time, 0);
   }

   // Create the MRI stepper
   void* arkode_mem = NULL;

   // explicit slow diffusion, explicit fast reactions
   arkode_mem = MRIStepCreate(ComputeRhsDiff, ComputeRhsReact,
                              time, nv_sol);

   // Set MRIStep options
   ARKStepSetUserData(arkode_mem, prob_data);
   MRIStepSetFixedStep(arkode_mem, 0.5, 0.5);
   MRIStepSetMaxNumSteps(arkode_mem, 500000);

   // Advance the solution in time
   Real tout = time + dtout; // first output time
   Real tret;                // return time
   for (int iout=0; iout < nout; iout++)
   {
      ier = MRIStepEvolve(arkode_mem, tout, nv_sol, &tret, ARK_NORMAL);
      if (ier < 0)
      {
         amrex::Print() << "Error in MRIStepEvolve" << std::endl;
         return;
      }

      // Get integration stats
      long nfs_evals, nff_evals;
      MRIStepGetNumRhsEvals(arkode_mem, &nfs_evals, &nff_evals);
      amrex::Print() << "t = " << std::setw(5) << tret
                     << "  slow evals " << std::setw(7) << nfs_evals
                     << "  fast evals " << std::setw(7) << nff_evals
                     << std::endl;

      // Write output
      if (plot_int > 0)
      {
         const std::string& pltfile = amrex::Concatenate("plt", iout+1, 5);
         MultiFab* sol = NV_MFAB(nv_sol);
         WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"}, *geom, tret,
                                  iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }
}
