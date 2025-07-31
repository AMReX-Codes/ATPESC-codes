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
 * Implementation file for 'intermediate' integration of 2D Advection-Diffusion
 * example problem using ARKODE. This treats the diffusion term of the problem
 * implicitly, but allows advection to be either implicit (via DIRK methods) or
 * explicit (via ARK-IMEX methods). The implicit portion is solved using a
 * Newton-Krylov method with unpreconditioned GMRES as the linear solver. This
 * program allows for either fixed time-step sizes or temporal adaptivity.
 * ---------------------------------------------------------------------------*/

#include "HandsOn.hpp"
#include <arkode/arkode_arkstep.h>
#include <sunlinsol/sunlinsol_spgmr.h>

using namespace amrex;

// -----------------------------------------------------------------------------
// Evolve the problem with ARKODE using an implicit or IMEX method
// -----------------------------------------------------------------------------

void ComputeSolution(N_Vector nv_sol, ProblemOpt* prob_opt,
                     ProblemData* prob_data)
{
  BL_PROFILE("ComputeSolution()");

  // Extract problem data and options
  Geometry* geom         = prob_data->geom;
  int       plot_int     = prob_opt->plot_int;
  int       arkode_order = prob_opt->arkode_order;
  int       nls_max_iter = prob_opt->nls_max_iter;
  int       ls_max_iter  = prob_opt->ls_max_iter;
  int       rhs_adv      = prob_opt->rhs_adv;
  Real      rtol         = prob_opt->rtol;
  Real      atol         = prob_opt->atol;
  Real      fixed_dt     = prob_opt->fixed_dt;
  Real      tfinal       = prob_opt->tfinal;
  Real      dtout        = prob_opt->dtout;
  int       max_steps    = prob_opt->max_steps;
  int       write_diag   = prob_opt->write_diag;

  // initial time, number of outputs, and error flag
  Real time = 0.0;
  int  nout = ceil(tfinal/dtout);
  int  ier  = 0;

  // Write a plotfile of the initial data
  if (plot_int > 0)
  {
    const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
    MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
    WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, time, 0);
  }

  // Create the SUNDIALS context
  SUNContext sunctx = nullptr;
  SUNContext_Create(ParallelContext::CommunicatorAll(), &sunctx);

  // Create the ARK stepper                    ***** UPDATED FROM HandsOn1 *****
  void* arkode_mem = nullptr;
  if (rhs_adv)
  {
    // explicit advection and implicit diffusion
    arkode_mem = ARKStepCreate(ComputeRhsAdv, ComputeRhsDiff, time, nv_sol,
                               sunctx);
  }
  else
  {
    // implicit advection and diffusion
    arkode_mem = ARKStepCreate(nullptr, ComputeRhsAdvDiff, time, nv_sol,
                               sunctx);
  }

  // Attach the user data structure
  ARKodeSetUserData(arkode_mem, prob_data);

  // Set the method order
  ARKodeSetOrder(arkode_mem, arkode_order);

  // Set the time step size or integration tolerances
  if (fixed_dt > 0.0)
    ARKodeSetFixedStep(arkode_mem, fixed_dt);
  else
    ARKodeSStolerances(arkode_mem, atol, rtol);

  // Set the max number of steps between outputs
  ARKodeSetMaxNumSteps(arkode_mem, max_steps);

  // Set logging file
  if (write_diag)
  {
    SUNLogger logger;
    ier = SUNContext_GetLogger(sunctx, &logger);
    if (ier != SUN_SUCCESS)
    {
      amrex::Print() << "Getting the logger failed" << std::endl;
      return;
    }
    SUNLogger_SetInfoFilename(logger, "HandsOn2.log");
  }

  // Create and attach GMRES linear solver     ***** UPDATED FROM HandsOn1 *****
  SUNLinearSolver LS = SUNLinSol_SPGMR(nv_sol, SUN_PREC_NONE, ls_max_iter, sunctx);
  ier = ARKodeSetLinearSolver(arkode_mem, LS, nullptr);
  if (ier != ARKLS_SUCCESS)
  {
    amrex::Print() << "Creating linear solver failed" << std::endl;
    return;
  }

  // Set max number of nonlinear iterations    ***** UPDATED FROM HandsOn1 *****
  ier = ARKodeSetMaxNonlinIters(arkode_mem, nls_max_iter);
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
    BL_PROFILE_VAR("ARKodeEvolve()", pevolve);
    ier = ARKodeEvolve(arkode_mem, tout, nv_sol, &tret, ARK_NORMAL);
    BL_PROFILE_VAR_STOP(pevolve);
    if (ier < 0)
    {
      amrex::Print() << "Error in ARKodeEvolve" << std::endl;
      return;
    }

    // Get integration stats                   ***** UPDATED FROM HandsOn1 *****
    long int nfe_evals, nfi_evals;
    ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe_evals);
    ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi_evals);
    if (nfe_evals > 0)
      amrex::Print() << "t = " << std::setw(5) << tret
                     << "  explicit evals = " << std::setw(7) << nfe_evals
                     << "  implicit evals = " << std::setw(7) << nfi_evals
                     << std::endl;
    else
      amrex::Print() << "t = " << std::setw(5) << tret
                     << "  RHS evals = " << std::setw(7) << nfi_evals
                     << std::endl;

    // Write output
    if (plot_int > 0)
    {
      const std::string& pltfile = amrex::Concatenate("plt", iout+1, 5);
      MultiFab* sol = amrex::sundials::N_VGetVectorPointer_MultiFab(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, tret, iout+1);
    }

    // Update output time
    tout += dtout;
    if (tout > tfinal) tout = tfinal;
  }

  // Output final solution statistics          ***** UPDATED FROM HandsOn1 *****
  long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nlcf, nni, ncfn, netf;
  nst = nst_a = nfe = nfi = nsetups = nli = nJv = nlcf = nni = ncfn = netf = 0;
  ARKodeGetNumSteps(arkode_mem, &nst);
  ARKodeGetNumStepAttempts(arkode_mem, &nst_a);
  ARKodeGetNumRhsEvals(arkode_mem, 0, &nfe);
  ARKodeGetNumRhsEvals(arkode_mem, 1, &nfi);
  ARKodeGetNumErrTestFails(arkode_mem, &netf);
  ARKodeGetNumNonlinSolvIters(arkode_mem, &nni);
  ARKodeGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
  ARKodeGetNumLinSolvSetups(arkode_mem, &nsetups);
  ARKodeGetNumLinIters(arkode_mem, &nli);
  ARKodeGetNumJtimesEvals(arkode_mem, &nJv);
  ARKodeGetNumLinConvFails(arkode_mem, &nlcf);

  amrex::Print() << "\nFinal Solver Statistics:\n"
                 << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n";
  if (nfe > 0)
    amrex::Print() << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n";
  else
    amrex::Print() << "   Total RHS evals = " << nfi << "\n";
  amrex::Print() << "   Total number of nonlinear iterations = " << nni << "\n"
                 << "   Total number of nonlinear solver convergence failures = " << ncfn << "\n"
                 << "   Total number of error test failures = " << netf << "\n";
  amrex::Print() << "   Total linear solver setups = " << nsetups << "\n"
                 << "   Total linear iterations = " << nli << "\n"
                 << "   Total number of Jacobian-vector products = " << nJv << "\n"
                 << "   Total number of linear solver convergence failures = " << nlcf << "\n";

  SUNContext_Free(&sunctx);

  return;
}


// -----------------------------------------------------------------------------
// Print help message with relevant options for Hands-on 2
// -----------------------------------------------------------------------------


void PrintHelp()
{
  amrex::Print()
    << std:: endl
    << "Usage: HandsOn2.exe [fname] [options]" << std::endl
    << "Options:" << std::endl
    << "  help=1" << std::endl
    << "    Print this help message and exit." << std::endl
    << "  plot_int=<int>" << std::endl
    << "    enable (1) or disable (0) plots [default=0]." << std::endl
    << "  arkode_order=<int>" << std::endl
    << "    ARKStep method order [default=4]." << std::endl
    << "  nls_max_iter=<int>" << std::endl
    << "    maximum number of nonlinear iterations [default=3]." << std::endl
    << "  ls_max_iter=<int>" << std::endl
    << "    maximum number of linear iterations [default=5]." << std::endl
    << "  rhs_adv=<int>" << std::endl
    << "    treat advection implicitly (0) or explicitly (1) [default=1]." << std::endl
    << "  fixed_dt=<float>" << std::endl
    << "    use a fixed time step size (if value > 0.0) [default=-1.0]." << std::endl
    << "  rtol=<float>" << std::endl
    << "    relative tolerance for time step adaptivity [default=1e-4]." << std::endl
    << "  atol=<float>" << std::endl
    << "    absolute tolerance for time step adaptivity [default=1e-9]." << std::endl
    << "  tfinal=<float>" << std::endl
    << "    final integration time [default=1e4]." << std::endl
    << "  dtout=<float>" << std::endl
    << "    time between outputs [default=tfinal]." << std::endl
    << "  max_steps=<int>" << std::endl
    << "    maximum number of internal steps between outputs [default=10000]." << std::endl
    << "  write_diag=<int>" << std::endl
    << "    output ARKStep time step adaptivity diagnostics to a file [default=1]." << std::endl
    << "  n_cell=<int>" << std::endl
    << "    number of cells on each side of the square domain [default=128]." << std::endl
    << "  max_grid_size=<int>" << std::endl
    << "    max size of boxes in box array [default=64]." << std::endl
    << "  advCoeffx=<float>" << std::endl
    << "    advection speed in the x-direction [default=5e-4]." << std::endl
    << "  advCoeffy=<float>" << std::endl
    << "    advection speed in the y-direction [default=2.5e-4]." << std::endl
    << "  diffCoeffx=<float>" << std::endl
    << "    diffusion coefficient in the x-direction [default=1e-6]." << std::endl
    << "  diffCoeffy=<float>" << std::endl
    << "    diffusion coefficient in the y-direction [default=1e-6]." << std::endl << std::endl
    << "If a file name 'fname' is provided, it will be parsed for each of the above" << std::endl
    << "options.  If an option is specified in both the input file and on the" << std::endl
    << "command line, then the command line option takes precedence." << std::endl << std::endl;
  return;
}


// -----------------------------------------------------------------------------
// Print relevant problem setup options for Hands-on 1
// -----------------------------------------------------------------------------


void PrintSetup(ProblemOpt& prob_opt, ProblemData& prob_data)
{
  // Ouput problem options and parameters
  amrex::Print()
    << "n_cell        = " << prob_data.n_cell        << std::endl
    << "max_grid_size = " << prob_data.max_grid_size << std::endl
    << "plot_int      = " << prob_opt.plot_int       << std::endl
    << "arkode_order  = " << prob_opt.arkode_order   << std::endl;
  if (prob_opt.rhs_adv)
    amrex::Print()
      << "ImEx treatment (implicit diffusion, explicit advection)" << std::endl;
  else
    amrex::Print()
      << "fully implicit treatment" << std::endl;
  if (prob_opt.fixed_dt > 0.0)
    amrex::Print()
      << "fixed_dt      = " << prob_opt.fixed_dt << std::endl;
  else
    amrex::Print()
      << "rtol          = " << prob_opt.rtol << std::endl
      << "atol          = " << prob_opt.atol << std::endl;
  amrex::Print()
    << "tfinal        = " << prob_opt.tfinal      << std::endl
    << "dtout         = " << prob_opt.dtout       << std::endl
    << "write_diag    = " << prob_opt.write_diag  << std::endl
    << "advCoeffx     = " << prob_data.advCoeffx  << std::endl
    << "advCoeffy     = " << prob_data.advCoeffy  << std::endl
    << "diffCoeffx    = " << prob_data.diffCoeffx << std::endl
    << "diffCoeffy    = " << prob_data.diffCoeffy << std::endl;
  amrex::Print()
    << "Newton nonlinear solver:" << std::endl
    << "  max_iter    = " << prob_opt.nls_max_iter << std::endl
    << "  ls_max_iter = " << prob_opt.ls_max_iter  << std::endl;
  return;
}
