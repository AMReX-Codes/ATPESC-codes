/*--------------------------------------------------------------------
  Time Integration and Nonlinear Solvers
  Hands-on Lessons with SUNDIALS + AMReX
  2019 Argonne Training Program in Extreme-Scale Computing

  Authors (alphabetical):
  David Gardner (gardner48@llnl.gov)
  John Loffeld (loffeld1@llnl.gov)
  Daniel Reynolds (reynolds@smu.edu)
  Donald Willcox (dewillcox@lbl.gov)

  --------------------------------------------------------------------
  Implementation file for 'simplest' ARKode integration of 2D
  Advection-Diffusion example problem.  This allows for either
  fixed-step or temporal adaptivity, but requires explicit time
  integration.
  --------------------------------------------------------------------*/

#include "HandsOn.hpp"
#include <arkode/arkode_arkstep.h>

using namespace amrex;

// -----------------------------------------------------------------------------
// Evolve the problem with ARKODE using an explicit method
// -----------------------------------------------------------------------------

void ComputeSolution(N_Vector nv_sol, ProblemOpt* prob_opt,
                     ProblemData* prob_data)
{
  BL_PROFILE("ComputeSolution()");

  // Extract problem data and options
  Geometry* geom         = prob_data->geom;
  int       plot_int     = prob_opt->plot_int;
  int       arkode_order = prob_opt->arkode_order;
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

  // Create the ARK stepper
  void* arkode_mem = ARKStepCreate(ComputeRhsAdvDiff, nullptr, time, nv_sol,
                                   *amrex::sundials::The_Sundials_Context());

  // Attach the user data structure to ARKStep
  ARKStepSetUserData(arkode_mem, prob_data);

  // Set the method order
  ARKStepSetOrder(arkode_mem, arkode_order);

  // Set the time step size or integration tolerances
  if (fixed_dt > 0.0)
    ARKStepSetFixedStep(arkode_mem, fixed_dt);
  else
    ARKStepSStolerances(arkode_mem, atol, rtol);

  // Set the max number of steps between outputs
  ARKStepSetMaxNumSteps(arkode_mem, max_steps);

  // Set file for writing ARKStep diagnostics
  FILE* diagfp = nullptr;
  if (write_diag)
  {
    diagfp = fopen("HandsOn1_diagnostics.txt", "w");
    ARKStepSetDiagnostics(arkode_mem, diagfp);
  }

  // Advance the solution in time
  Real tout = time + dtout; // first output time
  Real tret;                // return time
  for (int iout=0; iout < nout; iout++)
  {
    BL_PROFILE_VAR("ARKStepEvolve()", pevolve);
    ier = ARKStepEvolve(arkode_mem, tout, nv_sol, &tret, ARK_NORMAL);
    BL_PROFILE_VAR_STOP(pevolve);
    if (ier < 0)
    {
      amrex::Print() << "Error in ARKStepEvolve" << std::endl;
      return;
    }

    // Get integration stats
    long nfe_evals, nfi_evals;
    ARKStepGetNumRhsEvals(arkode_mem, &nfe_evals, &nfi_evals);
    amrex::Print() << "t = " << std::setw(5) << tret
                   << "  RHS evals = " << std::setw(7) << nfe_evals
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

  // Output final solution statistics
  long int nst, nst_a, nfe, nfi, netf;
  nst = nst_a = nfe = nfi = netf = 0;
  ARKStepGetNumSteps(arkode_mem, &nst);
  ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
  ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
  ARKStepGetNumErrTestFails(arkode_mem, &netf);
  amrex::Print() << "\nFinal Solver Statistics:\n"
                 << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n"
                 << "   Total RHS evals = " << nfe << "\n"
                 << "   Total number of error test failures = " << netf << "\n";

  // Close diagnostics file
  if (write_diag) fclose(diagfp);

  return;
}


// -----------------------------------------------------------------------------
// Print help message with relevant options for Hands-on 1
// -----------------------------------------------------------------------------


void PrintHelp()
{
  amrex::Print()
    << std:: endl
    << "Usage: HandsOn1.exe [fname] [options]" << std::endl
    << "Options:" << std::endl
    << "  help=1" << std::endl
    << "    Print this help message and exit." << std::endl
    << "  plot_int=<int>" << std::endl
    << "    enable (1) or disable (0) plots [default=0]." << std::endl
    << "  arkode_order=<int>" << std::endl
    << "    ARKStep method order [default=4]." << std::endl
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
  return;
}
