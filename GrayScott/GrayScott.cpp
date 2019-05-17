#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cvode/cvode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_mristep.h>
#include <sunlinsol/sunlinsol_spgmr.h>

#include "GrayScott.h"
#include "DiffOp.h"

#include "NVector_Multifab.h"

using namespace amrex;

int main(int argc, char* argv[])
{
   amrex::Initialize(argc,argv);
   amrex::Print() << "Hello world from AMReX version "
                  << amrex::Version() << "\n";

   DoProblem();

   amrex::Finalize();
   return 0;
}

void DoProblem()
{
   // What time is it now?  We'll use this to compute total run time.
   Real strt_time = amrex::second();


   // --- Problem setup ---
   int n_cell, max_grid_size, stepper, plot_int;
   Real tfinal, dtout;
   GrayScottProblem problem;
   ParseInputs(n_cell, max_grid_size, stepper, tfinal, dtout, plot_int,
               problem);

   // make BoxArray and Geometry
   BoxArray ba;
   Geometry geom;
   SetUpGeometry(ba, geom, problem, n_cell, max_grid_size);


   // --- Create solution vectors ---
   // How Boxes are distrubuted among MPI processes
   DistributionMapping dm(ba);

   // allocate the solution MultiFab
   int nGhost = 1;  // number of ghost cells for each array
   int nComp  = 2;  // number of components for each array
   MultiFab sol(ba, dm, nComp, nGhost);

   // build the flux MultiFabs
   Array<MultiFab, AMREX_SPACEDIM> flux;
   for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
   {
      // flux(dir) has one component, zero ghost cells, and is nodal in
      // direction dir
      BoxArray edge_ba = ba;
      edge_ba.surroundingNodes(dir);
      flux[dir].define(edge_ba, dm, 1, 0);
   }
   problem.flux = &flux;

   // create an N_Vector wrapper for the solution MultiFab
   sunindextype length = nComp * n_cell * n_cell;
   N_Vector nv_sol     = N_VMake_Multifab(length, &sol);

   // set the initial condition
   FillInitConds2D(sol, geom);


   // --- Time advance to end ---
   switch (stepper)
   {
   case 0:
      ComputeSolutionCV(nv_sol, &problem, tfinal, dtout, plot_int);
      break;
   case 1:
      ComputeSolutionARK(nv_sol, &problem, tfinal, dtout, plot_int);
      break;
   case 2:
      ComputeSolutionMRI(nv_sol, &problem, tfinal, dtout, plot_int);
      break;
   default:
      amrex::Print() << "Invalid stepper option" << std::endl;
      return;
   }


   // --- Print time and exit ---
   // Call the timer again and compute the maximum difference between the start
   // time and stop time over all processors
   Real stop_time = amrex::second() - strt_time;
   const int IOProc = ParallelDescriptor::IOProcessorNumber();
   ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

   // Tell the I/O Processor to write out the "run time"
   amrex::Print() << "Run time = " << stop_time << std::endl;
}

void ComputeDiffusion(MultiFab& sol, MultiFab& diffusion,
                      GrayScottProblem& problem)
{
   Geometry* geom = problem.geom;
   Array<MultiFab, AMREX_SPACEDIM>& flux = *(problem.flux);
   Real diffCoeffU = problem.diffCoeffU;
   Real diffCoeffV = problem.diffCoeffV;

   sol.FillBoundary(geom->periodicity());

   ComputeDiffFlux(sol, flux[0], flux[1], *geom, 0, diffCoeffU);
   ComputeDivergence(diffusion, flux[0], flux[1], *geom, 0);

   ComputeDiffFlux(sol, flux[0], flux[1], *geom, 1, diffCoeffV);
   ComputeDivergence(diffusion, flux[0], flux[1], *geom, 1);
}

int ComputeDiffusionNV(realtype t, N_Vector nv_sol, N_Vector nv_diffusion,
                       void* problem)
{
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* diffusion = NV_MFAB(nv_diffusion);
   GrayScottProblem *gs_problem = (GrayScottProblem*) problem;

   ComputeDiffusion(*sol, *diffusion, *gs_problem);

   return 0;
}

void ComputeReactions2D(MultiFab& sol, MultiFab& reactions,
                        GrayScottProblem& problem)
{
   Real A = problem.A;
   Real B = problem.B;

   for (MFIter mfi(sol); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& sol_fab = sol.array(mfi);
      Array4<Real> const& reactions_fab = reactions.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int j = lo.y; j <= hi.y; ++j) {
         for (int i = lo.x; i <= hi.x; ++i) {
            Real temp = sol_fab(i,j,0,0) * sol_fab(i,j,0,1) * sol_fab(i,j,0,1);
            reactions_fab(i,j,0,0) = A * (1.0 - sol_fab(i,j,0,0)) - temp;
            reactions_fab(i,j,0,1) = temp - (A + B) * sol_fab(i,j,0,1);
         }
      }
   }
}

int ComputeReactionsNV(realtype t, N_Vector nv_sol, N_Vector nv_reactions,
                       void* problem)
{
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* reactions = NV_MFAB(nv_reactions);
   GrayScottProblem *gs_problem = (GrayScottProblem*) problem;

   ComputeReactions2D(*sol, *reactions, *gs_problem);

   return 0;
}

void ComputeDiffusionReactions2D(MultiFab& sol, MultiFab& rhs,
                                 GrayScottProblem& problem)
{
   Geometry* geom = problem.geom;
   Array<MultiFab, AMREX_SPACEDIM>& flux = *(problem.flux);
   Real diffCoeffU = problem.diffCoeffU;
   Real diffCoeffV = problem.diffCoeffV;

   // Compute diffusion term
   sol.FillBoundary(geom->periodicity());

   ComputeDiffFlux(sol, flux[0], flux[1], *geom, 0, diffCoeffU);
   ComputeDivergence(rhs, flux[0], flux[1], *geom, 0);

   ComputeDiffFlux(sol, flux[0], flux[1], *geom, 1, diffCoeffV);
   ComputeDivergence(rhs, flux[0], flux[1], *geom, 1);

   // Compute reaction term
   Real A = problem.A;
   Real B = problem.B;

   for (MFIter mfi(sol); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& sol_fab = sol.array(mfi);
      Array4<Real> const& rhs_fab = rhs.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int j = lo.y; j <= hi.y; ++j) {
         for (int i = lo.x; i <= hi.x; ++i) {
            Real temp = sol_fab(i,j,0,0) * sol_fab(i,j,0,1) * sol_fab(i,j,0,1);
            rhs_fab(i,j,0,0) += A * (1.0 - sol_fab(i,j,0,0)) - temp;
            rhs_fab(i,j,0,1) += temp - (A + B) * sol_fab(i,j,0,1);
         }
      }
   }
}

int ComputeDiffusionReactionsNV(realtype t, N_Vector nv_sol, N_Vector nv_rhs,
                                void* problem)
{
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);
   GrayScottProblem *gs_problem = (GrayScottProblem*) problem;

   ComputeDiffusionReactions2D(*sol, *rhs, *gs_problem);

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
            fab(i,j,0,0) = 0.;
            fab(i,j,0,1) = 0.;

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

void ParseInputs(int& n_cell, int& max_grid_size, int& stepper, Real& tfinal,
                 Real& dtout, int& plot_int, GrayScottProblem& problem)
{
   Real diffCoeffU, diffCoeffV, A, B;

   // ParmParse is way of reading inputs from the inputs file
   ParmParse pp;

   // We need to get n_cell from the inputs file - this is the number of cells
   // on each side of a square domain.
   pp.get("n_cell", n_cell);

   // The domain is broken into boxes of size max_grid_size
   pp.get("max_grid_size", max_grid_size);

   // Default plot_int to -1, allow us to set it to something else in the inputs
   // file. If plot_int < 0 then no plot files will be written
   plot_int = -1;
   pp.query("plot_int", plot_int);

   // Specify which integration method to use (defaults to 0)
   // 0 = CVODE
   // 1 = ARKStep
   // 2 = MRIStep
   stepper = 0;
   pp.query("stepper", stepper);

   // Specify final time for integration (default to 1000)
   tfinal = 1.0e3;
   pp.query("tfinal", tfinal);

   // Specify output frequency (default to final time)
   dtout = tfinal;
   pp.query("dtout", dtout);

   // Get Gray-Scott problem coefficients
   pp.query("diffCoeffU", diffCoeffU);
   pp.query("diffCoeffV", diffCoeffV);
   pp.query("A", A);
   pp.query("B", B);

   amrex::Print() << "diffCoeffU = " << diffCoeffU;
   amrex::Print() << "  diffCoeffV = " << diffCoeffV;
   amrex::Print() << "  A = " << A;
   amrex::Print() << "  B = " << B << std::endl;

   problem.diffCoeffU = diffCoeffU;
   problem.diffCoeffV = diffCoeffV;
   problem.A = A;
   problem.B = B;
}

void SetUpGeometry(BoxArray& ba, Geometry& geom,
                   GrayScottProblem& problem,
                   int n_cell, int max_grid_size)
{
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

   problem.geom = &geom;
}

void ComputeSolutionMRI(N_Vector nv_sol, GrayScottProblem* problem,
                        Real tfinal, Real dtout, int plot_int)
{
   Geometry* geom = problem->geom;

   Real time = 0.0;                // time = starting time in the simulation
   int  ier  = 0;                  // error flag
   int  nout = ceil(tfinal/dtout); // number of outputs

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"},
                               *geom, time, 0);
   }

   // Create the MRI stepper
   void* arkode_mem = MRIStepCreate(ComputeDiffusionNV, ComputeReactionsNV,
                                    time, nv_sol);

   // Set MRIStep options
   MRIStepSetFixedStep(arkode_mem, 0.5, 0.5);
   MRIStepSetMaxNumSteps(arkode_mem, 500000);
   MRIStepSetUserData(arkode_mem, problem);

   // Advance the solution in time
   realtype tout = time + dtout; // first output time
   realtype tret;                // return time
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
         WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"},
                                  *geom, tret, iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }
}

void ComputeSolutionARK(N_Vector nv_sol, GrayScottProblem* problem,
                        Real tfinal, Real dtout, int plot_int)
{
   Geometry* geom = problem->geom;

   Real time = 0.0;                // time = starting time in the simulation
   int  ier  = 0;                  // error flag
   int  nout = ceil(tfinal/dtout); // number of outputs

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"},
                               *geom, time, 0);
   }

   // Create the ARK stepper

   // explicit diffusion, implicit reactions
   // void* arkode_mem = ARKStepCreate(ComputeDiffusionNV, ComputeReactionsNV,
   //                                  time, nv_sol);

   // explicit reactions, implicit diffusion
   void* arkode_mem = ARKStepCreate(ComputeReactionsNV, ComputeDiffusionNV,
                                    time, nv_sol);

   // Set ARKStep options
   ARKStepSetFixedStep(arkode_mem, 1.0);
   ARKStepSetMaxNumSteps(arkode_mem, 5000);
   ARKStepSetUserData(arkode_mem, problem);

   // Create and attach GMRES linear solver (without preconditioning)
   SUNLinearSolver LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, 100);
   ier = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
   if (ier != ARKLS_SUCCESS)
   {
      amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
      return;
   }

   // Advance the solution in time
   realtype tout = time + dtout; // first output time
   realtype tret;                // return time
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
         WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"},
                                  *geom, tret, iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }
}

void ComputeSolutionCV(N_Vector nv_sol, GrayScottProblem* problem,
                       Real tfinal, Real dtout, int plot_int)
{
   Geometry* geom = problem->geom;

   Real time = 0.0;                // time = starting time in the simulation
   int  ier  = 0;                  // error flag
   int  nout = ceil(tfinal/dtout); // number of outputs

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"},
                               *geom, time, 0);
   }

   // Create CVODE memory
   void* cvode_mem = CVodeCreate(CV_BDF);
   CVodeInit(cvode_mem, ComputeDiffusionReactionsNV, time, nv_sol);

   // Set CVODE options
   CVodeSStolerances(cvode_mem, 1.0e-4, 1.0e-9);
   CVodeSetMaxNumSteps(cvode_mem, 5000);
   CVodeSetUserData(cvode_mem, problem);

   // Create and attach GMRES linear solver (without preconditioning)
   SUNLinearSolver LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, 100);
   ier = CVodeSetLinearSolver(cvode_mem, LS, NULL);
   if (ier != CVLS_SUCCESS)
   {
      amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
      return;
   }

   // Advance the solution in time
   realtype tout = time + dtout; // first output time
   realtype tret;                // return time
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
         WriteSingleLevelPlotfile(pltfile, *sol, {"u", "v"},
                                  *geom, tret, iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }
}
