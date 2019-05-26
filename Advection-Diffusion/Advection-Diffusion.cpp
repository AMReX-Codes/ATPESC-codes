#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cvode/cvode.h>
#include <arkode/arkode_arkstep.h>
#include <sunlinsol/sunlinsol_spgmr.h>

#include "Advection-Diffusion.h"
#include "DiffOp.h"

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
   int nComp  = 1;  // number of components for each array
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
   Real advCoeffx = prob_data->advCoeffx;
   Real advCoeffy = prob_data->advCoeffy;

   // clear the RHS
   *rhs = 0.0;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // compute advection
   ComputeAdvectionUpwind(*sol, *rhs, *geom, 0, advCoeffx, advCoeffy);

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
   Real diffCoeff = prob_data->diffCoeff;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // clear the RHS
   *rhs = 0.0;

   // compute diffusion
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 0, diffCoeff);

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
   Real advCoeffx = prob_data->advCoeffx;
   Real advCoeffy = prob_data->advCoeffy;
   Real diffCoeff = prob_data->diffCoeff;

   // clear the RHS
   *rhs = 0.0;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // compute advection
   ComputeAdvectionUpwind(*sol, *rhs, *geom, 0, advCoeffx, advCoeffy);

   // compute diffusion
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 0, diffCoeff);

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
            Real x = prob_lo[0] + (((Real) i) + 0.5) * dx[0];

            if (x>0.25 && x<0.75 && y>0.25 && y<0.75)
            {
               fab(i,j,0,0) = sin(2 * M_PI * x - M_PI/2)
                  * sin(2 * M_PI * y - M_PI/2);
            }
            else
            {
               fab(i,j,0,0) = 0.0;
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

   // Default plot_int to -1, allow us to set it to something else in the inputs
   // file. If plot_int < 0 then no plot files will be written
   int plot_int = -1;
   pp.query("plot_int", plot_int);
   prob_opt.plot_int = plot_int;

   // Specify which integration method to use (defaults to CVODE)
   // 0 = CVODE
   // 1 = ARKStep
   int stepper = 0;
   pp.query("stepper", stepper);
   prob_opt.stepper = stepper;

   // Specify which RHS functions to use
   int rhs_adv  = 1;
   int rhs_diff = 1;
   pp.query("rhs_adv", rhs_adv);
   pp.query("rhs_diff", rhs_diff);
   prob_opt.rhs_adv  = rhs_adv;
   prob_opt.rhs_diff = rhs_diff;

   // Specify relative and absolute tolerances
   Real rtol = 1.0e-4;
   Real atol = 1.0e-9;
   pp.query("rtol", rtol);
   pp.query("atol", atol);
   prob_opt.rtol = rtol;
   prob_opt.atol = atol;

   // Specify final time for integration (default to 10,000)
   Real tfinal = 1.0e4;
   pp.query("tfinal", tfinal);
   prob_opt.tfinal = tfinal;

   // Specify output frequency (default to final time)
   Real dtout = tfinal;
   pp.query("dtout", dtout);
   prob_opt.dtout = dtout;

   // Advection coefficients
   Real advCoeffx = 5.0e-4;
   Real advCoeffy = 5.0e-4;
   pp.query("advCoeffx", advCoeffx);
   pp.query("advCoeffy", advCoeffy);
   prob_data.advCoeffx = advCoeffx;
   prob_data.advCoeffy = advCoeffy;

   // Diffusion coefficients
   Real diffCoeff = 2.0e-5;
   pp.query("diffCoeff", diffCoeff);
   prob_data.diffCoeff = diffCoeff;

   amrex::Print()
      << "n_cell        = " << n_cell        << std::endl
      << "max_grid_size = " << max_grid_size << std::endl
      << "plot_int      = " << plot_int      << std::endl
      << "stepper       = " << stepper       << std::endl
      << "rhs_adv       = " << rhs_adv       << std::endl
      << "rhs_diff      = " << rhs_diff      << std::endl
      << "rtol          = " << rtol          << std::endl
      << "atol          = " << atol          << std::endl
      << "tfinal        = " << tfinal        << std::endl
      << "dtout         = " << dtout         << std::endl;

   if (rhs_adv > 0)
      amrex::Print()
         << "advCoeffx     = " << advCoeffx << std::endl
         << "advCoeffy     = " << advCoeffy << std::endl;

   if (rhs_diff > 0)
      amrex::Print()
         << "diffCoeff     = " << diffCoeff << std::endl;

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
   RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                    {AMREX_D_DECL(1.0, 1.0, 1.0)});

   // This defines a Geometry object
   Vector<int> is_periodic(AMREX_SPACEDIM, 1);  // periodic in all direction
   geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());

   prob_data.geom = &geom;
}

void ComputeSolutionARK(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data)
{
   // Extract problem data and options
   Geometry* geom     = prob_data->geom;
   int       plot_int = prob_opt->plot_int;
   int       rhs_adv  = prob_opt->rhs_adv;
   int       rhs_diff = prob_opt->rhs_diff;
   Real      rtol     = prob_opt->rtol;
   Real      atol     = prob_opt->atol;
   Real      tfinal   = prob_opt->tfinal;
   Real      dtout    = prob_opt->dtout;

   // initial time, number of outputs, and error flag
   Real time = 0.0;
   int  nout = ceil(tfinal/dtout);
   int  ier  = 0;

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, time, 0);
   }

   // Create the ARK stepper
   void* arkode_mem = NULL;

   if (rhs_adv > 0 && rhs_diff > 0)
   {
      if (rhs_adv == 1 && rhs_diff == 1)
      {
         // explicit advection and diffusion
         arkode_mem = ARKStepCreate(ComputeRhsAdvDiff, NULL,
                                    time, nv_sol);
      }
      else if (rhs_adv == 1 && rhs_diff == 2)
      {
         // explicit advection and implicit diffusion
         arkode_mem = ARKStepCreate(ComputeRhsAdv, ComputeRhsDiff,
                                    time, nv_sol);
      }
      else if (rhs_adv == 2 && rhs_diff == 1)
      {
         // implicit advection and explicit diffusion
         arkode_mem = ARKStepCreate(ComputeRhsDiff, ComputeRhsAdv,
                                    time, nv_sol);
      }
      else if (rhs_adv == 2 && rhs_diff == 2)
      {
         // implicit advection and diffusion
         arkode_mem = ARKStepCreate(NULL, ComputeRhsAdvDiff,
                                    time, nv_sol);
      }
      else
      {
         amrex::Print() << "Invalid RHS options for ARKode" << std::endl;
         return;
      }
   }
   else if (rhs_adv > 0)
   {
      if (rhs_adv == 1)
      {
         // explicit advection
         arkode_mem = ARKStepCreate(ComputeRhsAdv, NULL,
                                    time, nv_sol);
      }
      else if (rhs_adv == 2)
      {
         // implicit advection
         arkode_mem = ARKStepCreate(NULL, ComputeRhsAdv,
                                    time, nv_sol);
      }
      else
      {
         amrex::Print() << "Invalid RHS options for ARKode" << std::endl;
         return;
      }
   }
   else if (rhs_diff > 0)
   {
      if (rhs_diff == 1)
      {
         // explicit diffusion
         arkode_mem = ARKStepCreate(ComputeRhsDiff, NULL,
                                    time, nv_sol);
      }
      else if (rhs_diff == 2)
      {
         // implicit diffusion
         arkode_mem = ARKStepCreate(NULL, ComputeRhsDiff,
                                    time, nv_sol);
      }
      else
      {
         amrex::Print() << "Invalid RHS options for ARKode" << std::endl;
         return;
      }
   }
   else
   {
      amrex::Print() << "Invalid RHS options for ARKode" << std::endl;
      return;
   }

   // Set ARKStep options
   ARKStepSStolerances(arkode_mem, atol, rtol);
   ARKStepSetUserData(arkode_mem, prob_data);

   // Create and attach GMRES linear solver (if necessary)
   if (rhs_adv == 2 || rhs_diff == 2)
   {
      SUNLinearSolver LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, 100);
      ier = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
      if (ier != ARKLS_SUCCESS)
      {
         amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
         return;
      }
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
         WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, tret, iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }
}

void ComputeSolutionCV(N_Vector nv_sol, ProblemOpt* prob_opt,
                       ProblemData* prob_data)
{
   // Extract problem data and options
   Geometry* geom     = prob_data->geom;
   int       plot_int = prob_opt->plot_int;
   int       rhs_adv  = prob_opt->rhs_adv;
   int       rhs_diff = prob_opt->rhs_diff;
   Real      rtol     = prob_opt->rtol;
   Real      atol     = prob_opt->atol;
   Real      tfinal   = prob_opt->tfinal;
   Real      dtout    = prob_opt->dtout;

   // initial time, number of outputs, and error flag
   Real time = 0.0;
   int  nout = ceil(tfinal/dtout);
   int  ier  = 0;

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, time, 0);
   }

   // Create CVODE memory
   void* cvode_mem = CVodeCreate(CV_BDF);

   if (rhs_adv == 1 && rhs_diff == 1)
   {
      // implicit Advection and Diffusion
      CVodeInit(cvode_mem, ComputeRhsAdvDiff, time, nv_sol);
   }
   else if (rhs_adv == 1)
   {
      // implicit Advection
      CVodeInit(cvode_mem, ComputeRhsAdv, time, nv_sol);
   }
   else if (rhs_diff == 1)
   {
      // implicit Diffusion
      CVodeInit(cvode_mem, ComputeRhsDiff, time, nv_sol);
   }
   else
   {
      amrex::Print() << "Invalid RHS options for CVODE" << std::endl;
      return;
   }

   // Set CVODE options
   CVodeSStolerances(cvode_mem, atol, rtol);
   CVodeSetUserData(cvode_mem, prob_data);

   // Create and attach GMRES linear solver
   SUNLinearSolver LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, 100);
   ier = CVodeSetLinearSolver(cvode_mem, LS, NULL);
   if (ier != CVLS_SUCCESS)
   {
      amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
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
         WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, tret, iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }
}
