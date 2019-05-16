#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <arkode/arkode_mristep.h>
#include <arkode/arkode_arkstep.h>
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
   // AMREX_SPACEDIM: number of dimensions
   int n_cell, max_grid_size, nsteps, plot_int;
   GrayScottProblem problem;
   ParseInputs(n_cell, max_grid_size, nsteps, plot_int, problem);

   // make BoxArray and Geometry
   BoxArray ba;
   Geometry geom;
   SetUpGeometry(ba, geom, problem, n_cell, max_grid_size);


   // --- Create solution vectors ---
   // How Boxes are distrubuted among MPI processes
   DistributionMapping dm(ba);

   // we allocate two soln multifabs; one will store the old state, the other
   // the new
   int nGhost = 1;  // number of ghost cells for each array
   int nComp  = 2;  // number of components for each array
   MultiFab sol_old(ba, dm, nComp, nGhost);
   MultiFab sol_new(ba, dm, nComp, nGhost);

   // build the flux multifabs
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

   sunindextype length = nComp * n_cell * n_cell;
   N_Vector nv_sol_old = N_VMake_Multifab(length, &sol_old);
   N_Vector nv_sol_new = N_VMake_Multifab(length, &sol_new);

#if (AMREX_SPACEDIM == 3)
   FillInitConds3D(sol_new, geom);
#else
   FillInitConds2D(sol_new, geom);
#endif


   // --- Time advance to end ---
   //ComputeSolution(sol_old, sol_new, problem, nsteps, plot_int);
   //ComputeSolutionNV(nv_sol_old, nv_sol_new, &problem, nsteps, plot_int);
   //ComputeSolutionMRI(nv_sol_old, nv_sol_new, &problem, nsteps, plot_int);
   ComputeSolutionAO(nv_sol_old, nv_sol_new, &problem, nsteps, plot_int);


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

int ComputeDiffusionNV(realtype t, N_Vector sol, N_Vector diffusion,
                       void* problem)
{
   MultiFab* sol_mf = NV_MFAB(sol);
   MultiFab* diffusion_mf = NV_MFAB(diffusion);
   GrayScottProblem *gs_problem = (GrayScottProblem*) problem;

   ComputeDiffusion(*sol_mf, *diffusion_mf, *gs_problem);

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

void ComputeReactions3D(MultiFab& sol, MultiFab& reactions,
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

      for (int k = lo.z; k <= hi.z; ++k) {
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               Real temp = sol_fab(i,j,0,0) * sol_fab(i,j,0,1) * sol_fab(i,j,0,1);
               reactions_fab(i,j,k,0) = A * (1.0 - sol_fab(i,j,k,0)) - temp;
               reactions_fab(i,j,k,1) = temp - (A + B) * sol_fab(i,j,k,1);
            }
         }
      }
   }
}

int ComputeReactionsNV(realtype t, N_Vector sol, N_Vector reactions,
                       void* problem)
{
   MultiFab* sol_mf = NV_MFAB(sol);
   MultiFab* reactions_mf = NV_MFAB(reactions);
   GrayScottProblem *gs_problem = (GrayScottProblem*) problem;

#if (AMREX_SPACEDIM == 3)
   ComputeReactions3D(*sol_mf, *reactions_mf, *gs_problem);
#else
   ComputeReactions2D(*sol_mf, *reactions_mf, *gs_problem);
#endif

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

void FillInitConds3D(MultiFab& sol, const Geometry& geom)
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
      for (int k = lo.z; k <= hi.z; ++k) {
         Real z = prob_lo[2] + (((Real) k) + 0.5) * dx[2];

         for (int j = lo.y; j <= hi.y; ++j) {
            Real y = prob_lo[1] + (((Real) j) + 0.5) * dx[1];

            for (int i = lo.x; i <= hi.x; ++i) {
               fab(i,j,k,0) = 0.;
               fab(i,j,k,1) = 0.;

               Real x = prob_lo[0] + (((Real) i) + 0.5) * dx[0];
               if (x>=1.14 && x<=1.33 && y>=1.14 && y<=1.33 && z>=1.14 && z<=1.33)
               {
                  fab(i,j,k,0) = 0.5;
                  fab(i,j,k,1) = 0.25;
               }
            }
         }
      }
   }
}

void ParseInputs(int& n_cell, int& max_grid_size,
                 int& nsteps, int& plot_int,
                 GrayScottProblem& problem)
{
   Real diffCoeffU, diffCoeffV, A, B;

   // ParmParse is way of reading inputs from the inputs file
   ParmParse pp;

   // We need to get n_cell from the inputs file - this is the number of cells
   // on each side of
   //   a square (or cubic) domain.
   pp.get("n_cell", n_cell);

   // The domain is broken into boxes of size max_grid_size
   pp.get("max_grid_size", max_grid_size);

   // Default plot_int to -1, allow us to set it to something else in the inputs
   // file. If plot_int < 0 then no plot files will be written
   plot_int = -1;
   pp.query("plot_int", plot_int);

   // Default nsteps to 10, allow us to set it to something else in the inputs
   // file
   nsteps = 10;
   pp.query("nsteps", nsteps);

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

void ComputeSolution(MultiFab& sol_old, MultiFab& sol_new,
                     GrayScottProblem& problem,
                     int nsteps, int plot_int)
{
   Geometry* geom = problem.geom;

   // time = starting time in the simulation
   Real time = 0.0;

   // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined
   // in the inputs file)
   if (plot_int > 0)
   {
      int n = 0;
      const std::string& pltfile = amrex::Concatenate("plt", n, 5);
      WriteSingleLevelPlotfile(pltfile, sol_new, {"u", "v"}, *geom, time, 0);
   }

   // allocate the RHS multifabs
   const BoxArray& ba = sol_old.boxArray();
   const DistributionMapping& dm = sol_old.DistributionMap();
   int nComp = sol_old.nComp();
   int nGhost = sol_old.nGrow();
   MultiFab RHS(ba, dm, nComp, 0);
   MultiFab diffusion(ba, dm, nComp, 0);

   // compute the time steps
   Real dt = 1.0;
   for (int n = 1; n <= nsteps; ++n)
   {
      MultiFab::Copy(sol_old, sol_new, 0, 0, nComp, 0);

      ComputeDiffusion(sol_old, diffusion, problem);
      ComputeReactions2D(sol_old, RHS, problem);

      // Euler time step
      MultiFab::Add(RHS, diffusion, 0, 0, nComp, 0);
      RHS.mult(dt, 0, nComp, 0);
      MultiFab::Add(sol_new, RHS, 0, 0, nComp, 0);
      time = time + dt;

      // Tell the I/O Processor to write out which step we're doing
      amrex::Print() << "Advanced step " << n << "\n";

      // Write a plotfile of the current data (plot_int was defined in the
      // inputs file)
      if (plot_int > 0 && n%plot_int == 0)
      {
         const std::string& pltfile = amrex::Concatenate("plt", n, 5);
         WriteSingleLevelPlotfile(pltfile, sol_new, {"u", "v"}, *geom, time, n);
      }
   }
}

void ComputeSolutionNV(N_Vector sol_old,
                       N_Vector sol_new,
                       GrayScottProblem* problem,
                       int nsteps, int plot_int)
{
   Geometry* geom = problem->geom;

   // time = starting time in the simulation
   Real time = 0.0;

   // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined
   // in the inputs file)
   if (plot_int > 0)
   {
      int n = 0;
      const std::string& pltfile = amrex::Concatenate("plt", n, 5);
      MultiFab* sol_new_mf = NV_MFAB(sol_new);
      WriteSingleLevelPlotfile(pltfile, *sol_new_mf, {"u", "v"},
                               *geom, time, 0);
   }

   // allocate the RHS multifabs
   N_Vector RHS = N_VClone_Multifab(sol_old);
   N_Vector diffusion = N_VClone_Multifab(sol_old);

   // compute the time steps
   Real dt = 1.0;
   for (int n = 1; n <= nsteps; ++n)
   {
      N_VScale(1.0, sol_new, sol_old);
      //MultiFab::Copy(sol_old, sol_new, 0, 0, nComp, 0);

      ComputeDiffusionNV(time, sol_old, diffusion, problem);
      ComputeReactionsNV(time, sol_old, RHS, problem);

      // Euler time step
      N_VLinearSum(1.0, RHS, 1.0, diffusion, RHS);
      N_VLinearSum(1.0, sol_new, dt, RHS, sol_new);
      time = time + dt;

      // Tell the I/O Processor to write out which step we're doing
      amrex::Print() << "Advanced step " << n << "\n";

      // Write a plotfile of the current data (plot_int was defined in the
      // inputs file)
      if (plot_int > 0 && n%plot_int == 0)
      {
         const std::string& pltfile = amrex::Concatenate("plt", n, 5);
         MultiFab* sol_new_mf = NV_MFAB(sol_new);
         WriteSingleLevelPlotfile(pltfile, *sol_new_mf, {"u", "v"},
                                  *geom, time, n);
      }
   }
}

void ComputeSolutionMRI(N_Vector sol_old,
                        N_Vector sol_new,
                        GrayScottProblem* problem,
                        int nsteps, int plot_int)
{
   Geometry* geom = problem->geom;

   Real time   = 0.0;          // time = starting time in the simulation
   Real tFinal = nsteps * 1.0;
   int  ier    = 0;            // error flag

   // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined
   // in the inputs file)
   if (plot_int > 0)
   {
      int n = 0;
      const std::string& pltfile = amrex::Concatenate("plt", n, 5);
      MultiFab* sol_new_mf = NV_MFAB(sol_new);
      WriteSingleLevelPlotfile(pltfile, *sol_new_mf, {"u", "v"},
                               *geom, time, 0);
   }

   // Copy new to old
   N_VScale(1.0, sol_new, sol_old);

   realtype tret;
   void* arkode_mem = MRIStepCreate(ComputeDiffusionNV, ComputeReactionsNV,
                                    time, sol_old);
   MRIStepSetFixedStep(arkode_mem, 0.5, 0.5);
   MRIStepSetMaxNumSteps(arkode_mem, 500000);
   MRIStepSetUserData(arkode_mem, problem);
   ier = MRIStepEvolve(arkode_mem, 1000.0, sol_new, &tret, ARK_NORMAL);
   if (ier < 0)
   {
      amrex::Print() << "Error in MRIStepEvolve" << std::endl;
      return;
   }

   long nfs_evals, nff_evals;
   MRIStepGetNumRhsEvals(arkode_mem, &nfs_evals, &nff_evals);
   amrex::Print() << nfs_evals << " slow evals "
                  << nff_evals << " fast evals" << std::endl;

   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", nsteps, 5);
      MultiFab* sol_new_mf = NV_MFAB(sol_new);
      WriteSingleLevelPlotfile(pltfile, *sol_new_mf, {"u", "v"},
                               *geom, nsteps, nsteps);
   }
}

void ComputeSolutionAO(N_Vector sol_old,
                       N_Vector sol_new,
                       GrayScottProblem* problem,
                       int nsteps, int plot_int)
{
   Geometry* geom = problem->geom;

   Real time   = 0.0;          // time = starting time in the simulation
   Real tFinal = nsteps * 1.0;
   int  ier    = 0;            // error flag

   // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined
   // in the inputs file)
   if (plot_int > 0)
   {
      int n = 0;
      const std::string& pltfile = amrex::Concatenate("plt", n, 5);
      MultiFab* sol_new_mf = NV_MFAB(sol_new);
      WriteSingleLevelPlotfile(pltfile, *sol_new_mf, {"u", "v"},
                               *geom, time, 0);
   }

   // Copy new to old
   N_VScale(1.0, sol_new, sol_old);

   realtype tret;
   //void* arkode_mem = ARKStepCreate(ComputeDiffusionNV, ComputeReactionsNV, time, sol_old);
   void* arkode_mem = ARKStepCreate(ComputeReactionsNV, ComputeDiffusionNV,
                                    time, sol_old);
   //ARKStepSStolerances(arkode_mem, 1.0, 1.0);
   ARKStepSetFixedStep(arkode_mem, 1.0);
   ARKStepSetMaxNumSteps(arkode_mem, 5000);
   ARKStepSetUserData(arkode_mem, problem);

   SUNLinearSolver LS = SUNLinSol_SPGMR(sol_new, PREC_NONE, 100);
   ier = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
   if (ier != ARKLS_SUCCESS)
   {
      amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
      return;
   }

   ier = ARKStepEvolve(arkode_mem, 1000.0, sol_new, &tret, ARK_NORMAL);
   if (ier < 0)
   {
      amrex::Print() << "Error in MRIStepEvolve" << std::endl;
      return;
   }
   //ARKStepEvolve(arkode_mem, 1e-6, sol_new, &tret, ARK_NORMAL);

   long nfe_evals, nfi_evals;
   ARKStepGetNumRhsEvals(arkode_mem, &nfe_evals, &nfi_evals);
   amrex::Print() << nfe_evals << " explicit evals "
                  << nfi_evals << " implicit evals" << std::endl;

   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", nsteps, 5);
      MultiFab* sol_new_mf = NV_MFAB(sol_new);
      WriteSingleLevelPlotfile(pltfile, *sol_new_mf, {"u", "v"},
                               *geom, nsteps, nsteps);
   }
}
