
#include <AMReX.H>
#include "MyTest.H"

#include <petsc.h>

amrex::Real TargetSolution(amrex::Real* coords);
PetscErrorCode FormFunctionGradient(Tao tao, Vec P, PetscReal *f, Vec G, void *ptr);

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    BL_PROFILE("main");
    MyTest mytest;

    PetscErrorCode     ierr;
    PetscReal          zero=0.0;
    int                nb, nl, nt;
    int                Nb, Nl, Nt;
    Vec                Plist[3];
    Vec                P;
    Tao                tao;
    PetscBool          flg;
    PetscMPIInt        size,rank;

    ierr = PetscInitialize(&argc, &argv, (char*)0, (char*)0); if (ierr) return ierr;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);

    // Create PETSc vectors for the bottom, left and top edge Dirichlet boundaries
    mytest.get_number_local_bcs(nb, nl, nt);
    mytest.get_number_global_bcs(Nb, Nl, Nt);

    // For debugging ...
    amrex::Print() << "num local: " << nb << " " << nl << " " << nt << "\n";
    amrex::Print() << "num global: " << Nb << " " << Nl << " " << Nt << "\n";

    ierr = VecCreateMPI(PETSC_COMM_WORLD, nb, Nb, &Plist[0]);CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, nl, Nl, &Plist[1]);CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, nt, Nt, &Plist[2]);CHKERRQ(ierr);

    // Combine the vectors for each edge into a single vector, initialize with zero
    ierr = VecCreateNest(PETSC_COMM_WORLD, 3, NULL, Plist, &P);CHKERRQ(ierr);
    ierr = VecSet(P, zero); CHKERRQ(ierr);

    // Set the target solution we want to recover on the right edge of the domain
    // u_target = 10.0 - y^2
    mytest.set_target_solution(TargetSolution);

    // Create and setup the TAO optimization algorithm
    ierr = TaoCreate(PETSC_COMM_WORLD, &tao); CHKERRQ(ierr);
    ierr = TaoSetType(tao, TAOBQNLS); CHKERRQ(ierr); // TAOBQNLS is a bound-constrained quasi-Newton linesearch alg,
    ierr = TaoSetInitialVector(tao, P); CHKERRQ(ierr);
    ierr = TaoSetObjectiveAndGradientRoutine(tao, FormFunctionGradient, &mytest); CHKERRQ(ierr);
    ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);

    // Start the optimization solution
    ierr = TaoSolve(tao); CHKERRQ(ierr);

    // Cleanup and exit
    ierr = VecDestroy(&P);CHKERRQ(ierr);
    ierr = PetscFinalize();
    amrex::Finalize();
}

/* -------------------------------------------------------------------- */

amrex::Real TargetSolution(amrex::Real* coords)
{
    amrex::Real y = coords[1];
    amrex::Real utarg = 10.0 - (y*y);
    return utarg;
}

/* -------------------------------------------------------------------- */
/*
    FormFunctionGradient - Evaluates the function, f(X), and gradient, G(X).

    Input Parameters:
.   tao  - the Tao context
.   X    - input vector
.   ptr  - optional user-defined context, as set by TaoSetFunctionGradient()

    Output Parameters:
.   G - vector containing the newly evaluated gradient
.   f - function value

    Note:
    Some optimization methods ask for the function and the gradient evaluation
    at the same time.  Evaluating both at once may be more efficient that
    evaluating each separately.
*/
PetscErrorCode FormFunctionGradient(Tao tao, Vec P, PetscReal *f, Vec G, void *ptr)
{
    MyTest            *mytest = (MyTest *) ptr;
    PetscInt          i;
    PetscInt          numNested;
    PetscErrorCode    ierr;
    PetscReal         ff=0;
    Vec               *Plist, *Glist;
    PetscInt          nb, nl, nt;
    const PetscScalar *pb, *pl, *pt;
    PetscScalar       *gb, *gl, *gt;

    PetscFunctionBeginUser;

    // Break up the vector of optimization variables into subvectors for bottom, left and top Dirichlet boundaries
    ierr = VecNestGetSubVecs(P, &numNested, &Plist);CHKERRQ(ierr);

    // Extract the subvector data as C arrays
    // NOTE: this is read-only to prevent accidental changes of the optimization variables
    ierr = VecGetArrayRead(Plist[0], &pb);CHKERRQ(ierr);
    ierr = VecGetArrayRead(Plist[1], &pl);CHKERRQ(ierr);
    ierr = VecGetArrayRead(Plist[2], &pt);CHKERRQ(ierr);

    ierr = VecGetLocalSize(Plist[0], &nb);CHKERRQ(ierr);
    ierr = VecGetLocalSize(Plist[1], &nl);CHKERRQ(ierr);
    ierr = VecGetLocalSize(Plist[2], &nt);CHKERRQ(ierr);

    // Set the boundary values from TAO into the AMReX solver
    mytest->update_boundary_values((int)nb, (const amrex::Real*)pb,
                                   (int)nl, (const amrex::Real*)pl,
                                   (int)nt, (const amrex::Real*)pt);

    // Extracted C arrays must be restored back into the parent PETSc vectors
    ierr = VecRestoreArrayRead(Plist[0], &pb);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(Plist[1], &pl);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(Plist[2], &pt);CHKERRQ(ierr);

    /// Solve the Laplace equations with the prescribed boundary values
    mytest->solve();

    // Compute the objective function using the AMReX solution
    // f = 0.5 * \int_0^1 (u(1, y) - u_targ)^2 dy
    // u_targ = 10.0 - y^2
    ff = mytest->calculate_obj_val();
    *f = ff;

    // Solve the adjoint problem
    // NOTE: Laplace equation is self-adjoint so we re-do the forward solution with a new RHS vector
    mytest->setup_adjoint_system(); 
    mytest->solve_adjoint_system(); 

    // The gradient copies the structure of the optimization variable with subvectors corresponding to 
    // the component of the gradient on each Dirichlet boundary edge
    ierr = VecNestGetSubVecs(G, &numNested, &Glist);CHKERRQ(ierr);

    // Extract the gradient data as C arrays
    // NOTE: this is read/write enabled so that we can copy gradient info from AMReX into the PETSc vector
    ierr = VecGetArray(Glist[0], &gb);CHKERRQ(ierr);
    ierr = VecGetArray(Glist[1], &gl);CHKERRQ(ierr);
    ierr = VecGetArray(Glist[2], &gt);CHKERRQ(ierr);

    // Calculate the gradient values using the adjoint computed above
    mytest->calculate_opt_gradient((int)nb, (amrex::Real*)gb,
                                  (int)nl, (amrex::Real*)gl,
                                  (int)nt, (amrex::Real*)gt);

    // Restore the arrays back into their parent vectors
    ierr = VecRestoreArray(Glist[0], &gb);CHKERRQ(ierr);
    ierr = VecRestoreArray(Glist[1], &gl);CHKERRQ(ierr);
    ierr = VecRestoreArray(Glist[2], &gt);CHKERRQ(ierr);

    // Perform other misc operations like visualization and I/O
    mytest->write_plotfile();
    mytest->update_counter();

    PetscFunctionReturn(0);
}
