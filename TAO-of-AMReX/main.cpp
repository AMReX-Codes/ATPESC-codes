
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

    /* ~~~~~ PETSc/TAO code begins here ~~~~~ */
    PetscErrorCode     ierr;
    PetscReal          zero=0.0;
    int                nb, nl, nt;
    int                Nb, Nl, Nt;
    Vec                *Plist;
    Vec                P;
    Tao                tao;
    PetscBool          flg;
    PetscMPIInt        size,rank;

    ierr = PetscInitialize(&argc, &argv, (char*)0, (char*)0); if (ierr) return ierr;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);

    mytest.get_number_local_bcs(nb, nl, nt);
    mytest.get_number_global_bcs(Nb, Nl, Nt);

    // For debugging ...
    amrex::Print() << "num local: " << nb << " " << nl << " " << nt << "\n";
    amrex::Print() << "num global: " << Nb << " " << Nl << " " << Nt << "\n";

    ierr = VecCreateMPI(PETSC_COMM_WORLD, nb, Nb, &Plist[0]);CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, nl, Nl, &Plist[1]);CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, nt, Nt, &Plist[2]);CHKERRQ(ierr);
    ierr = VecCreateNest(PETSC_COMM_WORLD, 3, NULL, Plist, &P);CHKERRQ(ierr);
    ierr = VecSet(P, zero); CHKERRQ(ierr);

    mytest.set_target_solution(TargetSolution);

    ierr = TaoCreate(PETSC_COMM_WORLD, &tao); CHKERRQ(ierr);
    ierr = TaoSetType(tao, TAOBQNLS); CHKERRQ(ierr);
    ierr = TaoSetInitialVector(tao, P); CHKERRQ(ierr);
    ierr = TaoSetObjectiveAndGradientRoutine(tao, FormFunctionGradient, &mytest); CHKERRQ(ierr);
    ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);
    ierr = TaoSolve(tao); CHKERRQ(ierr);

    ierr = PetscFinalize();

    amrex::Finalize();
}

/* -------------------------------------------------------------------- */

amrex::Real TargetSolution(amrex::Real* coords)
{
    amrex::Real x = coords[0];
    amrex::Real y = coords[1];
    amrex::Real z = coords[2];
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
    PetscErrorCode    ierr;
    PetscReal         ff=0;
    Vec               *Plist, *Glist;
    PetscInt          nb, nl, nt;
    const PetscScalar *pb, *pl, *pt;
    PetscScalar       *gb, *gl, *gt;

    PetscFunctionBeginUser;

    PetscInt          *numNested;
    *numNested = 3;

    /* 1) get boundary values from array x and set them into AMReX */
    ierr = VecNestGetSubVecs(P, numNested, &Plist);CHKERRQ(ierr);

    ierr = VecGetArrayRead(Plist[0], &pb);CHKERRQ(ierr);
    ierr = VecGetArrayRead(Plist[1], &pl);CHKERRQ(ierr);
    ierr = VecGetArrayRead(Plist[2], &pt);CHKERRQ(ierr);

    ierr = VecGetLocalSize(Plist[0], &nb);CHKERRQ(ierr);
    ierr = VecGetLocalSize(Plist[1], &nl);CHKERRQ(ierr);
    ierr = VecGetLocalSize(Plist[2], &nt);CHKERRQ(ierr);

    mytest->update_boundary_values((int)nb, (const amrex::Real*)pb,
                                   (int)nl, (const amrex::Real*)pl,
                                   (int)nt, (const amrex::Real*)pt);

    ierr = VecRestoreArrayRead(Plist[0], &pb);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(Plist[1], &pl);CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(Plist[2], &pt);CHKERRQ(ierr);

    /* 2) solve the Poisson problem with AMReX */
    mytest->solve(); // solve the poisson problem with the boundary values

    /* 3) compute the objective function value using solution from AMReX and store into ff */
    ff = mytest->calculate_obj_val();
    *f = ff;

    /* 4) solve the linearized adjoint problem with AMReX (in this case same as forward linearized problem) */
    mytest->setup_adjoint_system(); // calculate the adjoint system rhs
    mytest->solve_adjoint_system(); // solve for the adjoint

    /* 5) compute the gradient per the Overleaf document and store it into array g */
    ierr = VecNestGetSubVecs(G, numNested, &Glist);CHKERRQ(ierr);

    ierr = VecGetArray(Glist[0], &gb);CHKERRQ(ierr);
    ierr = VecGetArray(Glist[1], &gl);CHKERRQ(ierr);
    ierr = VecGetArray(Glist[2], &gt);CHKERRQ(ierr);

    mytest->calculate_opt_gradient((int)nb, (amrex::Real*)gb,
                                  (int)nl, (amrex::Real*)gl,
                                  (int)nt, (amrex::Real*)gt);

    ierr = VecRestoreArray(Glist[0], &gb);CHKERRQ(ierr);
    ierr = VecRestoreArray(Glist[1], &gl);CHKERRQ(ierr);
    ierr = VecRestoreArray(Glist[2], &gt);CHKERRQ(ierr);

    /* 6) Other misc operations like generating iterative plots can be done here */
    mytest->write_plotfile();
    mytest->update_counter();

    PetscFunctionReturn(0);
}
