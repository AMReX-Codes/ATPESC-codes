
#include <AMReX.H>
#include "MyTest.H"

#include <petsc.h>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    BL_PROFILE("main");
    MyTest mytest;
    // mytest.dfdp is a vector of Reals containing the derivative of f wrt BC values p:
    // {dfdp_x_lo, dfdp_x_hi, dfdp_y_lo, dfdp_y_hi, dfdp_z_lo, dfdp_z_hi}
    mytest.dfdp.resize(2*AMREX_SPACEDIM);

    /* ~~~~~ PETSc/TAO code begins here ~~~~~ */
    PetscErrorCode     ierr;
    PetscReal          zero=0.0;
    Vec                x;
    Tao                tao;
    PetscBool          flg;
    PetscMPIInt        size,rank;

    ierr = PetscInitialize(&argc, &argv, (char*)0, help); if (ierr) return ierr;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
    if (size >1) SETERRQ(PETSC_COMM_SELF, 1, "Incorrect number of processors");

    ierr = VecCreateSeq(PETSC_COMM_SELF, mytest.boundary.max_size(), &x);CHKERRQ(ierr); // vec sizing must come from AMReX
    ierr = TaoCreate(PETSC_COMM_SELF, &tao); CHKERRQ(ierr);
    ierr = TaoSetType(tao, TAOBQNLS); CHKERRQ(ierr);
    ierr = VecSet(x, zero); CHKERRQ(ierr); // starting point does not have to be zero
    ierr = TaoSetInitialVector(tao, x); CHKERRQ(ierr);
    ierr = TaoSetObjectiveAndGradientRoutine(tao, FormFunctionGradient, &mytest); CHKERRQ(ierr);
    ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);

    ierr = TaoSolve(tao); CHKERRQ(ierr);

    ierr = PetscFinalize();

    amrex::Finalize();
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
PetscErrorCode FormFunctionGradient(Tao tao, Vec X, PetscReal *f, Vec G, void *ptr)
{
    MyTest            *mytest = (MyTest *) ptr;
    PetscInt          i;
    PetscErrorCode    ierr;
    PetscReal         ff=0;
    PetscScalar       *gg;
    const PetscScalar *xx;

    PetscFunctionBeginUser;
    /* Get pointers to PETSc vector data */
    ierr = VecGetArrayRead(X, &xx); CHKERRQ(ierr);
    ierr = VecGetArray(G, &gg); CHKERRQ(ierr);

    /* ~~~~~~ AMReX code begins here ~~~~~~ */

    /* 1) get boundary values from array x and set them into AMReX */
    // get boundary values to set from TAO (these are available in the array xx)
    // TODO ...
    mytest.update_boundary_values(); // set the boundary values

    /* 2) solve the nonlinar problem with AMReX */
    mytest.solve(); // solve the poisson problem with the boundary values

    /* 3) compute the objective function value using solution from AMReX and store into ff */
    // TODO ...

    /* 4) solve the linearized adjoint problem with AMReX (in this case same as forward linearized problem) */
    mytest.setup_adjoint_system(); // calculate the adjoint system rhs
    mytest.solve_adjoint_system(); // solve for the adjoint

    /* 5) compute the gradient per the Overleaf document and store it into array g */
    mytest.calculate_opt_gradient(); // calculate the adjoint contribution to df/dp
    // pass the gradient to TAO (must be written into array gg)
    // TODO ...

    /* 6) Other misc operations like generating iterative plots can be done here */
    mytest.writePlotfile();

    /* ~~~~~~ AMReX code ends here ~~~~~~ */

    /* Restore PETSc vectors */
    ierr = VecRestoreArrayRead(X, &xx); CHKERRQ(ierr);
    ierr = VecRestoreArray(G, &gg); CHKERRQ(ierr);
    *f   = ff; // return the objective function value
    PetscFunctionReturn(0);
}
