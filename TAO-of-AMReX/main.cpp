
#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        MyTest mytest;

        // mytest.dfdp is a vector of Reals containing the derivative of f wrt BC values p:
        // {dfdp_x_lo, dfdp_x_hi, dfdp_y_lo, dfdp_y_hi, dfdp_z_lo, dfdp_z_hi}
        mytest.dfdp.resize(2*AMREX_SPACEDIM);

        // Presumably this block is a loop body with stopping conditions from TAO?
        {
            // get boundary values to set from TAO
            // TODO ...

            // AMReX operations to solve for dfdp ...
            mytest.update_boundary_values(); // set the boundary values
            mytest.solve(); // solve the poisson problem with the boundary values
            mytest.setup_adjoint_system(); // calculate the adjoint system rhs
            mytest.solve_adjoint_system(); // solve for the adjoint
            mytest.calculate_opt_gradient(); // calculate the adjoint contribution to df/dp

            // we now have dfdp so we can pass it to TAO
            // TODO ...
        }

        mytest.writePlotfile();
    }

    amrex::Finalize();
}
