
#include <AMReX.H>
#include "MyTest.H"

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    {
        BL_PROFILE("main");
        MyTest mytest;
        mytest.dfdp.resize(2*AMREX_SPACEDIM);

        {
            // get boundary values to set
            mytest.update_boundary_values(); // set the boundary values
            mytest.solve(); // solve the poisson problem with the boundary values
            mytest.setup_adjoint_system(); // calculate the adjoint system rhs
            mytest.solve_adjoint_system(); // solve for the adjoint
            mytest.calculate_opt_gradient(); // calculate the adjoint contribution to df/dp
        }
        
        mytest.writePlotfile();
    }

    amrex::Finalize();
}
