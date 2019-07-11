# ATPESC-codes

SUNDIALS+AMReX example codes for ATPESC 2019

## Building

To build the examples using the provided Makefile set the environment variables
`AMREX_INSTALL_DIR` and `SUNDIALS_INSTALL_DIR` to the directories for the AMReX
and Sundials installations respectively.  Also set the environment
variable `MPICXX` to the MPI C++ wrapper to use for compilation.  For
example, if AMReX was installed in `~/apps/amrex` and and Sundials in
`~/apps/sundials`, and the MPI C++ wrapper is just `mpicxx`, then the required 
 environment variables can be created with
```
export AMREX_INSTALL_DIR=~/apps/amrex
export SUNDIALS_INSTALL_DIR=~/apps/sundials
export MPICXX=mpicxx
```
for sh/bash/ksh/zsh shells or with
```
set AMREX_INSTALL_DIR ~/apps/amrex
set SUNDIALS_INSTALL_DIR ~/apps/sundials
set MPICXX mpicxx
```
for csh/tcsh shells. After the environment variables are set run `make` to build
the example executables, `Advection-Diffusion.exe` and `GrayScott.exe`.

If any of these environment variables are left unset, they default to
valid values for compilation on Cooley for the ATPESC 2019 workshop.
