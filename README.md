# ATPESC-codes

Example codes for ATPESC 2019

## Building

To build the examples using the provided Makefile set the environment variables
`AMREX_INSTALL_DIR` and `SUNDIALS_INSTALL_DIR` to the directories for the AMReX
and Sundials installations respectively. For example, if ARMeX was installed in
`~/apps/amrex` and and Sundials in `~/apps/sundials`, then the required
 environment variables can be created with
```
export AMREX_INSTALL_DIR=~/apps/amrex
export SUNDIALS_INSTALL_DIR=~/apps/sundials
```
for sh/bash/ksh shells or with
```
set AMREX_INSTALL_DIR ~/apps/amrex
set SUNDIALS_INSTALL_DIR ~/apps/sundials
```
for csh/tcsh shells. After the environment variables are set run `make` to build
the example executables, `Advection-Diffusion.exe` and `GrayScott.exe`.
