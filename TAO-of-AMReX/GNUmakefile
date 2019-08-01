DEBUG = FALSE

USE_MPI  = TRUE
USE_OMP  = FALSE

USE_HYPRE = TRUE
USE_PETSC = TRUE

HYPRE_DIR  ?= .
PETSC_HOME ?= .

COMP = gnu

DIM = 2 

AMREX_HOME ?= ../../..

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs 	:= Base Boundary AmrCore LinearSolvers/C_CellMG LinearSolvers/MLMG

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

