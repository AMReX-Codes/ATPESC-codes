DEBUG = FALSE

USE_MPI  = TRUE
USE_OMP  = FALSE

USE_HYPRE = FALSE
USE_PETSC = FALSE

HYPRE_DIR  ?= .
PETSC_HOME ?= $(PETSC_DIR)/$(PETSC_ARCH)

COMP = gnu

DIM = 2 

AMREX_HOME ?= /home/adener/amrex

INCLUDE_LOCATIONS += $(PETSC_DIR)/$(PETSC_ARCH)/include
INCLUDE_LOCATIONS += $(PETSC_DIR)/include
LIBRARY_LOCATIONS += $(PETSC_DIR)/$(PETSC_ARCH)/lib
LIBRARIES += -lpetsc

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs 	:= Base Boundary AmrCore LinearSolvers/C_CellMG LinearSolvers/MLMG

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

