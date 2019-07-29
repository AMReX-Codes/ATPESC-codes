# ------------------------------------------------------------------------------
# Makefile to build examples
# ------------------------------------------------------------------------------
# Note the following environment variables need to be set to in order to build:
# AMREX_INSTALL_DIR = path to AMReX installation
# MPICXX = mpicxx wrapper
#
# If any are unset, they assume the default values for compilation on Cooley
# ------------------------------------------------------------------------------

DIM = 2

# set default values for any unset environment variables
ifeq ($(AMREX_INSTALL_DIR),)
  AMREX_INSTALL_DIR = /projects/ATPESC2019/FASTMath/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/amrex-develop-fvkknjb363aqxpt6zu2tt6oojixnlnex
endif

# These will only apply if they are not set as environment variables
AMREX_INSTALL_DIR ?= /PATH/TO/AMREXLIB
HYPRE_DIR         ?= /PATH/TO/HYPRE

ifeq ($(MPICXX),)
  MPICXX = mpicxx
endif

CPPFLAGS = -Ishared -I$(AMREX_INSTALL_DIR)/include -I. -I$(HYPRE_DIR)/include
CXXFLAGS = -O2 -std=c++14
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib -L$(HYPRE_DIR)/lib -Wl,-rpath,$(AMREX_INSTALL_DIR)/lib 

LIBRARIES = -lamrex -lHYPRE

LIBRARIES += -lgfortran

main$(DIM)d.exe: main.cpp MyParticleContainer.cpp
	$(MPICXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

.PHONY: movie clean realclean pltclean

movie:
	ls -1 plt*/Header | tee movie.visit

clean:
	$(RM) *.o

realclean: clean
	$(RM) *~ *.exe

pltclean:
	$(RM) -rf plt*/ *.png *_diagnostics.txt __pycache__
