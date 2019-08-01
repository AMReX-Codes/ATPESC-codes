# ------------------------------------------------------------------------------
# Makefile to build examples
# ------------------------------------------------------------------------------
# Note the following environment variables need to be set to in order to build:
# AMREX_INSTALL_DIR = path to AMReX installation
# MPICXX = mpicxx wrapper
#
# If any are unset, they assume the default values for compilation on Cooley
# ------------------------------------------------------------------------------

DIM = 3

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

FC = gfortran

CPPFLAGS = -Ishared -I$(AMREX_INSTALL_DIR)/include -I. -I$(HYPRE_DIR)/include
CPPFLAGS = -I$(AMREX_INSTALL_DIR)/include -I. -I$(HYPRE_DIR)/include
CXXFLAGS = -O2 -std=c++14
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib -L$(HYPRE_DIR)/lib -Wl,-rpath,$(AMREX_INSTALL_DIR)/lib 

LIBRARIES = -lamrex -lHYPRE
LIBRARIES += -lgfortran

default: main$(DIM)d.exe

main$(DIM)d.exe: main.o MyParticleContainer.o writeEBsurface.o eb_to_pvf.o
	$(MPICXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

main.o: main.cpp MyParticleContainer.H writeEBsurface.H
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

MyParticleContainer.o: MyParticleContainer.cpp MyParticleContainer.H
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

writeEBsurface.o: writeEBsurface.cpp writeEBsurface.H
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

eb_to_pvf.o: eb_to_pvf.F90
	$(FC) -o $@ -c $(FFLAGS) $(CPPFLAGS) $<

.PHONY: clean realclean pltclean

clean:
	$(RM) *.o

realclean: clean
	$(RM) *~ *.exe

pltclean:
	$(RM) -rf plt*/ *.png eb_grid* eb*vtp

movie:
	ffmpeg -framerate 10 -pattern_type glob -i "p*.png" -r 5 -vf scale=1024:-1 pachinko.gif
