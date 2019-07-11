# ------------------------------------------------------------------------------
# Makefile to build examples
# ------------------------------------------------------------------------------
# Note the following environment variables need to be set to in order to build:
# AMREX_INSTALL_DIR = path to AMReX installation
# SUNDIALS_INSTALL_DIR = path to SUNDIALS installation
# MPICXX = mpicxx wrapper
#
# If any are unset, they assume the default values for compilation on Cooley
# ------------------------------------------------------------------------------

# set default values for any unset environment variables
ifeq ($(AMREX_INSTALL_DIR),)
  AMREX_INSTALL_DIR = /projects/ATPESC2019/FASTMath/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/amrex-develop-fvkknjb363aqxpt6zu2tt6oojixnlnex
endif
ifeq ($(SUNDIALS_INSTALL_DIR),)
  SUNDIALS_INSTALL_DIR = /projects/ATPESC2019/FASTMath/spack/opt/spack/linux-rhel7-x86_64/gcc-4.8.5/sundials-4.1.0-s7aljoyz5e53teghcs4miobyqszrphtt
endif
ifeq ($(MPICXX),)
  MPICXX = mpicxx
endif

CPPFLAGS = -Ishared -I$(AMREX_INSTALL_DIR)/include -I$(SUNDIALS_INSTALL_DIR)/include
CXXFLAGS = -O2 -std=c++11
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib -L$(SUNDIALS_INSTALL_DIR)/lib -L$(SUNDIALS_INSTALL_DIR)/lib64 -Wl,-rpath,$(AMREX_INSTALL_DIR)/lib,-rpath,$(SUNDIALS_INSTALL_DIR)/lib,-rpath,$(SUNDIALS_INSTALL_DIR)/lib64

LIBRARIES = -lamrex -lsundials_cvode -lsundials_arkode

LIBRARIES += -lgfortran

default: Advection-Diffusion.exe

Advection-Diffusion.exe: Advection-Diffusion.o shared/NVector_Multifab.o shared/Utilities.o
	$(MPICXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

Advection-Diffusion.o: Advection-Diffusion.cpp Advection-Diffusion.h
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

shared/NVector_Multifab.o: shared/NVector_Multifab.cpp shared/NVector_Multifab.h
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

shared/Utilities.o: shared/Utilities.cpp shared/Utilities.h
	$(MPICXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

.PHONY: movie clean realclean pltclean

movie:
	ls -1 plt*/Header | tee movie.visit

clean:
	$(RM) *.o shared/*.o

realclean: clean
	$(RM) *.exe

pltclean:
	$(RM) -rf plt*/
