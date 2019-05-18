# ------------------------------------------------------------------------------
# Makefile to create static library for MultiFAB NVector wrapper
# ------------------------------------------------------------------------------
# Note the following environment variables need to be set to in order to build:
# AMREX_INSTALL_DIR = path to AMReX installation
# SUNDIALS_INSTALL_DIR = path to SUNDIALS installation
# ------------------------------------------------------------------------------

CXX = mpicxx
FC = gfortran

CPPFLAGS = -I$(AMREX_INSTALL_DIR)/include -I$(SUNDIALS_INSTALL_DIR)/include
CXXFLAGS = -O2 -std=c++11
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib

LIBRARIES = -lamrex

LIBRARIES += -lgfortran

default: libs

NVector_Multifab.o: NVector_Multifab.cpp NVector_Multifab.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

libnvector_multifab.a: NVector_Multifab.o
	ar rsc libnvector_multifab.a NVector_Multifab.o

libs: libnvector_multifab.a

.PHONY: clean realclean

clean:
	$(RM) *.o

realclean: clean
	$(RM) *.a
