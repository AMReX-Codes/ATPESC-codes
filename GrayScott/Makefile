# ------------------------------------------------------------------------------
# Makefile to build 2D GrayScott (Reaction-Diffusion) example
# ------------------------------------------------------------------------------
# Note the following environment variables need to be set to in order to build:
# AMREX_INSTALL_DIR = path to AMReX installation
# SUNDIALS_INSTALL_DIR = path to SUNDIALS installation
# ------------------------------------------------------------------------------

CXX = mpicxx
FC = gfortran

CPPFLAGS = -I ../NVector_Multifab -I$(AMREX_INSTALL_DIR)/include -I$(SUNDIALS_INSTALL_DIR)/include
CXXFLAGS = -O2 -std=c++11
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib -L$(SUNDIALS_INSTALL_DIR)/lib64 -L../NVector_Multifab

LIBRARIES = -lamrex -lsundials_arkode -lnvector_multifab

LIBRARIES += -lgfortran

default: GrayScott.exe

GrayScott.exe: GrayScott.o DiffOp2D.o
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

GrayScott.o: GrayScott.cpp GrayScott.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

DiffOp2D.o: DiffOp2D.cpp DiffOp.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

.PHONY: movie clean realclean pltclean

movie:
	ls -1 plt*/Header | tee movie.visit

clean:
	$(RM) *.o

realclean: clean
	$(RM) *.exe

plotclean:
	$(RM) -rf plt*/
