# ------------------------------------------------------------------------------
# Makefile to build examples
# ------------------------------------------------------------------------------
# Note the following environment variables need to be set to in order to build:
# AMREX_INSTALL_DIR = path to AMReX installation
# SUNDIALS_INSTALL_DIR = path to SUNDIALS installation
# ------------------------------------------------------------------------------

ifeq ($(CXX),)
  CXX = mpicxx
endif
ifeq ($(FC),)
  FC = gfortran
endif

CPPFLAGS = -Ishared -I$(AMREX_INSTALL_DIR)/include -I$(SUNDIALS_INSTALL_DIR)/include
CXXFLAGS = -O2 -std=c++11
FFLAGS = -O2
LDFLAGS = -L$(AMREX_INSTALL_DIR)/lib -L$(SUNDIALS_INSTALL_DIR)/lib

LIBRARIES = -lamrex -lsundials_cvode -lsundials_arkode

LIBRARIES += -lgfortran

default: Advection-Diffusion.exe GrayScott.exe

Advection-Diffusion.exe: Advection-Diffusion/Advection-Diffusion.o shared/NVector_Multifab.o shared/DiffOp2D.o
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

Advection-Diffusion/Advection-Diffusion.o: Advection-Diffusion/Advection-Diffusion.cpp Advection-Diffusion/Advection-Diffusion.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

GrayScott.exe: GrayScott/GrayScott.o shared/NVector_Multifab.o shared/DiffOp2D.o shared/Reactions.o
	$(CXX) -o $@ $(CXXFLAGS) $^ $(LDFLAGS) $(LIBRARIES)

GrayScott/GrayScott.o: GrayScott/GrayScott.cpp GrayScott/GrayScott.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

shared/NVector_Multifab.o: shared/NVector_Multifab.cpp shared/NVector_Multifab.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

shared/DiffOp2D.o: shared/DiffOp2D.cpp shared/DiffOp.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

shared/Reactions.o: shared/Reactions.cpp shared/Reactions.h
	$(CXX) -o $@ -c $(CXXFLAGS) $(CPPFLAGS) $<

.PHONY: movie clean realclean pltclean

movie:
	ls -1 plt*/Header | tee movie.visit

clean:
	$(RM) Advection-Diffusion/*.o GrayScott/*.o shared/*.o

realclean: clean
	$(RM) *.exe

pltclean:
	$(RM) -rf plt*/
