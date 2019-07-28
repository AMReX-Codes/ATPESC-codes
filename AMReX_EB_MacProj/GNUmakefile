DEBUG = FALSE
USE_MPI  = TRUE
USE_OMP  = FALSE

COMP = gnu

DIM = 2

AMREX_HOME ?= ../../..

USE_EB = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package

Pdirs := AmrCore 
Pdirs += Base 
Pdirs += Boundary 
Pdirs += EB
Pdirs += LinearSolvers/MLMG
Pdirs += Particle

Ppack	+= $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

movie:
	./plot.py
	ffmpeg -framerate 20 -pattern_type glob -i "*.png" -r 15 -vf scale=1024:-1 macproj.gif
