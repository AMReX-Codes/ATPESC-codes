AMREX_LIBRARY_HOME ?= /home/dewillcox/dev-astro/amrex/libamrex

LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include

CFLAGS := -I$(INCDIR) -Werror=return-type -g -O3 -std=c++14
LFLAGS := -L$(LIBDIR) -lamrex -L/usr/lib/gcc/x86_64-linux-gnu/5/ -Wl,-Bsymbolic-functions -Wl,-z,relro -I/usr/include/mpich -I/usr/include/mpich -L/usr/lib/x86_64-linux-gnu -lmpichfort -lmpich -lgfortran -lquadmath

all:
	g++ -o main.exe main.cpp $(CFLAGS) $(LFLAGS)
