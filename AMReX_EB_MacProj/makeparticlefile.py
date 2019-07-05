#This will create a tracer particle input file based on the number
#of particles in each direction (NX,NY,NZ) and the physical domain
#boundaries (dimx,dimy,dimz). Amrex can then read this file in with
#the function TracerParticleContainer::InitFromAsciiFile()

import numpy as np

NX = 4
NY = 12
NZ = 12

dimx = [0,0.1]  #They'll start evenly close to the edge/source and flow through
dimy = [0,1]
dimz = [0,1]

dx = (dimx[1] - dimx[0]) / NX
dy = (dimy[1] - dimy[0]) / NY
dz = (dimz[1] - dimz[0]) / NZ

f = open("tracers_file", "w+")

for i in range(NX):
    for j in range(NY):
        for k in range(NZ):
            if (i==0 and j==0 and k==0):
                f.write("%d \r\n" % (NX*NY*NZ))
            else:
                f.write("%f %f %f \r\n" % ((i+1)*dx, (j+1)*dy, (k+1)*dz))
f.close()
                
