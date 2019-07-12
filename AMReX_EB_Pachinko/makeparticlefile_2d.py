#This will create a tracer particle input file based on the number
#of particles in each direction (NX,NY,NZ) and the physical domain
#boundaries (dimx,dimy,dimz). Amrex can then read this file in with
#the function TracerParticleContainer::InitFromAsciiFile()

import numpy as np

NX = 5
NY = 100

dimx = [0.01,0.1]  #They'll start evenly close to the edge/source and flow through
dimy = [0.01,0.99]

dx = (dimx[1] - dimx[0]) / NX
dy = (dimy[1] - dimy[0]) / NY

f = open("tracers_file_2d", "w")

for i in range(NX):
    for j in range(NY):
        if (i==0 and j==0):
            f.write("%d \r\n" % (NX*NY))
            f.write("%f %f \r\n" % ((i+1)*dx, (j+1)*dy))
        else:
            f.write("%f %f \r\n" % ((i+1)*dx, (j+1)*dy))
f.close()
                