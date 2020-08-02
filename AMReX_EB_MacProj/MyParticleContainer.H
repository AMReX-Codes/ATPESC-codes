#ifndef BL_MYPARTICLES_H_
#define BL_MYPARTICLES_H_

#include <AMReX_TracerParticles.H>

namespace InterpolationType {
    enum {NGP=0, CIC};
}

namespace PIdx {
    enum {BufferX=0, BufferY, BufferZ, Weight};
}

namespace amrex {

class MyParticleContainer
    : public ParticleContainer<1>
{
private:
    int m_number_particles_per_cell;

public:

    MyParticleContainer (ParGDBBase* gdb)
        : ParticleContainer<1>(gdb),
          m_number_particles_per_cell(0)
        {}

    MyParticleContainer (const Geometry            & geom,
                         const DistributionMapping & dmap,
                         const BoxArray            & ba)
        : ParticleContainer<1>(geom,dmap,ba),
          m_number_particles_per_cell(0)
        {}

    ~MyParticleContainer () {}

    void InitParticles(int nppc);

    void AdvectWithUmac (MultiFab* umac, int lev, Real dt);

    void DepositToMesh (MultiFab& phi, int interpolation=InterpolationType::NGP);
};

}

#endif

