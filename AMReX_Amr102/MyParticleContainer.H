#ifndef BL_MYPARTICLES_H_
#define BL_MYPARTICLES_H_

#include <AMReX_Particles.H>
#include <AMReX_ParticleMesh.H>
#include <AMReX_TracerParticle_mod_K.H>
#include <Indexing.H>

namespace Interpolation {
    enum {NGP=0, CIC};
}

namespace PIdx {
    enum {AMREX_D_DECL(BufferX=0, BufferY, BufferZ), Weight, NStructReal};
}

namespace amrex {

class MyParticleContainer
    : public ParticleContainer<PIdx::NStructReal>
{
private:
    int m_number_particles_per_cell;

public:

    MyParticleContainer (ParGDBBase* gdb)
        : ParticleContainer<PIdx::NStructReal>(gdb),
          m_number_particles_per_cell(0)
        {}

    MyParticleContainer (const Geometry            & geom,
                         const DistributionMapping & dmap,
                         const BoxArray            & ba)
        : ParticleContainer<PIdx::NStructReal>(geom,dmap,ba),
          m_number_particles_per_cell(0)
        {}

    ~MyParticleContainer () {}

    void InitParticles(int nppc, const MultiFab& phi, const MultiFab& ebvol);

    int NumParticlesPerCell() { return m_number_particles_per_cell; }

    void AdvectWithUmac (MultiFab* umac, int lev, Real dt);

    void DepositToMesh (MultiFab& phi, int interpolation=Interpolation::NGP);

    void InterpolateFromMesh (const MultiFab& phi, int interpolation=Interpolation::NGP);

    void RemoveCoveredParticles (const MultiFab& ebvol, int interpolation=Interpolation::NGP);
};

}

#endif

