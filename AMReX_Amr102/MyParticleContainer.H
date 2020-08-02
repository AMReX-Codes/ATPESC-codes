#ifndef BL_MYPARTICLES_H_
#define BL_MYPARTICLES_H_

#include <AMReX_TracerParticles.H>

namespace amrex {

class MyParticleContainer
    : public TracerParticleContainer
{
public:

    MyParticleContainer (ParGDBBase* gdb)
        : TracerParticleContainer(gdb)
        {}

    MyParticleContainer (const Geometry            & geom,
                         const DistributionMapping & dmap,
                         const BoxArray            & ba)
        : TracerParticleContainer(geom,dmap,ba), m_ba(ba), m_geom(geom)
        {}

    ~MyParticleContainer () {}

    void InitParticles (std::string initial_particle_file);

private:

     BoxArray m_ba;
     Geometry m_geom;

};

}

#endif

