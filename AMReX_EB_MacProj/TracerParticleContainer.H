#ifndef _TRACERPARTICLECONTAINER_H_
#define _TRACERPARTICLECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex
{

class TracerParticleContainer
  : public TracerParticleContainer<2*BL_SPACEDIM, 0>
{
public:
    
    TracerParticleContainer(const Geometry            & geom, 
			    const DistributionMapping & dmap,
			    const BoxArray            & ba);
    void InitParticles();

    void writeParticles(int n);
    

};}    
    


#endif 
