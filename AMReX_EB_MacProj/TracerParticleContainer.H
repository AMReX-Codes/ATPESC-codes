#ifndef _TRACERPARTICLECONTAINER_H_
#define _TRACERPARTICLECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex
{

class TracerParticleContainer
  : public ParticleContainer<2*BL_SPACEDIM, 0>
//: public TracerParticleContainer
{
public:
    
    ///
    /// This particle container fills a mask for quickly computing
    /// neighbor grids / tiles for a given particle
    ///
    TracerParticleContainer(const Geometry            & geom, 
			    const DistributionMapping & dmap,
			    const BoxArray            & ba);
			    //    int                         nneighbor);
  /*
    TracerParticleContainer(const Vector<Geometry>            & geom, 
                                  const Vector<DistributionMapping> & dmap,
                                  const Vector<BoxArray>            & ba,
                                  const Vector<int>                 & rr,
                                  int                               nneighbor);
  */
    ///
    /// Init one particle per cell with random velocities
    ///
    void InitParticles();

    ///
    /// Compute the short range forces on a tile's worth of particles.
    /// fillNeighbors must have already been called.
    ///
    //void computeForces();

    ///
    /// Compute the short range forces on a tile's worth of particles using
    /// the neighbor list instead of the N^2 approach.
    /// fillNeighbors must have already been called.
    ///
    //void computeForcesNL();

    ///
    /// Move the particles according to their forces, reflecting at domain boundaries
    ///
  //    void moveParticles(const Real dt);

    ///
    /// Save the particle data in an ASCII format
    ///
    void writeParticles(int n);
    
  // void getNeighborParticles(const int nGrow);
};}    
    


#endif // _NEIGHBORLISTPARTICLECONTAINER_H_
