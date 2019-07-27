
#include <MyParticleContainer.H>

namespace amrex {

Real 
MyParticleContainer::FindWinner (int n)
{
    BL_PROFILE("MyParticleContainer::FindWinner()");

    using ParticleType = MyParticleContainer::ParticleType;
    int nghost = 0;
    Real x = amrex::ReduceMax(*this, nghost,
       [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> Real
                                  { return p.pos(n); });

    ParallelDescriptor::ReduceRealMax(x);
    return x;
}

}
