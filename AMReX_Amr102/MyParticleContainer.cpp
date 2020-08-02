
#include <MyParticleContainer.H>

namespace amrex {

void
MyParticleContainer::InitParticles (std::string initial_particle_file)
{
    BL_PROFILE("MyParticleContainer::InitParticles()");

    InitFromAsciiFile(initial_particle_file,0);

    int lev = 0;
    auto& pmap = GetParticles(lev);
    for (auto& kv : pmap) {
       int grid = kv.first.first;
       auto& pbox = kv.second.GetArrayOfStructs();
       const int n = pbox.size();

       for (int i = 0; i < n; i++)
       {
            ParticleType& p = pbox[i];

            // We over-write the z-locations to make sure they're in the domain
            p.pos(2) = 0.5 * m_geom.ProbHi(2);
       }
    }
}

}
