
#include <MyParticleContainer.H>

namespace amrex {

void
MyParticleContainer::FindWinner (int lev, amrex::Real& x)
{
    BL_PROFILE("TracerParticleContainer::FindWinner()");
    BL_ASSERT(lev >= 0 && lev < GetParticles().size());

    x = -1.e10;
    auto& pmap = GetParticles(lev);
    for (auto& kv : pmap) 
    {
        int grid = kv.first.first;
        auto& pbox = kv.second.GetArrayOfStructs();
	const int n = pbox.size();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            x = std::max(x, p.m_rdata.pos[0]);

        }
    }
}

}

