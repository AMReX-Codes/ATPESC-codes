#ifndef _TRACERPARTICLECONTAINER_H_
#define _TRACERPARTICLECONTAINER_H_

#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_Particle.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

TracerParticleContainer::
TracerParticleContainer(const Geometry            & a_geom,
			const DistributionMapping & a_dmap,
			const BoxArray            & a_ba)
  : ParticleContainer<RealData::ncomps, IntData::ncomps> (a_geom, a_dmap, a_ba)
{}

void TracerParticleContainer::InitParticles() {

  
  BL_PROFILE("TracerParticleContainer::InitParticles");
      
  const int lev = 0;
  const Geometry& geom = Geom(lev);
  const Real* dx = geom.CellSize();
  const Real* plo = geom.ProbLo();
  
  std::mt19937 mt(0451);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  
  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    const Box& tile_box = mfi.tilebox();
    const RealBox tile_real_box { tile_box, dx, geom.ProbLo() };
    
    const int grid_id = mfi.index();
    const int tile_id = mfi.LocalTileIndex();
    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];
    
    const auto& boxlo = tile_box.smallEnd();
    ParticleType p;
    for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
      
      p.id() = ParticleType::NextID();
      p.cpu() = ParallelDescriptor::MyProc();
      
      p.pos(0) = tile_real_box.lo(0) + (iv[0]- boxlo[0] + 0.5)*dx[0];
      p.pos(1) = tile_real_box.lo(1) + (iv[1]- boxlo[1] + 0.5)*dx[1];
#if (BL_SPACEDIM == 3)
      p.pos(2) = tile_real_box.lo(2) + (iv[2]- boxlo[2] + 0.5)*dx[2];
#endif
      p.rdata(0) = dist(mt);
      p.rdata(1) = dist(mt);
#if (BL_SPACEDIM == 3)
      p.rdata(2) = dist(mt);
#endif
      
      p.rdata(BL_SPACEDIM)   = 0;
      p.rdata(BL_SPACEDIM+1) = 0;
#if (BL_SPACEDIM == 3)
      p.rdata(BL_SPACEDIM+2) = 0;
#endif
      
      particle_tile.push_back(p);
    }
  }
  
  Redistribute();
 
}





void TracerParticleContainer::writeParticles(int n)
{
    BL_PROFILE("TracerParticleContainer::writeParticles");
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}











  /*
      if ( ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber() ) {

	ParticleType p;

	p.id() = ParticleType::NextID();
	p.cpu() = ParallelDescriptor::MyProc();

	p.pos(0) = 0;
	p.pos(1) = 0.5;

	std::array<Real,2> attribs;
	attribs[Real vx] = 0.0;
	attribs[Real vy] = 0.0;

        // Add to level 0, grid 0, and tile 0
        // Redistribute() will move it to the proper place.  
	std::pair<int,int> key {0,0};
        auto& particle_tile = GetParticles(0)[key];

        particle_tile.push_back(p);
        particle_tile.push_back_real(attribs);

      }

      Redistribute();
      
  */


  /*
void
init_particles ()
{
  //if (level == 0)
  // {
      BL_ASSERT(TracerPC == nullptr);

      TracerPC.reset(new AmrTracerParticleContainer(parent));
      TracerPC->do_tiling = true;
      TracerPC->tile_size = IntVect(AMREX_D_DECL(1024000,4,4));

      AmrTracerParticleContainer::ParticleInitData pdata = {AMREX_D_DECL(0.0, 0.0, 0.0)};

      TracerPC->SetVerbose(0);
      TracerPC->InitOnePerCell(0.5, 0.5, 0.5, pdata);

      TracerPC->Redistribute();
      //    }


      */


#endif
