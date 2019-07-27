
#include <MyParticleContainer.H>

namespace amrex {

void 
MyParticleContainer::InitPachinko (std::string initial_tracer_file)
{
    BL_PROFILE("MyParticleContainer::InitPachinko()");

    InitFromAsciiFile(initial_tracer_file,0);

    int lev = 0;
    auto& pmap = GetParticles(lev);
    for (auto& kv : pmap) {
       int grid = kv.first.first;
       auto& pbox = kv.second.GetArrayOfStructs();
       const int n = pbox.size();

       for (int i = 0; i < n; i++)
       {
            ParticleType& p = pbox[i];

            p.rdata(PIdx::vx) =  0.0;
            p.rdata(PIdx::vy) = -1.0;
       }
    }
}

void 
MyParticleContainer::AdvectPachinko (Real dt)
{
    BL_PROFILE("MyParticleContainer::AdvectPachinko()");

    int lev = 0;

    Real xcent[14];
    Real ycent[14];

    xcent[0] = 0.3;
    ycent[0] = 0.3;

    xcent[1] = 0.6;
    ycent[1] = 0.3;

    xcent[2] = 0.9;
    ycent[2] = 0.3;

    xcent[3] = 0.15;
    ycent[3] = 0.7;

    xcent[4] = 0.45;
    ycent[4] = 0.7;

    xcent[5] = 0.75;
    ycent[5] = 0.7;

    xcent[6] = 1.05;
    ycent[6] = 0.7;

    xcent[7] = 0.3;
    ycent[7] = 1.1;

    xcent[8] = 0.6;
    ycent[8] = 1.1;

    xcent[9] = 0.9;
    ycent[9] = 1.1;

    xcent[10] = 0.15;
    ycent[10] = 1.5;

    xcent[11] = 0.45;
    ycent[11] = 1.5;

    xcent[12] = 0.75;
    ycent[12] = 1.5;

    xcent[13] = 1.05;
    ycent[13] = 1.5;

    // Particle radius
    Real prad = 0.02 ;

    // Obstracle radiu
    Real crad = 0.10 ;

    Real rad_squared = (prad+crad)*(prad+crad);

    Real grav = -10.;

    const auto prob_lo = Geom(0).ProbLoArray();
    const auto prob_hi = Geom(0).ProbHiArray();

    auto& pmap = GetParticles(lev);
    for (auto& kv : pmap) {
       int grid = kv.first.first;
       auto& pbox = kv.second.GetArrayOfStructs();
       const int n = pbox.size();

       // std::cout << "Number of particles: " << n << std::endl;

       for (int i = 0; i < n; i++)
       {
          ParticleType& p = pbox[i];

          // Let particles stop at the bottom
          if (p.pos(1) >= prob_lo[1] + 0.05) 
          {
            p.pos(0) += 0.5 * dt * p.rdata(PIdx::vx);
            p.pos(1) += 0.5 * dt * p.rdata(PIdx::vy);

            // Accleration under graivty
            p.rdata(PIdx::vy) += dt * grav;

            p.pos(0) += 0.5 * dt * p.rdata(PIdx::vx);
            p.pos(1) += 0.5 * dt * p.rdata(PIdx::vy);

            p.pos(0) += dt * p.rdata(PIdx::vx);
            p.pos(1) += dt * p.rdata(PIdx::vy);

            for (int ind = 0; ind < 14; ind++)
            {
               Real x_diff = p.pos(0) - xcent[ind];
               Real y_diff = p.pos(1) - ycent[ind];
               Real diff_sq = x_diff * x_diff + y_diff * y_diff;
               Real diff    = std::sqrt(diff_sq);

               if ( diff_sq < rad_squared )
               {

                 Real norm_x = x_diff / diff;
                 Real norm_y = y_diff / diff;

                 Real tang_x = -norm_y;
                 Real tang_y =  norm_x;

                 // Incoming velocity dot normal  = (norm_x, norm_y)
                 Real vel_norm = p.rdata(PIdx::vx) * norm_x + 
                                 p.rdata(PIdx::vy) * norm_y;

                 // Incoming velocity dot tangent = (x_tang, y_tang) = (-y_norm, x_norm)
                 Real vel_tang = p.rdata(PIdx::vx) * norm_y + 
                                 p.rdata(PIdx::vy) * norm_x;

                 // Original velocity was (vel_norm) * (norm_x, norm_y)
                 //                    +  (vel_tang) * (tang_x, tang_y)

                 // New      velocity  is MINUS (vel_norm) * (norm_x, norm_y)
                 //                          +  (vel_tang) * (tang_x, tang_y)

                 p.rdata(PIdx::vx) = -vel_norm * norm_x + vel_tang * tang_x;
                 p.rdata(PIdx::vy) = -vel_norm * norm_y + vel_tang * tang_y;

               }
            }

            // Bounce off left wall
            if (p.pos(0) < prob_lo[0]) 
            {
               p.pos(0) = -p.pos(0);
               p.rdata(PIdx::vx) = -p.rdata(PIdx::vx);
            }

            // Bounce off right wall
            if (p.pos(0) > prob_hi[0]) 
            {
               p.pos(0) = -p.pos(0);
               p.rdata(PIdx::vx) = -p.rdata(PIdx::vx);
            }

            // Stick to bottom wall
            // if (p.pos(1) < 0.05) 
            // {
            //    p.rdata(PIDx::vx) = 0.0;
            //    p.rdata(PIDx::vy) = 0.0;
            // }
         }
       }
    }
}

Real 
MyParticleContainer::FindWinner (int n)
{
    BL_PROFILE("MyParticleContainer::FindWinner()");

//  using ParticleType = MyParticleContainer::ParticleType;
    int nghost = 0;
    if (n != 0 && n != 1) amrex::Abort("Must specify 0 or 1 in FindWinner");

    Real x = amrex::ReduceMax(*this, nghost,
       [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> Real
                                  { return p.pos(n); });

    ParallelDescriptor::ReduceRealMax(x);
    return x;
}

}
