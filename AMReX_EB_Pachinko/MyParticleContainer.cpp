
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

    Real xcent[9];
    Real ycent[9];

    xcent[0] = 0.2;
    ycent[0] = 0.3;

    xcent[1] = 0.5;
    ycent[1] = 0.3;

    xcent[2] = 0.8;
    ycent[2] = 0.8;

    xcent[3] = 0.2;
    ycent[3] = 0.3;

    xcent[4] = 0.5;
    ycent[4] = 0.3;

    xcent[5] = 0.8;
    ycent[5] = 0.8;

    xcent[6] = 0.2;
    ycent[6] = 1.1;

    xcent[7] = 0.5;
    ycent[7] = 1.1;

    xcent[8] = 0.8;
    ycent[8] = 1.1;

    Real rad_squared = 0.1 * 0.1;

    auto& pmap = GetParticles(lev);
    for (auto& kv : pmap) {
       int grid = kv.first.first;
       auto& pbox = kv.second.GetArrayOfStructs();
       const int n = pbox.size();

       for (int i = 0; i < n; i++)
       {
            ParticleType& p = pbox[i];

            p.pos(0) += dt * p.rdata(PIdx::vx);
            p.pos(1) += dt * p.rdata(PIdx::vy);

            for (int ind = 0; ind < 9; ind++)
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

#if 0
                if (i == 0) 
                {
                   std::cout << "center " << xcent[ind] << " " << ycent[ind] << std::endl;
                   std::cout << "normal " << norm_x << " " << norm_y << std::endl;
                   std::cout << "Vel norm " << vel_norm << std::endl;
                   std::cout << "Vel tang " << vel_tang << std::endl;
                   std::cout << "Vel bef " << p.rdata(PIdx::vx) << " " << p.rdata(PIdx::vy) << std::endl;
                }
#endif

                 p.rdata(PIdx::vx) = -vel_norm * norm_x + vel_tang * tang_x;
                 p.rdata(PIdx::vy) = -vel_norm * norm_y + vel_tang * tang_y;

#if 0
                if (i == 0) 
                   std::cout << "Vel aft " << p.rdata(PIdx::vx) << " " << p.rdata(PIdx::vy) << "\n" << std::endl;
#endif

               }
            }
       }
    }
}

Real 
MyParticleContainer::FindWinner ()
{
    BL_PROFILE("MyParticleContainer::FindWinner()");

//  using ParticleType = MyParticleContainer::ParticleType;
    int nghost = 0;
    Real x = amrex::ReduceMax(*this, nghost,
       [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> Real
                                  { return p.pos(0); });

    ParallelDescriptor::ReduceRealMax(x);
    return x;
}

}
