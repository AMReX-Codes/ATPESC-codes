void ComputeAdvection(amrex::MultiFab& sol,
                      amrex::MultiFab& advection,
                      amrex::Geometry& geom,
                      int comp, amrex::Real advCoeff);

void ComputeDiffFlux(amrex::MultiFab& sol,
                     amrex::MultiFab& fx,
                     amrex::MultiFab& fy,
                     amrex::Geometry& geom,
                     int comp, amrex::Real diffCoeff);

void ComputeDivergence(amrex::MultiFab& div,
                       amrex::MultiFab& fx,
                       amrex::MultiFab& fy,
                       amrex::Geometry& geom,
                       int comp);
