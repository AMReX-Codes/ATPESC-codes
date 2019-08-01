#include <algorithm>
#include <writeEBsurface.H>

using namespace amrex;

extern "C" void eb_to_polygon ( const amrex::Real * dx,
                                const int * slo, const int * shi,
                                const void * flag,        const int * fglo, const int * fghi,
                                const amrex_real * bcent, const int * blo,  const int * bhi,
                                const amrex_real * ax,    const int * axlo, const int * axhi,
                                const amrex_real * ay,    const int * aylo, const int * ayhi,
                                const amrex_real * az,    const int * azlo, const int * azhi   );

extern "C" void write_eb_vtp ( int * myID );
extern "C" void write_pvtp ( int * nProcs );

void WriteEBSurface (const BoxArray & ba, const DistributionMapping & dmap, const Geometry & geom,
                     const EBFArrayBoxFactory & ebf) {

    const Real * dx = geom.CellSize();

    MultiFab mf_ba(ba, dmap, 1, 0, MFInfo(), ebf);

    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto & sfab    = static_cast<EBFArrayBox const &>(mf_ba[mfi]);
        const auto & my_flag = sfab.getEBCellFlagFab();

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered or
            my_flag.getType(bx) == FabType::regular) continue;

        std::array<const MultiCutFab *, AMREX_SPACEDIM> areafrac;
        const MultiCutFab * bndrycent;

        areafrac  =   ebf.getAreaFrac();
        bndrycent = &(ebf.getBndryCent());

        eb_to_polygon(dx, BL_TO_FORTRAN_BOX(bx),
                      BL_TO_FORTRAN_3D(my_flag),
                      BL_TO_FORTRAN_3D((* bndrycent)[mfi]),
                      BL_TO_FORTRAN_3D((* areafrac[0])[mfi]),
                      BL_TO_FORTRAN_3D((* areafrac[1])[mfi]),
                      BL_TO_FORTRAN_3D((* areafrac[2])[mfi]) );
    }

    int cpu = ParallelDescriptor::MyProc();
    int nProcs = ParallelDescriptor::NProcs();

    write_eb_vtp(& cpu);

    if(ParallelDescriptor::IOProcessor())
        write_pvtp(& nProcs);
}
