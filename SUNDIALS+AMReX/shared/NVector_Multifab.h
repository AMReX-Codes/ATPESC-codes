/*--------------------------------------------------------------------
  Time Integration and Nonlinear Solvers
  Hands-on Lessons with SUNDIALS + AMReX
  2019 Argonne Training Program in Extreme-Scale Computing

  Authors (alphabetical):
    David Gardner (gardner48@llnl.gov)
    John Loffeld (loffeld1@llnl.gov)
    Daniel Reynolds (reynolds@smu.edu)
    Donald Willcox (dewillcox@lbl.gov)

  --------------------------------------------------------------------
  Header file for N_Vector wrap of AMReX 'MultiFab' structure.
  --------------------------------------------------------------------*/

#ifndef NVECTOR_MULTIFAB_H
#define NVECTOR_MULTIFAB_H

#include <stdio.h>
#include <sundials/sundials_nvector.h>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

/* -----------------------------------------------------------------
 * Vector content structure
 * -----------------------------------------------------------------*/

struct _N_VectorContent_Multifab {
  sunindextype length;   /* vector length           */
  booleantype own_mf;    /* MultiFab ownership flag */
  amrex::MultiFab *mf;   /* wrapped MultiFab        */
};

typedef struct _N_VectorContent_Multifab *N_VectorContent_Multifab;

/* -----------------------------------------------------------------
 * Accessor Macros
 * -----------------------------------------------------------------*/

#define NV_CONTENT_M(v) ( (N_VectorContent_Multifab)(v->content) )

#define NV_LENGTH_M(v)  ( NV_CONTENT_M(v)->length )

#define NV_OWN_MF_M(v)  ( NV_CONTENT_M(v)->own_mf )

#define NV_MFAB(v)      ( NV_CONTENT_M(v)->mf )

/* -----------------------------------------------------------------
 * Exported functions
 * -----------------------------------------------------------------*/

N_Vector N_VNewEmpty_Multifab(sunindextype vec_length);
N_Vector N_VNew_Multifab(sunindextype vec_length,
                         const amrex::BoxArray &ba,
                         const amrex::DistributionMapping &dm,
                         sunindextype nComp,
                         sunindextype nGhost);
N_Vector N_VMake_Multifab(sunindextype vec_length,
                          amrex::MultiFab *mf);
sunindextype N_VGetLength_Multifab(N_Vector v);
N_Vector N_VCloneEmpty_Multifab(N_Vector w);
N_Vector N_VClone_Multifab(N_Vector w);
void N_VDestroy_Multifab(N_Vector v);
void N_VSpace_Multifab(N_Vector v, sunindextype *lrw,
                       sunindextype *liw);

/* standard vector operations */
void N_VLinearSum_Multifab(realtype a, N_Vector x,
                           realtype b, N_Vector y, N_Vector z);
void N_VConst_Multifab(realtype c, N_Vector z);
void N_VProd_Multifab(N_Vector x, N_Vector y, N_Vector z);
void N_VDiv_Multifab(N_Vector x, N_Vector y, N_Vector z);
void N_VScale_Multifab(realtype c, N_Vector x, N_Vector z);
void N_VAbs_Multifab(N_Vector x, N_Vector z);
void N_VInv_Multifab(N_Vector x, N_Vector z);
void N_VAddConst_Multifab(N_Vector x, realtype b, N_Vector z);
realtype N_VDotProd_Multifab(N_Vector x, N_Vector y);
realtype N_VMaxNorm_Multifab(N_Vector x);
realtype N_VWrmsNorm_Multifab(N_Vector x, N_Vector w);
realtype N_VWrmsNormMask_Multifab(N_Vector x, N_Vector w,
                                  N_Vector id);
realtype N_VMin_Multifab(N_Vector x);
realtype N_VWL2Norm_Multifab(N_Vector x, N_Vector w);
realtype N_VL1Norm_Multifab(N_Vector x);
void N_VCompare_Multifab(realtype c, N_Vector x, N_Vector z);
booleantype N_VInvTest_Multifab(N_Vector x, N_Vector z);
booleantype N_VConstrMask_Multifab(N_Vector c, N_Vector x,
                                   N_Vector m);
realtype N_VMinQuotient_Multifab(N_Vector num, N_Vector denom);

#ifdef __cplusplus
}
#endif

#endif
