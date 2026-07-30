# Script to configure the software environment on Polaris (as of July 30, 2026)
# to successfully build the SUNDIALS+AMReX hands-on codes

module use /soft/modulefiles
module use /eagle/ATPESC2026/usr/modulefiles
module load track-5-numerical
module load conda
conda activate base
export MPICH_GPU_SUPPORT_ENABLED=0
