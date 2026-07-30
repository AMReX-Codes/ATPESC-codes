#!/bin/bash

set -e
set -o pipefail

# Testing
build_type=RelWithDebInfo
#PKG_INSTALL_ROOT=${HOME}/local/ATPESC2026
#EX_INSTALL_ROOT=${HOME}/local/ATPESC2026

# ATPESC
PKG_INSTALL_ROOT=/eagle/ATPESC2026
EX_INSTALL_ROOT=/eagle/ATPESC2026

# ------------------------------------------------------------------------------

PKG_INSTALL_ROOT=${PKG_INSTALL_ROOT}/usr/MathPackages/time_integration_sundials
EX_INSTALL_ROOT=${EX_INSTALL_ROOT}/EXAMPLES/track-5-numerical/time_integration_sundials

# ------------------------------------------------------------------------------

# set sundials install directory name
install_dir=${PKG_INSTALL_ROOT}/sundials

if [ ! -d ${install_dir} ]; then

    build_dir=${PWD}/build/sundials
    rm -rf ${build_dir}

    # configure
    cmake \
        -S ${PWD}/sundials \
        -B ${build_dir} \
        -D CMAKE_BUILD_TYPE=${build_type} \
        -D CMAKE_INSTALL_PREFIX=${install_dir} \
        -D CMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -D EXAMPLES_ENABLE_C=OFF \
        -D EXAMPLES_ENABLE_CXX=OFF \
        -D EXAMPLES_ENABLE_CUDA=OFF \
        -D EXAMPLES_INSTALL=OFF \
        -D BUILD_STATIC_LIBS=ON \
        -D BUILD_SHARED_LIBS=OFF \
        -D SUNDIALS_LOGGING_LEVEL=3 \
        -D ENABLE_MPI=ON \
        -D MPI_C_COMPILER=$(which mpicc) \
        -D MPI_CXX_COMPILER=$(which mpicxx) \
        -D SUNDIALS_PRECISION=double \
        -D SUNDIALS_INDEX_SIZE=32 \
        -D ENABLE_CUDA=ON \
        -D CMAKE_CUDA_ARCHITECTURES=80 \
        2>&1 | tee sundials-build.log

    # build
    cmake \
        --build ${build_dir} --parallel 16 \
        2>&1 | tee -a sundials-build.log

    # install
    cmake \
        --install ${build_dir} \
        2>&1 | tee -a sundials-build.log
fi

# ------------------------------------------------------------------------------

# set AMReX install directory name
install_dir=${PKG_INSTALL_ROOT}/amrex/cpu

if [ ! -d ${install_dir} ]; then

    build_dir=${PWD}/build/amrex/cpu
    rm -rf ${build_dir}

    # configure
    cmake \
        -S ${PWD}/amrex \
        -B ${build_dir} \
        -D CMAKE_BUILD_TYPE=${build_type} \
        -D CMAKE_INSTALL_PREFIX=${install_dir} \
        -D AMReX_SPACEDIM=2 \
        -D AMReX_PLOTFILE_TOOLS=ON \
        -D AMReX_MPI=ON \
        -D MPI_C_COMPILER=$(which mpicc) \
        -D MPI_CXX_COMPILER=$(which mpicxx) \
        2>&1 | tee amrex-build-cpu.log

    # build
    cmake \
        --build ${build_dir} --parallel 16 \
        2>&1 | tee -a amrex-build-cpu.log

    # install
    cmake \
        --install ${build_dir} \
        2>&1 | tee -a amrex-build-cpu.log
fi

# ------------------------------------------------------------------------------

# set AMReX install directory name
install_dir=${PKG_INSTALL_ROOT}/amrex/gpu

if [ ! -d ${install_dir} ]; then

    build_dir=${PWD}/build/amrex/gpu
    rm -rf ${build_dir}

    # configure
    cmake \
        -S ${PWD}/amrex \
        -B ${build_dir} \
        -D CMAKE_BUILD_TYPE=${build_type} \
        -D CMAKE_INSTALL_PREFIX=${install_dir} \
        -D AMReX_SPACEDIM=2 \
        -D AMReX_TINY_PROFILE=ON \
        -D AMReX_PLOTFILE_TOOLS=ON \
        -D AMReX_MPI=ON \
        -D MPI_C_COMPILER=$(which mpicc) \
        -D MPI_CXX_COMPILER=$(which mpicxx) \
        -D AMReX_GPU_BACKEND=CUDA \
        -D CMAKE_CUDA_ARCHITECTURES=80 \
        -D AMReX_SUNDIALS=ON \
        -D SUNDIALS_DIR=${PKG_INSTALL_ROOT}/sundials/lib64/cmake/sundials \
        2>&1 | tee amrex-build-gpu.log

    # build
    cmake \
        --build ${build_dir} --parallel 16 \
        2>&1 | tee -a amrex-build-gpu.log

    # install
    cmake \
        --install ${build_dir} \
        2>&1 | tee -a amrex-build-gpu.log
fi

# ------------------------------------------------------------------------------

# set the build and install directory name
build_dir=${PWD}/build/atpesc_examples
install_dir=${EX_INSTALL_ROOT}

# remove old build directory
rm -rf ${build_dir}
rm -rf ${install_dir}

# configure build
cmake \
    -S ${PWD}/ATPESC-codes/SUNDIALS+AMReX \
    -B ${build_dir} \
    -D CMAKE_BUILD_TYPE=${build_type} \
    -D CMAKE_INSTALL_PREFIX=${install_dir} \
    -D MPI_C_COMPILER=$(which mpicc) \
    -D MPI_CXX_COMPILER=$(which mpicxx) \
    -D ENABLE_CUDA=ON \
    -D CMAKE_CUDA_ARCHITECTURES=80 \
    -D AMREX_ROOT=${PKG_INSTALL_ROOT}/amrex/gpu \
    -D SUNDIALS_ROOT=${PKG_INSTALL_ROOT}/sundials \
    -D INSTALL_SOURCE=OFF \
    2>&1 | tee atpesc-build.log

# build
cmake \
    --build ${build_dir} \
    2>&1 | tee -a atpesc-build.log

# install
cmake \
    --install ${build_dir} \
    2>&1 | tee -a atpesc-build.log

# ------------------------------------------------------------------------------

# Copy suntools
cp -r ${PWD}/sundials/suntools/src/suntools ${EX_INSTALL_ROOT}/.

# Copy amrex_fcompare
cp ${PKG_INSTALL_ROOT}/amrex/cpu/bin/amrex_fcompare ${EX_INSTALL_ROOT}/.

# Update group (just in case)
chgrp -R ATPESC_Instructors ${PKG_INSTALL_ROOT}
chgrp -R ATPESC_Instructors ${EX_INSTALL_ROOT}

# Update file permissons
chmod -R go=u ${PKG_INSTALL_ROOT}
chmod -R o-w ${PKG_INSTALL_ROOT}

chmod -R go=u ${EX_INSTALL_ROOT}
chmod -R o-w ${EX_INSTALL_ROOT}
