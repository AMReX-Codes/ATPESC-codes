#!/bin/bash

set -e
set -o pipefail

if [ ! -d ${PWD}/sundials ]; then
    git clone git@github.com:LLNL/sundials.git
    cd sundials
    git checkout v7.4.0
    cd -
fi

if [ ! -d ${PWD}/amrex ]; then
    git clone git@github.com:AMReX-Codes/amrex.git
    cd amrex
    git checkout 25.07
    cd -
fi

if [ ! -d ${PWD}/ATPESC-codes ]; then
    git clone git@github.com:AMReX-Codes/ATPESC-codes.git
fi
