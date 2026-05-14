#! /bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
VERSION_NUMBER=0.3.30

cd ${REPO_ROOT}/submodules

# clean
rm -rf OpenBLAS-${VERSION_NUMBER}

# acquire src
if ! wget https://github.com/OpenMathLib/OpenBLAS/archive/refs/tags/v${VERSION_NUMBER}.tar.gz; then
    exit 1
fi
tar -xvf v${VERSION_NUMBER}.tar.gz
rm v${VERSION_NUMBER}.tar.gz


cd ${REPO_ROOT}/submodules/OpenBLAS-${VERSION_NUMBER}
make clean    
make -j 8 TARGET=SANDYBRIDGE FC=gfortran CC=gcc NO_LAPACK=1