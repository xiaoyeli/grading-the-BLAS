#! /bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
VERSION_NUMBER=3.12.1

cd ${REPO_ROOT}/submodules

# clean
rm -rf lapack-${VERSION_NUMBER}

# acquire src
if ! wget https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v${VERSION_NUMBER}.tar.gz; then
    exit 1
fi
tar -xvf v${VERSION_NUMBER}.tar.gz
rm v${VERSION_NUMBER}.tar.gz


cd ${REPO_ROOT}/submodules/lapack-${VERSION_NUMBER}

make clean

cp INSTALL/make.inc.gfortran make.inc

# build tweaks for shared library + march=native
sed -i "/FFLAGS = -O2 -frecursive/c\FFLAGS = -O2 -frecursive -march=native -fPIC" make.inc
sed -i "s/librefblas.a/librefblas.so/" make.inc
sed -i "s|\$(AR) \$(ARFLAGS) \$@ \$^|gfortran -z muldefs -shared -Wl,-soname,librefblas.so -o \$@ \$^|" BLAS/SRC/Makefile
sed -i "/\$(RANLIB) \$@/d" BLAS/SRC/Makefile

# minimal build of LAPACK/BLAS
make -j 8 blaslib