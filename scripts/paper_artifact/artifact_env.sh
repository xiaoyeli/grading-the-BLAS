#!/bin/bash

set -v
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_NUM_THREADS=64
set +v

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

cd ${REPO_ROOT}/src

make -f Makefile.libs

NAIVE_TRIPLE_LOOP_SO_PATH=$(realpath libnaive_triple_loop_dgemm.so)
if [ -e "${NAIVE_TRIPLE_LOOP_SO_PATH}" ]
then
    echo "Found libnaive_triple_loop_dgemm.so at ${NAIVE_TRIPLE_LOOP_SO_PATH}"
    export NAIVE_TRIPLE_LOOP_SO_PATH
else
    echo "Warning: Could not find libnaive_triple_loop_dgemm.so"
fi

STRASSEN_SO_PATH=$(realpath libstrassen_dgemm.so)
if [ -e "${STRASSEN_SO_PATH}" ]
then
    echo "Found libstrassen_dgemm.so at ${STRASSEN_SO_PATH}"
    export STRASSEN_SO_PATH
else
    echo "Warning: Could not find libstrassen_dgemm.so"
fi

STRASSEN_ZEROING_SO_PATH=$(realpath libstrassen_zeroing_dgemm.so)
if [ -e "${STRASSEN_ZEROING_SO_PATH}" ]
then
    echo "Found libstrassen_zeroing_dgemm.so at ${STRASSEN_ZEROING_SO_PATH}"
    export STRASSEN_ZEROING_SO_PATH
else
    echo "Warning: Could not find libstrassen_zeroing_dgemm.so"
fi

cd ${REPO_ROOT}

NETLIBBLAS_LIB_PATH=$(find ${REPO_ROOT}/submodules -name "librefblas.a" | head -n 1 | xargs realpath)
if [ -e "${NETLIBBLAS_LIB_PATH}" ]
then
    echo "Found librefblas.a at ${NETLIBBLAS_LIB_PATH}"
    export NETLIBBLAS_LIB_PATH
else
    echo "Warning: Could not find librefblas.a"
fi

BLIS_SO_PATH=$(find ${REPO_ROOT}/submodules -name "libblis.so" | head -n 1 | xargs realpath)
if [ -e "${BLIS_SO_PATH}" ]
then
    echo "Found libblis.so at ${BLIS_SO_PATH}"
    export BLIS_SO_PATH
else
    echo "Warning: Could not find libblis.so"
fi

OPENBLAS_SO_PATH=$(find ${REPO_ROOT}/submodules -name "libopenblas.so" | head -n 1 | xargs realpath)
if [ -e "${OPENBLAS_SO_PATH}" ]
then
    echo "Found libopenblas.so at ${OPENBLAS_SO_PATH}"
    export OPENBLAS_SO_PATH
else
    echo "Warning: Could not find libopenblas.so"
fi

CUBLAS_SO_PATH=$(find ${CUDA_HOME}/../../math_libs -name "libcublas.so" | head -n 1 | xargs realpath)
if [ -e "${CUBLAS_SO_PATH}" ]
then
    echo "Found libcublas.so at ${CUBLAS_SO_PATH}"
    export CUBLAS_SO_PATH
else
    echo "Warning: Could not find libcublas.so"
fi

OZIMMU_SO_PATH=$(find ${REPO_ROOT}/submodules -name "libozimmu.so" | head -n 1 | xargs realpath)
if [ -e "${OZIMMU_SO_PATH}" ]
then
    echo "Found libozimmu.so at ${OZIMMU_SO_PATH}"
    export OZIMMU_SO_PATH
else
    echo "Warning: Could not find libozimmu.so"
fi

GEMMUL8_SO_PATH=$(find ${REPO_ROOT}/submodules -name "libgemmul8.so" | head -n 1 | xargs realpath)
if [ -e "${GEMMUL8_SO_PATH}" ]
then
    echo "Found libgemmul8.so at ${GEMMUL8_SO_PATH}"
    export GEMMUL8_SO_PATH
else
    echo "Warning: Could not find libgemmul8.so"
fi