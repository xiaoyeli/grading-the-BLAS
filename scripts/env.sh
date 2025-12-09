#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}

NETLIBBLAS_SO_PATH=$(find -name "librefblas.so" | head -n 1 | xargs realpath)
if [ -e "${NETLIBBLAS_SO_PATH}" ]
then
    echo "Found librefblas.so at ${NETLIBBLAS_SO_PATH}"
    NETLIBBLAS_LIB_PATH=$(dirname ${NETLIBBLAS_SO_PATH})
    export LIBRARY_PATH=${LIBRARY_PATH}:${NETLIBBLAS_LIB_PATH}
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${NETLIBBLAS_LIB_PATH}
else
    echo "Warning: Could not find librefblas.so"
fi

BLIS_SO_PATH=$(find -name "libblis.so" | head -n 1 | xargs realpath)
BLIS_H_PATH=$(find -name "blis.h" | head -n 1 | xargs realpath)
if [ -e "${BLIS_SO_PATH}" ]
then
    echo "Found libblis.so at ${BLIS_SO_PATH}"
    BLIS_LIB_PATH=$(dirname ${BLIS_SO_PATH})
    export LIBRARY_PATH=${LIBRARY_PATH}:${BLIS_LIB_PATH}
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BLIS_LIB_PATH}
else
    echo "Warning: Could not find libblis.so"
fi
if [ -e "${BLIS_H_PATH}" ]
then
    echo "Found blis.h at ${BLIS_H_PATH}"
    export CPATH=${CPATH}:$(dirname ${BLIS_H_PATH})
else
    echo "Warning: Could not find blis.h"
fi

OPENBLAS_SO_PATH=$(find -name "libopenblas.so" | head -n 1 | xargs realpath)
if [ -e "${OPENBLAS_SO_PATH}" ]
then
    echo "Found libopenblas.so at ${OPENBLAS_SO_PATH}"
    OPENBLAS_LIB_PATH=$(dirname ${OPENBLAS_SO_PATH})
    export LIBRARY_PATH=${LIBRARY_PATH}:${OPENBLAS_LIB_PATH}
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${OPENBLAS_LIB_PATH}
else
    echo "Warning: Could not find libopenblas.so"
fi

CUBLAS_SO_PATH=$(find ${CUDA_HOME}/../../math_libs -name "libcublas.so" | head -n 1 | xargs realpath)
if [ -e "${CUBLAS_SO_PATH}" ]
then
    echo "Found libcublas.so at ${CUBLAS_SO_PATH}"
    CUBLAS_LIB_PATH=$(dirname ${CUBLAS_SO_PATH})
    export LIBRARY_PATH=${LIBRARY_PATH}:${CUBLAS_LIB_PATH}
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CUBLAS_LIB_PATH}
else
    echo "Warning: Could not find libcublas.so"
fi

OZIMMU_SO_PATH=$(find -name "libozimmu.so" | head -n 1 | xargs realpath)
if [ -e "${OZIMMU_SO_PATH}" ]
then
    echo "Found libozimmu.so at ${OZIMMU_SO_PATH}"
    export OZIMMU_SO_PATH=${OZIMMU_SO_PATH}
else
    echo "Warning: Could not find libozimmu.so"
fi

set -v
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_NUM_THREADS=64
export OZIMMU_COMPUTE_MODE=fp64_int8_auto
export OZIMMU_AUTO_AVG_MANTISSA_LOSS_THRESHOLD=0
export OZIMMU_INTERCEPT_THRESHOLD_M=32
export OZIMMU_INTERCEPT_THRESHOLD_N=32
export OZIMMU_INTERCEPT_THRESHOLD_K=32
set +v