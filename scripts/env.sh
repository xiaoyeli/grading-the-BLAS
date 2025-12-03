#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}

PATH_OPENBLAS_LIB=$(find -name "libopenblas.so" | head -n 1 | xargs realpath)
if [ -e "${PATH_OPENBLAS_LIB}" ]
then
    echo "Found libopenblas.so at ${PATH_OPENBLAS_LIB}; adding to LD_LIBRARY_PATH"
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(dirname ${PATH_OPENBLAS_LIB})
else
    echo "Warning: Could not find libopenblas.so"
fi

PATH_REFBLAS_LIB=$(find -name "librefblas.so" | head -n 1 | xargs realpath)
if [ -e "${PATH_REFBLAS_LIB}" ]
then
    echo "Found librefblas.so at ${PATH_REFBLAS_LIB}; adding to LD_LIBRARY_PATH"
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(dirname ${PATH_REFBLAS_LIB})
else
    echo "Warning: Could not find librefblas.so"
fi

PATH_BLIS_LIB=$(find -name "libblis.so" | head -n 1 | xargs realpath)
if [ -e "${PATH_BLIS_LIB}" ]
then
    echo "Found libblis.so at ${PATH_BLIS_LIB}; adding to LD_LIBRARY_PATH"
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(dirname ${PATH_BLIS_LIB})
else
    echo "Warning: Could not find libblis.so"
fi

PATH_OZIMMU_LIB=$(find -name "libozimmu.so" | head -n 1 | xargs realpath)
if [ -e "${PATH_OZIMMU_LIB}" ]
then
    echo "Found libozimmu.so at ${PATH_OZIMMU_LIB}; adding to LD_LIBRARY_PATH"
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(dirname ${PATH_OZIMMU_LIB})
else
    echo "Warning: Could not find libozimmu.so"
fi

PATH_CUBLAS_LIB=$(find ${CUDA_HOME}/../../math_libs -name "libcublas.so" | head -n 1 | xargs realpath)
if [ -e "${PATH_CUBLAS_LIB}" ]
then
    echo "Found libcublas.so at ${PATH_CUBLAS_LIB}; adding to LIBRARY_PATH"
    export LIBRARY_PATH=${LIBRARY_PATH}:$(dirname ${PATH_CUBLAS_LIB})
else
    echo "Warning: Could not find libcublas.so"
fi

set -v
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_NUM_THREADS=64
export OZIMMU_COMPUTE_MODE=fp64_int8_3
set +v