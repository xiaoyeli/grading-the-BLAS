#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}

set -v
PATH_OPENBLAS_LIB=$(find -name "libopenblas.so" | head -n 1)
if [ -n "${PATH_OPENBLAS_LIB}" ]
then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(realpath $(dirname ${PATH_OPENBLAS_LIB}))
fi

PATH_REFBLAS_LIB=$(find -name "librefblas.so" | head -n 1)
if [ -n "${PATH_REFBLAS_LIB}" ]
then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(realpath $(dirname ${PATH_REFBLAS_LIB}))
fi

PATH_BLIS_LIB=$(find -name "libblis.so" | head -n 1)
if [ -n "${PATH_BLIS_LIB}" ]
then
    export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(realpath $(dirname ${PATH_BLIS_LIB}))
fi

export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_NUM_THREADS=64
set +v