#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}

set -v
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(realpath ./submodules/blis-2.0/lib/sandybridge/)
export OMP_NUM_THREADS=4
set +v