#! /bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}/submodules

# clean
rm -rf ozIMMU

# acquire
if ! git clone https://github.com/enp1s0/ozIMMU --recursive; then
    exit 1
fi

cd ozIMMU

# checkout specific commit
git checkout 08eea92

PATH_CUBLAS_LIB=$(find ${CUDA_HOME}/../../math_libs -name "libcublas.so" | head -n 1 | xargs realpath)
if [ -e "${PATH_CUBLAS_LIB}" ]
then
    export LIBRARY_PATH=${LIBRARY_PATH}:$(dirname ${PATH_CUBLAS_LIB})
else
    echo "Warning: Could not find libcublas.so; build may fail"
fi

# build
mkdir build
cd build
cmake ..
make -j4