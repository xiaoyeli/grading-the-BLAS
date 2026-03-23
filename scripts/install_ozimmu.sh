#! /bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
VERSION_HASH=08eea92

cd ${REPO_ROOT}/submodules

# clean
rm -rf ozIMMU-${VERSION_HASH}

# acquire
if ! git clone https://github.com/enp1s0/ozIMMU --recursive; then
    exit 1
fi

mv ozIMMU ozIMMU-${VERSION_HASH}
cd ozIMMU-${VERSION_HASH}

# checkout specific commit
git checkout ${VERSION_HASH}

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