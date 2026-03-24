#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}

DGEMM_IMPLEMENTATIONS=(
    NAIVE_TRIPLE_LOOP_EXTENDED_PRECISION
    OZIMMU
    OZIMMU_ADDB_SCALING
    GEMMUL8
    GEMMUL8_ADDB_SCALING
    STRASSEN
    STRASSEN_ADDB_SCALING
)

for impl in "${DGEMM_IMPLEMENTATIONS[@]}"
do
    echo "================== Testing ${impl}"
    make run MIN_DIM=64 MAX_DIM=64 DGEMM_IMPL=${impl} INPUT_TYPE=test2a
    make run MIN_DIM=64 MAX_DIM=64 DGEMM_IMPL=${impl} INPUT_TYPE=test2b
    make run MIN_DIM=64 MAX_DIM=64 DGEMM_IMPL=${impl} INPUT_TYPE=test2c
done