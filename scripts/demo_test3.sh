#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}

DGEMM_IMPLEMENTATIONS=(
    NAIVE_TRIPLE_LOOP_EXTENDED_PRECISION
    OZIMMU
    GEMMUL8
    GEMMUL8_DABD_SCALING
    STRASSEN
    STRASSEN_DABD_SCALING
)

for impl in "${DGEMM_IMPLEMENTATIONS[@]}"
do
    echo "================== Testing ${impl}"
    make run MIN_DIM=64 MAX_DIM=64 SCALING_FACTOR=500 DGEMM_IMPL=${impl} INPUT_TYPE=test3
done