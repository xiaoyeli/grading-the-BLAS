#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}

DGEMM_IMPLEMENTATIONS=(
    NAIVE_TRIPLE_LOOP_EXTENDED_PRECISION
    NAIVE_TRIPLE_LOOP_FP64
    OZIMMU
    GEMMUL8
    STRASSEN
    STRASSEN_DABD_SCALING
    STRASSEN_DABD_SCALING_AND_ZEROING
)

for impl in "${DGEMM_IMPLEMENTATIONS[@]}"
do
    echo "================== Testing ${impl}"
    make run MIN_DIM=32 MAX_DIM=32 DGEMM_IMPL=${impl} INPUT_TYPE=test1a
    make run MIN_DIM=32 MAX_DIM=32 DGEMM_IMPL=${impl} INPUT_TYPE=test1b
    make run MIN_DIM=32 MAX_DIM=32 SCALING_FACTOR=500 DGEMM_IMPL=${impl} INPUT_TYPE=test1c
done