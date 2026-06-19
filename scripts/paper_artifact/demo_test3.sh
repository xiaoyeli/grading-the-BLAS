#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
cd ${REPO_ROOT}

TEST_CONFIGS=(
    ./paper_artifact/config_files/test_configs/demo_3.mk
)

DGEMM_CONFIGS=(
    ./paper_artifact/config_files/dgemm_configs/ozaki1.mk
    ./paper_artifact/config_files/dgemm_configs/ozaki2_accu_mode.mk
    ./paper_artifact/config_files/dgemm_configs/ozaki2_fast_mode.mk
    ./paper_artifact/config_files/dgemm_configs/ozaki2_fast_mode_w_outer_scaling.mk
    ./paper_artifact/config_files/dgemm_configs/strassen.mk
    ./paper_artifact/config_files/dgemm_configs/strassen_w_outer_scaling.mk
)

for t in "${TEST_CONFIGS[@]}"
do
    for d in "${DGEMM_CONFIGS[@]}"
    do
        make run DGEMM_CONFIG_MK=${d} TEST_CONFIG_MK=${t}
    done
done