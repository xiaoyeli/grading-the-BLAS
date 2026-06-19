#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)
cd ${REPO_ROOT}

TEST_CONFIGS=(
    ./paper_artifact/config_files/test_configs/demo_4a.mk
    ./paper_artifact/config_files/test_configs/demo_4b.mk
)

DGEMM_CONFIGS=(
    ./paper_artifact/config_files/dgemm_configs/naive_triple_loop.mk
    ./paper_artifact/config_files/dgemm_configs/ozaki1.mk
    ./paper_artifact/config_files/dgemm_configs/ozaki2_accu_mode.mk
    ./paper_artifact/config_files/dgemm_configs/ozaki2_fast_mode.mk
    ./paper_artifact/config_files/dgemm_configs/strassen.mk
)

for t in "${TEST_CONFIGS[@]}"
do
    for d in "${DGEMM_CONFIGS[@]}"
    do
        make run DGEMM_CONFIG_MK=${d} TEST_CONFIG_MK=${t}
    done
done