#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}

CONFIG_FILES=(
    ./config_files/config_ozimmu.sh.template
    ./config_files/config_gemmul8_accu.sh.template
    ./config_files/config_gemmul8_fast.sh.template
    ./config_files/config_strassen.sh.template
)

for x in "${CONFIG_FILES[@]}"
do
    cp ${x} config.sh

    set -a
    source config.sh 5 5 500 500 100 test3 0 0 0 0 1
    set +a
    
    LIBRARY_PATH=${LIBRARY_PATH}:$(dirname ${SO_PATH}) LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(dirname ${SO_PATH}) make run
done

CONFIG_FILES=(
    ./config_files/config_gemmul8_fast.sh.template
    ./config_files/config_strassen.sh.template
)

for x in "${CONFIG_FILES[@]}"
do
    cp ${x} config.sh

    set -a
    source config.sh 5 5 500 500 100 test3 0 0 1 0 1
    set +a
    
    LIBRARY_PATH=${LIBRARY_PATH}:$(dirname ${SO_PATH}) LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(dirname ${SO_PATH}) make run
done