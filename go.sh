#!/bin/bash

set -a
source config.sh
set +a

export LIBRARY_PATH=${LIBRARY_PATH}:$(dirname ${SO_PATH})
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(dirname ${SO_PATH})

make -f Makefile.for_distribution run