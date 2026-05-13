#!/bin/bash

set -a
source config.sh 
set +a

LIBRARY_PATH=${LIBRARY_PATH}:$(dirname ${SO_PATH}) LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$(dirname ${SO_PATH}) make -f Makefile run