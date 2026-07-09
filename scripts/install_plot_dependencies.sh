#!/bin/bash

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
cd ${REPO_ROOT}

python3 -m venv .venv
source .venv/bin/activate
pip3 install "matplotlib>=3.5"