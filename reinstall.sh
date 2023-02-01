#!/usr/bin/env sh

# This script is used to rebuild COFFE is changes are made to the source files

set -eux

python3 setup.py build_ext -i
python3 -m pip install .

set +eux
