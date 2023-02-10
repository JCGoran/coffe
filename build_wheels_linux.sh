#!/usr/bin/env sh

set -eux
# script for building multiple manylinux images

# by default we build the manylinux2014 wheel
CIBW_BEFORE_ALL="sh build_wheels_linux_${1:-manylinux2014}.sh" \
    CIBW_MANYLINUX_X86_64_IMAGE="${1:-manylinux2014}" \
    cibuildwheel --platform linux

set +eux
