#!/usr/bin/env sh
set -eux

apt-get update
apt-get install -y libgsl-dev libfftw3-dev
git submodule update --init --recursive
git submodule foreach 'sh install.sh'

set +eux
