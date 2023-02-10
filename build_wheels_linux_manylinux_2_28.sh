#!/usr/bin/env sh
set -eux

yum install -y gsl-devel fftw-devel
git submodule update --init --recursive
git submodule foreach 'sh install.sh'

set +eux
