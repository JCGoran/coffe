#!/usr/bin/env sh
set -eux

# we need to build GSL from source as the one in the repos is too old (need
# version 2.1 for 2D interpolation)
GSL_VERSION="2.0"

curl -sL -vvv https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz --output libgsl.tar.gz
tar xf libgsl.tar.gz
(cd gsl-${GSL_VERSION} && ./autogen.sh && ./configure && make && make install)
yum install -y fftw-devel
git submodule update --init --recursive
git submodule foreach 'sh install.sh'

set +eux
