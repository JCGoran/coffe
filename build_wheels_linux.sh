#!/usr/bin/env sh

set -euxo
# script for building multiple manylinux images

# 1. manylinux_2_24
#CIBW_BEFORE_ALL="apt-get update && apt-get install -y libgsl-dev libfftw3-dev && git submodule update --init --recursive && git submodule foreach './install.sh'" \
#    CIBW_MANYLINUX_X86_64_IMAGE="manylinux_2_24" \
#    cibuildwheel --platform linux


# 2. manylinux2014
# we need to build GSL from source as the one in the repos is too old
GSL_VERSION="2.0"

CIBW_BEFORE_ALL="curl -sL -vvv https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz --output libgsl.tar.gz && tar xf libgsl.tar.gz && (cd gsl-${GSL_VERSION} && ./autogen.sh && ./configure && make && make install) && yum install -y fftw-devel && git submodule update --init --recursive && git submodule foreach './install.sh'" \
    CIBW_MANYLINUX_X86_64_IMAGE="manylinux2014" \
    cibuildwheel --platform linux


# 2. manylinux_2_28
#CIBW_BEFORE_ALL="yum install -y gsl-devel fftw-devel && git submodule update --init --recursive && git submodule foreach './install.sh'" \
#    CIBW_MANYLINUX_X86_64_IMAGE="manylinux_2_28" \
#    cibuildwheel --platform linux

set +euxo
