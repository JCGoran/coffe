#!/usr/bin/env sh

# THIS IS THE ONLY SCRIPT YOU NEED TO RUN TO INSTALL COFFE

set -eux

# in case of no conda env, we just install the Python-specific requirements
if [ -z "${CONDA_PREFIX-}" ] && [ -z "${CONDA_DEFAULT_ENV-}" ]
then
    python3 -m pip install Cython wheel
# otherwise, we need to install all of the requirements (compiler, GSL lib, etc.)
else
    conda install --channel conda-forge --file requirements.txt
fi

git submodule update --init --recursive
git submodule foreach './install.sh'
sh reinstall.sh

set +eux
