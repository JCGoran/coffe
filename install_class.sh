#!/usr/bin/env sh

set -e
# installer for CLASS when using a conda environment
CLASS_DIR="class_public"
CLASS_REMOTE_URL="https://github.com/JCGoran/class_public"
CLASS_BRANCH="feature/conda"

install_class(){
    if [ -z "${CONDA_PREFIX}" ]
    then
        printf 'You need to activate a conda environment using `conda activate [ENVIRONMENT]` before running this script\n'
        return 1
    fi

    if [ ! -d "${CLASS_DIR}" ]
    then
        printf 'Attempting to install CLASS in the current environment (%s)...\n' "${CONDA_DEFAULT_ENV}"
        printf 'Cloning to directory %s...\n' "${CLASS_DIR}"
        git clone --branch "${CLASS_BRANCH}" "${CLASS_REMOTE_URL}" "${CLASS_DIR}"
    fi

    make -C "${CLASS_DIR}" libclass.a && cp -a "${CLASS_DIR}/libclass.a" "${CONDA_PREFIX}/lib/" && cp -a "${CLASS_DIR}/include/"*.h "${CONDA_PREFIX}/include/"
    printf 'CLASS successfully installed\n'
}

install_class
set +e
