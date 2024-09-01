#!/usr/bin/env sh

set -ex

CLASS_INSTALL_DIR="${CLASS_INSTALL_DIR:-/opt/class_public_$(uname -m)}"
CUBA_INSTALL_DIR="${CUBA_INSTALL_DIR:-/opt/cuba_$(uname -m)}"
# the abs dir where this script is located (so we can call it from wherever)
script_dir="$(cd "$(dirname "$0")"; pwd -P)"

install_cuba(){
    if [ -e "${CUBA_INSTALL_DIR}/lib/libcuba.a" ]
    then
        printf "CUBA already installed at path %s\n" "${CUBA_INSTALL_DIR}"
        return 0
    fi
    cd "${script_dir}/../external/libcuba/"
    export MACOSX_DEPLOYMENT_TARGET='11.0'
    git clean -xdf .
    if [ "$(uname -s)" = 'Darwin' ]
    then
        CFLAGS="-fPIC -mmacosx-version-min=${MACOSX_DEPLOYMENT_TARGET}"
    else

        CFLAGS="-fPIC"
    fi
    autoreconf --install
    ./configure --prefix="${CUBA_INSTALL_DIR}" CFLAGS="${CFLAGS}"

    make install
    cd -
    printf 'CUBA installed in %s\n' "${CUBA_INSTALL_DIR}"
}



install_class(){
    if [ -e "${CLASS_INSTALL_DIR}/lib/libclass.a" ]
    then
        printf "CLASS already installed at path %s\n" "${CLASS_INSTALL_DIR}"
        return 0
    fi
    current_dir="${script_dir}/../external/class_public/"
    cd "${current_dir}"
    export MACOSX_DEPLOYMENT_TARGET='11.0'
    git clean -xdf .
    make libclass.a

    mkdir -p "${CLASS_INSTALL_DIR}/lib" "${CLASS_INSTALL_DIR}/include"
    cp -a "${current_dir}/libclass.a" "${CLASS_INSTALL_DIR}/lib/"
    cp -a "${current_dir}/include/"*.h "${CLASS_INSTALL_DIR}/include/"
    cd -
    printf 'CLASS installed in %s\n' "${CLASS_INSTALL_DIR}"
}

install_fftw(){
    if [ "$(uname -s)" = 'Darwin' ]
    then
        printf "Please use Conan to install FFTW\n"
        return 1
    fi
    yum install -y fftw-devel
}

install_gsl(){
    if [ "$(uname -s)" = 'Darwin' ]
    then
        printf "Please use Conan to install GSL\n"
        return 1
    fi
    yum install -y gsl-devel
}

install_libconfig(){
    if [ "$(uname -s)" = 'Darwin' ]
    then
        printf "Please use Conan to install libconfig\n"
        return 1
    fi
    yum install -y libconfig-devel
}

for arg in $@
do
    case "${arg}" in
        "gsl")
            install_gsl
            ;;
        "cuba")
            install_cuba
            ;;
        "class")
            install_class
            ;;
        "fftw")
            install_fftw
            ;;
        "libconfig")
            install_libconfig
            ;;
        *)
            ;;
    esac
done

set +ex
