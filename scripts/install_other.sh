#!/usr/bin/env sh

set -ex

CLASS_INSTALL_DIR="${CLASS_INSTALL_DIR:-/opt/class_public}"
CUBA_INSTALL_DIR="${CUBA_INSTALL_DIR:-/opt/cuba}"
# the abs dir where this script is located (so we can call it from wherever)
script_dir="$(cd "$(dirname "$0")"; pwd -P)"

install_cuba(){
    cd "${script_dir}/../external/libcuba/"
    autoreconf --install
    ./configure --prefix="${CUBA_INSTALL_DIR}" CFLAGS=-fPIC
    if [ -x "${CUBA_INSTALL_DIR}" ]
    then
        make install
    else
        printf 'Root access required to install CUBA to %s, please type in your password at the prompt below\n' "${CUBA_INSTALL_DIR}"
        sudo make install
    fi
    cd -
    printf 'CUBA installed in %s\n' "${CUBA_INSTALL_DIR}"
}



install_class(){
    current_dir="${script_dir}/../external/class_public/"
    cd "${current_dir}"
    make libclass.a

    if [ -x "${CLASS_INSTALL_DIR}" ]
    then
        mkdir -p "${CLASS_INSTALL_DIR}/lib" "${CLASS_INSTALL_DIR}/include"
        cp -a "${current_dir}/libclass.a" "${CLASS_INSTALL_DIR}/lib/"
        cp -a "${current_dir}/include/"*.h "${CLASS_INSTALL_DIR}/include/"
    else
        printf 'Root access required to install CLASS to %s, please type in your password at the prompt below\n' "${CLASS_INSTALL_DIR}"
        sudo mkdir -p "${CLASS_INSTALL_DIR}/lib" "${CLASS_INSTALL_DIR}/include"
        sudo cp -a "${current_dir}/libclass.a" "${CLASS_INSTALL_DIR}/lib/"
        sudo cp -a "${current_dir}/include/"*.h "${CLASS_INSTALL_DIR}/include/"
    fi
    cd -
    printf 'CLASS installed in %s\n' "${CLASS_INSTALL_DIR}"
}

install_fftw(){
    yum install -y fftw-devel
}

install_gsl(){
    yum install -y gsl-devel
}

install_libconfig(){
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
