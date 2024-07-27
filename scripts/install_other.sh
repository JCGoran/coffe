#!/usr/bin/env sh

set -ex

CLASS_INSTALL_DIR="/opt/class_public"
CLASS_REMOTE_URL="https://github.com/JCGoran/class_public"
CLASS_BRANCH="feature/conda"

CUBA_INSTALL_DIR="/opt/cuba"
CUBA_REMOTE_URL="https://github.com/JCGoran/libcuba"
CUBA_BRANCH="master"

install_cuba(){
    cuba_dir="$(mktemp -d)"
    git clone --depth 1 --branch "${CUBA_BRANCH}" "${CUBA_REMOTE_URL}" "${cuba_dir}"

    cd "${cuba_dir}"
    autoreconf --install
    ./configure --prefix="${CUBA_INSTALL_DIR}" CFLAGS=-fPIC
    make install
    cd -
    printf 'CUBA installed\n'
}



install_class(){
    class_dir="$(mktemp -d)"
    git clone --depth 1 --branch "${CLASS_BRANCH}" "${CLASS_REMOTE_URL}" "${class_dir}"

    cd "${class_dir}"
    make libclass.a

    if [ "$(uname)" = 'Darwin' ]
    then
        sudo mkdir -p "${CLASS_INSTALL_DIR}/lib" "${CLASS_INSTALL_DIR}/include"
        sudo cp -a "${class_dir}/libclass.a" "${CLASS_INSTALL_DIR}/lib/"
        sudo cp -a "${class_dir}/include/"*.h "${CLASS_INSTALL_DIR}/include/"

    else
        mkdir -p "${CLASS_INSTALL_DIR}/lib" "${CLASS_INSTALL_DIR}/include"
        cp -a "${class_dir}/libclass.a" "${CLASS_INSTALL_DIR}/lib/"
        cp -a "${class_dir}/include/"*.h "${CLASS_INSTALL_DIR}/include/"
    fi

    cd -
    printf 'CLASS installed\n'
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
