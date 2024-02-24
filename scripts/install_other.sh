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
    sudo mkdir -p "${CLASS_INSTALL_DIR}/lib" "${CLASS_INSTALL_DIR}/include"
    sudo cp -a "${class_dir}/libclass.a" "${CLASS_INSTALL_DIR}/lib/"
    sudo cp -a "${class_dir}/include/"*.h "${CLASS_INSTALL_DIR}/include/"
    cd -
    printf 'CLASS installed\n'
}

install_fftw(){
    yum install -y fftw-devel
}

install_gsl(){
    GSL_VERSION='2.0'
    curl -sL -vvv "https://ftp.gnu.org/gnu/gsl/gsl-${GSL_VERSION}.tar.gz" --output libgsl.tar.gz
    tar xf libgsl.tar.gz
    cd "gsl-${GSL_VERSION}"
    ./autogen.sh
    ./configure
    make -j 2 all
    make install
    cd -
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
        *)
            ;;
    esac
done

set +ex
