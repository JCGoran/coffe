ARG MANYLINUX_IMAGE=manylinux_2_28_x86_64
FROM quay.io/pypa/$MANYLINUX_IMAGE
LABEL author="Goran Jelic-Cizmek"

COPY install_other.sh .
RUN bash install_other.sh cuba gsl class fftw libconfig
