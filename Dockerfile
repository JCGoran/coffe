# NOTE: if you are building this on your system and are not pulling from Dockerhub,
# a lot or warnings may be displayed; to turn them off, use "docker build -q [PATH]",
# where [PATH] is the location of this Dockerfile;
# most are harmless and do not affect the outcome of the build
FROM gcc:8
MAINTAINER Goran Jelic-Cizmek "goran.jelic-cizmek@unige.ch"
# get all the standard required libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    autoconf-archive git \
    libconfig-dev \
    libgsl-dev \
    libfftw3-dev \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    wget && \
    apt-get clean && apt-get autoremove -y
# get CUBA and build it
WORKDIR /tmp/
RUN wget http://www.feynarts.de/cuba/Cuba-4.2.tar.gz
RUN tar xf /tmp/Cuba-4.2.tar.gz -C /tmp/
WORKDIR /tmp/Cuba-4.2/
RUN ./configure
RUN make -j lib
RUN make install
RUN rm -rf /tmp/Cuba*
# get CLASS and set it up
WORKDIR /tmp/
# CLASS has a small bug when using the static library,
# see https://github.com/lesgourg/class_public/issues/255 for details,
# which was fixed in the fork
RUN git clone -b master --single-branch https://github.com/lesgourg/class_public class
WORKDIR /tmp/class/
RUN make -j libclass.a
RUN cp libclass.a /usr/lib/x86_64-linux-gnu/
RUN cp include/*.h /usr/include/
# get COFFE and set it up
WORKDIR /
RUN git clone -q https://github.com/JCGoran/coffe
WORKDIR /coffe/
# this is necessary as long as the future branch isn't merged into master
RUN git checkout future
RUN autoreconf -i
RUN ./configure --enable-cuba --enable-class
RUN make -j coffe
# optional dependencies for making life easier
RUN pip3 install \
    wheel \
    numpy \
    scipy \
    matplotlib \
    jupyter
# installind dependencies for COFFE python "wrapper"
RUN pip3 install -r python/requirements.txt
RUN ln -s /coffe/coffe /usr/bin/coffe
