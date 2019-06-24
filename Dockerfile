# NOTE: if you are building this on your system and are not pulling from Dockerhub,
# a lot or warnings may be displayed; to turn them off, use "docker build -q [PATH]",
# where [PATH] is the location of this Dockerfile;
# most are harmless and do not affect the outcome of the build
FROM gcc:latest
MAINTAINER Goran Jelic-Cizmek "goran.jelic-cizmek@unige.ch"
# get all the standard required libraries
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    autoconf-archive git \
    libconfig-dev \
    libgsl-dev \
    libfftw3-dev \
    wget && \
    apt-get clean && apt-get autoremove -y
# get CUBA and build it
WORKDIR /tmp/
RUN wget http://www.feynarts.de/cuba/Cuba-4.2.tar.gz
RUN tar xf /tmp/Cuba-4.2.tar.gz -C /tmp/
WORKDIR /tmp/Cuba-4.2/
RUN ./configure
RUN make
RUN make install
RUN rm -rf /tmp/Cuba*
# get CLASS and set it up
WORKDIR /tmp/
# CLASS has a small bug when using the static library,
# see https://github.com/lesgourg/class_public/issues/255 for details,
# which was fixed in the fork
RUN git clone -b libfix --single-branch https://github.com/JCGoran/class_public class
WORKDIR /tmp/class/
RUN make libclass.a
RUN cp libclass.a /usr/lib/x86_64-linux-gnu/
RUN cp include/*.h /usr/include/
RUN rm -rf /tmp/class*
# get COFFE and set it up
WORKDIR /
RUN git clone -q https://github.com/JCGoran/coffe
WORKDIR /coffe/
# this is necessary as long as the future branch isn't merged into master
RUN git checkout future
RUN autoreconf -i
RUN ./configure --enable-cuba --enable-class --enable-doubleexp
RUN make
RUN ln -s /coffe/coffe /usr/bin/coffe
WORKDIR /data/
