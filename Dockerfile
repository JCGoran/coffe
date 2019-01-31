FROM gcc:latest
MAINTAINER Goran Jelic-Cizmek "goran.jelic-cizmek@unige.ch"
# get all the standard required libraries
RUN apt-get -qq update && apt-get install -yqq autoconf-archive git libconfig-dev libgsl-dev libfftw3-dev wget
WORKDIR /tmp/
# get CUBA and build it
RUN wget -q http://www.feynarts.de/cuba/Cuba-4.2.tar.gz
RUN tar xzf /tmp/Cuba-4.2.tar.gz -C /tmp/
WORKDIR /tmp/Cuba-4.2/
RUN ./configure > /dev/null 2>&1
RUN make > /dev/null 2>&1
RUN make install > /dev/null 2>&1
RUN rm -rf /tmp/Cuba*
WORKDIR /
# get COFFE and set it up
RUN git clone -q https://github.com/JCGoran/coffe
WORKDIR /coffe/
RUN autoreconf -i > /dev/null 2>&1
RUN ./configure --enable-cuba > /dev/null 2>&1
RUN make > /dev/null 2>&1
RUN ln -s /coffe/coffe /usr/bin/coffe
