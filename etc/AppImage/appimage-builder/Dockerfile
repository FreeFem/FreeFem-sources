FROM ubuntu:16.04

MAINTAINER Alexander Sashnov "sashnov@ngs.ru"

ENV DEBIAN_FRONTEND noninteractive

ENV HOME    /home/ubuntu

# WORKDIR directive will create the directory automatically
WORKDIR /home/ubuntu

RUN apt update && apt full-upgrade -y && apt install -y --no-install-recommends \
    autoconf \
    automake \
    autotools-dev \
    bison \
    ca-certificates \
    cmake \
    coinor-libipopt-dev \
    file \
    flex \
    freeglut3-dev \
    g++ \
    gcc \
    gdb \
    gfortran \
    ghostscript \
    git \
    gnuplot-qt \
    libgsl0-dev \
    libarpack2-dev \
    libfftw3-dev \
    libgmm++-dev \
    libhdf5-dev \
    liblapack-dev \
    libmumps-seq-dev \
    libnlopt-dev \
    libopenblas-dev \
    libscotch-dev \
    libsuitesparse-dev \
    libtet1.5-dev \
    locales \
    m4 \
    make \
    mpich \
    patch \
    pkg-config \
    python \
    sudo \
    unzip \
    valgrind \
    wget \
  && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# Reconfigure locale
RUN locale-gen en_US.UTF-8 && dpkg-reconfigure locales

# Add group & user + sudo
RUN groupadd -r ubuntu && useradd --create-home --gid ubuntu ubuntu && echo 'ubuntu ALL=NOPASSWD: ALL' > /etc/sudoers.d/ubuntu
