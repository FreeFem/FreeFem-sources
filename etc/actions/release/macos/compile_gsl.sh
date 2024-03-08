#!/bin/bash

set -e
set -x

PREFIX=$1
CC=$2
CXX=3

sudo mkdir -p "$PREFIX/pkg"
cd "$PREFIX/pkg"

wget https://ftp.gnu.org/gnu/gsl/gsl-2.7.1.tar.gz
tar zxf gsl-2.7.1.tar.gz
cd gsl-2.7.1
./configure --prefix="$PREFIX" CC="$CC" CXX="$CXX"
make -j2
make install
