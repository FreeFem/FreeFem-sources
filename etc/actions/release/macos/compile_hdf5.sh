#!/bin/bash

set -e
set -x

PREFIX=$1

sudo mkdir -p "$PREFIX/pkg"
cd "$PREFIX/pkg"

wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.14/hdf5-1.14.3/src/hdf5-1.14.3.tar.bz2
tar zxf hdf5-1.14.3.tar.bz2
cd hdf5-1.14.3
./configure --enable-cxx --prefix="$PREFIX" CC="$CC" CXX="$CXX"
make -j2
make install
