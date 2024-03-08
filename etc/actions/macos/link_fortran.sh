#!/bin/bash

set -e
set -x

ARCH=$1

LOCAL_PATH=/usr/local
if [ "$ARCH" = "ARM64" ]; then
	LOCAL_PATH=/opt/homebrew
fi

sudo ln -fs "$LOCAL_PATH/bin/gfortran-12" "$LOCAL_PATH/bin/gfortran"
sudo ln -fs "$LOCAL_PATH/bin/gcc-12" "$LOCAL_PATH/bin/gcc"
sudo ln -fs "$LOCAL_PATH/bin/g++-12" "$LOCAL_PATH/bin/g++"

# symlink dylib location for previous versions
sudo ln -fs "$LOCAL_PATH/opt/gcc/lib/gcc/12" "$LOCAL_PATH/opt/gcc/lib/gcc/11"
sudo ln -fs "$LOCAL_PATH/opt/gcc/lib/gcc/12" "$LOCAL_PATH/opt/gcc/lib/gcc/10"
