#!/bin/bash

set -e
set -x

sudo ln -fs /usr/local/bin/gfortran-12 /usr/local/bin/gfortran
sudo ln -fs /usr/local/bin/gcc-12 /usr/local/bin/gcc
sudo ln -fs /usr/local/bin/g++-12 /usr/local/bin/g++

# symlink dylib location for previous versions
sudo ln -fs /usr/local/opt/gcc/lib/gcc/12 /usr/local/opt/gcc/lib/gcc/11
sudo ln -fs /usr/local/opt/gcc/lib/gcc/12 /usr/local/opt/gcc/lib/gcc/10
