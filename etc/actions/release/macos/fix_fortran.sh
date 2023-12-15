#!/bin/bash

set -e
set -x

PACKAGE_DIR=$1
PREFIX=$2

FORTRAN_FIX="$HOME/fortran-fix"

mkdir "$FORTRAN_FIX"
find "$PACKAGE_DIR/$PREFIX" -name '*.dylib' >"$FORTRAN_FIX/dylib-ff.txt"

./bin/change-dylib-gfortran "$PREFIX/gnu" "$(cat "$FORTRAN_FIX/dylib-ff.txt")"
