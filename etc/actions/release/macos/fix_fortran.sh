#!/bin/bash

set -e
set -x

PACKAGE_DIR=$1
PREFIX=$2

FORTRAN_FIX="$HOME/fortran-fix"

mkdir "$FORTRAN_FIX"
find "$PACKAGE_DIR/$PREFIX" -name '*.dylib' >"$FORTRAN_FIX/dylib-ff.txt"

LIB_FORTRAN=$(grep gfortran $PREFIX/gnu/list-dylib-gfortran)
LIB_QUADMATH=$(grep libquadmath $PREFIX/gnu/list-dylib-gfortran)
LIB_FORTRAN_OLD_DIR=$(dirname "$LIB_FORTRAN")
LIB_FORTRAN_DYLIB=$(basename "$LIB_FORTRAN")
LIB_QUADMATH_DYLIB=$(basename "$LIB_QUADMATH")

./bin/change-dylib-gfortran "$LIB_FORTRAN_OLD_DIR" "$PREFIX/gnu" "$LIB_FORTRAN_DYLIB" "$LIB_QUADMATH_DYLIB" "$(cat "$FORTRAN_FIX/dylib-ff.txt")"
