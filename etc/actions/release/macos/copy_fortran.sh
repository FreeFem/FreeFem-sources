#!/bin/bash

set -e
set -x

PREFIX=$1

# Create path
mkdir -p "$PREFIX/gnu"

# Fortran check (get libraries)
FORTRAN_CHECK="$HOME/fortran-check"
mkdir "$FORTRAN_CHECK"
printf "program sizeofint\n  integer i\nend" >"$FORTRAN_CHECK/check.f90"
gfortran "$FORTRAN_CHECK/check.f90" -o "$FORTRAN_CHECK/check"
DYLD_PRINT_LIBRARIES=1 "$FORTRAN_CHECK/check" 2>&1 | awk '{print $NF}' | egrep -v '/usr/lib | "$FORTRAN_CHECK/check" | grep .dylib' >$PREFIX/gnu/list-dylib-gfortran
rm -rf "$FORTRAN_CHECK"

# Copy librairies
BREW_LIB_GFORTRAN=$(grep gfortran "$PREFIX/gnu/list-dylib-gfortran")

cd $(dirname "$BREW_LIB_GFORTRAN")
cp -f libgfortran.*.dylib libquadmath.*.dylib libgcc_s.*.dylib "$PREFIX/gnu/"

cd "$PREFIX/gnu"
ln -sf libgfortran.?.dylib libgfortran.dylib
ln -sf libquadmath.?.dylib libquadmath.dylib
