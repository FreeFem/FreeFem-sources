#!/bin/bash

## This script allows to install PETSc on vm2-2
## See ./README.md

echo "update ffPETSc(mpich)"
set -e

# change default  compiler
change_compiler=etc/jenkins/change_compiler/change_compiler-`uname -s`-`uname -r`-4_mpich.sh
test -f "$change_compiler" && echo  source file "$change_compiler"
test -f "$change_compiler" && cat  "$change_compiler"
test -f "$change_compiler" && source "$change_compiler"

if [ "$(uname)" == "Darwin" ]; then
  # in case where the OS type is Darwin
  PETSC_INSTALLDIR='/Users/Shared/mpich/ff-petsc'
elif [ "$(uname)" == "Linux" ]; then
  # in case where the OS type is Linux
PETSC_INSTALLDIR='/builds/Shared/mpich/ff-petsc'
fi

# configuration & build
autoreconf -i \
  && ./configure --enable-download --prefix=$PETSC_INSTALLDIR --enable-bemtool=no \
  && ./3rdparty/getall -a -o PETSc,Ipopt,NLopt,freeYams,FFTW,ARPACK,Gmm++,MMG3D,mshmet,MUMPS,htool \
  && ./etc/jenkins/blob/build_PETSc.sh \
  
if [ $? -eq 0 ]
then
  echo "ffpetsc update complete"
else
  echo "ffpetsc update complete"
  exit 1
fi
