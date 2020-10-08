#!/bin/bash

## This script allows to install PETSc on vm2-2
## See ./README.md

echo "update ffPETSc(openmpi)"
set -e

# change default  compiler
if [ "$(uname)" == "Darwin" ]; then
  # in case where the OS type is Darwin
  PETSC_INSTALLDIR='/Users/Shared/openmpi'
  change_compiler=etc/jenkins/change_compiler/change_compiler-`uname -s`-`uname -r`-4_openmpi.sh
elif [ "$(uname)" == "Linux" ]; then
  # in case where the OS type is Linux
  PETSC_INSTALLDIR='/builds/Shared/openmpi'
  change_compiler=etc/jenkins/change_compiler/change_compiler-`uname -s`-4_openmpi.sh
fi
echo try to source file  "$change_compiler"
test -f "$change_compiler" && echo  source file "$change_compiler"
test -f "$change_compiler" && cat  "$change_compiler"
test -f "$change_compiler" && source "$change_compiler"

# configuration & build
autoreconf -i \
  && ./configure --enable-download --prefix=$PETSC_INSTALLDIR \
  && ./3rdparty/getall -a -o PETSc,Ipopt,NLopt,freeYams,FFTW,Gmm++,MMG3D,mshmet,MUMPS,htool \
  && ./etc/jenkins/blob/build_PETSc.sh
  
if [ $? -eq 0 ]
then
  echo "ffpetsc update complete"
else
  echo "ffpetsc update fail"
  exit 1
fi
