#!/usr/bin/env bash

echo "Job 4 (openmpi)"

casejob=4_openmpi

# change default compiler
if [ "$(uname)" = "Darwin" ]
then
  # MacOS
  PETSC_DIR='/Users/Shared/openmpi/ff-petsc'
  change_compiler=etc/jenkins/change_compiler/change_compiler-$(uname -s)-$(uname -r)-$casejob.sh
elif [ "$(uname)" = "Linux" ]
then
  # Linux
  PETSC_DIR='/media/builds/shared/openmpi/ff-petsc'

  export MPIRUN=/usr/bin/mpirun.openmpi
  export MPICXX=/usr/bin/mpicxx.openmpi
  export MPIFC=/usr/bin/mpif90.openmpi
  export MPICC=/usr/bin/mpicc.openmpi
fi

petscVersion=$(grep "VERSION_GIT" "$PETSC_DIR/c/include/petscversion.h" | cut -c36-42 | cut -f1 -d "\"" )
echo "Installed PETSc version: $petscVersion"
echo "Try to source file $change_compiler"
test -f "$change_compiler" && echo "Source file $change_compiler"
test -f "$change_compiler" && cat "$change_compiler"
test -f "$change_compiler" && source "$change_compiler"

# configuration & build
autoreconf -i
./configure --enable-download --prefix="$WORKSPACE/install" \
  --with-petsc=$PETSC_DIR/r/lib \
  --with-petsc_complex=$PETSC_DIR/c/lib
./3rdparty/getall -a -o Ipopt,NLopt,freeYams,FFTW,Gmm++,MMG3D,mshmet,MUMPS
./etc/jenkins/blob/build.sh

if [ $? -eq 0 ]
then
  echo "Build process complete"
else
  echo "Build process failed"
  exit 1
fi

# check
./etc/jenkins/blob/check.sh

if [ $? -eq 0 ]
then
  echo "Check process complete"
else
  echo "Check process failed"
fi

# install
./etc/jenkins/blob/install.sh

if [ $? -eq 0 ]
then
  echo "Install process complete"
else
  echo "Install process failed"
  exit 1
fi

# uninstall
./etc/jenkins/blob/uninstall.sh

if [ $? -eq 0 ]
then
  echo "Uninstall process complete"
else
  echo "Uninstall process failed"
  exit 1
fi

# Jenkins tests results analyser
./etc/jenkins/resultForJenkins/resultForJenkins.sh

if [ $? -eq 0 ]
then
  echo "Jenkins process complete"
else
  echo "Jenkins process failed"
fi