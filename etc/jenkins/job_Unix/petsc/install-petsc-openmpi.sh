#!/usr/bin/env bash

if [ "$(uname)" = "Darwin" ]
then
    # TODO
    :
elif [ "$(uname)" = "Linux" ]
then
    PETSC_INSTALL_DIR='/media/builds/shared/openmpi'

    export MPIRUN=/usr/bin/mpirun.openmpi
    export MPICXX=/usr/bin/mpicxx.openmpi
    export MPIFC=/usr/bin/mpif90.openmpi
    export MPICC=/usr/bin/mpicc.openmpi
fi

PETScVersion=$(grep "VERSION=" 3rdparty/ff-petsc/Makefile | cut -c9-15)

# configuration & build
autoreconf -i
./configure --enable-download --prefix="$PETSC_INSTALL_DIR"
./3rdparty/getall -a -o PETSc
./etc/jenkins/blob/build_PETSc.sh
rm -rf "3rdparty/ff-petsc/petsc-$PETScVersion"

if [ $? -eq 0 ]
then
  echo "PETSc update complete"
else
  echo "PETSc update failed"
  rm -rf "$PETSC_INSTALL_DIR/ff-petsc"
fi
