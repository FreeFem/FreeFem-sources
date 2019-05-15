#!/bin/bash

# if sudo is needed
# launch ./build_PETSc.sh sudo

WITH_SUDO=""
if [ "$1" = "sudo" ]
then
	WITH_SUDO=sudo
fi

cd 3rdparty/ff-petsc \
	&& make petsc-slepc SUDO=${WITH_SUDO} \
	&& cd - \
	&& ./reconfigure
