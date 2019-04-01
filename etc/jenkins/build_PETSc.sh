#!/bin/sh

# if sudo is needed
# lanuch ./build_PETSc.sh sudo

cd 3rdparty/ff-petsc \
	&& make petsc-slepc SUDO=$1 \
	&& cd - \
	&& ./reconfigure
