#!/usr/bin/env bash

echo "Job PETsc"

freefemPETScVersion=$(grep "VERSION=" 3rdparty/ff-petsc/Makefile | cut -c9-15)

#################
# OpenMPI version
if [ "$(uname)" == "Darwin" ]
then
    # Macos
    openmpiPETScDir=''
    ## TODO
elif [ "$(uname)" == "Linux" ]
then
    # Linux
    openmpiPETScDir='/media/builds/shared/openmpi'
fi

major=$(grep "#define PETSC_VERSION_MAJOR" "$openmpiPETScDir/ff-petsc/c/include/petscversion.h" | cut -c34-35 )
minor=$(grep "#define PETSC_VERSION_MINOR" "$openmpiPETScDir/ff-petsc/c/include/petscversion.h" | cut -c34-36 )
subminor=$(grep "#define PETSC_VERSION_SUBMINOR" "$openmpiPETScDir/ff-petsc/c/include/petscversion.h" | cut -c34-35 )
openmpiPETScVersion=$major.$minor.$subminor

echo "Installer OpenMPI PETSc version: $openmpiPETScVersion"

# Check version
if [ "$freefemPETScVersion" = "$openmpiPETScVersion" ]
then
    echo "OpenMPI PETSc is up-to-date"
    echo "Nothing to do"
else
    echo "OpenMPI PETSc is outdated"
    echo "Update"
    rm -rf $openmpiPETScDir/ff-petsc
    ./etc/jenkins/job_Unix/petsc/install-petsc-openmpi.sh
fi

###############
# MPICH version
if [ "$(uname)" == "Darwin" ]
then
    # Macos
    mpichPETScDir=''
    ## TODO
elif [ "$(uname)" == "Linux" ]
then
    # Linux
    mpichPETScDir='/media/builds/shared/mpich'
fi

major=$(grep "#define PETSC_VERSION_MAJOR" "$mpichPETScDir/ff-petsc/c/include/petscversion.h" | cut -c34-35 )
minor=$(grep "#define PETSC_VERSION_MINOR" "$mpichPETScDir/ff-petsc/c/include/petscversion.h" | cut -c34-36 )
subminor=$(grep "#define PETSC_VERSION_SUBMINOR" "$mpichPETScDir/ff-petsc/c/include/petscversion.h" | cut -c34-35 )
mpichPETScVersion=$major.$minor.$subminor

echo "Installer MPICH PETSc version: $mpichPETScVersion"

# Check version
if [ "$freefemPETScVersion" = "$mpichPETScVersion" ]
then
    echo "MPICH PETSc is up-to-date"
    echo "Nothing to do"
else
    echo "MPICH PETSc is outdated"
    echo "Update"
    rm -rf $mpichPETScDir/ff-petsc
    ./etc/jenkins/job_Unix/petsc/install-petsc-mpich.sh
fi