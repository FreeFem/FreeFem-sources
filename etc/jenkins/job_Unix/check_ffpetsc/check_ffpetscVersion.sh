#!/bin/bash

releaseVersionffpetsc=$(grep "VERSION=" 3rdparty/ff-petsc/Makefile | cut -c9-15)
echo "release Version ffpetsc" $releaseVersionffpetsc

# compilation with openmpi
if [ "$(uname)" == "Darwin" ]; then
  # in case where the OS type is Darwin
  MAJOR=$(grep "#define PETSC_VERSION_MAJOR" /Users/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c34-35 )
  MINOR=$(grep "#define PETSC_VERSION_MINOR" /Users/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c34-36 )
  SUBMINOR=$(grep "#define PETSC_VERSION_SUBMINOR" /Users/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c34-35 )
  installedVersionffpetscO=$MAJOR.$MINOR.$SUBMINOR
  ffpetscDirectory=/Users/Shared/openmpi/
elif [ "$(uname)" == "Linux" ]; then
  # in case where the OS type is Linux
  MAJOR=$(grep "#define PETSC_VERSION_MAJOR" /Users/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c34-35 )
  MINOR=$(grep "#define PETSC_VERSION_MINOR" /Users/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c34-36 )
  SUBMINOR=$(grep "#define PETSC_VERSION_SUBMINOR" /Users/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c34-35 )
  installedVersionffpetscO=$MAJOR.$MINOR.$SUBMINOR
  ffpetscDirectory=/builds/Shared/openmpi/
fi
echo "installed Version ffpetsc openmpi" $installedVersionffpetscO

#  check the version
if [ "$releaseVersionffpetsc" == "$installedVersionffpetscO" ]; then
  echo "installed release version PETSc is up to date"
else
  # change default  compiler load in update_ffpetsc_*.sh
  echo "installed release version PETSc will be upgrated"
  rm -rf $ffpetscDirectory
  ./etc/jenkins/job_Unix/check_ffpetsc/update_ffpetsc_openmpi.sh
  echo " ************* upgrading openmpi ffpetsc success *************" \
  cd 3rdparty/ff-petsc/ && rm -rf petsc-$releaseVersionffpetsc && cd ../.. && make clean
fi

#define PETSC_VERSION_SUBMINOR 4
# compilation with mpich
if [ "$(uname)" == "Darwin" ]; then
  # in case where the OS type is Darwin
  MAJOR=$(grep "#define PETSC_VERSION_MAJOR" /Users/Shared/mpich/ff-petsc/c/include/petscversion.h | cut -c34-35 )
  MINOR=$(grep    "#define PETSC_VERSION_MINOR" /Users/Shared/mpich/ff-petsc/c/include/petscversion.h | cut -c34-36 )
  SUBMINOR=$(grep "#define PETSC_VERSION_SUBMINOR" /Users/Shared/mpich/ff-petsc/c/include/petscversion.h | cut -c34-35 )
  installedVersionffpetscM=$MAJOR.$MINOR.$SUBMINOR
  ffpetscDirectory=/Users/Shared/openmpi/
elif [ "$(uname)" == "Linux" ]; then
  # in case where the OS type is Linux
  MAJOR=$(grep "#define PETSC_VERSION_MAJOR" /Users/Shared/mpich/ff-petsc/c/include/petscversion.h | cut -c34-35 )
  MINOR=$(grep    "#define PETSC_VERSION_MINOR" /Users/Shared/mpich/ff-petsc/c/include/petscversion.h | cut -c34-36 )
  SUBMINOR=$(grep "#define PETSC_VERSION_SUBMINOR" /Users/Shared/mpich/ff-petsc/c/include/petscversion.h | cut -c34-35 )
  installedVersionffpetscM=$MAJOR.$MINOR.$SUBMINOR
  ffpetscDirectory=/builds/Shared/mpich/
fi
echo "installed Version ffpetsc mpich" $installedVersionffpetscM

#  check the version
if [ "$releaseVersionffpetsc" == "$installedVersionffpetscM" ]; then
  echo "installed release version PETSc is up to date"
else
  # change default  compiler load in update_ffpetsc_*.sh
  echo "installed release version PETSc will be upgrated" \
  rm -rf $ffpetscDirectory \
  && ./etc/jenkins/job_Unix/check_ffpetsc/update_ffpetsc_mpich.sh \
  echo " ************* upgrading mpich ffpetsc success *************"
fi
