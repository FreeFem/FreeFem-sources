#!/bin/bash

releaseVersionffpetsc=$(grep "VERSION=" 3rdparty/ff-petsc/Makefile | cut -c9-15)
echo "release Version ffpetsc" $releaseVersionffpetsc

# compilation with openmpi
if [ "$(uname)" == "Darwin" ]; then
  # in case where the OS type is Darwin
  installedVersionffpetscO=$(grep "VERSION_GIT" /Users/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c36-42 | cut -f1 -d "\"" )   #if version X.X.XX
  ffpetscDirectory=/Users/Shared/openmpi/
elif [ "$(uname)" == "Linux" ]; then
  # in case where the OS type is Linux
  installedVersionffpetscO=$(grep "VERSION_GIT" /builds/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c36-42 | cut -f1 -d "\"" )   #if version X.X.XX
  ffpetscDirectory=/builds/Shared/openmpi/
fi
echo "installed Version ffpetsc openmpi" $installedVersionffpetscO

#  check the version
if [ "$releaseVersionffpetsc" == "$installedVersionffpetscO" ]; then
        echo "installed release version PETSc is up to date"
else
        # change default  compiler load in update_ffpetsc_*.sh
        echo "installed release version PETSc will be upgrated" \
        rm -rf $ffpetscDirectory \
        && ./etc/jenkins/job_Unix/check_ffpetsc/update_ffpetsc_openmpi.sh \
        echo " ************* upgrading openmpi ffpetsc success *************" \
        && cd 3rdparty/ff-petsc/ && make -j4 clean && cd ../..
fi


# compilation with mpich
if [ "$(uname)" == "Darwin" ]; then
  # in case where the OS type is Darwin
  installedVersionffpetscM=$(grep "VERSION_GIT" /Users/Shared/mpich/ff-petsc/c/include/petscversion.h | cut -c36-42 | cut -f1 -d "\"" )   #if version X.X.XX
  ffpetscDirectory=/Users/Shared/mpich/
elif [ "$(uname)" == "Linux" ]; then
  # in case where the OS type is Linux
  installedVersionffpetscM=$(grep "VERSION_GIT" /builds/Shared/mpich/ff-petsc/c/include/petscversion.h | cut -c36-42 | cut -f1 -d "\"" )   #if version X.X.XX
  ffpetscDirectory=/builds/Shared/mpich/
fi
echo installed Version ffpetsc mpich "$installedVersionffpetscM"

#  check the version
if [ "$releaseVersionffpetsc" == "$installedVersionffpetscM" ]; then
        echo "installed release version PETSc is up to date"
else
        # change default  compiler load in update_ffpetsc_*.sh
        echo "installed release version PETSc will be upgrated" \
        rm -rf $ffpetscDirectory
        echo "************* upgrading mpich ffpetsc  *************" \
        && ./etc/jenkins/job_Unix/check_ffpetsc/update_ffpetsc_mpich.sh \
        echo " ************* upgrading mpich ffpetsc success *************" \
fi
