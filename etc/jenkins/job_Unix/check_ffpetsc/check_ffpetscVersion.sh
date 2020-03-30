#!/bin/bash

releaseVersionffpetsc=$(grep "VERSION=" 3rdparty/ff-petsc/Makefile | cut -c9-15)
echo "release Version ffpetsc" $releaseVersionffpetsc

if [ "$(uname)" == "Darwin" ]; then
  # in case where the OS type is Darwin
  installedVersionffpetsc=$(grep "VERSION=" /Users/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c36-42 | cut -f1 -d "\"" )   #if version X.X.XX
elif [ "$(uname)" == "Linux" ]; then
  # in case where the OS type is Linux
  installedVersionffpetsc=$(grep "VERSION=" /builds/Shared/openmpi/ff-petsc/c/include/petscversion.h | cut -c36-42 | cut -f1 -d "\"" )   #if version X.X.XX
fi
echo "installed Version ffpetsc" $installedVersionffpetsc

#  check the version
if [ "$releaseVersionffpetsc" == "$installedVersionffpetsc" ]; then
        echo "installed release version PETSc is up to date"
else
        # change default  compiler load in update_ffpetsc_*.sh
        echo "installed release version PETSc will be upgrated"
        echo "************* upgrading openmpi ffpetsc  *************" \
        ./etc/jenkins/job_Unix/check_ffpetsc/update_ffpetsc_openmpi.sh \
        echo " ************* upgrading mpich ffpetsc success *************" \
        cd 3rdparty/ff-petsc/ && make -j4 clean && cd ../.. \
        echo " ************* upgrading mpich ffpetsc  *************" \
    ./etc/jenkins/job_Unix/check_ffpetsc/update_ffpetsc_mpich.sh \
        echo " ************* upgrading mpich ffpetsc success *************"
fi

