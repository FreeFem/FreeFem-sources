#!/bin/bash


releaseVersionffpetsc=$(grep "VERSION=" 3rdparty/ff-petsc/Makefile) | cut -c9-15

if [ "$(uname)" == "Darwin" ]; then
  # in case where the OS type is Darwin
  installedVersionffpetsc=$(grep "VERSION=" /Users/Shared/openmpi/ff-petsc/c/include/petscversion.h) | cut -c36-42 | cut -f1 -d "\""    #if version X.X.XX
elif [ "$(uname)" == "Linux" ]; then
  # in case where the OS type is Linux
  installedVersionffpetsc=$(grep "VERSION=" /builds/Shared/openmpi/ff-petsc/c/include/petscversion.h) | cut -c36-42 | cut -f1 -d "\""    #if version X.X.XX
fi
#  check the version
if [ "$releaseVersionffpetsc" == "$installedVersionffpetsc" ]; then
	echo "installed release version PETSc is up to date"
	updatescript=0	
else
	echo "installed release version PETSc will be upgrated"
	updatescript=1
fi
# change default  compiler load in update_ffpetsc_*.sh
test $updatescript -gt 1 \
	&& echo "************* upgrading openmpi ffpetsc  *************" \
	&& ./etc/jenkins/job_Unix/check_ffpetsc/update_ffpetsc_openmpi.sh \
	&& echo " ************* upgrading mpich ffpetsc success *************" \
	&& cd 3rdparty/ff-petsc/ && make -j4 clean && cd ../.. \
	&& echo " ************* upgrading mpich ffpetsc  *************" \
    && ./etc/jenkins/job_Unix/check_ffpetsc/update_ffpetsc_mpich.sh \
	&& echo " ************* upgrading mpich ffpetsc success *************"
	

