#!/usr/bin/env bash

echo "Job 5 (openmpi)"
set -e

casejob=5_openmpi

# change default compiler
if [ "$(uname)" = "Darwin" ]
then
  # MacOS
  change_compiler=etc/jenkins/change_compiler/change_compiler-$(uname -s)-$(uname -r)-$casejob.sh
elif [ "$(uname)" == "Linux" ]
then
  # Linux
  change_compiler=etc/jenkins/change_compiler/change_compiler-$(uname -s)-$casejob.sh
fi

echo "Try to source file $change_compiler"
test -f "$change_compiler" && echo "Source file $change_compiler"
test -f "$change_compiler" && cat "$change_compiler"
test -f "$change_compiler" && source "$change_compiler"

# configuration & build
autoreconf -i
./configure --enable-download --enable-debug --prefix="$WORKSPACE/$JOB_NAME/install"
./3rdparty/getall -a -o PETSc,Ipopt,NLopt,freeYams,FFTW,Gmm++,MMG3D,mshmet,MUMPS
./etc/jenkins/blob/build_PETSc.sh
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
