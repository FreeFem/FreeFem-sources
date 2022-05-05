#!/usr/bin/env bash

echo "Job 3"

# change default compiler
if [ "$(uname)" = "Darwin" ]
then
  export CC=clang
  export CXX=clang++
  export FC=gfortran
  export F77=gfortran
fi

# configuration & build
autoreconf -i
./configure --enable-download --without-mpi --prefix="$WORKSPACE/$JOB_NAME/install"
./3rdparty/getall -a -o ARPACK,METIS,ParMETIS,ScaLAPACK,Scotch,SuiteSparse,SuperLU,mmg,parmmg,hpddm,bemtool,Boost,libpthread-google,TetGen,Ipopt,NLopt,freeYams,FFTW,Gmm++,MMG3D,mshmet,MUMPS,htool
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
