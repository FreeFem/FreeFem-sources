#!/usr/bin/env bash

## This job must be executed on VM2 machines
## See ./README.md

echo "Job 3"
casejob=3

set -e

# change default  compiler
change_compiler=etc/jenkins/change_compiler/change_compiler-`uname -s`-`uname -r`-$casejob.sh
echo try to source file  "$change_compiler"
test -f "$change_compiler" && echo  source file "$change_compiler"
test -f "$change_compiler" && cat  "$change_compiler"
test -f "$change_compiler" && source "$change_compiler"


# configuration & build
autoreconf -i \
  && ./configure --enable-download --without-mpi --prefix=/builds/workspace/freefem \
  && ./3rdparty/getall -a -o ARPACK,METIS,ParMETIS,ScaLAPACK,Scotch,SuiteSparse,SuperLU,mmg,parmmg,hpddm,bemtool,Boost,libpthread-google,TetGen,Ipopt,NLopt,freeYams,FFTW,Gmm++,MMG3D,mshmet,MUMPS,htool \
  && ./etc/jenkins/blob/build.sh

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

# visu for jenkins tests results analyser
./etc/jenkins/resultForJenkins/resultForJenkins.sh
