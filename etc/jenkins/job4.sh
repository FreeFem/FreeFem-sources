#!/bin/sh

## This job must be executed on VM2 machines
## See ./README.md

echo "Job 4"
casejob=4
# change default  compiler
change_compiler=etc/jenkins/change_compiler-`uname -s`-`uname -r`-$casejob.sh
echo try to source file  "$change_compiler"
test -f "$change_compiler" && echo  source file "$change_compiler"
test -f "$change_compiler" && cat  "$change_compiler"
test -f "$change_compiler" && source "$change_compiler"

if [ "$(uname)" == "Darwin" ]; then
  # in case where the OS type is Darwin
  PETSC_DIR='/Users/Shared/ff-petsc'
elif [ "$(uname)" == "Linux" ]; then
  # in case where the OS type is Linux
PETSC_DIR='/builds/Shared/ff-petsc'
fi

# configuration & build
autoreconf -i \
  && ./configure --enable-download --prefix=/builds/workspace/freefem \
  --with-petsc=$PETSC_DIR/real/lib \
  --with-petsc_complex=$PETSC_DIR/complex/lib \
  && ./3rdparty/getall -a \
  && chmod +x ./etc/jenkins/blob/build.sh && sh ./etc/jenkins/blob/build.sh

if [ $? -eq 0 ]
then
  echo "Build process complete"
else
  echo "Build process failed"
  exit 1
fi

# check
chmod +x ./etc/jenkins/blob/check.sh && sh ./etc/jenkins/blob/check.sh

if [ $? -eq 0 ]
then
  echo "Check process complete"
else
  echo "Check process failed"
  exit 1
fi

# install
chmod +x ./etc/jenkins/blob/install.sh && sh ./etc/jenkins/blob/install.sh

if [ $? -eq 0 ]
then
  echo "Install process complete"
else
  echo "Install process failed"
fi
