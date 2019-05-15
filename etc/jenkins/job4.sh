#!/bin/bash

## This job must be executed on VM2 machines
## See ./README.md

echo "Job 4"

# configuration & build
autoreconf -i \
  && ./configure --prefix=/builds/workspace/freefem \
  && chmod +x ./etc/jenkins/blob/build.sh && bash ./etc/jenkins/blob/build.sh

if [ $? -eq 0 ]
then
  echo "Build process complete"
else
  echo "Build process failed"
  exit 1
fi

# check
chmod +x ./etc/jenkins/blob/check.sh && bash ./etc/jenkins/blob/check.sh

if [ $? -eq 0 ]
then
  echo "Check process complete"
else
  echo "Check process failed"
fi

# install
chmod +x ./etc/jenkins/blob/install.sh && bash ./etc/jenkins/blob/install.sh

if [ $? -eq 0 ]
then
  echo "Install process complete"
else
  echo "Install process failed"
fi
