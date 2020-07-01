#!/bin/bash

./etc/jenkins/job_Unix/check_ffpetsc/job_ffpetsc.sh

if [ $? -eq 0 ]
then
  echo " ************* upgrading mpich ffpetsc success *************"
else
  echo " ************* upgrading mpich ffpetsc FAIL *************"
  exit 1
fi
