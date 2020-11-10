#!/usr/bin/env bash

## Jenkins variables
workspace=$WORKSPACE
job=$JOB_NAME

## Tests directories
logDirectory="${workspace}/log/${job}"
baseDirectory="examples"
declare -a directories
directories=("3d" "3dSurf" "bamg" "bug" "eigen" "examples" "ffddm" "hpddm" "misc" "mpi" "other" "plugin" "tutorial")

echo "Log will be put in ${logDirectory}"

## Create tree and copy log
for directory in "${directories[@]}"
do
    ## Create directory
    mkdir -p ${logDirectory}/${directory}
    ## Remove files, if any
    rm -f ${logDirectory}/${directory}/*
    ## Copy log
    cp -f ${baseDirectory}/${directory}/*.err ${logDirectory}/${directory}/ 2>/dev/null
done

## Exit with 0 if the last cp failed (just no error log)
exit 0
