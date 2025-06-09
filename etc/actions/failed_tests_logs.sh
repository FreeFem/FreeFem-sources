#!/usr/bin/env bash

FAILED_TESTS=( $(grep :test-result examples/*/*.trs | grep FAIL) )

for val in ${FAILED_TESTS[@]/#*FAIL/}
do
    IFS="::"
    read -ra path <<< $val
    printf  "\n******** ${path%.*}.log ********\n\n"
    tail -n 100 ${path%.*}.log
    unset IFS
done;
