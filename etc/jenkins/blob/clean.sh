#!/bin/bash

# if sudo is needed
# lanuch ./clean.sh sudo

WITH_SUDO=""
if [ "$1" = "sudo" ]
then
	WITH_SUDO=sudo
fi

make clean
${WITH_SUDO} rm -Rf /builds/freefem

if [ $? -eq 0 ]
then
	echo "Cleanup process succed"
else
	echo "Cleanup process failed"
fi
