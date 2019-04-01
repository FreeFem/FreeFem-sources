#!/bin/sh

# if sudo is needed
# lanuch ./clean.sh sudo

make clean
$SUDO rm -Rf /builds/freefem

if [ $? -eq 0 ]
then
	echo "Cleanup process succed"
else
	echo "Cleanup process failed"
fi
