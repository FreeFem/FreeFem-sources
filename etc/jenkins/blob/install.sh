#!/bin/bash

# if sudo is needed
# launch ./install.sh sudo

WITH_SUDO=""
if [ "$1" = "sudo" ]
then
	WITH_SUDO=sudo
fi

${WITH_SUDO} make install
