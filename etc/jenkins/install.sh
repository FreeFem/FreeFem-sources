#!/bin/sh

# if sudo is needed
# lanuch ./install.sh sudo

WITH_SUDO=""
if [ "$1" = "sudo" ]
then
	WITH_SUDO=sudo
fi

${WITH_SUDO} make install
