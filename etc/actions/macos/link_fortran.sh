#!/bin/bash

set -e
set -x

ARCH=$1

LOCAL_PATH=/usr/local
if [ "$ARCH" = "arm" ]; then
	LOCAL_PATH=/opt/homebrew
fi

sudo ln -fs "$LOCAL_PATH/bin/gfortran-14" "$LOCAL_PATH/bin/gfortran"
