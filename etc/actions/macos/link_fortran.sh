#!/bin/bash

set -e
set -x

ARCH=$1

LOCAL_PATH=/usr/local
if [ "$ARCH" = "arm" ]; then
	LOCAL_PATH=/opt/homebrew
fi

MAX_GFORTRAN_VERSIONS=$(printf '%s\n' $LOCAL_PATH/bin/gfortran-[0-9]* 2>/dev/null \
         | sort -V \
         | tail -n1)

sudo ln -fs "$MAX_GFORTRAN_VERSIONS" "$LOCAL_PATH/bin/gfortran"
