#!/bin/bash

set -x
set -u
set -e

## Parameters
TOKEN=$1
ORGANIZATION="FreeFem"
REPOSITORY="FreeFem-sources"
VERSION=`grep AC_INIT configure.ac | cut -d"," -f2`
RELEASE_TAG_NAME="v$VERSION"
distrib=`uname -s`-`uname -r`

DISTRIB="Ubuntu"

DEB_NAME="freefem_${VERSION}-1_amd64.deb"
GH_DEB_NAME="FreeFEM_${VERSION}_amd64.deb"

## DEB build
autoreconf -i
./configure --enable-download --enable-optim --enable-generic
./3rdparty/getall -a
make -j4
echo "FreeFEM: Finite Element Language" > description-pak
sudo checkinstall -D --install=no \
    --pkgname "freefem" --pkgrelease "1" \
    --pkgversion "${VERSION}" --pkglicense "LGPL-2+" \
    --pkgsource "https://github.com/FreeFem/FreeFem-sources" \
    --pkgaltsource "https://freefem.org/" \
    --maintainer "FreeFEM, Frédéric Hecht <frederic.hecht@sorbonne-universite.fr> "  --backup=no --default
    --requires= "libc6 \(\>= 2.27),g++ \(\>= 7), gcc \(\>= 7), gfortran \(\>= 7), libgsl-dev \(\>=2.4), libhdf5-dev \(\>=1.10.0), liblapack-dev \(\>= 3.7), libopenmpi-dev \(\>=2.1.1) ,libblas-dev \(\>= 3.7.1) "

## Rename DEB to include Ubuntu version
mv $DEB_NAME $GH_DEB_NAME

## Deploy in GitHub release
RELEASE=`curl 'https://api.github.com/repos/'$ORGANIZATION'/'$REPOSITORY'/releases/tags/'$RELEASE_TAG_NAME`
UPLOAD_URL=`printf "%s" "$RELEASE" | jq -r '.upload_url'`

if [ -x $UPLOAD_URL ]
then
	echo "Release does not exists"
	exit 1
else
  RESPONSE=`curl --data-binary "@$GH_DEB_NAME" -H "Authorization: token $TOKEN" -H "Content-Type: application/octet-stream" "$UPLOAD_URL=$GH_DEB_NAME"`
fi
