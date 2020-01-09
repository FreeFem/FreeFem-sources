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

if [ "$distrib" == "Linux-4.4.0-166-generic" ]; then
  # 16.04
DISTRIB="Ubuntu_16.04"
elif [ "$distrib" == "Linux-4.15.0-51-generic" ]; then
  # 18.04
DISTRIB="Ubuntu_18.04"
fi


DEB_NAME="freefem_${VERSION}-1_amd64.deb"
GH_DEB_NAME="FreeFEM_${VERSION}_${DISTRIB}_amd64.deb"

## DEB build
autoreconf -i
./configure --enable-download --enable-optim --enable-generic --prefix=/builds/workspace/freefem
./3rdparty/getall -a -o PETSc,Ipopt,NLopt,freeYams,FFTW,ARPACK,Gmm++,MMG3D,mshmet,MUMPS,htool
cd 3rdparty/ff-petsc && make petsc-slepc && cd -
./reconfigure
make -j4
echo "FreeFEM: Finite Element Language" > description-pak
sudo checkinstall -D --install=no \
    --pkgname "freefem" --pkgrelease "1" \
    --pkgversion "${VERSION}" --pkglicense "LGPL-2+" \
    --pkgsource "https://github.com/FreeFem/FreeFem-sources" \
    --pkgaltsource "https://freefem.org/" \
    --include=/builds/workspace/freefem/ff-petsc/ \
    --maintainer "FreeFEM" --backup=no --default

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
