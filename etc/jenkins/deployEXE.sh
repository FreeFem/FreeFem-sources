#!C:\msys64\usr\bin\bash.exe --login
source shell mingw64

set -x
set -u
set -e

## Parameters
TOKEN=$1
ORGANIZATION="FreeFem"
REPOSITORY="FreeFem-sources"
VERSION=`grep AC_INIT configure.ac | cut -d"," -f2 | tr - .`
RELEASE_TAG_NAME="v$VERSION"
EXE_NAME="Output/FreeFem++-${VERSION}-win7-64.exe"

## EXE build
autoreconf -i
./configure --enable-download --enable-optim --enable-generic
./3rdparty/getall -a
make -j4
cp AUTHORS readme/AUTHORS
touch readme/COPYING
make win32

## Deploy in GitHub release
RELEASE=`curl 'https://api.github.com/repos/'$ORGANIZATION'/'$REPOSITORY'/releases/tags/'$RELEASE_TAG_NAME`
UPLOAD_URL=`printf "%s" "$RELEASE" | jq -r '.upload_url'`

if [ -x $UPLOAD_URL ]
then
	echo "Release does not exists"
	exit 1
else
  RESPONSE=`curl --data-binary "@$EXE_NAME" -H "Authorization: token $TOKEN" -H "Content-Type: application/octet-stream" "$UPLOAD_URL=$EXE_NAME"`
fi
