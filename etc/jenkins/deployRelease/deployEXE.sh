#!C:\msys64\usr\bin\bash.exe --login
source shell mingw64

set -x
set -u
set -e

## Parameters
TOKEN=$1
ORGANIZATION="FreeFem"
REPOSITORY="FreeFem-sources"
VERSION=`grep AC_INIT configure.ac | cut -d"," -f2`
RELEASE_TAG_NAME="v$VERSION"
EXE_NAME="FreeFem++-${VERSION}-win64.exe"
GH_EXE_NAME="FreeFEM-${VERSION}-win7-64.exe"

## EXE build
autoreconf -i
./configure --enable-download --enable-optim --enable-generic
./3rdparty/getall -a -o PETSc,Ipopt,NLopt,freeYams,FFTW,Gmm++,MMG3D,mshmet,MUMPS,htool
## compile and install ff-petsc
cd 3rdparty/ff-petsc && make petsc-slepc && cd -
./reconfigure

make
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
    mv Output/$EXE_NAME $GH_EXE_NAME
    RESPONSE=`curl --data-binary "@$GH_EXE_NAME" -H "Authorization: token $TOKEN" -H "Content-Type: application/octet-stream" "$UPLOAD_URL=$GH_EXE_NAME"`
fi
