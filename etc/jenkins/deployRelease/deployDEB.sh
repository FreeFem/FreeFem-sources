#!/usr/bin/env bash

set -x
set -u
set -e

## Parameters
TOKEN=$1
ORGANIZATION="FreeFem"
REPOSITORY="FreeFem-sources"
VERSION=$(grep AC_INIT configure.ac | cut -d"," -f2 | cut -d"[" -f2 | cut -d"]" -f1)
RELEASE_TAG_NAME="v$VERSION"
OSRELEASE=$(lsb_release -r|awk '{print $2}')

BUILD_DIR="/usr/local"

DEB_NAME="freefem-${VERSION}-amd64-ubuntu${OSRELEASE}"
GH_DEB_NAME="FreeFEM-${VERSION}-amd64-ubnutu${OSRELEASE}.deb"

## DEB build
autoreconf -i
./configure --enable-download --enable-optim --enable-generic --prefix="$BUILD_DIR"
./3rdparty/getall -a -o PETSc,Ipopt,NLopt,freeYams,FFTW,Gmm++,MMG3D,mshmet,MUMPS
cd 3rdparty/ff-petsc && make petsc-slepc && cd -
./reconfigure
make -j"$(nproc)"
make -j"$(nproc)" install

## Create FreeFEM Debian package
mkdir "$DEB_NAME"
mkdir "$DEB_NAME/DEBIAN"
touch "$DEB_NAME/DEBIAN/control"
{
	echo "Package: freefem";
	echo "Version: $VERSION";
	echo "Section: custom";
	echo "Architecture: amd64";
	echo "Depends: libc6 (>= 2.31), g++ (>= 9.3), gcc (>= 9.3), gfortran (>= 9.3), libgsl-dev (>=2.5), libhdf5-dev (>=1.10.4), liblapack-dev (>= 3.9), libopenmpi-dev (>=4.0.3) ,freeglut3-dev (>= 2.8.1)";
	echo "Maintainer: FreeFEM, Frédéric Hecht <frederic.hecht@sorbonne-universite.fr>";
	echo "Description: FreeFEM, Finite Element Language software";
	echo "Homepage: https://freefem.org";
} >> "$DEB_NAME/DEBIAN/control"
mkdir -p "$DEB_NAME/usr/local/share/FreeFEM"
mkdir -p "$DEB_NAME/usr/local/bin"
mkdir -p "$DEB_NAME/usr/local/lib/ff++"
mkdir -p "$DEB_NAME/usr/share/doc/freefem"

cp -r "$BUILD_DIR/lib/ff++/$VERSION" "$DEB_NAME/usr/local/lib/ff++/$VERSION"
cp -r "$BUILD_DIR/ff-petsc/" "$DEB_NAME/usr/local/ff-petsc"
cp -r "$BUILD_DIR/share/FreeFEM/$VERSION" "$DEB_NAME/usr/local/share/FreeFEM/$VERSION"
cp -r "$BUILD_DIR/bin/FreeFem++" "$DEB_NAME/usr/local/bin/FreeFem++"
cp -r "$BUILD_DIR/bin/FreeFem++-mpi" "$DEB_NAME/usr/local/bin/FreeFem++-mpi"
cp -r "$BUILD_DIR/bin/FreeFem++-nw" "$DEB_NAME/usr/local/bin/FreeFem++-nw"
cp -r "$BUILD_DIR/bin/bamg" "$DEB_NAME/usr/local/bin/bamg"
cp -r "$BUILD_DIR/bin/cvmsh2" "$DEB_NAME/usr/local/bin/cvmsh2"
cp -r "$BUILD_DIR/bin/ff-c++" "$DEB_NAME/usr/local/bin/ff-c++"
cp -r "$BUILD_DIR/bin/ff-get-dep" "$DEB_NAME/usr/local/bin/ff-get-dep"
cp -r "$BUILD_DIR/bin/ff-mpirun" "$DEB_NAME/usr/local/bin/ff-mpirun"
cp -r "$BUILD_DIR/bin/ff-pkg-download" "$DEB_NAME/usr/local/bin/ff-pkg-download"
cp -r "$BUILD_DIR/bin/ffglut" "$DEB_NAME/usr/local/bin/ffglut"
cp -r "$BUILD_DIR/bin/ffmaster" "$DEB_NAME/usr/local/bin/ffmaster"
cp -r "$BUILD_DIR/bin/ffmedit" "$DEB_NAME/usr/local/bin/ffmedit"
cp AUTHORS "$DEB_NAME/usr/share/doc/freefem/AUTHOR"
cp README.md "$DEB_NAME/usr/share/doc/freefem/README.md"

dpkg-deb --build "$DEB_NAME/"
mv "$DEB_NAME.deb" "$GH_DEB_NAME"

## Deploy in GitHub release
RELEASE=$(curl "https://api.github.com/repos/$ORGANIZATION/$REPOSITORY/releases/tags/$RELEASE_TAG_NAME")
UPLOAD_URL=$(printf "%s" "$RELEASE" | jq -r '.upload_url')

if [ -x "$UPLOAD_URL" ]
then
	echo "Release does not exists"
	exit 1
else
	RESPONSE=$(curl --data-binary "@$GH_DEB_NAME" -H "Authorization: token $TOKEN" -H "Content-Type: application/octet-stream" "$UPLOAD_URL=$GH_DEB_NAME")
	echo "Github response:"
	echo "$RESPONSE"
fi

# clean the VM
rm -rf "$DEB_NAME"
rm -rf "$BUILD_DIR/ff++"
rm -rf "$BUILD_DIR/ff-petsc"
rm "$GH_DEB_NAME"
