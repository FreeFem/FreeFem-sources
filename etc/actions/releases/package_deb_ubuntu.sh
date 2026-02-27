#!/usr/bin/env bash

# Provided install tree (DESTDIR/usr/local) of FreeFEM, this scripts builds a DEB
# package of FreeFem and store it in GH_DEB_NAME
#
# Example for v4.16 on Ubuntu 24.04:
# FreeFEM-4.16-amd64-ubuntu24.04.deb

set -x
set -u
set -e

## Parameters
INSTALL_TREE=$1
VERSION=$2

REPOSITORY="FreeFem-sources"
OSRELEASE=$(lsb_release -r|awk '{print $2}')

DEB_NAME="freefem-${VERSION}-amd64-ubuntu${OSRELEASE}"
GH_DEB_NAME="FreeFEM-${VERSION}-amd64-ubuntu${OSRELEASE}.deb"

## Create FreeFEM Debian package
mkdir "$DEB_NAME"
mkdir "$DEB_NAME/DEBIAN"
touch "$DEB_NAME/DEBIAN/control"
{
	echo "Package: freefem";
	echo "Version: $VERSION";
	echo "Section: custom";
	echo "Architecture: amd64";
	echo "Depends: libc6 (>= 2.31), g++ (>= 9.3), gcc (>= 9.3), gfortran (>= 9.3), libgsl-dev (>=2.5), libhdf5-dev (>=1.10.4), liblapack-dev (>= 3.9), freeglut3-dev (>= 2.8.1)";
	echo "Maintainer: FreeFEM, Frédéric Hecht <frederic.hecht@sorbonne-universite.fr>";
	echo "Description: FreeFEM, Finite Element Language software";
	echo "Homepage: https://freefem.org"
} >> "$DEB_NAME/DEBIAN/control"
mkdir -p "$DEB_NAME/usr/local/share/FreeFEM"
mkdir -p "$DEB_NAME/usr/local/bin"
mkdir -p "$DEB_NAME/usr/local/lib/ff++"
mkdir -p "$DEB_NAME/usr/share/doc/freefem"

cp -r "$INSTALL_TREE/lib/ff++/$VERSION" "$DEB_NAME/usr/local/lib/ff++/$VERSION"
cp -r "$INSTALL_TREE/ff-petsc/" "$DEB_NAME/usr/local/ff-petsc"
cp -r "$INSTALL_TREE/share/FreeFEM/$VERSION" "$DEB_NAME/usr/local/share/FreeFEM/$VERSION"
cp -r "$INSTALL_TREE/bin/FreeFem++" "$DEB_NAME/usr/local/bin/FreeFem++"
cp -r "$INSTALL_TREE/bin/FreeFem++-mpi" "$DEB_NAME/usr/local/bin/FreeFem++-mpi"
cp -r "$INSTALL_TREE/bin/FreeFem++-nw" "$DEB_NAME/usr/local/bin/FreeFem++-nw"
cp -r "$INSTALL_TREE/bin/bamg" "$DEB_NAME/usr/local/bin/bamg"
cp -r "$INSTALL_TREE/bin/cvmsh2" "$DEB_NAME/usr/local/bin/cvmsh2"
cp -r "$INSTALL_TREE/bin/ff-c++" "$DEB_NAME/usr/local/bin/ff-c++"
cp -r "$INSTALL_TREE/bin/ff-get-dep" "$DEB_NAME/usr/local/bin/ff-get-dep"
cp -r "$INSTALL_TREE/bin/ff-mpirun" "$DEB_NAME/usr/local/bin/ff-mpirun"
cp -r "$INSTALL_TREE/bin/ff-pkg-download" "$DEB_NAME/usr/local/bin/ff-pkg-download"
cp -r "$INSTALL_TREE/bin/ffglut" "$DEB_NAME/usr/local/bin/ffglut"
cp -r "$INSTALL_TREE/bin/ffmaster" "$DEB_NAME/usr/local/bin/ffmaster"
cp -r "$INSTALL_TREE/bin/ffmedit" "$DEB_NAME/usr/local/bin/ffmedit"
cp AUTHORS "$DEB_NAME/usr/share/doc/freefem/AUTHOR"
cp README.md "$DEB_NAME/usr/share/doc/freefem/README.md"

dpkg-deb --build "$DEB_NAME/"
mv "$DEB_NAME.deb" "$GH_DEB_NAME"
