#!/bin/sh

set -e
set -x

VERSION=$1
DEB_FOLDER=$2
DEB_NAME=$3

mkdir -p "$DEB_FOLDER/DEBIAN"
touch "$DEB_FOLDER/DEBIAN/control"
{
	echo "Package: freefem"
	echo "Version: $VERSION"
	echo "Section: custom"
	echo "Architecture: amd64"
	echo "Depends: libc6 (>= 2.31), g++ (>= 9.3), gcc (>= 9.3), gfortran (>= 9.3), libgsl-dev (>=2.5), libhdf5-dev (>=1.10.4), liblapack-dev (>= 3.9), libopenmpi-dev (>=4.0.3) ,freeglut3-dev (>= 2.8.1)"
	echo "Maintainer: FreeFEM, Frédéric Hecht <frederic.hecht@sorbonne-universite.fr>"
	echo "Description: FreeFEM, Finite Element Language software"
	echo "Homepage: https://freefem.org"
} >>"$DEB_FOLDER/DEBIAN/control"
mkdir -p "$DEB_FOLDER/usr/local/share/FreeFEM"
mkdir -p "$DEB_FOLDER/usr/local/bin"
mkdir -p "$DEB_FOLDER/usr/local/lib/ff++"
mkdir -p "$DEB_FOLDER/usr/share/doc/freefem"

cp -r "/usr/local/lib/ff++/$VERSION" "$DEB_FOLDER/usr/local/lib/ff++/$VERSION"
cp -r "/usr/local/ff-petsc/" "$DEB_FOLDER/usr/local/ff-petsc"
cp -r "/usr/local/share/FreeFEM/$VERSION" "$DEB_FOLDER/usr/local/share/FreeFEM/$VERSION"
cp -r "/usr/local/bin/FreeFem++" "$DEB_FOLDER/usr/local/bin/FreeFem++"
cp -r "/usr/local/bin/FreeFem++-mpi" "$DEB_FOLDER/usr/local/bin/FreeFem++-mpi"
cp -r "/usr/local/bin/FreeFem++-nw" "$DEB_FOLDER/usr/local/bin/FreeFem++-nw"
cp -r "/usr/local/bin/bamg" "$DEB_FOLDER/usr/local/bin/bamg"
cp -r "/usr/local/bin/cvmsh2" "$DEB_FOLDER/usr/local/bin/cvmsh2"
cp -r "/usr/local/bin/ff-c++" "$DEB_FOLDER/usr/local/bin/ff-c++"
cp -r "/usr/local/bin/ff-get-dep" "$DEB_FOLDER/usr/local/bin/ff-get-dep"
cp -r "/usr/local/bin/ff-mpirun" "$DEB_FOLDER/usr/local/bin/ff-mpirun"
cp -r "/usr/local/bin/ff-pkg-download" "$DEB_FOLDER/usr/local/bin/ff-pkg-download"
cp -r "/usr/local/bin/ffglut" "$DEB_FOLDER/usr/local/bin/ffglut"
cp -r "/usr/local/bin/ffmaster" "$DEB_FOLDER/usr/local/bin/ffmaster"
cp -r "/usr/local/bin/ffmedit" "$DEB_FOLDER/usr/local/bin/ffmedit"
cp AUTHORS "$DEB_FOLDER/usr/share/doc/freefem/AUTHOR"
cp README.md "$DEB_FOLDER/usr/share/doc/freefem/README.md"

dpkg-deb --build "$DEB_FOLDER/"
mv "$DEB_FOLDER.deb" "$DEB_NAME"
