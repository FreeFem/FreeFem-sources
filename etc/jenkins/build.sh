#!/bin/sh

# Input arguments
# ./build.sh [PETSc (0/1), default 0] [sudo (0/1), default 0] [install (0), default 1]
withPETSc=0
withSudo=0
if [ "$#" -eq 0 ]
then
	echo "No arguments provided"
	echo "Run without petsc and sudo"
else
	if [ "$1" -eq 0 ] || [ "$1" -eq 1 ]
	then
		echo "PETSc argument: $1"
		withPETSc=$1
	fi

	if [ "$2" -eq 0 ] || [ "$2" -eq 1 ]
	then
		echo "Sudo argument: $2"
		withSudo=$2
	fi

	if [ "$3" -eq 0 ]
	then
		echo "Install argument: $3"
		withoutInstall=$3
	fi
fi

if [ $withSudo -eq 1 ]
then
	SUDO=sudo
else
	SUDO=""
fi

# Build process
if [ $withPETSc -eq 1 ]
then
	autoreconf -i \
		&& ./configure --enable-download --enable-optim --prefix=/builds/freefem \
		&& ./3rdparty/getall -a \
		&& cd 3rdparty/ff-petsc \
		&& make petsc-slepc SUDO=$SUDO \
		&& cd - \
		&& ./reconfigure \
		&& make -j4 \
		&& make check \
		&& if [ $withoutInstall ]; then echo "Without install"; else $SUDO make install; fi
else
	autoreconf -i \
		&& ./configure --enable-download --enable-optim --prefix=/builds/freefem \
		&& ./3rdparty/getall -a \
		&& make -j4 \
		&& make check \
		&& if [ $withoutInstall ]; then echo "Without install"; else $SUDO make install; fi
fi

if [ $? -eq 0 ]
then
	echo "Compilation process succed"
else
	exit 1
fi

# Clean process (not fail script if cleanup fails)
make clean
$SUDO rm -Rf /builds/freefem

if [ $? -eq 0 ]
then
	echo "Cleanup process succed"
else
	echo "Cleanup process failed"
fi

echo "That's all folks!"
