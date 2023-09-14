<!----------------------------------------------------------------------------------->
<!--- This file is part of FreeFEM.                                               --->
<!--- Laboratoire Jacques-Louis Lions                                             --->
<!--- Sorbonne Université, UMR 7598, Paris, F-75005 France                        --->
<!---                                                                             --->
<!--- FreeFEM is free software: you can redistribute it and/or modify             --->
<!--- it under the terms of the GNU Lesser General Public License as published by --->
<!--- the Free Software Foundation, either version 3 of the License, or           --->
<!--- (at your option) any later version.                                         --->
<!---                                                                             --->
<!--- FreeFEM is distributed in the hope that it will be useful,                  --->
<!--- but WITHOUT ANY WARRANTY; without even the implied warranty of              --->
<!--- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               --->
<!--- GNU Lesser General Public License for more details.                         --->
<!---                                                                             --->
<!--- You should have received a copy of the GNU Lesser General Public License    --->
<!--- along with FreeFEM.  If not, see <http://www.gnu.org/licenses/>.            --->
<!----------------------------------------------------------------------------------->

# How to compile FreeFem++ on MacOSX

##  Sep 2023

the mac version is compile wit the script 
 bin/compile-version-MacOS  gitbranch  path-where-to-compile configure-parameter
 
 configure-parameter compilation:
 generic for --enable-genenic
 debug for --enable-debug
 optim for --enable-optim
 
 and exact :
  bin/compile-version-MacOS master $HOME/tmp generic 
## Feb 2022

The Version 4.10 is compiled Under Monteray (Apple M1: arm64e) and on (Apple : x86 )

### for Apple M1 (arm64e)

Add with home brew 
1. install brew  see
[install Homebrew]{https://brew.sh/index_fr}

    
2. install this brew package to compile FreeFEM++

```bash
brew install \
aom				gettext			jpeg			libx11			pango \
autoconf		giflib			jpeg-xl			libxau			pcre \
automake		glib			libavif			libxcb			pixman \
bdw-gc			gmp				libcerf			libxdmcp		pkg-config \
bison			gnuplot			libevent		libxext			python@3.9 \
brotli			gnutls			libffi			libxrender		qt@5 \
ca-certificates		gobject-introspection		libidn2			lua			readline \
cairo			graphite2		libmpc			lzo			sqlite \
emacs			gsl				libnghttp2		m4			szip \
flex			guile			libpng			make			unbound \
fontconfig		harfbuzz		libpthread-stubs	mpdecimal		webp \
freetype		hdf5			libtasn1		mpfr			wget \
fribidi			icu4c			libtiff			nettle			xorgproto \
gcc				imath			libtool			openexr			xz \
gd				isl				libunistring	openssl@1.1		zip \
gdbm			jansson			libvmaf			p11-kit			zstd 
```

3. download the sources of  FreeFEM++ with git 

	```bash
	git clone git@github.com:FreeFem/FreeFem-sources.git ff++
	```
or the develop version 

	```bash
	git clone -b develop   git@github.com:FreeFem/FreeFem-sources.git ff++
	```
4. Compile FreeFem++ without petsc and MPI

	```bash
	cd ff++ 
	autoreconf -i 
	./configure --enable-summary --enable-download_arpack --prefix=/Applications/FreeFem++.app/Contents/ff-4.10/ --enable-download --enable-optim --enable-m64 F77=gfortran FC=gfortran CXXFLAGS=-Wno-undefined-var-template --enable-maintainer-mode
	make
	sudo make install
    ```
	
4. Compile FreeFem++ with PETSc and MPI

```bash
	cd ff++ 
	autoreconf -i 
	./configure --enable-summary --enable-download_arpack --prefix=/Applications/FreeFem++.app/Contents/ff-4.10/ --enable-download --enable-optim --enable-m64 F77=gfortran FC=gfortran CXXFLAGS=-Wno-undefined-var-template --enable-maintainer-mode
    cd 3rdparty/ff-petsc
	make 
```
    follow the output 
```bash
	cd - 
	./reconfigure 
	make -j 20 
	# check your version ...
	make -j 20  check
	# to install in /Applications/FreeFem++.app/Contents/ff-4.10/
	sudo make install -j 20 

```


## _october 2015_ (oldies)

The version 3.41 is compiled Under Yosemite (10.10.5),
xcode 7.0.1 + previous install of latex and with other brew install


Add installer Intel ifort compiler (ifort version 15.0.2 with mkl),
build mpich, libgsl, hdf5

The configure option used for all softwares are:
```bash
mpich-3.1.4: ./configure 'FC=ifort' 'F77=ifort' '-prefix=/usr/local/ff++/mpich'
gsl-1.16 ./configure 'CC=clang' '--prefix=/usr/local/ff++/mpich/'
hdf5-1.8.14 ./configure '--enable-cxx' 'CC=clang' 'CXX=clang++' '--prefix=/usr/local/ff++/mpich'
PETSc ./configure '--CFLAGS=-O2' '--COPTFLAGS=-O3' '--CXXFLAGS=-O2 -std=c++11' \
	'--CXXOPTFLAGS=-O3' '--FFLAGS=-O2' '--FOPTFLAGS=-O3' '--download-fftw' \
	'--download-hypre' '--download-metis' '--download-ml' '--download-mumps' \
	'--download-parmetis' '--download-ptscotch' '--download-scalapack' \
	'--download-suitesparse' '--download-superlu' '--prefix=/usr/local/ff++/mpich/petsc' \
	'--with-blas-lapack-lib=-L/opt/intel/mkl/lib -lmkl_intel_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_sequential -lm -lpthread' '--with-mpi-dir=/usr/local/ff++/mpich' 'PETSC_ARCH=arch-ff++'
FreeFem++: ./configure '--prefix=/usr/local/ff++/mpich/3.41' '--enable-download' '--enable-optim' '--enable-m64' \
	'F77=ifort' 'FC=ifort' '--enable-maintainer-mode' '--with-mkl=/opt/intel/mkl/lib' \
	'--with-petsc=/usr/local/ff++/mpich/petsc/lib/petsc/conf/petscvariables' '-with-hdf5=/usr/local/ff++/mpich/bin/h5cc' \
	'--with-gsl-prefix=/usr/local/ff++/mpich/' '--disable-pastix'
```

Please read INNOVATION file to see recent changes and news

Plus all previous install with brew install
```bash
499  brew install wget
500  brew install autoconf
501  brew install automake
504  brew install cmake
505  brew install gfortran
506  brew install gcc
540  brew install m4
585  brew install pkg-config
588  brew install gobject-introspection
```

LateX is install form [https://tug.org/mactex/mactex-download.html](https://tug.org/mactex/mactex-download.html).

## To compile a full version of FreeFem++ under MacOS

Under Yosemite (10.10.2, xcode 6.2 form scratch)

 * install Xcode, do clang to install command line
 * With brew install
	```bash
	499  brew install wget
	500  brew install autoconf
	501  brew install automake
	504  brew install cmake
	505  brew install gfortran
	506  brew install gcc
	540  brew install m4
	585  brew install pkg-config
	588  brew install gobject-introspection
	```
 * Install brew

Under Mavericks (10.9, xcode 5.0.2 form scratch)

 * Install Xcode 5.0.2, and the xcode command line tools,
 install Auxiliary Tools for Xcode (for PackageMaker)
 * Install Xcode command line
	```bash
	xcode-select --install
	```
 * Install gcc-4.9 form [http://hpc.sourceforge.net](http://hpc.sourceforge.net)
	```bash
	curl -O http://prdownloads.sourceforge.net/hpc/gfortran-4.9-bin.tar.gz?download
	sudo tar zxvf gfortran-4.9-bin.tar.gz -C /
	```
 * Install autoconf and automake (not in Xcode)

	I use the macport distribution form [http://www.macports.org](http://www.macports.org)
	```bash
	sudo port install autoconf
	sudo port install automake
	```
	Or with brew tool
 4) Install TeX from [ctan](http://mirrors.ctan.org/systems/mac/mactex/MacTeX.pkg)
 5) Install openmpi form the [source](http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.5.tar.bz2)
	```bash
	./configure 'CC=clang' 'CXX=clang++' 'FC=gfortran' 'F77=gfortran' --enable-ltdl-convenience
	make
	sudo make install
	```
 6) Install gsl
	```bash
	curl -O http://ftp.gnu.org/gnu/gsl/gsl-1.15.tar.gz
	tar zxvf gsl-1.15.tar.gz
	cd gsl-1.15.
	./configure CC=clang
	make
	sudo make install
	```
 7) Install git from the web [https://sourceforge.net/projects/git-osx-installer/](https://sourceforge.net/projects/git-osx-installer/)

 8) Download the sources
	```bash
	git clone git@github.com:FreeFem/FreeFem-sources.git
	```
 9) Compile FreeFem++
	```bash
	cd FreeFem-sources
	./configure '-with-suffix=macos-10.9' '-without-fltk' '--enable-download' '--enable-optim' 'MPIRUN=/usr/local/bin/mpirun' '--enable-m64' '--without-x' 'CC=clang' 'CXXFLAGS=-std=c++11' 'CXX=clang++' 'F77=/usr/local/bin/gfortran' 'FC=/usr/local/bin/gfortran' 'MPICXX=/usr/local/bin/mpic++' 'MPICC=/usr/local/bin/mpicc' 'MPIFC=/usr/local/bin/mpif90' '--enable-maintainer-mode'
	make
	sudo make install
	```
