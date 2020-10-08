# Ubuntu 16.04 - VM0
```
sudo apt install gcc g++ gfortran m4 patch git liblapack-dev flex bison
```

# Ubuntu 16.04 - VM1
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison
```

# Ubuntu 16.04 - VM2
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
```
[OpenMPI compilation](#openmpi)4.0.2
[MPICH compilation](#mpich)3.3.2

# Ubuntu 16.04 - VM2-2
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
local Petsc 3.12.4 installed on /Users/Shared (mpich and openmpi version)
```
[OpenMPI compilation](#openmpi)4.0.2
[MPICH compilation](#mpich)3.3.2

# Ubuntu 18.04 - VM0
```
sudo apt install gcc g++ gfortran m4 patch git liblapack-dev flex bison
```

# Ubuntu 18.04 - VM1
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison
```

# Ubuntu 18.04 - VM2
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
```
[OpenMPI compilation](#openmpi)4.0.2
[MPICH compilation](#mpich)3.3.2

# Ubuntu 18.04 - VM2-2
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
local Petsc 3.12.4 installed on /Users/Shared (mpich and openmpi version)
```
[OpenMPI compilation](#openmpi)4.0.3
[MPICH compilation](#mpich)3.3.2

# Ubuntu 19.10 - VM0
```
sudo apt install gcc g++ gfortran m4 patch git liblapack-dev flex bison
```

# Ubuntu 19.10 - VM1
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison
```

# Ubuntu 19.10 - VM2
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
```
[OpenMPI compilation](#openmpi)4.0.2
[MPICH compilation](#mpich)3.3.2

# Ubuntu 19.10 - VM2-2
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
local Petsc 3.12.4 installed on /builds/Shared/ (mpich and openmpi version)
```
[OpenMPI compilation](#openmpi)4.0.2
[MPICH compilation](#mpich)3.3.2
# MacOS 10.10.5 - VM0 (Xcode 7.2.1, gcc g++ gfortran 5.1)
```
brew install gcc@4.9 m4 git flex bison
```

# MacOS 10.10.5 - VM1 (Xcode 7.2.1, gcc g++ gfortran 5.1)
```
brew install gcc@4.9 m4 git flex bison suitesparse-4.5.4 hdf5
```

# MacOS 10.10.5 - VM2 (Xcode 7.2.1, gcc g++ gfortran 5.1)
```
brew install gcc@4.9 m4 git flex bison suitesparse-4.5.4 hdf5 cmake wget autoconf automake
```
[OpenMPI compilation](#openmpi) 2.1.6
[MPICH compilation](#mpich)3.3.2
# MacOS 10.10.5 - VM2-2 (Xcode 7.2.1, gcc g++ gfortran 5.1)
```
brew install gcc@4.9 m4 git flex bison suitesparse-4.5.4 hdf5 cmake wget autoconf automake
local Petsc 3.12.4 installed on /Users/Shared (mpich and openmpi version)
```
[OpenMPI compilation](#openmpi) 2.1.6
[MPICH compilation](#mpich)3.3.2
# MacOS 10.13 - VM0
```
brew install gcc m4 git flex bison
```

# MacOS 10.13 - VM1
```
brew install gcc m4 git flex bison suitesparse hdf5
```

# MacOS 10.13 - VM2
```
brew install gcc m4 git flex bison suitesparse hdf5 cmake wget autoconf automake
```
[OpenMPI compilation](#openmpi) 4.0.1
[MPICH compilation](#mpich)3.3.2

# MacOS 10.13 - VM2-2
```
brew install gcc m4 git flex bison suitesparse hdf5 cmake wget autoconf automake
local Petsc 3.12.4 installed on /Users/Shared (mpich and openmpi version)
```
[OpenMPI compilation](#openmpi)4.0.1
[MPICH compilation](#mpich)3.3.2

# MacOS 10.14 - VM0
```
brew install gcc m4 git flex bison
```

# MacOS 10.14 - VM1
```
brew install gcc m4 git flex bison suitesparse hdf5
```

# MacOS 10.14 - VM2
```
brew install gcc m4 git flex bison suitesparse hdf5 cmake wget autoconf automake
```
[OpenMPI compilation](#openmpi) 4.0.3
[MPICH compilation](#mpich)3.3.2

# MacOS 10.14 - VM2-2
```
brew install gcc m4 git flex bison suitesparse hdf5 cmake wget autoconf automake
local Petsc 3.12.4 installed on /Users/Shared (mpich and openmpi version)
```
[OpenMPI compilation](#openmpi)4.0.3
[MPICH compilation](#mpich)3.3.2

# Windows 7 and 10 - all VM
```
pacman -S autoconf automake-wrapper bash bash-completion \
    bison bsdcpio bsdtar bzip2 coreutils curl dash file filesystem \
    findutils flex gawk gcc gcc-fortran gcc-libs grep gzip inetutils info less lndir \
    make man-db git mingw-w64-x86_64-freeglut mingw-w64-x86_64-toolchain mingw-w64-x86_64-gsl mingw-w64-x86_64-hdf5 \
    mingw-w64-x86_64-openblas mintty msys2-keyring msys2-launcher-git \
    msys2-runtime ncurses pacman pacman-mirrors pactoys-git patch pax-git python \
    perl pkg-config pkgfile rebase sed tar tftp-hpa time tzcode unzip util-linux which \
    mingw-w64-x86_64-libmicroutils mingw-w64-x86_64-cmake
pacman -R mingw-w64-x86_64-gcc-ada mingw-w64-x86_64-gcc-objc mingw-w64-x86_64-gdb mingw-w64-x86_64-python
```
[MSMPI 10.1.2](https://www.microsoft.com/en-us/download/details.aspx?id=100593)


# OpenMPI

```
#MacOS 10.10.5
curl -L https://download.open-mpi.org/release/open-mpi/v2.1/openmpi-2.1.6.tar.gz--output openmpi-2.1.6.tar.gz
tar xf openmpi-2.1.6.tar.gz
cd openmpi-2.1.6/
./configure CC=clang CXX=clang++ --prefix=/usr/local  VM2 
./configure CC=gcc CXX=g++ --prefix=/usr/local (gcc g++ version5.1)  VM2-2
#MacOS 10.13.5
curl -L https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.1.tar.gz --output openmpi-4.0.1.tar.gz
tar xf openmpi-4.0.1.tar.gz
cd openmpi-4.0.1/
./configure CC=clang CXX=clang++ FC=gfortran-9 F77=gfortran-9 --prefix=/usr/local  VM2
./configure CC=gcc-9 CXX=g++-9 FC=gfortran-9 F77=gfortran-9 --prefix=/usr/local VM2-2
make -j4 all
make check
sudo make install
#Ubuntu 16.04 - Ubuntu 18.04 
openmpi 4.0.2 installed
# Ubuntu 19.10
openmpi 4.0.3 installed

```
# MPICH

```
curl -L http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz --output mpich-3.3.2.tar.gz
tar xf mpich-3.3.2.tar.gz
cd mpich-3.3.2/
#MacOS 10.10.5
./configure CC=clang CXX=clang++ --prefix=/usr/local/mpich3 VM2
./configure CC=gcc CXX=g++ --prefix=/usr/local/mpich3 (gcc g++ version5.1)  VM2-2
#MacOS 10.13.5
./configure CC=clang CXX=clang++ FC=gfortran-9 F77=gfortran-9 --prefix=/usr/local/mpich3 VM2  (update v3.3.2)
./configure CC=gcc-9 CXX=g++-9 FC=gfortran-9 F77=gfortran-9 --prefix=/usr/local/mpich3 VM2-2
make -j4 all
make check
sudo make install

#Ubuntu 16.04 - Ubuntu 18.04 - Ubuntu 19.10
wget http://www.mpich.org/static/downloads/3.3.1/mpich-3.3.1.tar.gz
tar xf mpich-3.3.1.tar.gz
cd mpich-3.3.1/
./configure --prefix=/usr/local/mpich3
make -j4
make -j4 install

```

Jenkins projets are [here](https://ci.inria.fr/freefem-dev/) and [here](https://ci.inria.fr/freefem-dev-2/) 


# Job 1
without mpi  (for macOS, compiled with gcc)
runs on VM0
# Job 2 
without-mpi  (for macOS, compiled with clang)
runs on VM1
# Job 3 sequential version of FreeFEM
with 3dparty, without-mpi  (for macOS, compiled with clang)
runs on VM2
# Job 4_openmpi Full version of FreeFEM compiled in release mode
with 3dparty, compilation in release mode (for macOS, compiled with gcc), uses a PETSc installed (curent version 3.12.4) 
runs on VM2-2
# Job 4_mpich Full version of FreeFEM compiled in release mode
with 3dparty, compilation in release mode (for macOS, compiled with gcc), uses a PETSc installed (curent version 3.12.42) 
runs on VM2-2
# Job 5_openmpi Full version of FreeFEM compiled in debug mode
with mpi with 3dparty,, compilation in debug mode (for macOS, compiled with clang) download and install PETSc at each build
runs on VM2
# Job 5_mpich Full version of FreeFEM compiled in debug mode
with mpi with 3dparty,, compilation in debug mode (for macOS, compiled with clang) download and install PETSc at each build
runs on VM2





