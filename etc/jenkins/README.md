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
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
```
# Ubuntu 16.04 - VM2-2
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
local Petsc 3.11.2 installed on /builds/Shared
```

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
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
```
# Ubuntu 18.04 - VM2-2
```
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libopenmpi-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
local Petsc 3.11.2 installed on /builds/Shared
```
# MacOS 10.10.5 - VM0 (Xcode 7.2.1, gcc g++ gfortran 4.9)
```
brew install gcc@4.9 m4 git flex bison
```

# MacOS 10.10.5 - VM1 (Xcode 7.2.1, gcc g++ gfortran 4.9)
```
brew install gcc@4.9 m4 git flex bison suitesparse-4.5.4 hdf5
```

[OpenMPI compilation](#openmpi) 2.1.6

# MacOS 10.10.5 - VM2 (Xcode 7.2.1, gcc g++ gfortran 4.9)
```
brew install gcc@4.9 m4 git flex bison suitesparse-4.5.4 hdf5 cmake wget autoconf automake
```
# MacOS 10.10.5 - VM2-2 (Xcode 7.2.1, gcc g++ gfortran 4.9)
```
brew install gcc@4.9 m4 git flex bison suitesparse-4.5.4 hdf5 cmake wget autoconf automake
local Petsc 3.11.2 installed on /Users/Shared
```

[OpenMPI compilation](#openmpi) 2.1.6

# MacOS 10.13 - VM0
```
brew install gcc m4 git flex bison
```

# MacOS 10.13 - VM1
```
brew install gcc m4 git flex bison suitesparse hdf5
```

[OpenMPI compilation](#openmpi)

# MacOS 10.13 - VM2
```
brew install gcc m4 git flex bison suitesparse hdf5 cmake wget autoconf automake
```
# MacOS 10.13 - VM2-2
```
brew install gcc m4 git flex bison suitesparse hdf5 cmake wget autoconf automake
local Petsc 3.11.2 installed on /Users/Shared
```

[OpenMPI compilation](#openmpi)



# Windows 7 - VM1
```
pacman -S autoconf automake-wrapper bash bash-completion \
    bison bsdcpio bsdtar bzip2 coreutils curl dash file filesystem \
    findutils flex gawk gcc gcc-fortran gcc-libs grep gzip inetutils info less lndir \
    make man-db git  mingw-w64-x86_64-toolchain mingw-w64-x86_64-hdf5 \
    mingw-w64-x86_64-arpack mingw-w64-x86_64-openblas mintty msys2-keyring msys2-launcher-git \
    msys2-runtime ncurses pacman pacman-mirrors pactoys-git patch pax-git python \
    perl pkg-config pkgfile rebase sed tar tftp-hpa time tzcode unzip util-linux which \
    mingw-w64-x86_64-libmicroutils  
```


# Windows 7 - VM3
```
pacman -S autoconf automake-wrapper bash bash-completion \
    bison bsdcpio bsdtar bzip2 coreutils curl dash file filesystem \
    findutils flex gawk gcc gcc-fortran gcc-libs grep gzip inetutils info less lndir \
    make man-db git mingw-w64-x86_64-freeglut mingw-w64-x86_64-tool-chain mingw-w64-x86_64-gsl mingw-w64-x86_64-hdf5 \
    mingw-w64-x86_64-openblas mintty msys2-keyring msys2-launcher-git \
    msys2-runtime ncurses pacman pacman-mirrors pactoys-git patch pax-git python \
    perl pkg-config pkgfile rebase sed tar tftp-hpa time tzcode unzip util-linux which \
    mingw-w64-x86_64-libmicroutils mingw-w64-x86_64-arpack 
```

[MSMPI 7](https://www.microsoft.com/en-us/download/details.aspx?id=49926)

# OpenMPI

```
curl -L https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.1.tar.gz --output openmpi-4.0.1.tar.gz
tar xf openmpi-4.0.1
cd openmpi-4.0.1/
#MacOS 10.10.5
./configure CC=clang CXX=clang++ --prefix=/usr/local  VM2 
./configure CC=gcc CXX=g++ --prefix=/usr/local (gcc g++ version5.1)  VM2-2
#MacOS 10.13.5
./configure CC=clang CXX=clang++ FC=gfortran-9 F77=gfortran-9 --prefix=/usr/local  VM2
./configure CC=gcc-9 CXX=g++-9 FC=gfortran-9 F77=gfortran-9 --prefix=/usr/local VM2-2
make -j4 all
make check
sudo make install
```

# Job 1
without mpi  (for macOS, compiled with gcc)
runs on VM0
# Job 2 
without-mpi  (for macOS, compiled with clang)
runs on VM1
# Job 3 sequential version of FreeFEM
with 3dparty, without-mpi  (for macOS, compiled with clang)
runs on VM2
# Job 4 Full version of FreeFEM compiled in release mode
with 3dparty, compilation in release mode (for macOS, compiled with gcc), PETSc installed (curent version 3.11.2) 
runs on VM2-2
# Job 5 Full version of FreeFEM compiled in debug mode
with mpi with 3dparty,, compilation in debug mode (for macOS, compiled with clang) download and install PETSc at each build
runs on VM2 every night

