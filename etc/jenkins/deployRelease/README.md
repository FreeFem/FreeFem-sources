```
Jenkins project: https://ci.inria.fr/freefem/
```

# Ubuntu system
minimal configuration:
g++ v7, gcc v7, gfortran v7, libgsl-dev v2.4, libhdf5-dev v1.10.0, liblapack-dev v3.7, libopenmpi-dev v2.1.1, freeglut3-dev v2.8.1

##FreeFEM_X.X-X_Ubuntu_amd64.deb
```
- release built on ubuntu 16.04 (Ubuntu Xenial 16.04 x64)
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
```

##FreeFEM_X.X-X_Ubuntu_withPETSc_amd64.deb
```
- release built on ubuntu 16.04 (Ubuntu Xenial 16.04 x64)
sudo apt install gcc g++ gfortran m4 patch git libblas-dev liblapack-dev libsuitesparse-dev libhdf5-dev libgsl-dev flex bison wget cmake autoconf automake autotools-dev
- contains the last PETSc version (binaries provide by FreeFEM compilation)  
```
# Windows system
##FreeFEM_X.X-X_-win64.exe
```
install msys2, git MS-MPI v10.1.2 inno
- release built on Windows 7 64bits
pacman -S autoconf automake-wrapper bash bash-completion \
    bison bsdcpio bsdtar bzip2 coreutils curl dash file filesystem \
    findutils flex gawk gcc gcc-fortran gcc-libs grep gzip inetutils info less lndir \
    make man-db git mingw-w64-x86_64-freeglut mingw-w64-x86_64-toolchain mingw-w64-x86_64-gsl mingw-w64-x86_64-hdf5 \
    mingw-w64-x86_64-openblas mintty msys2-keyring msys2-launcher-git \
    msys2-runtime ncurses pacman pacman-mirrors pactoys-git patch pax-git python \
    perl pkg-config pkgfile rebase sed tar tftp-hpa time tzcode unzip util-linux which \
    mingw-w64-x86_64-libmicroutils mingw-w64-x86_64-jq
pacman -R mingw-w64-x86_64-python mingw-w64-x86_64-gdb
```
# MacOS system
