# Ubuntu

**Current system:** Ubuntu 20.04

**Jenkins dashboard:** https://ci.inria.fr/freefem/job/FreeFEM-sources-deployDEB/

**Output:** FreeFEM-X.X.X-amd64-ubuntuY.Y.deb

Packages:

```
sudo apt-get install cpp freeglut3-dev g++ gcc gfortran \
    m4 make patch pkg-config wget python unzip \
    liblapack-dev libhdf5-dev libgsl-dev mpich \
    autoconf automake autotools-dev bison flex gdb git cmake
```

# Windows

**Current system:** Windows 10

**Jenkins dashboard:** https://ci.inria.fr/freefem/job/deployEXE/

**Output:** FreeFEM-X.X-X-win64.exe

Windows Install:

- [MSYS2](https://www.msys2.org/)
- [git](https://git-scm.com/downloads)
- [Microsoft MPI](https://www.microsoft.com/en-us/download/details.aspx?id=100593)
- [Inno setup](https://jrsoftware.org/isdl.php)

MSYS2 packages:

```
pacman -S autoconf automake-wrapper bison bsdcpio make git mingw-w64-x86_64-toolchain mingw-w64-x86_64-freeglut patch python  flex pkg-config pkgfile tar unzip mingw-w64-x86_64-cmake mingw-w64-x86_64-msmpi mingw-w64-x86_64-gsl mingw-w64-x86_64-jq
pacman -R mingw-w64-x86_64-python mingw-w64-x86_64-gdb mingw-w64-x86_64-gdb-multiarch

```

# MacOS system

Not currently build
