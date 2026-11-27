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

# Install FreeFem++ binaries

You can download FreeFem++ binaries for your platform on the [download page](https://doc.freefem.org/introduction/download.html#download). 
Once downloaded, double-click on the `.exe` installer and follow the instructions.

By default, FreeFem++ binaries are installed in `C:\Program Files (x86)\FreeFem++`.

You can either launch FreeFem++ from the launcher icon on your desktop
(sequential version), or from a Powershell if you want to use the parallel
version:

```shell
          cd "C:\Program Files (x86)\FreeFem++"
          & "C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" -n 4 .\FreeFem++-mpi.exe .\examples\mpi\MUMPS.edp
```

# Compile FreeFem++

FreeFem++ can be built on Windows using [MSYS2](https://www.msys2.org/). Use one
MSYS2 environment consistently, do not mix packages from UCRT64, MINGW64, or 
CLANG64.

## Select an MSYS2 environment

Open the MSYS2 shell for the environment you intend to use. UCRT64 is the
reference environment for the package examples below.

Other environments, like CLANG64, have not been tested but could work with
package names that may differ.

## Install build dependencies

```bash
# Update MSYS2
pacman -Syu
# Install dependencies
pacman -S \
  mingw-w64-ucrt-x86_64-toolchain \
  mingw-w64-ucrt-x86_64-freeglut \
  mingw-w64-ucrt-x86_64-gcc-fortran \
  mingw-w64-ucrt-x86_64-autotools \
  mingw-w64-ucrt-x86_64-openblas \
  mingw-w64-ucrt-x86_64-cmake \
  make m4 patch git flex bison unzip pkg-config pkgfile tar
```

### MPI support

MS-MPI is a prerequisite for building or running FreeFem++ with MPI. It can be 
downloaded from [here](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi).

For UCRT64 builds, also install the MS-MPI package from the UCRT64 repository:

```bash
pacman -S mingw-w64-ucrt-x86_64-msmpi
```

## Configure

From the FreeFem++ source directory, in the selected MSYS2 shell:

```bash
autoreconf -i
./configure --enable-download --enable-optim --enable-generic [--prefix=<install_prefix>]
```

This will automatically download missing FreeFem++ dependencies.
If you want PETSc/SLEPc's support, run:

```bash
./3rdparty/getall -a -o PETSc
cd 3rdparty/ff-petsc
make petsc-slepc
cd -
./reconfigure
```

## Build, check and install

```bash
make -j 4
make check
make install
```

## Run FreeFem++

FreeFem++ executables and their dlls will be located in the `/ucrt64/bin` (or in
`<install_prefix>/bin` if you specified an installation prefix).

From a Powershell (note the '\' instead of '/' in paths):

```shell
cd C:\msys64\ucrt64\bin
& "C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" -n 4 .\FreeFem++-mpi.exe ..\share\FreeFEM\4.17\examples\mpi\MUMPS.edp
```
