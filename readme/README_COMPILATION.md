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

# Compilation of FreeFem++ and bamg (mesh generator) under Unix, MacOs X or MinGW (Windows)

Read the following links, depending on your system for installation instructions:

* [FreeFem-doc](https://github.com/FreeFem/FreeFem-doc) -> Documentation website -> Installation
* [http://www.freefem.org/ff++/linux.php](http://www.freefem.org/ff++/linux.php) or [web/linux.php ](https://github.com/FreeFem/FreeFem-sources/blob/master/web/linux.php)
* [http://www.freefem.org/ff++/windows.php](http://www.freefem.org/ff++/windows.php) or [web/windows.php](https://github.com/FreeFem/FreeFem-sources/blob/master/web/windows.php)
* [http://www.freefem.org/ff++/macosx.php](http://www.freefem.org/ff++/macosx.php) or web/[macosx.php](https://github.com/FreeFem/FreeFem-sources/blob/master/web/linux.php)

Now (april 2018), the PDF doc is on GitHub on https://github.com/FreeFem/FreeFem-doc-pdf/

## To use a specific configuration file

Create file `config.param` and use the shell script:
```bash
./reconfigure
```
to reconfigure your system

## Examples of `config.param` files:

### Used for the precompiled MacOS version:
```bash
$ cat config.param
'--enable-download'
'--enable-optim'
'--enable-m64'
'F77=ifort'
'FC=ifort'
'--enable-maintainer-mode'
'--with-mkl=/opt/intel/mkl/lib'
'--with-petsc=/usr/local/ff++/mpich/petsc/lib/petsc/conf/petscvariables'
'-with-hdf5=/usr/local/ff++/mpich/bin/h5cc'
'--with-gsl-prefix=/usr/local/ff++/mpich/'
'--disable-pastix'
'--with-mpipath=/usr/local/ff++/mpich/bin/'
```

### Used for the MSWindows precompiled version:
```bash
$ cat /Volumes/C/msys64/home/hecht/ff++/config.param
'-with-glut=-lfreeglut -lglu32 -lopengl32 -lwinmm -lgdi32 -Wl,--subsystem,windows'
'-without-mpi'
'-with-blas=/home/hecht/64/bin/libopenblas.dll'
'-with-lapack=/home/hecht/64/mingw/bin/libopenblas.dll'
'--disable-hips'
'--disable-pastix'
'--enable-download'
'CXX=x86_64-w64-mingw32-g++'
'FC=x86_64-w64-mingw32-gfortran'
'F77=x86_64-w64-mingw32-gfortran'
'CC=x86_64-w64-mingw32-gcc'
```

### Used for the Ubuntu version

#### without MPI
```bash
$ cat config.param
'--enable-download'
'--without-mpi'
```

#### with MPI
```bash
$ cat config.param
'--enable-download'
```
