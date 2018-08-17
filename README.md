<!----------------------------------------------------------------------------------->
<!--- This file is part of FreeFem++.                                             --->
<!---                                                                             --->
<!--- FreeFem++ is free software: you can redistribute it and/or modify           --->
<!--- it under the terms of the GNU Lesser General Public License as published by --->
<!--- the Free Software Foundation, either version 3 of the License, or           --->
<!--- (at your option) any later version.                                         --->
<!---                                                                             --->
<!--- FreeFem++ is distributed in the hope that it will be useful,                --->
<!--- but WITHOUT ANY WARRANTY; without even the implied warranty of              --->
<!--- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               --->
<!--- GNU Lesser General Public License for more details.                         --->
<!---                                                                             --->
<!--- You should have received a copy of the GNU Lesser General Public License    --->
<!--- along with Foobar.  If not, see <http://www.gnu.org/licenses/>.             --->
<!----------------------------------------------------------------------------------->

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/710d25bb3c6040c19c3ff7c0f3201835)](https://www.codacy.com/app/sgarnotel/FreeFem-sources?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=FreeFem/FreeFem-sources&amp;utm_campaign=Badge_Grade)
[![Build Status](https://travis-ci.org/FreeFem/FreeFem-sources.svg?branch=master)](https://travis-ci.org/FreeFem/FreeFem-sources)
<a href="https://scan.coverity.com/projects/freefem-freefem-sources">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/15433/badge.svg"/>
</a>

| Branch | Ubuntu All<sup>1</sup> | Ubuntu No<sup>1</sup> | MacOSX All<sup>2</sup> | MacOSX No<sup>3</sup> | Windows 7 64bits<sup>4</sup> | Windows7 32bits<sup>4<sup> |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Develop | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-UbuntuAll)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-UbuntuAll/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-UbuntuNo)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-UbuntuNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-MacOSXAll)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-MacOSXAll/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-MacOSXNo)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-MacOSXNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-Windows7)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-Windows7) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-Windows7-32)](https://ci.inria.fr/freefem/job/FreeFem-source-devlop-Windows7-32) |
| Master | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-UbuntuAll)](https://ci.inria.fr/freefem/job/FreeFem-source-master-UbuntuAll/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-UbuntuNo)](https://ci.inria.fr/freefem/job/FreeFem-source-master-UbuntuNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-MacOSXAll)](https://ci.inria.fr/freefem/job/FreeFem-source-master-MacOSXAll/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-MacOSXNo)](https://ci.inria.fr/freefem/job/FreeFem-source-master-MacOSXNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-Windows7)](https://ci.inria.fr/freefem/job/FreeFem-source-master-Windows7) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-Windows7-32)](https://ci.inria.fr/freefem/job/FreeFem-source-master-Windows7-32) |

See [CI/CD Tools](#ci-cd-tools)

# FreeFem++ sources

FreeFem++ is a partial differential equation solver for non-linear multi-physics systems in 2D and 3D

Problems involving partial differential equations from several branches of physics such as fluid-structure interactions require interpolations of data on several meshes and their manipulation within one program.

FreeFem++ includes a fast interpolation algorithm and a language for the manipulation of data on multiple meshes. It is written in C++ and the FreeFem++ language is a C++ idiom.

## For users

The user documentation is available [here](https://github.com/FreeFem/FreeFem-doc) in Web format, or [here](https://github.com/FreeFem/FreeFem-doc-pdf/raw/master/freefem%2B%2Bdoc.pdf) in PDF format.

## For developers

All development efforts take place in the _develop_ branch (or in feature branches: cmake, geneo4PETSc, ... for specific projects)

**Do not commit on master branch !**

Have a look on the [Wiki](https://github.com/FreeFem/FreeFem-sources/wiki)!

## CI/CD Tools

All: all dependency packages are installed (computer with root access).<br/>
No : dependency packages are not installed (computer without root access).

<sup>1</sup>: Ubuntu 16.04 x86

<sup>2</sup>: mac OS X 10.13

<sup>3</sup>: mac OS X 10.9

<sup>4</sup>: Windows 7 + MSYS2 + MS MPI 7

__Executed commands:__

Automatic configuration:

```bash
autoreconf -i
```

Configuration:

```bash
./configure --enable-download --enable-optim --disable-pastix
```

Download:

```bash
./download/getall -a
```

PETSc:

```bash
cd download/ff-petsc
sed -i 's/--download-pastix //g' Makefile
make petsc-slepc
cd -
./reconfigure
```

Make:

```bash
make -j2
make check
```

Install:

```bash
(sudo) make install
```

See [CI/CD Tools Wiki](https://github.com/FreeFem/FreeFem-sources/wiki/CI-CD-Tools) for more informations.
