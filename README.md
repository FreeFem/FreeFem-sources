<!----------------------------------------------------------------------------------->
<!--- This file is part of FreeFEM.                                               --->
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

<details>
<summary> CI / CD tools </summary>

[Commit check](.github/CI.md)

[Release](.github/RELEASE.md)

Jenkins ([FreeFEM](https://ci.inria.fr/freefem/)):

_Master branch_

|                                                                                    Release                                                                                    |                                                                                 .pkg                                                                                  |                                                                                    AppImage                                                                                     |                                                                                 .deb                                                                                  |                                                                 .exe                                                                  |                                                                     Docker                                                                      |
| :---------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :-----------------------------------------------------------------------------------------------------------------------------------: | :---------------------------------------------------------------------------------------------------------------------------------------------: |
| [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFEM-sources-createRelease)](https://ci.inria.fr/freefem/view/Master/job/FreeFEM-sources-createRelease/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFEM-sources-deployPKG)](https://ci.inria.fr/freefem/view/Master/job/FreeFEM-sources-deployPKG/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFEM-sources-deployAppImage)](https://ci.inria.fr/freefem/view/Master/job/FreeFEM-sources-deployAppImage/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFEM-sources-deployDEB)](https://ci.inria.fr/freefem/view/Master/job/FreeFEM-sources-deployDEB/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=deployEXE)](https://ci.inria.fr/freefem/view/Master/job/deployEXE/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFEM-docker)](https://ci.inria.fr/freefem/view/Docker/job/FreeFEM-docker/) |

See [CI/CD Tools](#cicd-tools)

</details>

# FreeFEM sources

[FreeFEM](https://freefem.org) is a partial differential equation solver for non-linear multi-physics systems in 2D and 3D using the finite element method.

Problems involving partial differential equations from several branches of physics such as fluid-structure interactions require interpolations of data on several meshes and their manipulation within one program.

FreeFEM includes a fast interpolation algorithm and a language for the manipulation of data on multiple meshes. It is written in C++ and the FreeFEM language is a C++ idiom.

## For users

The user documentation is available [here](https://github.com/FreeFem/FreeFem-doc).

If you use FreeFEM for academic research, please use the following:

**BibTeX:**

```
@article{MR3043640,
  AUTHOR = {Hecht, F.},
  TITLE = {New development in FreeFem++},
  JOURNAL = {J. Numer. Math.},
  FJOURNAL = {Journal of Numerical Mathematics},
  VOLUME = {20}, YEAR = {2012},
  NUMBER = {3-4}, PAGES = {251--265},
  ISSN = {1570-2820},
  MRCLASS = {65Y15},
  MRNUMBER = {3043640},
  URL = {https://freefem.org/}
}
```

**APA:**

```
Hecht, F. (2012). New development in FreeFem++. Journal of numerical mathematics, 20(3-4), 251-266.
```

**ISO 690:**

```
HECHT, Frédéric. New development in FreeFem++. Journal of numerical mathematics, 2012, vol. 20, no 3-4, p. 251-266.
```

**MLA:**

```
Hecht, Frédéric. "New development in FreeFem++." Journal of numerical mathematics 20.3-4 (2012): 251-266.
```

## For developers

All development efforts take place in the _develop_ branch (or in feature branches: feature-cmake, geneo4PETSc, ... for specific projects)

**Do not commit on master branch!**

Have a look on the [Wiki](https://github.com/FreeFem/FreeFem-sources/wiki)!

## CI/CD Tools

### FreeFEM-dev

See [Jenkins configuration files](etc/jenkins)

### FreeFEM

All: all dependency packages are installed (computer with root access).<br/>
No: dependency packages are not installed (computer without root access).

<sup>1</sup>: Ubuntu 18.04 x86

<sup>2</sup>: macOS 10.13

<sup>3</sup>: macOS 10.9

<sup>4</sup>: Windows 7 + MSYS2 + MS MPI 7

**Executed commands:**

Automatic configuration:

```bash
autoreconf -i
```

Configuration:

```bash
./configure --enable-download --enable-optim
```

If you do not have administrator rights or do not want FreeFEM files scattered around on your machine, please use the `--prefix` option, e.g.:

```bash
./configure --enable-download --enable-optim --prefix=${HOME}/FreeFem-install
```

Download:

```bash
./3rdparty/getall -a
```

PETSc:

```bash
cd 3rdparty/ff-petsc
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
