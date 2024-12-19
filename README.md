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
<summary> CI </summary>

Each of the subsequent sections correspond to a workflow composed of a set of
jobs: Release/Debug + different plateform.  A workflow appeara as `failed` if
only one of its jobs has failed. To get more details about the failing
configuration, click on the corresponding badge below.

| master                                                                                                           | develop                                                                                                                         |
|:----------------------------------------------------------------------------------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------:|
| ![Minimal version](https://github.com/FreeFem/FreeFem-sources/actions/workflows/minimal.yml/badge.svg)           | ![Minimal version](https://github.com/FreeFem/FreeFem-sources/actions/workflows/minimal.yml/badge.svg?branch=develop)           |
| ![Sequential version](https://github.com/FreeFem/FreeFem-sources/actions/workflows/sequential.yml/badge.svg)     | ![Sequential version](https://github.com/FreeFem/FreeFem-sources/actions/workflows/sequential.yml/badge.svg?branch=develop)     |
| ![Full version OpenMPI](https://github.com/FreeFem/FreeFem-sources/actions/workflows/full-openmpi.yml/badge.svg) | ![Full version OpenMPI](https://github.com/FreeFem/FreeFem-sources/actions/workflows/full-openmpi.yml/badge.svg?branch=develop) |
| ![Full version MPICH](https://github.com/FreeFem/FreeFem-sources/actions/workflows/full-mpich.yml/badge.svg)     | ![Full version MPICH](https://github.com/FreeFem/FreeFem-sources/actions/workflows/full-mpich.yml/badge.svg?branch=develop)     |
| ![Full version MS-MPI](https://github.com/FreeFem/FreeFem-sources/actions/workflows/full-msmpi.yml/badge.svg)    | ![Full version MS-MPI](https://github.com/FreeFem/FreeFem-sources/actions/workflows/full-msmpi.yml/badge.svg?branch=develop)    |

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

See [Github Actions workflow files](.github/workflows)
