[![Codacy Badge](https://api.codacy.com/project/badge/Grade/710d25bb3c6040c19c3ff7c0f3201835)](https://www.codacy.com/app/sgarnotel/FreeFem-sources?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=FreeFem/FreeFem-sources&amp;utm_campaign=Badge_Grade)
[![Build Status](https://travis-ci.org/FreeFem/FreeFem-sources.svg?branch=master)](https://travis-ci.org/FreeFem/FreeFem-sources)
<a href="https://scan.coverity.com/projects/freefem-freefem-sources">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/15433/badge.svg"/>
</a>

| Branch | Ubuntu All | Ubuntu No | MacOSX All | MacOSX No | Windows All | Windows No |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Develop | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-UbuntuAll)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-UbuntuAll/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-UbuntuNo)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-UbuntuNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-MacOSXAll)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-MacOSXAll/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-MacOSXNo)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-MacOSXNo/) |  |  |
| Master | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-UbuntuAll)](https://ci.inria.fr/freefem/job/FreeFem-source-master-UbuntuAll/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-UbuntuNo)](https://ci.inria.fr/freefem/job/FreeFem-source-master-UbuntuNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-MacOSXAll)](https://ci.inria.fr/freefem/job/FreeFem-source-master-MacOSXAll/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-MacOSXNo)](https://ci.inria.fr/freefem/job/FreeFem-source-master-MacOSXNo/) |  |  ||

See [CI/CD Tools](#ci-cd-tools)

# FreeFem++ sources

## For users

The user documentation is available [here](https://github.com/FreeFem/FreeFem-doc) in Web format, or [here](https://github.com/FreeFem/FreeFem-doc-pdf/raw/master/freefem%2B%2Bdoc.pdf) in PDF format.

## For developers

All development take place in the develop branch (or in feature branches: cmake, geneo4PETSc, ... for specific projects)

**Do not commit on master branch !**

Have a look on the [Wiki](https://github.com/FreeFem/FreeFem-sources/wiki)!

## CI/CD Tools

All: all dependency packages are installed (computer with root access).<br/>
No : dependency packages are not installed (computer without root access).

__Executed commands:__

Automatic configuration:

```bash
autoreconf -i
```

Configuration:

```bash
./configure --enable-download --enable-optim --disable-pastix --prefix=/builds/freefem-source-*
```
`*` is `master` or `develop`

Download:

```bash
./download/getall -a
```

PETSc (only on Ubuntu):

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
