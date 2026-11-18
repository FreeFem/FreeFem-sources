############################################################################
# This file is part of FreeFEM.                                            #
#                                                                          #
# FreeFEM is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFEM is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          #
############################################################################
# SUMMARY : Makefile for downloaded arpack-ng
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

# Downloading and compiling arpack-ng
# -----------------------------------
# arpack-ng (https://github.com/opencollab/arpack-ng) is the maintained
# fork of the original 1996 ARPACK. It keeps the same F77 interface, so
# src/Eigen/arpackff.hpp works unchanged. The library is built with its
# own CMake build system into a static libarpack.a copied to ../lib,
# where configure.ac expects it (ARPACKLIB).

# Package selection (set by configure)
# DOWNLOAD_ARPACK : "install" when the downloaded library must be built
# ARPACKLIB       : final library path (empty when using a system arpack)

# Package information
PKG_VERSION=3.9.1
DIRPKG=../pkg
PKG_SRCDIR=ARPACK
PKG_TARBALL=$(DIRPKG)/arpack-ng-$(PKG_VERSION).tar.gz
PKG_BUILD=$(PKG_SRCDIR)/build
PKG_PREFIX=$(CURDIR)/prefix
CMAKE ?= cmake
CMAKE_GENERATOR=$(if $(WIN32DLLTARGET),-G "MSYS Makefiles")

PKGCOMMON_PACKTITLE=ARPACK
PKGCOMMON_PACKAGES=$(PKG_TARBALL)
include ../common.mak

# Common package goals
all-local:$(DOWNLOAD_ARPACK)
EXTRA_DIST=
compilepkg::$(ARPACKLIB)
.NOTPARALLEL:

# Download and extract
$(PKG_SRCDIR):$(PKG_TARBALL)
	rm -rf $(PKG_SRCDIR)
	mkdir -p $(PKG_SRCDIR)
	cd $(PKG_SRCDIR) && tar xzf ../$(PKG_TARBALL) --strip-components=1
	touch $(PKG_SRCDIR)

# Configure (CMake)
# BLAS/LAPACK come from FreeFEM's configure (@BLASLIBS@/@LAPACKLIBS@);
# pass them explicitly so CMake's FindBLAS/FindLAPACK cannot pick
# different ones than the rest of the build.
$(PKG_BUILD)/configure.done:$(PKG_SRCDIR) Makefile
	mkdir -p $(PKG_BUILD)
	rm -f $(PKG_BUILD)/CMakeCache.txt
	$(CMAKE) -S $(PKG_SRCDIR) -B $(PKG_BUILD) \
	 -DCMAKE_Fortran_COMPILER='$(F77)' \
	 -DCMAKE_Fortran_FLAGS='$(FFLAGS)' \
	 -DCMAKE_C_COMPILER="$(CC)" \
	 -DCMAKE_C_FLAGS='$(CFLAGS)' \
	 -DCMAKE_INSTALL_PREFIX="$(PKG_PREFIX)" \
	 -DCMAKE_INSTALL_LIBDIR=lib \
	 -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
	 -DBUILD_SHARED_LIBS=OFF \
	 -DTESTS=OFF \
	 -DEXAMPLES=OFF \
	 -DICB=OFF \
	 -DMPI=OFF \
	 -DBLAS_LIBRARIES='$(BLASLIBS)' \
	 -DLAPACK_LIBRARIES='$(LAPACKLIBS)' \
	 $(CMAKE_GENERATOR)
	touch $@

# Build and install (CMake)
$(PKG_BUILD)/install.done:$(PKG_BUILD)/configure.done
	$(CMAKE) --build $(PKG_BUILD)
	$(CMAKE) --install $(PKG_BUILD)
	touch $@

# Install into FreeFEM's 3rdparty tree
$(ARPACKLIB):$(PKG_BUILD)/install.done
	mkdir -p ../lib
	cp $(PKG_PREFIX)/lib/libarpack.a $(ARPACKLIB)

# Clean
clean-local::
	-rm -rf $(PKG_SRCDIR) $(PKG_BUILD) $(PKG_PREFIX)
	-rm -f ../lib/libarpack.a

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
