# ======================================================================
# Laboratoire Jacques-Louis Lions
# Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
# http://www.ljll.math.upmc.fr/lehyaric
# ======================================================================
# This file is part of Freefem++
#
# Freefem++ is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
#
# Freefem++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with Freefem++; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
# ======================================================================
# headeralh default=0 freefem make multipleauthors start=19/03/10 upmc

include Makefile.inc
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
dir3rdparty := $(abspath $(dir $(mkfile_path))/..)

DIR_INSTALL_REAL:=$(FF_prefix_petsc)/r
DIR_INSTALL_COMPLEX:=$(FF_prefix_petsc)/c

W_SUDO:=$(SHELL mkdir -p "$(DESTDIR)$(FF_prefix_petsc)" && test -w "$(DESTDIR)$(FF_prefix_petsc) || echo sudo)
ifeq ($(COMPILE_OPENBLAS),openblas)
PETSC_OPENBLAS := --download-f2cblaslapack
PETSC_OPENBLAS_C := --with-blaslapack-dir=$(DIR_INSTALL_REAL)
else
PETSC_OPENBLAS = --with-blaslapack-include='$(shell echo $(BLASINC) | sed "s/-I//g")' \
		--with-blaslapack-lib='$(BLAS_LAPACK_LIBS)'
PETSC_OPENBLAS_C = $(PETSC_OPENBLAS)
endif
PETSC_DOWNLOAD_C := --with-metis-dir=$(DIR_INSTALL_REAL) \
		--with-ptscotch-dir=$(DIR_INSTALL_REAL) \
		--with-suitesparse-include='$(SUITESPARSE_INCLUDE)' --with-suitesparse-lib='$(SUITESPARSE_LIB)' \
		--with-parmetis-dir=$(DIR_INSTALL_REAL) \
		--with-tetgen-dir=$(DIR_INSTALL_REAL)
ifeq ($(PETSC_PKG),yes)
pkg_dir=$(abspath $(dir $(mkfile_path))/pkg)
PETSC_DOWNLOAD := --download-scalapack=$(pkg_dir)/pkg-scalapack.tar.gz \
		--download-metis=$(pkg_dir)/pkg-metis.tar.gz \
		--download-ptscotch=$(pkg_dir)/scotch-v6.0.9.tar.gz \
		--download-mumps=$(pkg_dir)/pkg-mumps.tar.gz \
		--download-hypre=$(pkg_dir)/hypre.tar.gz \
		--download-parmetis=$(pkg_dir)/pkg-parmetis.tar.gz \
		--download-superlu=$(pkg_dir)/superlu.tar.gz \
		--download-suitesparse=$(pkg_dir)/SuiteSparse-5.7.1.tar.gz \
		--download-tetgen=$(pkg_dir)/tetgen1.5.1.tar.gz \
		--download-slepc=$(pkg_dir)/slepc.tar.gz \
		--download-hpddm=$(pkg_dir)/hpddm.tar.gz
ifeq ($(FFCMAKE),no)
PETSC_DOWNLOAD := $(PETSC_DOWNLOAD) --download-cmake=$(pkg_dir)/cmake.tar.gz
endif
PETSC_DOWNLOAD_C := $(PETSC_DOWNLOAD_C) --download-slepc=$(pkg_dir)/slepc.tar.gz
ifeq ($(WIN32DLLTARGET),)
PETSC_DOWNLOAD := $(PETSC_DOWNLOAD) --download-hpddm=$(pkg_dir)/hpddm.tar.gz
PETSC_DOWNLOAD_C := $(PETSC_DOWNLOAD_C) --download-hpddm=$(pkg_dir)/hpddm.tar.gz
endif
else
PETSC_DOWNLOAD := --download-metis --download-ptscotch \
		--download-hypre --download-parmetis \
		--download-superlu --download-suitesparse --download-tetgen \
		--download-slepc
# --download-slepc-commit=d3f7cd4a
ifeq ($(FFCMAKE),no)
PETSC_DOWNLOAD := $(PETSC_DOWNLOAD) --download-cmake
endif
PETSC_DOWNLOAD_C := $(PETSC_DOWNLOAD_C) --download-slepc
# --download-slepc-commit=d3f7cd4a
ifeq ($(WIN32DLLTARGET),)
PETSC_DOWNLOAD := $(PETSC_DOWNLOAD) --download-hpddm --download-hpddm-commit=ce6ce80
PETSC_DOWNLOAD_C := $(PETSC_DOWNLOAD_C) --download-hpddm --download-hpddm-commit=ce6ce80
endif
endif
ifeq ($(FF_generic_petsc), yes)
FLAGS_MTUNE := -mtune=generic
else
FLAGS_MTUNE := -mtune=native
endif
COMMON_FLAGS := MAKEFLAGS='' --with-debugging=0 COPTFLAGS='-O3 $(FLAGS_MTUNE)' CXXOPTFLAGS='-O3 $(FLAGS_MTUNE)' FOPTFLAGS='-O3 $(FLAGS_MTUNE)' --with-cxx-dialect=C++11 --with-ssl=0 --with-x=0 --with-fortran-bindings=0 --with-cudac=0
ifeq ($(WIN32DLLTARGET),)
##  Not on windows....
ifeq ($(MPICC)$(MPICXX)$(MPIFC),)
PETSC_DOWNLOAD += --with-cc='$(CC)' --with-cxx='$(CXX)' --download-mpich
PETSC_DOWNLOAD_C += --with-mpi-dir=$(DIR_INSTALL_REAL)
ifneq ($(FC),)
PETSC_DOWNLOAD += --with-fc='$(FC)' --download-scalapack --download-mumps
PETSC_DOWNLOAD_C += --with-scalapack-dir=$(DIR_INSTALL_REAL) --with-mumps-dir=$(DIR_INSTALL_REAL)
else
COMMON_FLAGS += --with-fc=0
endif
else
COMMON_FLAGS += --with-cc='$(MPICC)' --with-cxx='$(MPICXX)'
ifneq ($(MPIFC),)
COMMON_FLAGS += --with-fc='$(MPIFC)'
PETSC_DOWNLOAD += --download-scalapack --download-mumps
PETSC_DOWNLOAD_C += --with-scalapack-dir=$(DIR_INSTALL_REAL) --with-mumps-dir=$(DIR_INSTALL_REAL)
else
COMMON_FLAGS += --with-fc=0
endif
endif
else
##   On windows....
COMMON_FLAGS += --with-shared-libraries=0 \
	--with-cc='$(CC)' --with-cxx='$(CXX)' \
	CXXFLAGS='-fno-stack-protector' \
	CFLAGS='-fno-stack-protector' \
	--with-mpi-lib='$(MPI_LIB)' \
	--with-mpi-include='$(MPI_INC_DIR)' \
	--with-mpiexec='/C/Program\ Files/Microsoft\ MPI/Bin/mpiexec'
PETSC_DOWNLOAD += '--download-metis-cmake-arguments=-G "MSYS Makefiles"' \
	'--download-parmetis-cmake-arguments=-G "MSYS Makefiles"' \
	'--download-superlu-cmake-arguments=-G "MSYS Makefiles"' \
	'--download-hypre-configure-arguments=--build=x86_64-linux-gnu --host=x86_64-linux-gnu'
ifneq ($(FC),)
COMMON_FLAGS += --with-fc='$(FC)'
PETSC_DOWNLOAD += --download-scalapack --download-mumps
PETSC_DOWNLOAD_C += --with-scalapack-dir=$(DIR_INSTALL_REAL) --with-mumps-dir=$(DIR_INSTALL_REAL)
else
COMMON_FLAGS += --with-fc=0
endif
endif

FLAGS_CONF_PETSC_REAL = $(COMMON_FLAGS) --with-scalar-type=real $(PETSC_OPENBLAS) $(PETSC_DOWNLOAD)
FLAGS_CONF_PETSC_COMPLEX = $(COMMON_FLAGS) --with-scalar-type=complex $(PETSC_OPENBLAS_C) $(PETSC_DOWNLOAD_C)

all-local:
	@echo make ..
	@echo " # To build petsc and slepc do $ make petsc-slepc"
	@echo " # warning you have no write acces in $(DESTDIR)$(FF_prefix_petsc)"
	@echo " # need sudo if '$(W_SUDO)' == 'sudo' .."
	@echo " # run under root "
	@echo " # or create / add access to $(DESTDIR)$(FF_prefix_petsc) "
	@echo " # do for examples:"
	@echo " # PETSC_DOWNLOAD = $(PETSC_DOWNLOAD)"
	@echo "	sudo mkdir $(DESTDIR)$(FF_prefix_petsc)"
	@echo "	sudo chown -R $(LOGNAME) $(DESTDIR)$(FF_prefix_petsc)"
	@echo "	make petsc-slepc "

DIRPKG=../pkg
VERSION=3.13.4
VERSION_SLEPC=3.13.3
PACKAGE=../pkg/petsc-lite-$(VERSION).tar.gz
SRCDIR=petsc-$(VERSION)
MPI_DIR:=$(shell dirname $(MPI_INC_DIR))
W_MPICC:=$(strip $(shell which mpicc))
SUDO_WITH_ENV=$(if $(SUDO), $(SUDO) -E)
ifeq ($(LAPACKLIBS), $(BLASLIBS))
  BLAS_LAPACK_LIBS=$(LAPACKLIBS)
else
  BLAS_LAPACK_LIBS=$(LAPACKLIBS) $(BLASLIBS)
endif

PETSC_DIR:=$(shell pwd)/$(SRCDIR)
petsc-slepc: Makefile.inc dir-install WHERE-all
recompile: recompile-real recompile-complex
recompile-real:
	-rm */tag-make-real
	-rm */tag-install-real
	make WHERE
recompile-complex:
	-rm */tag-make-complex
	-rm */tag-install-complex
	@if [ -n "$(MPI_INCLUDE)" ]; then $(MAKE) -f Makefile.complex WHERE-complex ;\
	else echo " -- No PETSc, no MPI"; fi

dir-install: Makefile.inc
	-$(SUDO) mkdir -p $(DESTDIR)$(DIR_INSTALL_REAL)
	-$(SUDO) chown ${USER} $(DESTDIR)$(DIR_INSTALL_REAL)
	-$(SUDO) mkdir -p $(DESTDIR)$(DIR_INSTALL_COMPLEX)
	-$(SUDO) chown ${USER} $(DESTDIR)$(DIR_INSTALL_COMPLEX)
	@if test -w 	"$(DESTDIR)$(DIR_INSTALL_REAL)" -a -w 	"$(DESTDIR)$(DIR_INSTALL_COMPLEX)" ; then echo "no need of sudo"; echo > do-sudo; \
	else \
	 echo no write access in $(DESTDIR)$(DIR_INSTALL_REAL) and $(DESTDIR)$(DIR_INSTALL_COMPLEX); \
	 echo " do: make petsc-slepc SUDO=sudo "; \
	 echo " or add write access to 2 directory "; \
	 echo sudo >do-sudo; \
	fi
	test -e "$(DESTDIR)$(DIR_INSTALL_REAL)" -a -e "$(DESTDIR)$(DIR_INSTALL_COMPLEX)"
PETSc-real:
	@if [ -n "$(MPI_INCLUDE)" ]; then $(MAKE) WHERE ;\
	else echo " -- No PETSc, no MPI"; fi

PETSc-complex:
	@if [ -n "$(MPI_INCLUDE)" ]; then $(MAKE) -f Makefile.complex WHERE-complex ;\
	else echo " -- No PETSc, no MPI, no PETSc real"; fi

WHERE:../lib/WHERE.PETSc ../lib/WHERE.SLEPc

../lib/WHERE.PETSc:$(SRCDIR)/tag-install-real
	-mkdir ../lib
	echo 'petsc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib" -lpetsc' > $@
	echo 'petsc INCLUDE -I"$(DIR_INSTALL_REAL)/include"' >> $@

$(SRCDIR)/tag-make-real:$(SRCDIR)/tag-conf-real
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fr all
	touch $@


$(SRCDIR)/tag-install-real :$(SRCDIR)/tag-make-real
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fr install
	-test -x "`type -p otool`" && make changer
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fr check
	test -e $(DIR_INSTALL_REAL)/include/petsc.h
	test -e $(DIR_INSTALL_REAL)/lib/petsc/conf/petscvariables
	touch $@

Makefile.inc:../../config.status Makefile Makefile-PETSc.inc
	../../config.status --file="Makefile.inc:Makefile-PETSc.inc"

$(SRCDIR)/tag-conf-real:$(SRCDIR)/tag-tar
	cd $(SRCDIR) && ./configure --prefix=$(DIR_INSTALL_REAL) \
		$(FLAGS_CONF_PETSC_REAL) PETSC_ARCH=fr
	test -f $(SRCDIR)/fr/lib/petsc/conf/petscvariables
	touch $@
Make-petsc-download.mk:$(SRCDIR)/tag-install-real 
	test -e $(DIR_INSTALL_REAL)/lib/petsc/conf/petscvariables
	egrep 'PETSC_LIB_BASIC|SCALAPACK_|METIS_|MUMPS_|HPDDM_|TETGEN_|SUPERLU_|PTSCOTCH_|SUITESPARSE_' $(DIR_INSTALL_REAL)/lib/petsc/conf/petscvariables | sed 's/-I/ /g'|sort >$@
ifdef 	COMPLEX_CASE
#  version COMPLEX.....
$(SRCDIR)/tag-conf-complex:$(SRCDIR)/tag-tar Make-petsc-download.mk
	cd $(SRCDIR) && ./configure --prefix=$(DIR_INSTALL_COMPLEX) \
		$(FLAGS_CONF_PETSC_COMPLEX) PETSC_ARCH=fc
	test -f $(SRCDIR)/fc/lib/petsc/conf/petscvariables
	touch $@

$(SRCDIR)/tag-make-complex:$(SRCDIR)/tag-conf-complex Make-petsc-download.mk
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fc all
	touch $@

$(SRCDIR)/tag-install-complex :$(SRCDIR)/tag-make-complex
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fc install
	-test -x "`type -p otool`" && make changec
	test -e $(DIR_INSTALL_COMPLEX)/include/petsc.h
	touch $@
WHERE-complex:../lib/WHERE.PETSc-complex ../lib/WHERE.SLEPc-complex

../lib/WHERE.PETSc-complex:$(SRCDIR)/tag-install-complex
	test -e $(DIR_INSTALL_COMPLEX)/include/petsc.h
	-mkdir ../lib
	echo 'petsccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib" -lpetsc' > $@
	echo 'petsccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"' >> $@
	echo done $@
../lib/WHERE.SLEPc-complex:
	test -e $(DIR_INSTALL_COMPLEX)/include/slepc.h
	-mkdir ../lib
	echo 'slepccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib" -lslepc' > $@
	echo 'slepccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"' >> $@
	echo done $@

WHERE-all:Makefile WHERE WHERE-complex
	echo >$@
	test -z '$(SCALAPACK_LIB)' || echo scalapack LD $(SCALAPACK_LIB) >>$@
	test -z '$(SCALAPACK_INCLUDE)' || echo scalapack INCLUDE -I$(SCALAPACK_INCLUDE) >>$@
	test -z '$(METIS_LIB)' || echo metis LD $(METIS_LIB) >>$@
	test -z '$(METIS_INCLUDE)' || echo metis INCLUDE -I$(METIS_INCLUDE) >>$@
	test -z '$(MUMPS_LIB)' || echo mumps LD $(MUMPS_LIB) >>$@
	test -z '$(MUMPS_INCLUDE)' || echo mumps INCLUDE -I$(MUMPS_INCLUDE) >>$@
	test -z '$(SUPERLU_LIB)' || echo superlu LD $(SUPERLU_LIB) >>$@
	test -z '$(SUPERLU_INCLUDE)' || echo superlu INCLUDE -I$(SUPERLU_INCLUDE) >>$@
	test -z '$(PTSCOTCH_LIB)' || echo ptscotch LD $(PTSCOTCH_LIB) >>$@
	test -z '$(PTSCOTCH_INCLUDE)' || echo ptscotch INCLUDE -I$(PTSCOTCH_INCLUDE) >>$@
	test -z '$(TETGEN_LIB)' || echo tetgen LD $(TETGEN_LIB) >>$@
	test -z '$(TETGEN_INCLUDE)' || echo tetgen INCLUDE -I$(TETGEN_INCLUDE) >>$@
	test -z '$(PARMETIS_LIB)' || echo parmetis LD $(PARMETIS_LIB) >>$@
	test -z '$(PARMETIS_INCLUDE)' || echo parmetis INCLUDE -I$(PARMETIS_INCLUDE) >>$@
	test -z '$(HPDDM_LIB)' || echo hpddm LD $(HPDDM_LIB) >>$@
	test -z '$(HPDDM_INCLUDE)' || echo hpddm INCLUDE -I$(HPDDM_INCLUDE) >>$@
	echo 'petsc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib" -lpetsc' >>$@
	echo 'petsc INCLUDE -I"$(DIR_INSTALL_REAL)/include"' >>$@
	echo 'petsccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib" -lpetsc' >>$@
	echo 'petsccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"' >> $@
	echo 'slepc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib" -lslepc' >>$@
	echo 'slepc INCLUDE -I"$(DIR_INSTALL_REAL)/include"' >>$@
	echo 'slepccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib" -lslepc' >>$@
	echo 'slepccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"' >>$@

#
else
#
WHERE-complex:
	@if [ -n "$(MPI_INCLUDE)" ]; then $(MAKE) -f Makefile.complex WHERE-complex ;\
	else echo " -- No PETSc, no MPI"; fi

WHERE-all:Makefile Make-petsc-download.mk
	 $(MAKE) -f Makefile.complex $@
endif
$(SRCDIR)/tag-tar:$(PACKAGE)
	-tar xzf $(PACKAGE)
	patch -p1 < petsc-hpddm.patch
ifeq ($(WIN32DLLTARGET),)
	patch -p1 < petsc-metis.patch
endif
	touch $@
$(PACKAGE):
	../getall -o PETSc -a

changec:
	diri=$(DIR_INSTALL_COMPLEX)/lib ; \
	libpp=$$diri/libslepc.$(VERSION_SLEPC).dylib ; \
	opblibpp=`otool -L $$libpp | awk '/3rdparty/ {print $$1}'` ; \
	npblibpp=$$diri/`basename $$opblibpp` ; \
	echo -- $$opblibpp ; \
	echo -- $$npblibpp ; \
	test -n "$$opblibpp" && install_name_tool -change $$opblibpp $$npblibpp $$libpp ; \
	otool -L $$libpp | awk '/libpetsc/ {print "**",$$1}'
changer:
	diri=$(DIR_INSTALL_REAL)/lib ; \
	libpp=$$diri/libslepc.$(VERSION_SLEPC).dylib ; \
	opblibpp=`otool -L $$libpp | awk '/3rdparty/ {print $$1}'` ; \
	npblibpp=$$diri/`basename $$opblibpp` ; \
	echo -- $$opblibpp ; \
	echo -- $$npblibpp ; \
	test -n "$$opblibpp" && install_name_tool -change $$opblibpp $$npblibpp $$libpp ; \
	otool -L $$libpp | awk '/libpetsc/ {print "**",$$1}'

../lib/WHERE.SLEPc:
	-mkdir ../lib
	test -e $(DIR_INSTALL_REAL)/include/slepc.h
	echo 'slepc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib" -lslepc' > $@
	echo 'slepc INCLUDE -I"$(DIR_INSTALL_REAL)/include"' >> $@

clean:clean-local
clean-local:
	-cd $(SRCDIR) && $(MAKE) clean -C $(SRCDIR)
	-rm Makefile.inc FAIRE* ../lib/WHERE.PETSc* ../lib/WHERE.SLEPc*
	-rm -rf ../include/*PETSc*
	-rm -rf ../lib/lib*PETSc*
	-rm -rf $(SRCDIR)
	-rm -rf $(FF_prefix_petsc)
	-rm WHERE-all config.log *.done
	-if test -d $(FF_prefix_petsc) ; then echo " try of remove of $(FF_prefix_petsc) under sudo .."; sudo rm -rf $(FF_prefix_petsc) ; fi
	-rm do-sudo
echo:
	@echo SUDO: $(SUDO)
	@echo MTUNE: $(FLAGS_MTUNE)
	@echo "sudo using user env: $(SUDO_WITH_ENV)"
	@echo MPI_DIR: $(MPI_DIR)
	@echo " dir install real :" $(DIR_INSTALL_REAL)
	@echo " dir install complex :" $(DIR_INSTALL_COMPLEX)
	@echo " do-sudo auto ???:" $(SHELL cat do-sudo)
	@echo " BLAS_LAPACK_LIBS: $(BLAS_LAPACK_LIBS)"
	@echo " # need sudo if '$(W_SUDO)' == 'sudo' .."
	@echo " dir3rdparty: $(dir3rdparty)"
	-@otool -L $(DIR_INSTALL_COMPLEX)/lib/libslepc.$(VERSION_SLEPC).dylib | awk '/3rdparty/ {print $$1}'
	-@otool -L $(DIR_INSTALL_REAL)/lib/libslepc.$(VERSION_SLEPC).dylib | awk '/3rdparty/ {print $$1}'

install-destdir:
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fr install DESTDIR=$(DESTDIR)
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fc install DESTDIR=$(DESTDIR)

.NOTPARALLEL:
# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
