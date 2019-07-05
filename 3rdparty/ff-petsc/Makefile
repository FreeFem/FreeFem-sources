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
PETSC_ADDCONF?=
SPELC_ADDCONF?=

PYTHON2?=python2.7
include Makefile.inc
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
dir3rdparty := $(abspath $(dir $(mkfile_path))/..)
 
W_SUDO:=$(SHELL mkdir -p "$(DESTDIR)$(FF_prefix_petsc)" && test -w "$(DESTDIR)$(FF_prefix_petsc) || echo sudo)
W_SUDO1=$(SHELL ls	)
ifeq ($(FF_generic_petsc), yes)
  FLAGS_MTUNE := "-mtune=generic"
else
  FLAGS_MTUNE := "-mtune=native"
endif
ifeq ($(WIN32DLLTARGET),)
##   Not on  windows ....
  FFWINDOWS=false
  WITH_SLEPC=yes
  FLAGS_CONF_PETSC_REAL= 		--with-debugging=no \
		--with-cc='$(MPICC)' --with-cxx='$(MPICXX)' --with-fc='$(MPIFC)' \
		COPTFLAGS='-O3  $(FLAGS_MTUNE)' \
		CXXOPTFLAGS='-O3  $(FLAGS_MTUNE)' \
		FOPTFLAGS='-O3  $(FLAGS_MTUNE)' \
		CFLAGS=' $(FLAGS_MTUNE)' \
		CXXFLAGS='$(FLAGS_MTUNE)' \
		FFLAGS='$(FLAGS_MTUNE)' \
        --with-ssl=0 --with-x=0 --with-fortran-bindings=0 \
		--with-scalar-type=real $(WITH_MPI_DIR) $(PETSC_ADDCONF) \
		--with-blas-lapack-lib='$(BLAS_LAPACK_LIBS)' \
		--download-scalapack --download-metis --download-ptscotch \
		--download-mumps --download-ml --download-hypre --download-parmetis \
		--download-superlu --download-suitesparse --download-fftw --download-slepc 

FLAGS_CONF_PETSC_COMPLEX= 	--with-debugging=no \
		--with-cc='$(MPICC)' --with-cxx='$(MPICXX)' --with-fc='$(MPIFC)' \
		COPTFLAGS='-O3  $(FLAGS_MTUNE)' \
		CXXOPTFLAGS='-O3  $(FLAGS_MTUNE)' \
		FOPTFLAGS='-O3  $(FLAGS_MTUNE)' \
		CFLAGS=' $(FLAGS_MTUNE)' \
		CXXFLAGS='$(FLAGS_MTUNE)' \
		FFLAGS='$(FLAGS_MTUNE)' \
		--with-scalar-type=complex \
		--with-ssl=0 --with-x=0 --with-fortran-bindings=0 $(WITH_MPI_DIR) $(PETSC_ADDCONF) \
		--with-blas-lapack-lib='$(BLAS_LAPACK_LIBS)' \
		--with-scalapack-dir=$(DIR_INSTALL_REAL) \
		--with-metis-dir=$(DIR_INSTALL_REAL) \
		--with-ptscotch-dir=$(DIR_INSTALL_REAL) \
		--with-suitesparse-include='$(SUITESPARSE_INCLUDE)' --with-suitesparse-lib='$(SUITESPARSE_LIB)' \
		--with-mumps-dir=$(DIR_INSTALL_REAL) \
		--with-parmetis-dir=$(DIR_INSTALL_REAL) \
		--with-fftw-dir=$(DIR_INSTALL_REAL) \
		--download-slepc  

else
##   On  windows ....
  FFWINDOWS=true
  WITH_SLEPC=no
 FLAGS_CONF_PETSC_REAL= 	--with-python-exec=/usr/bin/python2 	-with-shared-libraries=0	--with-debugging=no \
	--with-cc='$(CC)' --with-cxx='$(CXX)' --with-fc=0 \
	COPTFLAGS='-O3  $(FLAGS_MTUNE) -D__int64=long\ long' \
	CXXOPTFLAGS='-O3  $(FLAGS_MTUNE) -D__int64=long\ long' \
	CFLAGS=' $(FLAGS_MTUNE) -D__int64=long\ long' \
	CXXFLAGS='$(FLAGS_MTUNE) -D__int64=long\ long' \
	--with-mpi-lib=['/C/Program\ Files\ \(x86\)/Microsoft\ SDKs/MPI/Lib/x64/msmpi.lib'] \
	--with-mpi-include=/C/Program\ Files\ \(x86\)/Microsoft\ SDKs/MPI/Include \
	--with-metis-lib=$(dir3rdparty)/lib/parmetis/libmetis.a \
	--with-metis-include=$(dir3rdparty)/include/parmetis \
	 --download-suitesparse \
	--with-ssl=0 --with-x=0 --with-debugging=0  \
	--with-mpiexec='/C/Program\ Files/Microsoft\ MPI/Bin/mpiexec' \
	--with-blas-lapack-lib='$(BLAS_LAPACK_LIBS)' 
  #--download-slepc 
  #	--with-blas-lapack-lib='/mingw64/lib/libopenblas.a' \

  
 FLAGS_CONF_PETSC_COMPLEX= --with-python-exec=/usr/bin/python2 	--with-shared-libraries=0	--with-debugging=no \
	--with-cc='$(CC)' --with-cxx='$(CXX)' --with-fc=0 \
	COPTFLAGS='-O3  $(FLAGS_MTUNE) -D__int64=long\ long' \
	CXXOPTFLAGS='-O3  $(FLAGS_MTUNE) -D__int64=long\ long' \
	CFLAGS=' $(FLAGS_MTUNE) -D__int64=long\ long' \
	CXXFLAGS='$(FLAGS_MTUNE) -D__int64=long\ long' \
	--with-mpi-lib=['/C/Program\ Files\ \(x86\)/Microsoft\ SDKs/MPI/Lib/x64/msmpi.lib'] \
	--with-mpi-include=/C/Program\ Files\ \(x86\)/Microsoft\ SDKs/MPI/Include \
	--with-metis-lib=$(dir3rdparty)/lib/parmetis/libmetis.a \
	--with-metis-include=$(dir3rdparty)/include/parmetis \
	 --download-suitesparse \
	--with-ssl=0 --with-x=0 --with-debugging=0  \
	--with-blas-lapack-lib='$(BLAS_LAPACK_LIBS)' \
	--with-mpiexec='/C/Program\ Files/Microsoft\ MPI/Bin/mpiexec' 
#	--download-slepc 
	
#	--with-suitesparse-include=$(dir3rdparty)/include \
#	--with-suitesparse-lib='[$(dir3rdparty)/lib/libumfpack.a,$(dir3rdparty)/lib/libcholmod.a,$(dir3rdparty)/lib/libcolamd.a,$(dir3rdparty)/lib/libamd.a,$(dir3rdparty)/lib/libsuitesparseconfig.a]' \	--with-suitesparse-lib='[$(dir3rdparty)/lib/libumfpack.a,$(dir3rdparty)/lib/libcolmod.a,$(dir3rdparty)/lib/libcholamd.a,$(dir3rdparty)/lib/libamd.a,$(dir3rdparty)/lib/libsuitesparseconfig.a]' \

   
endif
all-local:
	@echo make .. 
	@echo " # To build petsc and slepc do $ make petsc-slepc"
	@echo " # warning you have no write acces in $(DESTDIR)$(FF_prefix_petsc)"
	@echo " # need sudo if '$(W_SUDO)' == 'sudo'  .." $(W_SUDO1)
	@echo " # run under root "
	@echo " # or create / add access to $(DESTDIR)$(FF_prefix_petsc) "
	@echo " # do for examples:"
	@echo "	sudo mkdir  $(DESTDIR)$(FF_prefix_petsc)"
	@echo "	sudo chown -R $(LOGNAME)  $(DESTDIR)$(FF_prefix_petsc)"
	@echo "	make petsc-slepc "
# Downloading and compiling mumps
# -------------------------------

DIRPKG=../pkg
VERSION=3.11.2
VERSION_SLEPC=3.11.0
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
WITH_MPI_DIR=

ifeq ($(MPICC),$(W_MPICC))
  WITH_MPI_DIR=
else
 # WITH_MPI_DIR=--with-mpi-dir='$(MPI_DIR)'
endif

DIR_INSTALL_REAL:=$(FF_prefix_petsc)/real
DIR_INSTALL_COMPLEX:=$(FF_prefix_petsc)/complex
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
	 echo  no write access in $(DESTDIR)$(DIR_INSTALL_REAL) and $(DESTDIR)$(DIR_INSTALL_COMPLEX); \
	 echo " do:  make  petsc-slepc SUDO=sudo "; \
	 echo " or add write access to  2 directory "; \
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
	-test -x "`type -p otool`"  &&  make changer
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fr check
	test -e $(DIR_INSTALL_REAL)/include/petsc.h
	touch $@

Makefile.inc:../../config.status Makefile Makefile-PETSc.inc
	../../config.status --file="Makefile.inc:Makefile-PETSc.inc"

$(SRCDIR)/tag-conf-real:$(SRCDIR)/tag-tar
	cd $(SRCDIR) &&  $(PYTHON2) ./configure MAKEFLAGS='' \
		--prefix=$(DIR_INSTALL_REAL) \
		$(FLAGS_CONF_PETSC_REAL) PETSC_ARCH=fr
	test -f $(SRCDIR)/fr/lib/petsc/conf/petscvariables
	touch $@
Make-petsc-download.mk:$(SRCDIR)/tag-install-real $(DIR_INSTALL_REAL)/lib/petsc/conf/petscvariables
	egrep 'PETSC_LIB_BASIC|SCALAPACK_|METIS_|MUMPS_|SUPERLU_|FFTW_|PTSCOTCH_|SUITESPARSE_'   $(DIR_INSTALL_REAL)/lib/petsc/conf/petscvariables | sed 's/-I/ /g'|sort >$@
ifdef 	COMPLEX_CASE
#  version COMPLEX  .....
$(SRCDIR)/tag-conf-complex:$(SRCDIR)/tag-tar Make-petsc-download.mk
	cd $(SRCDIR) &&  $(PYTHON2) ./configure MAKEFLAGS='' \
		--prefix=$(DIR_INSTALL_COMPLEX) \
		$(FLAGS_CONF_PETSC_COMPLEX) PETSC_ARCH=fc
	test -f $(SRCDIR)/fc/lib/petsc/conf/petscvariables
	touch $@

$(SRCDIR)/tag-make-complex:$(SRCDIR)/tag-conf-complex Make-petsc-download.mk
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fc all
	touch $@

$(SRCDIR)/tag-install-complex :$(SRCDIR)/tag-make-complex
	cd $(SRCDIR) && $(MAKE)  PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fc install
	-test -x "`type -p otool`"  &&  make changec
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
ifeq ($(WITH_SLEPC),yes)	
	test -e $(DIR_INSTALL_COMPLEX)/include/slepc.h
	-mkdir ../lib
	echo 'slepccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib" -lslepc' > $@
	echo 'slepccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"' >> $@
	echo done $@
else
	echo no slepc sorry $@
endif
##  need include include Make-petsc-download.mk
WHERE-all:Makefile WHERE WHERE-complex
	echo >$@
	test -n  "$(SCALAPACK_LIB)"  && (  echo scalapack LD) $(SCALAPACK_LIB)  >>$@ )
	test -n  "$(SCALAPACK_INCLUDE)"  && (   echo scalapack INCLUDE -I$(SCALAPACK_INCLUDE)  >>$@ )
	test -n  "$(METIS_LIB)"  && (  echo metis LD $(METIS_LIB)  >>$@ )
	test -n  "$(METIS_INCLUDE)"  && (  echo metis INCLUDE -I$(METIS_INCLUDE)  >>$@ )
	test -n  "$(MUMPS_LIB)"  && (  echo mumps LD $(MUMPS_LIB)  >>$@ )
	test -n  "$(MUMPS_INCLUDE)"  && (  echo mumps INCLUDE -I$(MUMPS_INCLUDE)  >>$@ )
	test -n  "$(SUPERLU_LIB)"  && (  echo superlu LD $(SUPERLU_LIB)  >>$@ )
	test -n  "$(SUPERLU_INCLUDE)"  && (  echo superlu INCLUDE -I$(SUPERLU_INCLUDE)  >>$@ )
	test -n  "$(FFTW_LIB)"  && (  echo fftw LD $(FFTW_LIB)  >>$@ )
	test -n  "$(FFTW_INCLUDE)"  && (  echo fftw INCLUDE -I$(FFTW_INCLUDE)   >>$@ )
	test -n  "$(PTSCOTCH_LIB)"  && (  echo ptscotch LD $(PTSCOTCH_LIB)  >>$@ )
	test -n  "$(PTSCOTCH_INCLUDE)"  && (  echo ptscotch INCLUDE -I$(PTSCOTCH_INCLUDE)  >>$@ )
	test -n  "$(PARMETIS_LIB)"  && (  echo parmetis LD $(PARMETIS_LIB) >>$@ )
	test -n  "$(PARMETIS_INCLUDE)"  && (  echo parmetis INCLUDE -I$(PARMETIS_INCLUDE) )>0>$@
	echo 'petsc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib" -lpetsc' >>$@
	echo 'petsc INCLUDE -I"$(DIR_INSTALL_REAL)/include"'  >>$@
	echo 'petsccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib" -lpetsc' >>$@
	echo 'petsccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"'  >> $@
	test -n "$(WITH_SLEPC)" = yes && (  echo 'slepc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib" -lslepc'  >>$@ )
	test -n "$(WITH_SLEPC)" = yes && (  echo 'slepc INCLUDE -I"$(DIR_INSTALL_REAL)/include"'   >>$@ )
	test -n "$(WITH_SLEPC)" = yes && (  echo 'slepccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib" -lslepc'  >>$@ )
	test -n "$(WITH_SLEPC)" = yes && (  echo 'slepccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"'  )>> $@

#
else
#
WHERE-complex:
	@if [ -n "$(MPI_INCLUDE)" ]; then $(MAKE) -f Makefile.complex WHERE-complex ;\
	else echo " -- No PETSc, no MPI"; fi

WHERE-all:Makefile  Make-petsc-download.mk
	 $(MAKE) -f Makefile.complex $@
endif
$(SRCDIR)/tag-tar:$(PACKAGE)
	tar xvzf $(PACKAGE)
	patch -p1 < petsc-3.11.2.patch
	touch $@
$(PACKAGE):
	../getall -o PETSc -a

changec:
	diri=$(DIR_INSTALL_COMPLEX)/lib ; \
	libpp=$$diri/libslepc.$(VERSION_SLEPC).dylib ; \
    opblibpp=`otool -L $$libpp  | awk '/3rdparty/ {print $$1}'` ; \
	npblibpp=$$diri/`basename $$opblibpp` ; \
	echo -- $$opblibpp ; \
	echo -- $$npblibpp ; \
	test -n "$$opblibpp" &&  install_name_tool -change $$opblibpp $$npblibpp $$libpp ; \
	otool -L $$libpp  | awk '/libpetsc/ {print "**",$$1}'
changer:
	diri=$(DIR_INSTALL_REAL)/lib ; \
	libpp=$$diri/libslepc.$(VERSION_SLEPC).dylib ; \
    opblibpp=`otool -L $$libpp  | awk '/3rdparty/ {print $$1}'` ; \
	npblibpp=$$diri/`basename $$opblibpp` ; \
	echo -- $$opblibpp ; \
	echo -- $$npblibpp ; \
	test -n "$$opblibpp" &&  install_name_tool -change $$opblibpp $$npblibpp $$libpp ; \
	otool -L $$libpp  | awk '/libpetsc/ {print "**",$$1}'

../lib/WHERE.SLEPc:
ifeq ($(WITH_SLEPC),yes)	
	-mkdir ../lib
	test -e $(DIR_INSTALL_REAL)/include/slepc.h
	echo 'slepc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib" -lslepc' > $@
	echo 'slepc INCLUDE -I"$(DIR_INSTALL_REAL)/include"' >> $@
else
	echo "sorry no SLEPC" 
endif

clean:clean-local
clean-local:
	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR)
	-rm Makefile.inc FAIRE* ../lib/WHERE.PETSc* ../lib/WHERE.SLEPc*
	-rm -rf ../include/*PETSc*
	-rm -rf ../lib/lib*PETSc*
	-rm -rf $(SRCDIR)
	-rm -rf $(FF_prefix_petsc)
	-rm WHERE-all config.log *.done
	-if test -d $(FF_prefix_petsc) ; then echo " try of  remove of $(FF_prefix_petsc) under sudo .."; sudo rm -rf $(FF_prefix_petsc) ; fi
	-rm do-sudo
echo:
	@echo " with_dir_mpi #$(W_MPICC)# #$(MPICC)# ::: #$(WITH_MPI_DIR)#  "
	@echo SUDO: $(SUDO)
	@echo MTUNE: $(FLAGS_MTUNE)
	@echo "sudo using user env: $(SUDO_WITH_ENV)"
	@echo MPI_DIR: $(MPI_DIR)
	@echo " dir install real :" $(DIR_INSTALL_REAL)
	@echo " dir install complex :" $(DIR_INSTALL_COMPLEX)
	@echo " do-sudo auto ???:" $(SHELL  cat do-sudo)
	@echo " BLAS_LAPACK_LIBS: $(BLAS_LAPACK_LIBS)"
	@echo " # need sudo if '$(W_SUDO)' == 'sudo'  .." $(W_SUDO1)
	@echo " on WINDOWS: $(FFWINDOWS)"
	@echo " dir3rdparty: $(dir3rdparty)"
	-@otool -L $(DIR_INSTALL_COMPLEX)/lib/libslepc.$(VERSION_SLEPC).dylib  | awk '/3rdparty/ {print $$1}'
	-@otool -L $(DIR_INSTALL_REAL)/lib/libslepc.$(VERSION_SLEPC).dylib  | awk '/3rdparty/ {print $$1}'

install-destdir:
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fr install DESTDIR=$(DESTDIR)
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=fc install DESTDIR=$(DESTDIR)

.NOTPARALLEL:
# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
