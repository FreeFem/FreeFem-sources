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
W_SUDO:=$(SHELL mkdir -p "$(DESTDIR)$(prefix)/ff-petsc" && test -w "$(DESTDIR)$(prefix)/ff-petsc" || echo sudo)
W_SUDO1=$(SHELL ls	)

all-local:
	@echo " # To build petsc and slepc do $ make petsc-slepc"
	@echo " # warning you have no write acces in $(DESTDIR)$(prefix)/ff-petsc"
	@echo " # need sudo if '$(W_SUDO)' == 'sudo'  .." $(W_SUDO1)
	@echo " # run under root "
	@echo " # or create / add access to $(DESTDIR)$(prefix)/ff-petsc "
	@echo " # do for examples:"
	@echo "	sudo mkdir  $(DESTDIR)$(prefix)/ff-petsc"
	@echo "	sudo chown -R $(LOGNAME)  $(DESTDIR)$(prefix)/ff-petsc"
	@echo "	make petsc-slepc "
# Downloading and compiling mumps
# -------------------------------

DIRPKG=../pkg
VERSION=3.10.3
VERSION_SLEPC=3.10.1
# http://slepc.upv.es/download/download.php?filename=slepc-3.7.4.tar.gz
PACKAGE=../pkg/petsc-lite-$(VERSION).tar.gz
PACKAGE_SLEPC=../pkg/slepc-$(VERSION_SLEPC).tar.gz
SRCDIR=petsc-$(VERSION)
SRCDIR_SLEPC=slepc-$(VERSION_SLEPC)
MPI_DIR:=$(shell dirname $(MPI_INC_DIR))
W_MPICC:=$(strip $(shell which mpicc))
SUDO_WITH_ENV=$(if $(SUDO), $(SUDO) -E)
ifeq ($(LAPACKLIBS), $(BLASLIBS))
  BLAS_LAPACK_LIBS=$(LAPACKLIBS)
else
  BLAS_LAPACK_LIBS=$(LAPACKLIBS) $(BLASLIBS)
endif
ifeq ($(MPICC),$(W_MPICC))
  WITH_MPI_DIR=
else
  WITH_MPI_DIR=--with-mpi-dir='$(MPI_DIR)'
endif
DIR_INSTALL_REAL:=$(prefix)/ff-petsc/real
DIR_INSTALL_COMPLEX:=$(prefix)/ff-petsc/complex
PETSC_DIR:=$(PWD)/$(SRCDIR)
SLEPC_DIR:=$(PWD)/$(SRCDIR_SLEPC)
petsc-slepc: Makefile.inc dir-install WHERE-all
recompile: recompile-real recompile-complex
recompile-real:
	-rm */tag-make-real
	-rm */tag-install-real
	make WHERE
recompile-complex:
	-rm */tag-make-complex
	-rm */tag-install-complex
	@if [  -n "$(MPI_INCLUDE)" ]; then $(MAKE) -f Makefile.complex WHERE-complex ;\
	else echo " -- No  Petc , no MPI"; fi

dir-install: Makefile.inc
	-$(SUDO) mkdir -p $(DESTDIR)$(DIR_INSTALL_REAL)
	-$(SUDO) mkdir -p $(DESTDIR)$(DIR_INSTALL_COMPLEX)
	@if test -w 	"$(DESTDIR)$(DIR_INSTALL_REAL)" -a -w 	"$(DESTDIR)$(DIR_INSTALL_COMPLEX)" ; then echo "no need of sudo"; echo > do-sudo; \
	else \
	 echo  no write access in $(DESTDIR)$(DIR_INSTALL_REAL) and $(DESTDIR)$(DIR_INSTALL_COMPLEX); \
	 echo " do:  make  petsc-slepc SUDO=sudo "; \
	 echo " or add write access to  2 directory "; \
	 echo sudo >do-sudo; \
	fi
	test -e "$(DESTDIR)$(DIR_INSTALL_REAL)" -a -e "$(DESTDIR)$(DIR_INSTALL_COMPLEX)"
PETSc-real:
	@if [  -n "$(MPI_INCLUDE)" ]; then $(MAKE) WHERE ;\
	else echo " -- No  Petc , no MPI"; fi

PETSc-complex:
	@if [  -n "$(MPI_INCLUDE)" ]; then $(MAKE) -f Makefile.complex WHERE-complex ;\
	else echo " -- No  Petc , no MPI, no Petsc real"; fi

WHERE:../lib/WHERE.PETSc ../lib/WHERE.SLEPc

../lib/WHERE.PETSc:$(SRCDIR)/tag-install-real
	-mkdir ../lib
	echo 'petsc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib"   -lpetsc' > $@
	echo 'petsc INCLUDE -I"$(DIR_INSTALL_REAL)/include"'  >> $@

SLEPc-real:
	@if [  -n "$(MPI_INCLUDE)" ]; then $(MAKE) ../lib/WHERE.SLEPc ;\
	else echo " -- No  Petc , no MPI"; fi
SLEPc-complex:
	@if [  -n "$(MPI_INCLUDE)" ]; then $(MAKE) -f Makefile.complex ../lib/WHERE.SLEPc-complex ;\
	else echo " -- No  Petc , no MPI"; fi

$(SRCDIR)/tag-make-real:$(SRCDIR)/tag-conf-real
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-real all
	touch $@


$(SRCDIR)/tag-install-real :$(SRCDIR)/tag-make-real
	cd $(SRCDIR) && $(SUDO) $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-real install
	cd $(SRCDIR) && $(SUDO) $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-real check
	test -e $(DIR_INSTALL_REAL)/include/petsc.h
	touch $@

Makefile.inc:../../config.status Makefile Makefile-PETSc.inc
	../../config.status --file="Makefile.inc:Makefile-PETSc.inc"

$(SRCDIR)/tag-conf-real:$(SRCDIR)/tag-tar
	cd $(SRCDIR) &&  $(PYTHON2) ./configure \
		--prefix=$(DIR_INSTALL_REAL) \
		--with-debugging=no \
		--with-cc='$(CC)' --with-cxx='$(CXX)' --with-fc='$(FC)' \
		COPTFLAGS='-O3 -march=native -mtune=native' \
		CXXOPTFLAGS='-O3 -march=native -mtune=native' \
		FOPTFLAGS='-O3 -march=native -mtune=native' \
		--with-ssl=0 ---with-x=0 --with-fortran-bindings=0 \
	  '--with-scalar-type=real' \
	  $(WITH_MPI_DIR) $(PETSC_ADDCONF) \
		--with-blas-lapack-lib='$(BLAS_LAPACK_LIBS)' \
		--download-scalapack --download-metis --download-ptscotch \
		--download-mumps --download-ml --download-hypre --download-parmetis \
		--download-superlu --download-suitesparse --download-fftw \
		PETSC_ARCH=ff-real
	test -f $(SRCDIR)/ff-real/lib/petsc/conf/petscvariables
	touch $@
Make-petsc-download.mk:$(SRCDIR)/tag-install-real $(DIR_INSTALL_REAL)/lib/petsc/conf/petscvariables
	egrep 'SCALAPACK_|METIS_|MUMPS_|SUPERLU_|FFTW_|PTSCOTCH_|SUITESPARSE_'   $(DIR_INSTALL_REAL)/lib/petsc/conf/petscvariables | sed 's/-I/ /g'|sort >$@

ifdef 	MUMPS_LIB
#  version COMPLEX  .....
$(SRCDIR)/tag-conf-complex:$(SRCDIR)/tag-tar Make-petsc-download.mk
	cd $(SRCDIR) &&  $(PYTHON2) ./configure \
		--prefix=''$(DIR_INSTALL_COMPLEX)'' \
		--with-debugging=no \
		--with-cc='$(CC)' --with-cxx='$(CXX)' --with-fc='$(FC)' \		
		COPTFLAGS='-O3 -march=native -mtune=native' \
		CXXOPTFLAGS='-O3 -march=native -mtune=native' \
		FOPTFLAGS='-O3 -march=native -mtune=native' \
	  '--with-scalar-type=complex' \
		--with-ssl=0 ---with-x=0 --with-fortran-bindings=0 \
	  $(WITH_MPI_DIR) $(PETSC_ADDCONF)  \
		--with-blas-lapack-lib='$(BLAS_LAPACK_LIBS)' \
		--with-scalapack-include='$(SCALAPACK_INCLUDE)'  --with-scalapack-lib='$(SCALAPACK_LIB)' \
		--with-metis-include='$(METIS_INCLUDE)'  --with-metis-lib='$(METIS_LIB)' \
		--with-ptscotch-include='$(PTSCOTCH_INCLUDE)'  --with-ptscotch-lib='$(PTSCOTCH_LIB)' \
		--with-suitesparse-include='$(SUITESPARSE_INCLUDE)'  --with-suitesparse-lib='$(SUITESPARSE_LIB)' \
		--with-mumps-include='$(MUMPS_INCLUDE)'  --with-mumps-lib='$(MUMPS_LIB)' \
		--with-parmetis-include='$(PARMETIS_INCLUDE)'  --with-parmetis-lib='$(PARMETIS_LIB)' \
		--with-fftw-include='$(FFTW_INCLUDE)'  --with-fftw-lib='$(FFTW_LIB)' \
		PETSC_ARCH=ff-complex
	test -f $(SRCDIR)/ff-complex/lib/petsc/conf/petscvariables
	touch $@

$(SRCDIR)/tag-make-complex:$(SRCDIR)/tag-conf-complex Make-petsc-download.mk
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-complex all
	touch $@

$(SRCDIR)/tag-install-complex :$(SRCDIR)/tag-make-complex
	cd $(SRCDIR) && $(SUDO) $(MAKE)  PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-complex install
	test -e $(DIR_INSTALL_COMPLEX)/include/petsc.h
	touch $@
WHERE-complex:../lib/WHERE.PETSc-complex ../lib/WHERE.SLEPc-complex

../lib/WHERE.PETSc-complex:$(SRCDIR)/tag-install-complex
	test -e $(DIR_INSTALL_COMPLEX)/include/petsc.h
	-mkdir ../lib
	echo 'petsccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib"   -lpetsc' > $@
	echo 'petsccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"'  >> $@
	echo done $@
../lib/WHERE.SLEPc-complex:$(SRCDIR_SLEPC)/tag-install-complex
	test -e $(DIR_INSTALL_COMPLEX)/include/slepc.h
	-mkdir ../lib
	echo 'slepccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib"   -lslepc' > $@
	echo 'slepccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"'  >> $@
	echo done $@

##  need include include Make-petsc-download.mk
WHERE-all:Makefile WHERE WHERE-complex
	echo scalapack LD $(SCALAPACK_LIB) >$@
	echo scalapack INCLUDE -I$(CALAPACK_INCLUDE) >>$@
	echo metis LD $(METIS_LIB) >>$@
	echo metis INCLUDE -I$(METIS_INCLUDE) >>$@
	echo mumps LD $(MUMPS_LIB) >>$@
	echo mumps INCLUDE -I$(MUMPS_INCLUDE) >>$@
	echo superlu LD $(SUPERLU_LIB) >>$@
	echo superlu INCLUDE -I$(SUPERLU_INCLUDE) >>$@
	echo fftw LD $(FFTW_LIB) >>$@
	echo fftw INCLUDE -I$(FFTW_INCLUDE)  >>$@
	echo ptscotch LD $(PTSCOTCH_LIB) >>$@
	echo ptscotch INCLUDE -I$(PTSCOTCH_INCLUDE) >>$@
	echo parmetis LD $(PARMETIS_LIB)>>$@
	echo parmetis INCLUDE -I$(PARMETIS_INCLUDE)>>$@
	echo 'petsc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib"   -lpetsc' >> $@
	echo 'petsc INCLUDE -I"$(DIR_INSTALL_REAL)/include"'  >> $@
	echo 'petsccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib"   -lpetsc' >> $@
	echo 'petsccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"'  >> $@
	echo 'slepc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib"   -lslepc' >> $@
	echo 'slepc INCLUDE -I"$(DIR_INSTALL_REAL)/include"'  >> $@
	echo 'slepccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib"   -lslepc' >> $@
	echo 'slepccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"'  >> $@


$(SRCDIR_SLEPC)/tag-conf-complex:$(SRCDIR_SLEPC)/tag-tar
	export PETSC_DIR=$(PETSC_DIR) ; \
	export SLEPC_DIR=$(SLEPC_DIR) ; \
	export PETSC_ARCH=ff-complex ; \
	cd $(SRCDIR_SLEPC) &&   $(PYTHON2) ./configure  --prefix=$(DIR_INSTALL_COMPLEX) $(SPELC_ADDCONF)
	touch $@

$(SRCDIR_SLEPC)/tag-make-complex:$(SRCDIR_SLEPC)/tag-conf-complex
	export PETSC_DIR=$(PETSC_DIR) ; \
	export SLEPC_DIR=$(SLEPC_DIR) ; \
	export PETSC_ARCH=ff-complex ; \
	cd $(SRCDIR_SLEPC) &&   make
	touch $@


$(SRCDIR_SLEPC)/tag-install-complex :$(SRCDIR_SLEPC)/tag-make-complex
	export PETSC_DIR=$(PETSC_DIR) ; \
	export SLEPC_DIR=$(SLEPC_DIR) ; \
	export PETSC_ARCH=ff-complex ; \
	cd $(SRCDIR_SLEPC) &&   $(SUDO_WITH_ENV) make install
	test -e $(DIR_INSTALL_COMPLEX)/include/slepc.h
	cd $(SRCDIR_SLEPC) &&   make 'PETSC_DIR=$(PETSC_DIR)' 'SLEPC_DIR=$(SLEPC_DIR)' PETSC_ARCH=ff-complex check
	touch $@
#
else
#
WHERE-complex:
	@if [  -n "$(MPI_INCLUDE)" ]; then $(MAKE) -f Makefile.complex WHERE-complex ;\
	else echo " -- No  Petc , no MPI"; fi

WHERE-all:Makefile  Make-petsc-download.mk
	 $(MAKE) -f Makefile.complex $@
endif
$(SRCDIR)/tag-tar:$(PACKAGE)
	tar xvzf $(PACKAGE)
	touch $@
$(SRCDIR_SLEPC)/tag-tar:$(PACKAGE_SLEPC)
	rm -rf $(SRCDIR_SLEPC)
	tar xvzf $(PACKAGE_SLEPC)
#	mv `tar tfz $(PACKAGE_SLEPC) | head -1` $(SRCDIR_SLEPC)
#	cd $(SRCDIR_SLEPC);patch -p1 <../patch-slepc-3.8.3
	touch $@
$(PACKAGE):
	../getall -o PETSc -a
$(PACKAGE_SLEPC):
	../getall -o SLEPc -a

$(SRCDIR_SLEPC)/tag-conf-real:$(SRCDIR_SLEPC)/tag-tar
	export PETSC_DIR=$(PETSC_DIR) ; \
	export SLEPC_DIR=$(SLEPC_DIR) ; \
	export PETSC_ARCH=ff-real; \
	cd $(SRCDIR_SLEPC) &&   $(PYTHON2) ./configure --prefix=$(DIR_INSTALL_REAL) $(SPELC_ADDCONF)
	touch $@

$(SRCDIR_SLEPC)/tag-make-real:$(SRCDIR_SLEPC)/tag-conf-real
	export PETSC_DIR=$(PETSC_DIR) ; \
	export SLEPC_DIR=$(SLEPC_DIR) ; \
	export PETSC_ARCH=ff-real; \
	cd $(SRCDIR_SLEPC) &&   make
	touch $@


$(SRCDIR_SLEPC)/tag-install-real :$(SRCDIR_SLEPC)/tag-make-real
	export PETSC_DIR=$(DIR_INSTALL_REAL) ; \
	export SLEPC_DIR=$(SLEPC_DIR) ; \
	export PETSC_ARCH=ff-real; \
	cd $(SRCDIR_SLEPC) &&   $(SUDO_WITH_ENV) make install
	test -e $(DIR_INSTALL_REAL)/include/slepc.h
	cd $(SRCDIR_SLEPC) &&   make 'PETSC_DIR=$(PETSC_DIR)' 'SLEPC_DIR=$(SLEPC_DIR)' PETSC_ARCH=ff-real check

	touch $@

../lib/WHERE.SLEPc:$(SRCDIR_SLEPC)/tag-install-real
	-mkdir ../lib
	test -e $(DIR_INSTALL_REAL)/include/slepc.h
	echo 'slepc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib"   -lslepc' > $@
	echo 'slepc INCLUDE -I"$(DIR_INSTALL_REAL)/include"'  >> $@


clean:clean-local
clean-local:
	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR)
	-rm Makefile.inc FAIRE* ../lib/WHERE.PETSc* ../lib/WHERE.SLEPc*
	-rm -rf ../include/*PETSc*
	-rm -rf ../lib/lib*PETSc*
	-rm -rf $(SRCDIR)
	-rm -rf $(SRCDIR_SLEPC)
	-rm -rf $(prefix)/ff-petsc
	-rm WHERE-all config.log *.done
	-if test -d $(prefix)/ff-petsc ; then echo " try of  remove of $(prefix)/ff-petsc under sudo .."; sudo rm -rf $(prefix)/ff-petsc ; fi
	-rm do-sudo
echo:
	@echo " with_dir_mpi #$(W_MPICC)# #$(MPICC)# ::: #$(WITH_MPI_DIR)#  "
	@echo SUDO: $(SUDO)
	@echo "sudo using user env: $(SUDO_WITH_ENV)"
	@echo MPI_DIR: $(MPI_DIR)
	@echo " dir install real :" $(DIR_INSTALL_REAL)
	@echo " dir install complex :" $(DIR_INSTALL_COMPLEX)
	@echo " do-sudo auto ???:" $(SHELL  cat do-sudo)
	@echo " BLAS_LAPACK_LIBS: $(BLAS_LAPACK_LIBS)"
	@echo " # need sudo if '$(W_SUDO)' == 'sudo'  .." $(W_SUDO1)

install-destdir:
	cd $(SRCDIR) && $(SUDO) $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-real install  DESTDIR=$(DESTDIR)
	cd $(SRCDIR) && $(SUDO) $(MAKE)  PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-complex install DESTDIR=$(DESTDIR)
	cd $(SRCDIR_SLEPC) &&   $(SUDO) make install SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(DIR_INSTALL_COMPLEX) DESTDIR=$(DESTDIR)
	cd $(SRCDIR_SLEPC) &&   $(SUDO) make install SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(DIR_INSTALL_REAL) DESTDIR=$(DESTDIR)

.NOTPARALLEL:
# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
