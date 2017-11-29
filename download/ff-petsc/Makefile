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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with Freefem++; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
# ======================================================================
# headeralh default=0 freefem make multipleauthors start=19/03/10 upmc

include Makefile.inc

all-local:

# Downloading and compiling mumps
# -------------------------------

DIRPKG=../pkg
VERSION=3.8.0
VERSION_SLEPC=3.8.0
# http://slepc.upv.es/download/download.php?filename=slepc-3.7.4.tar.gz
PACKAGE=../pkg/petsc-lite-$(VERSION).tar.gz 
PACKAGE_SLEPC=../pkg/slepc-$(VERSION_SLEPC).tar.gz 
SRCDIR=petsc-$(VERSION)
SRCDIR_SLEPC=slepc-$(VERSION_SLEPC)
MPI_DIR:=$(shell dirname $(MPI_INC_DIR))
W_MPICC:=$(strip $(shell which mpicc))

ifeq ($(MPICC),$(W_MPICC))
  WITH_MPI_DIR=
else
  WITH_MPI_DIR="--with-mpi-dir='$(MPI_DIR)'"
endif	
DIR_INSTALL_REAL:=$(DESTDIR)$(prefix)/ff-petsc/real
DIR_INSTALL_COMPLEX:=$(DESTDIR)$(prefix)/ff-petsc/complex
PETSC_DIR:=$(PWD)/$(SRCDIR)
SLEPC_DIR:=$(PWD)/$(SRCDIR_SLEPC)
W_SUDO:=$(SHELL  mkdir -p "$(DESTDIR)$(prefix)/ff-petsc" && test -w "$(DESTDIR)$(prefix)/ff-petsc" || echo sudo)
petsc-slepc: WHERE-all

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
	cd $(SRCDIR) && $(W_SUDO) $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-real install  
	touch $@

Makefile.inc:../../config.status Makefile Makefile-PETSc.inc
	../../config.status --file="Makefile.inc:Makefile-PETSc.inc"

$(SRCDIR)/tag-conf-real:$(SRCDIR)/tag-tar
	cd $(SRCDIR) &&  python2 ./configure \
	--prefix=$(DIR_INSTALL_REAL) \
	--with-ssl=0 ---with-x=0 \
    '--with-scalar-type=real' \
    $(WITH_MPI_DIR) \
	--with-blas-lapack-lib='$(LAPACKLIBS) $(BLASLIBS)' \
	--download-scalapack --download-metis --download-ptscotch \
	--download-mumps --download-ml --download-hypre --download-parmetis \
	--download-superlu --download-pastix --download-suitesparse --download-fftw \
	PETSC_ARCH=ff-real
	touch $@
Make-petsc-download.mk:$(SRCDIR)/tag-install-real $(DIR_INSTALL_REAL)/lib/petsc/conf/petscvariables
	egrep 'SCALAPACK_|METIS_|MUMPS_|SUPERLU_|FFTW_|PTSCOTCH_|SUITESPARSE_'   $(DIR_INSTALL_REAL)/lib/petsc/conf/petscvariables | sed 's/-I/ /g'|sort >$@

ifdef 	MUMPS_LIB	
#  version COMPLEX  .....  
$(SRCDIR)/tag-conf-complex:$(SRCDIR)/tag-tar Make-petsc-download.mk
	cd $(SRCDIR) &&  python2 ./configure \
	--prefix=''$(DIR_INSTALL_COMPLEX)'' \
    '--with-scalar-type=complex' \
	--with-ssl=0 ---with-x=0 \
    $(WITH_MPI_DIR)  \
	--with-blas-lapack-lib='$(LAPACKLIBS) $(BLASLIBS)' \
	--with-scalapack-include='$(SCALAPACK_INCLUDE)'  --with-scalapack-lib='$(SCALAPACK_LIB)' \
	--with-metis-include='$(METIS_INCLUDE)'  --with-metis-lib='$(METIS_LIB)' \
	--with-ptscotch-include='$(PTSCOTCH_INCLUDE)'  --with-ptscotch-lib='$(PTSCOTCH_LIB)' \
	--with-suitesparse-include='$(SUITESPARSE_INCLUDE)'  --with-suitesparse-lib='$(SUITESPARSE_LIB)' \
	--with-mumps-include='$(MUMPS_INCLUDE)'  --with-mumps-lib='$(MUMPS_LIB)' \
	--with-parmetis-include='$(PARMETIS_INCLUDE)'  --with-parmetis-lib='$(PARMETIS_LIB)' \
	--with-fftw-include='$(FFTW_INCLUDE)'  --with-fftw-lib='$(FFTW_LIB)' \
	PETSC_ARCH=ff-complex

	touch $@

$(SRCDIR)/tag-make-complex:$(SRCDIR)/tag-conf-complex Make-petsc-download.mk 
	cd $(SRCDIR) && $(MAKE) PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-complex all
	touch $@

$(SRCDIR)/tag-install-complex :$(SRCDIR)/tag-make-complex 
	cd $(SRCDIR) && $(W_SUDO) $(MAKE)  PETSC_DIR=$(PETSC_DIR) PETSC_ARCH=ff-complex install
	touch $@
WHERE-complex:../lib/WHERE.PETSc-complex ../lib/WHERE.SLEPc-complex 

../lib/WHERE.PETSc-complex:$(SRCDIR)/tag-install-complex
	-mkdir ../lib
	echo 'petsccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib"   -lpetsc' > $@
	echo 'petsccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"'  >> $@
../lib/WHERE.SLEPc-complex:$(SRCDIR_SLEPC)/tag-install-complex
	-mkdir ../lib
	echo 'slepccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib"   -lslepc' > $@
	echo 'slepccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"'  >> $@


##  need include include Make-petsc-download.mk	
WHERE-all:Makefile WHERE-complex WHERE
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
	export PETSC_DIR=$(DIR_INSTALL_COMPLEX) ; \
	export SLEPC_DIR=$(SLEPC_DIR) ; \
	cd $(SRCDIR_SLEPC) &&   python2 ./configure  --prefix=$(DIR_INSTALL_COMPLEX) 
	touch $@

$(SRCDIR_SLEPC)/tag-make-complex:$(SRCDIR_SLEPC)/tag-conf-complex
	cd $(SRCDIR_SLEPC) &&   make SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(DIR_INSTALL_COMPLEX) 
	touch $@


$(SRCDIR_SLEPC)/tag-install-complex :$(SRCDIR_SLEPC)/tag-make-complex
	cd $(SRCDIR_SLEPC) &&   $(W_SUDO) make install SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(DIR_INSTALL_COMPLEX) 
	touch $@

../lib/WHERE.SLEPc-complex:$(SRCDIR_SLEPC)/tag-install-complex
	-mkdir ../lib
	echo 'slepccomplex LD -Wl,-rpath,"$(DIR_INSTALL_COMPLEX)/lib" -L"$(DIR_INSTALL_COMPLEX)/lib"   -lslepc' > $@
	echo 'slepccomplex INCLUDE -I"$(DIR_INSTALL_COMPLEX)/include"'  >> $@	
else	
WHERE-all:Makefile  Make-petsc-download.mk
	 $(MAKE) -f Makefile.complex $@ 
endif
$(SRCDIR)/tag-tar:$(PACKAGE)
	tar xvzf $(PACKAGE)
	touch $@
$(SRCDIR_SLEPC)/tag-tar:$(PACKAGE_SLEPC)
	tar xvzf $(PACKAGE_SLEPC)
#	patch -p1 <patch-slepc-3.7.4
	touch $@
$(PACKAGE):
	../getall -o PETSc -a
$(PACKAGE_SLEPC):
	../getall -o SLEPc -a

$(SRCDIR_SLEPC)/tag-conf-real:$(SRCDIR_SLEPC)/tag-tar
	export PETSC_DIR=$(DIR_INSTALL_REAL) ; \
	export SLEPC_DIR=$(SLEPC_DIR) ; \
	cd $(SRCDIR_SLEPC) &&   python2 ./configure  --prefix=$(DIR_INSTALL_REAL) 
	touch $@

$(SRCDIR_SLEPC)/tag-make-real:$(SRCDIR_SLEPC)/tag-conf-real
	cd $(SRCDIR_SLEPC) &&   make SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(DIR_INSTALL_REAL) 
	touch $@


$(SRCDIR_SLEPC)/tag-install-real :$(SRCDIR_SLEPC)/tag-make-real 
	cd $(SRCDIR_SLEPC) &&   $(W_SUDO) make install SLEPC_DIR=$(SLEPC_DIR) PETSC_DIR=$(DIR_INSTALL_REAL) 
	touch $@
../lib/WHERE.SLEPc:$(SRCDIR_SLEPC)/tag-install-real
	-mkdir ../lib
	echo 'slepc LD -Wl,-rpath,"$(DIR_INSTALL_REAL)/lib" -L"$(DIR_INSTALL_REAL)/lib"   -lslepc' > $@
	echo 'slepc INCLUDE -I"$(DIR_INSTALL_REAL)/include"'  >> $@


clean:clean-local
clean-local:
	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR) 
	-rm Makefile.inc FAIRE* ../lib/WHERE.PETSc* ../lib/WHERE.SLEPc* 
	-rm -rf ../include/*PETSc*
	-rm -rf ../lib/lib*PETSc* 
	-rm -rf $(SRCDIR)
	-rm -rf $(prefix)/ff-petsc
#	-rm $(PACKAGE1)
	-rm WHERE-all config.log *.done
	-if test -d $(prefix)/ff-petsc ; then echo " try of  remove of $(prefix)/ff-petsc under sudo .."; sudo -rf $(prefix)/ff-petsc ; fi

echo:
	echo " with_dir_mpi #$(W_MPICC)# #$(MPICC)# ::: #$(WITH_MPI_DIR)#  "
	echo W_SUDO: $(W_SUDO)
	echo MPI_DIR: $(MPI_DIR)
# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
