# ======================================================================
# Laboratoire Jacques-Louis Lions
# Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
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
# headeralh default=0 freefem make multipleauthors start=04/06/04 upmc brief="Makefile for downloaded ARPACK"

all-local:$(DOWNLOAD_ARPACK)
EXTRA_DIST=ARmake.m4 arpack-patch-lapack.tar.gz 

COMMON_PACKTITLE=ARPACK
include ../common.mak

# nothing specific to do for [[file:../common.mak::reinstall]]

reinstall::install

# Downloading and compiling ARPACK
# --------------------------------

# set in configure 
#ARPACKLIB=ARPACK/libarpack_ff++.a
DIRPKG=../pkg
ARPACK96_TAR_GZ=$(DIRPKG)/arpack96.tar.gz
PATCH_TAR_GZ=$(DIRPKG)/patch.tar.gz

COMMON_PACKAGES=$(ARPACK96_TAR_GZ) $(PATCH_TAR_GZ)

# FFCS: need to add $(RANLIB) under mingw64 to avoid "archive has no index" error

compile::$(ARPACKLIB)
$(ARPACKLIB): ARPACK/fait
	if [ -n '@FF_LAPACKdir@' ] ; then \
	$(F77) -c `echo $(FFLAGS)\ |sed -e s/-O.\*\ // ` ARPACK/LAPACK/dlamch.f -o ARPACK/LAPACK/dlamch.o; \
	fi; \
	cd ARPACK && make lib 
	if test -n '@FF_LAPACKdir@' ; then \
	$(AR) $(ARFLAGS) $(LAPACK_arpack_LIB)  ARPACK/SRC/*.o  ARPACK/UTIL/*.o ARPACK/LAPACK/*.o ;\
	$(RANLIB) $(LAPACK_arpack_LIB) ;\
	else \
	$(AR) $(ARFLAGS) $(ARPACKLIB)  ARPACK/SRC/*.o  ARPACK/UTIL/*.o ;\
	fi

ARPACK/fait: $(ARPACK96_TAR_GZ) $(PATCH_TAR_GZ) ARmake.m4 Makefile
	-rm -rf ARPACK
	gunzip -c $(ARPACK96_TAR_GZ) | tar xf -
	gunzip -c $(PATCH_TAR_GZ) | tar xf -
	gunzip -c arpack-patch-lapack.tar.gz | tar xf -
	for i in ARPACK/SRC/*.f ; do \
	    mv  $$i $$i.cpy; sed -e 's/, second/, secnd2/' -e 's/call *second/call secnd2/' <$$i.cpy >$$i;rm $$i.cpy; done 
	for i in  ARPACK/UTIL/second.f;  do  \
	   mv  $$i $$i.cpy; cat $$i.cpy| sed 's/ SECOND *(/ secnd2(/'|grep -v EXTERNAL  >$$i;rm $$i.cpy; done 
	m4  -DFF_BLASLIB="$(BLASLIB)" \
	    -DFF_ARPACKLIB="$(ARPACKLIB)" \
            -DFF_LAPACK_arpack_LIB="$(LAPACK_arpack_LIB)" \
	    -DFF_FC="@F77@" \
            -DFF_FFLAGS="@FFLAGS@" \
	    -DFF_LAPACKdir='@FF_LAPACKdir@' \
            -DFF_LDFLAGS="@LDFLAGS@" \
            -DFF_HOME=`pwd`/ARPACK \
	    -DFF_SECOND="@FF_SECOND@" \
	    -DFF_AR="@AR@" \
	    -DFF_ARFLAGS="@ARFLAGS@" \
	    -DFF_RANLIB="@RANLIB@" \
	    ARmake.m4 >ARPACK/ARmake.inc
	touch ARPACK/fait

clean-local::
	-rm -r ARPACK ../lib/libarpack.a

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
