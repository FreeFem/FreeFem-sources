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
# SUMMARY : Downloading and compiling extra BLAS libraries
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

DIRPKG=../pkg
BLAS_TGZ=$(DIRPKG)/blas-3.7.1.tgz
CBLAS_TGZ=$(DIRPKG)/cblas.tgz

PKGCOMMON_PACKTITLE=BLAS
PKGCOMMON_PACKAGES=$(BLAS_TGZ) $(CBLAS_TGZ) ../pkg/OpenBLAS.tar.gz
include ../common.mak

EXTRA_DIST= \
openblas.patches

# Downloading and compiling the Generic Blas
# ------------------------------------------

noinst_LIBRARIES=@DOWNLOADED_BLAS@
##compile::$(noinst_LIBRARIES)
EXTRA_LIBRARIES=libf77blas.a libcblas.a
pkgcompile::@DOWNLOADED_BLAS@

# List of files to compile (do not list them in *_SOURCES to prevent them from being included in distributions).

F77BLAS_SRC=BLAS/caxpy.f BLAS/crotg.f BLAS/dcopy.f BLAS/dsymv.f BLAS/lsame.f BLAS/sspmv.f BLAS/zaxpy.f BLAS/zhpr2.f	\
BLAS/ccopy.f BLAS/cscal.f BLAS/ddot.f BLAS/dsyr2.f BLAS/sasum.f BLAS/sspr2.f BLAS/zcopy.f BLAS/zhpr.f BLAS/cdotc.f	\
BLAS/csrot.f BLAS/dgbmv.f BLAS/dsyr2k.f BLAS/saxpy.f BLAS/sspr.f BLAS/zdotc.f BLAS/zrotg.f BLAS/cdotu.f BLAS/csscal.f	\
BLAS/dgemm.f BLAS/dsyr.f BLAS/scasum.f BLAS/sswap.f BLAS/zdotu.f BLAS/zscal.f BLAS/cgbmv.f BLAS/cswap.f BLAS/dgemv.f	\
BLAS/dsyrk.f BLAS/scnrm2.f BLAS/ssymm.f BLAS/zdrot.f BLAS/zswap.f BLAS/cgemm.f BLAS/csymm.f BLAS/dger.f BLAS/dtbmv.f	\
BLAS/scopy.f BLAS/ssymv.f BLAS/zdscal.f BLAS/zsymm.f BLAS/cgemv.f BLAS/csyr2k.f BLAS/dnrm2.f BLAS/dtbsv.f BLAS/sdot.f	\
BLAS/ssyr2.f BLAS/zgbmv.f BLAS/zsyr2k.f BLAS/cgerc.f BLAS/csyrk.f BLAS/drot.f BLAS/dtpmv.f BLAS/sdsdot.f BLAS/ssyr2k.f	\
BLAS/zgemm.f BLAS/zsyrk.f BLAS/cgeru.f BLAS/ctbmv.f BLAS/drotg.f BLAS/dtpsv.f BLAS/sgbmv.f BLAS/ssyr.f BLAS/zgemv.f	\
BLAS/ztbmv.f BLAS/chbmv.f BLAS/ctbsv.f BLAS/drotm.f BLAS/dtrmm.f BLAS/sgemm.f BLAS/ssyrk.f BLAS/zgerc.f BLAS/ztbsv.f	\
BLAS/chemm.f BLAS/ctpmv.f BLAS/drotmg.f BLAS/dtrmv.f BLAS/sgemv.f BLAS/stbmv.f BLAS/zgeru.f BLAS/ztpmv.f BLAS/chemv.f	\
BLAS/ctpsv.f BLAS/dsbmv.f BLAS/dtrsm.f BLAS/sger.f BLAS/stbsv.f BLAS/zhbmv.f BLAS/ztpsv.f BLAS/cher2.f BLAS/ctrmm.f	\
BLAS/dscal.f BLAS/dtrsv.f BLAS/snrm2.f BLAS/stpmv.f BLAS/zhemm.f BLAS/ztrmm.f BLAS/cher2k.f BLAS/ctrmv.f BLAS/dsdot.f	\
BLAS/dzasum.f BLAS/srot.f BLAS/stpsv.f BLAS/zhemv.f BLAS/ztrmv.f BLAS/cher.f BLAS/ctrsm.f BLAS/dspmv.f BLAS/dznrm2.f	\
BLAS/srotg.f BLAS/strmm.f BLAS/zher2.f BLAS/ztrsm.f BLAS/cherk.f BLAS/ctrsv.f BLAS/dspr2.f BLAS/icamax.f BLAS/srotm.f	\
BLAS/strmv.f BLAS/zher2k.f BLAS/ztrsv.f BLAS/chpmv.f BLAS/dasum.f BLAS/dspr.f BLAS/idamax.f BLAS/srotmg.f BLAS/strsm.f	\
BLAS/zher.f BLAS/chpr2.f BLAS/daxpy.f BLAS/dswap.f BLAS/isamax.f BLAS/ssbmv.f BLAS/strsv.f BLAS/zherk.f BLAS/chpr.f	\
BLAS/dcabs1.f BLAS/dsymm.f BLAS/izamax.f BLAS/sscal.f blas_xerbla.f BLAS/zhpmv.f

CBLAS_SRC=CBLAS/src/cblas_caxpy.c CBLAS/src/cblas_drot.c CBLAS/src/cblas_sgemm.c CBLAS/src/cblas_zher2.c		\
CBLAS/src/cblas_ccopy.c CBLAS/src/cblas_drotg.c CBLAS/src/cblas_sgemv.c CBLAS/src/cblas_zher2k.c			\
CBLAS/src/cblas_cdotc_sub.c CBLAS/src/cblas_drotm.c CBLAS/src/cblas_sger.c CBLAS/src/cblas_zher.c			\
CBLAS/src/cblas_cdotu_sub.c CBLAS/src/cblas_drotmg.c CBLAS/src/cblas_snrm2.c CBLAS/src/cblas_zherk.c			\
CBLAS/src/cblas_cgbmv.c CBLAS/src/cblas_dsbmv.c CBLAS/src/cblas_srot.c CBLAS/src/cblas_zhpmv.c CBLAS/src/cblas_cgemm.c	\
CBLAS/src/cblas_dscal.c CBLAS/src/cblas_srotg.c CBLAS/src/cblas_zhpr2.c CBLAS/src/cblas_cgemv.c CBLAS/src/cblas_dsdot.c	\
CBLAS/src/cblas_srotm.c CBLAS/src/cblas_zhpr.c CBLAS/src/cblas_cgerc.c CBLAS/src/cblas_dspmv.c CBLAS/src/cblas_srotmg.c	\
CBLAS/src/cblas_zscal.c CBLAS/src/cblas_cgeru.c CBLAS/src/cblas_dspr2.c CBLAS/src/cblas_ssbmv.c CBLAS/src/cblas_zswap.c	\
CBLAS/src/cblas_chbmv.c CBLAS/src/cblas_dspr.c CBLAS/src/cblas_sscal.c CBLAS/src/cblas_zsymm.c CBLAS/src/cblas_chemm.c	\
CBLAS/src/cblas_dswap.c CBLAS/src/cblas_sspmv.c CBLAS/src/cblas_zsyr2k.c CBLAS/src/cblas_chemv.c			\
CBLAS/src/cblas_dsymm.c CBLAS/src/cblas_sspr2.c CBLAS/src/cblas_zsyrk.c CBLAS/src/cblas_cher2.c CBLAS/src/cblas_dsymv.c	\
CBLAS/src/cblas_sspr.c CBLAS/src/cblas_ztbmv.c CBLAS/src/cblas_cher2k.c CBLAS/src/cblas_dsyr2.c CBLAS/src/cblas_sswap.c	\
CBLAS/src/cblas_ztbsv.c CBLAS/src/cblas_cher.c CBLAS/src/cblas_dsyr2k.c CBLAS/src/cblas_ssymm.c CBLAS/src/cblas_ztpmv.c	\
CBLAS/src/cblas_cherk.c CBLAS/src/cblas_dsyr.c CBLAS/src/cblas_ssymv.c CBLAS/src/cblas_ztpsv.c CBLAS/src/cblas_chpmv.c	\
CBLAS/src/cblas_dsyrk.c CBLAS/src/cblas_ssyr2.c CBLAS/src/cblas_ztrmm.c CBLAS/src/cblas_chpr2.c CBLAS/src/cblas_dtbmv.c	\
CBLAS/src/cblas_ssyr2k.c CBLAS/src/cblas_ztrmv.c CBLAS/src/cblas_chpr.c CBLAS/src/cblas_dtbsv.c CBLAS/src/cblas_ssyr.c	\
CBLAS/src/cblas_ztrsm.c CBLAS/src/cblas_cscal.c CBLAS/src/cblas_dtpmv.c CBLAS/src/cblas_ssyrk.c CBLAS/src/cblas_ztrsv.c	\
CBLAS/src/cblas_csscal.c CBLAS/src/cblas_dtpsv.c CBLAS/src/cblas_stbmv.c CBLAS/src/cdotcsub.f CBLAS/src/cblas_cswap.c	\
CBLAS/src/cblas_dtrmm.c CBLAS/src/cblas_stbsv.c CBLAS/src/cdotusub.f CBLAS/src/cblas_csymm.c CBLAS/src/cblas_dtrmv.c	\
CBLAS/src/cblas_stpmv.c CBLAS/src/dasumsub.f CBLAS/src/cblas_csyr2k.c CBLAS/src/cblas_dtrsm.c CBLAS/src/cblas_stpsv.c	\
CBLAS/src/ddotsub.f CBLAS/src/cblas_csyrk.c CBLAS/src/cblas_dtrsv.c CBLAS/src/cblas_strmm.c CBLAS/src/dnrm2sub.f	\
CBLAS/src/cblas_ctbmv.c CBLAS/src/cblas_dzasum.c CBLAS/src/cblas_strmv.c CBLAS/src/dsdotsub.f CBLAS/src/cblas_ctbsv.c	\
CBLAS/src/cblas_dznrm2.c CBLAS/src/cblas_strsm.c CBLAS/src/dzasumsub.f CBLAS/src/cblas_ctpmv.c CBLAS/src/cblas_f77.h	\
CBLAS/src/cblas_strsv.c CBLAS/src/dznrm2sub.f CBLAS/src/cblas_ctpsv.c CBLAS/src/cblas_globals.c				\
CBLAS/src/cblas_xerbla.c CBLAS/src/icamaxsub.f CBLAS/src/cblas_ctrmm.c CBLAS/src/cblas.h CBLAS/src/cblas_zaxpy.c	\
CBLAS/src/idamaxsub.f CBLAS/src/cblas_ctrmv.c CBLAS/src/cblas_icamax.c CBLAS/src/cblas_zcopy.c CBLAS/src/isamaxsub.f	\
CBLAS/src/cblas_ctrsm.c CBLAS/src/cblas_idamax.c CBLAS/src/cblas_zdotc_sub.c CBLAS/src/izamaxsub.f			\
CBLAS/src/cblas_ctrsv.c CBLAS/src/cblas_isamax.c CBLAS/src/cblas_zdotu_sub.c CBLAS/src/Makefile CBLAS/src/cblas_dasum.c	\
CBLAS/src/cblas_izamax.c CBLAS/src/cblas_zdscal.c CBLAS/src/sasumsub.f CBLAS/src/cblas_daxpy.c CBLAS/src/cblas_sasum.c	\
CBLAS/src/cblas_zgbmv.c CBLAS/src/scasumsub.f CBLAS/src/cblas_dcopy.c CBLAS/src/cblas_saxpy.c CBLAS/src/cblas_zgemm.c	\
CBLAS/src/scnrm2sub.f CBLAS/src/cblas_ddot.c CBLAS/src/cblas_scasum.c CBLAS/src/cblas_zgemv.c CBLAS/src/sdotsub.f	\
CBLAS/src/cblas_dgbmv.c CBLAS/src/cblas_scnrm2.c CBLAS/src/cblas_zgerc.c CBLAS/src/sdsdotsub.f CBLAS/src/cblas_dgemm.c	\
CBLAS/src/cblas_scopy.c CBLAS/src/cblas_zgeru.c CBLAS/src/snrm2sub.f CBLAS/src/cblas_dgemv.c CBLAS/src/cblas_sdot.c	\
CBLAS/src/cblas_zhbmv.c CBLAS/src/xerbla.c CBLAS/src/cblas_dger.c CBLAS/src/cblas_sdsdot.c CBLAS/src/cblas_zhemm.c	\
CBLAS/src/zdotcsub.f CBLAS/src/cblas_dnrm2.c CBLAS/src/cblas_sgbmv.c CBLAS/src/cblas_zhemv.c CBLAS/src/zdotusub.f

nodist_libf77blas_a_SOURCES= #$(F77BLAS_SRC)
nodist_libcblas_a_SOURCES= #$(CBLAS_SRC)
BUILT_SOURCES=@DOWNLOADED_BLAS_BUILT_SOURCES@

# -ICBLAS/include to find cblas.h
libcblas_a_CFLAGS=-DADD_ -ICBLAS/include

# "xerbla" exists in both BLAS and CBLAS. So we need to rename it to obtain two different object files.

BLAS:BLAS/fait

BLAS/fait:$(BLAS_TGZ)
	mkdir -p ../include ../lib
	pxerbla=`tar tf $(BLAS_TGZ) | grep xerbla.f`; \
	dirblas=`dirname $$pxerbla` ;\
	case $$dirblas in \
	BLAS) tar xvzf $(BLAS_TGZ) ;; \
	BLAS-*) tar xvzf $(BLAS_TGZ) ; mv $$dirblas BLAS ;; \
	.) mkdir -p BLAS; tar xvzf $(BLAS_TGZ) -c BLAS;; \
	esac
	cp BLAS/xerbla.f blas_xerbla.f
	touch BLAS/fait
$(F77BLAS_SRC): BLAS

CBLAS:CBLAS/fait

CBLAS/fait: $(CBLAS_TGZ)
	tar xvzf $(CBLAS_TGZ)
	cp CBLAS/include/*.h  CBLAS/src
	touch CBLAS/fait
$(CBLAS_SRC): CBLAS/fait

clean-local::
	-rm -r BLAS CBLAS blas_xerbla.f

# <<OpenBLAS>> ALH - 18/9/13 - Downloading and building the OpenBLAS
# ------------------------------------------------------------------

# to activate this, see [[file:../../configure.ac::OpenBLAS]]

all-local::@COMPILE_OPENBLAS@
pkgcompile::@COMPILE_OPENBLAS@
generic:

openblas:install.done ../lib/WHERE.blas

# The library may sometimes have a complex name with version number, just copy it into a standard location

../lib/WHERE.blas:install.done
	echo blas  LD -L@DIR@/lib  -lopenblas $(FLIBS)> $@
	echo blas INCLUDE -I@DIR@/include  >> $@
	echo blaslapack   LD -L@DIR@/lib  -lopenblas $(FLIBS) > $@
	echo blaslapack  INCLUDE -I@DIR@/include  >> $@


if FFCS_WINDOWS

# links2files is required for the MinGW compiler to understand where to find the library contents under Cygwin
links.done:openblas.done
	cd OpenBLAS && ../../../build/links2files
	touch $@
install.done:links.done
	cp OpenBLAS/libopenblas.a ../lib
	cd OpenBLAS && make BINARY=@SIZEOF_PTRINBIT@ CC=${CC} FC=${FC} NO_SHARED=1 DDYNAMIC_ARCH=1 PREFIX=../.. install
	touch $@
else
install.done:openblas.done
	cp OpenBLAS/libopenblas.a ../lib
	cd OpenBLAS && make BINARY=@SIZEOF_PTRINBIT@ CC=${CC} FC=${FC} DYNAMIC_ARCH=1 NO_SHARED=1 PREFIX=../.. install
	touch $@
endif

clean::
	-rm ../lib/libopenblas.a

openblas.done:openpatches.done
#
# DYNAMIC_ARCH=1 allows the Openblas to run fast on all the processor architectures that the FFCS users may have
#

	test  -d ../lib ||  mkdir  ../lib
	cd OpenBLAS && make BINARY=@SIZEOF_PTRINBIT@ "CC=$(CC)" "FC=$(FC)" DYNAMIC_ARCH=1 NO_SHARED=0  "FLIBS=$(FLIBS)" libs netlib
	touch $@

# The OpenBLAS directory is updated during the compilation, so the patching step should not depend on the directory date

openpatches.done:opendownload.done
	tar xvzf ../pkg/OpenBLAS.tar.gz
#
# ALH - 7/1/14 - The tar directory has changed names
#
	mv *OpenBLAS-* OpenBLAS
#	patch -u -p1 < openblas.patches
	touch $@

clean::
	-rm *.done
	-rm -r OpenBLAS

download::opendownload.done
opendownload.done:
	../getall -o OpenBLAS -a
	touch $@

# ALH - 6/11/13 - Since OpenBLAS is quite long to compile and it does not change very often, it's only cleaned as part of the specific
# target 'veryclean'.

veryclean::
	-rm -r *.done OpenBLAS xianyi-OpenBLAS-*

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
