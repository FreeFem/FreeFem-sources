# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

ATLASSRCLIB=ATLAS/lib/ff++/libatlas.a ATLAS/lib/ff++/libcblas.a	\
	ATLAS/lib/ff++/libf77blas.a
ATLASBINLIB=WinNT_PPRO256/libatlas.a

all-local:: $(DOWNLOAD_BLAS)
generic: 
atlas-source: $(SRCLIB)
atlas-binary: $(BINLIB)

# Downloading and compiling the Blas
# ----------------------------------

# ATLAS source
$(SRCLIB): ATLAS/Make.ff++
	cd ATLAS && make install arch=ff++

ATLAS/Make.ff++: atlas3.6.0.tar.bz2
	tar -x --bzip -f atlas3.6.0.tar.bz2
	@echo In the following questions, \"Architecture name\" must be ff++
	@echo "(Press ENTER)"
	@read
	cd ATLAS && make

atlas3.6.0.tar.bz2:
	wget -N http://heanet.dl.sourceforge.net/sourceforge/math-atlas/atlas3.6.0.tar.bz2

clean-local::
	-rm atlas*.tar.bz2
	-rm -r ATLAS

# ATLAS binary
$(BINLIB): atlas321_WinNT_PPRO.zip
	unzip atlas321_WinNT_PPRO.zip
	touch $(BINLIB)

atlas321_WinNT_PPRO.zip:
	wget -N http://www.netlib.org/atlas/archives/windows/atlas321_WinNT_PPRO.zip

clean-local::
	-rm atlas321_WinNT_PPRO.zip
	-rm -r WinNT_PPRO256

# Generic BLAS
noinst_LIBRARIES=@GENERIC_BLAS@
EXTRA_LIBRARIES=libf77blas.a libcblas.a

libf77blas_a_SOURCES=blas/caxpy.f blas/crotg.f blas/dcopy.f		\
blas/dsymv.f blas/lsame.f blas/sspmv.f blas/zaxpy.f blas/zhpr2.f	\
blas/ccopy.f blas/cscal.f blas/ddot.f blas/dsyr2.f blas/sasum.f		\
blas/sspr2.f blas/zcopy.f blas/zhpr.f blas/cdotc.f blas/csrot.f		\
blas/dgbmv.f blas/dsyr2k.f blas/saxpy.f blas/sspr.f blas/zdotc.f	\
blas/zrotg.f blas/cdotu.f blas/csscal.f blas/dgemm.f blas/dsyr.f	\
blas/scasum.f blas/sswap.f blas/zdotu.f blas/zscal.f blas/cgbmv.f	\
blas/cswap.f blas/dgemv.f blas/dsyrk.f blas/scnrm2.f blas/ssymm.f	\
blas/zdrot.f blas/zswap.f blas/cgemm.f blas/csymm.f blas/dger.f		\
blas/dtbmv.f blas/scopy.f blas/ssymv.f blas/zdscal.f blas/zsymm.f	\
blas/cgemv.f blas/csyr2k.f blas/dnrm2.f blas/dtbsv.f blas/sdot.f	\
blas/ssyr2.f blas/zgbmv.f blas/zsyr2k.f blas/cgerc.f blas/csyrk.f	\
blas/drot.f blas/dtpmv.f blas/sdsdot.f blas/ssyr2k.f blas/zgemm.f	\
blas/zsyrk.f blas/cgeru.f blas/ctbmv.f blas/drotg.f blas/dtpsv.f	\
blas/sgbmv.f blas/ssyr.f blas/zgemv.f blas/ztbmv.f blas/chbmv.f		\
blas/ctbsv.f blas/drotm.f blas/dtrmm.f blas/sgemm.f blas/ssyrk.f	\
blas/zgerc.f blas/ztbsv.f blas/chemm.f blas/ctpmv.f blas/drotmg.f	\
blas/dtrmv.f blas/sgemv.f blas/stbmv.f blas/zgeru.f blas/ztpmv.f	\
blas/chemv.f blas/ctpsv.f blas/dsbmv.f blas/dtrsm.f blas/sger.f		\
blas/stbsv.f blas/zhbmv.f blas/ztpsv.f blas/cher2.f blas/ctrmm.f	\
blas/dscal.f blas/dtrsv.f blas/snrm2.f blas/stpmv.f blas/zhemm.f	\
blas/ztrmm.f blas/cher2k.f blas/ctrmv.f blas/dsdot.f blas/dzasum.f	\
blas/srot.f blas/stpsv.f blas/zhemv.f blas/ztrmv.f blas/cher.f		\
blas/ctrsm.f blas/dspmv.f blas/dznrm2.f blas/srotg.f blas/strmm.f	\
blas/zher2.f blas/ztrsm.f blas/cherk.f blas/ctrsv.f blas/dspr2.f	\
blas/icamax.f blas/srotm.f blas/strmv.f blas/zher2k.f blas/ztrsv.f	\
blas/chpmv.f blas/dasum.f blas/dspr.f blas/idamax.f blas/srotmg.f	\
blas/strsm.f blas/zher.f blas/chpr2.f blas/daxpy.f blas/dswap.f		\
blas/isamax.f blas/ssbmv.f blas/strsv.f blas/zherk.f blas/chpr.f	\
blas/dcabs1.f blas/dsymm.f blas/izamax.f blas/sscal.f blas_xerbla.f	\
blas/zhpmv.f

libcblas_a_SOURCES=CBLAS/src/cblas_caxpy.c CBLAS/src/cblas_drot.c	\
CBLAS/src/cblas_sgemm.c CBLAS/src/cblas_zher2.c				\
CBLAS/src/cblas_ccopy.c CBLAS/src/cblas_drotg.c				\
CBLAS/src/cblas_sgemv.c CBLAS/src/cblas_zher2k.c			\
CBLAS/src/cblas_cdotc_sub.c CBLAS/src/cblas_drotm.c			\
CBLAS/src/cblas_sger.c CBLAS/src/cblas_zher.c				\
CBLAS/src/cblas_cdotu_sub.c CBLAS/src/cblas_drotmg.c			\
CBLAS/src/cblas_snrm2.c CBLAS/src/cblas_zherk.c				\
CBLAS/src/cblas_cgbmv.c CBLAS/src/cblas_dsbmv.c CBLAS/src/cblas_srot.c	\
CBLAS/src/cblas_zhpmv.c CBLAS/src/cblas_cgemm.c				\
CBLAS/src/cblas_dscal.c CBLAS/src/cblas_srotg.c				\
CBLAS/src/cblas_zhpr2.c CBLAS/src/cblas_cgemv.c				\
CBLAS/src/cblas_dsdot.c CBLAS/src/cblas_srotm.c CBLAS/src/cblas_zhpr.c	\
CBLAS/src/cblas_cgerc.c CBLAS/src/cblas_dspmv.c				\
CBLAS/src/cblas_srotmg.c CBLAS/src/cblas_zscal.c			\
CBLAS/src/cblas_cgeru.c CBLAS/src/cblas_dspr2.c				\
CBLAS/src/cblas_ssbmv.c CBLAS/src/cblas_zswap.c				\
CBLAS/src/cblas_chbmv.c CBLAS/src/cblas_dspr.c CBLAS/src/cblas_sscal.c	\
CBLAS/src/cblas_zsymm.c CBLAS/src/cblas_chemm.c				\
CBLAS/src/cblas_dswap.c CBLAS/src/cblas_sspmv.c				\
CBLAS/src/cblas_zsyr2k.c CBLAS/src/cblas_chemv.c			\
CBLAS/src/cblas_dsymm.c CBLAS/src/cblas_sspr2.c				\
CBLAS/src/cblas_zsyrk.c CBLAS/src/cblas_cher2.c				\
CBLAS/src/cblas_dsymv.c CBLAS/src/cblas_sspr.c CBLAS/src/cblas_ztbmv.c	\
CBLAS/src/cblas_cher2k.c CBLAS/src/cblas_dsyr2.c			\
CBLAS/src/cblas_sswap.c CBLAS/src/cblas_ztbsv.c CBLAS/src/cblas_cher.c	\
CBLAS/src/cblas_dsyr2k.c CBLAS/src/cblas_ssymm.c			\
CBLAS/src/cblas_ztpmv.c CBLAS/src/cblas_cherk.c CBLAS/src/cblas_dsyr.c	\
CBLAS/src/cblas_ssymv.c CBLAS/src/cblas_ztpsv.c				\
CBLAS/src/cblas_chpmv.c CBLAS/src/cblas_dsyrk.c				\
CBLAS/src/cblas_ssyr2.c CBLAS/src/cblas_ztrmm.c				\
CBLAS/src/cblas_chpr2.c CBLAS/src/cblas_dtbmv.c				\
CBLAS/src/cblas_ssyr2k.c CBLAS/src/cblas_ztrmv.c			\
CBLAS/src/cblas_chpr.c CBLAS/src/cblas_dtbsv.c CBLAS/src/cblas_ssyr.c	\
CBLAS/src/cblas_ztrsm.c CBLAS/src/cblas_cscal.c				\
CBLAS/src/cblas_dtpmv.c CBLAS/src/cblas_ssyrk.c				\
CBLAS/src/cblas_ztrsv.c CBLAS/src/cblas_csscal.c			\
CBLAS/src/cblas_dtpsv.c CBLAS/src/cblas_stbmv.c CBLAS/src/cdotcsub.f	\
CBLAS/src/cblas_cswap.c CBLAS/src/cblas_dtrmm.c				\
CBLAS/src/cblas_stbsv.c CBLAS/src/cdotusub.f CBLAS/src/cblas_csymm.c	\
CBLAS/src/cblas_dtrmv.c CBLAS/src/cblas_stpmv.c CBLAS/src/dasumsub.f	\
CBLAS/src/cblas_csyr2k.c CBLAS/src/cblas_dtrsm.c			\
CBLAS/src/cblas_stpsv.c CBLAS/src/ddotsub.f CBLAS/src/cblas_csyrk.c	\
CBLAS/src/cblas_dtrsv.c CBLAS/src/cblas_strmm.c CBLAS/src/dnrm2sub.f	\
CBLAS/src/cblas_ctbmv.c CBLAS/src/cblas_dzasum.c			\
CBLAS/src/cblas_strmv.c CBLAS/src/dsdotsub.f CBLAS/src/cblas_ctbsv.c	\
CBLAS/src/cblas_dznrm2.c CBLAS/src/cblas_strsm.c CBLAS/src/dzasumsub.f	\
CBLAS/src/cblas_ctpmv.c CBLAS/src/cblas_f77.h CBLAS/src/cblas_strsv.c	\
CBLAS/src/dznrm2sub.f CBLAS/src/cblas_ctpsv.c				\
CBLAS/src/cblas_globals.c CBLAS/src/cblas_xerbla.c			\
CBLAS/src/icamaxsub.f CBLAS/src/cblas_ctrmm.c CBLAS/src/cblas.h		\
CBLAS/src/cblas_zaxpy.c CBLAS/src/idamaxsub.f CBLAS/src/cblas_ctrmv.c	\
CBLAS/src/cblas_icamax.c CBLAS/src/cblas_zcopy.c CBLAS/src/isamaxsub.f	\
CBLAS/src/cblas_ctrsm.c CBLAS/src/cblas_idamax.c			\
CBLAS/src/cblas_zdotc_sub.c CBLAS/src/izamaxsub.f			\
CBLAS/src/cblas_ctrsv.c CBLAS/src/cblas_isamax.c			\
CBLAS/src/cblas_zdotu_sub.c CBLAS/src/Makefile CBLAS/src/cblas_dasum.c	\
CBLAS/src/cblas_izamax.c CBLAS/src/cblas_zdscal.c CBLAS/src/sasumsub.f	\
CBLAS/src/cblas_daxpy.c CBLAS/src/cblas_sasum.c				\
CBLAS/src/cblas_zgbmv.c CBLAS/src/scasumsub.f CBLAS/src/cblas_dcopy.c	\
CBLAS/src/cblas_saxpy.c CBLAS/src/cblas_zgemm.c CBLAS/src/scnrm2sub.f	\
CBLAS/src/cblas_ddot.c CBLAS/src/cblas_scasum.c				\
CBLAS/src/cblas_zgemv.c CBLAS/src/sdotsub.f CBLAS/src/cblas_dgbmv.c	\
CBLAS/src/cblas_scnrm2.c CBLAS/src/cblas_zgerc.c CBLAS/src/sdsdotsub.f	\
CBLAS/src/cblas_dgemm.c CBLAS/src/cblas_scopy.c				\
CBLAS/src/cblas_zgeru.c CBLAS/src/snrm2sub.f CBLAS/src/cblas_dgemv.c	\
CBLAS/src/cblas_sdot.c CBLAS/src/cblas_zhbmv.c CBLAS/src/xerbla.c	\
CBLAS/src/cblas_dger.c CBLAS/src/cblas_sdsdot.c				\
CBLAS/src/cblas_zhemm.c CBLAS/src/zdotcsub.f CBLAS/src/cblas_dnrm2.c	\
CBLAS/src/cblas_sgbmv.c CBLAS/src/cblas_zhemv.c CBLAS/src/zdotusub.f

BUILT_SOURCES=@GENERIC_BLAS_BUILT_SOURCES@

libcblas_a_CPPFLAGS=-DADD_

# "xerbla" exists in both BLAS and CBLAS. So we need to rename it to
# obtain two different object files.

BLAS: blas.tgz
	mkdir -p BLAS
	cd blas && tar xvzf ../blas.tgz
	cp BLAS/xerbla.f blas_xerbla.f

CBLAS: cblas.tgz
	tar xvzf cblas.tgz

blas.tgz:
	wget -N http://www.netlib.org/blas/blas.tgz

cblas.tgz:
	wget -N http://www.netlib.org/blas/blast-forum/cblas.tgz


clean-local::
	-rm -r BLAS CBLAS blas_xerbla.f
	-rm blas.tgz cblas.tgz
