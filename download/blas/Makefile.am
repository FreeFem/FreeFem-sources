# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

ATLASSRCLIB=ATLAS/lib/ff++/libatlas.a ATLAS/lib/ff++/libcblas.a	\
	ATLAS/lib/ff++/libf77blas.a
ATLASBINLIB=WinNT_PPRO256/libatlas.a

all-local:: $(DOWNLOAD_BLAS)
generic: 
atlas-source: $(ATLASSRCLIB)
atlas-binary: $(ATLASBINLIB)

# Downloading and compiling the Blas
# ----------------------------------

# ATLAS source
$(ATLASSRCLIB): ATLAS/Make.ff++
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
$(ATLASBINLIB): atlas321_WinNT_PPRO.zip
	unzip atlas321_WinNT_PPRO.zip
	touch $(ATLASBINLIB)

atlas321_WinNT_PPRO.zip:
	wget -N http://www.netlib.org/atlas/archives/windows/atlas321_WinNT_PPRO.zip

clean-local::
	-rm atlas321_WinNT_PPRO.zip
	-rm -r WinNT_PPRO256

# Generic BLAS
noinst_LIBRARIES=@GENERIC_BLAS@
EXTRA_LIBRARIES=libf77blas.a libcblas.a
libf77blas_a_SOURCES=
libcblas_a_SOURCES=
BUILT_SOURCES=@GENERIC_BLAS_BUILT_SOURCES@

libcblas_a_CPPFLAGS=-DADD_

# "xerbla" exists in both BLAS and CBLAS. So we need to rename it to
# obtain two different object files.

BLAS: blas.tgz
	mkdir -p BLAS
	cd BLAS && tar xvzf ../blas.tgz
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
