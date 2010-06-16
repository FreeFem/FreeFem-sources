# Downloading and compiling extra libraries
# -----------------------------------------

# $Id: Makefile.am,v 1.12 2009/10/07 02:54:53 hecht Exp $

SUBDIRS=blas arpack umfpack
EXTRA_DIST=f2c/Makefile f2c/f2c.h-int   f2c/fort77.sed	f2c/Makefile-MacOs tetgen/Makefile tetgen/tetgen1.4.2.patch superlu/Makefile \ 
./blacs/BLACS.patch \
./blacs/Bmake-blacs.inc \
./blacs/Makefile \
./f2c/f2c.h-int \
./f2c/fort77.sed \
./f2c/Makefile \
./f2c/Makefile-MacOs \
./f2c/tt.f \
./fftw/Makefile \
./fftw/Makefile.am \
./headers-sparsesolver.inc \
./hips/hips-1.2b-rc4-1.patch \
./hips/hips-1.2b-rc4-2.patch \
./hips/Makefile \
./hips/makefile-hips.inc \
./hips/SRC_SPKIT_makefile \
./hypre/ff-flags \
./hypre/Makefile \
./metis/Makefile \
./metis/Makefile-metis.in \
./metis/patch-metis \
./mmg3d/Makefile \
./mmg3d/makefile-mmg3d.inc \
./mshmet/Makefile \
./mshmet/Makefile-mshmet.inc \
./mshmet/mshmet.2009.09.16.patch \
./mshmet/mshmetlib-internal.h \
./mshmet/mshmetlib.c \
./mshmet/mshmetlib.h \
./mumps/Makefile \
./mumps/Makefile-mumps-4.9.2.inc \
./mumps/MUMPS_4.9.2.patch \
./parmetis/makefile \
./parmetis/Makefile-parmetis.in \
./parmetis/ParMetis-3.1.1-ok/stdheaders.h \
./parmetis/parmetis-3.1.1.patch \
./parms/Makefile \
./parms/makefile-parms.in \
./pastix/config-pastix-complex.in \
./pastix/config-pastix-real.in \
./pastix/Makefile \
./pastix/pastix_release_2200-blend.patch \
./pastix/pastix_release_2200-gen.patch \
./scalapack/Makefile \
./scalapack/SLmake-scalapack.inc \
./scotch/Makefile \
./scotch/Makefile-scotch.inc \
./scotch/scotch_5.1_esmumps.patch \
./superlu/make.inc \
./superlu/Makefile \
./superludist/make-superlu.inc \
./superludist/Makefile \
./superludist/superludist_2.3.patch \
./tetgen/Makefile \
./tetgen/tetgen1.4.2.patch \
./yams/Makefile \
./yams/makefile-yams.inc \
./yams/freeyams.2010.02.22.patch

all-recursive: lib include pkg

lib: 
	mkdir lib
include:
	mkdir  include
pkg:
	mkdir pkg

clean-local:
	-rm -rf  lib include 
	-rm */FAIT 
	for d in tetgen superlu fftw; do \
	if [ -d $$d ] ; then   ( cd "$$d";$(MAKE) clean ) ; fi; \
	done 
all-local: $(DOWNLOADCOMPILE) 

compile-pkg:
	-if [ -n "$(WGET)" ] ; then \
	for d in tetgen superlu fftw metis; do \
	 echo "(cd $$d;$(MAKE) ) " ;\
	 if [ -d $d ] ; then   ( cd "$$d";$(MAKE) clean ) ; fi; \	
	done ;\
	fi 