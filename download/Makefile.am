# Downloading and compiling extra libraries
# -----------------------------------------

# $Id: Makefile.am,v 1.16 2010/05/06 21:20:38 hecht Exp $

SUBDIRS=blas arpack umfpack  
EXTRA_DIST=f2c/Makefile f2c/f2c.h-int   f2c/fort77.sed	f2c/Makefile-MacOs tetgen/Makefile \
tetgen/tetgen1.4.2.patch superlu/Makefile \
fftw/Makefile metis/Makefile metis/Makefile-metis.in metis/patch-metis \
mmg3d/Makefile			mshmet/mshmetlib.h \
mmg3d/makefile-mmg3d.inc	yams/Makefile \
yams/freeyams.2010.02.22.patch \
mshmet/Makefile			yams/makefile-yams.inc \
mshmet/Makefile-mshmet.inc	yams/yamslib.c \
mshmet/mshmet.2009.09.16.patch	yams/yamslib.h \
mshmet/mshmetlib-internal.h	yams/yamslib_internal.h \
mshmet/mshmetlib.c \
mumps/MUMPS_4.9.2.patch			parmetis/makefile \
mumps/Makefile				parmetis/parmetis-3.1.1.patch \
mumps/Makefile-mumps-4.9.2.inc		pastix/Makefile \
blacs/BLACS.patch			pastix/config-pastix-complex.in \
blacs/Bmake-blacs.inc			pastix/config-pastix-real.in \
blacs/Makefile				pastix/pastix_release_2200-blend.patch \
pastix/pastix_long_complex.h		pastix/patch-pastix_long_complex.h \
hips/Makefile				pastix/pastix_release_2200-gen.patch \
hips/SRC_SPKIT_makefile			scalapack/Makefile \
hips/hips-1.2b-rc4.patch		scalapack/SLmake-scalapack.inc \
scotch/Makefile \
hips/makefile-hips.inc			scotch/Makefile-scotch.inc \
hypre/Makefile				scotch/scotch_5.1_esmumps.patch \
superludist/Makefile \
parmetis/Makefile-parmetis.in		superludist/make-superlu.inc \
superludist/superludist_2.3.patch  \
headers-sparsesolver.inc \
gmm/Makefile	gmm/cxxflags


MPI_SOFT=blacs parmetis  scalapack scotch superludist MUMPS pastix hypre hips	

all-recursive: bin lib include pkg

lib: 
	mkdir lib
bin: 
	mkdir bin
include:
	mkdir  include
pkg:
	mkdir pkg

clean-local:
	-rm -rf  lib include 
	-rm */FAIT 
	for d in tetgen superlu fftw; do \
	  ( cd "$$d";$(MAKE) clean ) ; \
	done 
all-local: $(DOWNLOADCOMPILE) 

compile-pkg:
	-if [ -n "$(WGET)" ] ; then \
	for d in tetgen superlu fftw metis yams mshmet mmg3d gmm ; do \
	if [ -d "$$d" ] ; then  echo "(cd $$d;$(MAKE) ) " ;   ( cd "$$d";$(MAKE) ) ; else echo ** missing "$$d" directory ; fi; \
	done ;\
	if [ -n "$(MPICC)" ] ; then \
	for d in $(MPI_SOFT) ; do \
	if [ -d "$$d" ] ; then  echo "(cd $$d;$(MAKE) ) " ;   ( cd "$$d";$(MAKE) ) ; else echo ** missing "$$d" directory ;  fi; \
	done ;\
	fi; \
	fi 