# Downloading and compiling extra libraries
# -----------------------------------------

# $Id: Makefile.am,v 1.16 2010/05/06 21:20:38 hecht Exp $

SUBDIRS=blas arpack umfpack  
EXTRA_DIST=f2c/Makefile f2c/f2c.h-int   f2c/fort77.sed	f2c/Makefile-MacOs tetgen/Makefile \
tetgen/tetgen1.4.2.patch superlu/Makefile \
fftw/Makefile metis/Makefile metis/Makefile-metis.in metis/patch-metis
MPI_SOFT=blacs parmetis  scalapack scotch superludist mumps pastix hypre hips	

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
	  ( cd "$$d";$(MAKE) clean ) ; \
	done 
all-local: $(DOWNLOADCOMPILE) 

compile-pkg:
	-if [ -n "$(WGET)" ] ; then \
	for d in tetgen superlu fftw metis yams mshmet mmg3d ; do \
	 echo "(cd $$d;$(MAKE) ) " ;\
	  ( cd "$$d";$(MAKE) ) ; \
	done ;\
	if [ -n "$(MPICC)" ] ; then \
	for d in $(MPI_SOFT) ; do \
	 echo "(cd $$d;$(MAKE) ) " ;\
	  ( cd "$$d";$(MAKE) ) ; \
	done ;\
	fi; \
	fi 