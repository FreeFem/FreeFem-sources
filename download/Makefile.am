# Downloading and compiling extra libraries
# -----------------------------------------

# $Id: Makefile.am,v 1.16 2010/05/06 21:20:38 hecht Exp $

SUBDIRS=blas arpack umfpack  
EXTRA_DIST= \
./nlopt/Make.inc ./nlopt/Makefile \
./blacs/BLACS.patch \
./blacs/Bmake-blacs.inc \
./blacs/Makefile \
./f2c/Makefile \
./f2c/Makefile-MacOs \
./f2c/f2c.h-int \
./f2c/fort77.sed \
./f2c/tt.f \
./fftw/Makefile \
./fftw/Makefile.am \
./fltk/Makefile.am \
./fltk/Makefile.in \
./gmm/Makefile \
./gmm/cxxflags \
./headers-sparsesolver.inc \
./hips/Makefile \
./hips/SRC_SPKIT_makefile \
./hips/hips-1.2b-rc4.patch \
./hips/makefile-hips.inc \
./hypre/Makefile \
./hypre/ff-flags.inc \
./metis/Makefile \
./metis/Makefile-metis.in \
./metis/patch-metis \
./mmg3d/Makefile \
./mmg3d/patch-mmg3dv4.diff \
./mshmet/Makefile \
./mshmet/Makefile-mshmet.inc \
./mshmet/mshmet.2011.03.06.patch \
./mshmet/mshmetlib-internal.h \
./mshmet/mshmetlib.c \
./mshmet/mshmetlib.h \
./mumps/Makefile \
./mumps-seq/Makefile-mumps-4.10.0.inc \
./mumps-seq/Makefile \
./mumps/Makefile-mumps-4.10.0.inc \
./parmetis/Makefile-parmetis.in \
./parmetis/ParMetis-3.1.1-ok/stdheaders.h \
./parmetis/makefile \
./parmetis/parmetis-3.1.1.patch \
./parms/Makefile \
./parms/makefile-parms.in \
./pastix/Makefile \
./pastix/config-pastix-complex.in \
./pastix/config-pastix-real.in \
./pastix/pastix_release_2200-blend.patch \
./pastix/pastix_release_2200-gen.patch \
./pastix/patch-pastix_long_complex.h \
./scalapack/Makefile \
./scalapack/SLmake-scalapack.inc \
./scotch/Makefile \
./scotch/Makefile-scotch.inc \
./scotch/scotch_5.1_esmumps.patch \
./superlu/Makefile \
./superlu/make.inc \
./superludist/Makefile \
./superludist/make-superlu.inc \
./superludist/superludist_2.3.patch \
./tetgen/Makefile \
./tetgen/tetgen1.4.2.patch \
./yams/Makefile \
./yams/freeyams.2011.09.23.patch \
./yams/makefile-yams.inc \
./yams/yamslib.c \
./yams/yamslib.h \
./yams/yamslib_internal.h \
ipopt/Makefile	\
ipopt/Makefile.inc.in



MPI_SOFT=blacs parmetis  scalapack scotch superludist mumps pastix  hips parms 
LIST_SOFT=tetgen superlu fftw metis mshmet yams mmg3d gmm  nlopt mumps-seq ipopt 
all-recursive: bin lib include pkg

lib: 
	mkdir lib
bin: 
	mkdir bin
include:
	mkdir  include
pkg:
	mkdir pkg



all-local:bin lib include pkg  $(DOWNLOADCOMPILE) 

WHERE-LD: tag-compile-pkg
	touch ../examples++-load/WHERE_LIBRARY-config ../examples++-load/WHERE_LIBRARY		
	grep LD ../examples++-load/WHERE_LIBRARY ../examples++-load/WHERE_LIBRARY-config >WHERE-LD

compile-dir: 
	@echo "\n\n ****** $(COMPILEDIR) ****** \n\n";
	-@if [ 0 -eq `egrep ':$(COMPILEDIR)' WHERE-LD | wc -l` ] ;then \
	  cd $(COMPILEDIR) && $(MAKE) $(DIRTARGET) ; \
	else \
	 echo $(COMPILEDIR) is in WHERE- files ;\
	fi 
compile-pkg: tag-compile-pkg WHERE-LD
tag-compile-pkg: bin lib include pkg FORCE
	-@if [ -n "$(WGET)" ] ; then \
	for d in $(LIST_SOFT) ; do $(MAKE) compile-dir COMPILEDIR=$$d ;done ;\
	if [ -n "$(MPICC)" ] ; then \
	for d in $(MPI_SOFT) ; do 	$(MAKE) compile-dir COMPILEDIR=$$d ; done;\
	fi;fi
	touch tag-compile-pkg
FORCE:

re-install: 
	$(MAKE) compile-pkg DIRTARGET=install
WHERE:
	$(MAKE) compile-pkg DIRTARGET=WHERE

clean-local:
	-rm -rf	 tag-* include lib bin 
	mkdir   include lib bin 	
	-rm */FAIT */FAIRE 
	-if [ -n "$(WGET)" ] ; then \
	for d in  $(LIST_SOFT)  ; do \
	if [ -d "$$d" ] ; then  echo "$(MAKE) clean -C $$d  " ;   $(MAKE) clean -C $$d  ; else echo ** missing "$$d" directory ; fi; \
	done ;\
	if [ -n "$(MPICC)" ] ; then \
	for d in $(MPI_SOFT) ; do \
	if [ -d "$$d" ] ; then  echo "$(MAKE) clean  -C $$d  " ;   $(MAKE) clean  -C $$d  ; else echo ** missing "$$d" directory ;  fi; \
	done ;\
	fi; \
	fi 
