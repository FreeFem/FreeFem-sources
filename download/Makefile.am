# Downloading and compiling extra libraries
# -----------------------------------------

# $Id: Makefile.am,v 1.16 2010/05/06 21:20:38 hecht Exp $

SUBDIRS=blas arpack umfpack  
EXTRA_DIST= \
./nlopt/Make.inc ./nlopt/Makefile \
./blacs/BLACS.patch \
./blacs/BLACS_gridinit_.c-return-values.patch \
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
./metis/metis-4.0_main_return.patch \
./mmg3d/Makefile \
./mmg3d/patch-mmg3dv4.diff \
./mshmet/Makefile \
./mshmet/Makefile-mshmet.inc \
./mshmet/mshmet.2011.03.06.patch \
./mshmet/mshmet.2012.04.25_i586.patch \
./mshmet/mshmetlib-internal.h \
./mshmet/mshmetlib.c \
./mshmet/mshmetlib.h \
./mumps/Makefile \
./mumps-seq/Makefile-mumps-4.10.0.inc \
./mumps-seq/Makefile \
./mumps/Makefile-mumps-4.10.0.inc \
./mumps/MUMPS_4.10.0.patch \
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
./scotch/Makefile.patch \
./scotch/scotch_5.1_esmumps.patch \
./superlu/Makefile \
./superlu/make.inc \
./superludist/Makefile \
./superludist/make-superlu.inc \
./superludist/superludist_2.3.patch \
./superludist_3.0-printf.patch \
./superludist_3.0-cast_warning.patch \
./superludist_3.0-return_values.patch \
./superludist_3.0-operation_undefined.patch \
./tetgen/Makefile \
./tetgen/tetgen1.4.2.patch \
./yams/Makefile \
./yams/freeyams.2012.02.05.patch \
./yams/freeyams.2012.02.05-return-values.patch \
./yams/makefile-yams.inc \
./yams/yamslib.c \
./yams/yamslib.h \
./yams/yamslib_internal.h \
ipopt/Makefile	\
ipopt/Makefile.inc.in \
ipopt/patch-IpBlas 

# FFCS: See [[file:../../../configure.ac::tools_problems_all_platforms]] for reasons why some tools may be deactivated

# <<no_more_blacs>> blacs is include in  scalapack 2.0.2 

MPI_SOFT=parmetis scalapack @FFCS_COMPILE_scotch@ @FFCS_COMPILE_superludist@ @FFCS_COMPILE_mumps@	\
	@FFCS_COMPILE_pastix@ @FFCS_COMPILE_hypre@ @FFCS_COMPILE_hips@ @FFCS_COMPILE_parms@
LIST_SOFT=tetgen superlu fftw metis @FFCS_COMPILE_mshmet@ @FFCS_COMPILE_yams@ @FFCS_COMPILE_mmg3d@ gmm nlopt	\
	@FFCS_COMPILE_mumps_seq@ @FFCS_COMPILE_ipopt@
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

# FFCS: need to stop at the first error to make sure that all libraries are correctly compiled
compile-dir: 
	@echo "\n\n ****** $(COMPILEDIR) ****** \n\n";
	@if [ 0 -eq `egrep ':$(COMPILEDIR)' WHERE-LD | wc -l` ] ;then \
	  cd $(COMPILEDIR) && $(MAKE) $(DIRTARGET) ; \
	else \
	 echo $(COMPILEDIR) is in WHERE- files ;\
	fi 
compile-pkg: tag-compile-pkg WHERE-LD

# FFCS: need to stop at the first error to make sure that all libraries are correctly compiled
tag-compile-pkg: bin lib include pkg FORCE
	@if [ -n "$(WGET)" ] ; then \
	for d in $(LIST_SOFT) ; do $(MAKE) compile-dir COMPILEDIR=$$d || exit 1;done ;\
	if [ -n "$(MPICC)" ] ; then \
	for d in $(MPI_SOFT) ; do 	$(MAKE) compile-dir COMPILEDIR=$$d || exit 1; done;\
	fi;fi
	touch tag-compile-pkg
FORCE:

re-install: 
	$(MAKE) compile-pkg DIRTARGET=install
WHERE:
	$(MAKE) compile-pkg DIRTARGET=WHERE

clean-local:
	-rm -rf	 tag-* include lib bin 
	-mkdir   include lib bin 	
	-rm */FAIT */FAIRE 
# 	FFCS - make sure that all directories are cleaned. Thisis especially important under Windows because there is no
# 	compilation dependencies control there (see
# 	[[file:c:/cygwin/home/alh/ffcs/dist/configure.ac::dependency_tracking]])
	for d in $(LIST_SOFT) $(MPI_SOFT) ; do $(MAKE) clean -C $$d ; done
