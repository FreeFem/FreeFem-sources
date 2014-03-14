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
./fftw/Makefile.am \
./fftw/Makefile.in \
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
./parmetis/makefile \
./parmetis/parmetis-3.1.1.patch \
./parms/Makefile \
./parms/makefile-parms.in \
./pastix/Makefile \
./pastix/config-pastix-complex.in \
./pastix/config-pastix-real.in \
./pastix/pastix_release_2200-blend.patch \
./pastix/patch-pastix_long_complex.h \
./pastix//all_macros.diff \
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
./superludist/superludist_3.0-printf.patch \
./superludist/superludist_3.0-cast_warning.patch \
./superludist/superludist_3.0-return_values.patch \
./superludist/superludist_3.0-operation_undefined.patch \
./tetgen/Makefile \
./tetgen/tetgen1.4.2.patch \
./tetgen/patches.win64 \
./yams/Makefile \
./yams/freeyams.2012.02.05.patch \
./yams/freeyams.2012.02.05-return-values.patch \
./yams/makefile-yams.inc \
./yams/yamslib.c \
./yams/yamslib.h \
./yams/yamslib_internal.h \
ipopt/Makefile	\
ipopt/Makefile.inc.in \
ipopt/patch-IpBlas  \
getall 

# FFCS: See [[file:../../../configure.ac::tools_problems_all_platforms]] for reasons why some tools may be deactivated

# <<no_more_blacs>> blacs is included in scalapack 2.0.2

MPI_SOFT=@TOOL_COMPILE_parmetis@ scalapack @TOOL_COMPILE_superludist@	\
	@TOOL_COMPILE_mumps@ @TOOL_COMPILE_pastix@			\
	@TOOL_COMPILE_hips@ @TOOL_COMPILE_parms@

LIST_SOFT=tetgen @TOOL_COMPILE_superlu@ fftw @TOOL_COMPILE_metis@	\
	@TOOL_COMPILE_scotch@ @TOOL_COMPILE_mshmet@			\
	@TOOL_COMPILE_yams@ @TOOL_COMPILE_mmg3d@ gmm			\
	@TOOL_COMPILE_nlopt@ @TOOL_COMPILE_mumps_seq@			\
	@TOOL_COMPILE_ipopt@

all-recursive: bin lib include pkg

lib: 
	mkdir lib
bin: 
	mkdir bin
include:
	mkdir  include
pkg:
	mkdir pkg

# ALH - /download/yams and /download/mshmet need /src/libMesh/libmesh.a but /src is compiled after /download, so we
# need to compile it now

lib/libMesh.a:lib include
	cd ../src/libMesh && $(MAKE) $(AM_MAKEFLAGS)
	test -f ../src/libMesh/libMesh.a
	mkdir -p include/libMesh
	cp ../src/libMesh/*h  include/libMesh
	echo libMesh  LD -L@DIR@/lib -lMesh  > lib/WHERE.libMesh
	echo libMesh INCLUDE -I@DIR@/include/libMesh >> lib/WHERE.libMesh
	cp ../src/libMesh/libMesh.a lib/libMesh.a

all-local:bin lib include lib/libMesh.a pkg  $(DOWNLOADCOMPILE) 

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
install-exec-local:	
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/lib
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/bin
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/include
	cp -rp lib $(DESTDIR)$(ff_prefix_dir)
	cp -rp include  $(DESTDIR)$(ff_prefix_dir)
	cp -rp bin  $(DESTDIR)$(ff_prefix_dir)

clean-local:
	-rm -rf	 tag-* include lib bin 
	-mkdir   include lib bin 	
	-rm */FAIT */FAIRE 
# 	FFCS - make sure that all directories are cleaned. Thisis especially important under Windows because there is no
# 	compilation dependencies control there (see
# 	[[file:c:/cygwin/home/alh/ffcs/dist/configure.ac::dependency_tracking]])
	for d in $(LIST_SOFT) $(MPI_SOFT) ; do $(MAKE) clean -C $$d ; done
