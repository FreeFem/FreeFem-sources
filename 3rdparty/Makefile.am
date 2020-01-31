# Downloading and compiling extra libraries
# -----------------------------------------

# $Id: Makefile.am,v 1.16 2010/05/06 21:20:38 hecht Exp $

SUBDIRS=blas arpack umfpack
EXTRA_DIST= \
./nlopt/Make.inc ./nlopt/Makefile \
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
./metis/Makefile \
./metis/Makefile-metis.in \
./metis/metis-5.1.patch \
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
./mumps-seq/Makefile-mumps-5.0.2.inc \
./mumps-seq/Makefile \
./mumps/Makefile-mumps-5.0.2.inc \
./parmetis/Makefile-parmetis.in \
./parmetis/makefile \
./parmetis/parmetis-4.0.3.patch \
./scalapack/Makefile \
./scalapack/SLmake-scalapack.inc \
./scotch/Makefile \
./scotch/Makefile-scotch.inc \
./scotch/scotch_6.0.4.patch \
./superlu/Makefile \
./superlu/make.inc \
./tetgen/Makefile \
./tetgen/patches.win64 \
./mmg/Makefile \
./parmmg/Makefile \
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
ff-petsc/Makefile \
ff-petsc/Makefile-PETSc.inc \
ff-petsc/Makefile.complex \
ff-petsc/petsc-3.11.2.patch \
getall

# FFCS: See [[file:../../../configure.ac::tools_problems_all_platforms]] for reasons why some tools may be deactivated


MPI_SOFT= @TOOL_COMPILE_scalapack@ @TOOL_COMPILE_parmetis@ @TOOL_COMPILE_parmmg@ @TOOL_COMPILE_mumps@

LIST_SOFT= @DOWNLOAD_FFTW@ @TOOL_COMPILE_tetgen@ @TOOL_COMPILE_metis@ @TOOL_COMPILE_superlu@ \
	@TOOL_COMPILE_scotch@ @TOOL_COMPILE_mshmet@			\
	@TOOL_COMPILE_yams@ @TOOL_COMPILE_mmg3d@ @TOOL_COMPILE_mmg@ @TOOL_COMPILE_gmm@	\
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
# ALH - /3rdparty/yams and /3rdparty/mshmet need /src/libMesh/libmesh.a but /src is compiled after /download, so we
# need to compile it now

lib/libMesh.a:lib include
	cd ../src/libMesh && $(MAKE) $(AM_MAKEFLAGS)
	test -f ../src/libMesh/libMesh.a
	mkdir -p include/libMesh
	cp ../src/libMesh/*h  include/libMesh
	echo libMesh  LD -L@DIR@/lib -lMesh  > lib/WHERE.libMesh
	echo libMesh INCLUDE -I@DIR@/include/libMesh >> lib/WHERE.libMesh
	cp ../src/libMesh/libMesh.a lib/libMesh.a

all-local:bin lib include lib/libMesh.a  pkg  $(DOWNLOADCOMPILE)  install-other

install-other: tag-install-other

tag-install-other:
	if test -n "$(TOOL_COMPILE_hpddm)" ; then  $(MAKE) install-hpddm; else true; fi
	if test -n "$(TOOL_COMPILE_htool)" ; then  $(MAKE) install-htool; else true; fi
	if test -n "$(TOOL_COMPILE_bemtool)" ; then  $(MAKE) install-bemtool; else true; fi
	if test -n "$(TOOL_COMPILE_boost)" ; then  $(MAKE) install-boost; else true; fi
	touch tag-install-other

WHERE-OTHER: lib/WHERE.hpddm lib/WHERE.htool lib/WHERE.bemtool lib/WHERE.boost
WHERE-LD: tag-compile-pkg  WHERE-OTHER
	touch ../plugin/seq/WHERE_LIBRARY-config ../plugin/seq/WHERE_LIBRARY
	-grep LD ../plugin/seq/WHERE_LIBRARY ../plugin/seq/WHERE_LIBRARY-config >WHERE-LD

# BEGIN HPDDM
install-hpddm: ./getall
	 if test -n "$(TOOL_COMPILE_hpddm)"; then   ./getall -o hpddm -a;  $(MAKE) include/hpddm/done.tag  lib/WHERE.hpddm; else true;fi
reinstall-hpddm:./getall
	-rm lib/WHERE.hpddm
	-rm ./pkg/hpddm.zip
	-rm -rf include/hpddm
	-test -n "$(TOOL_COMPILE_hpddm)" && $(MAKE) install-hpddm || true
include/hpddm/done.tag:./pkg/hpddm.zip
	-if test -f ./pkg/hpddm.zip ; then \
          cd include;rm -rf hpddm hpddm-*;  unzip ../pkg/hpddm.zip ; mv  hpddm-* hpddm ; touch hpddm/done.tag; \
         else true;  fi
lib/WHERE.hpddm: ./getall
	if test -d include/hpddm/include ; then \
	echo hpddm LD -L@DIR@/lib > $@ ;\
	echo hpddm  INCLUDE -I@DIR@/include/hpddm/include >> $@ ;\
	elif test -d include/hpddm/src ; then \
	echo hpddm LD -L@DIR@/lib > $@ ;\
	echo hpddm  INCLUDE -I@DIR@/include/hpddm/src >> $@ ;\
	else true; 	fi
#end HPDDM

# begin HTOOL
install-htool: ./getall
	 if test -n "$(TOOL_COMPILE_htool)"; then   ./getall -o htool -a;  $(MAKE) include/htool/done.tag  lib/WHERE.htool; else true;fi
reinstall-htool:./getall
	-rm lib/WHERE.htool
	-rm ./pkg/htool.zip
	-rm -rf include/htool
	-test -n "$(TOOL_COMPILE_htool)" && $(MAKE) install-htool || true
include/htool/done.tag:./pkg/htool.zip
	-if test -f ./pkg/htool.zip ; then \
          cd include;rm -rf htool htool-*;  unzip ../pkg/htool.zip ; mv  htool-* htool ; touch htool/done.tag; \
         else true;  fi
lib/WHERE.htool: ./getall
	if test -d include/htool/include ; then \
	echo htool LD -L@DIR@/lib > $@ ;\
	echo htool  INCLUDE -I@DIR@/include/htool/include >> $@ ;\
	else true; 	fi
# end HTOOL
# begin BEMTOOL
install-bemtool: ./getall
	 if test -n "$(TOOL_COMPILE_bemtool)"; then   ./getall -o bemtool -a;  $(MAKE) include/bemtool/done.tag  lib/WHERE.bemtool; else true;fi
reinstall-bemtool:./getall
	-rm lib/WHERE.bemtool
	-rm ./pkg/bemtool.zip
	-rm -rf include/BemTool
	-test -n "$(TOOL_COMPILE_bemtool)" && $(MAKE) install-bemtool || true
include/bemtool/done.tag:./pkg/bemtool.zip
	-if test -f ./pkg/bemtool.zip ; then \
          cd include;rm -rf BemTool BemTool-*;  unzip ../pkg/bemtool.zip ; mv BemTool-* BemTool ; touch BemTool/done.tag; \
         else true;  fi
lib/WHERE.bemtool: ./getall
	if test -d include/BemTool/ ; then \
	echo bemtool LD -L@DIR@/lib > $@ ;\
	echo bemtool  INCLUDE -I@DIR@/include/BemTool/ >> $@ ;\
	else true; 	fi
# end BEMTOOL

# begin BOOST
install-boost: ./getall
	 if test -n "$(TOOL_COMPILE_boost)"; then   ./getall -o Boost -a;  $(MAKE) boost/done.tag  lib/WHERE.boost; else true;fi
reinstall-boost:./getall
	-rm lib/WHERE.boost
	-rm ./pkg/boost.tar.gz
	-rm -rf boost
	-test -n "$(TOOL_COMPILE_boost)" && $(MAKE) install-boost || true
boost/done.tag:./pkg/boost.tar.gz
	if test -f ./pkg/boost.tar.gz ; then \
	rm -rf boost boost_*;   tar xvzf   pkg/boost.tar.gz  ; mv boost_for_bemtool boost; \
	else true;  fi	
	touch boost/done.tag
lib/WHERE.boost: ./getall
	if test -d boost/include ; then \
	echo boost  INCLUDE -I@DIR@/boost/include >> $@ ;\
	else true; 	fi
# end BOOST

#  to reinstall mpi under window afer clean 
reinstall-msmpi:
	-if test -f "$$MSMPI_INC"/mpif.h   ; then \
		echo " copy msmpi  in 3rdparty form $$MSMPI_INC and $$MSMPI_LIB64 or $$MSMPI_LIB32" ; \
		mkdir -p include/msmpi ;\
		mkdir -p lib/msmpi ;\
		cp "$$MSMPI_INC"/*.h include/msmpi ;\
		sed  's/INT_PTR_KIND()/@SIZEOF_PTR@/' <"$$MSMPI_INC"/mpif.h >include/msmpi/mpif.h ;\
		grep KIND include/msmpi/mpif.h; \
		test "@SIZEOF_PTR@" -eq 8 && cp "$$MSMPI_INC"/x64/*.h include/msmpi &&  cp "$$MSMPI_LIB64"/*.lib lib/msmpi ;\
		test "@SIZEOF_PTR@" -eq 4 && cp "$$MSMPI_INC"/x86/*.h include/msmpi &&  cp "$$MSMPI_LIB32"/*.lib lib/msmpi ;\
	fi
# FFCS: need to stop at the first error to make sure that all libraries are correctly compiled
compile-dir:
	@echo "\n\n ****** $(COMPILEDIR) ****** \n\n";
	@if [ 0 -eq `egrep ':$(COMPILEDIR) *LD' WHERE-LD | wc -l` ] ;then \
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
	$(MAKE) install-hpddm install-htool install-bemtool install-boost
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
	-rm -rf	 tag-* include lib bin boost
	-mkdir   include lib bin
	-rm */FAIT */FAIRE
# 	FFCS - make sure that all directories are cleaned. Thisis especially important under Windows because there is no
# 	compilation dependencies control there (see
# 	[[file:c:/cygwin/home/alh/ffcs/dist/configure.ac::dependency_tracking]])
	for d in $(LIST_SOFT) $(MPI_SOFT) ; do $(MAKE) clean -C $$d ; done
