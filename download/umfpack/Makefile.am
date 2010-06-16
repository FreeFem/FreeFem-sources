# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_UMFPACK)
EXTRA_DIST=UFconfig_mk.m4

# Downloading and compiling UMFPACK
# ---------------------------------

DIRPKG=../pkg
SUITESPARSE_TGZ=$(DIRPKG)/SuiteSparse-2.4.0.tar.gz

UMFPACKLIB=../lib/libumfpack.a 
AMDLIB=../lib/libamd.a

umfpack:$(UMFPACKLIB) $(AMDLIB)

# The 'lib' goal is replaced with 'libb' to avoid problems with the
# existing 'Lib' subdirectory on case-insensitive file systems
# the lib depend of the Makefile to force the reconstruction 
# if the parameter change

$(UMFPACKLIB): SuiteSparse/FF 
	cd SuiteSparse/UMFPACK/Source && make 
	-mkdir ../include ../lib
	cp -f  SuiteSparse/UMFPACK/Include/*.h ../include
	cp -f SuiteSparse/UFconfig/UFconfig.h  ../include
	cp  SuiteSparse/UMFPACK/Lib/libumfpack.a ../lib/libumfpack.a
	$(RANLIB) ../lib/libumfpack.a
$(AMDLIB): SuiteSparse/FF 
	cd SuiteSparse/AMD/Source && make 
	-mkdir ../include ../lib
	cp -f  SuiteSparse/AMD/Include/*.h ../include
	cp  SuiteSparse/AMD/Lib/libamd.a ../lib/libamd.a
	$(RANLIB) ../lib/libamd.a
UMFPACKv4.4:
	tar xzf UMFPACKv4.4.tar.gz
UMFPACKv4.4/AMD/Makefile2 UMFPACKv4.4/UMFPACK/Makefile2: UMFPACKv4.4  
	sed 's/lib:/libb:/' < UMFPACKv4.4/UMFPACK/`basename $@ 2`  >$@

SuiteSparse/DATE:$(SUITESPARSE_TGZ)
	tar zxvf $(SUITESPARSE_TGZ)
	touch SuiteSparse/DATE
SuiteSparse/FF:SuiteSparse/DATE  SuiteSparse/UFconfig/UFconfig.mk
	touch SuiteSparse/FF

SuiteSparse/UFconfig/UFconfig.mk:SuiteSparse/DATE Makefile UFconfig_mk.m4
	m4   -DFF_CC="$(CC)" \
             -DFF_CFLAGS="@CPPFLAGS@ @CFLAGS@ @BLASINC@" \
  	     -DFF_AR="@AR@" \
	     -DFF_ARFLAGS="@ARFLAGS@" \
	     -DFF_RANLIB="@RANLIB@" \
             -DFF_LIB="@BLASLIB@ @LIBS@" \
	     -DFF_CONFIG="@FF_UMFPACK_CONFIG@" \
              UFconfig_mk.m4 > UFconfig.mk
	-if  diff UFconfig.mk SuiteSparse/UFconfig/UFconfig.mk  2>&1 >/dev/null; then \
	    echo " same flags => no recompilation !  " ; \
	else \
	   echo "  recompile umfpack (some flags change) => clean umfpack" ;\
	   cp UFconfig.mk SuiteSparse/UFconfig/UFconfig.mk ; \
	   cd SuiteSparse/UMFPACK && make clean; \
	   cd ../AMD && make clean; \
	   cd ..;find . -name '*.exe'|xargs rm;  \
	fi; 

UMFPACK_Make.m4: Makefile UMFPACKv4.4
	m4   -DFF_CC="$(CC)" \
             -DFF_CFLAGS="@CPPFLAGS@ @CFLAGS@ @BLASINC@" \
             -DFF_LIB="@BLASLIB@ @LIBS@" \
	     -DFF_CONFIG="@FF_UMFPACK_CONFIG@" \
              UMFPACK_Make.m4 >Make.include
	-if  diff Make.include UMFPACKv4.4/UMFPACK/Make/Make.include 2>&1 >/dev/null; then \
	    echo " same flags => no recompilation !  " ; \
	else \
	   echo "  recompile umfpack (some flags change) => clean umfpack" ;\
	   cp Make.include  UMFPACKv4.4/UMFPACK/Make/Make.include; \
	   cd UMFPACKv4.4/UMFPACK && make clean; \
	   find . -name '*.exe'|xargs rm; \
	fi; 
UMFPACKv4.4.tar.gz:
	@WGET@ -N http://www.cise.ufl.edu/research/sparse/umfpack/v4.4/UMFPACKv4.4.tar.gz
$(SUITESPARSE_TGZ):
	cd $(DIRPKG);@WGET@ -N http://www.cise.ufl.edu/research/sparse/SuiteSparse/SuiteSparse-2.4.0.tar.gz

clean-local:
	-rm -rf UMFPACKv4.?.tar.gz UMFPACKv4.?
	-rm SuiteSparse*gz 
	-rm -rf SuiteSparse