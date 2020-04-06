############################################################################
# This file is part of FreeFEM.                                            #
#                                                                          #
# FreeFEM is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFEM is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          #
############################################################################
# SUMMARY : Downloading and compiling extra libraries for Umfpack
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

# see [[file:../../configure.ac::DOWNLOAD_UMFPACK]]
all-local: $(DOWNLOAD_UMFPACK)
EXTRA_DIST= SuiteSparse_config.mk.in 

# Downloading and compiling UMFPACK
# ---------------------------------

DIRPKG=../pkg
SUITESPARSE_TGZ=$(DIRPKG)/SuiteSparse-4.4.4.tar.gz

UMFPACKLIB=../lib/libumfpack.a 
AMDLIB=../lib/libamd.a
CAMDLIB=../lib/libcamd.a
CHOLMODLIB=../lib/libcholmod.a
COLAMDLIB=../lib/libcolamd.a
SUITESPARSECONFIGLIB=../lib/libsuitesparseconfig.a

# ALH - 4/9/13 - parallel make crashes ("pipe from processes is a directory"?). make sure that everything is run
# sequentially

umfpack:
	$(MAKE) MAKEFLAGS= $(AMDLIB)
	$(MAKE) MAKEFLAGS= $(UMFPACKLIB)
	$(MAKE) MAKEFLAGS= $(SUITESPARSECONFIGLIB)
#
# ALH - 19/9/13 - cholmod part built sequentially
#
	$(MAKE) MAKEFLAGS= $(CAMDLIB)
	$(MAKE) MAKEFLAGS= $(COLAMDLIB)
	$(MAKE) MAKEFLAGS= $(AMDLIB)
	$(MAKE) MAKEFLAGS= $(CHOLMODLIB)
	$(MAKE) MAKEFLAGS= $(SUITESPARSECONFIGLIB) 

# FFCS: SuiteSparse/*/Source are not able to compile in parallel from scratch ("pipe from processes is a
# directory"?). But specifying "make -j 1" is not enough (another error pops up). Use "$(MAKE) MAKEFLAGS=" instead
# (using "make MAKEFLAGS=" still produces an error on Cygwin).

$(SUITESPARSECONFIGLIB): SuiteSparse/FF 
	cd SuiteSparse/SuiteSparse_config && $(MAKE) MAKEFLAGS=
	mkdir -p ../include ../lib
	cp -f SuiteSparse/SuiteSparse_config/SuiteSparse_config.h  ../include
	cp  SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a ../lib/libsuitesparseconfig.a
	$(RANLIB) ../lib/libsuitesparseconfig.a

$(UMFPACKLIB): SuiteSparse/FF 
	cd SuiteSparse/UMFPACK/Lib && $(MAKE) MAKEFLAGS=
	mkdir -p ../include ../lib
	cp -f  SuiteSparse/UMFPACK/Include/*.h ../include
	cp -f SuiteSparse/SuiteSparse_config/SuiteSparse_config.h  ../include
	cp  SuiteSparse/UMFPACK/Lib/libumfpack.a ../lib/libumfpack.a
	$(RANLIB) ../lib/libumfpack.a

$(AMDLIB): SuiteSparse/FF 
	cd SuiteSparse/AMD/Lib && $(MAKE) MAKEFLAGS=
	mkdir -p ../include ../lib
	cp -f  SuiteSparse/AMD/Include/*.h ../include
	cp  SuiteSparse/AMD/Lib/libamd.a ../lib/libamd.a
	$(RANLIB) ../lib/libamd.a

$(CAMDLIB): SuiteSparse/FF 
	cd SuiteSparse/CAMD/Lib && $(MAKE) MAKEFLAGS=
	mkdir -p ../include ../lib
	cp -f  SuiteSparse/CAMD/Include/*.h ../include
	cp  SuiteSparse/CAMD/Lib/libcamd.a ../lib/libcamd.a
	$(RANLIB) ../lib/libcamd.a

$(CHOLMODLIB): SuiteSparse/FF 
	cd SuiteSparse/CHOLMOD/Lib && $(MAKE) MAKEFLAGS=
	mkdir -p ../include ../lib
	cp -f  SuiteSparse/CHOLMOD/Include/*.h ../include
	cp  SuiteSparse/CHOLMOD/Lib/libcholmod.a ../lib/libcholmod.a
	$(RANLIB) ../lib/libcholmod.a

$(COLAMDLIB): SuiteSparse/FF 
	cd SuiteSparse/COLAMD/Lib && $(MAKE) MAKEFLAGS=
	mkdir -p ../include ../lib
	cp -f  SuiteSparse/COLAMD/Include/*.h ../include
	cp  SuiteSparse/COLAMD/Lib/libcolamd.a ../lib/libcolamd.a
	$(RANLIB) ../lib/libcolamd.a

UMFPACKv4.4:
	tar xzf UMFPACKv4.4.tar.gz

# The 'lib' goal is replaced with 'libb' to avoid problems with the existing 'Lib' subdirectory on case-insensitive file
# systems the lib depend of the Makefile to force the reconstruction if the parameter change

UMFPACKv4.4/AMD/Makefile2 UMFPACKv4.4/UMFPACK/Makefile2: UMFPACKv4.4  
	sed 's/lib:/libb:/' < UMFPACKv4.4/UMFPACK/`basename $@ 2`  >$@

SuiteSparse/DATE:$(SUITESPARSE_TGZ)
	tar zxf $(SUITESPARSE_TGZ)
	touch SuiteSparse/DATE

SuiteSparse/FF:SuiteSparse/DATE  SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk
	touch SuiteSparse/FF

SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk:SuiteSparse/DATE Makefile SuiteSparse_config.mk.in
	../../config.status  --file="SuiteSparse_config.mk:SuiteSparse_config.mk.in"
	-if  diff SuiteSparse_config.mk SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk  2>&1 >/dev/null; then \
	    echo " same flags => no recompilation !  " ; \
	else \
	   echo "  recompile umfpack (some flags change) => clean umfpack colmod amd " ;\
	   mkdir -p SuiteSparse/SuiteSparse_config ;\
	   cp SuiteSparse_config.mk SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk ; \
	   (cd SuiteSparse/UMFPACK && make clean); \
	   (cd SuiteSparse/CHOLMOD && make clean); \
	   (cd SuiteSparse/COLAMD && make clean); \
	   (cd SuiteSparse/CAMD && make clean); \
	   (cd SuiteSparse/CAMD && make clean); \
	   (cd SuiteSparse/AMD && make clean); \
	   find . -name '*.exe'|xargs rm;  \
	fi; 
	touch SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk

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
	../getall -o SuiteSparse -a

clean-local:
	-rm $(AMDLIB) $(CAMDLIB) $(CHOLMODLIB) $(COLAMDLIB) $(SUITESPARSECONFIGLIB)
	-rm ../include/amd*.h 
	-rm ../include/umfpack*.h 
	-rm ../include/colmod*.h 
	-rm ../include/camd*.h 
	-rm ../include/amd*.h
	-rm ../include/SuiteSparse_config.h 
	-rm -rf UMFPACKv4.?.tar.gz UMFPACKv4.?
	-rm SuiteSparse*gz 
	-rm -rf SuiteSparse
#	-rm ../pkg/SuiteSparse-*
	-rm SuiteSparse_config.mk

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
