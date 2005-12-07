# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_UMFPACK)
EXTRA_DIST=UMFPACK_Make.m4

# Downloading and compiling UMFPACK
# ---------------------------------

UMFPACKLIB=UMFPACKv4.4/UMFPACK/Lib/libumfpack.a

umfpack: $(UMFPACKLIB)

# The 'lib' goal is replaced with 'libb' to avoid problems with the
# existing 'Lib' subdirectory on case-insensitive file systems
# the lib depend of the Makefile to force the reconstruction 
# if the parameter change

$(UMFPACKLIB): UMFPACKv4.4.tar.gz UMFPACK_Make.m4  Makefile UMFPACKv4.4/UMFPACK/Makefile2 UMFPACKv4.4/UMFPACK/Makefile2
	cd UMFPACKv4.4/UMFPACK && make -f Makefile2 libb		
UMFPACKv4.4:
	tar xzf UMFPACKv4.4.tar.gz
UMFPACKv4.4/AMD/Makefile2 UMFPACKv4.4/UMFPACK/Makefile2: UMFPACKv4.4  
	sed 's/lib:/libb:/' < UMFPACKv4.4/UMFPACK/`basename $@ 2`  >$@

UMFPACK_Make.m4: Makefile UMFPACKv4.4
	m4   -DFF_CC="$(CC)" \
             -DFF_CFLAGS="@CFLAGS@ @BLASINC@" \
             -DFF_LIB="@LIBS@ @BLASLIB@" \
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

clean-local:
	-rm -r UMFPACKv4.4.tar.gz UMFPACKv4.4
