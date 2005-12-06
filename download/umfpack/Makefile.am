# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_UMFPACK)
EXTRA_DIST=AMD_Make.include UMFPACK_Make.include

# Downloading and compiling UMFPACK
# ---------------------------------

UMFPACKLIB=UMFPACKv4.4/UMFPACK/Lib/libumfpack.a

umfpack: $(UMFPACKLIB)

# The 'lib' goal is replaced with 'libb' to avoid problems with the
# existing 'Lib' subdirectory on case-insensitive file systems

$(UMFPACKLIB): UMFPACKv4.4.tar.gz UMFPACK_Make.m4 UMFPACKv4.4
	sed 's/lib:/libb:/' < UMFPACKv4.4/UMFPACK/Makefile \
		> UMFPACKv4.4/UMFPACK/Makefile2
	cd UMFPACKv4.4/UMFPACK && make -f Makefile2 libb		
UMFPACKv4.4:
	tar xzf UMFPACKv4.4.tar.gz
	sed 's/lib:/libb:/' < UMFPACKv4.4/AMD/Makefile \
		> UMFPACKv4.4/AMD/Makefile2
UMFPACK_Make.m4: Makefile UMFPACKv4.4
	m4   -DFF_CC="$(CC)" -DFF_CFLAGS="@CFLAGS@" -DFF_LIB="@LIBS@" \
		 -DFF_CONFIG="@FF_UMFPACK_CONFIG@" \
              UMFPACK_Make.m4 >UMFPACKv4.4/UMFPACK/Make/Make.include
	-find . -name '*.[ao]'|xargs rm
UMFPACKv4.4.tar.gz:
	@WGET@ -N http://www.cise.ufl.edu/research/sparse/umfpack/v4.4/UMFPACKv4.4.tar.gz

clean-local:
	-rm -r UMFPACKv4.4.tar.gz UMFPACKv4.4
