# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_UMFPACK)

# Downloading and compiling UMFPACK
# ---------------------------------

UMFPACKLIB=UMFPACKv4.3/UMFPACK/Lib/libumfpack.a

umfpack: $(UMFPACKLIB)

# The 'lib' goal is replaced with 'libb' to avoid problems with
# existing 'Lib' subdir on case-insensitive file systems

$(UMFPACKLIB): UMFPACKv4.3.tar.gz
	tar xzf UMFPACKv4.3.tar.gz
	cp AMD_Make.include UMFPACKv4.3/AMD/Make/Make.include
	sed 's/lib:/libb:/' < UMFPACKv4.3/AMD/Makefile \
		> UMFPACKv4.3/AMD/Makefile2
	cd UMFPACKv4.3/AMD && make -f Makefile2 libb \
		CFLAGS="$(CFLAGS)" F77="$(F77)" \
		F77FLAGS="$(FFLAGS)" CC="$(CC)" LIB="$(BLASLIB) -lm"
	cp UMFPACK_Make.include UMFPACKv4.3/UMFPACK/Make/Make.include
	sed 's/lib:/libb:/' < UMFPACKv4.3/UMFPACK/Makefile \
		> UMFPACKv4.3/UMFPACK/Makefile2
	cd UMFPACKv4.3/UMFPACK && make -f Makefile2 libb \
		CFLAGS="$(CFLAGS)" F77="$(F77)" \
		F77FLAGS="$(FFLAGS)" CC="$(CC)" LIB="$(BLASLIB) -lm"

UMFPACKv4.3.tar.gz:
	wget -N http://www.cise.ufl.edu/research/sparse/umfpack/v4.3/UMFPACKv4.3.tar.gz

clean-local:
	-rm -r UMFPACKv4.3.tar.gz UMFPACKv4.3
