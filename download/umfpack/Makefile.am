# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_UMFPACK)

# Downloading and compiling UMFPACK
# ---------------------------------

UMFPACKLIB=UMFPACKv4.3/UMFPACK/Lib/libumfpack.a

umfpack: $(UMFPACKLIB)

$(UMFPACKLIB): UMFPACKv4.3.tar.gz
	tar xzf UMFPACKv4.3.tar.gz
	cp AMD_Make.include UMFPACKv4.3/AMD/Make/Make.include
	cd UMFPACKv4.3/AMD && make lib CFLAGS="$(CFLAGS)" F77="$(F77)" \
		F77FLAGS="$(FFLAGS)" CC="$(CC)" LIB="$(LIBS)"
	cp UMFPACK_Make.include UMFPACKv4.3/UMFPACK/Make/Make.include
	cd UMFPACKv4.3/UMFPACK && make lib CFLAGS="$(CFLAGS)" F77="$(F77)" \
		F77FLAGS="$(FFLAGS)" CC="$(CC)" LIB="$(LIBS)"

UMFPACKv4.3.tar.gz:
	wget -N http://www.cise.ufl.edu/research/sparse/umfpack/v4.3/UMFPACKv4.3.tar.gz

clean-local:
	-rm -r UMFPACK*.tar.gz
