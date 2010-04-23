# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include ../headers-sparsesolver.inc
include ff-flags


# Downloading and compiling hypre
# ------------------------------

DIRPKG=../pkg
SRCDIR=hypre-$(HYPRE_VERSION)
PACKAGE=$(DIRPKG)/hypre-$(HYPRE_VERSION).tar.gz
SERVER=https://computation.llnl.gov/casc/hypre/download
PKGDIR=../pkg
HYPRE_VERSION=2.4.0b
INSTALL=../..
SRCDIR2=hypre\-$(HYPRE_VERSION)

hypre: FAIT$(WIN32DLLTARGET)

FAITwin32-dll-target:
	echo "On Pure Win32 (to hard to compile) "

FAIT:	
	$(MAKE) FAIRE WHERE
	touch FAIT


#====== 
#
#
# Remark :: To test hypre examples we need to set parameter
#  --with-blas-lib-dirs= --with-blas-libs= 
#  --with-lapack-lib-dirs --with-lapack-lib
#======

FAIRE: $(SRCDIR) 
	cd $(SRCDIR)/src && ./configure  CXX="$(CXX)" CC="$(CC)" CFLAGS="$(CFLAGS)" CPP="$(CC) -E" CXXFLAGS="$(CXXFLAGS)" F77="$(F77)" FFLAGS="$(FFLAGS)" FLIBS="$(FLIBS)" LDFLAGS="$(LDFLAGS)" --with-MPI-include="$(FFMPIDIR)" --with-MPI-lib-dirs="$(FFMPIDIRLIBS)" --with-MPI-libs="$(FFMPILIB2)" # --with-lapack --with-blas  #// to add lapack to test
	# the following line is to avoid compilation of local blas and local lapack of hypre
	# begin patch
	cd $(SRCDIR)/src/config/; mv Makefile.config Makefile-tmp.config; 
	sed -e 's;${abs_top_builddir}/download/hypre/${SRCDIR}/src/lapack; ;'\
	<$(SRCDIR)/src/config/Makefile-tmp.config\
	>$(SRCDIR)/src/config/Makefile.config
	cd $(SRCDIR)/src/config/; mv Makefile.config Makefile-tmp.config; 
	sed -e 's;${abs_top_builddir}/download/hypre/${SRCDIR}/src/blas; ;'\
	<$(SRCDIR)/src/config/Makefile-tmp.config\
	>$(SRCDIR)/src/config/Makefile.config
	# end patch
	cd $(SRCDIR)/src && $(MAKE)
	mkdir -p $(SRCDIR)/$(INSTALL)/include/hypre
	mkdir -p $(SRCDIR)/$(INSTALL)/lib/hypre
	cp $(SRCDIR)/src/hypre/include/* $(SRCDIR)/$(INSTALL)/include/hypre
	cp $(SRCDIR)/src/hypre/lib/* $(SRCDIR)/$(INSTALL)/lib/hypre
	chmod 644 $(SRCDIR)/$(INSTALL)/include/hypre/HYPRE.h
	rm -rf $(SRCDIR)/src/hypre
WHERE:
	echo hypre LD '-L@DIR@/lib/hypre  -lHYPRE'  >$(SRCDIR)/$(INSTALL)/lib/WHERE.hypre
	echo hypre INCLUDE '-I@DIR@/include/hypre'>> $(SRCDIR)/$(INSTALL)/lib/WHERE.hypre


$(SRCDIR): $(PACKAGE) 
	tar xvzf $(PACKAGE)
	touch $(SRCDIR)

$(PACKAGE):
	-mkdir $(DIRPKG)
	cd $(PKGDIR); $(WGET)  $(SERVER)/`basename $(PACKAGE)`

clean: clean-local
	-rm ff-flags
clean-local:
	-rm -rf $(SRCDIR) FAIT 

ff-flags: ../Makefile
	grep 'abs_top_builddir *=' ../Makefile >> ff-flags
	grep 'CXX *=' ../Makefile >> ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'FLIBS *=' ../Makefile >> ff-flags
	grep 'CXXFLAGS *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'F77 *=' ../Makefile >> ff-flags
	grep 'FFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'FLIBS *=' ../Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags
	grep 'WIN32DLLTARGET *=' ../Makefile >>ff-flags

