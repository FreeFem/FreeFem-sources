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



#====== 
#
#
# Remark :: To test hypre examples we need to set parameter
#  --with-blas-lib-dirs= --with-blas-libs= 
#  --with-lapack-lib-dirs --with-lapack-lib
#======

FAIRE: $(SRCDIR)/FAIT

$(SRCDIR)/FAIT:$(SRCDIR)
	cd $(SRCDIR)/src && ./configure  CXX="$(CXX)" CC="$(CC)" CFLAGS="$(CFLAGS)" CPP="$(CC) -E" CXXFLAGS="$(CXXFLAGS)" F77="$(F77)" FFLAGS="$(FFLAGS)" FLIBS="$(FLIBS)" LDFLAGS="$(LDFLAGS)" --with-MPI-include="$(FFMPIDIR)" --with-MPI-lib-dirs="$(FFMPIDIRLIBS)" --with-MPI-libs="$(FFMPILIB2)" --prefix="${abs_top_builddir}/download"  --includedir="${abs_top_builddir}/download"/include/hypre  --libdir="${abs_top_builddir}/download"/lib/hypre 
# --with-lapack --with-blas  #// to add lapack to test
	# the following line is to avoid compilation of local blas and local lapack of hypre
	# begin patch
#	cd $(SRCDIR)/src/config/; mv Makefile.config Makefile-tmp.config; 
#	sed -e 's;${abs_top_builddir}/download/hypre/${SRCDIR}/src/lapack; ;'\
#	<$(SRCDIR)/src/config/Makefile-tmp.config\
#	>$(SRCDIR)/src/config/Makefile.config
#	cd $(SRCDIR)/src/config/; mv Makefile.config Makefile-tmp.config; 
#	sed -e 's;${abs_top_builddir}/download/hypre/${SRCDIR}/src/blas; ;'\
#	<$(SRCDIR)/src/config/Makefile-tmp.config\
#	>$(SRCDIR)/src/config/Makefile.config
#	# end patch
	cd $(SRCDIR)/src && $(MAKE)
	cd $(SRCDIR)/src && $(MAKE) install
	touch $(SRCDIR)/FAIT
install: 
	cd $(SRCDIR)/src && $(MAKE) install
WHERE:
	if [ -f $(SRCDIR)/FAIT ] ; then  \
	make install ; \
	echo hypre LD '-L@DIR@/lib/hypre  -lHYPRE'  >$(SRCDIR)/$(INSTALL)/lib/WHERE.hypre ;\
	echo hypre INCLUDE '-I@DIR@/include/hypre'>> $(SRCDIR)/$(INSTALL)/lib/WHERE.hypre ;\
	fi


$(SRCDIR): $(PACKAGE) 
	tar xvzf $(PACKAGE)
	touch $(SRCDIR)

$(PACKAGE):
	-mkdir $(DIRPKG)
	cd $(PKGDIR); $(WGET)   $(SERVER)/`basename $(PACKAGE)`

clean: clean-local
	-rm ff-flags
clean-local:
	-rm -rf $(SRCDIR) FAIT 

ff-flags: ../Makefile ff-flags.inc 
	../../config.status  --file="ff-flags:ff-flags.inc"

