# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include cxxflags



# Downloading and compiling FFTW
# ------------------------------

# FFTW information
#http://www.fftw.org/fftw-3.2.tar.gz
SRCDIR=fftw-$(FFTW_VERSION)
PACKAGE=fftw-$(FFTW_VERSION).tar.gz
SERVER=http://www.fftw.org
PKGDIR=../pkg
FFTW_VERSION=3.2
INSTALL=../..


fftw:FAIT

FAIT:	
	$(MAKE) FAIRE
	touch FAIT

FAIRE: $(SRCDIR) 
	cd $(SRCDIR) && ./configure  --disable-fortran  --prefix=`pwd`/$(INSTALL) CXX="$(CXX)" CC="$(CC)"  CXXFLAGS="$(CXXFLAGS)"
	cd $(SRCDIR) && make
	cd $(SRCDIR) && make install
	touch FAIT

$(SRCDIR): $(PKGDIR)/$(PACKAGE)
	zcat  $^ | tar xvf - 

$(PKGDIR)/$(PACKAGE):
	cd $(PKGDIR); $(WGET)  $(SERVER)/$(PACKAGE)

clean-local:
	-rm -rf fftw-* FAIT 

cxxflags: ../Makefile
	grep 'CXX *=' ../Makefile >cxxflags
	grep 'CC *=' ../Makefile >>cxxflags
	grep 'CXXFLAGS *=' ../Makefile >>cxxflags
	grep 'WGET *=' ../Makefile >>cxxflags
