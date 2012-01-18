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
FFTW_VERSION=3.2.1
INSTALL=../..
THIS=fftw3

fftw:FAIT
#   $(WIN32DLLTARGET)  

FAITwin32-dll-target:
	echo "On Pure Win32 (to hard to compile) "

FAIT:	
	$(MAKE) FAIRE
	touch FAIT 

FAIRE: $(SRCDIR) 
	cd $(SRCDIR) && ./configure --disable-dependency-tracking  --disable-fortran  --prefix=`pwd`/$(INSTALL) CXX="$(CXX)" CC="$(CC)"  CFLAGS="$(CFLAGS)" CPP='gcc -E' CXXFLAGS="$(CXXFLAGS)"
	cd $(SRCDIR) && make install
	touch FAIRE
install:  WHERE
	cd $(SRCDIR) && make install
WHERE: 
	-if [ FAIT ] ; then \
	echo $(THIS)  LD -L@DIR@/lib -l$(THIS)  >$(SRCDIR)/$(INSTALL)/lib/WHERE.$(THIS) ;\
	echo $(THIS) INCLUDE -I@DIR@/include   >> $(SRCDIR)/$(INSTALL)/lib/WHERE.$(THIS) ;\
	fi

$(SRCDIR): $(PKGDIR)/$(PACKAGE)
	gunzip -c   $^ | tar xvf - 

$(PKGDIR)/$(PACKAGE):
	cd $(PKGDIR); $(WGET)  $(SERVER)/$(PACKAGE)

clean: clean-local

clean-local:
	-rm -rf fftw-* FAIT FAIRE

cxxflags: ../Makefile
	grep 'CXX *=' ../Makefile >cxxflags
	grep 'CC *=' ../Makefile >>cxxflags
	grep 'CXXFLAGS *=' ../Makefile >>cxxflags
	grep 'CFLAGS *=' ../Makefile >>cxxflags
	grep 'WGET *=' ../Makefile >>cxxflags
	grep 'WIN32DLLTARGET *=' ../Makefile >>cxxflags

