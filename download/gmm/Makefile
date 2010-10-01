# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include cxxflags



# Downloading and compiling FFTW
# ------------------------------
PKG=gmm
# $(PKG) information
SRCDIR=$(PKG)-$(PKG_VERSION)
PACKAGE=$(PKG)-$(PKG_VERSION).tar.gz
SERVER=http://download.gna.org/getfem/stable
PKGDIR=../pkg
PKG_VERSION=4.1
INSTALL=../..


$(PKG):FAIT WHERE

FAITwin32-dll-target:
	echo "On Pure Win32 (to hard to compile) "

FAIT: 
	$(MAKE) install  WHERE
	touch FAIT

FAIRE: $(SRCDIR) 
	cd $(SRCDIR) && ./configure --disable-dependency-tracking   --prefix=`pwd`/$(INSTALL) CXX="$(CXX)" CC="$(CC)"  CFLAGS="$(CFLAGS)" CPP='gcc -E' CXXFLAGS="$(CXXFLAGS)"
	cd $(SRCDIR) && make
	touch FAIRE
install: FAIRE 
	cd $(SRCDIR) && make install


$(SRCDIR): $(PKGDIR)/$(PACKAGE)
	gunzip -c   $^ | tar xvf - 

$(PKGDIR)/$(PACKAGE):
	cd $(PKGDIR); $(WGET)  $(SERVER)/$(PACKAGE)

clean: clean-local

clean-local:
	-rm -rf $(PKG)-* FAIT FAIRE

cxxflags: ../Makefile
	grep 'CXX *=' ../Makefile >cxxflags
	grep 'CC *=' ../Makefile >>cxxflags
	grep 'CXXFLAGS *=' ../Makefile >>cxxflags
	grep 'CFLAGS *=' ../Makefile >>cxxflags
	grep 'WGET *=' ../Makefile >>cxxflags
	grep 'WIN32DLLTARGET *=' ../Makefile >>cxxflags

WHERE: 
	-if [ -f FAIT ] ; then \
	echo gmm INCLUDE '-I@DIR@/include'> $(SRCDIR)/$(INSTALL)/lib/WHERE.gmm ;\
	echo gmm LD '-L@DIR@/lib'>> $(SRCDIR)/$(INSTALL)/lib/WHERE.gmm ;\
	fi
