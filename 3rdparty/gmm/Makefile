# Downloading and compiling extra libraries
# -----------------------------------------

include cxxflags

# Downloading and compiling FFTW
# ------------------------------
PKG=gmm
# $(PKG) information
SRCDIR=$(PKG)-$(PKG_VERSION)
PACKAGE=$(PKG)-$(PKG_VERSION).tar.gz
SERVER=http://download.gna.org/getfem/stable
PKGDIR=../pkg
PKG_VERSION=4.2
INSTALL=../..


$(PKG):FAIT 

FAITwin32-dll-target:
	echo "On Pure Win32 (to hard to compile) "

FAIT: 
	$(MAKE) install  
	touch FAIT
	$(MAKE) WHERE  	

FAIRE: $(SRCDIR) 
	cd $(SRCDIR) && ./configure --disable-dependency-tracking   --prefix=`pwd`/$(INSTALL) CXX="$(CXX)" CC="$(CC)"  CFLAGS="$(CFLAGS)" CPP='gcc -E' CXXFLAGS="$(CXXFLAGS)"
	cd $(SRCDIR) && make
	touch FAIRE
install: FAIRE 
	cd $(SRCDIR) && make install


$(SRCDIR): $(PKGDIR)/$(PACKAGE)
	gunzip -c   $^ | tar xf - 

$(PKGDIR)/$(PACKAGE):
	../getall -o Gmm++ -a

clean: clean-local

clean-local:
	-rm -rf $(PKG)-* FAIT FAIRE $(SRCDIR) 
	-rm ../lib/WHERE.gmm

cxxflags: ../Makefile
	grep 'CXX *=' ../Makefile >cxxflags
	grep 'CC *=' ../Makefile >>cxxflags
	grep 'CXXFLAGS *=' ../Makefile >>cxxflags
	grep 'CFLAGS *=' ../Makefile >>cxxflags
	grep 'WGET *=' ../Makefile >>cxxflags
	grep 'WIN32DLLTARGET *=' ../Makefile >>cxxflags

WHERE: 
	-@if [ -f FAIT ] ; then \
	echo build ../lib/WHERE.gmm ;\
	echo gmm INCLUDE -I@DIR@/include  > ../lib/WHERE.gmm ;\
	echo gmm LD -L@DIR@/lib   >> ../lib/WHERE.gmm ;\
	fi
