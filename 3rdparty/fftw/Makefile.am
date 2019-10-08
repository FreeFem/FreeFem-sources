# Downloading and compiling extra libraries
# -----------------------------------------

all-local: $(DOWNLOAD_FFTW)

# Downloading and compiling FFTW
# ------------------------------

# FFTW information
FFTW_VERSION=3.3.8
SRCDIR=fftw-$(FFTW_VERSION)
PACKAGE=fftw-$(FFTW_VERSION).tar.gz
SERVER=http://www.fftw.org
THIS=fftw3

fftw: $(SRCDIR)/FAIT WHERE
	$(MAKE) WHERE

$(SRCDIR)/FAIT:$(SRCDIR)/FAIT-1
	cd $(SRCDIR) && ./configure --enable-shared --disable-dependency-tracking  --disable-fortran  --prefix=`cd ../..; pwd`  CXX="$(CXX)" CC="$(CC)"  CFLAGS="$(CFLAGS)" CPP='gcc -E' CXXFLAGS="$(CXXFLAGS)"
	cd $(SRCDIR) && make
	cd $(SRCDIR) && make install
	touch $(SRCDIR)/FAIT

$(SRCDIR)/FAIT-1:../pkg/$(PACKAGE)
	tar xvzf ../pkg/$(PACKAGE)
	touch $(SRCDIR)/FAIT-1

../pkg/$(PACKAGE):
	-mkdir ../pkg
	cd ../pkg;@WGET@ -N $(SERVER)/$(PACKAGE)

WHERE:$(SRCDIR)/FAIT
	-if [ $(SRCDIR)/FAIT ] ; then \
	echo $(THIS)  LD -L@DIR@/lib -l$(THIS)  >../lib/WHERE.$(THIS) ;\
	echo $(THIS) INCLUDE -I@DIR@/include   >> ../lib/WHERE.$(THIS) ;\
	fi


clean-local:
	-rm -rf fftw-* 
	-rm ../include/fftw3.f ../include/fftw3.f03 ../include/fftw3.h ../include/fftw3l.f03 ../include/fftw3q.f03	\
		../lib/libfftw3.a ../lib/libfftw3.la
