# Downloading and compiling extra libraries
# -----------------------------------------

all-local: $(DOWNLOAD_FFTW)

# Downloading and compiling FLTK
# ------------------------------

# FFTW information
SRCDIR=fltk-$(FFTW_VERSION)
PACKAGE=fltk-$(FFTW_VERSION).tar.gz
SERVER=http://www.fftw.org
INSTALL=install

fltk: $(SRCDIR)/$(INSTALL)

$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cd $(SRCDIR) && ./configure  --prefix=`pwd`/../.. CXX=$(CXX) CC=$(CC)
	cd $(SRCDIR) && make
	cd $(SRCDIR) && make install

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)

$(PACKAGE):
	@WGET@ -N $(SERVER)/$(PACKAGE)

clean-local:
	-rm -rf fftw-* 
	-rm ../include/fftw3.f ../include/fftw3.f03 ../include/fftw3.h ../include/fftw3l.f03 ../include/fftw3q.f03	\
		../lib/libfftw3.a ../lib/libfftw3.la
