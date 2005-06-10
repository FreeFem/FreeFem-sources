# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$


all-local:fftw

# Downloading and compiling FFTW
# ------------------------------

# FFTW information

SRCDIR=fftw-$(FFTW_VERSION)
PACKAGE=fftw-$(FFTW_VERSION).tar.gz
SERVER=http://www.fftw.org
INSTALL=../../install
FFTW_VERSION=3.0.1

fftw: $(SRCDIR)/$(INSTALL)

$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cd $(SRCDIR) && ./configure  --prefix=`pwd`/$(INSTALL) CXX=$(CXX) CC=$(CC) 
	cd $(SRCDIR) && make
	cd $(SRCDIR) && make install

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)

$(PACKAGE):
	wget  $(SERVER)/$(PACKAGE)

clean-local:
	-rm -rf fftw-*
