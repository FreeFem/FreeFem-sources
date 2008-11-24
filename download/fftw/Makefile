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
<<<<<<< Makefile
INSTALL=../../install
FFTW_VERSION=3.1.2
=======
INSTALL=../..
FFTW_VERSION=3.1.2
>>>>>>> 1.2

fftw: $(SRCDIR)/$(INSTALL)

$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cd $(SRCDIR) && ./configure  --prefix=`pwd`/$(INSTALL) CXX=$(CXX) CC=$(CC) 
	cd $(SRCDIR) && make
	cd $(SRCDIR) && make install

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)

$(PACKAGE):
	curl -O  $(SERVER)/$(PACKAGE)

clean-local:
	-rm -rf fftw-*
