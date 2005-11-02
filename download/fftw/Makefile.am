# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

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
	cd $(SRCDIR) && ./configure  --prefix=`pwd`/../$(INSTALL) CXX=$(CXX) CC=$(CC)
	cd $(SRCDIR) && make
	cd $(SRCDIR) && make install

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)

$(PACKAGE):
	wget -N $(SERVER)/$(PACKAGE)

clean-local:
	-rm -rf fftw-*
