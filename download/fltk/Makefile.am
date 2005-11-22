# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_FLTK)

# Downloading and compiling FLTK
# ------------------------------

# FLTK information
SRCDIR=fltk-$(FLTK_VERSION)
PACKAGE=fltk-$(FLTK_VERSION)-source.tar.bz2
SERVER=http://ovh.dl.sourceforge.net/sourceforge/fltk
INSTALL=install

fltk: $(SRCDIR)/$(INSTALL)

# Using system JPEG and PNG lib is a problem under Cygwin
$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cd $(SRCDIR) && ./configure --enable-threads --enable-localjpeg	\
		--enable-localpng $(FLTK_CONFIG_PARAM) --disable-gl	\
		@FLTK_XONLY@ --prefix=`pwd`/$(INSTALL) CXX=$(CXX) CC=$(CC)
	cd $(SRCDIR) && make
	cd $(SRCDIR) && make install

$(SRCDIR): $(PACKAGE)
	tar xvjf $(PACKAGE)

$(PACKAGE):
	@WGET@ -N $(SERVER)/$(PACKAGE)

clean-local:
	-rm -r fltk-*
