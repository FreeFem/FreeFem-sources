# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_FLTK)

# Downloading and compiling FLTK
# ------------------------------

# FLTK information
SRCDIR=fltk-1.1.5
PACKAGE=fltk-1.1.5-source.tar.bz2
SERVER=http://ovh.dl.sourceforge.net/sourceforge/fltk
INSTALL=install

fltk: $(SRCDIR)/$(INSTALL)

$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cd $(SRCDIR) && ./configure --enable-threads --prefix=`pwd`/$(INSTALL)
	cd $(SRCDIR) && make
	cd $(SRCDIR) && make install

$(SRCDIR): $(PACKAGE)
	tar xvjf $(PACKAGE)

$(PACKAGE):
	wget -N $(SERVER)/$(PACKAGE)

clean-local:
	-rm -r $(PACKAGE) $(SRCDIR)
