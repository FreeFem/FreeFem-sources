# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
all-local: mshmet

include ff-flags

# Downloading and compiling mshmet
# -------------------------------
# 
DIRPKG= ../pkg
SRCDIR= ./mshmet$(mshmet_VERSION)
#-$(mshmet_VERSION)
PACKAGE=$(DIRPKG)/mshmet$(mshmet_VERSION).tgz
SERVER=http://www.ann.jussieu.fr/~frey/ftp/archives/
INSTALL=../..
mshmet_VERSION=.2009.09.16


# ---------------------- 
#     mshmetlib

MSHMET_DIR = $(abs_top_builddir)/download/mshmet/mshmet$(mshmet_VERSION)
MSHMET_SRCDIRNOLIB = $(MSHMET_DIR)/sources
MSHMET_SRCDIR = $(MSHMET_DIR)/sourceslib
MSHMET_OBJDIR = $(MSHMET_DIR)/objects

mshmet: FAIRE WHERE

FAIRE:$(SRCDIR)/FAIT
	$(MAKE) $(SRCDIR)/$(INSTALL)  
	touch $(SRCDIR)/FAIT

$(SRCDIR)/FAIT:$(SRCDIR)
	cd $(MSHMET_DIR); make
	touch $(SRCDIR)/FAIT

install:$(SRCDIR)/FAIT 
	-mkdir $(SRCDIR)/$(INSTALL)/include/mshmet/
	cp $(MSHMET_SRCDIR)/mshmetlib.h  $(SRCDIR)/$(INSTALL)/include/mshmet/mshmetlib.h
	-mkdir $(SRCDIR)/$(INSTALL)/lib/mshmet/
	cp $(MSHMET_OBJDIR)/libmshmet.a  $(SRCDIR)/$(INSTALL)/lib/mshmet/libmshmet.a

WHERE: 
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	make install;  \
	echo mshmet  LD '-L@DIR@/lib -lmshmet'  >$(SRCDIR)/$(INSTALL)/lib/WHERE.mshmet ;\
	echo mshmet INCLUDE '-I@DIR@/include/mshmet'>> $(SRCDIR)/$(INSTALL)/lib/WHERE.mshmet ;\
	fi


$(SRCDIR): $(PACKAGE)
	-mkdir mshmet$(mshmet_VERSION)
	cd mshmet$(mshmet_VERSION); tar xvzf ../$(PACKAGE)
	-mkdir $(MSHMET_SRCDIR)
	cp $(MSHMET_SRCDIRNOLIB)/*.c $(MSHMET_SRCDIRNOLIB)/*.h $(MSHMET_SRCDIR)
	cp $(MSHMET_SRCDIRNOLIB)/compil.date $(MSHMET_SRCDIR)
	cp mshmetlib-internal.h mshmetlib.c mshmetlib.h $(MSHMET_SRCDIR)/
	cp Makefile-mshmet.inc $(MSHMET_DIR)/makefile
	cd $(MSHMET_SRCDIR); patch -p1 < ../../mshmet.2009.09.16.patch

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`

clean-local: 
	-make -C $(MSHMET_DIR) clean;
	-rm $(MSHMET_OBJDIR)/libmshmet.a

clean: clean-local
	-rm ff-flags
	-rm -r $(MSHMET_DIR)
	-rm $(SRCDIR)/$(INSTALL)/lib/mshmet/*.a
	-rm $(SRCDIR)/$(INSTALL)/include/mshmet/*.h
	-rm -r $(SRCDIR)/$(INSTALL)/lib/mshmet
	-rm -r $(SRCDIR)/$(INSTALL)/include/mshmet
	-rm -rf $(SRCDIR)
	-rm FAIT

ff-flags: ../Makefile 
	grep 'abs_top_builddir *=' ../Makefile >> ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags
	grep 'mshmet_VERSION *=' ./Makefile >> ff-flags


.PHONY: $(SRCDIR)/$(INSTALL)