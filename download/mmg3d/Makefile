# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
all-local: mmg3d

include ff-flags

# Downloading and compiling mmg3d
# -------------------------------
# 
DIRPKG= ../pkg
SRCDIR= ./mmg3dlib#-$(mmg3d_VERSION)
PACKAGE=$(DIRPKG)/mmg3dlib.tar.gz #-$(mmg3d_VERSION).tar.gz
INSTALL=../..
mmg3d_VERSION=
# mmg3d pas sur internet
SERVER=
WGET =

mmg3d: FAIRE

$(SRCDIR)/FAIT:
	$(MAKE) $(SRCDIR)/$(INSTALL)  
	touch FAIT

FAIRE: $(SRCDIR)/FAIT

$(SRCDIR)/$(INSTALL): $(PACKAGE) $(SRCDIR)
	cp makefile-mmg3d.inc $(SRCDIR)/makefile
	cd $(SRCDIR); $(MAKE) lib
	-mkdir $(SRCDIR)/$(INSTALL)/include/mmg3d/
	cp $(SRCDIR)/sources/libmmg3d.h  $(SRCDIR)/$(INSTALL)/include/mmg3d/libmmg3d.h
	-mkdir $(SRCDIR)/$(INSTALL)/lib/mmg3d/
	cp $(SRCDIR)/objects/libmmg3d.a  $(SRCDIR)/$(INSTALL)/lib/mmg3d/libmmg3d.a

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`


clean: 
	-rm ff-flags
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/mmg3d/
	-rm -rf $(SRCDIR)/$(INSTALL)/include/mmg3d/
	-rm -r $(SRCDIR)
	-rm FAIT

ff-flags: ../Makefile
	grep 'abs_top_builddir *=' ../Makefile >> ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags


.PHONY:$(SRCDIR)/$(INSTALL)