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
WGET=
LIBMMG3D=$(SRCDIR)/$(INSTALL)/lib/libmmg3d.a


mmg3d: 	
	-if [ -f $(PACKAGE) ] ;then \
	   $(MAKE) FAIRE; \
	else \
	   echo " file $(PACKAGE)  not found => no mmg3d tar file (download by hand) "; \
	fi

$(SRCDIR)/FAIT: $(SRCDIR)/tag-tar
	cp makefile-mmg3d.inc $(SRCDIR)/makefile
	$(MAKE) -C $(SRCDIR) lib
	-mkdir $(SRCDIR)/$(INSTALL)/include
	cp $(SRCDIR)/sources/libmmg3d.h  $(SRCDIR)/$(INSTALL)/include/libmmg3d.h
	-mkdir $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/objects/libmmg3d.a  $(LIBMMG3D)
	touch FAIT

FAIRE: $(SRCDIR)/FAIT

$(SRCDIR)/tag-tar: $(PACKAGE)
	tar xvzf $(PACKAGE)
	touch $(SRCDIR)/tag-tar

$(PACKAGE): 
	-mkdir $(DIRPKG);
	cd $(DIRPKG); $(WGET)   $(SERVER)/`basename $(PACKAGE)`


clean: 
	-rm ff-flags
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/libmmg3d.a
	-rm -rf $(SRCDIR)/$(INSTALL)/include/libmmg3d.a
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