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
	$(MAKE) -C $(SRCDIR) lib
	touch $(SRCDIR)/FAIT
install: $(SRCDIR)/FAIT
	-mkdir $(SRCDIR)/$(INSTALL)/include
	cp $(SRCDIR)/sources/libmmg3d.h  $(SRCDIR)/$(INSTALL)/include/libmmg3d.h
	-mkdir $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/objects/libmmg3d.a  $(LIBMMG3D)

WHERE: 
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	make install;  \
	echo mmg3d  LD '-L@DIR@/lib -lmmg3d'  >$(SRCDIR)/$(INSTALL)/lib/WHERE.mmg3d ;\
	echo mmg3d INCLUDE '-I@DIR@/include'>> $(SRCDIR)/$(INSTALL)/lib/WHERE.mmg3d ;\
	fi


FAIRE: $(SRCDIR)/FAIT install

$(SRCDIR)/tag-tar: $(PACKAGE)
	tar xvzf $(PACKAGE)
	touch $(SRCDIR)/tag-tar
	cp makefile-mmg3d.inc $(SRCDIR)/makefile


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