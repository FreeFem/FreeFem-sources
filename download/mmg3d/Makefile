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
ifeq "$(MMG3DLIB_TAR)" ""
OPT=1 
else
OPT=2
endif

mmg3d: mmg3d-$(OPT)

mmg3d-1: 	
	if [ -f $(PACKAGE) ] ;then \
	   $(MAKE) FAIRE-1; \
	else \
	   echo " file $(PACKAGE)  not found => no mmg3d tar file (download by hand) "; \
	fi
mmg3d-2: FAIT-2
FAIT-2: 	
	@if [ -f ../pkg/$(MMG3DLIB_TAR) ] ;then \
	   echo tar zxf ../pkg/$(MMG3DLIB_TAR)  -C .. ;\
	   tar zxf ../pkg/$(MMG3DLIB_TAR)  -C .. ;\
	   touch FAIT-2 ;\
	else \
	   echo " file ../pkg/$(MMG3DLIB_TAR)  not found => no mmg3d tar file (download by hand) "; \
	fi

FAIT-1: tag-tar
	$(MAKE) -C $(SRCDIR) lib
	touch FAIT-1
install:install-$(OPT) WHERE
install-1: FAIT-1
	-mkdir $(SRCDIR)/$(INSTALL)/include
	cp $(SRCDIR)/sources/libmmg3d.h  $(SRCDIR)/$(INSTALL)/include/libmmg3d.h
	-mkdir $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/objects/libmmg3d.a  $(LIBMMG3D)
install-2:

WHERE: 
	@-if [ -f FAIT-$(OPT) ] ; then \
	make install-$(OPT);  \
	echo mmg3d  LD '-L@DIR@/lib -lmmg3d'  >../lib/WHERE.mmg3d ;\
	echo mmg3d INCLUDE '-I@DIR@/include'>> ../lib/WHERE.mmg3d ;\
	echo build WHERE ./lib/WHERE.mmg3d ;\
	fi


FAIRE: FAIT install

tag-tar: $(PACKAGE)
	tar xvzf $(PACKAGE)
	touch tag-tar
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
	grep 'MMG3DLIB_TAR *=' ../Makefile >> ff-flags


.PHONY:$(SRCDIR)/$(INSTALL)