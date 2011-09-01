# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
all-local: mmg3d

include ff-flags

# Downloading and compiling mmg3d
# -------------------------------
# 
DIRPKG= ../pkg
SRCDIR= ./mmg3d4# lib#-$(mmg3d_VERSION)
PACKAGE=$(DIRPKG)/mmg3d4.0.tgz #mmg3dlib.tar.gz #-$(mmg3d_VERSION).tar.gz
INSTALL=..
mmg3d_VERSION=
# mmg3d pas sur internet
LIBMMG3D=$(INSTALL)/lib/libmmg3d-v4.a
OPT=4

mmg3d: build-4


build-4: 
	-make tag-tar
	if [ -d mmg3d4  ] ; then  make $(LIBMMG3D) install-4; else echo "Sorry no source " ; fi
$(LIBMMG3D):
	cd mmg3d4/build/sources/; $(CC) -c $(COPT) *.c;
	$(AR) $(ARFLAGS) $(LIBMMG3D) mmg3d4/build/sources/*.o
	$(CC) mmg3d4/build/sources/*.o -o ../bin/mmg3d
	touch build-4 FAIT-4
install-4:
	-mkdir ../include/mmg3d-v4/		
	cp mmg3d4/build/sources/*.h ../include/mmg3d-v4/
	touch FAIT-4


mmg3d: mmg3d-$(OPT)

mmg3d-1xxxx: 	
	if [ -f $(PACKAGE) ] ;then \
	   $(MAKE) FAIT-1; \
	else \
	   echo " file $(PACKAGE)  not found => no mmg3d tar file (download by hand) "; \
	fi
mmg3d-2xxx: FAIT-2xxx
FAIT-2xxx: 	
	@if [ -f ../pkg/$(MMG3DLIB_TAR) ] ;then \
	   echo tar zxf ../pkg/$(MMG3DLIB_TAR)  -C .. ;\
	   tar zxf ../pkg/$(MMG3DLIB_TAR)  -C .. ;\
	   touch FAIT-2 ;\
	else \
	   echo " file ../pkg/$(MMG3DLIB_TAR)  not found => no mmg3d tar file (download by hand) "; \
	fi

mmg3d-4:$(PACKAGE)

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET) http://www.math.u-bordeaux1.fr/~dobj/logiciels/download.php?file=`basename $(PACKAGE)`

FAIT-1: tag-tar
	$(MAKE) -C $(SRCDIR) lib
	touch FAIT-1
install:install-4 WHERE
install-1: FAIT-1
	-mkdir $(SRCDIR)/$(INSTALL)/include
	cp $(SRCDIR)/sources/libmmg3d.h  $(SRCDIR)/$(INSTALL)/include/libmmg3d.h
	-mkdir $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/objects/libmmg3d.a  $(LIBMMG3D)
install-2:

WHERE: 
	@-if [ -f FAIT-4 ] ; then \
	make install-4;  \
	echo mmg3d-v4  LD -L@DIR@/lib -lmmg3d-v4  >../lib/WHERE.mmg3d ;\
	echo mmg3d-v4 INCLUDE -I@DIR@/include/mmg3d-v4>> ../lib/WHERE.mmg3d ;\
	echo build WHERE ./lib/WHERE.mmg3d ;\
	fi


FAIRE: FAIT-4 install-4

tag-tar: $(PACKAGE)
	tar xvzf $(PACKAGE)
	touch mmg3d4/build/sources/dataff.h 	
	cd mmg3d4;patch -p1 <../patch-mmg3dv4.diff
	cat </dev/null >mmg3d4/build/sources/mmg3dConfig.h
	touch tag-tar
#	cp makefile-mmg3d.inc $(SRCDIR)/makefile


old-$(PACKAGE): 
	-mkdir $(DIRPKG);
	cd $(DIRPKG); $(WGET)   $(SERVER)/`basename $(PACKAGE)`


clean: 
	-rm ff-flags
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/libmmg3d*.a
	-rm -rf $(SRCDIR)/$(INSTALL)/include/mmg3d-v4
	-rm -r $(SRCDIR)
	-rm FAIT* mmg* flags-* build-4 tag-tar

ff-flags: ../Makefile
	grep 'abs_top_builddir *=' ../Makefile >> ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags
	grep 'MMG3DLIB_TAR *=' ../Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags


.PHONY:$(SRCDIR)/$(INSTALL)