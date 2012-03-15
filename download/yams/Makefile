# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
all-local: yams

include ff-flags

# Downloading and compiling yams
# -------------------------------
# 
DIRPKG= ../pkg
SRCDIR= ./freeyams$(yams_VERSION)
#-$(yams_VERSION)
PACKAGE=$(DIRPKG)/freeyams$(yams_VERSION).tgz
SERVER=http://www.ann.jussieu.fr/~frey/ftp/archives/
INSTALL=../..
yams_VERSION=.2012.02.05

# ---------------------- 
#     yamslib

YAMS_DIR = $(abs_top_builddir)/download/yams/$(SRCDIR)
YAMS_SRCDIRNOLIB = $(YAMS_DIR)/sources
YAMS_SRCDIR = $(YAMS_DIR)/sourceslib
YAMS_OBJDIR = $(YAMS_DIR)/objects

yams: FAIRE  



FAIRE: yamslib_internal.h yamslib.c yamslib.h  ../Makefile
	$(MAKE)  install  WHERE
	touch FAIRE

$(SRCDIR)/FAIT: $(SRCDIR)/PATCH  yamslib_internal.h yamslib.c yamslib.h 
	cp yamslib_internal.h yamslib.c yamslib.h $(YAMS_SRCDIR)/
	cp makefile-yams.inc $(YAMS_DIR)/makefile
	cd $(YAMS_DIR); make
	touch $(SRCDIR)/FAIT
install: $(SRCDIR)/FAIT
	sed  s/defines.h/freeyams_defines.h/ <$(YAMS_SRCDIR)/yamslib.h  >$(SRCDIR)/$(INSTALL)/include/freeyamslib.h
	cp $(YAMS_SRCDIR)/defines.h  $(SRCDIR)/$(INSTALL)/include/freeyams_defines.h
	-mkdir $(SRCDIR)/$(INSTALL)/lib	
	cp $(YAMS_OBJDIR)/libyams.a  $(SRCDIR)/$(INSTALL)/lib/libfreeyams.a
WHERE: 
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	make install;  \
	echo freeyams  LD -L@DIR@/lib -lfreeyams  >$(SRCDIR)/$(INSTALL)/lib/WHERE.freeyams ;\
	echo freeyams INCLUDE -I@DIR@/include>> $(SRCDIR)/$(INSTALL)/lib/WHERE.freeyams ;\
	fi


$(SRCDIR)/PATCH: $(PACKAGE) 
	-mkdir -p $(SRCDIR)
	cd $(SRCDIR); tar xvzf ../$(PACKAGE)
	-mkdir $(YAMS_SRCDIR)
	cp $(YAMS_SRCDIRNOLIB)/*.c $(YAMS_SRCDIRNOLIB)/*.h $(YAMS_SRCDIR)
	rm $(YAMS_SRCDIR)/memory.c 
	cp $(YAMS_SRCDIRNOLIB)/compil.date $(YAMS_SRCDIR)
	cd $(YAMS_SRCDIR); patch -p2 < ../../freeyams$(yams_VERSION).patch
	mv  $(YAMS_SRCDIR)/yams.c $(YAMS_SRCDIR)/..
	touch $(SRCDIR)/PATCH

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`

clean-local: 
	-rm FAIRE FAIT  $(SRCDIR)/FAIT
	-cd $(YAMS_DIR) &&  $(MAKE) -C $(YAMS_DIR)  clean
	-rm $(YAMS_OBJDIR)/libyams.a

clean: clean-local
	-rm ff-flags
	-rm $(SRCDIR)/$(INSTALL)/lib/libfreeyams.a
	-rm $(SRCDIR)/$(INSTALL)/include/*freeyams*.h
	-rm $(SRCDIR)/$(INSTALL)/lib/WHERE.freeyams
	-rm -rf $(YAMS_DIR)
	-rm -rf $(SRCDIR)
	-rm FAIT FAIRE

ff-flags: ../Makefile Makefile
	grep 'abs_top_builddir *=' ../Makefile > ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags
	grep 'yams_VERSION *=' ./Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags
	grep 'LIBS *=' ../Makefile >> ff-flags


.PHONY: $(SRCDIR)/$(INSTALL)