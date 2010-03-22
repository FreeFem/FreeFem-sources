# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
all-local: yams

include ff-flags

# Downloading and compiling yams
# -------------------------------
# 
DIRPKG= ../pkg
SRCDIR= ./yams2$(yams_VERSION)
#-$(yams_VERSION)
PACKAGE=$(DIRPKG)/yams2$(yams_VERSION).tgz
SERVER=
INSTALL=../..
yams_VERSION=.2010.02.22
WGET =

# ---------------------- 
#     yamslib

YAMS_DIR = $(abs_top_builddir)/download/yams/yams2$(yams_VERSION)
YAMS_SRCDIRNOLIB = $(YAMS_DIR)/sources
YAMS_SRCDIR = $(YAMS_DIR)/sourceslib
YAMS_OBJDIR = $(YAMS_DIR)/objects

yams: FAIRE

$(SRCDIR)/FAIT:
	$(MAKE) $(SRCDIR)/$(INSTALL)  
	touch $(SRCDIR)/FAIT

FAIRE: $(SRCDIR)/FAIT

$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cd $(YAMS_DIR); make
	-mkdir $(SRCDIR)/$(INSTALL)/include/yams/
	cp $(YAMS_SRCDIR)/yamslib.h  $(SRCDIR)/$(INSTALL)/include/yams/yamslib.h
	cp $(YAMS_SRCDIR)/defines.h  $(SRCDIR)/$(INSTALL)/include/yams/
	-mkdir $(SRCDIR)/$(INSTALL)/lib/yams/
	cp $(YAMS_OBJDIR)/libyams2.a  $(SRCDIR)/$(INSTALL)/lib/yams/libyams2.a

$(SRCDIR): $(PACKAGE)
	-mkdir yams2$(yams_VERSION)
	cd yams2$(yams_VERSION); tar xvzf ../$(PACKAGE)
	-mkdir $(YAMS_SRCDIR)
	cp $(YAMS_SRCDIRNOLIB)/*.c $(YAMS_SRCDIRNOLIB)/*.h $(YAMS_SRCDIR)
	cp $(YAMS_SRCDIRNOLIB)/compil.date $(YAMS_SRCDIR)
	cp yamslib_internal.h yamslib.c yamslib.h $(YAMS_SRCDIR)/
	cp makefile-yams.inc $(YAMS_DIR)/makefile
	cd $(YAMS_SRCDIR); patch -p1 < ../../yams2$(yams_VERSION).patch

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`

clean-local: 
	cd $(YAMS_DIR); make clean;
	-rm $(YAMS_OBJDIR)/libyams2.a

clean: clean-local
	-rm ff-flags
	-rm -r $(YAMS_DIR)
	-rm $(SRCDIR)/$(INSTALL)/lib/yams/*.a
	-rm $(SRCDIR)/$(INSTALL)/include/yams/*.h
	-rm -r $(SRCDIR)/$(INSTALL)/lib/yams
	-rm -r $(SRCDIR)/$(INSTALL)/include/yams
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
	grep 'yams_VERSION *=' ./Makefile >> ff-flags


.PHONY: $(SRCDIR)/$(INSTALL)