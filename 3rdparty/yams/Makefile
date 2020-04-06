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

YAMS_DIR = $(abs_top_builddir)/3rdparty/yams/$(SRCDIR)
YAMS_SRCDIRNOLIB = $(YAMS_DIR)/sources
YAMS_SRCDIR = $(YAMS_DIR)/sourceslib
YAMS_OBJDIR = $(YAMS_DIR)/objects

yams: FAIRE  



# FFCS - make sure that PATCH is done sequentially otherwise its error messages are drowned into other meaningless
# parallel compilation messages

FAIRE: $(SRCDIR)/PATCH yamslib_internal.h yamslib.c yamslib.h  ../Makefile
	$(MAKE)  install  WHERE
	touch FAIRE

$(SRCDIR)/FAIT: $(SRCDIR)/PATCH  yamslib_internal.h yamslib.c yamslib.h 
	cp yamslib_internal.h yamslib.c yamslib.h $(YAMS_SRCDIR)/
	cp makefile-yams.inc $(YAMS_DIR)/makefile
	cd $(YAMS_DIR); make
	touch $(SRCDIR)/FAIT
install: $(SRCDIR)/FAIT
	sed  s/defines.h/freeyams_defines.h/ <$(YAMS_SRCDIR)/yamslib.h  >../include/freeyamslib.h
	cp $(YAMS_SRCDIR)/defines.h  ../include/freeyams_defines.h
	-mkdir ../lib	
	cp $(YAMS_OBJDIR)/libyams.a  ../lib/libfreeyams.a

# FFCS - WHERE is made to depend on FAIT otherwise it may be built in parallel and not be activated because FAIT is not
# there yet

WHERE: $(SRCDIR)/FAIT
	echo freeyams  LD -L@DIR@/lib -lfreeyams  >$(SRCDIR)/$(INSTALL)/lib/WHERE.freeyams ;
	echo freeyams INCLUDE -I@DIR@/include>> $(SRCDIR)/$(INSTALL)/lib/WHERE.freeyams ;

$(SRCDIR)/PATCH: $(PACKAGE) 
	-mkdir -p $(SRCDIR)
	cd $(SRCDIR); tar xzf ../$(PACKAGE)
	-mkdir $(YAMS_SRCDIR)
	cp $(YAMS_SRCDIRNOLIB)/*.c $(YAMS_SRCDIRNOLIB)/*.h $(YAMS_SRCDIR)
	rm $(YAMS_SRCDIR)/memory.c 
	cp $(YAMS_SRCDIRNOLIB)/compil.date $(YAMS_SRCDIR)
	cd $(YAMS_SRCDIR) && \
	  patch -p2 < ../../freeyams$(yams_VERSION).patch && \
	  patch -p2 < ../../freeyams$(yams_VERSION)-return-values.patch
	mv  $(YAMS_SRCDIR)/yams.c $(YAMS_SRCDIR)/..
	touch $(SRCDIR)/PATCH

$(PACKAGE):
	../getall -o freeYams -a

# FFCS: only run make clean if cd to SRCDIR worked, otherwise this is infinite loop.

clean-local: 
	-rm FAIRE FAIT  $(SRCDIR)/FAIT
	-cd $(YAMS_DIR) &&  $(MAKE) -C $(YAMS_DIR)  clean
	-rm $(YAMS_OBJDIR)/libyams.a

# FFCS -simplifying all paths

 clean: clean-local
	-rm ff-flags
	-rm $(SRCDIR)/$(INSTALL)/lib/libfreeyams.a
	-rm $(SRCDIR)/$(INSTALL)/include/*freeyams*.h
	-rm $(SRCDIR)/$(INSTALL)/lib/WHERE.freeyams
	-rm ../lib/libfreeyams.a
	-rm ../include/*freeyams*.h
	-rm ../lib/WHERE.freeyams
	-rm -rf $(YAMS_DIR)
	-rm -rf $(SRCDIR)
#FH 	-rm $(PACKAGE)
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
