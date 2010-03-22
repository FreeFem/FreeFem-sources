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
SERVER=
INSTALL=../..
mmg3d_VERSION=
WGET =

# ---------------------- 
#     mmg3dlib

# working dirs
MMG3D_EXEDIR =     #$(HOME)/bin/
MMG3D_DIR    = $(abs_top_builddir)/download/mmg3d/mmg3dlib
MMG3D_SRCDIR = $(MMG3D_DIR)/sources
MMG3D_OBJDIR = $(MMG3D_DIR)/objects
MMG3D_ARCDIR = $(MMG3D_DIR)/archives
MMG3D_DIRDIR = $(MMG3D_DIR)/objects $(MMG3D_OBJDIR) $(MMG3D_ARCDIR)
MMG3D_INCDIR = -I$(MMG3D_DIR)/sources/ -I$(abs_top_builddir)/src/libMesh/
MMG3D_LDLDIR = 


# objects list
MMG3D_src    = $(wildcard $(MMG3D_SRCDIR)/*.c)
MMG3D_header = $(wildcard $(MMG3D_SRCDIR)/*.h)
MMG3D_objs   = $(patsubst $(MMG3D_SRCDIR)%,$(MMG3D_OBJDIR)%,$(MMG3D_src:.c=.o))
MMG3D_lib    = $(MMG3D_OBJDIR)/libmmg3d.a

## end makefile mmg3d

mmg3d: FAIRE

$(SRCDIR)/FAIT:
	$(MAKE) $(SRCDIR)/$(INSTALL)  
	touch $(SRCDIR)/FAIT

FAIRE: $(SRCDIR)/FAIT

$(SRCDIR)/$(INSTALL): $(SRCDIR) $(MMG3D_objs)
	$(AR) $(ARFLAGS) $(MMG3D_lib) $(MMG3D_OBJDIR)/*.o  
	$(RANLIB) $(MMG3D_lib)
	-mkdir $(SRCDIR)/$(INSTALL)/include/mmg3d/
	cp $(MMG3D_SRCDIR)/libmmg3d.h  $(SRCDIR)/$(INSTALL)/include/mmg3d/libmmg3d.h
	-mkdir $(SRCDIR)/$(INSTALL)/lib/mmg3d/
	cp $(MMG3D_OBJDIR)/libmmg3d.a  $(SRCDIR)/$(INSTALL)/lib/mmg3d/libmmg3d.a

## makefile mmg3d

$(MMG3D_OBJDIR)/%.o: $(MMG3D_SRCDIR)/%.c
	$(CC) $(OPT64) $(MMG3D_INCDIR) $(CFLAGS) -c $< -o $@

$(MMG3D_objs): $(PACKAGE) $(MMG3D_header)

## end makefile mmg3d

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`



clean: 
	-rm ff-flags
	-rm $(MMG3D_objs) libmmg3d.a
	-rm $(SRCDIR)/$(INSTALL)/lib/mmg3d/*.a
	-rm $(SRCDIR)/$(INSTALL)/include/mmg3d/*.h
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


.PHONY:$(SRCDIR)/$(INSTALL)