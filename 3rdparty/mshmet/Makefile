# Downloading and compiling extra libraries
# -----------------------------------------

all-local: mshmet

include ff-flags

# Now in distrib  and compiling mshmet
# -------------------------------
# 
DIRPKG= ../pkg
SRCDIR= ./src
#-$(mshmet_VERSION)
INSTALL=../..
mshmet_VERSION=.FH
FAIRE=$(SRCDIR)/FAIRE

# ---------------------- 
#     mshmetlib

MSHMET_DIR = $(abs_top_builddir)/3rdparty/mshmet/src
MSHMET_SRCDIR = $(MSHMET_DIR)/sourceslib
MSHMET_OBJDIR = $(MSHMET_DIR)/objects

mshmet: $(FAIRE) 

$(FAIRE):$(SRCDIR)/FAIT
	$(MAKE) WHERE $(SRCDIR)/$(INSTALL)  
	touch $(FAIRE)

$(SRCDIR)/FAIT:$(SRCDIR)/tag-tar
	cd $(MSHMET_DIR); $(MAKE)
	touch $(SRCDIR)/FAIT

# FFCS - libMesh is also required by yams, so we move all the rules to [[file:../../../Makefile.am]]
install:$(SRCDIR)/FAIT
	cp $(MSHMET_SRCDIR)/mshmetlib.h  $(SRCDIR)/$(INSTALL)/include/mshmetlib.h
	cp $(MSHMET_OBJDIR)/libmshmet.a  $(SRCDIR)/$(INSTALL)/lib/libmshmet.a

# FFCS - simplify makefile structure for automatic rebuilds
WHERE:install
	echo mshmet  LD -L@DIR@/lib -lmshmet  >$(SRCDIR)/$(INSTALL)/lib/WHERE.mshmet
	echo mshmet INCLUDE -I@DIR@/include>> $(SRCDIR)/$(INSTALL)/lib/WHERE.mshmet


$(SRCDIR)/tag-tar:
	touch $(SRCDIR)/tag-tar


clean:
	-rm ff-flags
#	FFCS - make sure that all directories are cleaned. Thisis especially important under Windows because there is no
#	compilation dependencies control there (see
#	[[file:c:/cygwin/home/alh/ffcs/dist/configure.ac::dependency_tracking]])
	-rm -rf mshmet.????.??.??
	-rm FAIT $(FAIRE) 

ff-flags: ../Makefile 
	grep 'abs_top_builddir *=' ../Makefile >> ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags


.PHONY: $(SRCDIR)/$(INSTALL)
