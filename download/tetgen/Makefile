# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

include cxxflags
all-local: tetgen

# Downloading and compiling Tetgen
# ------------------------------

# Tetgen information
DIRPKG=../pkg
SRCDIR=tetgen$(tetgen_VERSION)
PACKAGE=$(DIRPKG)/tetgen$(tetgen_VERSION).tar.gz
SERVER=http://www.wias-berlin.de/people/si
INSTALL=../..
tetgen_VERSION=1.4.2

tetgen: FAIT

FAIT:
	$(MAKE) FAIRE
	touch FAIT

FAIRE:$(SRCDIR)/$(INSTALL) cxxflags
	touch FAIT

$(SRCDIR)/$(INSTALL): $(SRCDIR) cxxflags
	cd $(SRCDIR);$(CXX) $(CXXFLAGS) -DSELF_CHECK  -DNDEBUG -DTETLIBRARY -c tetgen.cxx
	cd $(SRCDIR);$(CXX) $(CXXFLAGS) -DSELF_CHECK  -DNDEBUG -DTETLIBRARY -c predicates.cxx
	cd $(SRCDIR);$(AR) $(ARFLAGS)  $(INSTALL)/lib/libtet.a tetgen.o predicates.o
	cp $(SRCDIR)/tetgen.h $(SRCDIR)/$(INSTALL)/include 


$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)
	cd tetgen1.4.2;patch -p1 <../tetgen1.4.2.patch 

$(PACKAGE):
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`



clean-local:
	-rm -rf tetgen* FAIT
clean:
	-rm $(SRCDIR)/$(INSTALL)/lib/libtet.a
	-rm $(SRCDIR)/$(INSTALL)/include/tetgen.h
	-rm FAIT
	-rm -rf tetgem*
cxxflags: ../Makefile  Makefile
	grep 'CXX *=' ../Makefile >cxxflags
	grep 'CC *=' ../Makefile >>cxxflags
#  bug under ubuntu et cygwin si option optime  -O3 remove optim ???
	grep 'CXXFLAGS *=' ../Makefile | sed 's/ -O[0-9]* / /'  >>cxxflags
	grep 'WGET *=' ../Makefile >>cxxflags
	grep 'AR *=' ../Makefile >>cxxflags
	grep 'ARFLAGS *=' ../Makefile >>cxxflags
.PHONY:$(SRCDIR)/$(INSTALL)