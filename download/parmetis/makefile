# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Makefile.in

all-local: parmetis

# Downloading and compiling parmetis
# -------------------------------
# http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/ParMetis-3.1.1.tar.gz
# Metis information
DIRPKG=../pkg
SRCDIR=ParMetis-$(metis_VERSION)
PACKAGE=$(DIRPKG)/ParMetis-$(metis_VERSION).tar.gz
SERVER=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis
INSTALL=../..
metis_VERSION=3.1.1

parmetis: FAIRE

$(SRCDIR)/FAIT:
	$(MAKE) $(SRCDIR)/$(INSTALL)  
	touch $(SRCDIR)/FAIT

FAIRE: $(SRCDIR)/FAIT

Makefile.in: ../../config.status	Makefile
	grep 'abs_top_builddir *=' ../Makefile >> dirflags
	cat dirflags Makefile-parmetis.in > Makefile-parmetis-tmp.in
	../../config.status  --file="Makefile.in:Makefile-parmetis-tmp.in"
	case `uname` in *CYGWIN*)  cp Makefile.in Makefile.inn; sed "s/COPTIONS =/COPTIONS = -D__VC__/" <Makefile.inn > Makefile.in; rm Makefile.inn ;; esac 
	rm dirflags Makefile-parmetis-tmp.in

$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cp Makefile.in $(SRCDIR)/
	cd $(SRCDIR)/METISLib; make
	cd $(SRCDIR)/ParMETISLib; make 
	mkdir -p $(SRCDIR)/$(INSTALL)/include/parmetis
	cp $(SRCDIR)/METISLib/*.h $(SRCDIR)/$(INSTALL)/include/parmetis
	cp $(SRCDIR)/parmetis.h $(SRCDIR)/$(INSTALL)/include/parmetis
	mv $(SRCDIR)/$(INSTALL)/include/parmetis/metis.h $(SRCDIR)/$(INSTALL)/include/parmetis/metis-tmp.h
	sed -e 's;#include "../parmetis.h";#include "parmetis.h";'\
		-e 's;#include <defs.h>;#include "defs.h";'\
	<$(SRCDIR)/$(INSTALL)/include/parmetis/metis-tmp.h \
	>$(SRCDIR)/$(INSTALL)/include/parmetis/metis.h
	mkdir -p $(SRCDIR)/$(INSTALL)/lib/parmetis
	cp $(SRCDIR)/lib*metis.a $(SRCDIR)/$(INSTALL)/lib/parmetis/

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)
	cd $(SRCDIR)/ParMETISLib/; patch -p1 < ../../parmetis-3.1.1.patch
	touch $(SRCDIR)

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`


clean-local:
	-cd $(SRCDIR); make realclean
clean: clean-local
	-rm Makefile.in
	-rm $(SRCDIR)/$(INSTALL)/lib/parmetis/*.a
	-rm $(SRCDIR)/$(INSTALL)/include/parmetis/*.h
	-rm -rf $(SRCDIR)
.PHONY:$(SRCDIR)/$(INSTALL)