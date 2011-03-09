# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Makefile.in

all-local: metis

# Downloading and compiling Tetgen
# ------------------------------
# http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-4.0.tar.gz
# Metis information
DIRPKG=../pkg
SRCDIR=metis-$(metis_VERSION)
PACKAGE=$(DIRPKG)/metis-$(metis_VERSION).tar.gz
SERVER=http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis
INSTALL=../..
metis_VERSION=4.0

metis: FAIRE


$(SRCDIR)/FAIT:
	$(MAKE) install WHERE
	touch $(SRCDIR)/FAIT

install:$(SRCDIR)/tag-compile
#	$(MAKE) install
	mkdir -p $(SRCDIR)/$(INSTALL)/include/metis
	cp $(SRCDIR)/Lib/*.h $(SRCDIR)/$(INSTALL)/include/metis


FAIRE: $(SRCDIR)/FAIT 
	touch fAIRE

Makefile.in: ../../config.status Makefile-metis.in
	../../config.status  --file="Makefile.in:Makefile-metis.in"

$(SRCDIR)/tag-compile: $(SRCDIR)/tags
	cp Makefile.in $(SRCDIR)
	cd $(SRCDIR)/Lib;make
	-cd $(SRCDIR)/Programs;make 
	touch $(SRCDIR)/tag-compile



WHERE: 
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	echo metis LD '-L@DIR@/lib -lmetis'  >$(SRCDIR)/$(INSTALL)/lib/WHERE.metis ;\
	echo metis INCLUDE '-I@DIR@/include/metis'>> $(SRCDIR)/$(INSTALL)/lib/WHERE.metis ;\
	fi




$(SRCDIR)/tags: $(PACKAGE)
	tar xvzf $(PACKAGE)
	patch -p0  <patch-metis
	mv  $(SRCDIR)/Programs/Makefile $(SRCDIR)/Programs/Makefile-orig
	sed -e 's;BINDIR = ..;BINDIR = ../$(INSTALL)/bin;' \
            -e 's;../libmetis.a;../$(INSTALL)/lib/libmetis.a;' \
            -e 's;-L[.][.];-L../$(INSTALL)/lib;' \
          <$(SRCDIR)/Programs/Makefile-orig \
          >$(SRCDIR)/Programs/Makefile
	mv $(SRCDIR)/Lib/Makefile $(SRCDIR)/Lib/Makefile-orig
	sed   -e 's;../libmetis.a;../$(INSTALL)/lib/libmetis.a;' \
          <$(SRCDIR)/Lib/Makefile-orig \
          >$(SRCDIR)/Lib/Makefile
	touch $(SRCDIR)/tags

$(PACKAGE):
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`



clean-local:
	-cd $(SRCDIR) && $(MAKE) realclean  -C $(SRCDIR)
	-rm -rf metis*
clean: clean-local
	-rm $(SRCDIR)/$(INSTALL)/lib/libmetis.a
	-rm $(SRCDIR)/$(INSTALL)/include/metis.h
	-rm -rf $(SRCDIR)
	-rm FAIRE FAIT 