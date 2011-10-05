# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Makefile.in

all-local: metis

# Downloading and compiling Tetgen
# ------------------------------
# http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.0.1.tar.gz 
# http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-4.0.tar.gz
# Metis information
DIRPKG=../pkg
SRCDIR=metis-$(metis_VERSION)
PACKAGE=$(DIRPKG)/metis-$(metis_VERSION).tar.gz
SERVER=http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis
INSTALL=../..
metis_VERSION=5.0.1

metis: FAIRE_$(metis_VERSION)


$(SRCDIR)/FAIT:
	$(MAKE) install WHERE
	touch $(SRCDIR)/FAIT

install:$(SRCDIR)/tag-compile
#	$(MAKE) install
	mkdir -p $(SRCDIR)/$(INSTALL)/include 
	mkdir -p $(SRCDIR)/$(INSTALL)/lib
#	cp $(SRCDIR)/Lib/*.h $(SRCDIR)/$(INSTALL)/include/metis
	cd $(SRCDIR); make install 


FAIRE_$(metis_VERSION): $(SRCDIR)/FAIT 
	touch FAIRE_$(metis_VERSION)

Makefile.in: ../../config.status Makefile-metis.in
	../../config.status  --file="Makefile.in:Makefile-metis.in"

$(SRCDIR)/tag-compile: $(SRCDIR)/tags
#	cp Makefile.in $(SRCDIR)
#	cd $(SRCDIR)/Lib;make
#	-cd $(SRCDIR)/Programs;make 
	cd $(SRCDIR); make 
	cd $(SRCDIR); make install
	touch $(SRCDIR)/tag-compile



WHERE: 
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	echo metis LD -L@DIR@/lib -lmetis  >$(SRCDIR)/$(INSTALL)/lib/WHERE.metis ;\
	echo metis INCLUDE -I@DIR@/include >> $(SRCDIR)/$(INSTALL)/lib/WHERE.metis ;\
	fi




$(SRCDIR)/tags: $(PACKAGE)
	tar xvzf $(PACKAGE)
	cd $(SRCDIR); make config prefix=`pwd`/../../ cc='$(CC)'
#	patch -p0  <patch-metis	
#	mv  $(SRCDIR)/Programs/Makefile $(SRCDIR)/Programs/Makefile-orig
#	sed -e 's;BINDIR = ..;BINDIR = ../$(INSTALL)/bin;' \
#            -e 's;../libmetis.a;../$(INSTALL)/lib/libmetis.a;' \
#            -e 's;-L[.][.];-L../$(INSTALL)/lib;' \
#          <$(SRCDIR)/Programs/Makefile-orig \
#          >$(SRCDIR)/Programs/Makefile
#	mv $(SRCDIR)/Lib/Makefile $(SRCDIR)/Lib/Makefile-orig
#	sed   -e 's;../libmetis.a;../$(INSTALL)/lib/libmetis.a;' \
#          <$(SRCDIR)/Lib/Makefile-orig \
#          >$(SRCDIR)/Lib/Makefile
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
	-rm FAIRE* FAIT 