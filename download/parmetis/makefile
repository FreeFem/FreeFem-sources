# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Makefile.in

all-local: parmetis

# Downloading and compiling Tetgen
# ------------------------------
# http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.tar.gz
# Parmetis information
DIRPKG=../pkg
SRCDIR=ParMetis-$(parmetis_VERSION)
PACKAGE=$(DIRPKG)/ParMetis-$(parmetis_VERSION).tar.gz
SERVER=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/OLD
INSTALL=../..
parmetis_VERSION=3.1.1

parmetis: FAIRE

$(SRCDIR)/FAIT:$(SRCDIR)/tags Makefile.in
	-mkdir -p "../lib/parmetis"
	-mkdir -p "../include/parmetis"
	cp Makefile.in $(SRCDIR)
	cd $(SRCDIR)/METISLib; $(MAKE) 
	cd $(SRCDIR)/ParMETISLib;$(MAKE) 
	touch $(SRCDIR)/FAIT
install: 
#	cp $(SRCDIR)/libparmetis.a $(SRCDIR)/$(INSTALL)/lib
#	cp $(SRCDIR)/libmetis.a $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/parmetis.h ../include/parmetis
	cp $(SRCDIR)/METISLib/*.h ../include/parmetis
	sed  's;../parmetis.h;parmetis.h;' <$(SRCDIR)/METISLib/metis.h >../include/parmetis/metis.h
	-cd $(SRCDIR)/Programs;make 


WHERE:
	if [ -f $(SRCDIR)/FAIT ] ;then \
	make install ;\
	echo parmetis LD -L@DIR@/lib/parmetis -lparmetis -lmetis  >$(SRCDIR)/$(INSTALL)/lib/WHERE.parmetis ;\
	echo parmetis INCLUDE -I@DIR@/include/parmetis   >> $(SRCDIR)/$(INSTALL)/lib/WHERE.parmetis ;\
	fi

FAIRE: $(SRCDIR)/FAIT 
	$(MAKE) WHERE
	touch FAIRE

Makefile.in: ../../config.status Makefile-parmetis.in
	../../config.status  --file="Makefile.in:Makefile-parmetis.in"



$(SRCDIR)/tags: $(PACKAGE)
	tar xvzf $(PACKAGE)
	patch -p0  <parmetis-3.1.1.patch
	mv  $(SRCDIR)/Programs/Makefile $(SRCDIR)/Programs/Makefile-orig
	sed -e 's;BINDIR = .*$$;BINDIR = ../$(INSTALL)/bin;' \
            -e 's;../libparmetis.a;../$(INSTALL)/lib/parmetis/libparmetis.a;' \
            -e 's;../libmetis.a;../$(INSTALL)/lib/parmetis/libmetis.a;' \
            -e 's;-L[.][.];-L../$(INSTALL)/lib;' \
          <$(SRCDIR)/Programs/Makefile-orig \
          >$(SRCDIR)/Programs/Makefile
	mv $(SRCDIR)/ParMETISLib/Makefile $(SRCDIR)/ParMETISLib/Makefile-orig
	sed   -e 's;../libparmetis.a;../$(INSTALL)/lib/parmetis/libparmetis.a;' \
	      -e 's;../libmetis.a;../$(INSTALL)/lib/parmetis/libmetis.a;' \
          <$(SRCDIR)/ParMETISLib/Makefile-orig \
          >$(SRCDIR)/ParMETISLib/Makefile
	mv $(SRCDIR)/METISLib/Makefile $(SRCDIR)/METISLib/Makefile-orig
	sed   -e 's;../libparmetis.a;../$(INSTALL)/lib/parmetis/libparmetis.a;' \
	      -e 's;../libmetis.a;../$(INSTALL)/lib/parmetis/libmetis.a;' \
          <$(SRCDIR)/METISLib/Makefile-orig \
          >$(SRCDIR)/METISLib/Makefile
	touch $(SRCDIR)/tags

$(PACKAGE):
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`



clean-local:
	-cd $(SRCDIR) && $(MAKE) realclean  -C  $(SRCDIR)
	rm Makefile.in $(SRCDIR)/tags $(SRCDIR)/FAIT

clean: clean-local
	-rm $(SRCDIR)/$(INSTALL)/libparmetis.a
	-rm $(SRCDIR)/$(INSTALL)/libmetis.a
	-rm $(SRCDIR)/$(INSTALL)/include/parmetis.h
	-rm -rf  $(SRCDIR)/$(INSTALL)/include/metis
	-rm -rf $(SRCDIR)
.PHONY:$(SRCDIR)/$(INSTALL)