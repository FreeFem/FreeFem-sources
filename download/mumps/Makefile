# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Makefile.inc

all-local: mumps

# Downloading and compiling mumps
# ------------------------------
DIRPKG=../pkg
SRCDIR=MUMPS_$(VERSION)
PACKAGE1=$(DIRPKG)/MUMPS_$(VERSION).tar.gz
INSTALL=../..
VERSION=4.9.2

mumps: FAIRE

$(SRCDIR)/FAIT: $(SRCDIR)/tag-tar
	cp Makefile.inc $(SRCDIR)
	cd $(SRCDIR);make d z
	touch $(SRCDIR)/FAIT	
install: $(SRCDIR)/FAIT
	-mkdir -p $(SRCDIR)/$(INSTALL)/include
	cp $(SRCDIR)/include/*.h $(SRCDIR)/$(INSTALL)/include/
	-mkdir -p $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/lib/*.a $(SRCDIR)/$(INSTALL)/lib/

WHERE: 
	if [ -f $(SRCDIR)/FAIT  ] ;then \
	$(MAKE) install ; \
	echo "mumps LD '-L@DIR@/lib   -ldmumpsFREEFEM -lzmumpsFREEFEM  -lmumps_commonFREEFEM -lpordFREEFEM '"  >$(SRCDIR)/$(INSTALL)/lib/WHERE.mumps ;\
	echo "mumps INCLUDE '-I@DIR@/include'">> $(SRCDIR)/$(INSTALL)/lib/WHERE.mumps ;\
	fi


FAIRE:
	$(MAKE) install WHERE 
	touch FAIRE

Makefile.inc: ../../config.status	Makefile Makefile-mumps-$(VERSION).inc
	../../config.status  --file="Makefile.inc:Makefile-mumps-$(VERSION).inc"

$(SRCDIR)/$(INSTALL): $(SRCDIR)/tag-tar

$(SRCDIR)/tag-tar:
	tar xvzf $(PACKAGE1)
	patch -d MUMPS_4.9.2   -p 1    <MUMPS_4.9.2.patch
	touch $(SRCDIR)/tag-tar

clean-local:
	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR) 

clean: clean-local
	rm Makefile.inc FAIRE
	-rm -rf $(SRCDIR)/$(INSTALL)/include/*mumps*
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/lib*mumps* $(SRCDIR)/$(INSTALL)/lib/libpord*.a
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)