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
VERSION=4.10.0

mumps: FAIRE-$(VERSION)

$(SRCDIR)/FAIT: $(SRCDIR)/tag-tar
	cp Makefile.inc $(SRCDIR)
	cd $(SRCDIR);make d z
	touch $(SRCDIR)/FAIT	
install: $(SRCDIR)/FAIT
	-mkdir -p ../include/libseq
	cp $(SRCDIR)/include/*.h ../include/
	cp $(SRCDIR)/libseq/*.h  ../include/libseq
	-mkdir -p ../lib
	cp $(SRCDIR)/lib/*.a ../lib/
	cp $(SRCDIR)/libseq/libmpiseqFREEFEM.a ../lib/

WHERE: 
	if [ -f $(SRCDIR)/FAIT  ] ;then \
	$(MAKE) install ; \
	echo mumps LD -L@DIR@/lib   -ldmumpsFREEFEM -lzmumpsFREEFEM  -lmumps_commonFREEFEM -lpordFREEFEM  -lpthread  >../lib/WHERE.mumps ;\
	echo mumps INCLUDE -I@DIR@/include  >> ../lib/WHERE.mumps ;\
	fi


FAIRE-$(VERSION):
	$(MAKE) install WHERE 
	touch FAIRE-$(VERSION)

Makefile.inc: ../../config.status	Makefile Makefile-mumps-$(VERSION).inc
	../../config.status  --file="Makefile.inc:Makefile-mumps-$(VERSION).inc"

$(SRCDIR)/$(INSTALL): $(SRCDIR)/tag-tar

$(SRCDIR)/tag-tar:$(PACKAGE1)
	tar xvzf $(PACKAGE1)
	patch -d MUMPS_$(VERSION)   -p 1    <MUMPS_$(VERSION).patch
	touch $(SRCDIR)/tag-tar

$(PACKAGE1):
	cd `dirname $@`; $(WGET)   http://graal.ens-lyon.fr/MUMPS/`basename $(PACKAGE1)`

clean-local:
	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR) 
	rm Makefile.inc FAIRE* ../lib/WHERE.mumps
	-rm -rf ../include/*mumps*
	-rm -rf ../lib/lib*mumps* ../lib/libpord*.a ../lib/libmpiseq*.a
	-rm -rf $(SRCDIR)

clean: clean-local


.PHONY:$(SRCDIR)/$(INSTALL)