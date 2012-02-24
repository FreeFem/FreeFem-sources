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
	-mkdir -p $(SRCDIR)/$(INSTALL)/include/libseq
	cp $(SRCDIR)/include/*.h $(SRCDIR)/$(INSTALL)/include/libseq
	cp $(SRCDIR)/libseq/*.h  $(SRCDIR)/$(INSTALL)/include/libseq
	-mkdir -p $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/lib/*.a $(SRCDIR)/$(INSTALL)/lib/
	cp $(SRCDIR)/libseq/libmpiseqFREEFEM-SEQ.a $(SRCDIR)/$(INSTALL)/lib/

WHERE: 
	if [ -f $(SRCDIR)/FAIT  ] ;then \
	$(MAKE) install ; \
	echo mumps-seq LD -L@DIR@/lib   -ldmumpsFREEFEM-SEQ -lzmumpsFREEFEM-SEQ  -lmumps_commonFREEFEM-SEQ -lpordFREEFEM-SEQ   >$(SRCDIR)/$(INSTALL)/lib/WHERE.mumps ;\
	echo mumps-seq INCLUDE -I@DIR@/include/libseq  >> $(SRCDIR)/$(INSTALL)/lib/WHERE.mumpsseq ;\
	echo libseq LD -L@DIR@/lib  -lmpiseqFREEFEM-SEQ   >>$(SRCDIR)/$(INSTALL)/lib/WHERE.mumpsseq ;\
	echo libseq INCLUDE -I@DIR@/include/libseq>> $(SRCDIR)/$(INSTALL)/lib/WHERE.mumpsseq ;\
	fi


FAIRE-$(VERSION):
	$(MAKE) install WHERE 
	touch FAIRE-$(VERSION)

Makefile.inc: ../../config.status	Makefile Makefile-mumps-$(VERSION).inc
	../../config.status  --file="Makefile.inc:Makefile-mumps-$(VERSION).inc"

$(SRCDIR)/$(INSTALL): $(SRCDIR)/tag-tar

$(SRCDIR)/tag-tar:$(PACKAGE1)
	tar xvzf $(PACKAGE1)
#	patch -d MUMPS_4.9.2   -p 1    <MUMPS_4.9.2.patch
	touch $(SRCDIR)/tag-tar

$(PACKAGE1):
	cd `dirname $@`; $(WGET)   http://graal.ens-lyon.fr/MUMPS/`basename $(PACKAGE1)`
clean-local:
	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR) 

clean: clean-local
	rm Makefile.inc FAIRE* *~ 
	-rm -rf $(SRCDIR)/$(INSTALL)/include/*mumps*
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/lib*mumps* $(SRCDIR)/$(INSTALL)/lib/libpord*.a $(SRCDIR)/$(INSTALL)/lib/libmpiseq*.a
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)