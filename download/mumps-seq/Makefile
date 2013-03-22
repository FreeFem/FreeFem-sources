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
	cp $(SRCDIR)/include/*.h ../include/libseq
	cp $(SRCDIR)/libseq/*.h  ../include/libseq
	-mkdir -p ../lib
	cp $(SRCDIR)/lib/*.a ../lib/
	cp $(SRCDIR)/libseq/libmpiseqFREEFEM-SEQ.a ../lib/

WHERE: 
	if [ -f $(SRCDIR)/FAIT  ] ;then \
	$(MAKE) install ; \
	echo mumps-seq LD -L@DIR@/lib   -ldmumpsFREEFEM-SEQ -lzmumpsFREEFEM-SEQ  -lmumps_commonFREEFEM-SEQ -lpordFREEFEM-SEQ -lpthread    >../lib/WHERE.mumpsseq ;\
	echo mumps-seq INCLUDE -I@DIR@/include/libseq  >> ../lib/WHERE.mumpsseq ;\
	echo libseq LD -L@DIR@/lib  -lmpiseqFREEFEM-SEQ   >>../lib/WHERE.mumpsseq ;\
	echo libseq INCLUDE -I@DIR@/include/libseq>> ../lib/WHERE.mumpsseq ;\
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
	rm Makefile.inc FAIRE* *~ ../lib/WHERE.mumpsseq
	-rm -rf ../include/*mumps*
	-rm -rf ../lib/lib*mumps* ../lib/libpord*.a ../lib/libmpiseq*.a
	-rm -rf $(SRCDIR)

clean: clean-local


.PHONY:$(SRCDIR)/$(INSTALL)