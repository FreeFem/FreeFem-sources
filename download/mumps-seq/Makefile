# Downloading and compiling extra libraries
# -----------------------------------------

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
#
#	ALH - 22/5/12 - Mumps has difficulties compiling d & z in parallel
#
	cd $(SRCDIR) && make d
	cd $(SRCDIR) && make z
	touch $(SRCDIR)/FAIT	
install.done: $(SRCDIR)/FAIT
	-mkdir -p $(SRCDIR)/$(INSTALL)/include/libseq
	cp $(SRCDIR)/include/*.h $(SRCDIR)/$(INSTALL)/include/libseq
	cp $(SRCDIR)/libseq/*.h  $(SRCDIR)/$(INSTALL)/include/libseq
	-mkdir -p $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/lib/*.a $(SRCDIR)/$(INSTALL)/lib/
	cp $(SRCDIR)/libseq/libmpiseqFREEFEM-SEQ.a $(SRCDIR)/$(INSTALL)/lib/
	touch $@
clean::
	-rm *.done

# FFCS - install and WHERE need to be done sequentially, even in parallel builds
WHERE.done: install.done
	echo mumps-seq LD -L@DIR@/lib   -ldmumpsFREEFEM-SEQ -lzmumpsFREEFEM-SEQ  -lmumps_commonFREEFEM-SEQ -lpordFREEFEM-SEQ -lpthread >$(SRCDIR)/$(INSTALL)/lib/WHERE.mumpsseq ;
	echo mumps-seq INCLUDE -I@DIR@/include/libseq  >> $(SRCDIR)/$(INSTALL)/lib/WHERE.mumpsseq ;
	echo libseq LD -L@DIR@/lib  -lmpiseqFREEFEM-SEQ   >>$(SRCDIR)/$(INSTALL)/lib/WHERE.mumpsseq ;
	echo libseq INCLUDE -I@DIR@/include/libseq>> $(SRCDIR)/$(INSTALL)/lib/WHERE.mumpsseq ;
	touch $@

FAIRE-$(VERSION): install.done WHERE.done
	touch $@

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

# FFCS - $(SRCDIR)/$(INSTALL) is not valid if $(SRCDIR) is not there
clean:: clean-local
	-rm Makefile.inc FAIRE* *~ 
	-rm -rf ../include/*mumps*
	-rm -rf ../lib/lib*mumps* ../lib/libpord*.a ../lib/libmpiseq*.a ../lib/WHERE.mumpsseq
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)
