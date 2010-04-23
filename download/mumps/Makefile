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
	-mkdir -p $(SRCDIR)/$(INSTALL)/include
	cp $(SRCDIR)/include/*.h $(SRCDIR)/$(INSTALL)/include/
	-mkdir -p $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/lib/*.a $(SRCDIR)/$(INSTALL)/lib/
	touch $(SRCDIR)/FAIT
WHERE:
	echo "mumps LD '-L@DIR@//lib  -lpord -lmumps_common -ldmumps -lzmumps -lpord'"  >$(SRCDIR)/$(INSTALL)/lib/WHERE.mumps
	echo "mumps INCLUDE '-I@DIR@/include'">> $(SRCDIR)/$(INSTALL)/lib/WHERE.mumps


FAIRE: $(SRCDIR)/FAIT WHERE

Makefile.inc: ../../config.status	Makefile
	../../config.status  --file="Makefile.inc:Makefile-mumps-$(VERSION).inc"

$(SRCDIR)/$(INSTALL): $(SRCDIR)/tag-tar

$(SRCDIR)/tag-tar:
	tar xvzf $(PACKAGE1)
	patch -d MUMPS_4.9.2   -p 1    <MUMPS_4.9.2.patch
	touch $(SRCDIR)/tag-tar

clean-local:
	-cd $(SRCDIR); make clean

clean: clean-local
	rm Makefile.inc
	-rm -rf $(SRCDIR)/$(INSTALL)/include/*mumps*
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/lib*mumps* $(SRCDIR)/$(INSTALL)/lib/libpord*.a
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)