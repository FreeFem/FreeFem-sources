# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Makefile.in

all-local: parms

# Downloading and compiling parms
# -------------------------------
# http://www-users.cs.umm.edu/~saad/software/pARMS/pARMS_2.2.php
# Parms information
DIRPKG=../pkg
SRCDIR=pARMS_$(parms_VERSION)
PACKAGE=$(DIRPKG)/parms_$(parms_VERSION).tar.gz
SERVER= http://www-users.cs.umm.edu/~saad/software/pARMS/pARMS_2.2.php
INSTALL=../..
parms_VERSION=2.2

parms: FAIRE

$(SRCDIR)/FAIT:
	$(MAKE) $(SRCDIR)/$(INSTALL)  
	touch $(SRCDIR)/FAIT

FAIRE: $(SRCDIR)/FAIT

Makefile.in: ../../config.status	Makefile makefile-parms.in
	../../config.status  --file="Makefile.in:makefile-parms.in"
	case `uname` in *CYGWIN*)  cp Makefile.in Makefile.inn; sed "s/COPTIONS =/COPTIONS = -D__VC__/" <Makefile.inn > Makefile.in; rm Makefile.inn ;; esac 

$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cp Makefile.in $(SRCDIR)
	cd $(SRCDIR); $(MAKE) ./LIB/libparms$(DBG).a
	mkdir -p $(SRCDIR)/$(INSTALL)/include/parms
	cp $(SRCDIR)/INCLUDE/*.h $(SRCDIR)/$(INSTALL)/include/parms/
	cp $(SRCDIR)/LIB/*.a $(SRCDIR)/$(INSTALL)/lib/

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)
	touch $(SRCDIR)

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)   $(SERVER)/`basename $(PACKAGE)`

clean-local:
	-cd $(SRCDIR)/src; make cleanall

clean: clean-local
	-rm Makefile.in
	-rm -rf $(SRCDIR)
.PHONY:$(SRCDIR)/$(INSTALL)