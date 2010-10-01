# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Makefile.in

all-local: parms

# Downloading and compiling parms
# -------------------------------
# http://www-users.cs.umn.edu/~saad/software/pARMS/pARMS_2.2.php
# Parms information
DIRPKG=../pkg
SRCDIR=pARMS_$(parms_VERSION)
PACKAGE=$(DIRPKG)/parms_$(parms_VERSION).tar.gz
SERVER= http://www-users.cs.umn.edu/~saad/software/pARMS/pARMS_2.2.php
INSTALL=../..
parms_VERSION=2.2

parms: FAIRE

$(SRCDIR)/FAIT:
	$(MAKE) $(SRCDIR)/$(INSTALL)  
	touch $(SRCDIR)/FAIT

FAIRE: $(SRCDIR)/FAIT WHERE

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
	cd $(DIRPKG);$(WGET)   $(SERVER)

clean-local:
	-make cleanall -C $(SRCDIR)/src

clean: 
	-rm Makefile.in
	-rm -rf $(SRCDIR)

WHERE:
	echo parms LD '-L@DIR@/lib -lparms-O3'  >$(SRCDIR)/$(INSTALL)/lib/WHERE.parms
	echo parms INCLUDE '-I@DIR@/include/parms '>> $(SRCDIR)/$(INSTALL)/lib/WHERE.parms

.PHONY:$(SRCDIR)/$(INSTALL)

