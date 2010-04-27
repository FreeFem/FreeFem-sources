# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Makefile.inc

all-local: scotch

# Downloading and compiling scotch
# -------------------------------
# http://gforge.inria.fr/frs/download.php/23391/scotch_5.1.7_esmumps.tar.gz
# Scotch information
DIRPKG=../pkg
SRCDIR=scotch_$(scotch_VERSION)
PACKAGE=$(DIRPKG)/scotch_$(scotch_VERSION_LOCAL).tar.gz
SERVER=http://gforge.inria.fr/frs/download.php/23391
INSTALL=../..
scotch_VERSION=5.1_esmumps
scotch_VERSION_LOCAL=5.1.6_esmumps

scotch: FAIRE

$(SRCDIR)/AFAIRE:
	$(MAKE) install

FAIRE: $(SRCDIR)/FAIT WHERE

Makefile.inc: ../../config.status Makefile-scotch.inc Makefile
	../../config.status  --file="Makefile.inc:Makefile-scotch.inc"
	case `uname` in *CYGWIN*)  cp Makefile.inc Makefile.inn; sed "s/COPTIONS =/COPTIONS = -D__VC__/" <Makefile.inn > Makefile.inc; rm Makefile.inn ;; esac 

$(SRCDIR)/FAIT: $(SRCDIR)
	cp Makefile.inc $(SRCDIR)/src/
	cd $(SRCDIR)/src/; make scotch
	cd $(SRCDIR)/src/; make ptscotch
	touch $(SRCDIR)/FAIT

install:$(SRCDIR)/FAIT
	mkdir -p $(SRCDIR)/$(INSTALL)/include/scotch
	cp $(SRCDIR)/include/*.h $(SRCDIR)/$(INSTALL)/include/scotch/
	cp $(SRCDIR)/lib/*.a $(SRCDIR)/$(INSTALL)/lib/
WHERE:
	if [ -f  $(SRCDIR)/FAIT ] ; then \
	echo "ptscotch LD '-L@DIR@/lib  -lptesmumps -lptscotch -lptscotcherr '"  >$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch ;\
	echo "ptscotch INCLUDE '-I@DIR@/include/scotch'">> $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch ;\
	echo "scotch LD '-L@DIR@/lib   -lscotch -lscotcherr '"  >>$(SRCDIR)/$(INSTALL)/lib/WHERE.scotch ;\
	echo "scotch INCLUDE '-I@DIR@/include/scotch'">> $(SRCDIR)/$(INSTALL)/lib/WHERE.scotch ;\
	fi

$(SRCDIR): $(PACKAGE)
	tar xvzf $(PACKAGE)
	cd $(SRCDIR)/src/libscotch/; patch -p1 < ../../../scotch_5.1_esmumps.patch;
	touch $(SRCDIR)

$(PACKAGE):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)  --no-check-certificate  $(SERVER)/`basename $(PACKAGE)`

clean-local:
	-cd $(SRCDIR)/src; make realclean 
clean: clean-local
	-rm Makefile.inc
	-rm -rf $(SRCDIR)
.PHONY:$(SRCDIR)/$(INSTALL)