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

$(SRCDIR)/FAIT:
	$(MAKE) $(SRCDIR)/$(INSTALL)  
	touch $(SRCDIR)/FAIT

FAIRE: $(SRCDIR)/FAIT

Makefile.inc: ../../config.status	Makefile
	grep 'abs_top_builddir *=' ../Makefile >> dirflags
	cat dirflags  Makefile-mumps-$(VERSION).inc > Makefile-mumps-tmp.inc
	../../config.status  --file="Makefile.inc:Makefile-mumps-tmp.inc"
	case `uname` in *CYGWIN*)  cp Makefile.inc Makefile.incn; sed "s/COPTIONS =/COPTIONS = -D__VC__/" <Makefile.incn > Makefile.inc; rm Makefile.incn ;; esac
	rm dirflags
	rm Makefile-mumps-tmp.inc

$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cp Makefile.inc $(SRCDIR)
	cd $(SRCDIR);make d z
	mkdir -p $(SRCDIR)/$(INSTALL)/include/mumps
	cp $(SRCDIR)/include/*.h $(SRCDIR)/$(INSTALL)/include/mumps/
	mkdir -p $(SRCDIR)/$(INSTALL)/lib/mumps
	cp $(SRCDIR)/lib/*.a $(SRCDIR)/$(INSTALL)/lib/mumps/

$(SRCDIR):
	tar xvzf $(PACKAGE1)
	touch $(SRCDIR)

clean-local:
	-cd $(SRCDIR); make clean

clean: clean-local
	rm Makefile.inc
	-rm -rf $(SRCDIR)/$(INSTALL)/include/mumps/
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/mumps/
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)