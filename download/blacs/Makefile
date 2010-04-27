# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include Bmake.inc

all-local: blacs 

# Downloading and compiling BLACS
# ------------------------------
# http://www.netlib.org/blacs/
# Hips information
DIRPKG=../pkg
SRCDIR=BLACS
PACKAGE1=$(DIRPKG)/mpiblacs.tgz
PACKAGE2=$(DIRPKG)/mpiblacs-patch03.tgz
PACKAGE3=$(DIRPKG)/blacstester.tgz
SERVER=http://www.netlib.org/blacs/
INSTALL=../..

SYSERRORS = 

blacs: FAIRE


FAIRE: $(SRCDIR)/FAIT install

Bmake.inc: ../../config.status	Makefile Bmake-blacs.inc
	../../config.status  --file="Bmake.inc:Bmake-blacs.inc"


$(SRCDIR)/FAIT: $(SRCDIR)
	#cd $(SRCDIR)/INSTALL/;$(MAKE) xintface xsyserrors xcmpi_sane xfmpi_sane xtc_CsameF77 xtc_UseMpich
	cd $(SRCDIR);$(MAKE) mpi 
	touch $(SRCDIR)/FAIT

install: $(SRCDIR)/FAIT
	mkdir -p $(SRCDIR)/$(INSTALL)/lib/blacs
	cp $(SRCDIR)/LIB/*.a $(SRCDIR)/$(INSTALL)/lib/blacs/

$(SRCDIR): $(PACKAGE1) $(PACKAGE2) $(PACKAGE3)
	gunzip -c $(PACKAGE1) | tar xvf -
	gunzip -c $(PACKAGE2) | tar xvf -
	gunzip -c $(PACKAGE3) | tar xvf -
	cd $(SRCDIR)/SRC/MPI; patch -p1 < ../../../BLACS.patch
	cp Bmake.inc $(SRCDIR)
	mv $(SRCDIR)/SRC/MPI/Makefile $(SRCDIR)/SRC/MPI/Makefile.tmp
	sed -e 's;\.C;\.oo;g'\
	< $(SRCDIR)/SRC/MPI/Makefile.tmp \
	> $(SRCDIR)/SRC/MPI/Makefile
	touch $(SRCDIR)

$(PACKAGE1):
	-mkdir $(DIRPKG)
	cd $(DIRPKG);$(WGET)  $(SERVER)/`basename $(PACKAGE1)`

$(PACKAGE2):
	-mkdir $(DIRPKG)
	cd $(DIRPKG);$(WGET) $(SERVER)/`basename $(PACKAGE2)`

$(PACKAGE3):
	-mkdir $(DIRPKG)
	cd $(DIRPKG);$(WGET) $(SERVER)/`basename $(PACKAGE3)`

clean-local:
	-cd $(SRCDIR)/SRC/MPI; make clean 

clean: clean-local
	-rm Bmake.inc
	-rm $(SRCDIR)/$(INSTALL)/lib/blacs/*.a
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)  compile install