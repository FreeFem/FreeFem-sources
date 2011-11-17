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

blacs:FAIRE


FAIRE:$(SRCDIR)/FAIT
	$(MAKE) WHERE
	touch FAIRE

Bmake.inc: ../../config.status	Makefile Bmake-blacs.inc
	../../config.status  --file="Bmake.inc1:Bmake-blacs.inc"
	sed s/-DAdd__/-Df77IsF2C/ <Bmake.inc1 >Bmake.inc
	rm Bmake.inc1

$(SRCDIR)/FAIT: $(SRCDIR)
	#cd $(SRCDIR)/INSTALL/;$(MAKE) xintface xsyserrors xcmpi_sane xfmpi_sane xtc_CsameF77 xtc_UseMpich
	cd $(SRCDIR);$(MAKE) mpi 
	touch $(SRCDIR)/FAIT

WHERE: 
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	make install;  \
	echo blacs LD -L@DIR@/lib/blacs -lblacs_MPI-FREEFEM-0 	-lblacsF77init_MPI-FREEFEM-0	-lblacsCinit_MPI-FREEFEM-0  -lblacs_MPI-FREEFEM-0  >$(SRCDIR)/$(INSTALL)/lib/WHERE.blacs ;\
	fi

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
	-cd $(SRCDIR)/SRC/MPI &&  $(MAKE) clean -C $(SRCDIR)/SRC/MPI

clean: clean-local
	-rm Bmake.inc FAIRE
	-rm $(SRCDIR)/$(INSTALL)/lib/blacs/*.a
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)  compile install