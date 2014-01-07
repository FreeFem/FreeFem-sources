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

# FFCS: make sure that this makefile fails if the BLACS do not compile. But run links2files when the compilation fails because of
# a symbolic link (under Windows). Repeat it only by hand to avoid infinite loops!
$(SRCDIR)/FAIT: $(SRCDIR)
#	cd $(SRCDIR)/INSTALL/;$(MAKE) xintface xsyserrors xcmpi_sane xfmpi_sane xtc_CsameF77 xtc_UseMpich
	-cd $(SRCDIR) && $(MAKE) mpi
	../../../../build/links2files
	-cd $(SRCDIR) && $(MAKE) mpi
	../../../../build/links2files
	cd $(SRCDIR) && $(MAKE) mpi
	touch $(SRCDIR)/FAIT

WHERE: 
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	make install;  \
	echo blacs LD -L@DIR@/lib/blacs -lblacs_MPI-FREEFEM-0 	-lblacsF77init_MPI-FREEFEM-0	-lblacsCinit_MPI-FREEFEM-0  -lblacs_MPI-FREEFEM-0  >$(SRCDIR)/$(INSTALL)/lib/WHERE.blacs ;\
	fi

install: $(SRCDIR)/FAIT
	mkdir -p $(SRCDIR)/$(INSTALL)/lib/blacs
	cp $(SRCDIR)/LIB/*.a $(SRCDIR)/$(INSTALL)/lib/blacs/

# FFCS: some files in the untarred archives are symbolic links that need to be converted for MinGW compilers
$(SRCDIR): $(PACKAGE1) $(PACKAGE2) $(PACKAGE3)
	gunzip -c $(PACKAGE1) | tar xvf -
	gunzip -c $(PACKAGE2) | tar xvf -
	gunzip -c $(PACKAGE3) | tar xvf -
	patch -p0 < BLACS_gridinit_.c-return-values.patch
	cd $(SRCDIR)/SRC/MPI; patch -p1 < ../../../BLACS.patch
	cp Bmake.inc $(SRCDIR)
	mv $(SRCDIR)/SRC/MPI/Makefile $(SRCDIR)/SRC/MPI/Makefile.tmp
	sed -e 's;\.C;\.oo;g'  -e 's@: $$(MPIINCdir)/mpif.h@:@' \
	< $(SRCDIR)/SRC/MPI/Makefile.tmp \
	> $(SRCDIR)/SRC/MPI/Makefile
	../../../../build/links2files
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

# FFCS: only run make clean if cd to SRCDIR worked, otherwise this is an infinite loop.
clean-local:
	-cd $(SRCDIR)/SRC/MPI &&  $(MAKE) clean -C $(SRCDIR)/SRC/MPI
	-rm Bmake.inc FAIRE
	-rm $(SRCDIR)/$(INSTALL)/lib/blacs/*.a
	-rm -rf $(SRCDIR)

clean: clean-local


.PHONY:$(SRCDIR)/$(INSTALL)  compile install