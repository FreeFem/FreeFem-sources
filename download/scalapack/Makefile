# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$
include SLmake.inc

all-local: scalapack

# Downloading and compiling scalapack
# ------------------------------
# http://www.netlib.org/scalapack/
# Hips information
DIRPKG=../pkg
SRCDIR=scalapack-1.8.0
PACKAGE1=$(DIRPKG)/scalapack-1.8.0.tgz
SERVER=http://www.netlib.org/scalapack/
INSTALL=../..
DIR1 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/TOOLS/LAPACK/*.o
DIR2 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/TOOLS/*.o
DIR3 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/PBLAS/SRC/PBBLAS/*.o
DIR4 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/PBLAS/SRC/PTZBLAS/*.o
DIR5 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/PBLAS/SRC/PTOOLS/*.o
DIR6 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/PBLAS/SRC/*.o
DIR7 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/REDIST/SRC/*.o
DIR8 =  $(abs_top_builddir)/download/scalapack/$(SRCDIR)/SRC/*.o

scalapack: FAIRE

$(SRCDIR)/FAIT:$(SRCDIR)/tag-tar
	cp SLmake.inc $(SRCDIR)
	cd $(SRCDIR);$(MAKE) lib
	rm $(SRCDIR)/libscalapack.a
	$(ARCH) $(ARCHFLAGS) $(SRCDIR)/libscalapack.a $(DIR1) $(DIR2) $(DIR3) $(DIR4) $(DIR5) $(DIR6) $(DIR7) $(DIR8)
	$(RANLIB) $(SRCDIR)/libscalapack.a
	touch $(SRCDIR)/FAIT
install:
	mkdir -p $(SRCDIR)/$(INSTALL)/include
	cp $(SRCDIR)/SRC/*.h $(SRCDIR)/$(INSTALL)/include
	mkdir -p $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/libscalapack.a $(SRCDIR)/$(INSTALL)/lib

WHERE:
	if [ -f $(SRCDIR)/FAIT ] ; then \
	make install;  \
	echo scalapack LD -L@DIR@/lib  -lscalapack   >$(SRCDIR)/$(INSTALL)/lib/WHERE.scalapack ;\
	echo scalapack INCLUDE -I@DIR@/include  >> $(SRCDIR)/$(INSTALL)/lib/WHERE.scalapack ;\
	fi

FAIRE: $(SRCDIR)/FAIT
	$(MAKE) WHERE
	touch FAIRE 

SLmake.inc: ../../config.status	Makefile
	../../config.status  --file="SLmake.inc:SLmake-scalapack.inc"

$(SRCDIR)/tag-tar: $(PACKAGE1)
	gunzip -c $(PACKAGE1) | tar xvf -
	touch $(SRCDIR)/tag-tar

$(PACKAGE1):
	-mkdir $(DIRPKG)
	cd $(DIRPKG);$(WGET)  $(SERVER)/`basename $(PACKAGE1)`

clean-local:
	-cd $(SRCDIR) && $(MAKE) cleanlib  -C $(SRCDIR)

clean: clean-local
	rm SLmake.inc FAIRE FAIT
	-rm -rf $(SRCDIR)/$(INSTALL)/include/{pblas.h,pxsyevx.h,tools.h}
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/libscalapack.a
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)