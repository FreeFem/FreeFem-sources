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
#WGET  = curl -O 
#TOOLS
DIR1 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/TOOLS/LAPACK/*.o
DIR2 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/TOOLS/*.o
#PBLAS
DIR3 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/PBLAS/SRC/PBBLAS/*.o
DIR4 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/PBLAS/SRC/PTZBLAS/*.o
DIR5 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/PBLAS/SRC/PTOOLS/*.o
DIR6 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/PBLAS/SRC/*.o
#REDIST
DIR7 = $(abs_top_builddir)/download/scalapack/$(SRCDIR)/REDIST/SRC/*.o
#SRC LAPACK
DIR8 =  $(abs_top_builddir)/download/scalapack/$(SRCDIR)/SRC/*.o

scalapack: FAIRE

$(SRCDIR)/FAIT:
	$(MAKE) $(SRCDIR)/$(INSTALL)  
	touch $(SRCDIR)/FAIT

FAIRE: $(SRCDIR)/FAIT

SLmake.inc: ../../config.status	Makefile
	grep 'abs_top_builddir *=' ../Makefile >> dirflags
	cat dirflags SLmake-scalapack.inc > SLmake-tmp.inc
	../../config.status  --file="SLmake.inc:SLmake-tmp.inc"
	case `uname` in *CYGWIN*)  cp SLmake.inc SLmake.incn; sed "s/COPTIONS =/COPTIONS = -D__VC__/" <SLmake.incn > SLmake.inc; rm SLmake.incn ;; esac 
	rm dirflags
	rm SLmake-tmp.inc
$(SRCDIR)/$(INSTALL): $(SRCDIR)
	cp SLmake.inc $(SRCDIR)
	cd $(SRCDIR);$(MAKE) lib
	mkdir -p $(SRCDIR)/$(INSTALL)/include/scalapack
	cp $(SRCDIR)/SRC/*.h $(SRCDIR)/$(INSTALL)/include/scalapack/
	mkdir -p $(SRCDIR)/$(INSTALL)/lib/scalapack
	rm $(SRCDIR)/libscalapack.a
	$(ARCH) $(ARCHFLAGS) $(SRCDIR)/libscalapack.a $(DIR1) $(DIR2) $(DIR3) $(DIR4) $(DIR5) $(DIR6) $(DIR7) $(DIR8)
	$(RANLIB) $(SRCDIR)/libscalapack.a
	cp $(SRCDIR)/libscalapack.a $(SRCDIR)/$(INSTALL)/lib/scalapack/

$(SRCDIR): $(PACKAGE1) 
	gunzip -c $(PACKAGE1) | tar xvf -
	touch $(SRCDIR)

$(PACKAGE1):
	-mkdir $(DIRPKG)
	cd $(DIRPKG);$(WGET)  $(SERVER)/`basename $(PACKAGE1)`

clean-local:
	-cd $(SRCDIR); $(MAKE) cleanlib 

clean: clean-local
	rm SLmake.inc
	-rm -rf $(SRCDIR)/$(INSTALL)/include/scalapack/
	-rm -rf $(SRCDIR)/$(INSTALL)/lib/scalapack/
	-rm -rf $(SRCDIR)

.PHONY:$(SRCDIR)/$(INSTALL)