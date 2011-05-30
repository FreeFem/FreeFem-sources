# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: superludist

# Downloading and compiling SuperLU DIST
# ------------------------------

# SUPERLU information
#http://crd.lbl.gov/~xiaoye/SuperLU/superlu_3.1.tar.gz
DIRPKG=../pkg
SRCDIR=SuperLU_DIST_2.3
PACKAGE=superlu_dist_2.3.tar.gz
PACKAGE_PATH=$(DIRPKG)/$(PACKAGE)
SERVER=http://crd.lbl.gov/~xiaoye/SuperLU/
INSTALL=../..
SUPERLU_DIST_VERSION=2.3
DISTSUPERLULIB=../$(INSTALL)/lib/libsuperlu_dist_2.3.a
WGET=wget

include make.inc

superludist: FAIT

FAIT:	
	$(MAKE) FAIRE install WHERE
	touch FAIT


FAIRE: $(SRCDIR)/SRC/$(DSUPERLULIB)  
	touch FAIRE

make.inc: ../../config.status  Makefile
	grep 'abs_top_builddir *=' ../Makefile >> dirflags
	cat dirflags make-superlu.inc > make-superlu-tmp.inc
	rm dirflags
	../../config.status  --file="make.inc:make-superlu-tmp.inc"
	case 'uname' in *CYGWIN_) cp make.inc Makefile.inn; sed "s/COPTIONS =/COPTIONS = -D__VC__/" <Makefile.inn > Makefile.in; rm Makefile.inn ;; esac
	rm make-superlu-tmp.inc

$(SRCDIR)/SRC/$(DSUPERLULIB): $(SRCDIR)
	cp make.inc $(SRCDIR)/make.inc
	cd $(SRCDIR)/SRC; $(MAKE)
install:
	mkdir -p $(SRCDIR)/$(INSTALL)/include/superludist
	cp $(SRCDIR)/SRC/*.h  $(SRCDIR)/$(INSTALL)/include/superludist
WHERE:
	echo superlu_dist LD -L@DIR@/lib -lsuperlu_dist_2.3  >$(SRCDIR)/$(INSTALL)/lib/WHERE.superlu_dist
	echo superlu_dist INCLUDE -I@DIR@/include/superludist  >> $(SRCDIR)/$(INSTALL)/lib/WHERE.superlu_dist


$(SRCDIR): $(PACKAGE_PATH)
	tar xvzf $(PACKAGE_PATH)
	cd $(SRCDIR)/SRC/; patch -p1 < ../../superludist_2.3.patch
	touch $(SRCDIR)

$(PACKAGE_PATH):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);$(WGET)  $(SERVER)$(PACKAGE)

clean-local:
	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR)

clean: clean-local
	-rm -rf $(SRCDIR)
	-rm  $(SRCDIR)/$(INSTALL)/include/superludist/*.h
	-rm -r $(SRCDIR)/$(INSTALL)/include/superludist/
	-rm -r $(SRCDIR)/$(INSTALL)/lib/libsuperlu_dist*.a
	-rm -r make.inc
	-rm -rf FAIT
