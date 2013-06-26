# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: superludist

# Downloading and compiling SuperLU DIST
# ------------------------------

# SUPERLU information
#http://crd.lbl.gov/~xiaoye/SuperLU/superlu_3.1.tar.gz
#http://crd-legacy.lbl.gov/~xiaoye/SuperLU/superlu_dist_3.0.tar.gz
DIRPKG=../pkg
SRCDIR=SuperLU_DIST_3.0
PACKAGE=superlu_dist_3.0.tar.gz
PACKAGE_PATH=$(DIRPKG)/$(PACKAGE)
SERVER=http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
INSTALL=../..
SUPERLU_DIST_VERSION=3.0
DISTSUPERLULIB=../$(INSTALL)/lib/libsuperlu_dist_3.0.a
WGET=wget

include make.inc

superludist: FAIT

# FFCS: more dependencies for parallel builds
FAIT:$(SRCDIR)
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

# FFCS: we need an error if make breaks
$(SRCDIR)/SRC/$(DSUPERLULIB): $(SRCDIR)
	cp make.inc $(SRCDIR)/make.inc
	cd $(SRCDIR)/SRC && $(MAKE)
# FFCS: we need more dependencies for parallel builds
install:$(SRCDIR)/SRC/$(DSUPERLULIB)
	mkdir -p $(SRCDIR)/$(INSTALL)/include/superludist
	cp $(SRCDIR)/SRC/*.h  $(SRCDIR)/$(INSTALL)/include/superludist
WHERE:
	echo superlu_dist LD -L@DIR@/lib -lsuperlu_dist_3.0  >$(SRCDIR)/$(INSTALL)/lib/WHERE.superlu_dist
	echo superlu_dist INCLUDE -I@DIR@/include/superludist  >> $(SRCDIR)/$(INSTALL)/lib/WHERE.superlu_dist


$(SRCDIR): $(PACKAGE_PATH)
	tar xvzf $(PACKAGE_PATH)
	cd $(SRCDIR)/SRC/; patch -p1 < ../../superludist_2.3.patch
	cd $(SRCDIR)/SRC/; patch -p3 < ../../superludist_3.0-printf.patch
	cd $(SRCDIR)/SRC/; patch -p2 < ../../superludist_3.0-cast_warning.patch
	cd $(SRCDIR)/SRC/; patch -p2 < ../../superludist_3.0-return_values.patch
	cd $(SRCDIR)/SRC/; patch -p2 < ../../superludist_3.0-operation_undefined.patch
	touch $(SRCDIR)

$(PACKAGE_PATH):
	-mkdir $(DIRPKG);
	cd $(DIRPKG);../../../../build/download $(SERVER)$(PACKAGE) $(PACKAGE_PATH)

# FFCS - make sure that all directories are cleaned. Thisis especially important under Windows because there is no
# compilation dependencies control there (see [[file:c:/cygwin/home/alh/ffcs/dist/configure.ac::dependency_tracking]])

clean:
	-rm -rf SuperLU_DIST*
	-rm  ../include/superludist/*.h
	-rm -r ../include/superludist/
	-rm -r ../lib/libsuperlu_dist*.a
	-rm -r make.inc
	-rm -rf FAIT FAIRE
