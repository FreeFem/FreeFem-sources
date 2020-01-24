# Downloading and compiling extra libraries
# -----------------------------------------

all-local: parmmg

include ff-flags

# Downloading and compiling parmmg
# -------------------------------
#
DIRPKG= ../pkg
SRCDIR= .
PACKAGE=$(DIRPKG)/parmmg.zip
INSTALL=.
parmmg_VERSION=

parmmg:WHERE.done

FAIT.done:tag-tar
	-mkdir build
	cd build && cmake ../parmmg-sources \
	-DCMAKE_C_COMPILER=$(CC) \
	-DCMAKE_C_FLAGS="-O3 -fPIC" \
	-DDOWNLOAD_MMG=OFF \
	-DDOWNLOAD_METIS=OFF \
	-DMMG_DIR=../../mmg/mmg-sources \
	-DMMG_BUILDDIR=../../mmg/build \
	-DSCOTCH_DIR=$(scotch_dir) \
	-DMETIS_DIR=$(metis_dir) \
	-DCMAKE_BUILD_TYPE=Release && $(MAKE)
	touch FAIT.done

install.done:FAIT.done
	cp -r build/include/parmmg ../include/parmmg
	cp build/lib/libparmmg.a ../lib
	touch $@
clean::
	-rm install.done FAIT.done

parmmg:$(PACKAGE)

install:install.done WHERE.done

WHERE.done: install.done
	echo parmmg  LD -L@DIR@/lib -lparmmg  >../lib/WHERE.parmmg ;
	echo parmmg INCLUDE -I@DIR@/include/>> ../lib/WHERE.parmmg ;
	echo build WHERE ./lib/WHERE.parmmg ;
	touch $@
clean::
	-rm WHERE.done

FAIRE: FAIT.done install.done

tag-tar: $(PACKAGE)
	-rm -rf  ParMmg-*
	unzip $(PACKAGE) && mv ParMmg-* parmmg-sources
	touch tag-tar

$(PACKAGE):
	../getall -o parmmg -a

clean::
	-rm ff-flags
	-rm FAIT* 
	-rm -rf ParMmg-* parmmg-sources build
	-rm flags-* tag-tar*

ff-flags: ../Makefile Makefile
	grep 'abs_top_builddir *=' ../Makefile > ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags
	grep 'STD_LIBS *=' ../Makefile >> ff-flags
	grep 'CNOFLAGS *=' ../Makefile >> ff-flags
	grep 'FF_prefix_petsc *=' ../Makefile >> ff-flags
	-awk '/^scotch LD /&&/petsc/ {print "scotch_dir=$$(FF_prefix_petsc)/real" }' ../../plugin/seq/WHERE_LIBRARY-config >> ff-flags
	-awk '/^scotch LD /{print "scotch_dir=../../3rdparty/" }' ../lib/WHERE.scotch >> ff-flags
	-awk '/^metis LD /&&/petsc/ {print "metis_dir=$$(FF_prefix_petsc)/real" }' ../../plugin/seq/WHERE_LIBRARY-config >> ff-flags
	-awk '/^metis LD /{print "metis_dir=../../3rdparty/" }' ../lib/WHERE.metis >> ff-flags
	
.PHONY:$(SRCDIR)/$(INSTALL)
