# Downloading and compiling extra libraries
# -----------------------------------------

all-local: mmg

include ff-flags

# Downloading and compiling mmg
# -------------------------------
#
DIRPKG= ../pkg
SRCDIR= .
PACKAGE=$(DIRPKG)/mmg.zip
INSTALL=.
mmg_VERSION=

mmg:WHERE.done

FAIT.done:tag-tar
	-mkdir build
	cd build && cmake ../mmg-sources \
	-DCMAKE_C_COMPILER=$(CC) \
	-DCMAKE_C_FLAGS="-O3 -fPIC" \
	-DM_LIB="-lm" \
	-DCMAKE_BUILD_TYPE=Release && $(MAKE)
	touch FAIT.done

install.done:FAIT.done
	cp -r build/include/mmg ../include/mmg
	cp build/lib/libmmg.a ../lib
	touch $@
clean::
	-rm install.done FAIT.done

mmg:$(PACKAGE)

install:install.done WHERE.done

# FFCS - keep simplest makefile structure for automatic recompilations
WHERE.done: install.done
	echo mmg  LD -L@DIR@/lib -lmmg  >../lib/WHERE.mmg ;
	echo mmg INCLUDE -I@DIR@/include/>> ../lib/WHERE.mmg ;
	echo build WHERE ./lib/WHERE.mmg ;
	touch $@
clean::
	-rm WHERE.done

FAIRE: FAIT.done install.done

tag-tar: $(PACKAGE)
	-rm -rf  mmg-*
	unzip $(PACKAGE) && mv mmg-* mmg-sources
	patch -p1 <patch-mmg.txt
	touch tag-tar

$(PACKAGE):
	../getall -o mmg -a

clean::
	-rm ff-flags
	-rm FAIT* 
	-rm -rf mmg-* build
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
.PHONY:$(SRCDIR)/$(INSTALL)
