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

ifeq ($(WIN32DLLTARGET),)
CMAKE_GENERATOR := 
else
CMAKE_GENERATOR := -G "MSYS Makefiles"
endif

FAIT.done:tag-tar
	-mkdir build
	cd build && cmake ../mmg-sources \
	-DCMAKE_C_COMPILER=$(CC) \
	-DCMAKE_C_FLAGS="$(CFLAGS) -fPIC" \
	-DCMAKE_CXX_COMPILER=$(CXX) \
	-DCMAKE_CXX_FLAGS="$(CXXFLAGS) $(CXX11FLAGS) -fPIC" \
	-DM_LIB="-lm" -DUSE_VTK=OFF \
	-DSCOTCH_DIR=$(scotch_dir) \
	-DCMAKE_BUILD_TYPE=Release $(CMAKE_GENERATOR) && $(MAKE)
	touch FAIT.done

install.done:FAIT.done
	cp -r build/include/mmg ../include/mmg
	cp build/lib/libmmg.a ../lib
	touch $@
clean::
	-rm install.done FAIT.done

mmg:$(PACKAGE)

install:install.done WHERE.done

WHERE.done: install.done
	echo mmg LD -L@DIR@/lib -lmmg >../lib/WHERE.mmg ;
	echo mmg INCLUDE -I@DIR@/include/>> ../lib/WHERE.mmg ;
	echo build WHERE ./lib/WHERE.mmg ;
	touch $@
clean::
	-rm WHERE.done

FAIRE: FAIT.done install.done

tag-tar: $(PACKAGE)
	-rm -rf mmg-*
	unzip $(PACKAGE) && mv mmg-* mmg-sources
#	patch -p1 <patch-mmg 
	touch tag-tar

$(PACKAGE): FORCE
	../getall -o mmg -a

FORCE:

clean::
	-rm ff-flags
	-rm FAIT* 
	-rm -rf mmg-* build
	-rm flags-* tag-tar*
	-rm -rf ../include/mmg
	-rm ../lib/libmmg.a

ff-flags: ../Makefile Makefile
	grep 'abs_top_builddir *=' ../Makefile > ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'CXX *=' ../Makefile >> ff-flags
	grep 'CXXFLAGS *=' ../Makefile >> ff-flags
	grep 'CXX11FLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags
	grep 'STD_LIBS *=' ../Makefile >> ff-flags
	grep 'CNOFLAGS *=' ../Makefile >> ff-flags
	grep 'WIN32DLLTARGET *=' ../Makefile >> ff-flags
	grep 'FF_prefix_petsc_real *=' ../Makefile >> ff-flags
	-awk '/^scotch LD /&&/petsc/ {print "scotch_dir=$$(FF_prefix_petsc_real)" }' ../../plugin/seq/WHERE_LIBRARY-config >> ff-flags
	-awk '/^scotch LD /{print "scotch_dir=../../../3rdparty/" }' ../lib/WHERE.scotch >> ff-flags

.PHONY:$(SRCDIR)/$(INSTALL)
