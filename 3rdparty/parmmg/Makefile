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

ifeq ($(WIN32DLLTARGET),)
CMAKE_GUESS := 
else
CMAKE_GUESS := -DMPI_GUESS_LIBRARY_NAME=MSMPI -G "MSYS Makefiles"
endif

FAIT.done:tag-tar
	-mkdir build
	cd build && cmake ../parmmg-sources \
	-DCMAKE_C_COMPILER=$(CC) \
	-DCMAKE_C_FLAGS="$(CFLAGS) -fPIC" \
	-DCMAKE_CXX_COMPILER="$(CXX)" \
	-DCMAKE_CXX_FLAGS="$(CXXFLAGS) $(CXX11FLAGS) -fPIC" \
	-DCMAKE_Fortran_COMPILER="$(FC)" \
	-DMPI_C_COMPILER="$(MPICC)" \
	-DMPI_CXX_COMPILER="$(MPICXX)" \
	-DMPI_Fortran_COMPILER="$(MPIFC)" \
	-DMPIEXEC_EXECUTABLE="$(MPIRUN)" \
	-DDOWNLOAD_MMG=OFF \
	-DDOWNLOAD_METIS=OFF \
	-DMMG_DIR=../../mmg/mmg-sources \
	-DMMG_BUILDDIR=../../mmg/build \
	-DSCOTCH_DIR=$(scotch_dir) \
	-DMETIS_DIR=$(metis_dir) -DUSE_VTK=OFF \
	-DCMAKE_BUILD_TYPE=Release $(CMAKE_GUESS) && $(MAKE)
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
	echo parmmg LD -L@DIR@/lib -lparmmg >../lib/WHERE.parmmg ;
	echo parmmg INCLUDE -I@DIR@/include/>> ../lib/WHERE.parmmg ;
	echo build WHERE ./lib/WHERE.parmmg ;
	touch $@
clean::
	-rm WHERE.done

FAIRE: FAIT.done install.done

tag-tar: $(PACKAGE)
	-rm -rf ParMmg-* parmmg-sources
	unzip -q $(PACKAGE) && mv ParMmg-* parmmg-sources
	touch tag-tar

$(PACKAGE): FORCE
	../getall -o parmmg -a

FORCE:

clean::
	-rm ff-flags
	-rm FAIT* 
	-rm -rf ParMmg-* parmmg-sources build
	-rm flags-* tag-tar*
	-rm -rf ../include/parmmg
	-rm ../lib/libparmmg.a

ff-flags: ../Makefile Makefile
	grep 'abs_top_builddir *=' ../Makefile > ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'FC *=' ../Makefile >> ff-flags
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
	grep 'FF_prefix_petsc_real *=' ../Makefile >> ff-flags
	grep 'MPICC *=' ../Makefile >> ff-flags
	grep 'MPICXX *=' ../Makefile >> ff-flags
	grep 'MPIFC *=' ../Makefile >> ff-flags
	grep 'MPIRUN *=' ../Makefile >> ff-flags
	grep 'WIN32DLLTARGET *=' ../Makefile >> ff-flags
	-awk '/^scotch LD /&&/petsc/ {print "scotch_dir=$$(FF_prefix_petsc_real)" }' ../../plugin/seq/WHERE_LIBRARY-config >> ff-flags
	-awk '/^scotch LD /{print "scotch_dir=../../../3rdparty/" }' ../lib/WHERE.scotch >> ff-flags
	-awk '/^metis LD /&&/petsc/ {print "metis_dir=$$(FF_prefix_petsc_real)" }' ../../plugin/seq/WHERE_LIBRARY-config >> ff-flags
	-awk '/^metis LD /{print "metis_dir=../../../3rdparty/" }' ../lib/WHERE.metis >> ff-flags
	
.PHONY:$(SRCDIR)/$(INSTALL)
