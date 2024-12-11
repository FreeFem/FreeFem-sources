############################################################################
# This file is part of FreeFEM.                                            #
#                                                                          #
# FreeFEM is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFEM is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          #
############################################################################
# SUMMARY : Makefile for FreeFem++, adapted to Automake
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : Antoine Le Hyaric, Simon Garnotel
# E-MAIL  :
# add F. Hecht 
bin_PROGRAMS=bin/md2edp

SUBDIRS=3rdparty src \
	plugin  \
	examples 

EXTRA_DIST=README.md CHANGELOG.md AUTHORS\
	bin/keys.awk \
	bin/mpic++.in \
	bin/Install-MacOS.command \
	bin/FreeFem++-CoCoa.in \
	bin/save-in-files \
	bin/cleanregen.sh \
	bin/copysharedlibs.sh \
	bin/script/clean_all_installed_ffpp.sh \
	bin/script/PostInstall.m4 \
	bin/regtests.sh \
	bin/uninstall-ff++ \
	bin/build/cleancrlf \
	bin/build/orgindex \
	bin/build/download \
	bin/build/ff-md5 \
	bin/build/links2files \
	bin/modpath.iss \
	bin/Install-MacOS.command.in \
	bin/Searchdylib \
	bin/test-driver-ff \
	etc/tools/md2edp.cpp  \
	etc/tools/srcHeader.txt \
	etc/tools/.uncrustify \
	etc/debian/freefem++-doc.doc-base \
	etc/debian/freefem++.manpages \
	etc/debian/apt-ftparchive.conf \
	etc/debian/freefem++.docs \
	etc/debian/compat \
	etc/debian/freefem++.examples \
	etc/debian/changelog \
	etc/debian/freefem++-nw.install \
	etc/debian/freefem++-nw.prerm \
	etc/debian/freefem++-mpich.links \
	etc/debian/freefem++-doc.docs \
	etc/debian/freefem++-mpich.dirs \
	etc/debian/FreeFem++.1 \
	etc/debian/freefem++-mpich.install \
	etc/debian/freefem++-mpich.postinst \
	etc/debian/freefem++-nw.postinst \
	etc/debian/shlibs.local \
	etc/debian/rules \
	etc/debian/freefem++-mpich.prerm \
	etc/debian/copyright \
	etc/debian/freefem++-nw.dirs \
	etc/debian/CopyToServer.sh \
	etc/debian/control \
	etc/debian/freefem++-nw.links \
	etc/config/FreeFem++.scrpt-txt.in \
	etc/config/FreeFem++.app.in.tgz \
	etc/config/m4/acoptim.m4 \
	etc/config/m4/acmpi.m4 \
	etc/config/m4/ax_lib_hdf5.m4 \
	etc/config/m4/ax_lib_gsl.m4 \
	etc/config/m4/regtests.m4 \
	etc/config/m4/acmacros.m4 \
	etc/config/m4/ax_openmp.m4 \
	etc/config/m4/tags.rOZFHw \
	etc/config/m4/WindowsPackage.m4 \
	etc/FreeFem++.app/Contents \
	etc/FreeFem++.app/Contents/Resources \
	etc/FreeFem++.app/Contents/Resources/English.lproj \
	etc/FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib \
	etc/FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/classes.nib \
	etc/FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/NSAppleMenuImage.tiff \
	etc/FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/info.nib \
	etc/FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/objects.nib \
	etc/FreeFem++.app/Contents/Resources/FreeFem++.icns \
	etc/FreeFem++.app/Contents/Resources/script \
	etc/FreeFem++.app/Contents/Resources/Credits.rtf \
	etc/FreeFem++.app/Contents/Info.plist \
	etc/FreeFem++.app/Contents/PkgInfo \
	etc/logo/logo.ico \
	etc/logo/logo.bmp \
	etc/doxygen/mainpage.dox \
	etc/doxygen/Doxyfile \
	etc/jenkins \
	etc/jenkins/install.sh \
	etc/jenkins/deployDEB.sh \
	etc/jenkins/build.sh \
	etc/jenkins/check.sh \
	etc/jenkins/configure.sh \
	etc/jenkins/clean.sh \
	etc/jenkins/build_PETSc.sh \
	etc/ff/ff.xcodeproj \
	etc/ff/ff.xcodeproj/project.pbxproj  \
	readme/BUGS readme/HISTORY_BEFORE_2005 readme/README_MAC.md readme/TODO \
	readme/COPYRIGHT readme/INNOVATION readme/README_COMPILATION.md	readme/README_WINDOWS.md \
	readme/HISTORY readme/README readme/README_GIT.md readme/README_XCODE.md \
idp/CC.idp \
idp/DDM-Schwarz-macro.idp \
idp/DDM-funcs-v2.idp \
idp/ElasticLaw2d.idp \
idp/ExtractDofsonBorder.idp \
idp/Heat3d.idp \
idp/MPIGMRESmacro.idp \
idp/MPIplot.idp \
idp/MeshSurface.idp \
idp/ball-buildlayer.idp \
idp/buildmeshS.idp \
idp/cobrameshcavity.idp \
idp/cube.idp \
idp/ffddm.idp \
idp/ffddm_coarsemeshCS.idp \
idp/ffddm_functions.idp \
idp/ffddm_geneoCS.idp \
idp/ffddm_geneoCS_3rdlvl.idp \
idp/ffddm_geneoCS_saddlepoint.idp \
idp/ffddm_parameters.idp \
idp/ffddm_partitioning.idp \
idp/func-max.idp \
idp/getARGV.idp \
idp/gsl.idp \
idp/macro_ddm.idp \
idp/macro_ddm_substructuring.idp \
idp/mortar-msh.idp \
idp/movemeshsmooth.idp 


FF_MAC_PREFIX=FreeFem++v$(VERSION)$(ADD_PACKAGE_NAME)

FF_EXAMPLES_FILES=examples/tutorial/aile.msh \
	examples/tutorial/xyf \
	examples/3d/dodecaedre01.mesh \
	examples/3d/lac-leman-v4.msh \
	plugin/seq/load.link \
	examples/plugin/cube.msh \
	examples/plugin/fig.pgm \
	examples/plugin/lg.pgm \
	examples/mpi/regtests.sh \
	examples/misc/speedtest.sh \
        examples/*/*.md \
	examples/*/*.edp \
	examples/CheckAllEdp \
	examples/plugin/ch.pts \
	examples/plugin/g.gmesh

clean-local::
	-find . \( -name '*~' -or  -name ListOfUnAllocPtr.bin \) | xargs rm
	-rm examples/*/*.eps

# "dist" targets
clean-local::
	-rm -f freefem++-*.tar.gz freefem++-*.zip
	-rm -f Output/FreeFem++-*.exe

# Reduced compilation
# -------------------

unittests::
	-cd unit && make

quick:
	cd src/libMesh && $(MAKE) $(AM_MAKEFLAGS)
	cd src/lglib && $(MAKE) $(AM_MAKEFLAGS)
	cd src/fflib && $(MAKE) $(AM_MAKEFLAGS)
	cd src/std && $(MAKE) $(AM_MAKEFLAGS)

nw:
	cd src/libMesh && $(MAKE) $(AM_MAKEFLAGS)
	cd src/lglib && $(MAKE) $(AM_MAKEFLAGS)
	cd src/fflib && $(MAKE) $(AM_MAKEFLAGS)
	cd src/nw && $(MAKE) $(AM_MAKEFLAGS)

bamg:
	cd src/libMesh && $(MAKE) $(AM_MAKEFLAGS)
	cd src/lglib && $(MAKE) $(AM_MAKEFLAGS)
	cd src/fflib && $(MAKE) $(AM_MAKEFLAGS)
	cd src/bamg && $(MAKE) $(AM_MAKEFLAGS)

ide:
	cd src/libMesh && $(MAKE) $(AM_MAKEFLAGS)
	cd 3rdparty && $(MAKE) $(AM_MAKEFLAGS)
	cd src/lglib && $(MAKE) $(AM_MAKEFLAGS)
	cd src/fflib && $(MAKE) $(AM_MAKEFLAGS)
	cd src/ide && $(MAKE) $(AM_MAKEFLAGS) FreeFem++-cs$(EXEEXT)

bin/md2edp$(EXEEXT):etc/tools/md2edp.cpp 
	$(CXX) $(CXXFLAGS) etc/tools/md2edp.cpp -o bin/md2edp$(EXEEXT)

# Cleaning generated files which are stored in the CVS repository, to
# avoid technical CVS conflicts.
clean-gen:
	./cleanregen.sh

# Testing
# -------

# The standard automake goal "make check" is also valid. It just does
# not run any test that could prevent the user from working on its
# machine (because of unexpected windows opening right in the middle
# of the workspace).

visualcheck: all
	$(MAKE) $(AM_MAKEFLAGS) check VISUALCHECK=yes

# Windows package
# ---------------

# Windows package script (for Inno Setup). We extract version
# information from the Debian Changelog to get the package release
# number as well.
win32: WindowsPackage.iss FreeFEM-documentation.pdf
	rm -rf plugin/include-tmp
	cd plugin; mkdir include-tmp; cp -Lr seq/include/* */*.h include-tmp
	cp -r 3rdparty/include/. plugin/include-tmp/.
	/c/Program\ Files\ \(x86\)/Inno\ Setup\ 6/ISCC.exe "`cygpath.exe WindowsPackage.iss`"
WindowsPackage.iss: etc/config/m4/WindowsPackage.m4 configure.ac Makefile.am
	sed -i 's;/mingw64/ff-petsc;/c/Program\\ Files/FreeFem++/ff-petsc;g' plugin/seq/WHERE_LIBRARY-config
	m4 -DVERSION='$(VERSION)$(ADD_PACKAGE_NAME)' \
	-DMPIPROG='$(MPIPROG)' -DSIZEOFPTR='$(SIZEOF_PTRINBIT)' \
	-DHOSTOS='$(host_os)' -DFFPETSC='$(FF_prefix_petsc_win32)' \
	 etc/config/m4/WindowsPackage.m4  > WindowsPackage.iss
	echo loadpath += '"!\."'  >freefem++.pref
	echo includepath += '"!\idp"'  >>freefem++.pref

# Debian package
# --------------

deb:
	dpkg-buildpackage -rfakeroot
	@echo Now run CopyToServer.sh in debian subdirectory

# Build all versions
# ------------------


nativeX: $(FF_MAC_PREFIX)_MacOsX.tgz
	echo "done"
MacOsX: FreeFem++.app.tgz
	-rm -rf  OsXxx
	mkdir -p  OsXxx/Applications/
	make install DESTDIR="`pwd`/OsXxx"
	tar zxvf FreeFem++.app.tgz -C OsXxx/Applications/
	cd OsXxx;tar cvfz  ../$(FF_MAC_PREFIX)_MacOsX.tgz .
	rm -rf OsXxx
clean-local::
	-rm freefem++-$(VERSION).tar.gz
	-rm -rf FreeFem++v*_MacOS

# Native MacOS packaging
# ----------------------
bin/movefortranlib:etc/config/movefortranlib.in ./config.status
	./config.status --file=$@:etc/config/movefortranlib.in
	chmod a+x $@
bin/init-var.sh:Makefile
	echo "# buid form Makefile" >$@
	echo "bindir=$(bindir)" >>$@
	echo "ff_prefix=$(ff_prefix)" >>$@
	echo "version=$(VERSION)" >>$@
	echo "ffsrc=$(abs_top_srcdir)" >>$@
	echo "ff_prefix_dir_lib=$(ff_prefix_dir_lib)" >>$@
	echo "ff_prefix_dir_lib_mpi=$(ff_prefix_dir_lib_mpi)" >>$@
	echo "ff_data=$(pkgdatadir)/$(VERSION)" >>$@
	chmod a+x $@
	
	
	
install-exec-local:: bin/script/PostInstall.sh FreeFEM-documentation.pdf bin/md2edp$(EXEEXT)

	test `uname` != Darwin || $(MAKE) FreeFem++-CoCoa bin/script/PostInstall.sh
	$(mkinstalldirs) -m 755 "$(DESTDIR)$(pkgdatadir)/$(VERSION)"
	test `uname` != Darwin || $(mkinstalldirs) -m 755 "$(DESTDIR)/etc/paths.d/"
	test `uname` != Darwin || $(mkinstalldirs) -m 755 "$(DESTDIR)/usr/local/bin"
	make build-edp 
	tar cvf - $(FF_EXAMPLES_FILES)| (cd "$(DESTDIR)$(pkgdatadir)/$(VERSION)"; tar xvf -)
	$(INSTALL_SCRIPT) examples/CheckAll examples/CheckAllEdp "$(DESTDIR)$(pkgdatadir)/$(VERSION)"
	$(mkinstalldirs) -m 755 "$(DESTDIR)${bindir}"
	$(mkinstalldirs) -m 755 "$(DESTDIR)$(ff_prefix_dir)/idp"
	cd idp; for i in *.idp; do \
		if [ -f $$i ]; then $(INSTALL) -m 555 $$i "$(DESTDIR)$(ff_prefix_dir)/idp"; fi; done
	test ! -s FreeFEM-documentation.pdf || $(INSTALL_DATA) FreeFEM-documentation.pdf "$(DESTDIR)$(pkgdatadir)" 	
	test `uname` != Darwin || ( rm "$(DESTDIR)/usr/local/bin/FreeFem++-CoCoa"; $(INSTALL_SCRIPT) FreeFem++-CoCoa "$(DESTDIR)/usr/local/bin" )		
	test `uname` != Darwin || $(INSTALL_SCRIPT) FreeFem++-CoCoa "$(DESTDIR)${bindir}"
#    need sudo ...	
#winDssSSs	test `uname` != Darwin || echo "$(bindir)" > "$(DESTDIR)/etc/paths.d/FreeFem++"
#	test `uname` != Darwin || ( $(MAKE) FreeFem++.app.tgz; test -d "$(DESTDIR)/Applications" || mkdir -p "$(DESTDIR)/Applications"; tar zxf FreeFem++.app.tgz -C "$(DESTDIR)/Applications" )
#	test `uname` != Darwin || ( $(MAKE) bin/movefortranlib )

bin/script/PostInstall.sh:./Makefile bin/script/PostInstall.m4
	m4 "-DFF__FVER=$(PACKAGE_VERSION)" "-DFF_BINDIR=$(bindir)" "-DFF__DATADIR=$(pkgdatadir)" bin/script/PostInstall.m4 > bin/script/PostInstall.sh
	chmod a+x bin/script/PostInstall.sh

FreeFEM-documentation.pdf:Makefile
if DOWNLOAD
	-rm FreeFEM-documentation.pdf.md5
	-$(WGET) https://github.com/FreeFem/FreeFem-doc/raw/pdf/FreeFEM-documentation.pdf.md5
	-md5sum -c FreeFEM-documentation.pdf.md5 || (rm FreeFEM-documentation.pdf ; $(WGET) https://github.com/FreeFem/FreeFem-doc/raw/pdf/FreeFEM-documentation.pdf)
else
	touch FreeFEM-documentation.pdf
	echo "no downloaded documentation"
endif

ListFiles-natives: .FORCE
	test `uname` != Darwin || $(MAKE)  FreeFem++-CoCoa
	echo $(FF_EXAMPLES_FILES) >$@
	echo ./3rdparty/fftw/Makefile ./3rdparty/fftw/Makefile.am  >>$@
	find .  -name '*.[ei]dp' -o -name '*.h*' -o -name '*.cpp' -o -name '*.pgm'  |egrep  '[.]/examples/examples' >>$@
	find .   -name '*.h*' -o -name '*.cpp'  |egrep  '[.]/examples/examples' >>$@
List-agl-dylib: src/nw/FreeFem++
	otool -L src/nw/FreeFem++|egrep -v '/System/Library/|/usr/lib/'|awk '/.dylib/ {print $$1}' >$@
CheckMacLib.sh: src/nw/FreeFem++
	echo "for i in `otool -L src/nw/FreeFem++|egrep -v '/System/Library/|/usr/lib/'|awk '/.dylib/ {printf($$1.OFS)  }'` ; do test ! -f $$i && exit 1; done; exit 0"  >$@
	chmod a+x $@
CheckMPIMacLib.sh: src/mpi/FreeFem++-mpi
	echo "for i in `otool -L src/mpi/FreeFem++-mpi|egrep -v '/System/Library/|/usr/lib/'|awk '/.dylib/ {printf($$1.OFS)  }'` ; do test ! -f $$i && exit 1; done; exit 0"  >$@
	chmod a+x $@

$(FF_MAC_PREFIX)_Macos:documentation ListFiles-natives
	-mkdir $@
	cat ListFiles-natives|xargs tar chf - | (cd $@ ; tar xf - )
	/Developer/Tools/CpMac "FreeFem++(Carbon)" $@/FreeFem++
if DOWNLOAD
	cp   FreeFEM-documentation.pdf $@
endif

$(FF_MAC_PREFIX)_MacOsX: all documentation  ListFiles-natives List-agl-dylib
	-mkdir $@
if DOWNLOAD
	cp   FreeFEM-documentation.pdf $@
endif
	cat ListFiles-natives|xargs tar chf - | (cd $@ ; tar xf - )
	cd  $@ ; tar zxf ../src/agl/FreeFem++.app.tgz
	sed <src/agl/Info-plist.am >$@/FreeFem++.app/Contents/Info.plist \
           -e "s/@VVERSION@/$(VERSION)$(ADD_PACKAGE_NAME)/g" \
           -e "s/@DATE@/`date`/g"
	cp src/nw/FreeFem++ $@/FreeFem++.app/Contents/bin
	cp src/nw/ffglut $@/FreeFem++.app/Contents/bin
	cp src/medit/ffmedit $@/FreeFem++.app/Contents/bin
	cp examples/plugin/ff-get-dep  $@/FreeFem++.app/Contents/bin
	cp examples/plugin/ff-pkg-download  $@/FreeFem++.app/Contents/bin
	sed <examples/plugin/ff-c++ > $@/FreeFem++.app/Contents/bin/ff-c++ 	-e 's;FFAPPLI_INC;$@/FreeFem++.app/Contents/include;'
	chmod a+x $@/FreeFem++.app/Contents/bin/ff-c++
	-mkdir $@/FreeFem++.app/Contents/include
	-mkdir $@/FreeFem++.app/Contents/idp
	cp examples/plugin/include/* $@/FreeFem++.app/Contents/include
	cp examples/plugin/*.dylib  $@/FreeFem++.app/Contents/lib
	cp idp/*.idp  $@/FreeFem++.app/Contents/idp
	-if [ -s List-agl-dylib ]; then tar zchvf $@/OtherMacOsLib.tgz `cat List-agl-dylib`; fi;
	./config.status  --file=$@/Install-MacOS.command:Install-MacOS.command.in
	chmod a+rx $@/Install-MacOS.command
	-mkdir $@/FreeFem++.app/Contents/etc
	echo loadpath += \"./\"  >$@/FreeFem++.app/Contents/etc/freefem++.pref
	echo loadpath += \"$(ff_prefix_dir)/lib\"  >>$@/FreeFem++.app/Contents/etc/freefem++.pref
	echo includepath += \"$(ff_prefix_dir)/edp\"  >>$@/FreeFem++.app/Contents/etc/freefem++.pref
	echo includepath += \"$(ff_prefix_dir)/idp\"  >>$@/FreeFem++.app/Contents/etc/freefem++.pref


$(FF_MAC_PREFIX)_MacOsX.tgz: $(FF_MAC_PREFIX)_MacOsX
	tar zcvf $(FF_MAC_PREFIX)_MacOsX.tgz  $(FF_MAC_PREFIX)_MacOsX

# Linux binary-only package
# -------------------------

# Include kernel and libc version in static package name
PACKAGE_NAME=FreeFem++v$(VERSION)_linux-$(KERNEL_VERSION)_$(LIBC_VERSION)$(OPTIM_TYPE)

linux-package: $(PACKAGE_NAME).tgz

# No direct dependency to "all" to be able to debug the packaging
# procedure on its own.

$(PACKAGE_NAME):  ListFiles-natives
	cat ListFiles-natives|xargs tar cfh - | (cd $@ ; tar xf - )
	-mkdir $@
	cp src/std/FreeFem++ $@
	./copysharedlibs.sh src/std/FreeFem++ $@
	cp src/nw/FreeFem++-nw $@
	./copysharedlibs.sh src/nw/FreeFem++-nw $@
	cp src/ide/FreeFem++-cs $@
	./copysharedlibs.sh src/ide/FreeFem++-cs $@
	cp src/ide/FreeFem++-server $@
	./copysharedlibs.sh src/ide/FreeFem++-server $@
	cp src/ide/FreeFem++-client $@
	./copysharedlibs.sh src/ide/FreeFem++-client $@
	-cp src/glx/FreeFem++-glx $@
	-./copysharedlibs.sh src/glx/FreeFem++-glx $@
	-cp src/mpi/FreeFem++-mpi $@
	-./copysharedlibs.sh src/mpi/FreeFem++-mpi $@


$(PACKAGE_NAME).tgz: $(PACKAGE_NAME)
	tar cvzf $@ $<

clean-local::
	-rm -r $(PACKAGE_NAME) $(PACKAGE_NAME).tgz
autofiles:AutoGeneratedFile.tar.gz
LIST_REQUIRE_AM= compile config.guess config.sub depcomp install-sh  missing  test-driver
LIST_AM_UNUSED=ar-lib  mdate-sh ylwrap py-compile texinfo.tex
LIST_GENERATE_FILE= configure config.h.in \
   Makefile.in \
    3rdparty/Makefile.in \
    3rdparty/blas/Makefile.in \
    3rdparty/arpack/Makefile.in \
    3rdparty/umfpack/Makefile.in \
    3rdparty/fftw/Makefile.in \
    src/Makefile.in \
    src/bamglib/Makefile.in \
    src/Graphics/Makefile.in \
    src/ffgraphics/Makefile.in \
    src/ffgraphics/client/Makefile.in \
    src/ffgraphics/server/Makefile.in \
    src/femlib/Makefile.in \
    src/Algo/Makefile.in \
    src/lglib/Makefile.in \
    src/fflib/Makefile.in \
    src/nw/Makefile.in \
    src/mpi/Makefile.in \
    src/bamg/Makefile.in \
    src/libMesh/Makefile.in \
    src/medit/Makefile.in \
    src/bin-win32/Makefile.in \
    plugin/Makefile.in \
	plugin/seq/Makefile.in \
	plugin/mpi/Makefile.in \
	examples/Makefile.in \
    examples/plugin/Makefile.in \
    examples/tutorial/Makefile.in \
    examples/misc/Makefile.in \
    examples/mpi/Makefile.in \
    examples/bem/Makefile.in \
    examples/hpddm/Makefile.in \
    examples/ffddm/Makefile.in \
    examples/eigen/Makefile.in \
    examples/examples/Makefile.in \
    examples/bug/Makefile.in \
    examples/3d/Makefile.in \
    examples/3dSurf/Makefile.in \
    examples/3dCurve/Makefile.in

LIST_GENERATE_FILE_AM= \
  Makefile.am \
    3rdparty/Makefile.am \
    3rdparty/blas/Makefile.am \
    3rdparty/arpack/Makefile.am \
    3rdparty/umfpack/Makefile.am \
    3rdparty/fftw/Makefile.am \
    src/Makefile.am \
    src/bamglib/Makefile.am \
    src/Graphics/Makefile.am \
    src/ffgraphics/Makefile.am \
    src/ffgraphics/client/Makefile.am \
    src/ffgraphics/server/Makefile.am \
    src/femlib/Makefile.am \
    src/Algo/Makefile.am \
    src/lglib/Makefile.am \
    src/fflib/Makefile.am \
    src/libMesh/Makefile.am \
    src/nw/Makefile.am \
    src/mpi/Makefile.am \
    src/bamg/Makefile.am \
    src/medit/Makefile.am \
    src/bin-win32/Makefile.am \
    plugin/Makefile.am \
	plugin/seq/Makefile.am \
	plugin/mpi/Makefile.am \
	examples/Makefile.am \
    examples/plugin/Makefile.am \
    examples/tutorial/Makefile.am \
    examples/misc/Makefile.am \
    examples/mpi/Makefile.am \
    examples/bem/Makefile.am \
    examples/hpddm/Makefile.am \
    examples/ffddm/Makefile.am \
    examples/eigen/Makefile.am \
    examples/examples/Makefile.am \
    examples/bug/Makefile.am \
    examples/3d/Makefile.am \
    examples/3dSurf/Makefile.am \
    examples/3dCurve/Makefile.am


#$(LIST_GENERATE_FILE):$(LIST_GENERATE_FILE_AM)
#	@echo "WARNING the configure file is older than  configure build flies"
#	@echo "Rebuild configure: do one the three case"
#	@echo " if you have autoconf    # make conf"
#	@echo " or without              #  make conf-without-autoconf"
#	@echo " or by pass this problem #  make conf-touch"
conf:
	autoreconf
	./config.status --recheck
conf-without-autoconf:
	tar zxvf AutoGeneratedFile.tar.gz
	./config.status --recheck
conf-touch:
	touch $(LIST_GENERATE_FILE)
AutoGeneratedFile.tar.gz:$(LIST_GENERATE_FILE)  configure.ac $(LIST_REQUIRE_AM)
	tar cvfz $@  $(LIST_GENERATE_FILE) $(LIST_REQUIRE_AM)

FreeFem++-CoCoa:bin/FreeFem++-CoCoa.in ./config.status
		./config.status  --file="FreeFem++-CoCoa:bin/FreeFem++-CoCoa.in"
FreeFem++.app.tgz:./config.status etc/config/FreeFem++.app.in.tgz etc/config/FreeFem++.scrpt-txt.in
	tar zxf etc/config/FreeFem++.app.in.tgz
	./config.status --file=FreeFem++.scrpt-txt:etc/config/FreeFem++.scrpt-txt.in
	-rm FreeFem++.app/Contents/Resources/Scripts/main.scpt
	 osacompile <FreeFem++.scrpt-txt -o FreeFem++.app/Contents/Resources/Scripts/main.scpt
	tar zcf FreeFem++.app.tgz FreeFem++.app

# ALH - build a robust TAGS file for a whole directory tree (more robust than the automake TAGS directive because of the
# multiple source trees downloaded within FF)

findtags:
	find . -type f -name '*.?pp'|grep -v lg.tab|xargs etags
clean-local::
	-rm TAGS

build-edp:bin/md2edp
	make bin/md2edp `echo  examples/*/*.md ' ' | sed 's/.md /.edp /g'`
##  New rule to build .edp file for md file
.md.edp: 
	bin/md2edp $<

##  trick to check if the configure configure.ac Makefile.in .. are correct in time
Makefile: @MAINTAINER_MODE_FALSE@ configure  $(LIST_GENERATE_FILE_AM)
	echo "configure Makefile:$(LIST_GENERATE_FILE_AM)"  >Makefile-for-Checkam
	echo '@false'| tr '@' '\t'  >> Makefile-for-Checkam
	if $(MAKE) -f Makefile-for-Checkam 2>&1 >/dev/null ; \
	 then ok=1 ; \
	else \
	$(MAKE) -f Makefile-for-Checkam ; \
	echo "Need to rebuid configure and Makefile files " ;\
	echo "Do autoreconf -i or tar zxvf tar zxvf AutoGeneratedFile.tar.gz" ;\
	echo " and ./reconfigure" ;\
	ok=0; \
	fi ;  \
	test $$ok -eq 1
.FORCE:
