############################################################################
# This file is part of FreeFem++.                                          #
#                                                                          #
# FreeFem++ is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFem++ is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFem++. If not, see <http://www.gnu.org/licenses/>.        #
############################################################################
# SUMMARY : Makefile for FreeFem++, adapted to Automake
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : Antoine Le Hyaric, Simon Garnotel
# E-MAIL  :

SUBDIRS=3rdparty src \
	plugin/seq plugin/mpi \
	plugin/hpddm plugin/ffddm \
	examples/examples examples/misc \
	examples/tutorial examples/other \
	examples/plugin examples/mpi \
	examples/eigen examples/bug \
	examples/3d examples/hpddm examples/ffddm

EXTRA_DIST=README.md CHANGELOG.md

FF_MAC_PREFIX=FreeFem++v$(VERSION)$(ADD_PACKAGE_NAME)

FF_EXAMPLES_FILES=examples/tutorial/aile.msh \
	examples/tutorial/xyf \
	examples/3d/dodecaedre01.mesh \
	examples/3d/lac-leman-v4.msh \
	examples/load/load.link \
	examples/load/cube.msh \
	examples/load/fig.pgm	\
	examples/plugin/lg.pgm \
	examples/mpi/regtests.sh \
	examples/other/speedtest.sh \
	examples/*/*.edp \
	examples/CheckAllEdp \
	examples/plugin/ch.pts \
	examples/plugin/g.gmesh

clean-local::
	-find . \( -name '*~' -or  -name ListOfUnAllocPtr.bin \) | xargs rm
	-rm examples/*/*.eps

# "dist" targets
clean-local::
	-rm freefem++-*.tar.gz freefem++-*.zip
	-rm Output/FreeFem++-*.exe

# Reduced compilation
# -------------------

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

speedtest: all
	cd examples/other && $(MAKE) $(AM_MAKEFLAGS) check

# Windows package
# ---------------

# Windows package script (for Inno Setup). We extract version
# information from the Debian Changelog to get the package release
# number as well.
win32:WindowsPackage.iss
	test ! -d examples/plugin/include-tmp  || rm -r examples/plugin/include-tmp
	cd examples/plugin;tar zxvf include.tar.gz;mkdir include-tmp;cp -Lr *.h include/* include-tmp
	/c/Program\ Files\ \(x86\)/Inno\ Setup\ 5/ISCC.exe  "`cygpath.exe WindowsPackage.iss`"
WindowsPackage.iss: m4/WindowsPackage.m4 configure.ac Makefile.am
	m4 -DVERSION='$(VERSION)$(ADD_PACKAGE_NAME)' \
	-DMPIPROG='$(MPIPROG)' -DSIZEOFPTR='$(SIZEOF_PTRINBIT)' \
	-DHOSTOS='$(host_os)'  m4/WindowsPackage.m4 > WindowsPackage.iss
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
install-exec-local::
	test `uname` != Darwin || $(MAKE) FreeFem++-CoCoa  bin/script/PostInstall.sh
	$(mkinstalldirs) -m 755 $(DESTDIR)$(pkgdatadir)/$(VERSION)
	test `uname` != Darwin || $(mkinstalldirs) -m 755 $(DESTDIR)/etc/paths.d/
	test `uname` != Darwin || $(mkinstalldirs) -m 755 $(DESTDIR)/usr/local/bin
	tar cvf - $(FF_EXAMPLES_FILES)| (cd $(DESTDIR)$(pkgdatadir)/$(VERSION); tar xvf -)
	$(INSTALL_SCRIPT) examples/CheckAll examples/CheckAllEdp $(DESTDIR)$(pkgdatadir)/$(VERSION)
	$(mkinstalldirs) -m 755 $(DESTDIR)${bindir}
	test `uname` != Darwin || $(INSTALL_SCRIPT) FreeFem++-CoCoa  $(DESTDIR)${bindir}
	test `uname` != Darwin ||  echo $(bindir) >$(DESTDIR)/etc/paths.d/FreeFem++
	test `uname` != Darwin || $(INSTALL_SCRIPT) FreeFem++-CoCoa  $(DESTDIR)/usr/local/bin
	test `uname` != Darwin || ( $(MAKE) FreeFem++.app.tgz  ; test -d $(DESTDIR)/Applications ||  mkdir $(DESTDIR)/Applications ; tar zxf FreeFem++.app.tgz  -C $(DESTDIR)/Applications )

script/PostInstall.sh:./Makefile bin/script/PostInstall.m4
	m4 bin/script/PostInstall.m4 "-DFF__FVER=$(PACKAGE_VERSION)" "-DFF_BINDIR=$(bindir)" "-DFF__DATADIR=$(pkgdatadir)" > bin/script/PostInstall.sh
	chmod a+x bin/script/PostInstall.sh



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
.FORCE:

$(FF_MAC_PREFIX)_Macos:documentation ListFiles-natives
	-mkdir $@
	cat ListFiles-natives|xargs tar chf - | (cd $@ ; tar xf - )
	/Developer/Tools/CpMac "FreeFem++(Carbon)" $@/FreeFem++
	cp   freefem++-doc.pdf $@

$(FF_MAC_PREFIX)_MacOsX: all documentation  ListFiles-natives List-agl-dylib
	-mkdir $@
	cp   freefem++doc.pdf $@
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
	plugin/hpddm/Makefile.in \
	plugin/ffddm/Makefile.in \
	examples/Makefile.in \
    examples/plugin/Makefile.in \
    examples/tutorial/Makefile.in \
    examples/misc/Makefile.in \
    examples/mpi/Makefile.in \
    examples/hpddm/Makefile.in \
    examples/ffddm/Makefile.in \
    examples/eigen/Makefile.in \
    examples/examples/Makefile.in \
    examples/bug/Makefile.in \
    examples/other/Makefile.in \
    examples/3d/Makefile.in 


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
    src/femlib/Makefile.am \
    src/Algo/Makefile.am \
    src/lglib/Makefile.am \
    src/fflib/Makefile.am \
    src/nw/Makefile.am \
    src/mpi/Makefile.am \
    src/bamg/Makefile.am \
    src/libMesh/Makefile.am \
    src/medit/Makefile.am \
    src/bin-win32/Makefile.am \
    plugin/Makefile.am \
	plugin/seq/Makefile.am \
	plugin/mpi/Makefile.am \
	plugin/hpddm/Makefile.am \
	plugin/ffddm/Makefile.am \
	examples/Makefile.am \
    examples/plugin/Makefile.am \
    examples/tutorial/Makefile.am \
    examples/misc/Makefile.am \
    examples/mpi/Makefile.am \
    examples/hpddm/Makefile.am \
    examples/ffddm/Makefile.am \
    examples/eigen/Makefile.am \
    examples/examples/Makefile.am \
    examples/bug/Makefile.am \
    examples/other/Makefile.am \
    examples/3d/Makefile.am 


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
FreeFem++.app.tgz:./config.status etc/config/FreeFem++.app.in.tgz
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
