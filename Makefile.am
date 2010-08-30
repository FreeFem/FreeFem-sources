# Makefile for FreeFem++, adapted to Automake
# -------------------------------------------

# Adaptation to Automake: Antoine Le Hyaric - LJLL Paris 6 -
# lehyaric@ann.jussieu.fr - 13/5/04

# $Id$

SUBDIRS=download src DOC examples++-tutorial examples++	\
	examples++-eigen examples++-load examples++-mpi	\
	examples++-bug examples++-chapt3 examples++-other \
	examples++-3d

EXTRA_DIST=regtests.sh config-wrapper.h  \
FreeFem++.mcp HISTORY HISTORY_BEFORE_2005  BUGS TODO regtests.m4 \
WindowsPackage.m4 README README_ARPACK  README_CVS README_CW  README_WINDOWS \
logo.ico copysharedlibs.sh COPYRIGHT edp.nedit INNOVATION reconfigure \
mode-mi-edp.zip aclocal.m4  acmacros.m4   acoptim.m4    \
freefem++.spec crimson-freefem++.zip \
Install-MacOS.command.in \
examples-bamg/NACA012/naca.awk  examples-bamg/quadloop/dotest.pl \
examples-bamg/square/*_g.* examples-bamg/square/do* examples-bamg/NACA012/[adp]* \
examples-bamg/test/dotest*.pl 0ldUserReadMe.txt CheckAllEdp


FF_MAC_PREFIX=FreeFem++v$(VERSION)$(ADD_PACKAGE_NAME)

FF_EXAMPLES_FILES = COPYRIGHT HISTORY HISTORY_BEFORE_2005 README BUGS TODO INSTALL INSTALL-MacOSX INNOVATION \
mode-mi-edp.zip \
examples++-tutorial/aile.msh examples++-tutorial/xyf \
examples++-3d/dodecaedre01.mesh  \
examples++-3d/lac-leman-v4.msh \
examples++-load/load.link \
./examples++-mpi/regtests.sh ./examples++-other/speedtest.sh  \
./download/fftw/Makefile ./download/fftw/Makefile.am \
examples++*/*.[ie]dp  CheckAllEdp


# Creates a file named "ChangeLog" containing the chronology of all
# modifications to the source files. Needs "cvs2cl" to be installed.

changelog:
	cvs2cl

# History before 2005 is stored in the file "HISTORY_BEFORE_2005"
all-local:: @HISTORY@

# try and avoid running this under fakeroot (otherwise we may face
# problems trying to connect to CVS as pseudo-root). Debian packaging
# uses fakeroot.
history:
	if test "$$FAKED_MODE" = ""; then\
		cvs2cl --hide-filenames -l "-d>01/01/2005" --file HISTORY;\
		fi

documentation:
	cd DOC && $(MAKE) $(AM_MAKEFLAGS) documentation

DOC/freefem++doc.pdf:
	cd $(@D) && $(MAKE) $(AM_MAKEFLAGS) $(@F)

clean-local::
	-find . \( -name '*~' -or  -name ListOfUnAllocPtr.bin \) |xargs rm 
	-rm examples*/*.eps 

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
	cd download && $(MAKE) $(AM_MAKEFLAGS) 
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
	cd examples++-other && $(MAKE) $(AM_MAKEFLAGS) check

# Windows package
# ---------------

# Windows package script (for Inno Setup). We extract version
# information from the Debian Changelog to get the package release
# number as well.
win32:WindowsPackage.iss
	cd examples++-load;tar zxvf include.tar.gz;mkdir include-tmp;cp -Lr include/* include-tmp
WindowsPackage.iss: WindowsPackage.m4 configure.ac Makefile.am
	m4 -DVERSION=$(VERSION)$(ADD_PACKAGE_NAME) \
		WindowsPackage.m4 > WindowsPackage.iss

# Debian package
# --------------

deb:
	dpkg-buildpackage -rfakeroot
	@echo Now run CopyToServer.sh in debian subdirectory

# Build all versions
# ------------------


nativeX: $(FF_MAC_PREFIX)_MacOsX.tgz 
	echo "done"

clean-local::
	-rm freefem++-$(VERSION).tar.gz
	-rm -rf FreeFem++v*_MacOS

# Native MacOS packaging
# ----------------------
install-exec-local:: 
	$(mkinstalldirs) -m 755 $(DESTDIR)$(pkgdatadir)/$(VERSION)
	tar cvf - $(FF_EXAMPLES_FILES)| (cd $(DESTDIR)$(pkgdatadir)/$(VERSION); tar xvf -)

ListFiles-natives: .FORCE 
#	echo COPYRIGHT HISTORY HISTORY_BEFORE_2005 README BUGS TODO INSTALL INSTALL-MacOSX INNOVATION >$@	
#	echo mode-mi-edp.zip >> $@
	echo $(FF_EXAMPLES_FILES) >$@
#	echo examples++-tutorial/aile.msh examples++-tutorial/xyf >> $@
#	echo examples++-3d/dodecaedre01.mesh  >> $@
#	echo examples++-load/load.link >> $@
#	echo ./examples++-mpi/regtests.sh ./examples++-other/speedtest.sh >> $@
	echo ./download/fftw/Makefile ./download/fftw/Makefile.am  >>$@
	find .  -name '*.[ei]dp' -o -name '*.h*' -o -name '*.cpp'  |egrep  '[.]/examples++' >>$@
	find .   -name '*.h*' -o -name '*.cpp'  |egrep  '[.]/examples++' >>$@
List-agl-dylib: src/nw/FreeFem++
	otool -L src/nw/FreeFem++|egrep -v '/System/Library/|/usr/lib/'|awk '/.dylib/ {print $$1}' >$@
.FORCE:

$(FF_MAC_PREFIX)_Macos:documentation ListFiles-natives
	-mkdir $@ 
	cat ListFiles-natives|xargs tar chf - | (cd $@ ; tar xf - )  
	/Developer/Tools/CpMac "FreeFem++(Carbon)" $@/FreeFem++ 
	cp   DOC/freefem++-doc.pdf $@

$(FF_MAC_PREFIX)_MacOsX: all documentation  ListFiles-natives List-agl-dylib
	-mkdir $@ 
	cp   DOC/freefem++doc.pdf $@
	cat ListFiles-natives|xargs tar chf - | (cd $@ ; tar xf - )  
	cd  $@ ; tar zxf ../src/agl/FreeFem++.app.tgz
	sed <src/agl/Info-plist.am >$@/FreeFem++.app/Contents/Info.plist \
           -e "s/@VVERSION@/$(VERSION)$(ADD_PACKAGE_NAME)/g" \
           -e "s/@DATE@/`date`/g"  
	cp src/nw/FreeFem++ $@/FreeFem++.app/Contents/bin
	cp src/nw/ffglut $@/FreeFem++.app/Contents/bin
	cp src/medit/ffmedit $@/FreeFem++.app/Contents/bin
	cp examples++-load/ff-get-dep  $@/FreeFem++.app/Contents/bin
	cp examples++-load/ff-pkg-download  $@/FreeFem++.app/Contents/bin
	sed <examples++-load/ff-c++ > $@/FreeFem++.app/Contents/bin/ff-c++ 	-e 's;FFAPPLI_INC;$@/FreeFem++.app/Contents/include;' 
	chmod a+x $@/FreeFem++.app/Contents/bin/ff-c++
	-mkdir $@/FreeFem++.app/Contents/include 
	-mkdir $@/FreeFem++.app/Contents/idp 
	cp examples++-load/include/* $@/FreeFem++.app/Contents/include
	cp examples++-load/*.dylib  $@/FreeFem++.app/Contents/lib
	cp examples++-*/*.idp  $@/FreeFem++.app/Contents/idp       
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
