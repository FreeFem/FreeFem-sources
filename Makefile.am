# Makefile for FreeFem++, adapted to Automake
# -------------------------------------------

# Adaptation to Automake: Antoine Le Hyaric - LJLL Paris 6 -
# lehyaric@ann.jussieu.fr - 13/5/04

# $Id$

SUBDIRS=download src DOC examples++-tutorial examples++	\
	examples++-eigen examples++-load examples++-mpi	\
	examples++-bug examples++-chapt3 examples++-other

EXTRA_DIST=regtests.sh config-wrapper.h arpack/arpack++/include/*.h \
FreeFem++.*mcp HISTORY HISTORY_BEFORE_2005  BUGS TODO regtests.m4 \
WindowsPackage.m4 README README_ARPACK  README_CVS README_CW  README_WINDOWS \
logo.ico copysharedlibs.sh COPYRIGHT edp.nedit INNOVATION reconfigure \
mode-mi-edp.zip aclocal.m4              acmacros.m4             acoptim.m4    \
examples-bamg/square/*_g.* examples-bamg/square/do*  \
examples-bamg/NACA012/[adp]* examples-bamg/NACA012/naca.awk \
 examples-bamg/quadloop/dotest.pl examples-bamg/test/dotest*.pl \
freefem++.spec

FF_MAC_PREFIX=FreeFem++v$(VERSION)$(ADD_PACKAGE_NAME)

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
	cd DOC && make documentation

DOC/freefem++doc.pdf:
	cd doc && make freefem++doc.pdf

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
	cd src/lglib && make
	cd src/fflib && make
	cd src/std && make

agl:
	cd src/lglib && make
	cd src/fflib && make
	cd src/agl&& make

nw:
	cd src/lglib && make
	cd src/fflib && make
	cd src/nw&& make

bamg:
	cd src/lglib && make
	cd src/fflib && make
	cd src/bamg && make

ide:
	cd download && make
	cd src/lglib && make
	cd src/fflib && make
	cd src/ide && make FreeFem++-cs$(EXEEXT)

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
	make check VISUALCHECK=yes

speedtest: all
	cd examples++-other && make check

# Windows package
# ---------------

# Windows package script (for Inno Setup). We extract version
# information from the Debian Changelog to get the package release
# number as well.
win32:WindowsPackage.iss
	cd examples++-load;tar zxvf include.tar.gz;cp -L include/* include-tmp
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

native:$(FF_MAC_PREFIX)_MacOS $(FF_MAC_PREFIX)_MacOsX.tgz DOC/freefem++doc.pdf
	echo "done"
native9:$(FF_MAC_PREFIX)_MacOS 
	echo "done"
nativeX: $(FF_MAC_PREFIX)_MacOsX.tgz 
	echo "done"

clean-local::
	-rm freefem++-$(VERSION).tar.gz
	-rm -rf FreeFem++v*_MacOS

# Native MacOS packaging
# ----------------------
ListFiles-natives: .FORCE 
	echo COPYRIGHT HISTORY HISTORY_BEFORE_2005 README BUGS TODO INSTALL INSTALL-MacOSX INNOVATION >$@	
	echo mode-mi-edp.zip >> $@
	echo examples++-tutorial/aile.msh examples++-tutorial/xyf >> $@
	echo examples++-load/load.link >> $@
	echo ./examples++-mpi/regtests.sh ./examples++-other/speedtest.sh >> $@
	echo ./download/fftw/Makefile ./download/fftw/Makefile.am  >>$@
	find .  -name '*.edp' -o -name '*.h*' -o -name '*.cpp'  |egrep  '[.]/examples++' >>$@
List-agl-dylib: src/agl/FreeFem++-agl
	otool -L src/agl/FreeFem++-agl|egrep -v '/System/Library/|/usr/lib/'|awk '/.dylib/ {print $$1}' >$@
.FORCE:

$(FF_MAC_PREFIX)_Macos:documentation ListFiles-natives
	-mkdir $@ 
	cat ListFiles-natives|xargs tar cf - | (cd $@ ; tar xf - )  
	/Developer/Tools/CpMac "FreeFem++(Carbon)" $@/FreeFem++ 
	cp   DOC/freefem++-doc.pdf $@

$(FF_MAC_PREFIX)_MacOsX: all documentation  ListFiles-natives List-agl-dylib
	-mkdir $@ 
	cp   DOC/freefem++doc.pdf $@
	cat ListFiles-natives|xargs tar cf - | (cd $@ ; tar xf - )  
	cd  $@ ; tar zxf ../src/agl/FreeFem++.app.tgz
	sed <src/agl/Info-plist.am >$@/FreeFem++.app/Contents/Info.plist \
           -e "s/@VVERSION@/$(VERSION)$(ADD_PACKAGE_NAME)/g" \
           -e "s/@DATE@/`date`/g"  
	-mkdir -p  $@/FreeFem++.app/Contents/FreeFem++.app/Contents
	cp $@/FreeFem++.app/Contents/Info.plist $@/FreeFem++.app/Contents/FreeFem++.app/Contents/Info.plist
	cp src/agl/FreeFem++-agl $@/FreeFem++.app/Contents/FreeFem++.app/FreeFem++
	cp src/FreeFem++-CoCoa $@
	-if [ -s List-agl-dylib ]; then tar zchvf $@/OtherMacOsLib.tgz `cat List-agl-dylib`; fi; 
	cp -p Install-MacOS.command $@

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
	cat ListFiles-natives|xargs tar cf - | (cd $@ ; tar xf - )  
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
