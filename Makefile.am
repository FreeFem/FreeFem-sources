# Makefile for FreeFem++, adapted to Automake
# -------------------------------------------

# Adaptation to Automake: Antoine Le Hyaric - LJLL Paris 6 -
# lehyaric@ann.jussieu.fr - 13/5/04

# $Id$

SUBDIRS=download src DOC examples++-tutorial examples++	\
	examples++-eigen examples++-load examples++-mpi	\
	examples++-bug examples++-other

EXTRA_DIST=regtests.sh config-wrapper.h arpack/arpack++/include/*.h		  \
FreeFem++.*mcp HISTORY HISTORY_BEFORE_2005 README BUGS TODO regtests.m4



# Creates a file named "ChangeLog" containing the chronology of all
# modifications to the source files. Needs "cvs2cl" to be installed.

changelog:
	cvs2cl

# History before 2005 is stored in the file "HISTORY_BEFORE_2005"
all-local:: @HISTORY@
history:
	cvs2cl --hide-filenames -l "-d>01/01/2005" --file HISTORY

documentation:
	cd DOC && make documentation

clean-local::
	-find . \( -name '*~' -or  -name ListOfUnAllocPtr.bin \) |xargs rm 
	-rm examples*/*.eps 

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

ide:
	cd download && make
	cd src/lglib && make
	cd src/fflib && make
	cd src/ide && make FreeFem++-cs$(EXEEXT) FreeFem++-cs-server$(EXEEXT)

ide-client:
	cd download && make
	cd src/lglib && make
	cd src/fflib && make
	cd src/ide && make FreeFem++-cs$(EXEEXT)

# Cleaning generated files which are stored in the CVS repository, to
# avoid technical CVS conflicts.
clean-gen:
	-rm configure
	-find . -name Makefile.in -exec rm {} \;
	-rm examples++/regtests.edp

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

WindowsPackage.iss: WindowsPackage.m4 configure.ac Makefile.am
	m4 -DVERSION=$(VERSION)$(OPTIM_TYPE) \
		WindowsPackage.m4 > WindowsPackage.iss

# Debian package
# --------------

deb:
	dpkg-buildpackage -rfakeroot
	@echo Now run CopyToServer.sh in debian subdirectory

# Build all versions
# ------------------

native:FreeFem++v$(VERSION)_MacOS FreeFem++v$(VERSION)_MacOsX.tgz
	echo "done"
native9:FreeFem++v$(VERSION)_MacOS 
	echo "done"
nativeX: FreeFem++v$(VERSION)_MacOsX.tgz
	echo "done"

clean-local::
	-rm freefem++-$(VERSION).tar.gz
	-rm -rf FreeFem++v*_MacOS

# Native MacOS packaging
# ----------------------

FreeFem++v$(VERSION)_MacOS:documentation	
	-mkdir $@  $@/examples++ $@/examples++-tutorial $@/examples++-bug $@/examples++-eigen
	/Developer/Tools/CpMac "FreeFem++(Carbon)" $@/FreeFem++ 
	cp  COPYRIGHT HISTORY HISTORY_BEFORE_2005 README BUGS TODO  $@ 
	cp  examples++/*.edp   $@/examples++
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf  examples++-tutorial/*.edp   $@/examples++-tutorial
	-cp  examples++-bug/*.edp   $@/examples++-bug
	cp   DOC/manual-full.pdf $@

FreeFem++v$(VERSION)_MacOsX: all documentation
	-mkdir $@  $@/examples++ $@/examples++-tutorial $@/examples++-bug $@/examples++-eigen $@/examples++-load
	cp COPYRIGHT HISTORY HISTORY_BEFORE_2005 README BUGS TODO INSTALL-MacOSX  $@
	cp  examples++/*.edp   $@/examples++
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf  examples++-tutorial/*.edp   $@/examples++-tutorial
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	-cp  examples++-bug/*.edp   $@/examples++-bug
	cp   DOC/manual-full.pdf $@
	cd  FreeFem++v$(VERSION)_MacOsX ; tar zxf ../src/agl/FreeFem++.app.tgz
	sed <src/agl/Info-plist.am >$@/FreeFem++.app/Contents/Info.plist \
           -e "s/@VVERSION@/$(VERSION)/g" \
           -e "s/@DATE@/`date`/g"  
	-mkdir -p  $@/FreeFem++.app/Contents/FreeFem++.app/Contents
	cp $@/FreeFem++.app/Contents/Info.plist $@/FreeFem++.app/Contents/FreeFem++.app/Contents/Info.plist
	cp src/agl/FreeFem++-agl $@/FreeFem++.app/Contents/FreeFem++.app/FreeFem++
	cp src/FreeFem++-CoCoa $@

FreeFem++v$(VERSION)_MacOsX.tgz: FreeFem++v$(VERSION)_MacOsX
	tar zcvf FreeFem++v$(VERSION)_MacOsX.tgz  FreeFem++v$(VERSION)_MacOsX

# Linux binary-only package
# -------------------------

# Include kernel and libc version in static package name
PACKAGE_NAME=FreeFem++v$(VERSION)_linux-$(KERNEL_VERSION)_$(LIBC_VERSION)$(OPTIM_TYPE)

linux-package: $(PACKAGE_NAME).tgz
$(PACKAGE_NAME):all
	-mkdir $@  $@/examples++ $@/examples++-tutorial $@/examples++-bug \
		$@/examples++-eigen $@/examples++-load
	cp src/std/FreeFem++ $@
	./copysharedlibs.sh src/std/FreeFem++ $@
	cp src/nw/FreeFem++-nw $@
	./copysharedlibs.sh src/nw/FreeFem++-nw $@
	cp src/ide/FreeFem++-cs $@
	./copysharedlibs.sh src/ide/FreeFem++-cs $@
	cp src/ide/FreeFem++-cs-server $@
	./copysharedlibs.sh src/ide/FreeFem++-cs-server $@
	-cp src/glx/FreeFem++-glx $@
	-./copysharedlibs.sh src/glx/FreeFem++-glx $@
	-cp src/mpi/FreeFem++-mpi $@
	-./copysharedlibs.sh src/mpi/FreeFem++-mpi $@
	cp  COPYRIGHT HISTORY HISTORY_BEFORE_2005 README BUGS TODO  $@ 
	cp  examples++/*.edp   $@/examples++
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf \
		examples++-tutorial/*.edp   $@/examples++-tutorial
	-cp  examples++-bug/*.edp   $@/examples++-bug

$(PACKAGE_NAME).tgz: $(PACKAGE_NAME)
	tar cvzf $@ $<

clean-local::
	-rm $(PACKAGE_NAME) $(PACKAGE_NAME).tgz
