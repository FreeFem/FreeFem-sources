# Makefile for FreeFem++, adapted to Automake
# -------------------------------------------

# Adaptation to Automake: Antoine Le Hyaric - LJLL Paris 6 -
# lehyaric@ann.jussieu.fr - 13/5/04

# $Id$

SUBDIRS=download src DOC examples++-tutorial examples++	\
	examples++-eigen examples++-load examples++-mpi	\
	examples++-bug examples++-other

EXTRA_DIST=regtests.sh config-wrapper.h arpack/arpack++/include/*.h		  \
FreeFem++.*mcp HISTORY README BUGS TODO regtests.m4				  \
FreeFem++.app.tgz 


# Creates a file named "ChangeLog" containing the chronology of all
# modifications to the source files. Needs "cvs2cl" to be installed.

changelog:
	cvs2cl -w -S --fsf

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
	m4 -DVERSION=$(VERSION)-$(FREEFEM_PACKAGING)$(OPTIM_TYPE) \
		WindowsPackage.m4 > WindowsPackage.iss

# Debian package
# --------------

deb:
	dpkg-buildpackage -rfakeroot
	@echo Now run CopyToServer.sh in debian subdirectory

# Build all versions
# ------------------

native:FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOS FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOsX.tgz
	echo "done"
native9:FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOS 
	echo "done"
nativeX: FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOsX.tgz
	echo "done"

clean-local::
	-rm freefem++-$(VERSION)-$(FREEFEM_PACKAGING).tar.gz
	-rm -rf FreeFem++v*_MacOS

# Native MacOS packaging
# ----------------------

FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOS:documentation	
	-mkdir $@  $@/examples++ $@/examples++-tutorial $@/examples++-bug $@/examples++-eigen
	/Developer/Tools/CpMac "FreeFem++(Carbon)" $@/FreeFem++ 
	cp  COPYRIGHT HISTORY README BUGS TODO  $@ 
	cp  examples++/*.edp   $@/examples++
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf  examples++-tutorial/*.edp   $@/examples++-tutorial
	-cp  examples++-bug/*.edp   $@/examples++-bug
	cp   DOC/manual-full.pdf $@

FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOsX: all documentation
	-mkdir $@  $@/examples++ $@/examples++-tutorial $@/examples++-bug $@/examples++-eigen $@/examples++-load
	cp COPYRIGHT HISTORY README BUGS TODO INSTALL-MacOSX  $@
	cp  examples++/*.edp   $@/examples++
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf  examples++-tutorial/*.edp   $@/examples++-tutorial
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	-cp  examples++-bug/*.edp   $@/examples++-bug
	cp   DOC/manual-full.pdf $@
	cd  FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOsX ; tar zxf ../FreeFem++.app.tgz
	sed <Info-plist.am >$@/FreeFem++.app/Contents/Info.plist \
           -e "s/@VVERSION@/$(VERSION)-$(FREEFEM_PACKAGING)/g" \
           -e "s/@DATE@/`date`/g"  
	-mkdir -p  $@/FreeFem++.app/Contents/FreeFem++.app/Contents
	cp $@/FreeFem++.app/Contents/Info.plist $@/FreeFem++.app/Contents/FreeFem++.app/Contents/Info.plist
	cp src/agl/FreeFem++-agl $@/FreeFem++.app/Contents/FreeFem++.app
	cp src/FreeFem++-CoCoa $@

FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOsX.tgz: FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOsX
	tar zcvf FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOsX.tgz  FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_MacOsX

# Linux binary-only package
# -------------------------

# Include kernel and libc version in static package name
PACKAGE_NAME=FreeFem++v$(VERSION)-$(FREEFEM_PACKAGING)_linux-$(KERNEL_VERSION)_$(LIBC_VERSION)$(OPTIM_TYPE)

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
	cp  COPYRIGHT HISTORY README BUGS TODO  $@ 
	cp  examples++/*.edp   $@/examples++
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf \
		examples++-tutorial/*.edp   $@/examples++-tutorial
	-cp  examples++-bug/*.edp   $@/examples++-bug

$(PACKAGE_NAME).tgz: $(PACKAGE_NAME)
	tar cvzf $@ $<

clean-local::
	-rm $(PACKAGE_NAME) $(PACKAGE_NAME).tgz
