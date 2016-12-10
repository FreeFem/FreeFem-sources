# Makefile for FreeFem++, adapted to Automake
# -------------------------------------------

# Adaptation to Automake: Antoine Le Hyaric - LJLL Paris 6 -
# lehyaric@ann.jussieu.fr - 13/5/04

SUBDIRS=download src  examples++-load  examples++-mpi \
	 examples++-tutorial examples++	\
	examples++-eigen	\
	examples++-bug examples++-chapt3 examples++-other \
	examples++-3d  examples++-hpddm  DOC 

EXTRA_DIST=test-driver-ff regtests.sh config-wrapper.h  \
FreeFem++.mcp HISTORY HISTORY_BEFORE_2005  BUGS TODO regtests.m4 \
INSTALL INSTALL-MacOSX INNOVATION \
WindowsPackage.m4 README README_ARPACK  README_WINDOWS README_MAC \
logo.ico copysharedlibs.sh COPYRIGHT edp.nedit  reconfigure \
mode-mi-edp.zip aclocal.m4  acmacros.m4   acoptim.m4 ax_lib_hdf5.m4  \
freefem++.spec crimson-freefem++.zip \
Install-MacOS.command.in \
examples-bamg/NACA012/naca.awk  examples-bamg/quadloop/dotest.pl \
examples-bamg/square/*_g.* examples-bamg/square/do* examples-bamg/NACA012/[adp]* \
examples-bamg/test/dotest*.pl 0ldUserReadMe.txt CheckAllEdp CheckAll  \
WHERE_LIBRARY-mkl FreeFem++-CoCoa.in  \
uninstall-ff++ \
./build/cleancrlf		./build/download		./build/links2files \
./build/orgindex ./build/ff-md5  script/PostInstall.m4 \
web/download.php	web/index.php		web/macosx.php		web/menu.php		web/phpfiles.php \
web/ffstyle.css		web/linux.php		web/main.php		web/news.php		web/windows.php 


FF_MAC_PREFIX=FreeFem++v$(VERSION)$(ADD_PACKAGE_NAME)

FF_EXAMPLES_FILES = COPYRIGHT HISTORY HISTORY_BEFORE_2005 README README_WINDOWS README_MAC BUGS TODO INSTALL INSTALL-MacOSX INNOVATION \
mode-mi-edp.zip \
examples++-tutorial/aile.msh examples++-tutorial/xyf \
examples++-3d/dodecaedre01.mesh  \
examples++-3d/lac-leman-v4.msh \
examples++-load/load.link \
examples++-load/cube.msh \
examples++-load/fig.pgm	examples++-load/lg.pgm \
./examples++-mpi/regtests.sh ./examples++-other/speedtest.sh  \
./download/fftw/Makefile ./download/fftw/Makefile.am \
examples++*/*.[ie]dp  CheckAllEdp \
examples++-load/ch.pts 


# Creates a file named "ChangeLog" containing the chronology of all
# modifications to the source files. Needs "cvs2cl" to be installed.

changelog:
	cvs2cl

# try and avoid running this under fakeroot (otherwise we may face
# problems trying to connect to CVS as pseudo-root). Debian packaging
# uses fakeroot.
history:
	if test "$$FAKED_MODE" = ""; then\
		hg log >HISTORY ;\
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
	test ! -d examples++-load/include-tmp  || rm -r examples++-load/include-tmp 
	cd examples++-load;tar zxvf include.tar.gz;mkdir include-tmp;cp -Lr *.h include/* include-tmp
WindowsPackage.iss: WindowsPackage.m4 configure.ac Makefile.am
	m4 -DVERSION='$(VERSION)$(ADD_PACKAGE_NAME)' \
	-DMPIPROG='$(MPIPROG)' -DSIZEOFPTR='$(SIZEOF_PTRINBIT)' \
	-DHOSTOS='$(host_os)'  WindowsPackage.m4 > WindowsPackage.iss

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
	test `uname` != Darwin || $(MAKE) FreeFem++-CoCoa  script/PostInstall.sh 
	$(mkinstalldirs) -m 755 $(DESTDIR)$(pkgdatadir)/$(VERSION)
	test `uname` != Darwin || $(mkinstalldirs) -m 755 $(DESTDIR)/etc/paths.d/
	test `uname` != Darwin || $(mkinstalldirs) -m 755 $(DESTDIR)/usr/local/bin
	tar cvf - $(FF_EXAMPLES_FILES)| (cd $(DESTDIR)$(pkgdatadir)/$(VERSION); tar xvf -)
	$(INSTALL_SCRIPT) CheckAll CheckAllEdp $(DESTDIR)$(pkgdatadir)/$(VERSION)
	$(mkinstalldirs) -m 755 $(DESTDIR)${bindir}
	test `uname` != Darwin || $(INSTALL_SCRIPT) FreeFem++-CoCoa  $(DESTDIR)${bindir}
	test `uname` != Darwin ||  echo $(bindir) >$(DESTDIR)/etc/paths.d/FreeFem++
	test `uname` != Darwin || $(INSTALL_SCRIPT) FreeFem++-CoCoa  $(DESTDIR)/usr/local/bin
	test `uname` != Darwin || ( $(MAKE) FreeFem++.app.tgz  ; test -d $(DESTDIR)/Applications ||  mkdir $(DESTDIR)/Applications ; tar zxf FreeFem++.app.tgz  -C $(DESTDIR)/Applications )

script/PostInstall.sh:./Makefile script/PostInstall.m4
	m4 script/PostInstall.m4 "-DFF__FVER=$(PACKAGE_VERSION)" "-DFF_BINDIR=$(bindir)" "-DFF__DATADIR=$(pkgdatadir)" >script/PostInstall.sh
	chmod a+x script/PostInstall.sh



ListFiles-natives: .FORCE 
	test `uname` != Darwin || $(MAKE)  FreeFem++-CoCoa
	echo $(FF_EXAMPLES_FILES) >$@
	echo ./download/fftw/Makefile ./download/fftw/Makefile.am  >>$@
	find .  -name '*.[ei]dp' -o -name '*.h*' -o -name '*.cpp' -o -name '*.pgm'  |egrep  '[.]/examples++' >>$@
	find .   -name '*.h*' -o -name '*.cpp'  |egrep  '[.]/examples++' >>$@
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
autofiles:AutoGeneratedFile.tar.gz

LIST_GENERATE_FILE= \
DOC/Makefile.in			download/umfpack/Makefile.in	examples++-tutorial/Makefile.in	src/femlib/Makefile.in \
Makefile.in			examples++-3d/Makefile.in	examples++/Makefile.in		src/fflib/Makefile.in \
config.h.in			examples++-bug/Makefile.in	src/Algo/Makefile.in		src/lglib/Makefile.in \
configure			examples++-chapt3/Makefile.in	src/Graphics/Makefile.in	src/libMesh/Makefile.in \
download/Makefile.in		examples++-eigen/Makefile.in	src/Makefile.in			src/medit/Makefile.in \
download/arpack/Makefile.in	examples++-load/Makefile.in	src/bamg/Makefile.in		src/mpi/Makefile.in \
download/blas/Makefile.in	examples++-mpi/Makefile.in	src/bamglib/Makefile.in		src/nw/Makefile.in \
download/fftw/Makefile.in	examples++-other/Makefile.in	src/bin-win32/Makefile.in     	examples++-hpddm/Makefile.in 

LIST_GENERATE_FILE_AM= \
DOC/Makefile.am			download/umfpack/Makefile.am	examples++-tutorial/Makefile.am	src/femlib/Makefile.am \
Makefile.am			examples++-3d/Makefile.am	examples++/Makefile.am		src/fflib/Makefile.am \
                                examples++-bug/Makefile.am	src/Algo/Makefile.am		src/lglib/Makefile.am \
configure.ac			examples++-chapt3/Makefile.am	src/Graphics/Makefile.am	src/libMesh/Makefile.am \
download/Makefile.am		examples++-eigen/Makefile.am	src/Makefile.am			src/medit/Makefile.am \
download/arpack/Makefile.am	examples++-load/Makefile.am	src/bamg/Makefile.am		src/mpi/Makefile.am \
download/blas/Makefile.am	examples++-mpi/Makefile.am	src/bamglib/Makefile.am		src/nw/Makefile.am \
download/fftw/Makefile.am	examples++-other/Makefile.am	src/bin-win32/Makefile.am	examples++-hpddm/Makefile.am

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
AutoGeneratedFile.tar.gz:$(LIST_GENERATE_FILE)  configure.ac
	tar cvfz $@  $(LIST_GENERATE_FILE)

FreeFem++-CoCoa:FreeFem++-CoCoa.in ./config.status
		./config.status  --file="FreeFem++-CoCoa:FreeFem++-CoCoa.in"
FreeFem++.app.tgz:./config.status
	tar zxf FreeFem++.app.tgz
	./config.status --file=FreeFem++.scrpt-txt:FreeFem++.scrpt-txt.in
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


