# Makefile for FreeFem++, adapted to Automake
# -------------------------------------------

# Adaptation to Automake: Antoine Le Hyaric - LJLL Paris 6 -
# lehyaric@ann.jussieu.fr - 13/5/04

# $Id$

SUBDIRS=download src DOC examples++ examples++-tutorial		\
	examples++-eigen examples++-load examples++-mpi

EXTRA_DIST=regtests.pl arpack/arpack++/include/*.h src/*/*.h			  \
FreeFem++.*mcp.sit ffpc*mcp.zip HISTORY README BUGS TODO			  \
FreeFem++.app/Contents/Info.plist						  \
FreeFem++.app/Contents/MacOS/FreeFem++ FreeFem++.app/Contents/PkgInfo		  \
FreeFem++.app/Contents/Resources/Contents/Info.plist				  \
FreeFem++.app/Contents/Resources/Credits.rtf					  \
FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/classes.nib		  \
FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/info.nib		  \
FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/NSAppleMenuImage.tiff \
FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/objects.nib		  \
FreeFem++.app/Contents/Resources/script

# Creates a file named "ChangeLog" containing the chronology of all
# modifications to the source files. Needs "cvs2cl" to be installed.

changelog:
	cvs2cl -w -S --fsf

clean-local::
	-find . \( -name '*~' -or  -name ListOfUnAllocPtr.bin \) |xargs rm 
	-rm examples*/*.eps 

# Testing
# -------

# Standard automake goal "make check" is also valid. It just does not
# run any test that could prevent the user from working on its machine
# (because of unexpected windows opening right in the middle of the
# workspace).

visualtest: all
	src/std/FreeFem++ examples++-tutorial/all.edp

speedtest: all
	src/nw/FreeFem++-nw examples++-other/lap3-cpu.edp \
		|grep -E -- '-- lap (cholesky|CG|UMFPACK) +[0-9]+x[0-9]+'

# Build all versions
# ------------------

native:FreeFem++v$(VERSION)_MacOS FreeFem++v$(VERSION)_MacOsX.tgz
	echo "done"

clean-local::
	-rm freefem++-$(VERSION).tar.gz
	-rm -rf FreeFem++v*_MacOS

# Native MacOS packaging
# ----------------------

FreeFem++v$(VERSION)_MacOS:
	-mkdir $@  $@/examples++ $@/examples++-tutorial $@/examples++-bug $@/examples++-eigen
	/Developer/Tools/CpMac FreeFem++ $@ 
	cp  COPYRIGHT HISTORY README BUGS TODO  $@ 
	cp  examples++/*.edp   $@/examples++
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf  examples++-tutorial/*.edp   $@/examples++-tutorial
	-cp  examples++-bug/*.edp   $@/examples++-bug
	cp  DOC/manual.ps.gz DOC/manual.pdf $@

FreeFem++v$(VERSION)_MacOsX:
	-mkdir $@  $@/examples++ $@/examples++-tutorial $@/examples++-bug $@/examples++-eigen $@/examples++-load
	cp COPYRIGHT HISTORY README BUGS TODO INSTALL-MacOSX  $@
	cp  examples++/*.edp   $@/examples++
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf  examples++-tutorial/*.edp   $@/examples++-tutorial
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	-cp  examples++-bug/*.edp   $@/examples++-bug
	cp  DOC/manual.ps.gz DOC/manual.pdf $@
	cp -r FreeFem++.app FreeFem++v$(VERSION)_MacOsX/FreeFem++.app
	cp script/FreeFem++-CoCoa $@

FreeFem++v$(VERSION)_MacOsX.tgz: FreeFem++v$(VERSION)_MacOsX
	tar zcvf FreeFem++v$(VERSION)_MacOsX.tgz FreeFem++v$(VERSION)_MacOsX
