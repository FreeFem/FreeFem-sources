# Makefile for FreeFem++, adapted to Automake
# -------------------------------------------

# Adaptation to Automake: Antoine Le Hyaric - LJLL Paris 6 -
# lehyaric@ann.jussieu.fr - 13/5/04

# $Id$

SUBDIRS=src DOC examples++ examples++-tutorial examples++-eigen	\
examples++-load examples++-mpi examples++-bug

EXTRA_DIST=regtests.pl arpack/arpack++/include/*.h src/*/*.h			  \
FreeFem++.*mcp.sit ffpc*mcp.zip COPYRIGHT HISTORY README BUGS TODO		  \
FreeFem++.app/Contents/Info.plist						  \
FreeFem++.app/Contents/MacOS/FreeFem++ FreeFem++.app/Contents/PkgInfo		  \
FreeFem++.app/Contents/Resources/Contents/Info.plist				  \
FreeFem++.app/Contents/Resources/Credits.rtf					  \
FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/classes.nib		  \
FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/info.nib		  \
FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/NSAppleMenuImage.tiff \
FreeFem++.app/Contents/Resources/English.lproj/MainMenu.nib/objects.nib		  \
FreeFem++.app/Contents/Resources/script

clean-local::
	-find . \( -name '*~' -or  -name ListOfUnAllocPtr.bin \) |xargs rm 
	-rm examples*/*.eps 

# Build all versions
# ------------------

native:FreeFem++v$(VERSION)_Win.zip FreeFem++v$(VERSION)_MacOS FreeFem++v$(VERSION)_MacOsX.tgz
	echo "done"

clean-local::
	-rm freefem++-$(VERSION).tar.gz
	-rm -rf FreeFem++v*_MacOS
	-rm -rf FreeFem++v*_Win

# Native Windows version packaging
# --------------------------------

FreeFem++v$(VERSION)_Win:pc/FreeFem++.exe
	-mkdir $@  $@/examples++ $@/examples++-tutorial $@/examples++-bug $@/examples++-eigen
	cp  pc/FreeFem++.exe $@
	cp COPYRIGHT HISTORY README BUGS TODO  $@ 
	cp  examples++/*.edp   $@/examples++
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf examples++-tutorial/*.edp   $@/examples++-tutorial
	-cp  examples++-bug/*.edp   $@/examples++-bug
	cp  DOC/manual.ps.gz DOC/manual.pdf $@
zip:FreeFem++v$(VERSION)_Win.zip

FreeFem++v$(VERSION)_Win.zip:FreeFem++v$(VERSION)_Win
	-rm FreeFem++v$(VERSION)_Win.zip
	zip -9l FreeFem++v$(VERSION)_Win.zip `find FreeFem++v$(VERSION)_Win -type f | egrep -v 'exe$$|pdf$$|gz$$'`
	zip -9 FreeFem++v$(VERSION)_Win.zip `find FreeFem++v$(VERSION)_Win -type f | egrep  'exe$$|pdf$$|gz$$'`

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

# WWW access
# ----------

www=rascasse.inria.fr:www_gamma/Gamma/cdrom/ftp/freefem/
ftp=rascasse.inria.fr:/ftp_gamma/freefem/
ftpk=baobab.ann.jussieu.fr:public_html/ftp/freefem/

www:  FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit FreeFem++v$(VERSION)_MacOsX.tgz
	scp  HISTORY DOC/manual.pdf DOC/manual.ps.gz freefem++-$(VERSION).tar.gz FreeFem++v$(VERSION)_MacOsX.tgz FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit  $(ftpk)
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_Win.zip  freefem++.zip"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_MacOS.sit freefem++.sit"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf freefem++-$(VERSION).tar.gz freefem++.tar.gz"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf freefem++-$(VERSION).tar.gz freefem++.tgz"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_MacOsX.tgz freefem++_MacOsX.tgz"
	scp  HISTORY  baobab.ann.jussieu.fr:www/.
	ssh  baobab.ann.jussieu.fr "cd www/.; ./.update++ $(VERSION) <ff++.htmx >freefem++.htm"
	scp  HISTORY DOC/manual.pdf DOC/manual.ps.gz freefem++-$(VERSION).tar.gz FreeFem++v$(VERSION)_MacOsX.tgz FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit  $(ftp)
	scp HISTORY DOC/manual.pdf DOC/manual.ps.gz  $(www)
	scp FreeFem++v$(VERSION)_Win.zip $(www)/freefem++.zip
	scp freefem++-$(VERSION).tar.gz $(www)/freefem++.tar.gz
	scp FreeFem++v$(VERSION)_MacOS.sit $(www)/freefem++.sit
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION)_Win.zip  freefem++.zip"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION)_MacOS.sit freefem++.sit"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf freefem++-$(VERSION).tar.gz freefem++.tar.gz"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf freefem++-$(VERSION).tar.gz freefem++.tgz"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION)_MacOsX.tgz freefem++_MacOsX.tgz"
	ssh rascasse.inria.fr "cd ~/public_html/; ./.update++ $(VERSION) <ff++.htmx >freefem++.htm"

get-arpack:
	-mkdir ../arpack;
	cd ../arpack; \
	wget  http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz ;\
	wget  http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz  ;\
	wget  http://www.caam.rice.edu/software/ARPACK/ARPACK++/arpack++.tar.gz ; \
	gunzip -c arpack96.tar.gz | tar xvf - ;\
	gunzip -c patch.tar.gz | tar xvf - ;\
	gunzip -c arpack++.tar.gz| tar xvf - 
