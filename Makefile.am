# Makefile for FreeFem++, adapted to Automake
# -------------------------------------------

# Adaptation to Automake: Antoine Le Hyaric - LJLL Paris 6 -
# lehyaric@ann.jussieu.fr - 13/5/04

# $Id$

SUBDIRS=src DOC examples++

www=rascasse.inria.fr:www_gamma/Gamma/cdrom/ftp/freefem/
ftp=rascasse.inria.fr:/ftp_gamma/freefem/
ftpk=baobab.ann.jussieu.fr:public_html/ftp/freefem/

allex:examples++-tutorial/all.edp examples++/all.edp examples++-eigen/all.edp

clean-local::
	-find . \( -name '*~' -or  -name ListOfUnAllocPtr.bin \) |xargs rm 
	-rm examples*/*.eps 

versions:manuals  allex FreeFem++v$(VERSION).tar.gz FreeFem++v$(VERSION)_Win.zip FreeFem++v$(VERSION)_MacOS FreeFem++v$(VERSION)_MacOsX.tgz
	echo "done"
tgz: FreeFem++v$(VERSION).tar.gz 
	echo "done"
clean-version:
	-rm FreeFem++v$(VERSION).tar FreeFem++v$(VERSION).tar.gz 
	-rm -rf FreeFem++v*_MacOS
	-rm -rf FreeFem++v*_Win
FreeFem++v$(VERSION).tar:  FORCE
	mkdir FreeFem++v$(VERSION)
	(                                                             \
	echo COPYRIGHT HISTORY BUGS  Makefile* README README_CW README_ARPACK TODO *.mcp.*             ; \
	echo src/Makefile-* src/FreeFem++-CoCoa; \
	find src -type f -name '*pp' -o -name '*.[yhr]'             ; \
	echo DOC/plots/*.eps DOC/Makefile DOC/*sty DOC/*.tex  DOC/manual.ps.gz DOC/manual.pdf   ; \
	echo FreeFem++.*mcp.sit   ffpc*mcp.zip ;\
        find FreeFem++.app -type f|egrep -v '.DS_Store|Resources/FreeFem++'; \
	echo examples++/*.edp                                       ; \
	echo examples++-tutorial/aile.msh examples++-tutorial/xyf  examples++-tutorial/*.edp ; \
	echo examples++-eigen/*.edp                                   ; \
	echo examples++-bug/*.edp                                   ; \
	echo examples++-mpi/*.edp                                   ; \
	echo examples++-load/*.edp examples++-load/*pp examples++-load/*.link ; \
	echo arpack/arpack++/include ; \
	) | xargs tar cf - | (cd FreeFem++v$(VERSION); tar xf -)
	tar cvf $@ FreeFem++v$(VERSION)
	rm -rf FreeFem++v$(VERSION)
%.gz:%
	gzip -9f $*

FreeFem++v$(VERSION)_Win:pc/FreeFem++.exe FORCE
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
FreeFem++v$(VERSION)_MacOS: FORCE
	-mkdir $@  $@/examples++ $@/examples++-tutorial $@/examples++-bug $@/examples++-eigen
	/Developer/Tools/CpMac FreeFem++ $@ 
	cp  COPYRIGHT HISTORY README BUGS TODO  $@ 
	cp  examples++/*.edp   $@/examples++
	cp  examples++-eigen/*.edp   $@/examples++-eigen
	cp  examples++-tutorial/aile.msh examples++-tutorial/xyf  examples++-tutorial/*.edp   $@/examples++-tutorial
	-cp  examples++-bug/*.edp   $@/examples++-bug
	cp  DOC/manual.ps.gz DOC/manual.pdf $@



FreeFem++v$(VERSION)_MacOsX: FORCE
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
	tar zcvf FreeFem++v$(VERSION)_MacOsX.tgz FreeFem++v$(VERSION)_MacOsX;


%/all.edp:% FORCE 
	(cd $*; echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
         for i in *`ls *.edp|grep -v ^all.edp$$` ; do  \
	echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
	echo "verbosity=verbosityy;" ; \
        echo \{ include \"$$i\"\;\}\; ;\
	echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@
lex.o: lg.tab.hpp
lg.tab.o:lg.tab.hpp 

www:  FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit FreeFem++v$(VERSION)_MacOsX.tgz
	scp  HISTORY DOC/manual.pdf DOC/manual.ps.gz FreeFem++v$(VERSION).tar.gz FreeFem++v$(VERSION)_MacOsX.tgz FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit  $(ftpk)
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_Win.zip  freefem++.zip"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_MacOS.sit freefem++.sit"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION).tar.gz freefem++.tar.gz"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION).tar.gz freefem++.tgz"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_MacOsX.tgz freefem++_MacOsX.tgz"
	scp  HISTORY  baobab.ann.jussieu.fr:www/.
	ssh  baobab.ann.jussieu.fr "cd www/.; ./.update++ $(VERSION) <ff++.htmx >freefem++.htm"
	scp  HISTORY DOC/manual.pdf DOC/manual.ps.gz FreeFem++v$(VERSION).tar.gz FreeFem++v$(VERSION)_MacOsX.tgz FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit  $(ftp)
	scp HISTORY DOC/manual.pdf DOC/manual.ps.gz  $(www)
	scp FreeFem++v$(VERSION)_Win.zip $(www)/freefem++.zip
	scp FreeFem++v$(VERSION).tar.gz $(www)/freefem++.tar.gz
	scp FreeFem++v$(VERSION)_MacOS.sit $(www)/freefem++.sit
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION)_Win.zip  freefem++.zip"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION)_MacOS.sit freefem++.sit"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION).tar.gz freefem++.tar.gz"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION).tar.gz freefem++.tgz"
	ssh rascasse.inria.fr "cd /ftp_gamma/freefem;ln -sf FreeFem++v$(VERSION)_MacOsX.tgz freefem++_MacOsX.tgz"
	ssh rascasse.inria.fr "cd ~/public_html/; ./.update++ $(VERSION) <ff++.htmx >freefem++.htm"
Makefile-$(HOSTTYPE):
	echo " ------------------------------------ " ; \
	echo "Sorry the the $(src)/Makefile-$HOSTYPE do not exist" ;\
	echo "You can build  this file from existing  file :  $(src)/Makefile-macintosh" ;\
	echo " ------------------------------------ " 

install:   $(BIN_DIR)/.
	   @for i in $(TARGETS) ; do \
	     if [ -f $(COMPILE_DIR)/$$i ] ; then \
		cp  $(COMPILE_DIR)/$$i  $(BIN_DIR)/$$i; \
		chmod 755 $(BIN_DIR)/$$i; \
	     	echo " Install $(BIN_DIR)/$$i"  ; fi; \
	     if [ -f $(COMPILE_DIR)-g/$$i ] ; then \
		cp  $(COMPILE_DIR)-g/$$i $(BIN_DIR)/$$i-g ;\
		chmod 755 $(BIN_DIR)/$$i; \
		echo " Install $(BIN_DIR)/$$i-g"  ; fi; \
	   done; \
           if [ \( -f $(COMPILE_DIR)/FreeFem++-agl \) -a \( -n "$(BUILD_COCOA)"  \)  ] ; then \
		echo " Install $(BIN_DIR)/../FreeFem++.app/Contents/MacOS/FreeFem++"; \
		cp $(COMPILE_DIR)/FreeFem++-agl $(BIN_DIR)/../FreeFem++.app/Contents/MacOS/FreeFem++ ;\
		echo " Install FreeFem++.app/Resources/MacOS/FreeFem++-agl"; \
		cp $(COMPILE_DIR)/FreeFem++-agl ./FreeFem++.app/Contents/Resources/FreeFem++ ;\
                echo Install /Applications/FreeFem++.app; \
                (tar cf - FreeFem++.app| (cd /Applications;tar xf -)) ;\
                echo Install FreeFem++-CoCoa; \
                cp src/FreeFem++-CoCoa  $(BIN_DIR);\
	   fi ;
	echo "End of Installation "

get-arpack:
	-mkdir ../arpack;
	cd ../arpack; \
	wget  http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz ;\
	wget  http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz  ;\
	wget  http://www.caam.rice.edu/software/ARPACK/ARPACK++/arpack++.tar.gz ; \
	gunzip -c arpack96.tar.gz | tar xvf - ;\
	gunzip -c patch.tar.gz | tar xvf - ;\
	gunzip -c arpack++.tar.gz| tar xvf - 
FORCE:

