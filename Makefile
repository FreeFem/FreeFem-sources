# -----  not change after this line ------------
TARGETS=FreeFem++ FreeFem++-x11 FreeFem++-mpi  FreeFem++-nw FreeFem++-glx
# --------------------
src?=./src
include $(src)/Makefile-$(HOSTTYPE)
CXXFLAGS =  $(OPTFLAGS)  $(FFFLAGS) $(includedir)  $(INCARPACKPP)
CXXMPIFLAGS=   $(CXXFLAGS) $(MPIFLAGS) -DPARALLELE $(INCARPACKPP)
LIBS=$(LIBUMFPACK) $(LIBARPACK) $(LIBF77)  $(LIBLOCAL)
VERSION=1.38
www=rascasse.inria.fr:/net/saupe/www_gamma/Gamma/cdrom/ftp/freefem/
ftp=rascasse.inria.fr:/ftp_gamma/freefem/
ftpk=baobab.ann.jussieu.fr:public_html/ftp/freefem/
OBJETS= $(EIGEN) AFunction.o   lex.o  lgfem.o  lgmesh.o \
 Drawing.o    \
 gibbs.o CheckPtr.o fem.o QuadratureFormular.o FESpace.o \
 Element_RT.o mshptg.o FQuadTree.o \
 QuadTree.o R2.o Meshio.o Mesh2.o Metric.o \
 BamgFreeFem.o MeshDraw.o MeshGeom.o MeshQuad.o SetOfE4.o \
 MeshRead.o MeshWrite.o problem.o strversionnumber.o InitFunct.o \
 lgalgo.o Element_P2h.o load.o


VPATH= $(src) $(src)/femlib $(src)/bamglib $(src)/Graphics $(src)/mpi $(src)/Algo $(src)/Eigen
includedir= -I$(INCLUDEX11)    -I$(src) -I$(src)/femlib -I$(src)/bamglib -I$(src)/Graphics -I$(src)/mpi -I$(src)/Algo $(IUMFPACK)

.SUFFIXES:.cpp .o

all: $(COMPILE_DIR) FORCE 
	$(MAKE) -C $(COMPILE_DIR) allc -f`pwd`/Makefile OPTFLAGS="$(OOPTFLAGS)" FFFLAGS="$(OFFFLAGS)" src=`pwd`/src
all-g: $(COMPILE_DIR)-g FORCE 
	$(MAKE) -C $(COMPILE_DIR)-g allc -f`pwd`/Makefile OPTFLAGS="$(GOPTFLAGS)" FFFLAGS="$(GFFFLAGS)" src=`pwd`/src

world: all all-g glx glx-g x11 x11-g $(BUILD_COCOA) $(BUILD_MPI)

$(COMPILE_DIR):
	-mkdir $@
$(COMPILE_DIR)-g:
	-mkdir $@

allc:  FreeFem++ FreeFem++-nw 


mpi: $(COMPILE_DIR) FORCE
	$(MAKE) -C $(COMPILE_DIR) FreeFem++-mpi -f`pwd`/Makefile OPTFLAGS="$(OOPTFLAGS)" FFFLAGS="$(OFFFLAGS)" src=`pwd`/src
mpi-g: $(COMPILE_DIR)-g FORCE
	$(MAKE) -C $(COMPILE_DIR)-g FreeFem++-mpi -f`pwd`/Makefile OPTFLAGS="$(GOPTFLAGS)" FFFLAGS="$(GFFFLAGS)" src=`pwd`/src

agl: $(COMPILE_DIR)-g FORCE
	$(MAKE) -C $(COMPILE_DIR) FreeFem++-agl -f`pwd`/Makefile OPTFLAGS="$(OOPTFLAGS)" FFFLAGS="$(OFFFLAGS)" src=`pwd`/src
agl-g:$(COMPILE_DIR)-g FORCE
	$(MAKE) -C $(COMPILE_DIR)-g FreeFem++-agl -f`pwd`/Makefile OPTFLAGS="$(GOPTFLAGS)" FFFLAGS="$(GFFFLAGS)" src=`pwd`/src
glx-g:$(COMPILE_DIR)-g FORCE
	$(MAKE) -C $(COMPILE_DIR)-g FreeFem++-glx -f`pwd`/Makefile OPTFLAGS="$(GOPTFLAGS)" FFFLAGS="$(GFFFLAGS)" src=`pwd`/src
glx:$(COMPILE_DIR)-g FORCE
	$(MAKE) -C $(COMPILE_DIR) FreeFem++-glx -f`pwd`/Makefile OPTFLAGS="$(OOPTFLAGS)" FFFLAGS="$(OFFFLAGS)" src=`pwd`/src
x11-g:$(COMPILE_DIR)-g FORCE
	$(MAKE) -C $(COMPILE_DIR)-g FreeFem++-x11 -f`pwd`/Makefile OPTFLAGS="$(GOPTFLAGS)" FFFLAGS="$(GFFFLAGS)" src=`pwd`/src
x11:$(COMPILE_DIR)-g FORCE
	$(MAKE) -C $(COMPILE_DIR) FreeFem++-x11 -f`pwd`/Makefile OPTFLAGS="$(OOPTFLAGS)" FFFLAGS="$(OFFFLAGS)" src=`pwd`/src


allex:examples++-tutorial/all.edp examples++/all.edp examples++-eigen/all.edp

yy:$(src)/lg.tab.hpp  $(src)/lg.tab.cpp $(src)/versionnumber.hpp
	echo ""
manuals:FORCE
	cd DOC; make all
FreeFem++:%: lg.tab.o   $(OBJETS) Xrgraph.o
	$(CXX) $(CXXFLAGS) $<  -o $@   $(OBJETS) Xrgraph.o $(LIBX11) $(LIBS)
FreeFem++-nw:%: lg.tab.o   $(OBJETS) sansrgraph.o
	$(CXX) $(CXXFLAGS) $<  -o $@   $(OBJETS)  sansrgraph.o $(LIBS)
FreeFem++-mpi: lg.tab.cpp  parallelempi.cpp $(OBJETS) sansrgraph.o  
	$(CXXMPI) $(CXXMPIFLAGS)   $(src)/lg.tab.cpp  $(src)/mpi/parallelempi.cpp  -o $@   $(OBJETS)  sansrgraph.o $(LIBS) $(MPILIBS)
FreeFem++-agl: lg.tab.cpp  macglrgraf.cpp  $(OBJETS) macglrgraf.o  
	$(CXX) $(CXXFLAGS)   $(src)/lg.tab.cpp  -o $@   $(OBJETS)   macglrgraf.o $(LIBGL) $(LIBS)
FreeFem++-glx: lg.tab.cpp  xglrgraf.cpp  $(OBJETS) xglrgraf.o  
	$(CXX) $(CXXFLAGS)   $(src)/lg.tab.cpp  -o $@   $(OBJETS)   xglrgraf.o $(LIBGLX) $(LIBS)

FreeFem++-x11:%: lg.tab.o   $(OBJETS) Xrgraph.o
	$(CXX) $(CXXFLAGS) $<  -o $@   $(OBJETS) Xrgraph.o $(LIBX11) $(LIBS)
%.o : %.c

clean: 
	-rm *.o  #* core
	-rm lg.*
	-find . \( -name '*~' -or  -name ListOfUnAllocPtr.bin \) |xargs rm 
	-rm examples*/*.eps 
$(src)/lg.tab.hpp  $(src)/lg.tab.cpp: lg.y
	-rm $(src)/lg.tab.hpp  $(src)/lg.tab.cpp lg.tab.cpp.h
	bison -dtv -p lg  $< -o lg.tab.cpp
	-mv lg.tab.cpp.h  lg.tab.hpp # pour un  pb entre des versions de bison 
	mv lg.tab.[hc]pp $(src)  
.cpp.o:
	$(CXX)  -c $(CXXFLAGS) $< 
$(src)/versionnumber.hpp: FORCE
	(echo "#define VersionFreeFempp " '$(VERSION)';echo "#define VersionFreeFemDate " '"'`date`'"') >$@
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


tyty: $(src)/. $(src)/*/. 
	egrep '^# *include' src/*/*.h src/*/*.hpp src/*.hpp src/*/*.cpp src/*.cpp | grep -v alloca.h| grep '"' \
	    | awk -F'[:"]' ' \
	/.hpp:#/ { nn=split($$1,bbb,"[/]");c=bbb[nn];d=$$(NF-1); l[c] = l[c] " " d  ;} \
	/.h:#/ { nn=split($$1,bbb,"[/]");c=bbb[nn];d=$$(NF-1); l[c] = l[c] " " d  ;} \
	/.cpp:#/  {nn=split($$1,bbb,"[/.]");c=bbb[nn-1];d=$$(NF-1);print c ".o:" d,l[d] }' 
    
$(src)/MakeDepend: $(src)/. $(src)/*/. 
	egrep '^# *include' src/*/*.h src/*/*.hpp src/*.hpp src/*/*.cpp src/*.cpp | grep -v alloca.h| grep '"' \
	    | awk -F'[:"]' ' \
	/.hpp:#/ { nn=split($$1,bbb,"[/]");c=bbb[nn];d=$$(NF-1); l[c] = l[c] " " d  ;} \
	/.h:#/ { nn=split($$1,bbb,"[/]");c=bbb[nn];d=$$(NF-1); l[c] = l[c] " " d  ;} \
	/.cpp:#/  {nn=split($$1,bbb,"[/.]");c=bbb[nn-1];d=$$(NF-1);print c ".o:" d,l[d] } \
	' >$@

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

include $(src)/MakeDepend


www:  FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit FreeFem++v$(VERSION)_MacOsX.tgz
	scp  HISTORY DOC/manual.pdf DOC/manual.ps.gz FreeFem++v$(VERSION).tar.gz FreeFem++v$(VERSION)_MacOsX.tgz FreeFem++v$(VERSION)_Win.zip  FreeFem++v$(VERSION)_MacOS.sit  $(ftpk)
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_Win.zip  freefem++.zip"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_MacOS.sit freefem++.sit"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION).tar.gz freefem++.tar.gz"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION).tar.gz freefem++.tgz"
	ssh  baobab.ann.jussieu.fr "cd www/ftp/freefem;ln -sf FreeFem++v$(VERSION)_MacOsX.tgz freefem++_MacOsX.tgz"
	scp  HISTORY  baobab.ann.jussieu.fr:/www/.
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
	   fi ; \
           
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

