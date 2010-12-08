# $Id$

all-local: @LOAD_COMPILE@
TESTS=@LOAD_TESTS@
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src REGEDP=regtests.edp		\
	EXEEXT=$(EXEEXT) X11PROG=$(X11PROG) GLXPROG=$(GLXPROG)		\
	AGLPROG=$(AGLPROG) MPIPROG=$(MPIPROG) IDEPROG="$(IDEPROG)"	\
	VISUALCHECK=$(VISUALCHECK)
EXTRA_DIST=ff-get-dep.in ff-pkg-download.in load.link.in  include.tar.gz \
 BernadiRaugel.cpp Morley.cpp dfft.cpp  \
 funcTemplate.cpp mat_dervieux.cpp mat_psi.cpp myfunction.cpp \
 mat_dervieux.cpp convect_dervieux.edp funcTemplate.cpp \
 regtests.m4  makeref.edp ref.edp myType.edp addNewType.cpp \
  qf11to25.cpp  Element_P3.cpp  Element_P3.hpp  Element_P4.cpp  Element_P4.hpp \
  plot-fb-P3.edp  plot-fb-P4.edp  LaplaceP3.edp   LaplaceP4.edp \
  Element_P3dc.cpp  Element_P3dc.hpp Element_P4dc.cpp  Element_P4dc.hpp \
  Element_PkEdge.cpp \
  SuperLu.cpp SuperLu.edp \
  NewSolver.cpp \
  LapDG3.edp              convectchacon.edp       plot-fb-P3dc.edp \
  LapDG4.edp              dfft.edp                plot-fb-P4dc.edp \
  LaplaceP3.edp                  ref.edp \
  LaplaceP4.edp           funcTemplate.edp        regtests.edp \
  NSP2BRP0.edp            load.edp                testFE-P3.edp \
  all.edp                 makeref.edp             testFE-P3dc.edp \
  bilapMorley.edp         myType.edp              testFE.edp \
  convect_dervieux.edp    testFEMorley.edp \
  testFE-PkEdge.edp \
  splitmesh3.cpp  splitmesh3.edp  splitmesh6.cpp  splitmesh6.edp \
  README_SuperLU \
  tetgen.cpp\
  msh3.cpp msh3.hpp\
  ffrandom.cpp ffrandom.edp  \
  medit.cpp  \
  bmo.cpp lgbmo.cpp bmo.hpp \
 fflapack.cpp	lapack.cpp	clapack.h lapack.edp \
 fig.pgm ppm2rnm.cpp ppmimg.h  ppm2rnm.edp \
 DxWriter.cpp  metis.cpp pcm.cpp pcm.hpp pcm2rnm.cpp \
 UMFPACK64.cpp iovtk.cpp gmsh.cpp MetricKuate.cpp  \
 ffnewuoa.cpp	ffnewuoa.edp	newuoa.f Element_P1dc1.cpp ilut.cpp ilut.edp \
 freeyams.cpp mmg3d.cpp  mshmet.cpp \
 BinaryIO.cpp	

LIST_COMPILE=myfunction.$(DYLIB_SUFFIX)  BernadiRaugel.$(DYLIB_SUFFIX) \
         Morley.$(DYLIB_SUFFIX) funcTemplate.$(DYLIB_SUFFIX) addNewType.$(DYLIB_SUFFIX) \
         qf11to25.$(DYLIB_SUFFIX)  Element_P3.$(DYLIB_SUFFIX)    Element_P4.$(DYLIB_SUFFIX)  \
	 Element_P3dc.$(DYLIB_SUFFIX)    Element_P4dc.$(DYLIB_SUFFIX) 	\
	 Element_PkEdge.$(DYLIB_SUFFIX)  msh3.$(DYLIB_SUFFIX) \
	splitmesh3.$(DYLIB_SUFFIX)  splitmesh6.$(DYLIB_SUFFIX) \
	ffrandom.$(DYLIB_SUFFIX)  \
	medit.$(DYLIB_SUFFIX)  \
	mat_dervieux.$(DYLIB_SUFFIX) lgbmo.$(DYLIB_SUFFIX) mat_psi.$(DYLIB_SUFFIX)\
	ppm2rnm.$(DYLIB_SUFFIX) DxWriter.$(DYLIB_SUFFIX) \
	pcm2rnm.$(DYLIB_SUFFIX) $(DYLIB_OTHER_COMPILE)  \
	iovtk.$(DYLIB_SUFFIX) 	 gmsh.$(DYLIB_SUFFIX) MetricKuate.$(DYLIB_SUFFIX) \
	Element_P1dc1.$(DYLIB_SUFFIX) BinaryIO.$(DYLIB_SUFFIX)

LIST_COMPILE_PKG= 	tetgen.$(DYLIB_SUFFIX) SuperLu.$(DYLIB_SUFFIX) dfft.$(DYLIB_SUFFIX) \
	metis.$(DYLIB_SUFFIX) UMFPACK64.$(DYLIB_SUFFIX) NewSolver.$(DYLIB_SUFFIX) \
	lapack.$(DYLIB_SUFFIX)  fflapack.$(DYLIB_SUFFIX) ffnewuoa.$(DYLIB_SUFFIX) ilut.$(DYLIB_SUFFIX) \
	freeyams.$(DYLIB_SUFFIX) mmg3d.$(DYLIB_SUFFIX) mshmet.$(DYLIB_SUFFIX)

bin_PROGRAMS=



load_compile:ff-c++ freefem++.pref WHERE_LIBRARY-download  include  $(LIST_COMPILE)  $(LIST_COMPILE_PKG)
	echo " finish build list $(DYLIB_SUFFIX)"



.cpp.$(DYLIB_SUFFIX): ff-c++
	 ./ff-c++ -auto  $< 

include: include.tar.gz
	rm -rf include 
	gunzip -c include.tar.gz| tar xvf -
	-rm -f include/._*
	touch include 
#	-cp -rf ../download/include/* include/.

Ref: makeref.edp
	../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../regtests.m4
	m4 regtests.m4 > makeref.edp

all-local: all.edp regtests.edp load_compile

all.edp: 
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# To check the scripts against their reference values
regtests.edp: regtests.m4 ../regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp

$(LIST_COMPILE):ff-c++ include Makefile 

clean-local:
	-rm *.o  load.link WHERE_LIBRARY-download ff-get-dep ff-c++ ff-pkg-download \
	  $(LIST_COMPILE) $(LIST_COMPILE_PKG)  regtests.edp makeref.edp
	-rm -rf include
ff-c++:load.link.in  Makefile load.link WHERE_LIBRARY-download ff-get-dep
	../config.status  --file=ff-c++:load.link.in
	chmod a+x ff-c++
load.link:load.link.in  Makefile
	../config.status  --file=load.link:load.link.in
	chmod a+x load.link
ff-pkg-download:ff-pkg-download.in Makefile
	../config.status  --file=$@:$@.in
	chmod a+x $@
	cp $@ ../download/bin
ff-get-dep:ff-get-dep.in Makefile
	../config.status  --file=$@:$@.in
	chmod a+x $@
	-if [ -d ../download/bin ] ;then cp $@ ../download/bin; fi
WHERE_LIBRARY-download:ff-pkg-download FORCE  WHERE_LIBRARY-config WHERE_LIBRARY
	@./ff-pkg-download >$@-new
	@diff $@-new $@  || (cp $@-new $@ ;touch WHERE_LIBRARY)
FORCE: ;
freefem++.pref:
	echo loadpath = \"./\" >freefem++.pref
install-exec-local:: load_compile
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/lib
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/include
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/etc
	$(INSTALL_DATA)  clapack.h ppmimg.h include/* $(DESTDIR)$(ff_prefix_dir)/include
	$(INSTALL)  -m 555 $(LIST_COMPILE) $(DESTDIR)$(ff_prefix_dir)/lib
	echo loadpath += \"./\"  >$(DESTDIR)$(ff_prefix_dir)/etc/freefem++.pref
	echo loadpath += \"$(ff_prefix_dir)/lib\"  >>$(DESTDIR)$(ff_prefix_dir)/etc/freefem++.pref
	echo includepath += \"$(ff_prefix_dir)/idp\"  >>$(DESTDIR)$(ff_prefix_dir)/etc/freefem++.pref
	$(INSTALL_SCRIPT) ff-c++  $(DESTDIR)${bindir}
	$(INSTALL_SCRIPT) ff-pkg-download  $(DESTDIR)${bindir}
	$(INSTALL_SCRIPT) ff-get-dep  $(DESTDIR)${bindir}
	-for i in $(LIST_COMPILE_PKG); do \
	 if [ -f $$i ] ; then 	$(INSTALL)  -m 555 $$i $(DESTDIR)$(ff_prefix_dir)/lib; fi; done




