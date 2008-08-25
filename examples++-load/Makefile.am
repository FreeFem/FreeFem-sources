# $Id$

all-local: @LOAD_COMPILE@
TESTS=@LOAD_TESTS@
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src REGEDP=regtests.edp		\
	EXEEXT=$(EXEEXT) X11PROG=$(X11PROG) GLXPROG=$(GLXPROG)		\
	AGLPROG=$(AGLPROG) MPIPROG=$(MPIPROG) IDEPROG="$(IDEPROG)"	\
	VISUALCHECK=$(VISUALCHECK)
EXTRA_DIST=load.link dfft.link  include.tar.gz \
 BernadiRaugel.cpp Morley.cpp dfft.cpp  \
 funcTemplate.cpp mat_dervieux.cpp mat_psi.cpp myfunction.cpp \
 mat_dervieux.cpp convect_dervieux.edp funcTemplate.cpp AFunction_ext.hpp \
 regtests.m4  makeref.edp ref.edp myType.edp addNewType.cpp \
  qf11to25.cpp  Element_P3.cpp  Element_P3.hpp  Element_P4.cpp  Element_P4.hpp \
  plot-fb-P3.edp  plot-fb-P4.edp  LaplaceP3.edp   LaplaceP4.edp \
  Element_P3dc.cpp  Element_P3dc.hpp Element_P4dc.cpp  Element_P4dc.hpp \
  Element_PkEdge.cpp \
  SuperLu.cpp SuperLu.edp \
  NewSolver.cpp \
  LapDG3.edp              convectchacon.edp       plot-fb-P3dc.edp \
  LapDG4.edp              dfft.edp                plot-fb-P4dc.edp \
  LaplaceP3.edp           ff-mt19937ar.edp        ref.edp \
  LaplaceP4.edp           funcTemplate.edp        regtests.edp \
  NSP2BRP0.edp            load.edp                testFE-P3.edp \
  all.edp                 makeref.edp             testFE-P3dc.edp \
  bilapMorley.edp         myType.edp              testFE.edp \
  convect_dervieux.edp    testFEMorley.edp \
  testFE-PkEdge.edp \
  splitmesh3.cpp  splitmesh3.edp  splitmesh6.cpp  splitmesh6.edp \
  SuperLU_3.0-include-ff.tar.gz REAME_SuperLU \
  glumesh2D.cpp glumesh2D.edp \
  glumesh3D.cpp glumesh3D.edp buildlayer.cpp tetgen.cpp\
  LayerMesh.hpp TransfoMesh_v2.hpp LayerMesh.cpp  TransfoMesh_v2.cpp\
  ffrandom.cpp ffrandom.edp  \
  popen.cpp 


LIST_COMPILE=myfunction.$(DYLIB_SUFFIX)  BernadiRaugel.$(DYLIB_SUFFIX) \
         Morley.$(DYLIB_SUFFIX) funcTemplate.$(DYLIB_SUFFIX) addNewType.$(DYLIB_SUFFIX) \
         qf11to25.$(DYLIB_SUFFIX)  Element_P3.$(DYLIB_SUFFIX)    Element_P4.$(DYLIB_SUFFIX)  \
	 Element_P3dc.$(DYLIB_SUFFIX)    Element_P4dc.$(DYLIB_SUFFIX) 	\
	 Element_PkEdge.$(DYLIB_SUFFIX)   \
	glumesh2D.$(DYLIB_SUFFIX)  glumesh3D.$(DYLIB_SUFFIX) buildlayer.$(DYLIB_SUFFIX) tetgen.$(DYLIB_SUFFIX)\
	splitmesh3.$(DYLIB_SUFFIX)  splitmesh6.$(DYLIB_SUFFIX) \
	ffrandom.$(DYLIB_SUFFIX)  \
	popen.$(DYLIB_SUFFIX)  \
	$(DYLIB_OTHER_COMPILE)

load_compile:include  $(LIST_COMPILE)


.cpp.$(DYLIB_SUFFIX): 
	 ./load.link CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS) $(CPPFLAGS)"  $< 
glumesh2D.$(DYLIB_SUFFIX): TransfoMesh_v2.cpp  LayerMesh.cpp glumesh2D.cpp 
	 ./load.link  CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS) $(CPPFLAGS)" TransfoMesh_v2.cpp glumesh2D.cpp
glumesh3D.$(DYLIB_SUFFIX): TransfoMesh_v2.cpp  LayerMesh.cpp glumesh3D.cpp 
	 ./load.link  CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS) $(CPPFLAGS)" TransfoMesh_v2.cpp glumesh3D.cpp
buildlayer.$(DYLIB_SUFFIX): TransfoMesh_v2.cpp  LayerMesh.cpp buildlayer.cpp 
	 ./load.link  CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS) $(CPPFLAGS)" TransfoMesh_v2.cpp LayerMesh.cpp  buildlayer.cpp
tetgen.$(DYLIB_SUFFIX): TransfoMesh_v2.cpp LayerMesh.cpp $(TETGEN_LOADLINK)
	 ./load.link  CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS) $(CPPFLAGS)" -I../download/include -L../download/lib -ltet TransfoMesh_v2.cpp  LayerMesh.cpp tetgen.cpp

include: include.tar.gz
	rm -rf include 
	gunzip -c include.tar.gz| tar xvf -
	-rm -f include/._*
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

$(LIST_COMPILE):include Makefile

clean-local:
	-rm *.o  $(LIST_COMPILE)  regtests.edp makeref.edp
	-rm -rf include
ff-c++:load.link
	sed <load.link >ff-c++ -e 's;INC=";INC="-I$(ff_prefix_dir)/include;'
	chmod a+x ff-c++

install-exec-local:: load_compile
	mkdir -p $(ff_prefix_dir)/lib
	mkdir -p $(ff_prefix_dir)/include
	mkdir -p $(ff_prefix_dir)/etc
	mkdir -p $(ff_prefix_dir)/edp
	$(INSTALL_DATA) include/* $(ff_prefix_dir)/include
	$(INSTALL) -m 555 $(LIST_COMPILE) $(ff_prefix_dir)/lib
	echo loadpath += \"$(ff_prefix_dir)/lib\"  >$(ff_prefix_dir)/etc/freefem++.pref
	echo includepath += \"$(ff_prefix_dir)/edp\"  >>$(ff_prefix_dir)/etc/freefem++.pref
	$(INSTALL_SCRIPT) ff-c++  ${bindir}



