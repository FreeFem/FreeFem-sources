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
 mat_dervieux.cpp convect_dervieux.edp funcTemplate.cpp AFunction_ext.hpp funcTemplate.edp \
 NSP2BRP0.edp bilapMorley.edp convect_dervieux.edp convectchacon.edp \
 dfft.edp  funcTemplate.edp load.edp testFE.edp testFEMorley.edp \
  regtests.edp    regtests.m4  makeref.edp ref.edp myType.edp addNewType.cpp

LIST_COMPILE=myfunction.$(DYLIB_SUFFIX)  BernadiRaugel.$(DYLIB_SUFFIX) \
         Morley.$(DYLIB_SUFFIX) funcTemplate.$(DYLIB_SUFFIX) addNewType.$(DYLIB_SUFFIX)

load_compile:include  $(LIST_COMPILE)

.cpp.$(DYLIB_SUFFIX): include
	./load.link $< 

include: include.tar.gz
	tar zxvf include.tar.gz
	-rm include/._*
Ref: makeref.edp
	../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../regtests.m4
	m4 regtests.m4 > makeref.edp

# To check the scripts against their reference values
regtests.edp: regtests.m4 ../regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp

clean-local:
	-rm myfunction.o myfunction.so myfunction.dll
