# $Id$

all-local: @LOAD_COMPILE@
TESTS=@LOAD_TESTS@
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src REGEDP=load.edp		\
	EXEEXT=$(EXEEXT) X11PROG=$(X11PROG) GLXPROG=$(GLXPROG)		\
	AGLPROG=$(AGLPROG) MPIPROG=$(MPIPROG) IDEPROG="$(IDEPROG)"	\
	VISUALCHECK=$(VISUALCHECK)
EXTRA_DIST=load.edp load.link myfunction.cpp  dfft.cpp dfft.edp dfft.link include.tar.gz

myfunction.so: myfunction.cpp ../src/fflib/AFunction.hpp ../config.h include
	CXX=$(CXX) ./load.link myfunction
include: include.tar.gz
	tar zxvf include.tar.gz