# $Id$

all-local: myfunction.so
TESTS=../regtests.pl
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src REGEDP=load.edp	\
	EXEEXT=$(EXEEXT) X11PROG=$(X11PROG) GLXPROG=$(GLXPROG)	\
	AGLPROG=$(AGLPROG) MPIPROG=$(MPIPROG) VISUALCHECK=$(VISUALCHECK)
EXTRA_DIST=load.edp load.link myfunction.cpp

myfunction.so: myfunction.cpp
	./load.link myfunction
