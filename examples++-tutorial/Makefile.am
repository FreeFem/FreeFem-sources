# $Id$

all-local: all.edp regtests.edp
TESTS=../regtests.sh
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src REGEDP=regtests.edp	\
	EXEEXT=$(EXEEXT) X11PROG=$(X11PROG) GLXPROG=$(GLXPROG)	\
	AGLPROG=$(AGLPROG) MPIPROG=$(MPIPROG) VISUALCHECK=$(VISUALCHECK)
EXTRA_DIST=*.edp aile.msh xyf all.edp regtests.edp

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|refumf|refnoumf)\.edp$$'` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# To create a new set of reference values in "ref.edp"
Ref: makeref.edp
	../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4
	m4 -DUMFPACKLIB="@UMFPACKLIB@" regtests.m4 > makeref.edp

# To check the scripts against their reference values
regtests.edp: regtests.m4
	m4 -DUMFPACKLIB="@UMFPACKLIB@" -DASSERT regtests.m4 > regtests.edp
