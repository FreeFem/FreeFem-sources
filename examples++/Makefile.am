# $Id$

all-local: all.edp regtests.edp
TESTS=../regtests.pl
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src REGEDP=regtests.edp	\
	EXEEXT=$(EXEEXT) X11PROG=$(X11PROG) GLXPROG=$(GLXPROG)	\
	AGLPROG=$(AGLPROG) MPIPROG=$(MPIPROG)
EXTRA_DIST=*.edp all.edp regtests.edp

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests)\.edp$$'` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# Special test script, without any troublesome window-opening scripts
regtests.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(aaa-adp|testadp|all|regtests)\.edp$$'` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@
