# $Id$

TESTS=../regtests.pl
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src			\
	EXEEXT=$(EXEEXT) MPIPROG=$(MPIPROG)
EXTRA_DIST=*.edp

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E ^all\.edp$$` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@
