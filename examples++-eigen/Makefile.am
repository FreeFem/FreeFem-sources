# $Id$

TESTS=../regtests.pl
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src FF_PROG_EXT=$(FF_PROG_EXT) \
	EXEEXT=$(EXEEXT)
EXTRA_DIST=*.edp all.edp

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v ^all.edp$$` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@
