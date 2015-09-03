# $Id$
# F. hecht july 2014  add  LapEigenValueFunc.edp example 
all-local: all.edp freefem++.pref

LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw SKIP=$(SKIP_TESTS_EIGEN)

TESTS=BeamEigenValue.edp BeamEigenValueperio.edp Lap3dEigenValue.edp  LapEigenValue.edp LapnosymEigenValue.edp Stokes-eigen.edp VP-Steklov-Poincare.edp neuman.edp  WGM-sphere.edp free-cyl-axi.edp LapEigenValueFunc.edp
#  LapComplexEigenValue.edp  (buggus) 
EXTRA_DIST=*.edp all.edp

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

freefem++.pref:
	echo loadpath = \"../examples++-load/\" >freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref

