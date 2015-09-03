# $Id$

all-local: all.edp regtests.edp freefem++.pref
TESTS=BlackScholes2D.edp NSNewton.edp NSprojection.edp condensor.edp convects.edp heatex.edp lame.edp membrane.edp membranerror.edp muwave.edp optimcontrol.edp potential.edp schwarz.edp sound.edp stokes.edp test1.edp testbed.edp thermal.edp thermic.edp
# if .. no arpack ... 
XFAIL_TESTS=sound.edp

EXTRA_DIST=*.edp  all.edp regtests.edp regtests.m4 ref.edp
LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# To create a new set of reference values in "ref.edp"
Ref: makeref.edp
	../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../regtests.m4
	m4 regtests.m4 > makeref.edp

# To check the scripts against their reference values
regtests.edp: regtests.m4 ../regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp
freefem++.pref:
	echo loadpath = \"../examples++-load/\" >freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref

FORCE: