# $Id$

all-local: all.edp regtests.edp freefem++.pref

TESTS=aadaptation.edp aalapacien.edp aalaplace-nc.edp aamove.edp aaRT.edp arrayoFVh.edp bilap.edp D2.edp demo1.edp demo.edp funct.edp lapacienprecon.edp lap_mat.edp NSP1P1b.edp NSP1P1.edp NSP1P2.edp parareal.edp Richard.edp teste.edp testFE.edp wafer-heating-laser-axi.edp renumbering.edp
LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw

EXTRA_DIST=*.edp all.edp regtests.edp regtests.m4 ref.edp

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# To create a new set of reference values in "ref.edp"
clean-local:
	-rm *.fg *.eps *~ *.ps  Th?.faces  Th?.points Th?.BB Th?.mesh.gmsh   ListOf?llocPtr*.bin ffglut*.ppm 	toto.txt
Ref: makeref.edp
	../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../regtests.m4
	m4 regtests.m4 > makeref.edp
freefem++.pref:
	echo loadpath = \"../examples++-load/\" >freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref

# To check the scripts against their reference values
regtests.edp: regtests.m4 ../regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp
FORCE:
