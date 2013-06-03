# $Id$

all-local: all.edp regtests.edp  freefem++.pref
TESTS=adapt.edp adaptindicatorP1.edp adaptindicatorP2.edp algo.edp array.edp a_tutorial.edp beam.edp  calculus.edp cavity.edp convect2.edp convect-apt.edp convect.edp dumptable.edp ex-vf.edp FE.edp fluidStructAdapt.edp fluidStruct.edp freeboundary.edp freeboundary-weak.edp LapDG2.edp Laplace.edp LaplaceP1bis.edp LaplaceP1.edp LaplaceP1P2h.edp LaplaceRT.edp mesh.edp movemesh.edp nolinear-elas.edp NSUzawaCahouetChabart.edp onde.edp periodic4.edp Periodic.edp plot.edp readmesh.edp region.edp saverestore.edp schwarz-gc.edp schwarz-no-overlap.edp schwarz-overlap.edp sparse-matrix.edp sparse-cmatrix.edp StokesUzawa.edp tablefunction.edp
LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw


LIST_IDP=*.idp
EXTRA_DIST=*.edp *.idp aile.msh xyf all.edp regtests.edp regtests.m4 ref.edp

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# To create a new set of reference values in "ref.edp"
Ref: makeref.edp freefem++.pref
	../src/nw/FreeFem++-nw makeref.edp
install-exec-local:: 
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/idp
	$(INSTALL)  -m 555  $(LIST_IDP)  $(DESTDIR)$(ff_prefix_dir)/idp


makeref.edp: regtests.m4 ../regtests.m4
	m4 regtests.m4 > makeref.edp
freefem++.pref:
	echo loadpath = \"../examples++-load/\" >freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref

# To check the scripts against their reference values
regtests.edp: regtests.m4 ../regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp
FORCE: