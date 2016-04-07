# $Id$  

all-local: all.edp regtests.edp  freefem++.pref

TESTS=3d-Leman.edp ArrayFE-3d.edp EqPoisson.edp Lac.edp Laplace-Adapt-3d.edp Laplace-Adapt-aniso-3d.edp Laplace3d.edp LaplaceRT-3d.edp NSI3d-carac.edp NSI3d.edp Period-Poisson-cube-ballon.edp Poisson-cube-ballon.edp Poisson.edp Poisson3d.edp Stokes.edp TruncLac.edp  beam-3d.edp cone.edp convect-3d.edp cube-period.edp cylinder-3d.edp cylinder.edp fallingspheres.edp first.edp meditddm.edp p.edp periodic-3d.edp pyramide.edp  refinesphere.edp  schwarz-nm-3d.edp sphere2.edp sphere6.edp tetgencube.edp tetgenholeregion.edp intlevelset3d.edp
XFAIL_TESTS=Laplace-Adapt-aniso-3d.edp fallingspheres.edp Laplace-Adapt-3d.edp Period-Poisson-cube-ballon.edp Poisson-cube-ballon.edp cylinder.edp refinesphere.edp tetgencube.edp tetgenholeregion.edp Poisson3d.edp  schwarz-nm-3d.edp 
LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw

LIST_IDP=MeshSurface.idp cube.idp

EXTRA_DIST=*.edp  all.edp regtests.edp regtests.m4 ref.edp dodecaedre01.mesh \
  lac-leman-v4.msh $(LIST_IDP)



all.edp:  Makefile
	(echo "NoGraphicWindow=true;NoUseOfWait=true;int verbosityy=verbosity;int MEM1234=storageused();"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do  \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;searchMethod=0;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << " mem leak = " <<storageused() - MEM1234 << endl << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# To create a new set of reference values in "ref.edp"
Ref: makeref.edp freefem++.pref
	export  FF_LOADPATH=../examples++-load/.;../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../regtests.m4
	m4 regtests.m4 > makeref.edp

freefem++.pref:
	echo loadpath = \"../examples++-load/\" >freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref
3d-Leman.edp:freefem++.pref
install-exec-local:: 
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/idp
	$(INSTALL)  -m 555  $(LIST_IDP)  $(DESTDIR)$(ff_prefix_dir)/idp



# To check the scripts against their reference values
regtests.edp: regtests.m4 ../regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp
FORCE:
