# $Id$  

all-local:   freefem++.pref

TESTS=Helmholtz-2d-Cobracavity.edp		Maxwell-3d-simple.edp			diffusion-3d-intermediate.edp		diffusion-3d-simple.edp \
Helmholtz-2d-marmousi.edp		Maxwell_Cobracavity.edp			diffusion-3d-minimal-ddm.edp		diffusion-3d-test.edp \
Helmholtz-2d-simple.edp			Navier-2d-marmousi2.edp			diffusion-3d-minimal-direct.edp		elasticity-3d-simple.edp \
Helmholtz-3d-overthrust.edp		Navier-3d-overthrust.edp		diffusion-3d-obsolete.edp		elasticity-3d-thirdlevelgeneo.edp \
Helmholtz-3d-semiopencavity.edp		Richards-2d.edp				diffusion-3d-scalability.edp		natural_convection.edp \
Helmholtz-3d-simple.edp			diffusion-2d-thirdlevelgeneo.edp	diffusion-3d-sequential.edp		natural_convection_3D_obstacle.edp

LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" SKIP=$(SKIP_TESTS_MPI)

LIST_IDP=cobrameshcavity.idp			ffddm_EuanCS.idp			ffddm_geneoCS_3rdlvl.idp		ffddm_partitioning.idp \
cobrameshcavity2d.idp			ffddm_functions.idp			ffddm_globals.idp			semiopen-cylindrical-cavity-mesh.idp \
ffddm.idp				ffddm_geneoCS.idp			ffddm_parameters.idp

EXTRA_DIST=*.edp  





freefem++.pref:Makefile.am
	echo loadpath = \"../examples++-mpi/ ../examples++-hpddm/\" >freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref
	
$(TESTS)::freefem++.pref
install-exec-local:: 
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/idp
	$(INSTALL)  -m 555  $(LIST_IDP)  $(DESTDIR)$(ff_prefix_dir)/idp

FORCE:
