# $Id$  

all-local:   freefem++.pref

TESTS=	Maxwell-3d-simple.edp		Maxwell_Cobracavity.edp			diffusion-3d-simple.edp \
diffusion-3d-minimal-ddm.edp	 \
Helmholtz-2d-simple.edp			diffusion-3d-minimal-direct.edp		elasticity-3d-simple.edp \
elasticity-3d-thirdlevelgeneo.edp \
Richards-2d.edp					natural_convection.edp \
Helmholtz-3d-simple.edp			diffusion-2d-thirdlevelgeneo.edp	natural_convection_3D_obstacle.edp \
Helmholtz-2d-HPDDM-BGMRES.edp

LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" SKIP=$(SKIP_TESTS_MPI)

LIST_IDP= cobrameshcavity.idp			ffddm_EuanCS.idp			ffddm_geneoCS_3rdlvl.idp		ffddm_partitioning.idp \
ffddm_functions.idp \
ffddm.idp				ffddm_geneoCS.idp			ffddm_parameters.idp

EXTRA_DIST=*.edp  





freefem++.pref:Makefile.am
	echo loadpath = \"../examples++-mpi/\" >freefem++.pref
	echo loadpath += \"../examples++-hpddm/\" >>freefem++.pref
	echo loadpath += \"../examples++-load/\" >>freefem++.pref
	echo includepath = \"../examples++-mpi/\" >>freefem++.pref
	echo includepath += \"../examples++-hpddm/\" >>freefem++.pref
	echo includepath += \"../examples++-3d/\" >>freefem++.pref
	echo includepath += \"../examples++-tutorial/\" >>freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref
	
$(TESTS)::freefem++.pref
install-exec-local:: 
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/idp
	$(INSTALL)  -m 555  $(LIST_IDP)  $(DESTDIR)$(ff_prefix_dir)/idp

FORCE:
