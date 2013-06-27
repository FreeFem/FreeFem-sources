# $Id$
TESTS=DDM-Schwarz-Lame-2d.edp DDM-Schwarz-Lame-3d.edp DDM-Schwarz-Lap-2dd.edp DDM-Schwarz-Lap-3d.edp DDM-Schwarz-Stokes-2d.edp LaplaceRT-3d-matrix-mumps.edp MPICGLap.edp MPIGMRES2D.edp MPIGMRES3D.edp MUMPS.edp NSI3d-carac-mumps.edp NSI3d-carac.edp Stokes-v1-matrix-mumps.edp Stokes-v1-matrix-superludist.edp Stokes-v2-matrix-mumps.edp Stokes-v3-matrix-mumps.edp VG.edp beam-3d-matrix-superludist.edp cavityNewtow-MUMPS.edp chaleur3D-hips.edp chaleur3D-mumps.edp chaleur3D-superludist.edp cmaes-mpi-VarIneq.edp essai.edp  mortar-DN-4-mpi.edp testsolver_MUMPS.edp testsolver_SuperLU_DIST.edp testsolver_pastix.edp
XFAILS_TESTS=$(TESTS)

LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4" SKIP=$(SKIP_TESTS_MPI)

LIST_IDP= AddLayer2d.idp		DDM-Schwarz-macro.idp	MPIGMRESmacro.idp	getARGV.idp \
AddLayer3d.idp		DDM-funcs-v2.idp	MPIplot.idp		mortar-msh.idp Heat3d.idp

EXTRA_DIST=*.edp *.idp  regtests.sh  \
MPICG.cpp \
MUMPS_FreeFem.cpp \
complex_SuperLU_DIST_FreeFem.cpp \
complex_pastix_FreeFem.cpp \
dSuperLU_DIST.cpp \
ffsuperludistoption-1.hpp \
ffsuperludistoption.hpp \
generaldefs.h \
hips_FreeFem.cpp \
hypre_FreeFem.cpp \
interfacepastix.cpp \
parms_FreeFem.cpp \
real_SuperLU_DIST_FreeFem.cpp \
real_pastix_FreeFem.cpp mpi-cmaes.cpp \
ffmumps_fileparam.txt		ffpastix_iparm_dparm.txt	ffsuperlu_dist_fileparam.txt \
MUMPS.cpp   \
MUMPS.edp  dmatrix.hpp              


LIST_COMPILE=MUMPS_FreeFem.$(DYLIB_SUFFIX) hips_FreeFem.$(DYLIB_SUFFIX) \
 interfacepastix.$(DYLIB_SUFFIX) \
complex_SuperLU_DIST_FreeFem.$(DYLIB_SUFFIX) real_SuperLU_DIST_FreeFem.$(DYLIB_SUFFIX) \
complex_pastix_FreeFem.$(DYLIB_SUFFIX) real_pastix_FreeFem.$(DYLIB_SUFFIX) \
dSuperLU_DIST.$(DYLIB_SUFFIX) MPICG.$(DYLIB_SUFFIX) parms_FreeFem.$(DYLIB_SUFFIX) \
mpi-cmaes.$(DYLIB_SUFFIX) MUMPS.$(DYLIB_SUFFIX)  

all-local:$(LIST_COMPILE) freefem++.pref
	echo "finish compile load mpi solver !"

.cpp.$(DYLIB_SUFFIX): ../examples++-load/ff-c++ ../examples++-load/WHERE_LIBRARY-download
	- ../examples++-load/ff-c++ -auto  $< 
install-exec-local:: 
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/lib
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/include
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/etc
	-for i in $(LIST_COMPILE); do \
	 if [ -f $$i ] ; then 	$(INSTALL)  -m 555 $$i $(DESTDIR)$(ff_prefix_dir)/lib; fi; done
	$(INSTALL)  -m 555  $(LIST_IDP)  $(DESTDIR)$(ff_prefix_dir)/idp


freefem++.pref:
	echo loadpath = \"../examples++-load/\" >freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref
DDM-Schwarz-Lame-3d.edp:freefem++.pref: