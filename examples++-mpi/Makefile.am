# $Id$
TESTS=DDM-Schwarz-Lame-2d.edp DDM-Schwarz-Lame-3d.edp DDM-Schwarz-Lap-2dd.edp DDM-Schwarz-Lap-3d.edp DDM-Schwarz-Stokes-2d.edp LaplaceRT-3d-matrix-mumps.edp MPICGLap.edp MPIGMRES2D.edp MPIGMRES3D.edp MUMPS.edp NSI3d-carac-mumps.edp NSI3d-carac.edp Stokes-v1-matrix-mumps.edp Stokes-v1-matrix-superludist.edp Stokes-v2-matrix-mumps.edp Stokes-v3-matrix-mumps.edp VG.edp beam-3d-matrix-superludist.edp cavityNewtow-MUMPS.edp chaleur3D-hips.edp chaleur3D-mumps.edp chaleur3D-superludist.edp cmaes-mpi-VarIneq.edp essai.edp  mortar-DN-4-mpi.edp testsolver_MUMPS.edp testsolver_SuperLU_DIST.edp testsolver_pastix.edp
XFAIL_TESTS=$(TESTS)

LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" SKIP=$(SKIP_TESTS_MPI)

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

# FFCS - list modified to disable some downloaded tools depending on the platform (see
# [[file:../../../configure.ac::tools_problems_all_platforms]] for reasons why some tools may be
# deactivated).

LIST_COMPILE=@FFCS_DYLIB_mumps@ @FFCS_DYLIB_hips@ @FFCS_DYLIB_superludist@ @FFCS_DYLIB_pastix@ MPICG.$(DYLIB_SUFFIX)	\
@FFCS_DYLIB_parms@ mpi-cmaes.$(DYLIB_SUFFIX)

if FFCS_MPIOK

# FFCS - do not set loadpath to be able to run an external version of FF on the examples in this directory with
# [[../../mkffref]]

all-local:$(LIST_COMPILE)
	@echo Warning missing mpi plugin: `for i in $(LIST_COMPILE); do if test ! -s $i ; then j=1; echo "$i," ;fi; done`
	echo "finish compile load mpi solver !"
else
all-local:
endif

.cpp.$(DYLIB_SUFFIX): ../examples++-load/ff-c++ ../examples++-load/WHERE_LIBRARY-download
#	FFCS needs an error exit code to make sure that all libraries are correctly compiled
#
#	FFCS on Windows inserts all MPI options (include, lib, ...) in the $MPICXX compiler script, instead of "$CXX $MPI_xxx",
#	but it needs '-mpi' to do that
#
###	-../examples++-load/ff-c++ -auto  $< 
	../examples++-load/ff-c++ -auto -mpi $<

install-exec-local:: 
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/lib
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/include
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/etc
	-for i in $(LIST_COMPILE); do \
	 if [ -f $$i ] ; then 	$(INSTALL)  -m 555 $$i $(DESTDIR)$(ff_prefix_dir)/lib; fi; done
	$(INSTALL)  -m 555  $(LIST_IDP)  $(DESTDIR)$(ff_prefix_dir)/idp


freefem++.pref:Makefile
	echo loadpath = \"../examples++-load/\" >freefem++.pref
	echo includepath = \"../examples++-3d/\" >>freefem++.pref
	echo includepath += \"../examples++-tutorial/\" >>freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref

# FFCS - cleaning is useful sometimes
clean-local::
	-rm *.$(DYLIB_SUFFIX)
	-rm *.$(OBJEXT)
