# ======================================================================
# Laboratoire Jacques-Louis Lions
# Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
# http://www.ljll.math.upmc.fr/lehyaric
# ======================================================================
# This file is part of Freefem++
# 
# Freefem++ is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
# 
# Freefem++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with Freefem++; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
# ======================================================================
# headeralh default=0 freefem make multipleauthors start=19/03/10 upmc

TESTS=  DDM-Schwarz-Lame-2d.edp DDM-Schwarz-Lame-3d.edp DDM-Schwarz-Lap-2dd.edp DDM-Schwarz-Lap-3d.edp	\
	DDM-Schwarz-Stokes-2d.edp LaplaceRT-3d-matrix-mumps.edp MPICGLap.edp MPIGMRES2D.edp MPIGMRES3D.edp MUMPS.edp	\
	NSI3d-carac-mumps.edp NSI3d-carac.edp Stokes-v1-matrix-mumps.edp Stokes-v1-matrix-superludist.edp		\
	Stokes-v2-matrix-mumps.edp Stokes-v3-matrix-mumps.edp VG.edp beam-3d-matrix-superludist.edp			\
	cavityNewtow-MUMPS.edp chaleur3D-hips.edp chaleur3D-mumps.edp chaleur3D-superludist.edp cmaes-mpi-VarIneq.edp	\
	essai.edp mortar-DN-4-mpi.edp testsolver_MUMPS.edp testsolver_SuperLU_DIST.edp 

XFAIL_TESTS=$(TESTS)

LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" SKIP=$(SKIP_TESTS_MPI)

LIST_IDP= AddLayer2d.idp		DDM-Schwarz-macro.idp	MPIGMRESmacro.idp	getARGV.idp \
AddLayer3d.idp		DDM-funcs-v2.idp	MPIplot.idp		mortar-msh.idp Heat3d.idp func-max.idp 

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
MUMPS.edp  dmatrix.hpp  ffmumps_fileparam.txt		ffpastix_iparm_dparm.txt	ffsuperlu_dist_fileparam.txt

# FFCS - list modified to disable some downloaded tools depending on the platform (see
# [[file:../../../configure.ac::tools_problems_all_platforms]] for reasons why some tools may be
# deactivated).

LIST_COMPILE=@TOOL_DYLIB_mumps@ @TOOL_DYLIB_hips@ @TOOL_DYLIB_superludist@ @TOOL_DYLIB_pastix@ MPICG.$(DYLIB_SUFFIX)	\
@TOOL_DYLIB_parms@ mpi-cmaes.$(DYLIB_SUFFIX)

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
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/lib/mpi
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/include
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/etc
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/idp
	-for i in $(LIST_COMPILE); do \
	 if [ -f $$i ] ; then 	$(INSTALL)  -m 555 $$i $(DESTDIR)$(ff_prefix_dir)/lib/mpi; fi; done
	$(INSTALL)  -m 555  $(LIST_IDP)  $(DESTDIR)$(ff_prefix_dir)/idp


freefem++.pref:Makefile
	echo 'loadpath = "./;../examples++-load/"' >freefem++.pref
	echo 'includepath = "../examples++-3d/"' >>freefem++.pref
	echo 'includepath += "../examples++-tutorial/"' >>freefem++.pref
#	echo loadpath += \"./\" >>freefem++.pref
# force the build for freefem++.pref  FH.. 
DDM-Schwarz-Lame-2d.edp:freefem++.pref
# FFCS - cleaning is useful sometimes
clean-local::
	-rm $(LIST_COMPILE)
	-rm *.$(OBJEXT)

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
