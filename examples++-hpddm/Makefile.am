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

TESTS=  diffusion-3d.edp	elasticity-3d.edp	schwarz.edp

XFAIL_TESTS=$(TESTS)

LOG_DRIVER=$(SHELL) $(top_srcdir)/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" SKIP=$(SKIP_TESTS_MPI)

LIST_IDP=  \
getARGV.idp			schwarz-skel.idp \
additional.idp		macro_3d.idp		schwarz-additional_macro.idp	 \
argv.idp			macro.idp			schwarz-add-macro.idp		schwarz-macro.idp

EXTRA_DIST= $(TESTS) $(LIST_IDP) ./hpddm.cpp ./skeleton.cpp ./ff-hpddm.hpp ./removeDOF.cpp ./mymsh3.cpp ./PETSc.hpp ./init.hpp \
include-hpddm/ARPACK.hpp		include-hpddm/MKL_PARDISO.hpp		include-hpddm/dmatrix.hpp		include-hpddm/option_impl.hpp \
include-hpddm/BDD.hpp			include-hpddm/MUMPS.hpp			include-hpddm/eigensolver.hpp		include-hpddm/preconditioner.hpp \
include-hpddm/BLAS.hpp			include-hpddm/PaStiX.hpp		include-hpddm/enum.hpp			include-hpddm/schur.hpp \
include-hpddm/FETI.hpp			include-hpddm/ScaLAPACK.hpp		include-hpddm/iterative.hpp		include-hpddm/schwarz.hpp \
include-hpddm/HPDDM.hpp			include-hpddm/SuiteSparse.hpp		include-hpddm/matrix.hpp		include-hpddm/subdomain.hpp \
include-hpddm/Hypre.hpp			include-hpddm/coarse_operator.hpp	include-hpddm/operator.hpp		include-hpddm/wrapper.hpp \
include-hpddm/LAPACK.hpp		include-hpddm/coarse_operator_impl.hpp	include-hpddm/option.hpp


# FFCS - list modified to disable some downloaded tools depending on the platform (see
# [[file:../../../configure.ac::tools_problems_all_platforms]] for reasons why some tools may be
# deactivated).

LIST_COMPILE=@TOOL_DYLIB_hpddm@ 

if FFCS_MPIOK

# FFCS - do not set loadpath to be able to run an external version of FF on the examples in this directory with
# [[../../mkffref]]

all-local:$(LIST_COMPILE)
	@echo Warning missing mpi plugin: `for i in $(LIST_COMPILE); do if test ! -s $i ; then j=1; echo "$i," ;fi; done`
	echo "finish compile load hpddm  solver !"
else
all-local:
endif

.cpp.$(DYLIB_SUFFIX): ../examples++-load/ff-c++ ../examples++-load/WHERE_LIBRARY-download ../examples++-load/WHERE_LIBRARY-config
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
	-for i in $(LIST_COMPILE); do \
	 if [ -f $$i ] ; then 	$(INSTALL)  -m 555 $$i $(DESTDIR)$(ff_prefix_dir)/lib/mpi; fi; done
	$(INSTALL)  -m 555  $(LIST_IDP)  $(DESTDIR)$(ff_prefix_dir)/idp


freefem++.pref:Makefile
	echo loadpath = \"../examples++-mpi/\" >freefem++.pref
	echo loadpath += \"../examples++-load/\" >>freefem++.pref
	echo includepath = \"../examples++-3d/\" >>freefem++.pref
	echo includepath += \"../examples++-tutorial/\" >>freefem++.pref
	echo loadpath += \"./\" >>freefem++.pref
# force the build for freefem++.pref  FH.. 
$(TEST) :freefem++.pref
# FFCS - cleaning is useful sometimes
clean-local::
	-rm $(LIST_COMPILE)
	-rm *.$(OBJEXT)

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
