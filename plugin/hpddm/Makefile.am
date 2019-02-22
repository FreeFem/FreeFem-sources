############################################################################
# This file is part of FreeFem++.                                          #
#                                                                          #
# FreeFem++ is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFem++ is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFem++. If not, see <http://www.gnu.org/licenses/>.        #
############################################################################
# SUMMARY : Makefile for FreeFem++
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

# FFCS - list modified to disable some downloaded tools depending on the platform (see
# [[file:../../../configure.ac::tools_problems_all_platforms]] for reasons why some tools may be
# deactivated).

LIST_COMPILE=@TOOL_DYLIB_hpddm@

if FFCS_MPIOK

# FFCS - do not set loadpath to be able to run an external version of FF on the examples in this directory with
# [[../../mkffref]]

all-local:$(LIST_COMPILE)
	@echo Warning missing mpi plugin: `for i in $(LIST_COMPILE); do if test ! -s $i ; then j=1; echo "$i," ;fi; done`
	echo "finish compile load hpddm solver !"
else
all-local:
endif

.cpp.$(DYLIB_SUFFIX): ../seq/ff-c++ ../seq/WHERE_LIBRARY-download ../seq/WHERE_LIBRARY-config
#	FFCS needs an error exit code to make sure that all libraries are correctly compiled
#
#	FFCS on Windows inserts all MPI options (include, lib, ...) in the $MPICXX compiler script, instead of "$CXX $MPI_xxx",
#	but it needs '-mpi' to do that
#
###	-../seq/ff-c++ -auto $<
	../seq/ff-c++ -auto -mpi $<

install-exec-local::
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/lib/mpi
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/include
	$(mkinstalldirs) -m 755 $(DESTDIR)$(ff_prefix_dir)/etc
	-for i in $(LIST_COMPILE); do \
	 if [ -f $$i ] ; then 	$(INSTALL) -m 555 $$i $(DESTDIR)$(ff_prefix_dir)/lib/mpi; fi; done


# force the build for freefem++.pref FH.. (do not works 2016)
#$(TEST) :freefem++.pref
# FFCS - cleaning is useful sometimes
clean-local::
	-rm $(LIST_COMPILE)  hpddm.dylib hpddm_substructuring.dylib PETSc.dylib PETSc-complex.dylib SLEPc-complex.dylib SLEPc.dylib
	-cd ../mpi/.; rm $(LIST_COMPILE)  hpddm.dylib hpddm_substructuring.dylib PETSc.dylib PETSc-complex.dylib SLEPc-complex.dylib SLEPc.dylib
	-rm *.$(OBJEXT)

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
