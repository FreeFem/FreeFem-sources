############################################################################
# This file is part of FreeFEM.                                            #
#                                                                          #
# FreeFEM is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFEM is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          #
############################################################################
# SUMMARY : Makefile for FreeFEM
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...
.NOTPARALLEL:
LOG_DRIVER = $(SHELL) $(top_srcdir)/bin/test-driver-ff
TESTS_ENVIRONMENT = TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" LIST_CONDITIONAL="$(LIST_CONDITIONAL)"

all-local: freefem++.pref


TESTS_MPI_BEMTOOL = Helmholtz_circle_Dirichlet_all_direct.edp \
	Helmholtz_Cobracavity.edp \
	Maxwell_cube_EFIE.edp

TESTS_MPI_BEMTOOL_MUMPS = Helmholtz-2d-FEM-BEM-coupling-MUMPS.edp \
	Helmholtz-2d-FEM-BEM-coupling-MUMPS-composite.edp

TESTS_MPI_BEMTOOL_GSL = Helmholtz_circle_Dirichlet_simple.edp

if MPI
if BEMTOOL
CONDITIONAL_MPI_BEMTOOL = $(TESTS_MPI_BEMTOOL)

if MUMPS
CONDITIONAL_MPI_BEMTOOL_MUMPS = $(TESTS_MPI_BEMTOOL_MUMPS)
endif MUMPS

if METIS
CONDITIONAL_MPI_BEMTOOL_METIS = $(TESTS_MPI_BEMTOOL_METIS)
endif METIS

if GSL
CONDITIONAL_MPI_BEMTOOL_GSL = $(TESTS_MPI_BEMTOOL_GSL)
endif GSL

endif BEMTOOL
endif MPI

TESTS = $(TESTS_MPI_BEMTOOL) \
	$(TESTS_MPI_BEMTOOL_MUMPS) \
	$(TESTS_MPI_BEMTOOL_GSL) \
	$(TESTS_MPI_BEMTOOL_METIS)

LIST_CONDITIONAL = $(CONDITIONAL_MPI_BEMTOOL) \
	$(CONDITIONAL_MPI_BEMTOOL_MUMPS) \
	$(CONDITIONAL_MPI_BEMTOOL_GSL) \
	$(CONDITIONAL_MPI_BEMTOOL_METIS)

XFAIL_TESTS = $(TESTS)

EXTRA_DIST = *.edp

freefem++.pref: Makefile
	echo "loadpath = ../../plugin/mpi/" > freefem++.pref
	echo "loadpath += ../../plugin/seq/" >>freefem++.pref
	echo "includepath = ../../idp/" >> freefem++.pref

clean-local::
