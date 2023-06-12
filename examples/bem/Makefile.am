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


TESTS_MPI_BEMTOOL = 

if MPI
CONDITIONAL_MPI =  $(TESTS_MPI)




if BEMTOOL
CONDITIONAL_MPI_BEMTOOL = $(TESTS_MPI_BEMTOOL)
if MUMPS
CONDITIONAL_MPI_BEMTOOL_MUMPS = $(TESTS_MPI_BEMTOOL_MUMPS)
endif MUMPS
if METIS
CONDITIONAL_MPI_BEMTOOL_METIS = $(TESTS_MPI_BEMTOOL_METIS)
endif METIS
endif BEMTOOL
endif MPI

TESTS = $(TESTS_MPI) \
	$(TESTS_MPI_BEMTOOL_METIS) \
	$(TESTS_MPI_BEMTOOL_MUMPS) \
	$(TESTS_MPI_BEMTOOL)

LIST_CONDITIONAL = $(CONDITIONAL_MPI) \
	$(CONDITIONAL_MPI_BEMTOOL_METIS) \
	$(CONDITIONAL_MPI_BEMTOOL_MUMPS) \
	$(CONDITIONAL_MPI_BEMTOOL)
	

XFAIL_TESTS = $(TESTS)

EXTRA_DIST = *.edp

freefem++.pref: Makefile
	echo "loadpath = ../../plugin/mpi/" > freefem++.pref
	echo "loadpath += ../../plugin/seq/" >>freefem++.pref
	echo "includepath = ../../idp/" >> freefem++.pref

clean-local::
