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
TESTS_ENVIRONMENT = TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw"LIST_CONDITIONAL="$(LIST_CONDITIONAL)"

all-local: freefem++.pref


TESTS_MPI = diffusion-2d-thirdlevelgeneo.edp \
	diffusion-3d-minimal-ddm.edp\
	diffusion-3d-minimal-direct.edp \
	diffusion-3d-simple.edp \
	elasticity-3d-simple.edp \
	elasticity-3d-thirdlevelgeneo.edp \
	Helmholtz-2d-HPDDM-BGMRES.edp \
	Helmholtz-2d-simple.edp \
	Helmholtz-3d-simple.edp \
	Maxwell-3d-simple.edp \
	natural_convection_3D_obstacle.edp \
	natural_convection.edp \
	Richards-2d.edp \
	heat-torus-3d-surf.edp


TESTS_METIS = Maxwell_Cobracavity.edp

if METIS
CONDITIONAL_METIS = $(TESTS_METIS)
endif METIS

TESTS = $(TESTS_MPI) \
	$(TESTS_METIS)

LIST_CONDITIONAL = $(CONDITIONAL_OTHER) \
	$(CONDITIONAL_METIS)


XFAIL_TESTS = $(TESTS)

EXTRA_DIST = *.edp

freefem++.pref: Makefile.am
	echo loadpath = \"../../plugin/mpi/\" > freefem++.pref
	echo loadpath += \"../../plugin/seq/\" >> freefem++.pref
	echo includepath = \"../../idp/\" >> freefem++.pref

$(TESTS)::freefem++.pref

FORCE:
