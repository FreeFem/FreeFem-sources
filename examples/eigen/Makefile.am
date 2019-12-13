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

LOG_DRIVER = $(SHELL) $(top_srcdir)/bin/test-driver-ff
TESTS_ENVIRONMENT = TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw SKIP=$(SKIP_TESTS_EIGEN) LIST_CONDITIONAL="$(LIST_CONDITIONAL)"

all-local: all.edp freefem++.pref


TESTS_ARPACK = BeamEigenValue.edp \
	BeamEigenValueperio.edp \
	condition-number.edp \
	free-cyl-axi.edp \
	Lap3dEigenValue.edp \
	LapComplexEigenValue.edp \
	LapEigenValue.edp \
	LapEigenValueFunc.edp \
	LapEigenValueFuncV2.edp \
	LapnosymComplexEigenValue.edp \
	LapnosymEigenValue.edp \
	neuman.edp \
	Stokes-eigen.edp \
	VP-Steklov-Poincare.edp \
	WGM-sphere.edp \
	LapEigenBeltrami.edp 

if ARPACK
CONDITIONAL_ARPACK = $(TESTS_ARPACK)
endif ARPACK

TESTS =	$(TESTS_ARPACK)

LIST_CONDITIONAL = $(CONDITIONAL_ARPACK)

XFAIL_TESTS = LapnosymEigenValue.edp

EXTRA_DIST = *.edp all.edp

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

freefem++.pref:Makefile
	echo includepath = \"../../idp/\" > freefem++.pref
	echo loadpath = \"../../plugin/seq/\" >> freefem++.pref