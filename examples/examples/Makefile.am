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
TESTS_ENVIRONMENT = TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw LIST_CONDITIONAL="$(LIST_CONDITIONAL)"

all-local: all.edp regtests.edp freefem++.pref

TESTS_OTHER = BlackScholes2D.md \
	condensor.md \
	convects.md \
	Heat.md \
	heatex.md \
	lame.md \
	membrane.md \
	membranerror.edp \
	muwave.md \
	NSNewton.md \
	NSprojection.md \
	optimcontrol.md \
	potential.md \
	schwarz.md \
	stokes.md \
	stokes_composite.md \
	stokes_periodic_composite.md \
	test1.edp \
	thermal.md \
	thermic.md \
	dist-projection.md

TESTS_ARPARCK = sound.md


if ARPACK
CONDITIONAL_ARPACK = $(TESTS_ARPACK)
endif ARPACK

TESTS = $(TESTS_OTHER) \
	$(TESTS_ARPACK)
	
LIST_CONDITIONAL = $(TESTS_OTHER) \
	$(TESTS_ARPACK)



# if .. no arpack or too long
XFAIL_TESTS = NSprojection.md

EXTRA_DIST = *.edp all.edp regtests.edp regtests.m4 ref.edp

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# To create a new set of reference values in "ref.edp"
Ref: makeref.edp
	../../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 regtests.m4 > makeref.edp

# To check the scripts against their reference values
regtests.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp

freefem++.pref:Makefile
	echo includepath = \"../../idp/\" > freefem++.pref
	echo loadpath = \"../../plugin/seq/\" >> freefem++.pref

FORCE:
