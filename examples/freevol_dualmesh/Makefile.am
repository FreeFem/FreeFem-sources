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

all-local: all.edp freefem++.pref

# the sequential examples of FreeVol on the dual mesh (make check); the parallel examples *_ddm.edp are run with
# ff-mpirun -np 4 (not in make check); the long ones are shortened by a line // PARAM -Tf ... (bin/test-driver-ff)
TESTS = freevol_advectiondm.edp \
	freevol_advectiondm2d_periodic.edp \
	freevol_advectiondm1d.edp \
	freevol_advectiondm1d_periodic.edp \
	freevol_advectiondm1d_adapt.edp \
	freevol_advectiondm2d_adapt.edp \
	freevol_advectiondm3d.edp \
	freevol_advectiondm3d_periodic.edp \
	freevol_advectiondmS.edp \
	freevol_advectiondmS_periodic.edp \
	freevol_eulerdm.edp \
	freevol_eulerdm2d_periodic.edp \
	freevol_eulerdm1d.edp \
	freevol_eulerdm1d_periodic.edp \
	freevol_eulerdm3d.edp \
	freevol_eulerdm3d_periodic.edp \
	freevol_eulerdmS_periodic.edp \
	freevol_gasdmS.edp \
	freevol_igddm.edp \
	freevol_igddm2d_periodic.edp \
	freevol_igddm1d.edp \
	freevol_igddm1d_periodic.edp \
	freevol_igddm3d.edp \
	freevol_igddm3d_periodic.edp \
	freevol_igddmS_periodic.edp \
	freevol_swtopodm.edp \
	freevol_swtopodm2d_periodic.edp \
	freevol_swtopodm2d_adapt.edp \
	freevol_swtopodm1d.edp \
	freevol_swtopodm1d_periodic.edp \
	freevol_swtopodm1d_adapt.edp \
	freevol_swtopodm3d.edp \
	freevol_swtopodm3d_periodic.edp \
	freevol_swtopodmS.edp \
	freevol_swtopodmS_periodic.edp \
	test_P0dual.edp \
	test_P0dualS.edp \
	test_meshD.edp \
	testdualmesh.edp

LIST_CONDITIONAL = $(TESTS)

EXTRA_DIST = *.edp all.edp

all.edp: Makefile
	(echo "NoGraphicWindow=true;NoUseOfWait=true;int verbosityy=verbosity;int MEM1234=storageused();"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;searchMethod=0;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << " mem leak = " <<storageused() - MEM1234 << endl << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@


freefem++.pref:Makefile
	echo includepath = \"../../idp/\" > freefem++.pref
	echo loadpath = \"../../plugin/seq/\" >> freefem++.pref
	echo loadpath += \"../../plugin/mpi/\" >> freefem++.pref

FORCE:
