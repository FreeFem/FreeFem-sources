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

# the sequential examples of FreeVol (make check); the parallel examples *_ddm.edp are run with
# ff-mpirun -np 4 (not in make check); the long ones are shortened by a line // PARAM -Tf ... (bin/test-driver-ff)
TESTS = freevol_advection1d.edp \
	freevol_advection1d_periodic.edp \
	freevol_advection2d.edp \
	freevol_advection3d.edp \
	freevol_advection3d_periodic.edp \
	freevol_advectionS.edp \
	freevol_advectionS_periodic.edp \
	freevol_advection1d_adapt.edp \
	freevol_advection2d_adapt.edp \
	freevol_advection2d_periodic.edp \
	freevol_bbm1d.edp \
	freevol_bbm2d.edp \
	freevol_bbm2d_circle.edp \
	freevol_euler1d.edp \
	freevol_euler1d_periodic.edp \
	freevol_euler2d.edp \
	freevol_euler2d_periodic.edp \
	freevol_euler3d.edp \
	freevol_euler3d_periodic.edp \
	freevol_igd1d.edp \
	freevol_igd1d_periodic.edp \
	freevol_igd2d.edp \
	freevol_igd2d_periodic.edp \
	freevol_igd3d.edp \
	freevol_igd3d_periodic.edp \
	freevol_sw2d_circle.edp \
	freevol_swtopo1d.edp \
	freevol_swtopo1d_periodic.edp \
	freevol_swtopo2d.edp \
	freevol_swtopo2d_periodic.edp \
	freevol_swtopo3d.edp \
	freevol_swtopo3d_periodic.edp \
	freevol_swtopoS.edp \
	freevol_swtopoS_periodic.edp

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
