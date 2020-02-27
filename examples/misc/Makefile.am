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

LOG_DRIVER=$(SHELL) $(top_srcdir)/bin/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw

all-local: all.edp regtests.edp freefem++.pref

TESTS=aaa-adp.edp aadaptation.edp \
	aalapacien.edp aalaplace-nc.edp \
	aamove.edp aaRT.edp arrayoFVh.edp \
	bilap.edp ccc-adp.edp \
	D2.edp demo.edp \
	demo1.edp FE-medit.edp \
	funct.edp include.edp \
	lap_mat.edp lapacienprecon.edp \
	NSP1P1.edp NSP1P1b.edp NSP1P1bl.edp \
	NSP1P2.edp parareal.edp \
	renumbering.edp Richard.edp \
	Stokes-P1nc-P0+P1.edp testadp.edp \
	teste.edp testFE.edp \
	wafer-heating-laser-axi.edp  torture-search2d.edp \
	lap3-cpu.edp
# if .. no arpack or too long
XFAIL_TESTS=NSP1P1.edp

EXTRA_DIST=*.edp all.edp regtests.edp regtests.m4 ref.edp torture.msh speedtest.sh

all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# To create a new set of reference values in "ref.edp"
clean-local:
	-rm *.fg *.eps *~ *.ps Th?.faces Th?.points Th?.BB Th?.mesh.gmsh ListOf?llocPtr*.bin ffglut*.ppm 	toto.txt

Ref: makeref.edp
	../../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 regtests.m4 > makeref.edp

freefem++.pref:Makefile
	echo includepath = \"../../idp/\" > freefem++.pref
	echo loadpath = \"../../plugin/seq/\" >> freefem++.pref

# To check the scripts against their reference values
regtests.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp

FORCE:
