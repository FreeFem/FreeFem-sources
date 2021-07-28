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

TESTS=forall.edp adapt.edp adaptindicatorP1.edp adaptindicatorP2.edp algo.edp array.edp a_tutorial.edp beam.edp  calculus.edp cavity.edp convect2.edp convect-apt.edp convect.edp dumptable.edp ex-vf.edp FE.edp fluidStructAdapt.edp fluidStruct.edp freeboundary.edp freeboundary-weak.edp LapDG2.edp Laplace.edp LaplaceP1bis.edp LaplaceP1.edp LaplaceP1P2h.edp LaplaceRT.edp mesh.edp movemesh.edp nolinear-elas.edp NSUzawaCahouetChabart.edp onde.edp periodic4.edp Periodic.edp plot.edp readmesh.edp region.edp saverestore.edp schwarz-gc.edp schwarz-no-overlap.edp schwarz-overlap.edp sparse-matrix.edp sparse-cmatrix.edp StokesUzawa.edp tablefunction.edp intlevelset.edp mesh-internal.edp minlenedge.edp AdjointSolve.edp

EXTRA_DIST=*.edp aile.msh xyf all.edp regtests.edp regtests.m4 ref.edp


all.edp:
	(echo "NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$' | LC_ALL=C sort` ; do \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

clean-local:
	-rm *.fg *.eps *~ *.ps *.mesh *.mesh.gmsh f.txt	u.txt A.matrix	mm.matrix mm.b sphere-a.faces	sphere-a.points	sphere.faces	sphere.points toto.Th		toto.am_fmt	toto.dbg.gmsh toto.Th.gmsh	toto.dbg	toto.msh plot.gp ListOfAllocPtr-8.bin ffglut*.ppm th.msh Th1.msh		Th12.msh	th.msh 	emptymesh-2.msh	 thermic.dat

# To create a new set of reference values in "ref.edp"
Ref: makeref.edp freefem++.pref
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
