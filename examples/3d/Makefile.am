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

TESTS_OTHER = 3d-leman.md \
	ArrayFE-3d.edp \
	beam-3d.md \
	cone.md \
	convect-3d.md \
	cube-period.md \
	cylinder-3d.md \
	EqPoisson.md \
	extract-boundary3d.md \
	first.md \
	intlevelset3d.md \
	Lac.md \
	Laplace3d.md \
	LaplaceRT-3d.md \
	meditddm.edp \
	NSI3d-carac.md \
	NSI3d.edp \
	p.edp \
	periodic-3d.md \
	Poisson.edp \
	pyramide.md \
	sphere2.md \
	sphere6.md \
	Stokes.md \
	TruncLac.edp \
	crack-3d.md \
	periodic3.edp \
	LapDG3d1.md \
	LapDG3d.md \
	Connectivity-3d.md \
	Elasticity-simple-support-BC.md

TESTS_METIS = schwarz-nm-3d.md

TESTS_TETGEN = cylinder.edp \
	Poisson3d.md \
	Poisson-cube-ballon.md \
	refinesphere.md \
	tetgencube.edp \
	tetgenholeregion.md

TESTS_TETGEN_MMG = Sphere-Isocahedron.md Period-Poisson-cube-ballon.md \
	sphereincube.md

TESTS_TETGEN_MMG_MSHMET = Laplace-Adapt-aniso-3d.md

TESTS_TETGEN_GSL = bottle.md

TESTS_TETGEN_MMG3D_MMG = fallingspheres.md

TESTS_TETGEN_MSHMET = Laplace-Adapt-3d.md

if METIS
CONDITIONAL_METIS = $(TESTS_METIS)
endif METIS

if TETGEN
CONDITIONAL_TETGEN = $(TESTS_TETGEN)
if MMG
CONDITIONAL_TETGEN_MMG = $(TESTS_TETGEN_MMG)
if MSHMET
CONDITIONAL_TETGEN_MMG_MSHMET = $(TESTS_TETGEN_MMG_MSHMET)
endif MSHMET
endif MMG

if GSL
CONDITIONAL_TETGEN_GSL = $(TESTS_TETGEN_GSL)
endif GSL
if MMG3D
if MMG
CONDITIONAL_TETGEN_MMG3D_MMG = $(TESTS_TETGEN_MMG3D_MMG)
endif MMG
endif MMG3D
if MSHMET
CONDITIONAL_TETGEN_MSHMET = $(TESTS_TETGEN_MSHMET)
endif MSHMET

endif TETGEN

TESTS = $(TESTS_OTHER) \
	$(TESTS_METIS) \
	$(TESTS_TETGEN) \
	$(TESTS_TETGEN_FREEYAMS) \
	$(TESTS_TETGEN_FREEYAMS_MMG3D_MSHMET) \
	$(TESTS_TETGEN_GSL)  \
	$(TESTS_TETGEN_MMG3D_MMG) \
	$(TESTS_TETGEN_MSHMET)

LIST_CONDITIONAL = $(TESTS_OTHER) \
	$(CONDITIONAL_METIS) \
	$(CONDITIONAL_TETGEN) \
	$(CONDITIONAL_TETGEN_MMG) \
	$(CONDITIONAL_TETGEN_MMG_MSHMET) \
	$(CONDITIONAL_TETGEN_GSL) \
	$(CONDITIONAL_TETGEN_MMG3D_MMG) \
	$(CONDITIONAL_TETGEN_MSHMET)

XFAIL_TESTS =

EXTRA_DIST = *.edp \
	all.edp \
	regtests.edp \
	regtests.m4 \
	ref.edp \
	dodecaedre01.mesh \
	lac-leman-v4.msh

all.edp: Makefile
	(echo "NoGraphicWindow=true;NoUseOfWait=true;int verbosityy=verbosity;int MEM1234=storageused();"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;searchMethod=0;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << " mem leak = " <<storageused() - MEM1234 << endl << "------------------------------------------------------------------------------ " << endl;' ;\
	done) > $@

# To create a new set of reference values in "ref.edp"
Ref: makeref.edp freefem++.pref
	export FF_LOADPATH=../../plugin/seq/.;../../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 regtests.m4 > makeref.edp

freefem++.pref:Makefile
	echo includepath = \"../../idp/\" > freefem++.pref
	echo loadpath = \"../../plugin/seq/\" >> freefem++.pref

3d-leman.md: freefem++.pref

# To check the scripts against their reference values
regtests.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp
echo:
	echo $(LIST_CONDITIONAL)
	echo $(LIST)
	
FORCE:
