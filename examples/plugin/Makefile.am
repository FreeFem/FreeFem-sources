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
TESTS_ENVIRONMENT = TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw LIST_SKIP="$(LIST_SKIP)"

all-local: freefem++.pref all.edp regtests.edp

if FFTW3
CONDITIONAL_FFTW3 = dfft.edp
else 
SKIP_FFTW3 = dfft.edp
endif !FFTW3

if FORTRAN
CONDITIONAL_FORTRAN = ffnewuoa.edp
else
SKIP_FORTRAN = ffnewuoa.edp
endif !FORTRAN

if GMM
CONDITIONAL_GMM = ilut.edp
else
SKIP_GMM = ilut.edp
endif !GMM

if GSL
CONDITIONAL_GSL = gsl.edp
else
SKIP_GSL = gsl.edp
endif !GSL

if HDF5
CONDITIONAL_HDF5 = iohd5-beam-2d.edp \
	iohd5-beam-3d.edp
else
SKIP_HDF5 = iohd5-beam-2d.edp \
	iohd5-beam-3d.edp
endif !HDF5

if IPOPT
CONDITIONAL_IPOPT = IpoptLap.edp \
	IpoptMinSurf.edp \
	IpoptMinSurfVol.edp \
	IpoptVI.edp \
	IpoptVI2.edp \
	IpoptTest.edp
else
SKIP_IPOPT = IpoptLap.edp \
	IpoptMinSurf.edp \
	IpoptMinSurfVol.edp \
	IpoptVI.edp \
	IpoptVI2.edp \
	IpoptTest.edp
endif !IPOPT

if METIS
CONDITIONAL_METIS = metis.edp \
	schwarz-nm.edp
if SCOTCH
CONDITIONAL_METIS_SCOTCH = scotch.edp
else
SKIP_METIS_SCOTCH = scotch.edp
endif
else
SKIP_METIS = metis.edp schwarz-nm.edp scotch.edp
endif !METIS

if MKL
CONDITIONAL_MKL = PARDISO.edp
else
SKIP_MKL = PARDISO.edp
endif

if MSHMET
CONDITIONAL_MSHMET = aniso.edp
else
SKIP_MSHMET = aniso.edp
endif !MSHMET

if MUMPS_seq
CONDITIONAL_MUMPS_seq = LapMUMPS_seq.edp
else
SKIP_MUMPS_seq = LapMUMPS_seq.edp
endif !MUMPS_seq

if MUMPS_seq
if SUPERLU
CONDITIONAL_MUMPS_SUPERLU = lap-solvers.edp
else
SKIP_MUMPS_SUPERLU = lap-solvers.edp
endif !SUPERLU
else
SKIP_MUMPS_SUPERLU = lap-solvers.edp
endif !MUMPS_seq

if NLOPT
CONDITIONAL_NLOPT = VarIneq2.edp
else
SKIP_NLOPT = VarIneq2.edp
endif !NLOPT

if SUPERLU
CONDITIONAL_SUPERLU = SuperLU.edp
else
SKIP_SUPERLU = SuperLU.edp
endif !SUPERLU

if TETGEN
CONDITIONAL_TETGEN = buildlayermesh.edp \
	checkglumeshcube.edp \
	convexehull3d.edp \
	tetgenholeregion_rugby.edp \
	ttestio.edp
else
SKIP_TETGEN = buildlayermesh.edp \
	checkglumeshcube.edp \
	convexehull3d.edp \
	tetgenholeregion_rugby.edp \
	ttestio.edp
endif !TETGEN

if UMFPACK
CONDITIONAL_UMFPACK = LapLNewSolver.edp \
	LapNewSolver.edp \
	LapUmfpack64.edp
else
SKIP_UMFPACK = LapLNewSolver.edp \
	LapNewSolver.edp \
	LapUmfpack64.edp
endif !UMFPACK

TESTS = APk-AdaptEpsDeltaPk.edp \
	APk-ExplicitPkTest.edp \
	APk-FreeFemQA.edp \
	APk-MetricPk.edp \
	bfstream.edp \
	bilapHCT.edp \
	bilapMorley.edp \
	bilapP3-hct-like.edp \
	bmo.edp \
	ClosePoints.edp \
	cmaes-oven.edp \
	cmaes-VarIneq.edp \
	convect_dervieux.edp \
	convectchacon-3d.edp \
	convectchacon.edp \
	cube.edp \
	distance.edp \
	distance2.edp \
	distance3.edp \
	Element_QF.edp \
	ffrandom.edp \
	funcTemplate.edp \
	glumesh3D.edp \
	IncompleteCholesky.edp \
	isoline.edp \
	lame-TD-NSS.edp \
	lapack.edp \
	LapDG3.edp \
	LapDG4.edp \
	LaplaceP2pnc.edp \
	LaplaceP3-3d.edp \
	LaplaceP3.edp \
	LaplaceP4.edp \
	LaplaceRT1.edp \
	LaplaceRT13d.edp \
	LaplaceRT2.edp \
	layer.edp \
	Leman-mesh.edp \
	load.edp \
	meditddm.edp \
	MetricKuate.edp \
	myfunction2.edp \
	myType.edp \
	NS_P2BR_P0.edp \
	pipe.edp \
	plot-fb-P3.edp \
	plot-fb-P3dc.edp \
	plot-fb-P4.edp \
	plot-fb-P4dc.edp \
	plotfb.edp \
	ppm2rnm.edp \
	provadxw.edp \
	ref.edp \
	RT2-ff.edp \
	SaveHB_sample_z.edp \
	SaveHB_sample.edp \
	shell.edp \
	splitedges.edp \
	splitmesh3.edp \
	splitmesh6.edp \
	Stokes-surface-tension-axi.edp \
	test-ElementMixte.edp \
	testdist.edp \
	testFE_P2BR.edp \
	testFE-P2Bulle3.edp \
	testFE-P2pnc.edp \
	testFE-P3.edp \
	testFE-P3dc.edp \
	testFE-P4.edp \
	testFE-P4dc.edp \
	testFE-PkEdge.edp \
	testFEHCT.edp \
	testFEMorley.edp \
	testp1dcnc.edp \
	thresholdings.edp \
	Schur-Complement.edp \
	$(CONDITIONAL_FFTW3) \
	$(CONDITIONAL_FORTRAN) \
	$(CONDITIONAL_GMM) \
	$(CONDITIONAL_GSL) \
	$(CONDITIONAL_IPOPT) \
	$(CONDITIONAL_METIS) \
	$(CONDITIONAL_METIS_SCOTCH) \
	$(CONDITIONAL_MKL) \
	$(CONDITIONAL_MSHMET) \
	$(CONDITIONAL_MUMPS_seq) \
	$(CONDITIONAL_MUMPS_SUPERLU) \
	$(CONDITIONAL_NLOPT) \
	$(CONDITIONAL_SUPERLU) \
	$(CONDITIONAL_TETGEN) \
	$(CONDITIONAL_UMFPACK) \
    $(SKIP_FFTW3) \
	$(SKIP_FORTRAN) \
	$(SKIP_GMM) \
	$(SKIP_GSL) \
	$(SKIP_IPOPT) \
	$(SKIP_METIS) \
	$(SKIP_METIS_SCOTCH) \
	$(SKIP_MKL) \
	$(SKIP_MSHMET) \
	$(SKIP_MUMPS_seq) \
	$(SKIP_MUMPS_SUPERLU) \
	$(SKIP_NLOPT) \
	$(SKIP_SUPERLU) \
	$(SKIP_TETGEN) \
	$(SKIP_UMFPACK)	
LIST_SKIP= 	$(SKIP_FFTW3) \
	$(SKIP_FORTRAN) \
	$(SKIP_GMM) \
	$(SKIP_GSL) \
	$(SKIP_IPOPT) \
	$(SKIP_METIS) \
	$(SKIP_METIS_SCOTCH) \
	$(SKIP_MKL) \
	$(SKIP_MSHMET) \
	$(SKIP_MUMPS_seq) \
	$(SKIP_MUMPS_SUPERLU) \
	$(SKIP_NLOPT) \
	$(SKIP_SUPERLU) \
	$(SKIP_TETGEN) \
	$(SKIP_UMFPACK)

	
# all test can fail must but clean ????  FH...

XFAIL_TESTS = $(TESTS)

EXTRA_DIST = *.edp all.edp regtests.edp regtests.m4 cube.msh  ch.pts fig.pgm lg.pgm

freefem++.pref: Makefile
	echo loadpath = "../../plugin/seq/" > freefem++.pref
	echo includepath = "../../idp/" >> freefem++.pref

clean-local::
	-rm -r include
	-rm include.done

Ref: makeref.edp
	../../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 regtests.m4 > makeref.edp

all.edp:
	mkdir o
	@test -f ../../src/bin-win32/FreeFem++.exe && ffpp=../../src/bin-win32/FreeFem++.exe || ffpp=../../src/nw/FreeFem++; \
	(echo "NoGraphicWindow=true;NoUseOfWait=true;int verbosityy=verbosity;"; \
	for i in *`ls *.edp|grep -v -E '^(all|regtests|makeref|ref)\.edp$$'` ; do  \
	    if $$ffpp -ns -jc -v 0 $$i  1>/dev/null 2>/dev/null ; then \
		echo ' cout << "--------- file : '$$i' --------------------------------------------------------" << endl;' ;\
		echo "verbosity=verbosityy;" ; \
		echo \{ include \"$$i\"\;\}\; ;\
		echo ' cout << "------------------------------------------------------------------------------ " << endl;' ;\
	    else  \
	        echo ' cout << "--------- error compile  file : '$$i' ---------------------------------------------" << endl;'; \
	      fi; done) > $@
	@awk '/error compile/ {print "Error freefem++ file ," $$6}' $@
skip:
	echo  skip examples:  "$(LIST_SKIP)"
# To check the scripts against their reference values
regtests.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp

clean-local::
