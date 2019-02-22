############################################################################
# This file is part of FreeFem++.                                          #
#                                                                          #
# FreeFem++ is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFem++ is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFem++. If not, see <http://www.gnu.org/licenses/>.        #
############################################################################
# SUMMARY : Makefile for FreeFem++
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

LOG_DRIVER=$(SHELL) $(top_srcdir)/bin/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPP) FLAGS_FFPP=-nw

all-local: all.edp regtests.edp

TESTS=APk-AdaptEpsDeltaPk.edp APk-ExplicitPkTest.edp APk-FreeFemQA.edp APk-MetricPk.edp IPOTest.edp IpOptMinSurf.edp IpoptLap.edp IpoptMinSurfVol.edp IpoptVI.edp IpoptVI2.edp LapDG3.edp LapDG4.edp LapLNewSolver.edp LapMUMPS_seq.edp LapNewSolver.edp LapUmfpack64.edp LaplaceP3.edp LaplaceP4.edp LaplaceRT1.edp Leman-mesh.edp MetricKuate.edp NS_P2BR_P0.edp PARDISO.edp SuperLU.edp VarIneq2.edp  bilapMorley.edp bmo.edp buildlayermesh.edp checkglumeshcube.edp cmaes-VarIneq.edp cmaes-oven.edp convect_dervieux.edp convectchacon-3d.edp convectchacon.edp convexehull3d.edp dfft.edp ffnewuoa.edp ffrandom.edp funcTemplate.edp glumesh3D.edp gsl.edp  ilut.edp isoline.edp lame-TD-NSS.edp lap-solvers.edp lapack.edp layer.edp load.edp meditddm.edp metis.edp mshmet.edp myType.edp myfunction2.edp pipe.edp plot-fb-P3.edp plot-fb-P3dc.edp plot-fb-P4.edp plot-fb-P4dc.edp plotfb.edp ppm2rnm.edp provadxw.edp ref.edp refinesphere.edp  schwarz-nm.edp scotch.edp shell.edp splitedges.edp splitmesh3.edp splitmesh6.edp test-ElementMixte.edp testFE-P3.edp testFE-P3dc.edp testFE-P4.edp testFE-P4dc.edp testFE-PkEdge.edp testFE_P2BR.edp testFEMorley.edp tetgencube.edp  tetgenholeregion_rugby.edp thresholdings.edp ttestio.edp bfstream.edp bilapP3-hct-like.edp ClosePoints.edp cube.edp bilapHCT.edp	testFEHCT.edp Stokes-surface-tension-axi.edp Element_QF.edp distance.edp	distance2.edp	distance3.edp testdist.edp testp1dcnc.edp testFE-P2Bulle3.edp LaplaceRT2.edp LaplaceP2pnc.edp testFE-P2pnc.edp LaplaceP3-3d.edp LaplaceRT13d.edp SaveHB_sample.edp SaveHB_sample_z.edp IncompleteCholesky.edp test_LaplacianSurf.edp
# all test can fail must but clean ????  FH...

XFAIL_TESTS=$(TESTS)

EXTRA_DIST=*.edp all.edp regtests.edp regtests.m4

clean-local::
	-rm -r include
	-rm include.done

Ref: makeref.edp
	../../src/nw/FreeFem++-nw makeref.edp

makeref.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 regtests.m4 > makeref.edp

all.edp:
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

# To check the scripts against their reference values
regtests.edp: regtests.m4 ../../etc/config/m4/regtests.m4
	m4 -DASSERT regtests.m4 > regtests.edp

clean-local::
