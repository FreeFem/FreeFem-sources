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
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" SKIP=$(SKIP_TESTS_MPI)

all-local: freefem++.pref

TESTS=  DDM-Schwarz-Lame-2d.edp DDM-Schwarz-Lame-3d.edp DDM-Schwarz-Lap-2dd.edp DDM-Schwarz-Lap-3d.edp	\
	DDM-Schwarz-Stokes-2d.edp LaplaceRT-3d-matrix-mumps.edp MPICGLap.edp MPIGMRES2D.edp MPIGMRES3D.edp MUMPS.edp	\
	NSI3d-carac-mumps.edp NSI3d-carac.edp Stokes-v1-matrix-mumps.edp 	\
	Stokes-v2-matrix-mumps.edp Stokes-v3-matrix-mumps.edp VG.edp 	\
	cavityNewtow-MUMPS.edp  cmaes-mpi-VarIneq.edp	\
	essai.edp mortar-DN-4-mpi.edp testsolver_MUMPS.edp  parmetis.edp parmetis-3d.edp \
	LapMPIEigenValue.edp htool.edp

XFAIL_TESTS=$(TESTS)

EXTRA_DIST=*.edp

freefem++.pref:Makefile
	echo "loadpath = ../../plugin/mpi/" >> freefem++.pref
	echo "loadpath += ../../plugin/seq/" > freefem++.pref
	echo "includepath = ../../idp/" >> freefem++.pref

clean-local::
