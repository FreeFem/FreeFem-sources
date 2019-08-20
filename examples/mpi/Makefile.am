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
.NOTPARALLEL:
LOG_DRIVER = $(SHELL) $(top_srcdir)/bin/test-driver-ff
TESTS_ENVIRONMENT = TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" SKIP=$(SKIP_TESTS_MPI)

all-local: freefem++.pref

if MPI
CONDITIONAL_MPI =  \
	cmaes-mpi-VarIneq.edp \
	essai.edp \
	mortar-DN-4-mpi.edp \
	MPICGLap.edp

if HTOOL
CONDITIONAL_MPI_HTOOL = htool.edp
endif HTOOL

if METIS
CONDITIONAL_MPI_METIS = DDM-Schwarz-Lame-2d.edp \
	DDM-Schwarz-Lame-3d.edp \
	DDM-Schwarz-Lap-2dd.edp \
	DDM-Schwarz-Lap-3d.edp \
	DDM-Schwarz-Stokes-2d.edp \
	MPIGMRES2D.edp \
	MPIGMRES3D.edp

if PARMETIS
CONDITIONAL_MPI_METIS_PARMETIS = parmetis-3d.edp \
	parmetis.edp
endif PARMETIS
endif METIS

if MUMPS
CONDITIONAL_MPI_MUMPS = cavityNewtow-MUMPS.edp \
	LaplaceRT-3d-matrix-mumps.edp \
	MUMPS.edp \
	NSI3d-carac-mumps.edp \
	Stokes-v1-matrix-mumps.edp \
	Stokes-v2-matrix-mumps.edp \
	Stokes-v3-matrix-mumps.edp \
	testsolver_MUMPS.edp \
	VG.edp

if ARPACK
CONDITIONAL_MPI_MUMPS_ARPACK = LapMPIEigenValue.edp
endif ARPACK
endif MUMPS
endif MPI

TESTS = $(CONDITIONAL_MPI) \
	$(CONDITIONAL_MPI_HTOOL) \
	$(CONDITIONAL_MPI_METIS) \
	$(CONDITIONAL_MPI_METIS_PARMETIS) \
	$(CONDITIONAL_MPI_MUMPS) \
	$(CONDITIONAL_MPI_MUMPS_ARPACK)

XFAIL_TESTS = $(TESTS)

EXTRA_DIST = *.edp

freefem++.pref: Makefile
	echo "loadpath = ../../plugin/mpi/" > freefem++.pref
	echo "loadpath += ../../plugin/seq/" >>freefem++.pref
	echo "includepath = ../../idp/" >> freefem++.pref

clean-local::
