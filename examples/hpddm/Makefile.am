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
# SUMMARY : Makefile for FreeFem++
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

LOG_DRIVER=$(SHELL) $(top_srcdir)/bin/test-driver-ff
TESTS_ENVIRONMENT=TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" SKIP=$(SKIP_TESTS_MPI)

all-local: freefem++.pref

TESTS=diffusion-2d-PETSc.edp diffusion-2d-PETSc-complex.edp diffusion-3d-PETSc.edp \
	elasticity-2d-PETSc.edp elasticity-3d-PETSc.edp \
	stokes-2d-PETSc.edp stokes-3d-PETSc.edp \
	heat-2d-PETSc.edp diffusion-periodic-2d-PETSc.edp diffusion-periodic-2d.edp diffusion-periodic-balanced-2d-PETSc.edp helmholtz-2d-PETSc.edp heat-2d-io.edp \
	laplace-2d-SLEPc.edp laplace-complex-2d-SLEPc.edp stokes-2d-PETSc-fieldsplit.edp stokes-3d-io.edp \
	laplace-2d-spherical-harmonics-SLEPc.edp schrodinger-1d-harmonic-oscillator-SLEPc.edp schrodinger-2d-harmonic-oscillator-SLEPc.edp \
	laplace-2d-torus-SLEPc.edp schrodinger-1d-square-well-SLEPc.edp transpose-solve-PETSc.edp stokes-3d-PETSc-fieldsplit.edp \
	PETSc-block.edp diffusion-2d.edp diffusion-3d.edp diffusion-3d-simple.edp schrodinger-2d-axial-well-SLEPc.edp \
	diffusion-2d-substructuring.edp diffusion-2d-substructuring-PETSc.edp \
	elasticity-2d.edp Newton-adaptmesh-2d-PETSc.edp Newton-vi-adaptmesh-2d-PETSc.edp \
	stokes-3d.edp stokes-2d.edp \
	heat-2d.edp heat-3d.edp Newton-adaptmesh-2d-PETSc.edp \
	helmholtz-2d.edp maxwell-3d.edp \
	iterative.edp withPartitioning.edp \
	elasticity-2d-substructuring.edp elasticity-block.edp \
	diffusion-2d-substructuring-withPartitioning.edp \
	elasticity-3d.edp elasticity-3d-simple.edp laplace-lagrange-PETSc.edp \
	elasticity-3d-SNES-PETSc.edp vi-2d-PETSc.edp Newton-vi-2d-PETSc.edp Newton-2d-PETSc.edp neo-Hookean-2d-PETSc.edp bratu-2d-PETSc.edp \
	natural-convection-2d-PETSc-fieldsplit.edp

XFAIL_TESTS=$(TESTS)

EXTRA_DIST=*.edp

freefem++.pref:Makefile
	echo loadpath = \"../../plugin/mpi/\" > freefem++.pref
	echo loadpath += \"../../plugin/seq/\" >> freefem++.pref
	echo includepath = \"../../idp/\" >> freefem++.pref

clean-local:
	rm -rf sol_heat_2d_io* sol_stokes_3d_io*
