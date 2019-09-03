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
TESTS_ENVIRONMENT = TEST_FFPP=$(TEST_FFPPMPI) FLAGS_FFPP="-np 4 -nw" LIST_CONDITIONAL="$(LIST_CONDITIONAL)"

all-local: freefem++.pref

TESTS_MPI = withPartitioning.edp

TESTS_HPDDM = diffusion-2d.edp \
	diffusion-2d-mg.edp \
	diffusion-2d-substructuring.edp \
	diffusion-2d-substructuring-withPartitioning.edp \
	diffusion-3d.edp \
	diffusion-3d-simple.edp \
	diffusion-periodic-2d.edp \
	elasticity-2d.edp \
	elasticity-2d-substructuring.edp \
	elasticity-3d.edp \
	elasticity-3d-simple.edp \
	elasticity-block.edp \
	heat-2d.edp \
	heat-2d-io.edp \
	heat-3d.edp \
	helmholtz-2d.edp \
	helmholtz-2d-mg.edp \
	iterative.edp \
	maxwell-3d.edp \
	stokes-2d.edp \
	stokes-3d.edp \
	stokes-3d-io.edp

TESTS_HPDDM_PETSC = diffusion-2d-substructuring-PETSc.edp

TESTS_PETSC = bratu-2d-PETSc.edp \
	diffusion-2d-PETSc.edp \
	diffusion-3d-PETSc.edp \
	diffusion-periodic-2d-PETSc.edp \
	diffusion-periodic-balanced-2d-PETSc.edp \
	elasticity-2d-PETSc.edp \
	elasticity-3d-PETSc.edp \
	elasticity-3d-SNES-PETSc.edp \
	heat-2d-PETSc.edp \
	laplace-lagrange-PETSc.edp \
	natural-convection-2d-PETSc-fieldsplit.edp \
	neo-Hookean-2d-PETSc.edp \
	Newton-2d-PETSc.edp \
	Newton-adaptmesh-2d-PETSc.edp \
	Newton-vi-2d-PETSc.edp \
	Newton-vi-adaptmesh-2d-PETSc.edp \
	PETSc-block.edp \
	laplace-RT-2d-PETSc.edp \
	stokes-2d-PETSc.edp \
	stokes-2d-PETSc-fieldsplit.edp \
	stokes-2d-block-PETSc.edp \
	stokes-3d-PETSc.edp \
	stokes-3d-PETSc-fieldsplit.edp \
	transpose-solve-PETSc.edp \
	vi-2d-PETSc.edp \
	orego.edp \
	heat-2d-TS-PETSc.edp \
	heat-2d-TS-RHS-PETSc.edp \
	advection-2d-TS-PETSc.edp \
	toy-Tao-PETSc.edp \
	minimal-surface-2d-Tao-PETSc.edp \
	Schur-complement-PETSc.edp

TESTS_PETSC_SLEPC = laplace-2d-SLEPc.edp \
	laplace-2d-spherical-harmonics-SLEPc.edp \
	laplace-2d-torus-SLEPc.edp \
	schrodinger-1d-harmonic-oscillator-SLEPc.edp \
	schrodinger-1d-square-well-SLEPc.edp \
	schrodinger-2d-axial-well-SLEPc.edp \
	schrodinger-2d-harmonic-oscillator-SLEPc.edp


TESTS_PETSCCOMPLEX = diffusion-2d-PETSc-complex.edp \
	helmholtz-2d-PETSc.edp

TESTS_PETSCCOMPLEX_SLEPCCOMPLEX = laplace-complex-2d-SLEPc.edp

if MPI
CONDITIONAL_MPI = $(TESTS_MPI)
endif MPI

if HPDDM
CONDITIONAL_HPDDM = $(TESTS_HPDDM)
if PETSC
CONDITIONAL_HPDDM_PETSC = $(TESTS_HPDDM_PETSC)
endif PETSC
endif HPDDM


if PETSC
CONDITIONAL_PETSC = $(TESTS_PETSC)
if SLEPC
CONDITIONAL_PETSC_SLEPC = $(TESTS_PETSC_SLEPC)
endif SLEPC
endif PETSC

if PETSCCOMPLEX
CONDITIONAL_PETSCCOMPLEX = $(TESTS_PETSCCOMPLEX)
if SLEPCCOMPLEX
CONDITIONAL_PETSCCOMPLEX_SLEPCCOMPLEX = $(TESTS_PETSCCOMPLEX_SLEPCCOMPLEX)
endif SLEPCCOMPLEX
endif PETSCCOMPLEX

TESTS = $(TESTS_MPI) \
	$(TESTS_HPDDM) \
	$(TESTS_HPDDM_PETSC) \
	$(TESTS_PETSC) \
	$(TESTS_PETSC_SLEPC) \
	$(TESTS_PETSCCOMPLEX) \
	$(TESTS_PETSCCOMPLEX_SLEPCCOMPLEX)

LIST_CONDITIONAL = $(CONDITIONAL_MPI) \
	$(CONDITIONAL_HPDDM) \
	$(CONDITIONAL_HPDDM_PETSC) \
	$(CONDITIONAL_PETSC) \
	$(CONDITIONAL_PETSC_SLEPC) \
	$(CONDITIONAL_PETSCCOMPLEX) \
	$(CONDITIONAL_PETSCCOMPLEX_SLEPCCOMPLEX)

XFAIL_TESTS = $(TESTS)

EXTRA_DIST = *.edp

freefem++.pref: Makefile
	echo loadpath = \"../../plugin/mpi/\" > freefem++.pref
	echo loadpath += \"../../plugin/seq/\" >> freefem++.pref
	echo includepath = \"../../idp/\" >> freefem++.pref

clean-local:
	rm -rf sol_heat_2d_io* sol_stokes_3d_io*
