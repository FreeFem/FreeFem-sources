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

TESTS_MPI = withPartitioning.edp buildRecursive.edp createPartition.edp reconstructDmesh.edp redistributeDmesh.edp convect.edp

TESTS_HPDDM = diffusion-substructuring-2d.edp \
	diffusion-substructuring-withPartitioning-2d.edp \
	elasticity-2d.edp \
	elasticity-substructuring-2d.edp \
	elasticity-block.edp \
	heat-2d.edp \
	heat-io-2d.edp \
	heat-3d.edp \
	helmholtz-2d.edp \
	helmholtz-mg-2d.edp \
	iterative.edp \
	maxwell-3d.edp \
	heat-torus-3d-surf.edp

TESTS_HPDDM_ARPACK = diffusion-2d.edp \
	diffusion-mg-2d.edp \
	diffusion-3d.edp \
	diffusion-simple-3d.edp \
	diffusion-periodic-2d.edp \
	elasticity-3d.edp \
	elasticity-simple-3d.edp \
	stokes-2d.edp \
	stokes-3d.edp \
	stokes-io-3d.edp

TESTS_PETSC = bratu-2d-PETSc.edp \
	diffusion-2d-PETSc.edp \
	diffusion-3d-PETSc.edp \
	diffusion-periodic-2d-PETSc.edp \
	diffusion-periodic-balanced-2d-PETSc.edp \
	elasticity-2d-PETSc.edp \
	elasticity-3d-PETSc.edp \
	elasticity-SNES-3d-PETSc.edp \
	heat-2d-PETSc.edp \
	laplace-lagrange-PETSc.edp \
	natural-convection-fieldsplit-2d-PETSc.edp \
	neo-Hookean-2d-PETSc.edp \
	newton-2d-PETSc.edp \
	newton-adaptmesh-2d-PETSc.edp \
	newton-vi-2d-PETSc.edp \
	newton-vi-adaptmesh-2d-PETSc.edp \
	block-PETSc.edp \
	laplace-RT-2d-PETSc.edp \
	stokes-2d-PETSc.edp \
	stokes-fieldsplit-2d-PETSc.edp \
	stokes-block-2d-PETSc.edp \
	MatLoad-PETSc.edp \
	stokes-3d-PETSc.edp \
	stokes-fieldsplit-3d-PETSc.edp \
	transpose-solve-PETSc.edp \
	bratu-hpddm-2d-PETSc.edp \
	vi-2d-PETSc.edp \
	orego-Tao-PETSc.edp \
	heat-TS-2d-PETSc.edp \
	heat-TS-RHS-2d-PETSc.edp \
	advection-TS-2d-PETSc.edp \
	toy-Tao-PETSc.edp \
	minimal-surface-Tao-2d-PETSc.edp \
	Schur-complement-PETSc.edp \
	maxwell-2d-PETSc.edp \
	maxwell-3d-PETSc.edp \
	laplace-adapt-3d-PETSc.edp \
	diffusion-mg-2d-PETSc.edp \
	diffusion-mg-3d-PETSc.edp \
	save-load-Dmesh.edp \
	navier-stokes-2d-PETSc.edp \
	oseen-2d-PETSc.edp \
	transfer.edp \
	DMPlex-PETSc.edp \
	PtAP-2d-PETSc.edp \
	restriction-2d-PETSc.edp \
	function-PETSc.edp \
	bilaplace-2d-PETSc.edp
# diffusion-substructuring-2d-PETSc.edp <- broken, problem with BDDC

TESTS_PETSC_PARMMG = distributed-parmmg.edp laplace-adapt-dist-3d-PETSc.edp

TESTS_PETSC_SLEPC = laplace-2d-SLEPc.edp \
	laplace-spherical-harmonics-2d-SLEPc.edp \
	laplace-torus-2d-SLEPc.edp \
	schrodinger-harmonic-oscillator-1d-SLEPc.edp \
	schrodinger-square-well-1d-SLEPc.edp \
	schrodinger-axial-well-2d-SLEPc.edp \
	schrodinger-harmonic-oscillator-2d-SLEPc.edp \
	laplace-beltrami-3d-surf-SLEPc.edp \
	laplace-beltrami-3d-line-SLEPc.edp \
	stokes-2d-SLEPc.edp \
	mf-2d-SLEPc.edp

TESTS_PETSCCOMPLEX = diffusion-2d-PETSc-complex.edp \
	helmholtz-2d-PETSc-complex.edp \
	helmholtz-mg-2d-PETSc-complex.edp \
	maxwell-mg-3d-PETSc-complex.edp

TESTS_PETSCCOMPLEX_SLEPCCOMPLEX = laplace-2d-SLEPc-complex.edp \
	navier-stokes-2d-SLEPc-complex.edp helmholtz-2d-SLEPc-complex.edp \
	nonlinear-2d-SLEPc-complex.edp

TESTS_PETSCCOMPLEX_BEMTOOL_BOOST_HTOOL = helmholtz-3d-surf-PETSc-complex.edp helmholtz-3d-line-PETSc-complex.edp helmholtz-coupled-2d-PETSc-complex.edp helmholtz-dense-3d-line-PETSc-complex.edp

if MPI
CONDITIONAL_MPI = $(TESTS_MPI)

if HPDDM
CONDITIONAL_HPDDM = $(TESTS_HPDDM)
if ARPACK
CONDITIONAL_HPDDM_ARPACK = $(TESTS_HPDDM_ARPACK)
endif ARPACK
endif HPDDM

if PETSC
CONDITIONAL_PETSC = $(TESTS_PETSC)
if PARMMG
CONDITIONAL_PETSC_PARMMG = $(TESTS_PETSC_PARMMG)
endif PARMMG
if SLEPC
CONDITIONAL_PETSC_SLEPC = $(TESTS_PETSC_SLEPC)
endif SLEPC
endif PETSC

if PETSCCOMPLEX
CONDITIONAL_PETSCCOMPLEX = $(TESTS_PETSCCOMPLEX)
if SLEPCCOMPLEX
CONDITIONAL_PETSCCOMPLEX_SLEPCCOMPLEX = $(TESTS_PETSCCOMPLEX_SLEPCCOMPLEX)
endif SLEPCCOMPLEX
if BEMTOOL
if BOOST
if HTOOL
CONDITIONAL_PETSCCOMPLEX_BEMTOOL_BOOST_HTOOL = $(TESTS_PETSCCOMPLEX_BEMTOOL_BOOST_HTOOL)
endif HTOOL
endif BOOST
endif BEMTOOL
endif PETSCCOMPLEX

endif MPI

TESTS = $(TESTS_MPI) \
	$(TESTS_HPDDM) \
	$(TESTS_HPDDM_ARPACK) \
	$(TESTS_PETSC) \
	$(TESTS_PETSC_PARMMG) \
	$(TESTS_PETSC_SLEPC) \
	$(TESTS_PETSCCOMPLEX) \
	$(TESTS_PETSCCOMPLEX_SLEPCCOMPLEX) \
	$(TESTS_PETSCCOMPLEX_BEMTOOL_BOOST_HTOOL)

LIST_CONDITIONAL = $(CONDITIONAL_MPI) \
	$(CONDITIONAL_HPDDM) \
	$(CONDITIONAL_HPDDM_ARPACK) \
	$(CONDITIONAL_PETSC) \
	$(CONDITIONAL_PETSC_PARMMG) \
	$(CONDITIONAL_PETSC_SLEPC) \
	$(CONDITIONAL_PETSCCOMPLEX) \
	$(CONDITIONAL_PETSCCOMPLEX_SLEPCCOMPLEX) \
	$(CONDITIONAL_PETSCCOMPLEX_BEMTOOL_BOOST_HTOOL)

XFAIL_TESTS = $(TESTS)

EXTRA_DIST = *.edp cube.msh rectangle.msh

freefem++.pref: Makefile
	echo loadpath = \"../../plugin/mpi/\" > freefem++.pref
	echo loadpath += \"../../plugin/seq/\" >> freefem++.pref
	echo includepath = \"../../idp/\" >> freefem++.pref

clean-local:
	rm -rf heat-io-2d_* stokes-io-3d_* laplace-adapt-3d_* diffusion-mg-3d_* maxwell-mg-3d_* createPartition_* helmholtz-3d-surf.vtu helmholtz-3d-line.vtu dump-2d_* dump-3d_* navier-stokes-2d_* laplace-adapt-dist-3d_* transfer_2d_* transfer_3d_* stokes-mat-* stokes-rhs-* oseen-2d_*
