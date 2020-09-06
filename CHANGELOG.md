<!--USE THIS TEMPLATE TO COMPLETE THE CHANGELOG-->
<!--
## [Version number] - YYYY-MM-DD
### Added
-

### Changed
-

### Deprecated
-

### Removed
-

### Fixed
-

### Security
-
-->

# Changelog
All notable changes to this project will be documented in this file.

## [Unreleased]
### Added

### Changed
- change the language definition to use tyep as a construction function with named arguments for bem plugin
- PETSc version 3.14.0
- ARPACK compiled by SLEPc

### Deprecated

### Removed
- plugins thresholdings, symmetrizeCSR, and fflapack and associed example

### Fixed
- problem compilation with gfortran-10 of arpack and mumps (add -fallow-argument-mismatch flags- 

## [4.7]
### Added

- new way to build matrix beetween 2d Finite element 2d and Curve finite element to do mortar (Thank to Axel ) , see first example `examples/tutorial/mortar-DN-4-v4.5.edp`
- add `Ns` normal vector  in R^3 on meshS (normal of the surface) of current point (to day Ns of [x,y,0] plan  is [0,0,-1])  no be compatibe to exterior normal. 
- add `Tl` tangent vector in R^3 on meshL (tangent vector of the line/curve) of current point 
- compile ffmaster / ffslave example under windows (thanks to johann@ifado.de)
- Boolean parameter `spiltpbedge` in `buildmesh` to split in to edge with two boundary vertices
- interface to PETSc DMPlex, see `examples/hpddm/DMPlex-PETSc.edp`
- function `MatDestroy`
- function `MatPtAP` and `transferMat` for parallel interpolation between non-matching grids, see `examples/hpddm/PtAP-2d-PETSc.edp` or `examples/hpddm/diffusion-mg-2d-PETSc.edp`
- preliminary interface to `SVDSolve` from SLEPc to compute singular value decompositions, see `examples/hpddm/mf-2d-SLEPc.edp` or `examples/hpddm/helmholtz-2d-SLEPc-complex.edp`
- preliminary interface to `NEPSolve` from SLEPc to solve nonlinear eigenvalue problems, see `examples/hpddm/nonlinear-2d-SLEPc-complex.edp`
- `transpose` parameter when constructing a `Mat` for defining a matrix-free transposed operation
- interface to `PetscMemoryGetCurrentUsage`
- add P2b, RT0, RT1 surface FE (P2bS, RT0S, RT1S))
- add operator interpolate (2d->3d surface)
- add operator x = A'\*b; where x, b are array and A 2 dim array (full matrix) and generate an error in case of b'\*A or b'\*A expression
- function `MatLoad` to load a PETSc `Mat` from disk, see `examples/hpddm/MatLoad-PETSc.edp`
- possibility to assemble a symmetric `HMatrix<complex>` and to densify a `HMatrix<complex>` into a `Mat<complex>`

### Changed
- moved Htool to its new GitHub location
- ScaLAPACK and MUMPS are not compiled by PETSc anymore if there is no Fortran compiler
- MPICH is compiled by PETSc if no MPI is detected during configure, see https://community.freefem.org/t/feature-request-use-download-mpich-on-ubuntu/407
- PETSc version 3.13.5
- force `--with-cudac=0` in `make petsc-slepc`, see https://github.com/FreeFem/FreeFem-sources/issues/141
- change DSL keyword P1dc3dL->P1dcL and P1dc3dS->P1dcS
- rename `view`, `hasType`, `changeSchur` to respectively `ObjectView`, `HasType`, and `ChangeSchur`

### Deprecated
- rename `changeNumbering`, `globalNumbering`, `originalNumbering`, `changeOperator`, `destroyRecycling`, and `attachCoarseOperator` to respectively `ChangeNumbering`, `GlobalNumbering`, `OriginalNumbering`, `ChangeOperator`, `DestroyRecycling`, and `AttachCoarseOperator`
- `Nt` the normal vector of the current (wrong on meshL) use `Ns` pr `Tl` 
### Removed
- `augmentation` routine from the PETSc plugin
- `MPIF77` variable

### Fixed
- lot of mistake in MeshL element add a example o ckeck lot of thing `tutomesh1d.edp`
- fixed problem of change of mesh when rebuild 2d mesh with builmesh, .... (Thank to P. Jovilet to points this problem)
- missing METIS library when using SuiteSparse compiled by PETSc
- missing `-fno-stack-protector` when building PETSc on Windows, see https://community.freefem.org/t/error-loading-complex-petsc-slepc-library/370
- fixed ffglut for the plotting of FE array solution
- fixed  ffglut bug on MacOS Catalina , draw inn only half windows screen (Apple Bug ???)
- correct P0VF  finite element 
- `abs` function of array

## [4.6]

### Added
- new algorithm for searching element contening a point (more safe) in mesh of type mesh3, meshS, or meshL.
- new function `hasType` to know if a PETSc component has been installed, e.g., `hasType("PC", "hypre")`
- eigenvalue problems on linear elements, cf. `examples/eigen/LapEigen1DBeltrami.edp` or `examples/hpddm/laplace-beltrami-3d-line-SLEPc.edp`
- `--download-cmake` in PETSc configure if there is no CMake available
- flags `--with-[slepc|slepccomplex]-include` and `--with-[slepc|slepccomplex]-ldflags` for when SLEPc has been built outside of FreeFEM or PETSc
- interface to `KSPSetResidualHistory` and `KSPGetIterationNumber`
- interface to `mpiWaitAll`
- new function extract, allows to build a curve mesh from a 2d mesh (can extract a labeled boundary, apply a geometric transformation)
- ffglut can plot a vectorial FE function in surface 3d
- distributed ParMmg interface, cf. `examples/hpddm/distributed-parmmg.edp` or `examples/hpddm/laplace-adapt-dist-3d-PETSc.edp`
- new parallel interpolator on non-matching meshes, cf. `examples/hpddm/transfer.edp`
- ability to solve problems in single precision or with 64 bit integers
- tool to read data form vtk file only in 3d (cf. plugin iovtk a first example `examples/plugin/iovtk.edp`)
- tool to read/wrile ply file of meshL, mesh3, MeshS : Polygon File Format / Stanford Triangle Format do  `load "ioply"`
     see `examples//3dSurf/operatorsOnMeshS.edp`

### Changed
- new `tgv` values: -10 => zero row, -20 => zero row/column
- Windows binary now shipped with PETSc/SLEPc
- BEM examples are now in `examples/mpi`
- plot border type is now in 3d (border 2d and 3d)
- PETSc version 3.13.0

### Deprecated

### Fixed
- `--enable-download_package` may now be used to download a single package, e.g., `--enable-download_metis`
- compilation of PETSc under Windows
- compilation of plugins when using static libraries
- correct detection problem in FE type when use a vectorial FE
- macro concatenation with spaces in arguments
- correct bug in `plugin/seq/Schur-Complement.cpp`
- correct ambigity bug in `plugin/seq/bfstream.cpp` (reading real or integer)
- compilation fo  plugin libff-mmap-semaphore.c under windows

## [4.5]

### Added
- for windows version: rename under mpi `MUMPS` in `MUMPS_mpi` and in sequentiel in `MUMPS_seq` due to confict between seq. and mpi version so all MUMPS load become `MUMPS_seq` or `MUMPS_mpi`in all examples
- correct link edition with fortran mpi under windows juste use the msmpi (just use `libmsmpi.dll` )
- new `mmg` and `parmmg` (parallel mmg) plugins interfacing mmg5 and parmmg libraries, to replace `mmg3d-v4.0` and `freeyams` (Thanks to P-H Tournier)
- a true 3d anisotropic mesh adaptation `examples/3d/Laplace-Adapt-aniso-3d.edp`
- an example to extract surface mesh from isovalue in `examples/3dSurf/Pinochio.edp`
- function `f.eatspace` to reach eof on istream file which return false in case of EOF.
- function `f.length` to get the istream file length
- Interface to `PetscLogStagePush()`/`PetscLogStagePop()`
- Ability to directly assemble a `Mat` using a `varf`
- New `bem` plugin for the Boundary Element Method (using htool and BemTool libraries)
- New DSL for BEM (varfbem see examples/bem)
- add int0d to apply Neumann BC (curve FE), differential operators (dx,dy,...), compute an 1d integral
- add P1dc FE for Border FEM (possible to define a new FE with plugin)
- PETSc as a subdomain solver for HPDDM

### Changed
- correct ffglut (bug in case of changing number of nb isovalue)
- PETSc version 3.12.4
- Change the point search triangle algorithm to be sure in any case (in test)
- Sline operator renamed to segment
- In square3, segment, movemesh functions: geometry transformation can now be [X] or [X,Y] or [X,Y,Z] according to the minimal shape element dim
- PETSc now download OpenBLAS if there is no BLAS found by FreeFEM configure

### Deprecated
- freeyams plugin
- mmg3d-v4.0 plugin

### Fixed
- fix plot for curve mesh


## [4.4-3]
### Added
- Preliminary support for symmetric distributed PETSc matrices (MATMPISBAIJ instead of MATMPIAIJ)
- Interface to AMS, Hiptmair--Xu preconditioner for problems in H(curl), see maxwell-3d-PETSc.edp
- FEM on curve 3D (in test)
- P0, P1, P2 curve 3D FE (scalar for the moment)
- i/o medit and vtk format for curve FE
- checkMesh() function, allow to remove multiple vertices, elements and border elements (argument: precisvertice(double),removeduplicate(bool))
- possible to build a curve mesh from a surface, ThS = buildBdMesh(ThS) and define this new mesh by meshL ThL= ThS.Gamma
- can extract a border part of a meshL (meshL ThL = extract(ThL,label=llabs))
- Support for optimized boundary conditions with PETSc, see helmholtz-2d-PETSc-complex.edp
- buildmeshL() function: build meshL from borders
- `mpiCommSelf` keyword

### Changed
- function buildSurface(...) renamed by buildBdMesh(...)
- line3(...) renamed by SLine(...)

### Removed
- FFTW is not compiled by PETSc anymore
- Spurious outputs in TetGen plugin
- curve3 type -> border
- hypre examples since it is not downloaded by FreeFEM for many months (use PETSc instead)
- `dscalprod` routine from HPDDM and PETSc plugins, use `A(u, v)` with `A` a `Mat` or a `schwarz` object
- `export` function for `macro_ddm.idp`, use `savevtk` as in the sequential iovtk plugin

### Fixed
- plotMPI function for plotting 3D solutions, problem with serialize
- variable mes in clean_mesh function
- correct bug verflow in plugin iohdf5
- correct problem with buffer iostrean function (buffer must be out of range )
- correct i/o vtk and by defaut write at binary format
- fix an overflow in RT13d FE
- problem with auto-build of border mesh
## [4.4-2]
### Added
- add matrix and array tools (FH)
   ```
    matrix A=eye(10);
	real[int,int] af = eye(10,10);
    real[int,int] a(10,10);
    int[int] I=[1,3,6];
    real[int] d = a.diag ; // get the diag of full matrix (no copy)
	real[int] dI= d(I); // init a array from renumbering array
    real[int] c= a(:,1)(I); // init a array from renumbering array
    real[int] aa= a.asarray; //  view full the matrice  as an array (no copy)
	a(2:5,3:7).diag= 200;
	a.diag += 100;
	```
- adding of a global variable `lockOrientation` to allows the building of mesh without checking the orientation elements (AF)
- add plugin tool to build matrix edge/P1 with sign `mat_edgeP1`	(FH)
- new examples `diffusion-2d-mg.edp` and `helmholtz-2d-mg.edp` showing how to use user-defined coarse corrections
- support for nonzero scalars in PETSc block matrices
- simpler constructor for sequential HPDDM matrices (no need for the restriction array and the partition of unity)
- array of `Mat` and `schwarz` types
- add mpi meshS (serialize object)

### Changed
- correct mistake in mpirank in case of broadcast with comm (thank tp PHT)
- update fftw to v3.3.8 and openblas v0.3.6
- in movemesh23 correct the argument label -> region to change label
- new implementation for the moving mesh functions, new arguments: boolean cleanmesh, removemultiple, rebuildborder
- new PETSc version 3.12
- templatize movemesh, setMesh functions
- add conditional tests in make check

### Fixed
- spurious output in PARDISO
- fix problem in ffglut (AF)
- detect hdf5 and gsl if no enable-download

### Security

## [4.4]
### Added
- interface to `TSSolve`, DAE/ODE solvers from PETSc
- interface to `TaoSolve`, Toolkit for Advance Optimization from PETSc
- simpler constructor for sequential PETSc matrices (no need for the restriction array and the partition of unity)
- some unit tests

### Changed
- PETSc version 3.11.3
- replaced custom implementations (`RNM::real`, `RNM::norm2`, and `Fem2D::norm`) by C++11 functions
- API of the macro `plotMPI`
- switched to inexact coarse operators in HPDDM by default
- RHS and solution vectors permuted in `IterativeMethod` and `DDM`
- `.mesh` are now saved using version 2 (which stores floating-point scalars in double precision)

### Removed
- legacy linear solver interfaces using the old matrix type
- dot products using CBLAS because of errors at link time
- Newtow function (bad name)

### Fixed
- assertion failure with some 3D meshes when doing `trunc(Th, true)` (thanks to F. Feppon)
- compile error when plotting arrays of vectorial functions

## [4.2.1]
### Added
- nested fieldsplit example `examples/hpddm/natural-convection-2d-PETSc-fieldsplit.edp`
- `int[int][int] array;` is now supported (a size was previously needed, i.e., `array(0);`)
- check selectivity during `make check`, depending on available 3rd party librairies
- new CI/CD tools for `develop` branch
- new gestion of mesh3 - meshS coupling
- square3, buildSurface... operators for meshS

### Changed
- SLEPc is now directly downloaded by PETSc with `--download-slepc`
- HPDDM and PETSc API have been simplified, instead of an `int[int]` and an `int[int][int]`, only a single `int[int][int]` is now needed
- `build` macros for HPDDM and PETSc have been simplified to follow the above API change, two parameters have been permuted as well to match the HPDDM and PETSc constructors
- PETSc version 3.11.2 and HPDDM with multilevel GenEO

### Removed
- old interfaces that were not maintained anymore (pARMS, PaStiX, hips) and that are available through PETSc
- spurious outputs when destroying some meshes
- old surface msh3 type, replaced by meshS

### Fixed
- multiple segmentation faults when using unitialized values (thanks to G. Sadaka)
- nested fieldsplits in the PETSc interface
- memory leaks in `SNESSolve` (nonlinear PETSc solvers)
- bug fix of `Cofactor` function
- various bug fixes on surface mesh

## [4.1]
### Fixed
- missing conj operation is some hermitian operation on complex sparse matrix like  A+c*B', A*B' thanks to P-H Tournier
- writing CheckAllEdp to be compatible with new tree
- fix eps in trunc in case of very anisotrope mesh, Thank G. Sadaka

### Added
- CMake, thanks to [https://github.com/cdoucet](https://github.com/cdoucet)
- Surface finite element, thanks to [AFourmont](https://github.com/AFourmont)
- AppImage generation, thanks to [Alexander Sashnov](https://github.com/asashnov)

### Changed
- PETSc/SLEPc version 3.11

## [4.0]
- correct bug in RT1Ortho and RT2Ortho  2d in the computation of derivative (2018-01-30, Thank to Bryan.Bosworth@colorado.edu)
- uniformize 2d/3d in element, method   EdgeOrientation(e)  now return +1/-1
- change all the sparse matrix structure
- remove all map matrix jan 2019

### Added
- surface finite element (in progress)
- nElementonB (version 2d and 3d of nTonEge)
- area  ( same the lenEdge) in 3d
- add function labels to get the array of label of a  mesh
- add function regions to get the array of label of a  mesh
- correct big bug in 	toRarray,toZarray, toCarray  transform [ ... ] array to int[int], real[int], complex[int]  

## [3.62] - 2018-08-31
### Added
- add  x0=true/false, add veps=eps in solver parameters to initialazed of not the the CG , GMRES algo
  with 0 or previous value and veps is to get the absolue tolerance
- A tool of solve adjoint matrix A with only one single LU decomposition with LU, UMFPACK, GMRES  
    `u[]=A'^-1*b;`  
- Add plugin to save matrix in Harwell-Boeing format (see Harwell-Boeing format)

### Fixed
- Fix bug in `trunc` (2d) in case of very fine mesh (eps too small)

## [3.61] - 2018-06-20
### Added
- Add name parameter `kerneln=`, `kernelt=`, `kerneldim=` for dissection solver
- Add option in method `toClose` function in `fquatree` to get the nearst point (for intersect meshes)
- Add missing file `curvature.edp`
- Add `imax`, `jmax`, `imin`, `jmin` to get index of row or column of the min or max coefficient<br/>
  We have: `A(A.imin,A.jmin) = A.min`
- Add cosmetics in macro (macro name, macro line...)

### Changed
- Pass to PETSc/SLEPc version 3.8.4/3.8.3

### Fixed
- Fix launchff.exe bug under windows 64 to choose a filescrip if no parameter
- Fix the label definition in case of `intalledges` in 2d
- Fix mpi_comm with MUMPS (very rare)

## 3.60 - 2018-04-13
### Changed
- The main distribution is now on Github

[Unreleased]: https://github.com/FreeFem/FreeFem-sources/compare/v4.7..develop
[4.7]: https://github.com/FreeFem/FreeFem-sources/compare/v4.6...v4.7
[4.6]: https://github.com/FreeFem/FreeFem-sources/compare/v4.5...v4.6
[4.5]: https://github.com/FreeFem/FreeFem-sources/compare/v4.4-3...v4.5
[4.4-3]: https://github.com/FreeFem/FreeFem-sources/compare/v4.4-2...v4.4-3
[4.4-2]: https://github.com/FreeFem/FreeFem-sources/compare/v4.4...v4.4-2
[4.4]: https://github.com/FreeFem/FreeFem-sources/compare/v4.2.1...v4.4
[4.2.1]: https://github.com/FreeFem/FreeFem-sources/compare/v4.0...v4.2.1
[4.1]: https://github.com/FreeFem/FreeFem-sources/compare/v4.0...v4.1
[4.0]: https://github.com/FreeFem/FreeFem-sources/compare/3.62...v4.0
[3.62]: https://github.com/FreeFem/FreeFem-sources/compare/3.61...3.62
[3.61]: https://github.com/FreeFem/FreeFem-sources/compare/v3.60...3.61
