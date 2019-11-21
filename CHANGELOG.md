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
- Preliminary support for symmetric distributed PETSc matrices (MATMPISBAIJ instead of MATMPIAIJ)
- Interface to AMS, Hiptmair--Xu preconditioner for problemis in H(curl), see maxwell-3d-PETSc.edp

### Changed
-

### Removed
- FFTW is not compiled by PETSc anymore
- Spurious outputs in TetGen plugin

### Fixed
-  plotMPI function for plotting 3D solutions, problem with serialize
-  variable mes in clean_mesh function
-  correct bug verflow in plugin iohdf5
- correct problem with buffer iostrean function (buffer must be out of range )

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

[Unreleased]: https://github.com/FreeFem/FreeFem-sources/compare/v4.4-2...develop
[4.4-2]: https://github.com/FreeFem/FreeFem-sources/compare/v4.4...v4.4-2
[4.4]: https://github.com/FreeFem/FreeFem-sources/compare/v4.2.1...v4.4
[4.2.1]: https://github.com/FreeFem/FreeFem-sources/compare/v4.0...v4.2.1
[4.1]: https://github.com/FreeFem/FreeFem-sources/compare/v4.0...v4.1
[4.0]: https://github.com/FreeFem/FreeFem-sources/compare/3.62...v4.0
[3.62]: https://github.com/FreeFem/FreeFem-sources/compare/3.61...3.62
[3.61]: https://github.com/FreeFem/FreeFem-sources/compare/v3.60...3.61
