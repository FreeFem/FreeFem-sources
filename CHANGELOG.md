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

## [Unreleased] 4.1
### Fixed
-- missing conj operation is some hermitian operation on complex sparse matrix like  A+c*B', A*B' thanks to P-H Tournier
-- writing CheckAllEdp to be compatible with new tree
-- fix eps in trunc in case of very anisotrope mesh, Thank G. Sadaka

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
- Add cosmetics in macro (macro name, macro line, ...)

### Changed
- Pass to PETSc/SLEPc version 3.8.4/3.8.3

### Fixed
- Fix launchff.exe bug under windows 64 to choose a filescrip if no parameter
- Fix the label definition in case of `intalledges` in 2d
- Fix mpi_comm with MUMPS (very rare)

## 3.60 - 2018-04-13
### Changed
- The main distribution is now on Github

[Unreleased]: https://github.com/FreeFem/FreeFem-sources/compare/v4.0...develop
[4.0]: https://github.com/FreeFem/FreeFem-sources/compare/3.62...v4.0
[3.62]: https://github.com/FreeFem/FreeFem-sources/compare/3.61...3.62
[3.61]: https://github.com/FreeFem/FreeFem-sources/compare/v3.60...3.61
