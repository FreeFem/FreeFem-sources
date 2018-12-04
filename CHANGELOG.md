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
## [Version 4] 2018-12-01

- change all the sparse amtrice structure 

### Added 
- surface finite elemet (in progress)	

- nElementonB (version 2d and 3d of nTonEge)
- area  ( same the lenEdge) in 3d 

## [Unreleased] 2018-08-31

### Added

- add  x0=true/false, add veps=eps in solver parameters to initiaze of not the the CG , GMRES algo 
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

[Unreleased]: https://github.com/FreeFem/FreeFem-sources/compare/3.61...master
[3.61]: https://github.com/FreeFem/FreeFem-sources/compare/v3.60...3.61
