## Documentation
You can find the PETSc and SLEPc documentation [here](https://doc.freefem.org/documentation/petsc/index.html).

## Examples
### Linear problems
| Filename                                                                                                                                                      | Comments (preconditioners, numerical schemes)                                 |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------|
| [diffusion-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-2d-PETSc.edp)                                       | Distributed LU/Cholesky, domain decomposition and multigrid methods           |
| [diffusion-2d-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-2d-PETSc-complex.edp)                       | &nbsp;                                                                        |
| [heat-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/heat-2d-PETSc.edp)                                                 | Transient diffusion equation, same as above                                   |
| [diffusion-periodic-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-periodic-2d-PETSc.edp)                     | Periodic boundary conditions, multigrid methods                               |
| [diffusion-periodic-balanced-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-periodic-balanced-2d-PETSc.edp)   | Better load balancing than above example                                      |
| [diffusion-substructuring-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-substructuring-2d-PETSc.edp)         | Balancing Domain Decomposition with Constraints                               |
| [diffusion-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-3d-PETSc.edp)                                       | Three-dimensional problem, domain decomposition and multigrid methods         |
| [diffusion-mg-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-mg-3d-PETSc.edp)                                 | Geometric multigrid methods                                                   |
| [helmholtz-2d-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-2d-PETSc-complex.edp)                       | Domain decomposition methods with optimized boundary conditions               |
| [helmholtz-mg-2d-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-mg-2d-PETSc-complex.edp)                 | Geometric multigrid methods                                                   |
| [laplace-RT-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-RT-2d-PETSc.edp)                                     | Vectorial two-dimensional problem with a block preconditioner (fieldsplit)    |
| [laplace-adapt-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-adapt-3d-PETSc.edp)                               | Three-dimensional problem with *h* adaptivity, multigrid methods              |
| [laplace-lagrange-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-lagrange-PETSc.edp)                               | Laplace equation with constraints and a block preconditioner (fieldsplit)     |
| [elasticity-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/elasticity-2d-PETSc.edp)                                     | Vectorial problem, domain decomposition (GenEO) and multigrid methods         |
| [elasticity-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/elasticity-3d-PETSc.edp)                                     | &nbsp;                                                                        |
| [stokes-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-2d-PETSc.edp)                                             | Distributed LU/Cholesky                                                       |
| [stokes-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-3d-PETSc.edp)                                             | &nbsp;                                                                        |
| [stokes-block-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-block-2d-PETSc.edp)                                 | Stokes equation defined as a block system with four matrices (fieldsplit)     |
| [stokes-fieldsplit-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-fieldsplit-2d-PETSc.edp)                       | Block preconditioner (fieldsplit)                                             |
| [stokes-fieldsplit-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-fieldsplit-3d-PETSc.edp)                       | &nbsp;                                                                        |
| [maxwell-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/maxwell-2d-PETSc.edp)                                           | Direct LU/Cholesky                                                            |
| [maxwell-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/maxwell-3d-PETSc.edp)                                           | Multigrid method                                                              |
| [maxwell-mg-3d-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/maxwell-mg-3d-PETSc-complex.edp)                     | Two-grid preconditioner                                                       |
| [helmholtz-3d-surf-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-3d-surf-PETSc-complex.edp)             | BEM with hierarchical matrices from Htool                                     |
| [helmholtz-3d-line-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-3d-line-PETSc-complex.edp)             |                                                                               |

### Nonlinear problems
| Filename                                                                                                                                                          | Comments (preconditioners, numerical schemes)                                 |
|-------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------|
| [bratu-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/bratu-2d-PETSc.edp)                                                   | &nbsp;                                                                        |
| [bratu-hpddm-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/bratu-hpddm-2d-PETSc.edp)                                       | GenEO with reused coarse spaces                                               |
| [newton-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/newton-2d-PETSc.edp)                                                 | &nbsp;                                                                        |
| [newton-adaptmesh-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/newton-adaptmesh-2d-PETSc.edp)                             | Newton method and *h* adaptivity                                              |
| [newton-vi-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/newton-vi-2d-PETSc.edp)                                           | Newton method and a variational inequality                                    |
| [newton-vi-adaptmesh-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/newton-vi-adaptmesh-2d-PETSc.edp)                       | Newton method, *h* adaptivity, and a variational inequality                   |
| [elasticity-SNES-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/elasticity-SNES-3d-PETSc.edp)                               | Linear elasiticty with a Newton method                                        |
| [neo-Hookean-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/neo-Hookean-2d-PETSc.edp)                                       | Nonlinear elasticity                                                          |
| [navier-stokes-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/navier-stokes-2d-PETSc.edp)                                   | Steady-state Navier--Stokes equation for linear stability analysis            |
| [natural-convection-fieldsplit-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/natural-convection-fieldsplit-2d-PETSc.edp)   | Newton method and *h* adaptivity                                              |
| [vi-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/vi-2d-PETSc.edp)                                                         | Variational inequalities                                                      |

### Time steppers and optimizers
| Filename                                                                                                                                      | Comments (preconditioners, numerical schemes) |
|-----------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------|
| [advection-TS-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/advection-TS-2d-PETSc.edp)                 | Implicit and explicit schemes                 |
| [heat-TS-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/heat-TS-2d-PETSc.edp)                           | &nbsp;                                        |
| [heat-TS-RHS-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/heat-TS-RHS-2d-PETSc.edp)                   | &nbsp;                                        |
| [minimal-surface-Tao-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/minimal-surface-Tao-2d-PETSc.edp)   | Minimal surface problem                       |
| [orego-Tao-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/orego-Tao-PETSc.edp)                             | &nbsp;                                        |
| [toy-Tao-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/toy-Tao-PETSc.edp)                                 | &nbsp;                                        |

### Eigenvalue problems
| Filename                                                                                                                                                              | Comments (preconditioners, numerical schemes) |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------|
| [laplace-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-2d-SLEPc.edp)                                                   | &nbsp;                                                |
| [laplace-spherical-harmonics-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-spherical-harmonics-2d-SLEPc.edp)           | &nbsp;                                                |
| [laplace-torus-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-torus-2d-SLEPc.edp)                                       | &nbsp;                                                |
| [schrodinger-axial-well-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-axial-well-2d-SLEPc.edp)                     | &nbsp;                                                |
| [schrodinger-harmonic-oscillator-1d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-harmonic-oscillator-1d-SLEPc.edp)   | &nbsp;                                                |
| [schrodinger-harmonic-oscillator-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-harmonic-oscillator-2d-SLEPc.edp)   | &nbsp;                                                |
| [schrodinger-square-well-1d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-square-well-1d-SLEPc.edp)                   | &nbsp;                                                |
| [laplace-2d-SLEPc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-2d-SLEPc-complex.edp)                                   | &nbsp;                                                |
| [laplace-beltrami-3d-surf-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-beltrami-3d-surf-SLEPc.edp)                       | Eigenvalue problem on a surface                       |
| [laplace-beltrami-3d-line-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-beltrami-3d-line-SLEPc.edp)                       | Eigenvalue problem on a curve                         |
| [navier-stokes-2d-SLEPc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/navier-stokes-2d-SLEPc-complex.edp)                       | Linear stability analysis of Navier--Stokes equations |

### Miscellaneous
| Filename                                                                                                                          | Comments (preconditioners, numerical schemes) |
|-----------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| [transpose-solve-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/transpose-solve-PETSc.edp)     | Solving a transposed system                                         |
| [Schur-complement-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/Schur-complement-PETSc.edp)   | Computing an exact Schur complement                                 |
| [block-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/block-PETSc.edp)                         | &nbsp;                                                              |
| [heat-io-2d.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/heat-io-2d.edp)                           | Automatic ParaView animation output                                 |
| [stokes-io-3d.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-io-3d.edp)                       | &nbsp;                                                              |
| [buildRecursive.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/buildRecursive.edp)                   | Recursive mesh partitioning (for geometric multigrid)               |
| [withPartitioning.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/withPartitioning.edp)               | Connectivity construction with a user-supplied partitioning         |
| [createPartition.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/createPartition.edp)                 | Creation of different partitions of unity using the same DD         |
| [save-load-Dmesh.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/save-load-Dmesh.edp)                 | Saving and loading a distributed mesh for restarting a computation  |
| [transfer.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/transfer.edp)                               | Parallel interpolation of finite element functions                  |

### Reproducible science
| Article                                                                                                                                                       | Source code                                                   |
|---------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------|
| [Augmented Lagrangian Preconditioner for Large-Scale Hydrodynamic Stability Analysis](https://www.sciencedirect.com/science/article/pii/S0045782519301914)    | [GitHub repository](https://github.com/prj-/moulin2019al)     |
| [A Multilevel Schwarz Preconditioner Based on a Hierarchy of Robust Coarse Spaces](https://hal.archives-ouvertes.fr/hal-02151184/document)                    | [GitHub repository](https://github.com/prj-/aldaas2019multi)  |
