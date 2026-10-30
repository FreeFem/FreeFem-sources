## Documentation
You can find the PETSc and SLEPc documentation [here](https://doc.freefem.org/documentation/petsc/index.html).

## Examples
### Linear problems
| Filename                                                                                                                                                    | Comments (preconditioners, numerical schemes)                                 |
|---------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------|
| [diffusion-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-2d-PETSc.edp)                                 | Distributed LU/Cholesky, domain decomposition and multigrid methods               |
| [diffusion-2d-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-2d-PETSc-complex.edp)                 | &nbsp;                                                                            |
| [heat-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/heat-2d-PETSc.edp)                                           | Transient diffusion equation, same as above                                       |
| [diffusion-periodic-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-periodic-2d-PETSc.edp)               | Periodic boundary conditions, multigrid methods                                   |
| [diffusion-substructuring-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-substructuring-2d-PETSc.edp)   | Balancing Domain Decomposition with Constraints                                   |
| [diffusion-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-3d-PETSc.edp)                                 | Three-dimensional problem, domain decomposition and multigrid methods             |
| [diffusion-mg-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-mg-2d-PETSc.edp)                           | Geometric non-nested multigrid methods                                            |
| [diffusion-mg-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/diffusion-mg-3d-PETSc.edp)                           | Geometric nested multigrid methods                                                |
| [helmholtz-2d-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-2d-PETSc-complex.edp)                 | Domain decomposition methods with optimized boundary conditions                   |
| [helmholtz-mg-2d-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-mg-2d-PETSc-complex.edp)           | Geometric multigrid methods                                                       |
| [laplace-RT-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-RT-2d-PETSc.edp)                               | Vectorial two-dimensional problem with a block preconditioner (fieldsplit)        |
| [laplace-RT-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-RT-3d-PETSc.edp)                               | Vectorial three-dimensional problem with a block preconditioner (fieldsplit)      |
| [laplace-adapt-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-adapt-3d-PETSc.edp)                         | Three-dimensional problem with *h* adaptivity, multigrid methods using Mmg        |
| [laplace-adapt-dist-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-adapt-dist-3d-PETSc.edp)               | Three-dimensional problem with fully-distributed *h* adaptivity using ParMmg      |
| [laplace-lagrange-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-lagrange-PETSc.edp)                         | Laplace equation with constraints and a block preconditioner (fieldsplit)         |
| [elasticity-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/elasticity-2d-PETSc.edp)                               | Vectorial problem, domain decomposition (GenEO) and multigrid methods             |
| [elasticity-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/elasticity-3d-PETSc.edp)                               | &nbsp;                                                                            |
| [elasticity-block-hpddm-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/elasticity-block-hpddm-2d-PETSc.edp)       | Domain decomposition (GenEO for saddle-point problems)                            |
| [stokes-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-2d-PETSc.edp)                                       | Distributed LU/Cholesky                                                           |
| [stokes-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-3d-PETSc.edp)                                       | &nbsp;                                                                            |
| [stokes-block-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-block-2d-PETSc.edp)                           | Stokes equation defined as a block system with four matrices (fieldsplit)         |
| [stokes-fieldsplit-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-fieldsplit-2d-PETSc.edp)                 | Block preconditioner (fieldsplit)                                                 |
| [stokes-fieldsplit-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-fieldsplit-3d-PETSc.edp)                 | &nbsp;                                                                            |
| [stokes-block-hpddm-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-block-hpddm-2d-PETSc.edp)               | Domain decomposition (GenEO for saddle-point problems)                            |
| [oseen-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/oseen-2d-PETSc.edp)                                         | Oseen problem preconditioned by Pressure Convection--Diffusion (PCD)              |
| [maxwell-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/maxwell-2d-PETSc.edp)                                     | Direct LU/Cholesky                                                                |
| [maxwell-3d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/maxwell-3d-PETSc.edp)                                     | Multigrid method                                                                  |
| [maxwell-mg-3d-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/maxwell-mg-3d-PETSc-complex.edp)               | Two-grid preconditioner                                                           |
| [helmholtz-3d-surf-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-3d-surf-PETSc-complex.edp)       | BEM with hierarchical matrices from Htool on surfaces                             |
| [helmholtz-3d-line-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-3d-line-PETSc-complex.edp)       | BEM with hierarchical matrices from Htool on lines                                |
| [helmholtz-coupled-2d-PETSc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-coupled-2d-PETSc-complex.edp) | FEM coupled with BEM in a MatNest using a block preconditioner (fieldsplit)       |
| [PtAP-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/PtAP-2d-PETSc.edp)                                           | Parallel interpolation on non-matching grids + P' A P operation                   |
| [restriction-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/restriction-2d-PETSc.edp)                             | Coupled problems with one being defined on a restriction of the mesh of the other |
| [bilaplace-2d-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/bilaplace-2d-PETSc.edp)                                 | Bilaplacian using Morley finite element                                           |

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
| [orego-TS-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/orego-TS-PETSc.edp)                               | &nbsp;                                        |
| [toy-Tao-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/toy-Tao-PETSc.edp)                                 | &nbsp;                                        |

### Eigenvalue problems
| Filename                                                                                                                                                              | Comments (preconditioners, numerical schemes) |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| [laplace-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-2d-SLEPc.edp)                                                   | &nbsp;                                                                                          |
| [laplace-spherical-harmonics-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-spherical-harmonics-2d-SLEPc.edp)           | &nbsp;                                                                                          |
| [laplace-torus-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-torus-2d-SLEPc.edp)                                       | &nbsp;                                                                                          |
| [schrodinger-axial-well-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-axial-well-2d-SLEPc.edp)                     | &nbsp;                                                                                          |
| [schrodinger-harmonic-oscillator-1d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-harmonic-oscillator-1d-SLEPc.edp)   | &nbsp;                                                                                          |
| [schrodinger-harmonic-oscillator-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-harmonic-oscillator-2d-SLEPc.edp)   | &nbsp;                                                                                          |
| [schrodinger-square-well-1d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/schrodinger-square-well-1d-SLEPc.edp)                   | &nbsp;                                                                                          |
| [laplace-2d-SLEPc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-2d-SLEPc-complex.edp)                                   | &nbsp;                                                                                          |
| [laplace-beltrami-3d-surf-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-beltrami-3d-surf-SLEPc.edp)                       | Eigenvalue problem on a surface                                                                 |
| [laplace-beltrami-3d-line-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/laplace-beltrami-3d-line-SLEPc.edp)                       | Eigenvalue problem on a curve                                                                   |
| [navier-stokes-2d-SLEPc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/navier-stokes-2d-SLEPc-complex.edp)                       | Linear stability analysis of Navier--Stokes equations                                           |
| [stokes-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/stokes-2d-SLEPc.edp)                                                     | Inf-sup constant of Stokes equations                                                            |
| [helmholtz-2d-SLEPc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/helmholtz-2d-SLEPc-complex.edp)                               | SVD of the discretized Helmholtz equation                                                       |
| [mf-2d-SLEPc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/mf-2d-SLEPc.edp)                                                              | Matrix-free SVD with user-provided matrix--vector and matrix transpose--vector multiplications |
| [nonlinear-2d-SLEPc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/nonlinear-2d-SLEPc-complex.edp)                                | Nonlinear eigenvalue problem Ax = sqrt(lambda)x                                                |
| [blasius-stability-1d-SLEPc-complex.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/blasius-stability-1d-SLEPc-complex.edp)                | Polynomial eigenvalue problem (A2 lambda^2 x + A1 lambda + A0)x = 0                            |

### Miscellaneous
| Filename                                                                                                                          | Comments (preconditioners, numerical schemes) |
|-----------------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| [transpose-solve-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/transpose-solve-PETSc.edp)     | Solving a transposed system                                              |
| [Schur-complement-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/Schur-complement-PETSc.edp)   | Computing an exact Schur complement                                      |
| [block-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/block-PETSc.edp)                         | &nbsp;                                                                   |
| [buildRecursive.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/buildRecursive.edp)                   | Recursive mesh partitioning (for geometric multigrid)                    |
| [withPartitioning.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/withPartitioning.edp)               | Connectivity construction with a user-supplied partitioning              |
| [PartitionCreate.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/PartitionCreate.edp)                 | Creation of different partitions of unity using the same DD              |
| [Dmesh-Save-Load.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/Dmesh-Save-Load.edp)                 | Saving and loading a distributed mesh for restarting a computation       |
| [transfer.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/transfer.edp)                               | Parallel interpolation of finite element functions                       |
| [DmeshReconstruct.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/DmeshReconstruct.edp)               | Construction of a distributed mesh from sequential nonoverlapping meshes |
| [distributed-parmmg.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/distributed-parmmg.edp)           | Distributed ParMmg interface for mesh adaptation                         |
| [DmeshRedistribute.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/DmeshRedistribute.edp)             | Gather and scatter a Dmesh on different MPI communicators                |
| [DMPlex-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/DMPlex-PETSc.edp)                       | Load and partition a mesh using a DM in two and three dimensions         |
| [function-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/function-PETSc.edp)                   | Custom plugin showing how to directly access PETSc in C++                |
| [convect.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/convect.edp)                                 | Calling the convert FreeFEM function on a distributed mesh               |
| [MatLoad-PETSc.edp](https://github.com/FreeFem/FreeFem-sources/tree/develop/examples/hpddm/MatLoad-PETSc.edp)                     | Loading a distributed matrix/vector with various formats                 |

### Reproducible science
| Article                                                                                                                                                                                             | Source code                                                   |
|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------|
| [Augmented Lagrangian preconditioner for large-scale hydrodynamic stability analysis](https://www.sciencedirect.com/science/article/pii/S0045782519301914)                                          | [GitHub repository](https://github.com/prj-/moulin2019al)     |
| [A multilevel Schwarz preconditioner based on a hierarchy of robust coarse spaces](https://hal.archives-ouvertes.fr/hal-02151184/document)                                                          | [GitHub repository](https://github.com/prj-/aldaas2019multi)  |
| [KSPHPDDM and PCHPDDM: extending PETSc with advanced Krylov methods and robust multilevel overlapping Schwarz preconditioners](https://www.sciencedirect.com/science/article/pii/S0898122121000055) | [GitHub repository](https://github.com/prj-/jolivet2020petsc) |

## PCHPDDM with subdomains shared by MPI processes

[schwarz-2d-PETSc.edp](schwarz-2d-PETSc.edp) exercises PCGASM with shared
subdomains and PCASM with several subdomains per process.
[schwarz-neumann-2d-PETSc.edp](schwarz-neumann-2d-PETSc.edp) builds distributed
Neumann auxiliary matrices, checks them against sequential references, solves
a manufactured reaction-diffusion problem with PCHPDDM, and plots the solution.

The PCHPDDM example requires these matching development branches:

| Project | Repository and branch | Tested commit |
| --- | --- | --- |
| FreeFEM | this repository, `jolivet/pchpddm-with-gasm` | this branch |
| PETSc | [indiano68/petsc](https://gitlab.com/indiano68/petsc), `ci-test` | `99ef6b11cb3cc2132929c8460223fd2384e5b7ba` |
| HPDDM | [prj-/hpddm-dev](https://github.com/prj-/hpddm-dev), `fabrizzi/multprocpersub` | `4f83feb526fb22484bcbb76cc30765a65e01ff01` |

Build PETSc with real, double-precision scalars, shared libraries, SLEPc, MUMPS,
and the HPDDM revision above. Use the same MPI installation for PETSc, FreeFEM,
and the launcher. For example, add these options to the PETSc configuration
appropriate for the target machine:

```sh
./configure --with-scalar-type=real --with-precision=double \
  --with-shared-libraries=1 --download-slepc --download-mumps \
  --download-scalapack --download-metis --download-parmetis \
  --download-hpddm=https://github.com/prj-/hpddm-dev/archive/4f83feb526fb22484bcbb76cc30765a65e01ff01.tar.gz
```

After building PETSc, configure and rebuild FreeFEM against that installation.
For an in-place PETSc build, export `PETSC_DIR` and the target's `PETSC_ARCH`,
then add this option to the usual FreeFEM configuration:

```sh
./configure --with-petsc="${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables"
```

For a prefix-installed PETSc, use its `lib/petsc/conf/petscvariables` instead.
FreeFEM takes the matching HPDDM include directory from this file; rebuild its
PETSc plugin after changing either dependency. The examples load `PETSc` and
`macro_ddm.idp` through FreeFEM's normal search paths, with no architecture name,
library suffix, or installation directory embedded in their source.

From this directory, run four subdomains on eight processes:

```sh
ff-mpirun -np 8 -nw schwarz-neumann-2d-PETSc.edp \
  -global 24 -subdomains 4 -hpddm_setup -ksp_view \
  -pc_hpddm_levels_1_eps_view \
  -pc_hpddm_levels_1_sub_mat_mumps_icntl_4 2 \
  -fglut solution.ffglut
```

The example enables inertia with an absolute eigenvalue threshold of `0.1`,
MUMPS Cholesky/LDLT, and `-pc_hpddm_levels_1_st_share_sub_ksp`. There is no fixed
eigenvector count. The views show the negative-pivot counts, the requested
numbers of eigenvectors, and whether the subdomain KSP is shared.

Change `-np` and `-subdomains` independently: tested pairs are `8/4`, `4/4`,
`4/2`, `4/3`, and `1/1`. PCHPDDM infers the grouping from the auxiliary matrices,
including unequal group sizes; do not also specify
`-pc_hpddm_total_subdomains`. The Neumann example requires one active subdomain
per rank, so the subdomain count must not exceed the process count. In serial,
PETSc uses a single level and does not need the eigensolver options.

`-fglut` enables plotting with `-nw` and saves a replayable graphics stream.
The final plot is also saved as `schwarz-neumann-solution.eps` in the current
directory; choose another relative filename with `-plot solution.eps`.
Without `-nw`, the usual FreeFEM graphical viewer displays the solution.
