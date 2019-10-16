# What is ffddm ?

In the acronym `ffddm`, `ff` stands for  FreeFEM and `ddm` for domain decomposition methods. The idea behind ffddm is to simplify the use of parallel solvers in FreeFEM: distributed direct methods and domain decomposition methods.  

Parallelism is an important issue because, since about 2004, the clock speed of cores stagnates at 2-3 GHz. The increase in performance is almost entirely due to the increase in the number of cores per processor. All major processor vendors are producing multicore chips and now every machine is a parallel machine. Waiting for the next generation machine does not guarantee anymore a better performance of a software. To keep doubling performance parallelism must double. It implies a huge effort in algorithmic development.

Thanks to `ffddm`, FreeFEM users have access to high-level functionalities for specifying and solving their finite element problems in parallel. The first task handled by `ffddm` is the data distribution among the processors. This is done via an overlapping domain decomposition and a related distributed linear algebra. Then, solving a linear system is possible either via an interface to the parallel [MUMPS](http://mumps.enseeiht.fr/) solver or by using domain decomposition methods as preconditioners to the GMRES Krylov method. The `ffddm` framework makes it easy to use scalable Schwarz methods enhanced by a coarse space correction built either from a coarse mesh or a [GenEO](https://link.springer.com/article/10.1007%2Fs00211-013-0576-y#page-1) (Generalized Eigenvalue in the Overlap) coarse space, see also the book [An Introduction to Domain Decomposition Methods: algorithms, theory, and parallel implementation](http://bookstore.siam.org/ot144/). State-of-the-art three level methods are also implemented in `ffddm`.  

The `ffddm` framework is entirely written in the FreeFEM language and the '.idp' scripts can be found [here](https://github.com/FreeFem/FreeFem-sources/blob/develop/idp) ('ffddm*.idp' files). It makes it also a very good tool for learning and prototyping domain decomposition methods without compromising efficiency.

 `ffddm` can also act as a wrapper for the [HPDDM](https://github.com/hpddm/hpddm) library. HPDDM is an efficient implementation of various domain decomposition methods and a variety of Krylov subspace algorithms, with advanced block and recycling methods for solving sequences of linear systems with multiple right-hand sides: GMRES and Block GMRES, CG, Block CG, and Breakdown-Free Block CG, GCRO-DR and Block GCRO-DR. For more details on how to use HPDDM within `ffddm`, see [the ffddm documentation](https://doc.freefem.org/documentation/ffddm/documentation.html#ffddmdocumentationhpddmffddm).

# Getting Started

## Minimal Example

```cpp
macro dimension 2// EOM            // 2D or 3D
include "ffddm.idp"
mesh Th = square(50,50);    // global mesh
// Step 1: Decompose the mesh
ffddmbuildDmesh( P , Th , mpiCommWorld )
// Step 2: Define your finite element
macro def(u)  u // EOM
macro init(u) u // EOM
ffddmbuildDfespace( P , P , real , def , init , P2 )
// Step 3: Define your problem
macro grad(u) [dx(u), dy(u)] // EOM
macro Varf(varfName, meshName, VhName)
    varf varfName(u,v) = int2d(meshName)(grad(u)'* grad(v)) + int2d(meshName)(1*v)
                       + on(1, u = 0);  // EOM
ffddmsetupOperator( P , P , Varf )
PVhi ui, bi;
ffddmbuildrhs( P , Varf , bi[] )
// Step 4: Define the one level DD preconditioner
ffddmsetupPrecond( P , Varf )
// Step 5: Define the two-level GenEO Coarse Space
ffddmgeneosetup( P , Varf )
// Step 6: Solve the linear system with GMRES
PVhi x0i = 0;
ui[] = PfGMRES(x0i[], bi[], 1.e-6, 200, "right");
ffddmplot(P, ui, "u")
Pwritesummary
```

This example solves a Laplace problem in 2D in parallel with a two-level GenEO domain decomposition method. To try this example, just copy and paste the script above in a file 'test.edp' and run it on 2 cores with

```console
$ ff-mpirun -np 2 test.edp -wg
```

## Documentation

You can find the `ffddm` documentation [here](https://doc.freefem.org/documentation/ffddm/index.html).

## Examples

| File name                                                                                                                                       | $M^{-1}_1$ | $M^{-1}_2$  | inexact CS  | comments                                           |
| ----------------------------------------------------------------------------------------------------------------------------------------------- | ---------- | ----------- | ----------- | -------------------------------------------------- |
| [diffusion-3d-minimal-direct.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/diffusion-3d-minimal-direct.edp)       | &nbsp;     | &nbsp;      | &nbsp;      | direct solver <br> MUMPS                           |
| [diffusion-3d-minimal-ddm.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/diffusion-3d-minimal-ddm.edp)             | RAS        | GenEO       | &nbsp;      |                                                    |
| [diffusion-3d-simple.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/diffusion-3d-simple.edp)                       | RAS        | GenEO       | &nbsp;      | comparison with <br> direct solver                 |
| [diffusion-2d-thirdlevelgeneo.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/diffusion-2d-thirdlevelgeneo.edp)     | RAS        | GenEO       | RAS + GenEO | &nbsp;                                             |
| [elasticity-3d-simple.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/elasticity-3d-simple.edp)                     | RAS        | GenEO       | &nbsp;      | &nbsp;                                             |
| [elasticity-3d-thirdlevelgeneo.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/elasticity-3d-thirdlevelgeneo.edp)   | RAS        | GenEO       | RAS + GenEO | &nbsp;                                             |
| [Helmholtz-2d-simple.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-2d-simple.edp)                       | ORAS       | Coarse Mesh | &nbsp;      | &nbsp;                                             |
| [Helmholtz-2d-marmousi.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-2d-marmousi.edp)                   | ORAS       | Coarse Mesh | &nbsp;      | &nbsp;                                             |
| [Helmholtz-3d-simple.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-3d-simple.edp)                       | ORAS       | Coarse Mesh | &nbsp;      |                                                    |
| [Helmholtz-3d-overthrust.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-3d-overthrust.edp)               | ORAS       | &nbsp;      | &nbsp;      | &nbsp;                                             |
| [Helmholtz-2d-HPDDM-BGMRES.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Helmholtz-2d-HPDDM-BGMRES.edp)           | ORAS       | &nbsp;      | &nbsp;      | multi-rhs <br> Block GMRES <br> with HPDDM         |
| [Navier-2d-marmousi2.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Navier-2d-marmousi2.edp)                       | ORAS       | Coarse Mesh | &nbsp;      | &nbsp;                                             |
| [Maxwell-3d-simple.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Maxwell-3d-simple.edp)                           | ORAS       | Coarse Mesh | &nbsp;      | &nbsp;                                             |
| [Maxwell_Cobracavity.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Maxwell_Cobracavity.edp)                       | ORAS       | Coarse Mesh | ORAS        | &nbsp;                                             |
| [natural_convection.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/natural_convection.edp)                         | ORAS       | Coarse Mesh |             | nonlinear                                          |
| [natural_convection_3D_obstacle.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/natural_convection_3D_obstacle.edp) | ORAS       | Coarse Mesh | &nbsp;      | nonlinear                                          |
| [Richards-2d.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/Richards-2d.edp)                                       | RAS        | &nbsp;      | &nbsp;      | nonlinear <br> time dependent <br> mesh adaptation |
| [heat-torus-3d-surf.edp](https://github.com/FreeFem/FreeFem-sources/blob/develop/examples/ffddm/heat-torus-3d-surf.edp)                         | RAS        | GenEO       | &nbsp;      | 3d surface <br> time dependent                     |
