# What is ffddm ?

In the acronym `ffddm`, `ff` stands for  FreeFem++ and `ddm` for domain decomposition methods. The idea behind ffddm is to simplify the use of parallel solvers in FreeFem++: distributed direct methods and domain decomposition methods.  

Parallelism is an important issue because, since about 2004, the clock speed of cores stagnates at 2-3 GHz. The increase in performance is almost entirely due to the increase in the number of cores per processor. All major processor vendors are producing multicore chips and now every machine is a parallel machine. Waiting for the next generation machine does not guarantee anymore a better performance of a software. To keep doubling performance parallelism must double. It implies a huge effort in algorithmic development. 

Thanks to `ffddm`, FreeFem++ users have access to high-level functionalities for specifying and solving their finite element problems in parallel. The first task handled by `ffddm` is the data distribution among the processors. This is done via an overlapping domain decomposition and a related distributed linear algebra. Then, solving a linear system is possible either via an interface to the parallel [MUMPS](http://mumps.enseeiht.fr/) solver or by using domain decomposition methods as preconditioners to the GMRES Krylov method. The `ffddm` framework makes it easy to use scalable Schwarz methods enhanced by a coarse space correction built either from a coarse mesh or a [GenEO](https://link.springer.com/article/10.1007%2Fs00211-013-0576-y#page-1) (Generalized Eigenvalue in the Overlap) coarse space, see also the book [An Introduction to Domain Decomposition Methods: algorithms, theory, and parallel implementation](http://bookstore.siam.org/ot144/). State-of-the-art three level methods are also implemented in `ffddm`.  

The `ffddm` framework is entirely written in the FreeFem++ language and the '.idp' scripts can be found [here](https://github.com/FreeFem/FreeFem-sources/tree/master/examples%2B%2B-ffddm). It makes it also a very good tool for learning and prototyping domain decomposition methods without compromising efficiency.

 `ffddm` can also act as a wrapper for the [HPDDM](https://github.com/hpddm/hpddm) library. HPDDM is an efficient implementation of various domain decomposition methods and a variety of Krylov subspace algorithms, with advanced block and recycling methods for solving sequences of linear systems with multiple right-hand sides: GMRES and Block GMRES, CG, Block CG, and Breakdown-Free Block CG, GCRO-DR and Block GCRO-DR. For more details on how to use HPDDM within `ffddm`, see [here METTRE LE LIEN](DOC)
 

# Getting Started

[Tutorial slides](doc/tutorial-slides.html)

[Introduction to DD](doc/introddm.md)

[Complete Documentation](doc/documentation.md)

### Minimal Example

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
// Step 5: Solve the linear system with GMRES
PVhi x0i = 0;
ui[] = PfGMRES(x0i[], bi[], 1.e-6, 200, "right");
ffddmplot(P, ui, "u")
Pwritesummary
```

### Examples

nonlinear or time dependent ?

| File name                          | d   | $M^{-1}_1$ | $M^{-1}_2$  | inexact CS |
| ---------------------------------- | --- | ---------- | ----------- | ---------- |
| Helmholtz-2d-marmousi.edp          | 2D  | ORAS       | Coarse Mesh | X          |
| Helmholtz-2d-simple.edp            | 2D  | ORAS       | Coarse Mesh | X          |
| Helmholtz-3d-overthrust.edp        | 3D  | ORAS       | X           | X          |
| Helmholtz-3d-simple.edp            | 3D  | ORAS       | Coarse Mesh | X          |
| Maxwell_Cobracavity.edp            | 2D  | ORAS       | Coarse Mesh | ORAS       |
| Navier-2d-marmousi2.edp            | 2D  | ORAS       | Coarse Mesh | X          |
| Richards-2d.edp                    | 2D  |            |             |            |
| diffusion-2d-thirdlevelgeneo.edp   | 2D  |            |             |            |
| diffusion-3d-intermediate.edp      | 3D  |            |             |            |
| diffusion-3d-minimal-ddm.edp       | 3D  |            |             |            |
| diffusion-3d-minimal-direct.edp    | 3D  |            |             |            |
| diffusion-3d-simple.edp            | 3D  |            |             |            |
| elasticity-3d-simple.edp           | 3D  |            |             |            |
| elasticity-3d-thirdlevelgeneo.edp  | 3D  |            |             |            |
| natural_convection.edp             | 2D  |            |             |            |
| natural_convection_3D_obstacle.edp | 3D  |            |             |            |

