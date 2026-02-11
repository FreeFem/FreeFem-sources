---
name: coupledphysics
category: Coupled Physics
layout: mpi
---

## Coupled-physics problem in 2D

In this example, we will solve a coupled-physics problem with two domains in 2D using domain decomposition, Krylov solvers and LU. The setup consists of two domains that share an interface. On each domain a reaction-diffusion equation is defined,

$$
-\nabla u_1 + \eta u_1 = f_1\text{ on }\Omega_1\\
-\nabla u_2 + \eta u_2 = f_2\text{ on }\Omega_2
$$

where $u_1$ and $u_2$ correspond to the local solution on the two domains, $\Omega_1$ and $\Omega_2$, and $f_1$ and $f_2$ are the two corresponding right hand sides. For this example, we assume that $f_1 = f_2 = f$. In addition to the local equations on the two domains, we also define some cross-interface conditions for how the two domains interact,

$$
\left(\frac{1}{\alpha}p_1 + u_1\right) + \left(\frac{1}{\alpha}p_2 - u_2\right) = 0 \text{ on } \Gamma_1\\
\left(\frac{1}{\alpha}p_12+ u_2\right) + \left(\frac{1}{\alpha}p_1 - u_1\right) = 0 \text{ on } \Gamma_2
$$

with $\alpha$ a fixed constant Robin parameter, and $p_1$ and $p_2$ are the solutions defined on the interface $\Gamma_1$ and $\Gamma_2$ respectively. The $\alpha$ value is calculated according to the formula,

$$
\alpha = \left[\left(\left(\frac{\pi}{L}\right)^2 + 1\right) \left(\left(\frac{\pi}{h}\right)^2 + 1\right)\right]^{\frac{1}{4}}
$$

### Weak form

We turn this into its weak form by using the usual test function approach, resulting in the definition of the local-to-domain problem,

$$
\mathscr{A}_i = \begin{bmatrix}L_{ii} & B_{ii}\\B_{ii}^T & C_{ii}\end{bmatrix}
$$

for $i=1,2$ as given by

$$
L_{ii} \leftarrow \int\int_{\Omega_i}(\nabla u_i \nabla v_i + \eta u_i v_i)\\
B_{ii} \leftarrow -\int_{\Gamma_i}p_i v_i \\
B_{ii}^T \leftarrow -\int_{\Gamma_i} u_i q_i \\
C_{ii} \leftarrow -\frac{1}{\alpha}\int_{\Gamma_i}u_i v_i.
$$

Similarly, the weak form for the cross-interface matrices between $\Omega_1$ and $\Omega_2$ are defined as

$$
B_{12} \leftarrow - \int_{\Gamma_1} p_2v_1\\
B_{21} \leftarrow - \int_{\Gamma_2} p_1v_2\\
C_{12} \leftarrow - \frac{1}{\alpha} \int_{\Gamma_1} p_2q_1\\
C_{21} \leftarrow - \frac{1}{\alpha} \int_{\Gamma_2} p_1q_2.
$$

### Preparing constants and building the mesh

First, we include the necessary headers and define some global constants and compute $\alpha$,

~~~freefem
include "getARGV.idp"
load "PETSc"
// Controlling the number of dof's
int nn = getARGV("-n", 200);
// domain 1
func Pku1 = P1; func Pkp1 = P1;
// domain 2
func Pku2 = P1; func Pkp2 = P1;
// alpha
real tmp1 = pow(pi/2,2) + 1;
real tmp2 = pow((nn-1)*pi, 2) + 1;
real tmp3 = tmp1*tmp2;
real alpha = pow(tmp3, 0.25);
~~~

Next, we build a mesh (a star in this example) and split it into two halves along the line $x=0$. We also extract the interface meshes along the split:

~~~freefem
border star(t=0, 2*pi){
    x = 1.0*(1 + 0.3*sin(5*t))*cos(t);
    y = 1.0*(1 + 0.3*sin(5*t))*sin(t);
}
mesh Th = buildmesh(star(nn));
mesh Th1 = trunc(Th, x<=0, split=1, label=123);
mesh Th2 = trunc(Th, x>=0, split=1, label=123);
int[int] labs = [123];
meshL ThI1 = extract(Th1, label=labs);
meshL ThI2 = extract(Th2, label=labs);
~~~

### Domain decomposition

In order to use domain decomposition, we will take advantage of PETSc's capabilities through the `macro_ddm` script in FreeFEM. However, the datatype `mesh` requires `macro_ddm` to be set up for a 2D volume mesh, whereas the `meshL` datatype requires a 3D line mesh setup. To accommodate this, we will scope the relevant parts of the code doing the decomposition. To prepare for this, we need to prepare some variables that will get defined inside of those scoped regions but are also needed outside of them:

~~~freefem
Mat L11, L22;
Mat C11, C22;
Mat interpI12, interpI21;

fespace Ph1(Th1, P0);
fespace Ph2(Th2, P0);
Ph1 part1;
Ph2 part2;
~~~

The `fespace`s defined will hold the decomposition which will allow us to reuse the decomposition of the surface mesh in the line mesh along the interface. Next, we decompose our surface meshes using a custom partitioning performed with `metis`, storing the partitioning in `part1` and `part2`:

~~~freefem
{
    macro dimension()2// EOM
    macro meshtype()V// EOM
    macro partitioner()metis//
    include "macro_ddm.idp"
    partitionerSeq(part1[], Th1, mpisize);
    partitionerSeq(part2[], Th2, mpisize);
    partitionerPar(part1[], Th1, mpiCommWorld, mpisize);
    partitionerPar(part2[], Th2, mpiCommWorld, mpisize);
    macro Th1UserPartitioning()part1[]// EOM
    macro Th2UserPartitioning()part2[]// EOM
    MatCreate(Th1, L11, Pku1);
    MatCreate(Th2, L22, Pku2);
}
~~~

Following this, `L11` and `L22` are distributed PETSc matrices, according to the decomposition, that we will later fill with the relevant data. We now will interpolate the decomposition along the interface in order to have matching decompositions along the interface meshes:

~~~freefem
fespace PhI1(ThI1, P0);
fespace PhI2(ThI2, P0);
PhI1 partI1 = part1;
PhI2 partI2 = part2;
~~~

There are also certain operations that are not defined if a local subdomain does not include any part of the interface,

~~~freefem
int[int] labs1 = labels(Th1);
int[int] labs2 = labels(Th2);
int hasI1 = 0, hasI2 = 0;
for(int i = 0; i < labs1.n; i++) if(labs1[i] == 123) hasI1 = 1;
for(int i = 0; i < labs2.n; i++) if(labs2[i] == 123) hasI2 = 1;
~~~

We can now decompose the interface meshes in a similar way as we did with the surface meshes. Since we already have a partitioning, however, we don't need to create a new one,

~~~freefem
{
    macro dimension()3// EOM
    macro meshtype()L// EOM
    include "macro_ddm.idp"
    macro ThI1UserPartitioning()partI1[]// EOM
    macro ThI2UserPartitioning()partI2[]// EOM
    MatCreate(ThI1, C11, Pkp1);
    MatCreate(ThI2, C22, Pkp2);
    MatInterpolate(ThI1, Pkp1, C11, ThI2, Pkp2, C22, interpI21);
    MatInterpolate(ThI2, Pkp2, C22, ThI1, Pkp1, C11, interpI12);
}
~~~

In this scoped section, we also set up two interpolation matrices in order to move from the interface mesh in one of $\Omega_1$ or $\Omega_2$ to the other. Since the decompositions of the two domains is not assumed to match along the interface, we will need to set up our cross-interface matrices in one of the domains and project it across the interface as we will see later.

With all the meshes decomposed, we can now set up our distributed finite element spaces,

~~~freefem
fespace Vh1(Th1, Pku1);
fespace Vh2(Th2, Pku2);
fespace VhI1(ThI1, Pkp1);
fespace VhI2(ThI2, Pkp2);
~~~

### Problem specification

We define our problem, including the corresponding right-hand side, in the usual FreeFEM manner,

~~~freefem
func f = x^2*y*100;
macro grad(u) [dx(u),dy(u)] // EOM
varf ProbL1(u, v) =  int2d(Th1)(grad(u)'*grad(v) + 1e-1*u*v) - int2d(Th1)(f * v);
varf ProbB1(u, v) =  int1d(ThI1)(u*v);
varf ProbC1(u, v) = -int1d(ThI1)(1./alpha * u*v);
varf ProbL2(u, v) =  int2d(Th2)(grad(u)'*grad(v) + 1e-1*u*v) - int2d(Th2)(f * v);
varf ProbB2(u, v) =  int1d(ThI2)(u*v);
varf ProbC2(u, v) = -int1d(ThI2)(1./alpha * u*v);
real[int] rhs(Vh1.ndof + Vh2.ndof + VhI1.ndof + VhI2.ndof);
rhs = 0;
rhs(0:Vh1.ndof-1) = ProbL1(0,Vh1);
rhs(Vh1.ndof+VhI1.ndof:Vh1.ndof+VhI1.ndof+Vh2.ndof-1) = ProbL2(0,Vh2);
~~~

For simplicity, we define the same problem on both domains and use the same $f$ for both, however, this can be anything here.

### Constructing distributed matrices

We are now ready to set up and compose the various distributed block matrices that will make up the full system definition. We start with the matrices defining the problem local to $\Omega_1$ and $\Omega_2$,

~~~freefem
// L11 and L22
matrix mL11 = ProbL1(Vh1, Vh1, solver=GMRES);
matrix mL22 = ProbL2(Vh2, Vh2, solver=GMRES);
L11 = mL11;
L22 = mL22;
// B11 and B22
matrix mB11(Vh1.ndof, VhI1.ndof);
matrix mB22(Vh2.ndof, VhI2.ndof);
if(hasI1) mB11 = ProbB1(VhI1, Vh1, solver=GMRES);
if(hasI2) mB22 = ProbB2(VhI2, Vh2, solver=GMRES);
Mat B11(L11, C11);
Mat B22(L22, C22);
B11 = mB11;
B22 = mB22;
// C11 and C22
matrix mC11(VhI1.ndof,VhI1.ndof);
matrix mC22(VhI2.ndof,VhI2.ndof);
if(hasI1) mC11 = ProbC1(VhI1, VhI1, solver=GMRES);
if(hasI2) mC22 = ProbC2(VhI2, VhI2, solver=GMRES);
C11 = mC11;
C22 = mC22;
~~~

Note that the value of the `solver` parameter is completely irrelevant as we will solve the system through PETSc later-on. However, the default solver in distributed FreeFEM includes collective MPI operations, which does not allow us to only conditionally fill `B11`, `B22`, `C11`, and `C22`, as otherwise the execution would be stalled as soon as one of the subdomains does not include any part of the boundary.

The cross-interface matrices will be set up on the source domain only, as we cannot assume that we have matching decompositions on the other side of the interface:

~~~freefem
// BT12 and BT21
matrix mB12(Vh1.ndof, VhI1.ndof), mB21(Vh2.ndof, VhI2.ndof);
if(hasI1) mB12 = ProbB1(VhI1, Vh1, solver=GMRES);
if(hasI2) mB21 = ProbB2(VhI2, Vh2, solver=GMRES);
matrix mBT21 = -mB21';
matrix mBT12 = -mB12';
Mat BT12(C11, L11);
Mat BT21(C22, L22);
BT12 = mBT12;
BT21 = mBT21;
// C12 and C21
matrix mC12(VhI1.ndof,VhI1.ndof), mC21(VhI2.ndof, VhI2.ndof);
if(hasI1) mC12 = ProbC1(VhI1, VhI1, solver=GMRES);
if(hasI2) mC21 = ProbC2(VhI2, VhI2, solver=GMRES);
Mat C12(C11, C11);
Mat C21(C22, C22);
C12 = mC12;
C21 = mC21;
~~~

In order to project the target domain across the interface onto the right subdomains of that domain, we make use of the interpolation matrices set up initially. PETSc allows us to define a block matrix constructed out of block matrices and block matrix products, which allows us to “fix” the target domains of the cross-interface block matrices while defining the full system block matrix:

~~~freefem
Mat A = [[     L11      ,      B11     ,       0       ,       0      ],
         [     B11'     ,      C11     , interpI12*BT21, interpI12*C21],
         [      0       ,       0      ,      L22      ,      B22     ],
         [interpI21*BT12, interpI21*C12,      B22'     ,      C22     ]];
~~~

### Solving the system

For solving the problem we use PETSc's `fgmres` solver on the global system preconditioned by a solve of the two 2x2 diagonal blocks, which define each a saddle point problem local to domain $\Omega_1$ and $\Omega_2$ respectively. These two saddle point systems are solved with PETSc's `gmres` solver preconditioned in turn by `lu`:

~~~freefem
set(A, sparams=" -ksp_type fgmres" +
               " -ksp_pc_side right" +
               " -pc_type fieldsplit" +
               " -pc_fieldsplit_type additive" +
               " -pc_fieldsplit_0_fields 0,1" +
               " -pc_fieldsplit_1_fields 2,3" +
               " -fieldsplit_0_ksp_type preonly" +
               " -fieldsplit_0_pc_type lu" +
               " -fieldsplit_1_ksp_type preonly" +
               " -fieldsplit_1_pc_type lu" +
               " -ksp_max_it 100" +
               " -ksp_gmres_restart 100" +
               " -ksp_monitor_true_residual" +
               " -ksp_converged_reason" +
               " -ksp_gmres_modifiedgramschmidt");
~~~

The `fieldsplit` preconditioner in PETSc allows us to define a block-Jacobi type splitting, with each split containing the fields corresponding to the respective 2x2 block matrix. We are now ready to solve the problem and plot the solution:

~~~freefem
real[int] sol = A^-1 * rhs;
macro dimension()2// EOM
macro meshtype()V// EOM
include "macro_ddm.idp"
macro def(u)u//
Vh1 u1;
Vh2 u2;
u1[] = sol(0:Vh1.ndof-1);
u2[] = sol(Vh1.ndof+VhI1.ndof:Vh1.ndof+VhI1.ndof+Vh2.ndof-1);
plotMPI(Th1, u1, Pku1, def, real, cmm="Global velocity (u1)");
plotMPI(Th2, u2, Pku2, def, real, cmm="Global velocity (u2)");
~~~

For 16 MPI processes, a possible combined plot is shown below:

![][_solution]

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/coupledphysics/solution.png
