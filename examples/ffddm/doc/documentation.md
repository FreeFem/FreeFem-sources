<a name="test"></a>

## Minimal example

```cpp
macro dimension 3// EOM            // 2D or 3D

include "ffddm.idp"

int[int] LL = [2,2, 1,2, 2,2];
mesh3 ThGlobal = cube(10, 10, 10, [x, y, z], label = LL);      // global mesh

macro grad(u) [dx(u), dy(u), dz(u)]// EOM    // three-dimensional gradient

macro Varf(varfName, meshName, VhName)
    varf varfName(u,v) = int3d(meshName)(grad(u)''* grad(v)) + int3d(meshName)(v) + on(1, u = 1.0);
// EOM
       
// Domain decomposition
ffddmbuildDmesh( Lap , ThGlobal , mpiCommWorld )

macro def(i)i// EOM                         // scalar field definition
macro init(i)i// EOM                        // scalar field initialization
ffddmbuildDfespace( Lap , Lap , real , def , init , P1 )

ffddmsetupOperator( Lap ,Lap , Varf )

real[int] rhsi(0);
ffddmbuildrhs( Lap , Varf , rhsi )

LapVhi def(ui);

//Direct solve
ui[] = Lapdirectsolve(rhsi);

Lapwritesummary

ffddmplot(Lap,ui,"u");
```

## Overlapping mesh decomposition

```cpp
ffddmbuildDmesh(pr,Th,comm)
```

decomposes the mesh **Th** into overlapping submeshes. The mesh will be distributed over the mpi ranks of communicator **comm**. This will create and expose variables whose names will be prefixed by **pr**, see below (# is the concatenation operator).
The way the initial mesh **Th** is partitioned depends on the value of [ffddmpartitioner](parameters.md#global-parameters).  
The size of the overlap between subdomains (its width in terms of number of mesh elements) is given by [ffddmoverlap](parameters.md#global-parameters).  
The level of refinement of the resulting submeshes with respect to the input mesh **Th** is given by [ffddmsplit](parameters.md#global-parameters).  
If [ffddmexclude](parameters.md#global-parameters) $\neq 0$, the first [ffddmpCS](parameters.md#global-parameters) mpi ranks of **comm** will be excluded from the spatial domain decomposition, in order to dedicate them later to the coarse problem (for two-level preconditioners).  
The label of the new border of the submeshes (the interface between the subdomains) is given by [ffddminterfacelabel](parameters.md#global-parameters).

***defines***:

- `int pr#npart` number of subdomains for this decomposition; should be equal to mpiSize(**comm**) - [ffddmexclude](parameters.md#global-parameters) * [ffddmpCS](parameters.md#global-parameters)
- `meshN[int] pr#aTh` array (size `pr#npart`) of local meshes of the subdomains. In the standard parallel case, only the local mesh for this mpi rank `pr#aTh[mpiRank(pr#commddm)]` is defined (unless this mpi rank is excluded from the spatial domain decomposition, i. e. `prmesh#excluded` = 1, see below). In the sequential case, all local meshes are defined.
- `meshN pr#Thi` the local mesh of the subdomain for this mpi rank, i. e. `pr#aTh[mpiRank(pr#commddm)]` in the parallel case
- `int pr#numberIntersection` the number of neighbors for this mpi rank
- `int[int] pr#arrayIntersection` the list of neighbor ranks in `pr#commddm` for this mpi rank
- `int pr#pCS` equal to [ffddmpCS](parameters.md#global-parameters)
- `int pr#exclude` equal to [ffddmexclude](parameters.md#global-parameters)
- `int pr#excluded` *true* if [ffddmexclude](parameters.md#global-parameters) is *true* ($\neq 0$) and mpiRank(**comm**) < `pr#pCS`. In this case, this mpi rank will be excluded from the spatial domain decomposition and will only work on the coarse problem.
- `mpiComm pr#commddm` mpi communicator for ranks participating in the spatial domain decomposition (ranks 0 to `pr#npart`-1 in **comm** if `pr#exclude` is *false*, ranks `pr#pCS` to `pr#pCS`+`pr#npart`-1 otherwise)
- `mpiComm pr#commCS` mpi communicator for ranks participating in the assembly and resolution of the coarse problem for two-level preconditioners (ranks 0 to `pr#pCS` - 1 in **comm**)
- `mpiComm pr#commself` self mpi communicator (this mpi rank only), used for factorizing local matrices

***For advanced users***:

- `int pr#binexactCS` 
- `int pr#inexactCSsplit`
- `int pr#isincomm`
- `meshN[int] pr#aThborder`

***Remark for sequential use*** (see [``-seqddm``](parameters.md#command-line-arguments)):
- `meshN[int] pr#aTh` array (size `pr#npart`) of local meshes of the subdomains

#int pr#binexactgeneoCS

#fespace pr#VhiP1(pr#Thi,P1);

#pr#VhiP1[int] pr#partitionIntersectionbasei(0);

#meshN pr#Thglob = minimalMesh;

#matrix[int] pr#RihP1(pr#npart);
#pr#VhiP1[int] pr#DP1(pr#npart);

#NewMacro pr#mpicomm()comm EndMacro

***depends on***:
- [ffddmpartitioner](parameters.md#global-parameters)
- [ffddmpCS](parameters.md#global-parameters)
- [ffddmexclude](parameters.md#global-parameters)
- [ffddmoverlap](parameters.md#global-parameters)
- [ffddmsplit](parameters.md#global-parameters)
- [ffddminterfacelabel](parameters.md#global-parameters)

***see also***:

## Local finite element spaces

```cpp
ffddmbuildDfespace(pr,prmesh,scalar,def,init,Pk)
```

builds the local finite element spaces and associated distributed operators on top of the mesh decomposition **prmesh**. This will create and expose variables whose names will be prefixed by **pr**, see below. It is assumed that **[`ffddmbuildDmesh`](#overlapping-mesh-decomposition)** has already been called with prefix **prmesh** in order to build the mesh decomposition.  
The local finite element spaces of type **Pk** (where **Pk** is the type of finite element: P1, [P2,P2,P1], ...) are defined on the local meshes of the subdomains based on the mesh decomposition previously created with prefix **prmesh**.  
**scalar** determines the type of data for this finite element: *real* or *complex*.  
Two macros, **def** and **init**, are needed: **def** specifies how to define a finite element function in the finite element space **Pk**, and **init** specifies how to interpolate a scalar function onto the (possibly multiple) components of **Pk**. Two examples are given below:  

For scalar P2 finite elements and complex-valued problems:
```cpp
macro def(u) u// EOM
macro init(u) u// EOM
ffddmbuildDfespace(myFEprefix,mymeshprefix,complex,def,init,P2)
```

For vectorial [P2,P2,P1] finite elements and real-valued problems:
```cpp
macro def(u) [u, u#B, u#C]// EOM
macro init(u) [u, u, u]// EOM
ffddmbuildDfespace(myFEprefix,mymeshprefix,real,def,init,[P2,P2,P1])
```

In practice, this builds the necessary distributed operators associated to the finite element space: the local partition of unity functions $(D_i)_{i=1,...,N}$ (see `pr#Dk` and `pr#Dih` below) as well as the function `pr#update` (see below) which synchronizes local vectors $(u_i)_{i=1,...,N}$ between neighboring subdomains, performing the equivalent of $u_i = R_i (\sum_{j=1}^N R_j^T u_j)$ or $u_i = R_i (\sum_{j=1}^N R_j^T D_j u_j)$ in a distributed parallel environment.  
`pr#scalprod` (see below) performs the parallel scalar product for vectors defined on this finite element.

***defines***:

- `pr#prmesh` macro, saves the parent prefix **prmesh** of the mesh decomposition
- `pr#K` macro, saves the type of data **scalar** for this finite element space (*real* or *complex*)
- `func pr#fPk` saves the type of finite element **Pk**, e.g. *P1, [P2,P2,P1], ...*
- `fespace pr#Vhi` the local finite element space for this mpi rank, defined on the local mesh `prmesh#Thi`
- `int pr#Ndofglobal` the total number of degrees of freedom $n$ for this finite element discretization
- `pr#mdef` macro, saves the macro **def** giving the definition of a finite element function in the finite element space **Pk**
- `pr#minit` macro, saves the macro **init** specifying how to interpolate a scalar function onto the (possibly multiple) components of a finite element function of **Pk**. This is used to create the local partition of unity function in `pr#Vhi`, by interpolating the local P1 partition of unity function onto the components of `pr#Vhi`. For non Lagrange finite element spaces (e.g. *RT0*, *Edge03d*, ...), see **[`ffddmbuildDfespaceEdge`](#local-finite-element-spaces)**. 
- `pr#K[int][int] pr#Dk` array (size `prmesh#npart`) of local partition of unity vectors in the subdomains, equivalent to $(D_i)_{i=1,...,N}$. In the standard parallel case, only the local partition of unity vector for this mpi rank `pr#Dk[mpiRank(prmesh#commddm)]` is defined (unless this mpi rank is excluded from the spatial domain decomposition, i. e. `prmesh#excluded` = 1). In the sequential case, all local partition of unity vectors are defined.
- `matrix<pr#K>[int] pr#Dih` array (size `prmesh#npart`) similar to `pr#Dk` but in *matrix* form, allowing for easier *matrix*-*matrix* multiplications. `pr#Dih[i]` is a diagonal matrix, with the diagonal equal to `pr#Dk[i]`.
- `fespace pr#Vhglob` the global finite element space defined on the global mesh `prmesh#Thglob`. Defined only if [``-noGlob``](parameters.md#command-line-arguments) is not used.
- `matrix<pr#K>[int] pr#Rih` array (size `prmesh#npart`) of restriction matrices from the global finite element space to the local finite element spaces on the local submeshes of the subdomains. In the standard parallel case, only the restriction matrix for this mpi rank `pr#Rih[mpiRank(prmesh#commddm)]` is defined (unless this mpi rank is excluded from the spatial domain decomposition, i. e. `prmesh#excluded` = 1). In the sequential case, all restriction matrices are defined. The restriction matrices `pr#Rih` are defined only if [``-noGlob``](parameters.md#command-line-arguments) is not used.
- `func int pr#update(scalar[int] ui, bool scale)` The function `pr#update` synchronizes the local vector *ui* between subdomains by exchanging the values of *ui* shared with neighboring subdomains (in the overlap region) using point-to-point MPI communications. If *scale* is *true*, *ui* is multiplied by the local partition of unity beforehand. This is equivalent to
$u_i = R_i (\sum_{j=1}^N R_j^T u_j)$ when *scale* is *false* and $u_i = R_i (\sum_{j=1}^N R_j^T D_j u_j)$ when *scale* is *true*.
- `func scalar pr#scalprod(scalar[int] ai, scalar[int] bi)` The function `pr#scalprod` computes the global scalar product of two vectors whose local restriction to the subdomain of this mpi rank are *ai* and *bi*. The result is computed as $\sum_{j=1}^N (D_j a_j, b_j)$.

***Remark:***


***For advanced users***:

matrix<pr#K>[int] pr#restrictionIntersection(0);

NewMacro pr#mdefpart udefpart EndMacro

NewMacro pr#minitpart uinitpart EndMacro

func pr#fPkP0 = mPkP0;

pr#K[int][int] pr#rcv(0);
pr#K[int][int] pr#snd(0);

***depends on***:

***see also***:

- **[`ffddmbuildDfespaceEdge`](#local-finite-element-spaces)**

## Define the problem to solve

```cpp
ffddmsetupOperator(pr,prfe,Varf)
```

builds the distributed operator associated to the variational problem given by **Varf**, on top of the distributed finite element space **prfe**. This will create and expose variables whose names will be prefixed by **pr**, see below. It is assumed that **[`ffddmbuildDfespace`](#local-finite-element-spaces)** has already been called with prefix **prfe** in order to define the distributed finite element space.  
In practice, this builds the so-called local 'Dirichlet' matrices $A_i = R_i A R_i^T$, the restrictions of the global operator $A$ to the subdomains (see `pr#aRd`below). The matrices correspond to the discretization of the bilinear form given by the macro **Varf**, which represents the abstract variational form of the problem. These matrices are then used to implement the action of the global operator $A$ on a local vector (the parallel matrix-vector product with $A$), see `pr#A` below.  
At this point, we already have the necessary data to be able to solve the problem with a parallel direct solver (*MUMPS*), which is the purpose of the function `pr#directsolve` (see below). See **[`ffddmbuildrhs`](#define-the-problem-to-solve)** for building the right-hand side.  

The macro **Varf** is required to have three parameters: the name of the variational form, the mesh, and the finite element space. The variational form given in this 'abstract' format will then be used by *ffddm* to assemble the discrete operators by setting the appropriate mesh and finite element space as parameters. An example is given below:

```cpp
macro myVarf(varfName, meshName, VhName)
    varf varfName(u,v) = int3d(meshName)(grad(u)''* grad(v)) + on(1, u = 1.0);
// EOM

ffddmsetupOperator(myprefix,myFEprefix,myVarf)

```
***Remark*** In this simple example, the third parameter *VhName* is not used. However, for more complex cases such as non-linear or time dependent problems where the problem depends on a solution computed at a previous step, it is useful to know for which discrete finite element space the variational form is being used. See for example TODO  

***defines***:

- `pr#prfe` macro, saves the parent prefix **prfe** of the finite element space
- `int pr#verbosity` the level of verbosity for this problem, initialized with the value of [ffddmverbosity](parameters.md#global-parameters)
- `pr#writesummary` macro, prints a summary of timings for this problem, such as the time spent to assemble local matrices or solve the linear system.
- `matrix<prfe#K> pr#Aglobal` the global matrix $A$ corresponding to the discretization of the variational form given by the macro **Varf** on the global finite element space `prfe#Vhglob`. Defined only in the sequential case.
- `matrix<prfe#K>[int] pr#aRd` array (size `prfe#prmesh#npart`) of so-called local 'Dirichlet' matrices in the subdomains; these are the restrictions of the global operator to the subdomains, equivalent to $A_i = R_i A R_i^T$ with $A$ the global matrix corresponding to the discretization of the variational form given by the macro **Varf** on the global finite element space. In the standard parallel case, only the local matrix for this mpi rank `pr#aRd[mpiRank(prmesh#commddm)]` is defined (unless this mpi rank is excluded from the spatial domain decomposition, i. e. `prmesh#excluded` = 1). In the sequential case, all local matrices are defined.
- `func prfe#K[int] pr#A(prfe#K[int] &ui)` The function `pr#A` computes the parallel matrix-vector product, i.e. the action of the global operator $A$ on the local vector $u_i$. The computation is equivalent to $R_i (\sum_{j=1}^N R_j^T D_j A_j u_j)$ and is performed in parallel using local matrices `pr#aRd` and the function `prfe#update`. In the sequential case, the global matrix **pr#Aglobal** is used instead.
- `func prfe#K[int] pr#directsolve(prfe#K[int]& rhsi)` The function `pr#directsolve` allows to solve the linear system $A x = b$ in parallel using the parallel direct solver *MUMPS*. The matrix is given to *MUMPS* in distributed form through the local matrices `pr#aRd`. The input *rhsi* is given as a distributed vector (*rhsi* is the restriction of the global right-hand side $b$ to the subdomain of this mpi rank, see **[`ffddmbuildrhs`](#define-the-problem-to-solve))** and the returned vector is local as well.

NewMacro pr#plot(u,s)

***For advanced users***:

NewMacro pr#fromVhi(ui,VhName,res)

***depends on***:

- [ffddmverbosity](parameters.md#global-parameters)

```cpp
ffddmbuildrhs(pr,Varfrhs,rhs)
```

builds the right-hand side associated to the variational form given by **Varfrhs** for the problem corresponding to prefix **pr**.  The resulting right-hand side vector **rhs** corresponds to the discretization of the abstract linear form given by the macro **Varfrhs** (see **[`ffddmsetupOperator`](#define-the-problem-to-solve)** for more details on how to define the abstract variational form as a macro).  
The input vector **rhs** is resized and contains the resulting local right-hand side $R_i b$, the restriction of the global right-hand side $b$ to the subdomain of this mpi rank. In the sequential case, the global right-hand side vector $b$ is assembled instead.  
An example is given below:

```cpp
macro myVarfrhs(varfName, meshName, VhName)
    varf varfName(u,v) = intN(meshName)(v) + on(1, u = 1.0);
// EOM

real[int] rhsi(0);
ffddmbuildrhs(myprefix,myVarfrhs,rhsi)
```

## One level preconditioners

```cpp
ffddmsetupPrecond(pr,VarfPrec)
```

builds the one level preconditioner for problem **pr**. This will create and expose variables whose names will be prefixed by **pr**, see below. It is assumed that **[`ffddmsetupOperator`](#define-the-problem-to-solve)** has already been called with prefix **pr** in order to define the problem to solve.  
In practice, this builds and performs the factorization of the local matrices used in the one level preconditioner. The local matrices depend on the choice of [ffddmprecond](parameters.md#global-parameters) and **VarfPrec**, see `pr#aR`below.

***defines***:

- `string pr#prec` equal to [ffddmprecond](parameters.md#global-parameters). Sets the type of one level preconditioner $M^{-1}_1$ to be used: "asm" (*Additive Schwarz*), "ras" (*Restricted Additive Schwarz*), "oras" (*Optimized Restricted Additive Schwarz*), "soras" (*Symmetric Optimized Restricted Additive Schwarz*) or "none" (no preconditioner).

- `matrix<pr#prfe#K>[int] pr#aR` array (size `prfe#prmesh#npart`) of local matrices used for the one level preconditioner. Each mpi rank of the spatial domain decomposition performs the $LU$ (or $LDL^T$) factorization of the local matrix corresponding to its subdomain using the direct solver *MUMPS*.

   - If **VarfPrec** is not a previously defined macro (just put *null* for example), the matrices `pr#aR` are set to be equal to the so-called local 'Dirichlet' matrices `pr#aRd` (see **[`ffddmsetupOperator`](#define-the-problem-to-solve)**). This is for the classical ASM preconditioner $M^{-1}_1 = M^{-1}_{\text{ASM}} = \sum_{i=1}^N R_i^T A_i^{-1} R_i$ or classical RAS preconditioner $M^{-1}_1 = M^{-1}_{\text{RAS}} = \sum_{i=1}^N R_i^T D_i A_i^{-1} R_i$ (it is assumed that [ffddmprecond](parameters.md#global-parameters) is equal to "asm" or "ras").
   - If **VarfPrec** is a macro, it is assumed that **VarfPrec** defines an abstract bilinear form (see **[`ffddmsetupOperator`](#define-the-problem-to-solve)** for more details on how to define the abstract variational form as a macro).
    - If [ffddmprecond](parameters.md#global-parameters) is equal to "asm" or "ras", the matrices `pr#aR` will be assembled as local 'Dirichlet' matrices in the same manner as `pr#aRd`, but using the bilinear form defined by **VarfPrec** instead. This defines the ASM preconditioner as $M^{-1}_1 = M^{-1}_{\text{ASM}} = \sum_{i=1}^N R_i^T {(A_i^{\text{Prec}})}^{-1} R_i$ and the RAS preconditioner as $M^{-1}_1 = M^{-1}_{\text{RAS}} = \sum_{i=1}^N R_i^T D_i {(A_i^{\text{Prec}})}^{-1} R_i$, where $A_i^{\text{Prec}} = R_i A^{\text{Prec}} R_i^T$.
    - If [ffddmprecond](parameters.md#global-parameters) is equal to "oras" or "soras", the matrices `pr#aR` will correspond to the discretization of the variational form **VarfPrec** in the subdomains $\Omega_i$. In particular, various boundary conditions can be imposed at the interface between subdomains (corresponding to mesh boundary of label [ffddminterfacelabel](parameters.md#global-parameters) set by the parent call to **[`ffddmbuildDmesh`](#overlapping-mesh-decomposition)**), such as Optimized Robin boundary conditions. We note the ORAS preconditioner as $M^{-1}_1 = M^{-1}_{\text{ORAS}} = \sum_{i=1}^N R_i^T D_i {(B_i^{\text{Prec}})}^{-1} R_i$ and the SORAS preconditioner as $M^{-1}_1 = M^{-1}_{\text{SORAS}} = \sum_{i=1}^N R_i^T D_i {(B_i^{\text{Prec}})}^{-1} D_i R_i$.
- `func pr#prfe#K[int] pr#PREC1(pr#prfe#K[int] &ui)` The function `pr#PREC1` computes the parallel application of the one level preconditioner $M^{-1}_1$, i.e. the action of $M^{-1}_1$ on the local vector $u_i$. In the sequential case, it computes the action of $M^{-1}_1$ on a global vector. The action of the inverse of local matrices `pr#aRd` is computed by forward-backward substitution using their $LU$ (or $LDL^T$) decomposition.
- `func pr#prfe#K[int] pr#PREC(pr#prfe#K[int] &ui)` The function `pr#PREC` corresponds to the action of the preconditioner $M^{-1}$ for problem **pr**. It coincides with the one level preconditioner `pr#PREC1` after the call to **[`ffddmsetupPrecond`](#one-level-preconditioners)**. If a second level is subsequently added (see the next section about **[Two level preconditioners](#two-level-preconditioners)**), it will then coincide with the two level preconditioner $M^{-1}_2$ (see `pr#PREC2level`).
- `func pr#prfe#K[int] pr#fGMRES(pr#prfe#K[int]& x0i, pr#prfe#K[int]& bi, real eps, int nbiter, string sprec)` The function `pr#fGMRES` allows to solve the linear system $A x = b$ in parallel using the flexible GMRES method preconditioned by $M^{-1}$. The action of the global operator $A$ is given by `pr#A`, the action of the preconditioner $M^{-1}$ is given by `pr#PREC` and the scalar products are computed by `pr#scalprod`. More details are given in the section **[Solving the linear system](#solving-the-linear-system)**.

***For advanced users***:

NewMacro pr#localmacroaug pr#prfe#prmesh#buildAug EndMacro  
IFMACRO(pr#localmacroaug,1)  
matrix<pr#prfe#K> pr#CSinterp;  
ENDIFMACRO 

## Two level preconditioners

The main ingredient of a two level preconditioner is the so-called 'coarse space' matrix $Z$.  
$Z$ is a rectangular matrix of size $n \times n_c$, where usually $n_c \ll n$.  
$Z$ is used to build the 'coarse space operator' $E = Z^T A Z$, a square matrix of size $n_c \times n_c$. We can then define the 'coarse space correction operator' $Q = Z E^{-1} Z^T = Z (Z^T A Z)^{-1} Z^T$, which can then be used to enrich the one level preconditioner through a correction formula. The simplest one is the *additive* coarse correction: $M^{-1}_2 = M^{-1}_1 + Q$. See `pr#corr` below for all other available correction formulas.

There are multiple ways to define a relevant coarse space $Z$ for different classes of problems. **[`ffddmgeneosetup`](#building-the-geneo-coarse-space)** defines a coarse space correction operator by building the GenEO coarse space, while **[`ffddmEuansetup`](#building-the-coarse-space-from-a-coarse-mesh)** builds the coarse space using a coarse mesh.  
After a call to either **[`ffddmgeneosetup`](#building-the-geneo-coarse-space)** or **[`ffddmEuansetup`](#building-the-coarse-space-from-a-coarse-mesh)**, the following variables and functions are set up:

- `int pr#ncoarsespace` the size of the coarse space $n_c$.
- `string pr#corr` initialized with the value of [ffddmcorrection](parameters.md#global-parameters). Specifies the type of coarse correction formula to use for the two level preconditioner. The possible values are:  
$$
\begin{align}
&&\text{"AD"}:&&\textit{Additive}, \quad &M^{-1} = M^{-1}_2 = \phantom{(I - Q A) }M^{-1}_1\phantom{ (I - A Q)} + Q\\
&&\text{"BNN"}:&&\textit{Balancing Neumann-Neumann}, \quad &M^{-1} = M^{-1}_2 = (I - Q A) M^{-1}_1 (I - A Q) + Q\\
&&\text{"ADEF1"}:&&\textit{Adapted Deflation Variant 1}, \quad &M^{-1} = M^{-1}_2 = \phantom{(I - Q A) }M^{-1}_1 (I - A Q) + Q\\
&&\text{"ADEF2"}:&&\textit{Adapted Deflation Variant 2}, \quad &M^{-1} = M^{-1}_2 = (I - Q A) M^{-1}_1\phantom{ (I - A Q)} + Q\\
&&\text{"RBNN1"}:&&\textit{Reduced Balancing Variant 1}, \quad &M^{-1} = M^{-1}_2 = (I - Q A) M^{-1}_1 (I - A Q)\\
&&\text{"RBNN2"}:&&\textit{Reduced Balancing Variant 2}, \quad &M^{-1} = M^{-1}_2 = (I - Q A) M^{-1}_1\phantom{ (I - A Q)}\\
&&\text{"none"}:&&\textit{no coarse correction}, \quad &M^{-1} = M^{-1}_2 = \phantom{(I - Q A) }M^{-1}_1\phantom{ (I - A Q)}\\
\end{align}
$$
Note that *AD*, *ADEF1* and *RBNN2* only require one application of $Q$, while *BNN*, *ADEF2* and *RBNN1* require two. The default coarse correction is *ADEF1*, which is cheaper and generally as robust as *BNN* or *ADEF2*.

- `func pr#prfe#K[int] pr#Q(pr#prfe#K[int] &ui)` The function `pr#Q` computes the parallel application of the coarse correction operator $Q$, i.e. the action of $Q = Z E^{-1} Z^T$ on the local vector $u_i$. In the sequential case, it computes the action of $Q$ on a global vector.  
The implementation differs depending on the method used to build the coarse space (with GenEO or using a coarse mesh), but the idea is the same: the action of the transpose of the distributed operator $Z$ on the distributed vector $u_i$ is computed in parallel, with the contribution of all subdomains being gathered in a vector of size $n_c$ in the mpi process of rank 0. The action of the inverse of the coarse space operator $E$ is then computed by forward-backward substitution using its $LU$ (or $LDL^T$) decomposition previously computed by the first `pr#prfe#prmesh#pCS` ranks of the mpi communicator. The result is then sent back to all subdomains to perform the last application of $Z$ and obtain the resulting local vector in each subdomain.
- `func pr#prfe#K[int] pr#PREC2level(pr#prfe#K[int] &ui)` The function `pr#PREC2level` computes the parallel application of the two level preconditioner $M^{-1}_2$, i.e. the action of $M^{-1}_2$ on the local vector $u_i$. In the sequential case, it computes the action of $M^{-1}_2$ on a global vector. The two level preconditioner depends on the choice of the coarse correction formula which is determined by `pr#corr`, see above.


***For advanced users***:

int pr#bEuan = 0;

### Building the GenEO coarse space

```cpp
ffddmgeneosetup(pr,Varf)
```

This builds the GenEO coarse space for problem **pr**. This will create and expose variables whose names will be prefixed by **pr**, see below. It is assumed that **[`ffddmsetupPrecond`](#one-level-preconditioners)** has already been called for prefix **pr** in order to define the one level preconditioner for problem **pr**. The GenEO coarse space is $Z = (R_i^T D_i V_{i,k})^{i=1,...,N}_{\lambda_{i,k} \ge \tau}$, where $V_{i,k}$ are eigenvectors corresponding to eigenvalues $\lambda_{i,k}$ of the following local generalized eigenvalue problem in subdomain $i$:

$D_i A_i D_i V_{i,k} = \lambda_{i,k} A_i^{\text{Neu}} V_{i,k}$,  
 
where $A_i^{\text{Neu}}$ is the local Neumann matrix of subdomain $i$ (with Neumann boundary conditions at the subdomain interface).  
In practice, this builds and factorizes the local Neumann matrices $A_i^{\text{Neu}}$ corresponding to the abstract bilinear form given by the macro **Varf** (see **[`ffddmsetupOperator`](#define-the-problem-to-solve)** for more details on how to define the abstract variational form as a macro). In the GenEO method, the abstract bilinear form **Varf** is assumed to be the same as the one used to define the problem **pr** through the previous call to **[`ffddmsetupOperator`](#define-the-problem-to-solve)**. The local generalized eigenvalue problem is then solved in each subdomain to find the eigenvectors $V_{i,k}$ corresponding to the largest eigenvalues $\lambda_{i,k}$ (see `pr#Z` below). The number of computed eigenvectors $\nu$ is given by [ffddmnu](parameters.md#global-parameters). The eigenvectors selected to enter $Z$ correspond to eigenvalues $\lambda_{i,k}$ larger than $\tau$, where the threshold parameter $\tau$ is given by [ffddmtau](parameters.md#global-parameters). If [ffddmtau](parameters.md#global-parameters) $= 0$, all [ffddmnu](parameters.md#global-parameters) eigenvectors are selected.   Finally, the coarse space operator $E = Z^T A Z$ is assembled and factorized (see `pr#E` below).

***defines***:

- `pr#prfe#K[int][int] pr#Z` array of local eigenvectors $Z_{i,k} = D_i V_{i,k}$ obtained by solving the local generalized eigenvalue problem above in the subdomain of this mpi rank using *Arpack*. The number of computed eigenvectors $\nu$ is given by [ffddmnu](parameters.md#global-parameters). The eigenvectors selected to enter $Z$ correspond to eigenvalues $\lambda_{i,k}$ larger than $\tau$, where the threshold parameter $\tau$ is given by [ffddmtau](parameters.md#global-parameters). If [ffddmtau](parameters.md#global-parameters) $= 0$, all [ffddmnu](parameters.md#global-parameters) eigenvectors are selected.

- `matrix<pr#prfe#K> pr#E` the coarse space operator $E = Z^T A Z$. The matrix `pr#E` is assembled in parallel and is factorized by the parallel direct solver *MUMPS* using the first `pr#prfe#prmesh#pCS` ranks of the mpi communicator, with mpi rank 0 as the master process. The number of mpi processes dedicated to the coarse problem is set by the underlying mesh decomposition of problem **pr**, which also specifies if these mpi ranks are excluded from the spatial decomposition or not. These parameters are set by [ffddmpCS](parameters.md#global-parameters) and [ffddmexclude](parameters.md#global-parameters) when calling **[`ffddmbuildDmesh`](#overlapping-mesh-decomposition)** (see **[`ffddmbuildDmesh`](#overlapping-mesh-decomposition)** for more details).

***For advanced users***:

int pr#si;

pr#sizelg(pr#prfe#prmesh#npart), pr#offseti(pr#prfe#prmesh#npart);

int[int] pr#sizelgworld(mpiSize(pr#prfe#prmesh#mpicomm)), pr#offsetiworld(mpiSize(pr#prfe#prmesh#mpicomm));

matrix<pr#prfe#K> pr#matN;

### Building the coarse space from a coarse mesh

```cpp
ffddmEuansetup(pr,Thc,VarfEprec,VarfA)
```

builds the coarse space for problem **pr** from a coarse problem which corresponds to the discretization of a variational form on a coarser mesh **Thc** of $\Omega$. This will create and expose variables whose names will be prefixed by **pr**, see below. It is assumed that **[`ffddmsetupPrecond`](#one-level-preconditioners)** has already been called for prefix **pr** in order to define the one level preconditioner for problem **pr**. The abstract variational form for the coarse problem can differ from the original problem **pr** and is given by macro **VarfEprec** (see **[`ffddmsetupOperator`](#define-the-problem-to-solve)** for more details on how to define the abstract variational form as a macro).  
The coarse space $Z$ corresponds to the interpolation operator from the 

***defines***:

matrix<pr#prfe#K> pr#AglobEprec;

matrix<pr#prfe#K>[int] pr#aRdEprec(pr#prfe#prmesh#npart);

func pr#prfe#K[int] pr#AEprec(pr#prfe#K[int] &x)

matrix<pr#prfe#K> pr#Zeuan

matrix<pr#prfe#K> pr#Zeuani

matrix<pr#prfe#K> pr#Eeuan;

## Solving the linear system

```cpp
func pr#prfe#K[int] pr#fGMRES(pr#prfe#K[int]& x0i, pr#prfe#K[int]& bi, real eps, int nbiter, string sprec)
```
