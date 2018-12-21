

# Domain Decomposition (DD)
When the size of a three dimensional problem is large (whatever it means), it is necessary to distribute data among several processors especially for solving linear systems. A natural way is to do it via domain decomposition.


## Mesh Decomposition 
<!--- INSERER FIGURE? The overlap is the minimum width of the overlap between sub-meshes.
implicit  global renvoie a implicit en fait --->
The starting point is a collection of $N$ sub-meshes $(Th_i)_{i=1}^N$ that together form a global mesh

$$Th:= \cup_{i=1}^N Th_i\,.$$
 This induces a natural decomposition of the global finite element space $Vh$ on $Th$ into $N$ local finite element spaces $(Vh_i)_{i=1}^N$ each of them defined on $Th_i$.  

 **Note** By global, we mean that the corresponding structure can be refered to in the code (most often only) by its local values. In computer science term, it corresponds to a distributed data where each piece of data is stored by a MPI process.


## Distributed Linear Algebra
The domain decomposition induces a natural decomposition of the set of the global degrees of freedom (d.o.f.) ${\mathcal N}$ of the finite element space $Vh$ into the $N$ subsets of d.o.f.'s $({\mathcal N})_{i=1}^N$ each associated with the local finite element space $Vh_i$. We have thus 

$${\mathcal N} = \cup_{i=1}^N {\mathcal N}_i\,, $$
but with duplications of some of the d.o.f.'s.

 Associated with this decomposition of the set of d.o.f.'s ${\mathcal N}$, a *distributed vector* is a collection of local vectors $({\mathbf V_i}_{1\le i\le N})$ so that the values on the duplicated d.o.f.'s are the same.


**Note**  In mathematical terms, it can be described as follows for a real valued problem. For a real value problem, simply replace $\R$ with $\C$. Let $R_i$ be the restriction operator from $\R^{\#{\mathcal N}}$ to $\R^{\#{\mathcal N}_i}$, where $\#{\mathcal N}_i$ denotes the number of elements of ${\mathcal N}_i$. A collection of local vectors $({\mathbf V}_i)_{1\le i\le N}\in \Pi_{i=1}^N \R^{\#{\mathcal N}_i}$ is a distributed vector iff there exists a global vector ${\mathbf V}\in\R^{\#{\mathcal N}}$ such that for all subset $1\le i\le N$, we have:

$$
{\mathbf V}_i = R_i\,{\mathbf V}\,.
$$
We will also say that the collection of local vectors $({\mathbf V}_i)_{1\le i\le N}$ is consistent.

## Partition of Unity Matrices (POUM)
 Let $(D_i)_{1\le i \le N}$ be square diagonal matrices of size $\#{\mathcal N}_i$ which form a partition of unity in the sense that:

$$
  Id_{} = \sum_{i=1}^N R_i^T\,D_i\,R_i\text{ in }\R^{\#{\mathcal N}\times \#{\mathcal N}} \,.
$$
For instance if a degree of freedom is shared by $k$ subdomains defining the corresponding entry of the diagonal matrix $D$ to be $1/k$ yields partition of unity matrices. The matrices $R_i$ and $D_i$ are the heart of distributed linear algebra.

### Distributed scalar product
For two global vectors ${\mathbf U}$ and ${\mathbf V}$ of size $\#{\mathcal N}$, the formula for the scalar product ${\mathbf V}^T\,{\mathbf U}=({\mathbf U},\,{\mathbf V})$ in terms of their distributed vector counterparts is:

$$({\mathbf U}, {\mathbf V}) = \left({\mathbf U}, \sum_{i=1}^N R_i^T D_i R_i {\mathbf V}\right) = \sum_{i=1}^N(R_i {\mathbf U}, D_i R_i {\mathbf V})
=\sum_{i=1}^N\left({\mathbf U}_i, D_i {\mathbf V}_i\right)\,.
$$

Local scalar products are performed concurrently. Thus, the implementation is parallel except for the sum which corresponds to a MPI\_Reduce call across the $N$ MPI processes. Note also that the implementation relies on the knowledge of a partition of unity so that the FreeFem++ syntax is `dscalprod(Di,u,v)` or equivalently  `pr#scalprod(u,v)` where `pr` is a user defined prefix that refers to the domain decomposition and thus implicitely also to the partition of unity.

### Update
From a collection of local vectors $({\mathbf U}_i)_{1\le i \le N}$, it is possible ensure consistency of the duplicated data and thus creating a distributed vector $({\mathbf V}_i)_{1\le i \le N}$ by calling the function `pr#update(Ui, TRUE)` where `pr` is a user defined prefix that refers to the domain decomposition. This function performs the following operation for all $1\le i \le N$:

$$
 {\mathbf V}_i \leftarrow R_i\, \sum_{j=1}^N R_j^T D_j {\mathbf U}_j
$$

**Note** The implementation corresponds to

$$
 {\mathbf V}_i := R_i \sum_{j=1}^N R_j^T D_j {\mathbf U}_j = D_i {\mathbf U}_i + \sum_{j\in \mathcal{O}(i)} R_i\,R_j^T\,D_j {\mathbf U}_j
$$
where $\mathcal{O}(i)$ is the set of neighbors of subdomain $i$. Therefore, the matrix vector product is computed in three steps:
- concurrent computing of $D_j {\mathbf U}_j$ for all $1\le j\le N$;
- neighbor to neighbor MPI-communications ($R_i\,R_j^T$) ;
- concurrent sum of neighbor contributions.


## Distributed Matrix and Vector resulting from a variational formulation   
The discretization of a variational formulation on the global mesh $Th$ yields a global matrix $A$ and a global right hand side $\mathbf{RHS}$. Thanks to the sparsity of finite element matrices for partial differential equations and thanks to the overlap between subdomains, the knowledge of the local matrix $R_i A R_i^T$ on each subdomain $1\le i\le N$ is sufficient to perform the matrix-vector product $A\times \mathbf{U}$ for any global vector $\mathbf{U}$. Once the problem has been set up by a call to `ffddmsetupOperator(myprefix,myFEprefix,myVarf)`, the matrix-vector product is performed by calling the function `pr#A(Ui)` where `pr` is a user defined prefix that refers to the problem at hand which itself implicitly refers to the triplet (domain decomposition, finite element, variational formulation). See more on problem defintion in this [documentation](./documentation.md#Define-the-problem-to-solve) and more on distributed linear algebra in chapter 8 of ["An Introduction to Domain Decomposition Methods: algorithms, theory and parallel implementation" SIAM 2015](http://bookstore.siam.org/ot144/) or in the [appendix](#Appendix) below (A FAIRE OU NON? ).


## Distributed Linear Solvers
In many cases, we are interested in the solution of the problem  in terms of the vector of d.o.f.'s $\mathbf{X}$ that satisfies:

$$A\, \mathbf{X} = \mathbf{RHS}\,.$$

### Distributed Direct Solvers

A COMPLETER

### Schwarz methods

We consider the solve of the equation $A\, \mathbf{X} = \mathbf{RHS}$ by a flexible GMRES method preconditioned by domain decomposition methods.

#### Restricted Additive Schwarz (RAS)
The RAS preconditioner reads:

$$
M^{-1}_{RAS} := \sum_{j=1}^N R_j^T D_j (R_j\, A\,R_j^T)^{-1} R_j\,.
$$

Let $A_{i}$ denote the local matrix $(R_i\, A\,R_i^T)$. The application of the operator $M^{-1}_{RAS}$ to a distributed right hand side $(\mathbf{RHS}_i)_{i=1}^N$ consists in computing:

$$
R_i\, \sum_{j=1}^N R_j^T\,D_j\, A_{j}^{-1}\,\, \mathbf{ RHS}_j
= D_i\, A_{i}^{-1}\, \mathbf{ RHS}_i + \sum_{j\in \mathcal{O}(i)} (R_i\,R_j^T)\,D_j\, A_{j}^{-1}\, \mathbf{ RHS}_j\,.
$$

This task is performed by first solving concurrently on all subdomains a linear system for ${\mathbf Y}_j$ for all $1\le j \le N$:

$$
A_{j}\, {\mathbf Y}_j = \mathbf{RHS}_j\,.
$$

Each local vector ${\mathbf Y}_j$ is weighted by the partition of unity matrix $D_j$. Then data transfers between neighboring subdomains implement the $R_i\,R_j^T\,D_j\,{\mathbf Y}_j$ formula. The contribution from neighboring subdomains are summed locally. This pattern is very similar to that of the [update](#Update) procedure.

#### Two level methods
The RAS method is called a one-level method in the sense that sub-domains only interact with their direct neighbors. For some problems such as Darcy problems or static elasticiy problems and when the number of subdomains is large, such one-level methods may suffer from a slow convergence. The fix is to add to the preconditioner an auxiliary coarse problem that couples all subdomains at each iteration and is inexpensive to calculate. We consider two ways to build this coarse problem, see below [Coarse Mesh](#coarse-mesh) and [GenEO](#geneo)


##### Coarse Mesh
A first possibility is to discretize the problem on a coarse mesh, following the same principle as multi-grid methods. For 3-D problems, a coarsening of the mesh size by a factor 2, reduces by a factor $2^3=8$ the size of the coarse problem which is then easier to solve by a direct method.

##### GenEO
For highly heterogeneous or anisotropic problems, two level methods based on coarse meshes might fail and a more sophisticated construction must be used. A provable robust coarse space called GenEO is built by first solving the following local generalized eigenvalue problem in parallel for each subdomain $1\le i\le N$, where $A_i^{\text{Neu}}$ denotes the local matrix resulting from the variational formulation:

$$
D_i A_i D_i\, V_{i,k} = \lambda_{i,k}\, A_i^{\text{Neu}} \,V_{i,k}
$$
The eigenvectors selected to enter the coarse space correspond to eigenvalues $\lambda_{i,k} \ge \tau$, where the threshold parameter $\tau$ is user-defined. The precise formulas are given in this [documentation](./documentation.md#building-the-geneo-coarse-space). From a mathematical point of view, it has been proved that for a symmetric positive definite matrix $A$, the spectrum of the  preconditioned by the two-level method with a GenEO coarse space lies in the interval $[\displaystyle \frac{1}{1+k_1\,\tau} , k_0 ]$.

**Note** A heuristic that justifies this construction is as follows. We first introduce the Additive Schwarz method (ASM) which can be seen as a symmetrized variant of the RAS preconditioner:

$$
	M_{ASM}^{-1} := \sum_{j=1}^N R_j^T A_j^{-1} R_j\,.
$$
It can be proved that the lower bound for the eigenvalue of $M_{ASM}^{-1}\,A$ is close to zero (which is bad for convergence) whereas the upper bound depends only on the number of neigbors of a subdomain (which is good for convergence).

Second, we also introduce the following preconditioner $M^{-1}_{NN}$:

$$
	M^{-1}_{NN} := \sum_{1\le j\le N} D_i\,(A_j^{\text{Neu}})^{-1} D_j\,.
$$
We have a very good lower bound for the preconditioned operator $M^{-1}_{NN}\,A$ that does not depend on the number of subdomains  but only on the maximum multiplicity of intersections $k_1$  (which is good for convergence). But the upper bound for this preconditioner is very large (which is bad for convergence).

Now, if we compare formulas for $M^{-1}_{NN}$ and $M^{-1}_{ASM}$, we may suspect that vectors $\mathbf{V}_{ik}$ for which $D_i\, (A_i^{\text{Neu}})^{-1}\,D_i\,\mathbf{V}_{ik}$ and $A_{i}^{-1}\,\mathbf{V}_{ik}$ have very different values are responsible for the slow convergence and should contribute to the coarse space. This is a way to interpret the above generalized eigenvalue problem which controls the lower bound of the two-level preconditioned system.

# exemple sequentiel
Est ce sa place ici?????
puis detail de comment passer du code sequentiel au code parallele . Doc auto contenu sans regarder les fichiers sources


# Appendix

Reprendre le livre

## Distributed linear algebra




@book{Dolean:2015:IDDSiam,
	Author = {Dolean, Victorita and Jolivet, Pierre and Nataf, Fr\'ed\'eric},
	Date-Added = {2015-12-24 16:49:03 +0000},
	Date-Modified = {2016-01-29 10:15:37 +0000},
	Publisher = {SIAM},
	Title = {

bibliography: bookddm.bib
