[`test`](documentation.md#test)

## Command-line arguments

- `-ffddm_verbosity N`, the level of verbosity of **ffddm**, see [ffddmverbosity](#global-parameters) (default 3).
- `-seqddm N` use **ffddm** in sequential mode, with N the number of subdomains.
- `-noGlob` if present, do not define any global quantity (such as saving the global mesh for plotting or building the global restriction matrices). Cannot be used in sequential mode or with plotting.
- `-ffddm_partitioner N` specifies how to partition the initial domain, see [ffddmpartitioner](#global-parameters) (default 1, *metis*).
- `-ffddm_overlap N` specifies the width of the overlap region between subdomains, see [ffddmoverlap](#global-parameters) (default 1). 
- `-ffddm_master_p N`, number of master processes for the coarse problem (for two level preconditioners), see [ffddmpCS](#global-parameters) (default 1).
- `-ffddm_master_exclude 0|1` exclude master processes from the domain decomposition, see [ffddmexclude](#global-parameters) (default 0).
- `ffddm_split N`, level of refinement of the local submeshes with respect to the initial global mesh, see [ffddmsplit](#global-parameters) (default 1).
- `-ffddm_schwarz_method S`, specifies the type of one level preconditioner $M^{-1}_1$: "asm" (*Additive Schwarz*), "ras" (*Restricted Additive Schwarz*), "oras" (*Optimized Restricted Additive Schwarz*), "soras" (*Symmetric Optimized Restricted Additive Schwarz*) or "none" (no preconditioner), see [ffddmprecond](#global-parameters) (default "ras").
- `-ffddm_geneo_nu N`, number of local eigenvectors to compute in each subdomain when solving the local generalized eigenvalue problem for the GenEO method, see [ffddmnu](#global-parameters) (default 20).
- `-ffddm_geneo_threshold R`, threshold parameter for selecting local eigenvectors when solving the local generalized eigenvalue problems for the GenEO method, see [ffddmtau](#global-parameters) (default 0.5). If the command-line parameter **-ffddm_geneo_nu N** is used, then [ffddmtau](#global-parameters) is initialized to 0.
- `-ffddm_schwarz_coarse_correction S`, specifies the coarse correction formula to use for the two level preconditioner: "AD" (*Additive*), "BNN" (*Balancing Neumann-Neumann*), "ADEF1" (*Adapted Deflation Variant 1*), "ADEF2" (*Adapted Deflation Variant 2*), "RBNN1" (*Reduced Balancing Variant 1*), "RBNN2" (*Reduced Balancing Variant 2*) or "none" (no coarse correction), see [ffddmcorrection](#global-parameters) (default "ADEF1").

## Global parameters

- `ffddmverbosity` initialized by command-line argument **-ffddm_verbosity N**, specifies the level of verbosity of **ffddm** (default 3).
- `ffddmpartitioner` initialized by command-line argument **-ffddm_partitioner N**, specifies how to partition the initial domain:
	- N=0: user-defined partition through the definition of a macro, see **[`ffddmbuildDmesh`](documentation.md#overlapping-mesh-decomposition)**
	- N=1: use the automatic graph partitioner *metis* (default)
	- N=2: use the automatic graph partitioner *scotch*
- `ffddmoverlap` initialized by command-line argument **-ffddm_overlap N**, specifies the number of layers of mesh elements in the overlap region between subdomains N >= 1 (default 1). ***Remark*** The actual width of the overlap region between subdomains is 2N, since each subdomain is extended by N layers of elements in a symmetric way.
- `ffddminterfacelabel` the label of the new border of the subdomain meshes (the interface between the subdomains) (default 10). Used for imposing problem-dependent boundary conditions at the interface between subdomains for the preconditioner, for example optimized Robin boundary conditions (see ORAS).
- `ffddmpCS` initialized by command-line argument **-ffddm_master_p N**, number of mpi processes used for the assembly and resolution of the coarse problem for two level preconditioners (default 1).
- `ffddmexclude` initialized by command-line argument **-ffddm_master_exclude**, 0 or 1 (default 0). If true, mpi ranks participating in the assembly and resolution of the coarse problem for two level preconditioners will be excluded from the spatial domain decomposition and will only work on the coarse problem.
- `ffddmsplit` initialized by command-line argument **ffddm_split N**, level of refinement of the local submeshes with respect to the initial global mesh (default 1). This is useful for large problems, where we want to avoid working with a very large global mesh. The idea is to start from a coarser global mesh, and generate finer local meshes in parallel during the mesh decomposition step in order to reach the desired level of refinement for the subdomains. For example, calling **[`ffddmbuildDmesh`](documentation.md#overlapping-mesh-decomposition)** with [ffddmsplit](#global-parameters) = 3 will generate local submeshes where each mesh element of the initial mesh is split into $3^d$ elements.
- `ffddmprecond` initialized by command-line argument **-ffddm_schwarz_method S**, specifies the type of one level preconditioner $M^{-1}_1$ to build when calling **[`ffddmsetupPrecond`](documentation.md#one-level-preconditioners)**: "asm" (*Additive Schwarz*), "ras" (*Restricted Additive Schwarz*), "oras" (*Optimized Restricted Additive Schwarz*), "soras" (*Symmetric Optimized Restricted Additive Schwarz*) or "none" (no preconditioner). Default is "ras". See **[`ffddmsetupPrecond`](documentation.md#one-level-preconditioners)** for more details.
- `ffddmnu` initialized by command-line argument **-ffddm_geneo_nu N**, number of local eigenvectors to compute in each subdomain when solving the local generalized eigenvalue problem for the GenEO method (default 20). See **[`ffddmgeneosetup`](documentation.md#building-the-geneo-coarse-space)** for more details.
- `ffddmtau` initialized by command-line argument **-ffddm_geneo_threshold R**, threshold parameter for selecting local eigenvectors when solving the local generalized eigenvalue problems for the GenEO method (default 0.5). If the command-line parameter **-ffddm_geneo_nu N** is used, then [ffddmtau](#global-parameters) is initialized to 0. See **[`ffddmgeneosetup`](documentation.md#building-the-geneo-coarse-space)** for more details.
- `ffddmcorrection` initialized by command-line argument **-ffddm_schwarz_coarse_correction S**, specifies the coarse correction formula to use for the two level preconditioner: "AD" (*Additive*), "BNN" (*Balancing Neumann-Neumann*), "ADEF1" (*Adapted Deflation Variant 1*), "ADEF2" (*Adapted Deflation Variant 2*), "RBNN1" (*Reduced Balancing Variant 1*), "RBNN2" (*Reduced Balancing Variant 2*) or "none" (no coarse correction). Default is "ADEF1". See the section about **[Two level preconditioners](documentation.md#two-level-preconditioners)** for more details.


