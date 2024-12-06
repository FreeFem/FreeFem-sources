---
name: Lac
category: mesh
folder: 3d-examples
---
# Solution of the Radiative Transfer equations in the Chamonix valley

### Companion script of 2 articles published in JCP
	Golse, F., Hecht, F., Pironneau, O., Tournier, P. H., & Smets, D. (2022). Radiative Transfer For Variable 3D Atmospheres.
 arXiv preprint arXiv:2208.06410.  
 https://doi.org/10.1016/j.jcp.2022.111864  
 https://doi.org/10.1016/j.jcp.2023.112531

This script also illustrates the use of the Htool Hierarchical Matrix library to compress a dense matrix stemming from a custom user-defined operator defined in a user C++ FreeFEM plugin. The script uses freefem-mpi to run in parallel.

### About Htool
Htool (https://github.com/htool-ddm/htool) is designed for electromagnetism, acoustics, etc, by boundary element methods.  
Interface with Htool for a new problem, with a new kernel like here, can be done by defining and interfacing a Generator class deriving from the virtual class VirtualGenerator of Htool. This is done here in the freefem plugin `RadiativeTransfer_htool.cpp`.

The file can be compiled by typing in a terminal window

	ff-c++ -auto RadiativeTransfer_htool.cpp

### About RadiativeTransfer_htool.cpp
The "Generator" class contains the description of the kernel; it is created through a user-defined constructor, for example here GenSE which takes for arguments a mesh3, a meshS, two vectors and a new user-defined type KappaGrid (also defined in the plugin); thence, to use it:

	Generator GenSE(Th3,onGamma,ThS,seeface,kappag);

The Generator can then be passed to Htool through the Build() function, which returns the compressed H-Matrix built by Htool from the user-defined copy_submatrix (or get_coef) routine  of the Generator computing the matrix coefficients:

	HMatrix HSE = Build(GenSE,VhS,Vh,eta=100,eps=1e-2,minclustersize=10);

### To launch this script:
for example:

~~~bash
ff-mpirun -np 8 chamonix3.edp -ns -wg -Kmax 1 -fullmesh -n 15
~~~

which means run with 8 cores, no script output, with graphics, with one level of $\kappa$ (grey case) and a mesh refinement defined by n=15.

In order to run the test cases from the JCP article you can select  the chamonix smaller mesh or you can download the original mesh at https://freefem.org/misc/chamonix_100K_surf.meshb,
 load the file in the readmesh command, and run with argument '-fullmesh':
~~~bash
	wget https://freefem.org/misc/chamonix_100K_surf.meshb
~~~
then run with
~~~bash
	 ff-mpirun -np 8 chamonix3.edp -ns -wg -Kmax 9
~~~

### Remark
You are reading this Markdown file which contains explainations and the FreeFEM script. When running Markdown files, FreeFEM selects the code lines enclosed in Markdown FreeFEM code blocks.

For a user friendly environment install the **vscode-FreeFEM** extension in the **VSCode** app. If your code is "not trusted" do Ctrl+Maj+P and look for "trust workspace settings".

## The problem solved
$$
\begin{aligned}
{}&\omega\cdot\nabla I_\nu+{\rho\kappa} I_\nu={\rho\kappa}(1-a)B_\nu(T)+{\rho\kappa} a J_\nu\,,\forall\omega\in S^2,~\forall \nu\in (\nu_{min},\nu_{max})
\\
&I_\nu({x},\omega)\!=\!Q_\nu({x},\omega)\,,\quad\forall\omega\cdot n({x})<0\,,\,\,\forall {x}\in\partial\Omega.
\\
&u\nabla T -\nabla\cdot(\kappa_\theta\nabla T)+\int_{\nu_{min}}^{\nu_{max}}\rho\kappa(1-a)(J_\nu-B_\nu(T))d\nu\,=0,\quad  \forall x\in\Omega\subset R^3,
\\ & \hbox{ where } J_\nu:=\tfrac1{4\pi}\int_{{S^2}}I_\nu d\omega\,,
\end{aligned}
$$
$B_\nu=\nu^3/(e^{-\frac\nu T}-1)$ is the normalized Planck function with $T= T_{Kelvin}/4798$.  The frequencies are also normalized: the physical frequencies are divided by $10^{14}$.

The unknowns are the radiation intensity $(\nu,\omega,x)\mapsto I_\nu(\omega,x)$ and the temperature $x\mapsto T(x)$, $\omega\in S^2$ are the directions of the rays in the unit sphere $S^2$, and $\nu$ is its frequency.  The density $\rho$ is a given function of $x$, usually linear, exponential or constant. To define your own you must select the "constant" predefined type; $\kappa$ is a given function of $\nu$ read from the file $\texttt{kappafile}$. The boundary conditions are applied at points and rays entering the domain ($n$ is the outer normal).

The program handles $u\ne 0$, $\kappa_\theta\ne 0$, but we document only the zero case.

Two cases are of interest:
  1. Infrared radiations enter $\Omega$ from the ground and  no rays enter from the upper atmosphere.
  2. Sunrays enter $\Omega$ from the the top $z=Z$, crosses the atmosphere undisturbed and reflects on the ground. So here too one assumes that the rays are emitted by the ground and no rays come from the top.
  The part of $\partial\Omega$ where the rays enter is called $\Sigma$.

#### Reformulation of the first PDE
If $S$ is the right-hand side of the first PDE, the method of characteristics gives
$$
\begin{aligned}
J_\nu(x)&:=\frac1{4\pi}\int_{S^2} I(x,\omega){d}\omega
=\frac1{4\pi}\int_{S^2} I({x}_\Sigma({x},\omega))e^{-\int_0^{\tau_{{x},\omega}}\kappa({x}-\omega s){d} s}{d}\omega
\cr& + \frac1{4\pi}\int_{S^2}\int_0^{\tau_{{x},\omega}} e^{-\int_0^s\kappa({x}-\omega s'){d} s'} S({x}-\omega s){d} s{d}\omega
%\cr&
=: S^E_\nu({x}) + {\mathcal J}[S]({x}),
\end{aligned}
$$
where
$$
\begin{aligned}&
 S^E_\nu({x}):=
 \frac{1}{4\pi}\int_\Gamma  Q_\nu({x}',\frac{{x}'-{x}}{|{x}'-{x}|})
 \frac{[({x}'-{x})\cdot n({x}')]_+}{|{x}'-{x}|^3} e^{-\int_{[{x},{x}']}\kappa}{d}\Gamma({x}'),
 \cr&
{\mathcal J}_\nu[S]({x}):=
 \frac1{4\pi}\int_{\Omega} S({x}')\frac{e^{-\int_{[{x},{x}']}\kappa}}{|{x}'-{x}|^2}{d}{x}'.
 \end{aligned}
$$

### The iterative method

If $T$  and an old $J_\nu$ are given, then we can compute $S$ and then compute a new $J_\nu$ for all $\nu,x$ and then finally update $T(x)$ from the integral equation of the problem.
The integral equation is solved by a Newton method at every $x$.

To compute the two integrals the integrands $S$ and $Q$ are approximated by  continuous piecewise polynomials on a tetraedral (resp. triangular) mesh of $\Omega$ (resp. $\Sigma$). Then, one computes once and for all $S_\nu^{ij}$ which is $S_\nu^E(x)$ when $Q_\nu$ is replaced by the finite element hat function associated with vertex $x^i$ and similarly for ${\mathcal J}_\nu$, producing ${\mathcal J}_\nu^{ij}$ .  Notice that $S_\nu$ is a rectangular matrix while $ {\mathcal J}_\nu$ is a square matrix.

These matrices are stored as ${\mathcal H}$-matrices to reduce the computing time to $O(n\log n)$ where $n$ is the number of vertices in the mesh.

## The freefem script
The following libraries are used

~~~freefem
// NBPROC 8
// PARAM -ns -wg -Kmax 9
load "bem"  // boundary element methods for electromagnetism (Htool)
load "RadiativeTransfer_htool"  // contains the functions specific to this application program in C++
load "qf11to25"  // needed for better quadrature rules of integrals
load "medit"  // for 3D graphics
load "iovtk"  // interface to vtk
load "mmg"  // for mesh adaptivity
load "shell"  // interface with unix terminal
load "gsl" // gnu scientific library
load "ffrandom" // for random functions
srandomdev(); // for a random def of kappa and a

include "getARGV.idp" // to specify parameters in the terminal command

verbosity=0;  // minimal
~~~

There are two types of parameters to define: physical quantities and algorithmic parameters.

### Physical Parameters
- The light source from the ground is
  $$
  Q_\nu(x,\omega) =Q_0(\beta+(1-\beta){\bf 1}_{z<s_{al}}) \cos(\omega\cdot n)B_\nu(T_{sun})\cos\gamma
  $$
- $\gamma$ is the angle of incidence of the light source (to study the difference between morning and evening)
- $Q_0=2.10^{-5}$ is its intensity (to understand why the sun light comes from the ground, see the paper)
- $T_{sun}=1.209$ is the scaled sun temperature (5800K before rescaling)
- $ (\beta+(1-\beta){\bf 1}_{z<s_{al}})$ is a factor to account for snow on the ground above altitude $s_{al}$
- The scattering coefficient depends on the presence of a cloud, assumed to be in a band $x\in[1.5-\sqrt{0.5},1.5+\sqrt{0.5}]$ between altitude $cloud_m$ and $cloud_M$.  It is possible to have a random density/scattering in the cloud.
 Here the scattering coefficient is $a_{scat}{\bf 1}_{z\in(cloud_m,cloud_M)}{\bf 1}_{(x-1.5)^2<0.5}$.

### Algorithmic Parameters
- `Newton` is the maximum number of Newton iterations for the temperature equation
- `heat`: the temperature equation is the heat PDE with thermal diffusion
- `Niter` is the number of global iterations $T\to I\to T$

~~~freefem
/********************* SCRIPT PARAMETERS *********************/
int dc = 1, 						// samples every dc points in the nu integrals (speedup with dc>1, danger, do not change)
Newton = 40,						// max Newton iter for the temp eq
n = getARGV("-n", 13),				// controls the number of vertices in the mesh
Niter = 7,  						// outer loop for convergence
Kmax = getARGV("-Kmax",1), 		    // 1 if kappa constant // in [2,10] otherwise
gamma = getARGV("-gamma", 0), 		// incidence of sunlight (0 is vertical, 1 is morning and -1 evening)
heat = 0, 							// if 1, temperature eq is solved
nosnow = getARGV("-nosnow", 0), 	// 0 or 1 if no snow
cloud = getARGV("-cloud",0), 		// 0 or 1 if cloud
altkappa = getARGV("-altkappa", 0);	// if 1 (resp-1) Kappa is changed in a range for CO2 (resp CH4)

real sal = 10*nosnow+0.25, beta=0.3, 	// snow altitude limit and snow albedo
	Tsun=1.209, Q0=2.0e-5, sigma=sqr(sqr(pi))/15; // the Sunlight

real ascat =0, cloudm=0.3, cloudM=0.7; // scattering in altitude [zscalm,zscalM]
~~~

### Names of I/O files

There are two input files and several output files.

One input file, `geminitransmittance.txt`, defines $\nu\mapsto\kappa_\nu$ pointwise, the other defines `Thf` the mesh and altitude of the ground (see `readmeshS`).

All output file paths are relative to the `basedir`; their names reflect the test case. Hence output files from one case do not overwrite outputs of another cases if the cases are different.

~~~freefem
/********************* OUTPUT FILES *********************/
string basedir("chamonix-output/");
if(mpirank==0) mkdir(basedir);
string tempestr("noon");
if (gamma == 1) tempestr = "morning";
if (gamma == -1) tempestr = "evening";
if (nosnow) tempestr = tempestr + "nosnow";
if (cloud) tempestr = tempestr + "cloud";
if (ascat>0)  tempestr = tempestr + "Scatt";
if (heat) tempestr = tempestr + "heat";
if (Kmax > 1) tempestr = tempestr + "K";  // res file non grey have a xx.K
if (Kmax > 1 && altkappa>0) tempestr = tempestr + "K";  // res file have a xx.KK
if (Kmax > 1 && altkappa<0) tempestr = tempestr + "KK"; // res file have a xx.KKK
int mem = storageused(); // memory used at this stage
~~~

### Read the mesh which defines the ground

The file which defines the mesh must be next to this .md file.  For the format of a surface mesh see the freefem documentation.
It contains a triangulation where the vertices have 3 coodinates (x,y,z).

To construct the 3d mesh above `Thf` it is convenient to use an auxiliary mesh which is `Thf` where all vertices are at altitude zero.

~~~freefem
/********************* BUILD THE MESH *********************/
mesh3 Th3;
mesh Th2;
if (mpirank == 0)
{
	// read the surface mesh of the Chamonix valley:
	meshS Thf = readmeshS("chamonix_surf.meshb");
	//meshS Thf = readmeshS("chamonix_100K_surf.meshb"); // uncomment to use the original fine mesh ; you can download it at https://freefem.org/misc/chamonix_100K_surf.meshb
	// build the volume mesh above the valley:
	real[int] bb(6);
	boundingbox(Thf, bb);
	int nn = 5*n;

	Th2 = movemesh(Thf,transfo=[x,y]); // flatten the mesh, will be the input of buildlayers
	fespace Vhf(Thf, P1);
	fespace Vh2(Th2, P1);
	Vhf zzf = z;
	Vh2 z2;
	z2[] = zzf[]; // get the altitude as a P1 function on the flattened 2D mesh
	
	bool fullmesh = usedARGV("-fullmesh") != -1; // use the original fine adapted surface mesh ; beware, expensive case !
	if (!fullmesh) Th2 = square(nn, nn, [bb[0]+x*(bb[1]-bb[0]),bb[2]+y*(bb[3]-bb[2])]);
	z2 = z2;

	real zmin = z2[].min, zmax = 1;
	plot(z2);

	int[int] l23 = [1,1,2,2,3,3,4,4], lup = [0,6], lbot = [0,5];
	real hx = (bb[1]-bb[0])/nn;
	int nnz = fullmesh ? 5 : (zmax-zmin)/hx;

	func cc = (zmax-z2)/(zmax-zmin);
	Th3 = buildlayers(Th2, coef=cc, nnz, zbound=[z2,zmax],
		labelmid=l23, labelup=lup, labeldown=lbot);
	
	int[int] reqt=[5];
	
	if (fullmesh) Th3 = mmg3d(Th3,requiredTriangle=reqt,hmin=hx,hmax=hx,hgrad=1.5,mem=8192);
}
~~~

### Define the finite element spaces

The finite element space of degree 1 on the volume mesh (tetrahedra) is called `Vh`.
The finite element space of degree 1 on the surface mesh (triangles) is called `VhS`.
It contains all radiative boundaries of the `Th` (in this case, the ground).

~~~freefem
broadcast(processor(0),Th3);
broadcast(processor(0),Th2);

real[int] bounds(6); // Domain bounds
boundingbox(Th3, bounds);
if (mpirank == 0) cout<< "x,y,z bounds: " << bounds << endl;

int[int] labRT = [5]; // emitting boundary
meshS ThS = extract(Th3, label=labRT, angle=pi); // extract as a surface mesh
if (mpirank == 0) cout << "Number of Vertices = " << Th3.nv << endl;
fespace Vh(Th3, P1);
fespace VhS(ThS, P1);
~~~

### Define the absorption $(\nu,x)\mapsto\rho(x)\kappa_\nu$.

The absorption is the product of the density $\rho(x)$ and a function $\nu\mapsto\kappa_\nu$.

For $\kappa_\nu$ we rely on a file read from the Gemini experiment web site. The number of lines in the file must be specified (more convenient than an EOF test). The file has 2 columns, one for $\nu[k]$ and one for the transmittance $1-\kappa[k]$ which is capped in $(0,1)$. The last 200 values are entered analytically because the data from Gemini are not complete.

To account for greenhouse gases we change the gemini data in given $\nu$-ranges, according to the value of `altkappa`.

~~~freefem
/********************* READ KAPPA[nu] *********************/
int Nnu = 483; // number of frequencies, will be augmented by 200
real[int] nu(200+Nnu), kappa(200+Nnu);
ifstream kappafile("geminitransmittance.txt");
int j =- 1;
while (j++ < Nnu-1) {
	real auxx, kappaux, nuj;
	kappafile >> nuj >> kappaux >> auxx >> auxx >> auxx;
	kappa[200+j] = 1 - min(1.,max(0.,kappaux)); // necessary
	nu[200+j] = nuj;
	if (altkappa<0 && nuj > 3/4. && nuj < 3/3.) kappa[200+j]=1;
	if (altkappa>0 && nuj > 3/18. && nuj < 3/14.) kappa[200+j]=1;
}
for (int j=0; j<200; j++) { nu[j] = 0.01+(nu[200]-0.01)*j/200.; kappa[j] = kappa[200]; }
Nnu += 200;
real kappa0=0;
for (int j=0; j<Nnu; j++) kappa0 += kappa[j]/Nnu;  // mean kappa is used in Newton steps
// absorption is kappa[nu]*rhof(x,y,z), but integrals are sampled by kappaApprox[K]*rhof(x,y,z)
real[int] kappaApprox(Kmax); // kappa(nu,x) = kappaApprox(K(nu))*rhof(x)  cf Lebesgue integral
for (int K=0; K<Kmax; K++) kappaApprox[K] = 1-real(K)/Kmax; // grey case: kappaApprox[0]=1

int[int] kappaInK(Kmax); // detect if kappa in (kappaApprox[K],kappaApprox[K-1])
kappaInK[0] = (Kmax==1); // used for the grey case
for (int K=1; K<Kmax; K++) {
	kappaInK[K]=0;
	for (int k=0; k<Nnu-dc-1; k+=dc) // every dc points
	if ( (kappa[k]>kappaApprox[K]) && (kappa[k]<=kappaApprox[K-1]) )
		kappaInK[K] = 1;
	if (mpirank == 0)
		cout << K << "K kappaApprox " << kappaApprox[K] << " " << kappaInK[K] << " " << kappaApprox[K-1] <<endl;
}
~~~

The last block above is a piecewise constant interpolation of $\nu\mapsto\kappa_\nu$ in `K` levels. `kappaApprox[K]` is the value of the `K`$^{th}$ level.
- if `kappaInK[K]` $=0$ then there exists a `k` such that `kappa[k]` $\in [$
`kappaApprox[K-1]`, `kappaApprox[K]`$]$.
- When needed `kappa[k]`  can be approximated by `kappaApprox[K-1]` where `K` is such that `kappaApprox[K-1]` $\leq$ `kappa[k]` $<$ `kappaApprox[K]`.

Now the x-dependency of the absorption is defined. It is the product of the density `rhof` and the change of density in the cloud.

The **scattering** coefficient is  defined below in `scatf`.

~~~freefem
/********************* DEFINE the x-dependence of kappa, called rho *********************/

/************************ Generation of cloud density *********************************/
fespace Vh2(Th2,P1);
Vh2 xi;
if(mpirank==0) {
	real sigma=0.85, mu=2.7;
	real musig2=mu+sigma*sigma;
	real aux = 1-(mu*mu+2*sigma*sigma*log(sigma))/musig2/musig2;
	gslrng ffng;
	gslrngset(ffng,random());

	for(int i=0;i<Th2.nv;i++)
		xi[][i]=0.125*musig2*(1- sqrt(aux+pow(gslrangaussian(ffng,1)*sigma/musig2,2))); //average(xi)=0.5
}
broadcast(processor(0),xi[]);

/*****************************************/
func rhof = 0.5*(1-z/2) * (1 + (Kmax > 1)) * (1 + cloud * xi*(z>cloudm)*(z<cloudM)*(sqr(x-1.5)<0.5));//(sqr(x-1.5)+sqr(y+1.5)<0.5)) );
xi = (1 + cloud * (1+xi*(sqr(x-1.5)<0.5))); // only to plot
Vh rhox = rhof; // rhox is the P1 function equal to rhof
int[int] fforderr = [1];
if (mpirank == 0) savevtk(basedir+tempestr+"_volkappa.vtk",Th3,rhox,order=fforderr);
if(mpirank==0) {plot(xi,value=1); plot(rhox,value=1,wait=0);} //medit("Cloud",Th3,rhox);}

func scatf = ascat*(z>cloudm)*(z<cloudM)*(sqr(x-1.5)<0.5);

Vh scatx=scatf; // rhox,scatx is the P1 function equal to rhof, scatf
real surf = int3d(Th3)(1.);
if (mpirank == 0) cout << "arithmetic mean of kappa(nu)*rhof = " << kappa0*int3d(Th3)(rhox)/surf << endl << endl;
~~~

The wind velocity could be read from files as a result of a Navier-Stokes simulation or defined by hand.

~~~freefem
/********************* DEFINE WIND *********************/
Vh u1, u2, u3, dT, dTh;
/*
{
ifstream velocity1("maillageCham/velocity1.txt");
ifstream velocity2("maillageCham/velocity2.txt");
ifstream velocity3("maillageCham/velocity3.txt");
velocity1 >> u1[];
velocity2 >> u2[];
velocity3 >> u3[]; 
if (mpirank==0) plot(u2, wait=1);
}
*/
u2 = 1*z; u1 = 0; u3 = 0;
~~~

### Surface integral vectors

Some integrals on the surface are computed.  For instance, `seeface1` is a vector on the vertices `i` of `ThSof` an integral on `ThS` -computed with mass lumping- of the vertical component of the normal `Ns.z` + `gamma*(Ns.x -..)` times `Q0S` times the $P^1$- hat function associated to vertex `i`.

~~~freefem
/********************* Light intensity is x cos(angle with normal)^-. Stored in seeface *********************/
VhS Q0S = Q0*(beta+(1-beta)*(z<sal)); // Light intensity modified by snow
varf vseeface1(u,v) = int2d(ThS,qft=qf1pTlump)(Q0S*(Ns.z+gamma*(Ns.x-0.2*Ns.y))*v/sqrt(1+0.04*gamma*gamma));
real[int] seeface1 = vseeface1(0,VhS);
varf vlump(u,v) = int2d(ThS,qft=qf1pTlump)(v);
real[int] blump = vlump(0,VhS);
seeface1 ./= blump;
for [i, bi : seeface1] bi = -min(bi,0.);

varf vseeface2(u,v) = int2d(ThS,qft=qf1pTlump)(Q0S*(-1+gamma*(0-0.2*0))*v/sqrt(1+0.04*gamma*gamma));
real[int] seeface2 = vseeface2(0,VhS);
seeface2 ./= blump;
for [i, bi : seeface2] bi = -min(bi,0.);

real[int] seeface = [seeface1, seeface2];
~~~

### Interpolation matrix

The rectangular matrix `R` allows to find where on surface S is a point given in the 3D FE space.
`[I,J,K]` refers to a compact definition of a matrix storing the row, column indices and coefficients when `K(I,J)` $\neq 0$. `onGamma` contains the index of the surface vertex corresponding to the same vertex of index `J[i]` in the volume mesh.

###
~~~freefem
/********************* Interpolation operator from Th3 to ThS *********************/
Vh onGamma;
{
       onGamma[] = -1;
       matrix R = interpolate(VhS,Vh);
       int[int] I(1),J(1);
       real[int] K(1);
       [I,J,K] = R;
       for (int i=0; i<I.n; i++)
               if (K[i] > 0.5) onGamma[][J[i]] = I[i];
}

/********************* Define truncated volume mesh Th3t (without the surface) *********************/
// Th3t will be used to remove surface-surface interactions from the boundary operator,
// as they will be treated separately.

fespace Ph0(Th3,P0);
varf vGamma(u,v) = int3d(Th3)(u*v);
matrix RGamma = vGamma(Vh, Ph0);
Ph0 tGamma;
Vh onG = (onGamma == -1 ? 0 : 1);
tGamma[] = RGamma*onG[];
mesh3 Th3t = trunc(Th3, tGamma == 0);
fespace Vht(Th3t,P1);

matrix Rt = interpolate(Vht,Vh); // restriction matrix from Th3 to Th3t
~~~

### Computation of the ${\mathcal H}$-matrices

For each level `K` there are two ${\mathcal H}$-matrices, one called `HVolume[K]` and the other called `HSE`. `HVolume[K]` needs to be stored and applied at each iteration, while `HSE` is computed and used on the fly for each `K` to compute `SE[K]`, the discrete counterpart of $ S^E_\nu$.  
The dynamic plugin `RadiativeTransfer_htool.dylib` contains the functions necessary for their construction.  
As the elements contain an integral of $\kappa$ on a segment, there is an interpolation of $\kappa$ on a uniform square grid to speed up the computations (see the type `KappaGrid`).  
In the end,
$$
\text{HSE}(i,j) = \int_\Sigma w^j(x')((x^i-x').n) \exp(-\int_{(x^i,x')}{\bf 1}_{\kappa(\nu,x'')\approx \kappa_{approx}[K]}d x'')/|x^i-x'|^3 d\Sigma(x'),
$$
where $x^i$ is the vertex $i$ of the volume mesh.

Similarly,
$$
\text{HVolume[K]}(i,j) = \int_\Omega w^j(x') \kappa(nu,x') \exp(-\int_{(x^i,x')}\kappa_{approx}[K]{\bf 1}_{\kappa(nu,x'')}d x'')/|x^i-x'|^2 dx'
$$

~~~freefem
/********************* BUILD H-Matrices *********************/
verbosity = 1;
HMatrix[int] HVolume(Kmax); // The Volume integral matrices
Vh[int] SE(Kmax); // The surface integral vector operator

Vht SEt;
VhS ones = 1;

for (int K=min(Kmax-1,1); K<Kmax; K++)
if (kappaInK[K]) {
	if (mpirank == 0) cout<<"building surface matrix K= "<<K<<endl;
	KappaGrid kappag(bounds,0.01,kappaApprox[K-(K>0)]*rhof);

	// We use the truncated volume mesh Th3t (without the surface) for the boundary operator,
	// as we use an analytic formula for surface-surface interactions
	Generator GenSE(Th3t,ThS,seeface,kappag);
	HMatrix HSE = Build(GenSE,VhS,Vht,eta=100,eps=1e-2,minclustersize=10);
	if (mpirank == 0) cout << HSE.infos << endl; 
	display(HSE);
	SEt[] = HSE*ones[]; // SE = \int_∑ ((x-x').n) exp(-\int_(x,x')kappa(nu,x'))/|x-x'|^3
	SE[K][] = Rt'*SEt[]; // extend the result from truncated Vht to Vh
	// use an analytic formula for surface-surface interactions:
	for [i,bi:SE[K][]] if (onGamma[][i] > -0.5) bi = 0.25*seeface1[int(onGamma[][i]+0.5)];

	if (mpirank == 0) cout<<"Building volume matrix K= "<<K<<endl;

	Generator GenVolume(Th3,kappag);
	HVolume[K] = Build(GenVolume,Vh,Vh,eta=100,eps=1e-2,minclustersize=10);
	if (mpirank == 0) cout << HVolume[K].infos << endl; 
	display(HVolume[K]); // = \int_Omega kappa(nu,x') exp(-\int_(x,x')kappa(nu,x'))/|x-x'|^2
}
~~~

### Non constant thermal diffusion

~~~freefem
/********************* OUTER LOOP *********************/ 
Vh T = 0.01; // initial value of the temperature for iterations
Vh res, sT4=sqr(sqr(T)), J=sigma*sT4;
// This is when kappatheta>0
real kappatheta=0.02;
problem heatedp(dT, dTh) = int3d(Th3) (0.02*(dx(dT)*dx(dTh)+dy(dT)*dy(dTh)+dz(dT)*dz(dTh))
		+ (u1*dx(dT)+u2*dy(dT)+u3*dz(dT))*dTh + 15*dT*dTh)
	+ int3d(Th3) (kappatheta*(dx(T)*dx(dTh)+dy(T)*dy(dTh)+dz(T)*dz(dTh)));
~~~

### The grey case

If $\kappa$ is independent of $\nu$ Stefan's law can be used because $\int_0^\infty B_\nu(T)d\nu=\sigma T^4$ where $\sigma$ is the Stefan constant.
In this case $T$ can be found explicitly from the temperature equation.

~~~freefem
if (Kmax == 1) { /********************* GREY CASE kappaApprox[0] = 1 *********************/
	for (int niter=0; niter<Niter; niter++) {
		sT4 = sigma*T*T*T*T*(1-scatx)+scatx*J; // results are independent of scatx
		res[] = HVolume[0]*sT4[];
		J[] = SE[0][];
		J = J*sigma*sqr(sqr(Tsun));
		J[] += res[];
		T = sqrt(sqrt(abs(J/sigma)));
		if (heat) {
			heatedp;
			T = T + dT;
		}
		sT4 = sT4-J;
		if (mpirank == 0) {
			//plot(J, cmm="J at iter "+ niter, value=1, fill=1);
			plot(sT4, cmm="sT4 at iter "+ niter , value=1, fill=1);
			//plot(T, cmm="T at iter "+ niter, value=1, fill=1);
		}
	}
}
~~~

### The general case

~~~freefem
else { /********************* General Nonlinear loop *********************/
	// The method:
	// Let 4 pi J(nu,x) = int_om I dw = \int_∑(x') exp(-kappa(nu,x')|x'-x|)(angle(N,x'-x))^2Bnu(nu,Tsun)d∑ 
	// 							+ \int_O(x') exp(-kappa(nu,x')|x'-x|)/|x'-x|^2 Bnu(nu,T(x'))dx'
	// get T(x) from : \int_0^∞ kappa(nu,x) Bnu(nu,T(x)) dnu = \int_0^∞ kappa(nu,x)J(nu,x)dnu =: J by Newton
	// nu-integrals as double sum: \int_0^∞ f(kappa(nu),nu)=\sum_K=1..Kmax \sum_k f(nu_k,kappa=K/Kmax)dnu_k
	// = \sum_K=1..Kmax H(kappa)*\sum_k F(K/Kmax,nu)dnu_k   when f=H(kappa)*F(kappa,nu)

	// Note: function Bnu is defined in the C++ plugin instead in order to speed up computations. The function corresponds to:
	// func real Bnu(real nu, real T) {return nu*nu*nu/(exp(nu/T)-1);}

	real[int] BnuTsun(Nnu);
	for (int k=0; k<Nnu-1 ;k++) BnuTsun[k] = Bnu(nu[k],Tsun); // aim is to speed up below
	Vh aux,aux2,aux3, bux, pres, paux = 0;
	Vh[int] Jold(Kmax); for (int k=0; k<Kmax;k++) Jold[k]=J;

	real[int] compint(Vh.ndof);
	varf vint(u,v) = int3d(Th3)(v);
	compint = vint(0,Vh);

	// cut expensive loops into MPI slices for parallel computation
	int[int] sizes(mpisize);
	int[int] offsets(mpisize);
	offsets = 0;
	int chunksize = Vh.ndof%mpisize == 0 ? Vh.ndof/mpisize : Vh.ndof/mpisize + 1;

	for (int i = 0; i < mpisize; i++)
	if (i != mpisize - 1) {
		sizes[i] = chunksize;
		offsets[i] = i*chunksize;
	}
	else {
		sizes[i] = Vh.ndof - i*chunksize;
		offsets[i] = i*chunksize;
	}

	int ii = offsets[mpirank];
	real[int] chunkbuff(sizes[mpirank]), chunkbuff2(sizes[mpirank]), chunkbuff3(sizes[mpirank]);
~~~

## The main iteration loop

Again here two integrals are computed for each level $\kappa_{approx}[K]$ of $\kappa$.
Notice that all computations are done in parallel on chunks of global vectors like  `res[]`, `aux[]`. Thence at some points the chunks need to be gathered to the global vectors, when used as input for ${\mathcal H}$-matrix vector multiplies.

~~~freefem
	for (int niter=0; niter<Niter; niter++){ // Nonlinear loop
		bux=0; J=0;
		for (int K=1; K<Kmax; K++)
		if (kappaInK[K]) {
			real auxs = 0;
			chunkbuff2 = 0;
			chunkbuff3 = 0;
			for (int k=0; k<Nnu-dc-1; k+=dc)  // every dc points, dc can be changed
			if ( (kappa[k]>kappaApprox[K]) && (kappa[k]<=kappaApprox[K-1]) ) {    
				real dnu = nu[k+dc] - nu[k];
				for [i,bi:chunkbuff] bi = rhox[](i+ii)*Bnu(nu[k],T[][i+ii]);
				chunkbuff3 += dnu*kappaApprox[K-1]*chunkbuff ; // = int_{kappa[k]=K} kappa Bnu(nu,T(x)) dnu
				if(ascat){
					for [i,bi:chunkbuff] bi = rhox[](i+ii)*Jold[K][][i+ii];
					chunkbuff2 += dnu*kappaApprox[K-1]*chunkbuff; //  = int_{kappa[k]=K} kappa Jold_nu dnu
				}
				auxs += BnuTsun[k]*dnu; // auxs = Q0*int_{kappa(nu)=K} Bnu(nu,Tsun) dnu
			}
			if(ascat) mpiAllgatherv(chunkbuff2, aux2[], mpiCommWorld, sizes, offsets); // aux2= int_K (kappa Jold) dnu
			mpiAllgatherv(chunkbuff3, aux3[], mpiCommWorld, sizes, offsets); // aux3=int_K (kappa Bnu(nu,T(x)))dnu
			if(ascat) aux = (1-scatx)*aux3 + scatx*aux2; // aux = int_K ( (1-scat)Bnu + scat Jold )kappa dnu
				 else aux=aux3;
			res[] = HVolume[K]*aux[];
			for [i,bi:chunkbuff] bi = res[][i+ii] + SE[K][][i+ii]*rhox[][i+ii]*kappaApprox[K-1]*auxs;
			mpiAllgatherv(chunkbuff, res[], mpiCommWorld, sizes, offsets);
			bux[] += aux3[]; // bux = \sum_K  \int_{kappa(nu)=K} kappa_K Bnu(nu[k],T(x))dnu
			J[] += res[]; // J(x) =  \sum_K  int_{kappa(nu)=K} kappa_K exp(-kappa|x'-x|)/|x'-x|^2 Bnu(nu,T(x)) dnu
			Jold[K][] = res[];
		}
		bux = bux - J;
		if (mpirank == 0) {
			//plot(J, cmm="J at n="+ niter ,value=1,fill=1);
			plot(bux, cmm="Bnu-J at n = "+ niter ,value=1, fill=1); // J = bux => converged
		}
~~~

### Quasi-Newton iterations

The problem is of type $F(T)=0$ for each $x$.  The derivative $F'$ is computed as if it was the grey case, i.e. $(\sigma T^4)'=4\sigma T^3$. Hence
$$
T^{m+1} = T ^m - \frac12\frac{bux-J}{4\sigma\rho\kappa_0 {T^m}^3}
$$
where $F(T^m)=bux-J$.  The $\frac12$ is to slow down the iterations.

~~~freefem
		int newton = 0; // Newton iterations for T

		for [i,bi:chunkbuff] bi = sqrt(sqrt(abs(J[][i+ii]/sigma/rhox[][i+ii]/kappa0))); // guess to start Newton
		mpiAllgatherv(chunkbuff, T[], mpiCommWorld, sizes, offsets);

		pres = abs(bux/J-1);
		while (pres[]'*compint > 0.0001*surf && newton < Newton) {
			if (newton>0) { // compute bux with new T
				chunkbuff2 = 0; 
				for (int k=0; k<Nnu-dc-1; k+=dc) { 
					real dnu = nu[k+dc]-nu[k];
					for [i,bi:chunkbuff] bi = rhox[](i+ii)*kappa[k]*Bnu(nu[k],T[][i+ii]);
					chunkbuff2 += dnu*chunkbuff;
				}
				mpiAllgatherv(chunkbuff2, bux[], mpiCommWorld, sizes, offsets);
			}
			pres = bux/J-1;
			if (mpirank == 0) {
				if (Kmax == 1) {
					paux = rhox*kappaApprox[0]*sigma*T^4; 
					cout << niter << " <- niter, error ->" << int3d(Th3)(bux/paux-1)/surf << endl;
				} 
				else 
					cout << niter << " <- niter, newton ->" << newton << " error = " << pres[]'*compint <<endl;
			}
			pres = abs(bux/J-1);
			// bux=paux; T=sqrt(sqrt(J/(sigma*rhox*kappaApprox[0]))); // error check in the grey case
			// T = T - 0.5*(bux-J)/(4*sigma*rhox*kappa0*T*T*T);
			// Newton relaxation parameter 0.5
			for [i,bi:chunkbuff] bi = T[][i+ii] - 0.5*(bux[][i+ii]-J[][i+ii])/(4*sigma*rhox[][i+ii]*kappa0*T[][i+ii]*T[][i+ii]*T[][i+ii]);
			mpiAllgatherv(chunkbuff, T[], mpiCommWorld, sizes, offsets);
			if (newton++ > Newton) bux=J;
		}

		if (heat) {
			heatedp; T=T + dT;
		}

		plot(T, cmm="T at iter = "+ niter, value=1, fill=1);
	}
	// End of non linear loop
~~~

## Graphics and other outputs

### Display of $\nu\mapsto J_\nu$

When $K=9$ (which is the maximum we did) $\nu\mapsto J_\nu$ is printed at point $(x,y,z)=(0.3,1.5,-1.5)$.

~~~freefem
	/********************* Display J_nu(x,y,z) *********************/
	if (Kmax == 9) {
		string lightstr(basedir+tempestr+"light.txt");
		ofstream myfile2(lightstr);
		for (int k=0; k<Nnu-dc-1; k+=dc) { // every dc points, can be changed
			int K = -1; 
			for (int K1=1; K1<Kmax; K1++)
			if ( (kappa[k]>kappaApprox[K1]) && (kappa[k]<=kappaApprox[K1-1]) ) K = K1;
			if (K > 0) {
				for [i,bi:chunkbuff] bi = Bnu(nu[k],T[][i+ii]);
				mpiAllgatherv(chunkbuff, res[], mpiCommWorld, sizes, offsets);
				bux[] = HVolume[K]*res[];
				for [i,bi:chunkbuff] bi = bux[][i+ii] + BnuTsun[k]*SE[K][][i+ii];
				mpiAllgatherv(chunkbuff, J[], mpiCommWorld, sizes, offsets);
				if( k>380 && k<389 ) Jold[k-380]=1e5*J; // for graphics
				if (mpirank == 0) myfile2 << 3/nu[k] << " " << 1e5*J(0.3,1.5,-1.5) << " " << kappa[k]
					<< " " << int(kappa[k]*10+0.5)/10. << " " << kappaApprox[K-1] << " " << K-1 << " " << k << endl;
			}
		}
	}
}
int tmem = storagetotal(); // max memory used at this stage
~~~

### Display of the temperature in Celcius degrees

- 1D plot of $z\mapsto T[x_0,y_0,z]$ at $x_0=1.5,y_0=-1.5$ which is approximately at the city of Chamonix.

~~~freefem
/********************* GRAPHICS *********************/
real hx=bounds[1]-bounds[0], hy=bounds[3]-bounds[2], hz=bounds[5]-bounds[4]; // domain size
mesh Tp = square(10*n, 5*n, [bounds[2]+hx*x, bounds[4]+hz*y]);
fespace Vp(Tp,P1);
Vp Temp = 4798*T(1.5,x,y)-273;
plot(Temp, value=1, fill=1);

if (mpirank == 0) {
	ofstream myfile(basedir+tempestr+"tempe.txt");
	real[int] XX(2000), YY(2000);
	int kk = 0;
	for (real Z=bounds[4]; Z<bounds[5];Z+=0.01) {
		real Tx = 4798*T(1.5,-1.5,Z)-273;
		XX[kk] = Z; YY[kk++] = Tx;
		cout << Z << "  " << Tx << endl;
		myfile << Z << "  " << Tx << endl;
	}
	XX.resize(kk); YY.resize(kk);
//	plot([XX,YY], cmm = "T in ["+YY.min+","+YY.max+"]");
}
~~~

  Two files in the vtk format are written for $T$:
  - One for $T$ in Celcius on the volume mesh
  -  The other for $T$ in Celcius on the surface mesh

~~~freefem
VhS Ts = 4798*T-273, dTs=4798*dT;
//	plot(u1); plot(u2); plot(u3);
meshS ThS2 = movemesh(ThS,[x,y,2*z]); // to enhance graphics
fespace VhS2(ThS2,P1);
VhS2 Ts2;
Ts2[]= dTs[];
Ts2[] = Ts[];
plot(Ts2,fill=1,value=1);
//if (mpirank == 0) medit("Temperature",ThS2,Ts2);
	T= 4798*T-273;
//if (mpirank == 0) medit("Temperature",Th3,T);
int[int] fforder = [1];
if (mpirank == 0) savevtk(basedir+tempestr+"_vol.vtk",Th3,T,order=fforder);
if (mpirank == 0) savevtk(basedir+tempestr+".vtk",ThS2,Ts2,order=fforder);


/*
if (mpirank == 0){
	ofstream f("solS.dat");
	f << Ts2[];
	ofstream g("solV.dat");
	g << T[];
}
*/

/*
if (mpirank == 0){
	ifstream f("solS.dat");
	VhS2 Ts2c;
	f >> Ts2c[];
	Vh Tc;
	ifstream g("solV.dat");
	g >> Tc[];

	Tc[] -= T[];
	Ts2c[] -= Ts2[];
	plot(Ts2c,fill=1,value=1,cmm="diff");
	cout << "relative L2 error surf = " << Ts2c[].l2 / Ts2[].l2 << ", vol = " << Tc[].l2 / T[].l2 << endl;
	savevtk(basedir+tempestr+"_vol_diff.vtk",Th3,Tc,order=fforder);
	savevtk(basedir+tempestr+"_diff.vtk",ThS2,Ts2c,order=fforder);
 }
*/

if (mpirank == 0) cout << "file name: " << basedir+tempestr << endl;
if (mpirank == 0) cout << "Memory used (Go) = " << (((tmem-mem)/1000.)/1000.)/1000. << endl;
~~~

|The temperature in a cutting plane in the middle of the domain |
|-----------------------|
|![][_coupe3D]          |

[_coupe3D]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/mpi/chamonix/coupe3D.png