---
name: Periodic-Poisson-cube-ballon
category: math
layout: 3d
---

## Construction of a periodic domain with a cube and a ball

First a surface mesh is made for a cube plus a ball inside; $\texttt{SurfaceHex}$ triangulate the surface of a cube and $\texttt{Sphere}$ triangulate the surface of a sphere.  Then $\texttt{tetg}$ is invoked to make a 3d triangulation of the volume delimited by these 2 surfaces.
~~~freefem
verbosity=1;

load "tetgen"
load "medit"
include "MeshSurface.idp"

mesh3 Th;
try {
  Th=readmesh3("Th-hex-sph.mesh");
 }
catch(...)
  { 
    real hs = 0.2;  // mesh size on sphere 
    int[int]  NN=[11,9,10];
    real [int,int]  BB=[[-1.1,1.1],[-.9,.9],[-1,1]];
    int [int,int]  LL=[[1,2],[3,4],[5,6]];
    
    ////////////////////////////////
    meshS ThHS = SurfaceHex(NN,BB,LL,1)+Sphere(0.5,hs,7,1); // "gluing" surface meshs to tolat boundary meshes
    real voltet=(hs^3)/6.;
    cout << " voltet = " << voltet << endl;
    real[int] domaine = [0,0,0,1,voltet,0,0,0.7,2,voltet];
    
    Th = tetg(ThHS,switch="pqaAYY",nbofregions=2,regionlist=domaine);    
    // Tetrahelize the interior of the cube with tetgen
    medit("tetg",Th,wait=1);
 //   savemesh(Th,"Th-hex-sph.mesh");
    // FFCS: testing 3d plots
    plot(Th,wait=1);
 }
~~~
The periodic $P^1$ finite element space Vh is defined
~~~freefem
fespace Ph(Th,P0);
verbosity=50;
fespace Vh(Th,P1,periodic=[[3,x,z],[4,x,z],[1,y,z],[2,y,z],[5,x,y],[6,x,y]]);// back and front
verbosity=1;
Ph reg=region;

cout << "  centre = " << reg(0,0,0) << endl;
cout << " exterieur = " << reg(0,0,0.7) << endl;
~~~
An elliptic problem is set with discontinuous coefficients of 100 in region 1 and 2 in region 2.
~~~freefem
macro Grad(u) [dx(u),dy(u),dz(u)] // EOM

Vh uh,vh;
real x0=0.3,y0=0.4,z0=06;
func f= sin(x*2*pi+x0)*sin(y*2*pi+y0)*sin(z*2*pi+z0);
real gn = 1.;
real cf= 1;
problem P(uh,vh,solver=sparsesolver)=
     int3d(Th,1)( Grad(uh)'*Grad(vh)*100) +  int3d(Th,2)( Grad(uh)'*Grad(vh)*2)
  + int3d(Th) (vh*f)
  ; 
  
  P;

plot(uh,wait=1, nbiso=6);
medit("   uh ",Th, uh,wait=1); 
~~~

|The solution and a cut into the domain (do "b" then "m" then f1 then f2 in the medit window) to see in the volume |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Period-Poisson-cube-ballon/solution.png
