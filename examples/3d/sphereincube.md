---
name: sphereincube
category: mesh
layout: 3d
---

##  Build the mesh of a cube with a sphere inside it

The surface mesh of the cube is built using $\texttt{SurfaceHex}$ and that of the sphere with $\texttt{Sphere}$.  Then $\texttt{tetg}$ builds a volumic mesh from the sum of both surface meshes.

~~~freefem
verbosity=1;
load "tetgen"
load "medit"
include "MeshSurface.idp"

real hs = 0.1; // mesh size
real Lex=1,Lin =0.5;
int nLex = 2*Lex/hs;

int[int]  Nex=[nLex,nLex,nLex];
real [int,int]  Bex=[[-Lex,Lex],[-Lex,Lex],[-Lex,Lex]];
int [int,int]  Lbex=[[1,1],[1,1],[1,1]] ;
int  Ls=2 ;

////////////////////////////////
meshS ThHS = SurfaceHex(Nex,Bex,Lbex,1)+Sphere(Lin,hs,Ls,1); // "gluing" surface meshs to total boundary meshes
real voltet=(hs^3)/6.;
cout << " voltet = " << voltet << endl;
real[int] domaine = [0,0,0,2,voltet,0,0,(Lex+Lin)/2,1,voltet];

mesh3 Th = tetg(ThHS,switch="pqaAYY",regionlist=domaine);
Th= trunc(Th,region==1);
// Tetrahelize the interior of the cube with tetgen
// medit("tetg",Th,wait=1);
// savemesh(Th,"Th.mesh");
plot(Th);
~~~

|The 3d mesh             |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/sphereincube/solution.png