---
name: fallingspheres
category: mesh
folder: 3d
---

## Test of mmg3d for moving objects in a mesh
The mesh around 2 spheres in a box is built, starting with the surface mesh of the 3 shapes
~~~freefem
load "tetgen" 
load "medit" 
load "mmg3d-v4.0"
include "MeshSurface.idp"

// build mesh of a box (311)  wit 2 holes  (300,310)

real hs = 0.8; 
int[int]  N=[4/hs,8/hs,11.5/hs];
real [int,int]  B=[[-2,2],[-2,6],[-10,1.5]];
int [int,int]  L=[[311,311],[311,311],[311,311]];
meshS ThH = SurfaceHex(N,B,L,1);
meshS ThSg = Sphere(1,hs,300,-1); // "gluing" surface meshs to tolat boundary meshes
meshS ThSd = Sphere(1,hs,310,-1); 

plot (ThH,ThSg,ThSd);
ThSd=movemesh(ThSd,[x,4+y,z]);

meshS ThHS=ThH+ThSg+ThSd;
medit("ThHS", ThHS); // strike the esc key to continue

~~~
The surface mesh ThHS is input for the TetGen  function $\texttt{tetg}$.

~~~freefem
real voltet=(hs^3)/6.;
cout << " voltet = " << voltet << endl;
real[int] domaine = [0,0,-4,1,voltet];
real [int] holes=[0,0,0,0,4,0];
mesh3 Th = tetg(ThHS,switch="pqaAAYYQ",regionlist=domaine,holelist=holes);    
medit("Box-With-two-Ball",Th);
~~~
You may check that Th is a volumic mesh by using the option key f1 and f2 within medit window.

Now the 2 spheres will be moved vertically  with velocity vit. The rest of the mesh is moved according to uh which is solution of the Laplacian with uh=vit on the sheres and 0 on the box.
~~~freefem
verbosity=0;
real[int] vit=[0,0,-0.3];
func zero = 0.;
func dep  = vit[2];

fespace Vh(Th,P1); 
macro Grad(u) [dx(u),dy(u),dz(u)] //

Vh uh,vh; //  to compute the displacemnt field 
problem Lap(uh,vh,solver=CG) = int3d(Th)(Grad(uh)'*Grad(vh))  //') for emacs
				  + on(310,300,uh=dep) +on(311,uh=0.); 

for(int it=0; it<29; it++){ 
  cout<<"  ITERATION       "<<it<<endl;
  Lap;
  plot(Th,uh, wait=1); // hit return to compute the next time step
  Th=mmg3dv4(Th,opt="-O 9",displacement=[zero,zero,uh]); 
 }
~~~

| The solution           |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/fallingspheres/solution.png
