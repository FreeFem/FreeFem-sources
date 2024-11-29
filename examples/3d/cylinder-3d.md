---
name: cylinder-3d
category: mesh
folder: 3d
---
## Build volumic or surface meshes of a 3d cylinder

This is a very simple example to show how to build a 3d mesh of a cylinder by extruding the 2d mesh of a circle
~~~freefem
load "medit"// contains the interface with medit graphic display
load "tetgen"// used by the second part of this program

border C(t=0,2*pi) { x = cos(t); y=sin(t); label=1;}
mesh Baseh = buildmesh(C(20));
plot(Baseh,wait=0);// press return to go to the next plot

int[int] rup=[0,1],  rdown=[0,2], rmid=[1,3];
func zmin= 1;
func zmax= 10;
int nlayer=100;
mesh3 Th=buildlayers(Baseh,nlayer,
  coef= 1.,
  zbound=[zmin,zmax],
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);
// To see the volumic mesh use the f1 then f2 control functions on the keyboard
medit("Cyl1",Th,wait=1); // press escape in the medit window to go to the next plot

plot(Th,cmm="Cyl1");
~~~

| The solution is on the figure below |
|------------------------|
|![][_cylinder]          |


## Build a 2d mesh of a 3d cylinder

Here we wrap above the mesh of a 2d disk the mesh of a square and close it with the mesh of the same 2d disk but at the top of the cylinder.  We need to call $\texttt{Tetgen}$ functions.
Here the axis of the cylinder is Ox
~~~freefem
int nx=10; // nomber of points along the axis
int nth=50; // number of points on the cercle
real xmin=1.,xmax=3.; // bottom and top of the cylinder

// The mesh of the disk
mesh Thcercle = buildmesh(C(nth));

// This mesh will be wrapped around and above the bottom cercle
mesh Thcarre=square(nx,nth,[xmin+x*(xmax-xmin),2*pi*y]);
~~~
We use $\texttt{movemesh23}$ to transform the mesh of the square into a lateral circular mesh; the same function is used to translate the mesh of the disk at $\texttt{xmin}$ and $\texttt{xmax}$ perpendicular to Ox.
~~~freefem
cout<<" Build a lateral mesh\n";
func f1 = x;
func f2 = cos(y);
func f3 = sin(y);
meshS Thsurf1=movemesh23(Thcarre,transfo=[f1,f2,f3],orientation=-1);
meshS Thsurf2=movemesh23(Thcercle,transfo=[xmin,x,y],orientation=-1);
meshS Thsurf3=movemesh23(Thcercle,transfo=[xmax,x,y],orientation=1);
~~~
Finally the 3 meshes are glued
~~~freefem
meshS Thsurf=Thsurf1+Thsurf2+Thsurf3; // To see that the mesh is surfacic only use the f1 then f2 control functions. To close the window hit the escape key
medit("Cyl2",Thsurf,wait=1);
~~~
If a volumic mesh is desired from the surface mesh (to fill it inside with tetraedra, use the $\texttt{Tetgen}$ function as follows
~~~freefem
real voltet= ( ( (2*pi)/50 )^3 )/6.;
cout << "  voltet = " << voltet << endl;
real[int] domaine = [1.5,0.,0.,1,voltet];
mesh3 Th2=tetg(Thsurf,switch="pqaaAAYYQ",nbofregions=1,regionlist=domaine);

//savemesh(Th,"cyl.mesh");
medit("Cyl3",Th2,wait=1);
// FFCS: testing 3d plots
plot(Th2,cmm="Cyl3");
~~~

| This mesh is surfacic |
|-----------------------|
|![][_surface]          |

| This mesh is volumic |
|----------------------|
|![][_volume]          |

[_cylinder]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/cylinder-3d/cylinder.png

[_surface]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/cylinder-3d/surface.png

[_volume]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/cylinder-3d/volume.png
