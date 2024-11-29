---
name: sphere6
category: mesh
layout: 3d
---

##  Build a 3d mesh for a sphere by mapping the surface mesh of a cube

The mesh of the cube is built by gluing the 6 meshes of the faces and then use the mapping
$$
(x,y,z)\mapsto\frac{(x,y,z)}{\sqrt{x^2+y^2+z^2}}.
$$
~~~freefem
load "medit" 
real sqrt2=sqrt(2.);
real onesqrt2=sqrt2/2.;

mesh TS= square(10,10);
TS=adaptmesh(TS,sqrt(1+x*x+y*y),err=0.003,periodic=[[4,y],[1,x],[2,y],[3,x]]);
TS=adaptmesh(TS,sqrt(1+x*x+y*y),err=0.003,periodic=[[4,y],[1,x],[2,y],[3,x]]);
TS=TS+movemesh(TS,[-x,y])+movemesh(TS,[x,-y])+movemesh(TS,[-x,-y]);//  build symetrique mesh
plot(TS,wait=1); 
int orientation=1;
func f = 1;
int[int]  ref=[0,1]; 
meshS Thx0 = movemesh23(TS,transfo=[-f,x,y],orientation=-orientation,label=ref);
meshS Thx1 = movemesh23(TS,transfo=[+f,x,y],orientation=+orientation,label=ref);
meshS Thy0 = movemesh23(TS,transfo=[x,-f,y],orientation=+orientation,label=ref);
meshS Thy1 = movemesh23(TS,transfo=[x,+f,y],orientation=-orientation,label=ref);
meshS Thz0 = movemesh23(TS,transfo=[x,y,-f],orientation=-orientation,label=ref);
meshS Thz1 = movemesh23(TS,transfo=[x,y,+f],orientation=+orientation,label=ref);
meshS Tcube= Thx0+Thx1+Thy0+Thy1+Thz0+Thz1;
//savemesh(Tcube,"T.mesh");
//exec("ffmedit T.mesh");
//medit("Tcube",Tcube);
plot(Tcube,wait=1);
func R = sqrt(x*x+y*y+z*z); 
meshS Th = movemesh(Tcube,[x/R,y/R,z/R]);
plot(Th,wait=1);
//savemesh(Th,"T.mesh");
//exec("ffmedit T.mesh");
//medit("Th",Th);

plot(Th);
~~~

|The 3d mesh             |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/sphere6/SOLUTION.png
