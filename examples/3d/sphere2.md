---
name: sphere2
category: mesh
layout: 3d
---

##  Build a 3d mesh for a sphere by glueing the meshes of two half spheres

~~~freefem
load "medit"
border BC(t=0,2*pi){ x=cos(t);y=sin(t);label=1;}
mesh TC= buildmesh(BC(100));
func r
 = (1e-5+abs(1-square(x)-square(y)))^0.5;
real cc=100;
//TC=adaptmesh(TC,[cc*x/r,0,cc*y/r],IsMetric=1);
plot(TC,wait=1);
meshS Thup=movemesh23(TC,transfo=[x,y,sqrt(abs(1-square(x)-square(y)))]);
meshS Thdown=movemesh23(TC,transfo=[x,y,-sqrt(abs(1-square(x)-square(y)))], orientation=-1);
verbosity=10;
meshS Th= Thup+Thdown;

bool usemedit=false;
if(usemedit)
  {
    medit("Thup",Thup,wait=1);
    medit("Thdown",Thdown,wait=1);
    medit("Th",Th,wait=1);
  }
//savemesh(Th,"toto.mesh");
plot(Thup);
plot(Thdown);
plot(Th);
~~~

|The 3d mesh             |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/sphere2/solution.png