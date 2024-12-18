---
name: tetgenholeregion
category: mesh
layout: 3d
---

## Various Experiments with TetGen

First the surface mesh of a sphere is built by mapping the mesh of a square.

The first step is to obtain an adapted periodic mesh for a rectangle

~~~freefem
load "tetgen"
load "medit"

mesh Th=square(10,20,[x*pi-pi/2,2*y*pi]);
//  a parametrization of a sphere 
func f1 =cos(x)*cos(y);
func f2 =cos(x)*sin(y);
func f3 = sin(x);
//  partiel derivative of the parametrization DF
func f1x=sin(x)*cos(y);   
func f1y=-cos(x)*sin(y);
func f2x=-sin(x)*sin(y);
func f2y=cos(x)*cos(y);
func f3x=cos(x);
func f3y=0;
// $  M = DF^t DF $
func m11=f1x^2+f2x^2+f3x^2;
func m21=f1x*f1y+f2x*f2y+f3x*f3y;
func m22=f1y^2+f2y^2+f3y^2;

func perio=[[4,y],[2,y],[1,x],[3,x]];  
real hh=0.1;
real vv= 1/square(hh);
verbosity=2;
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
//medit("squaremesh",Th,wait=1);
plot(Th);
~~~

Next, $\texttt{movemesh23}$ will be used to map a 2d mesh to a 3d shape

~~~freefem
verbosity=2;
// construction of the surface of spheres
real Rmin  = 1.;
func f1min = Rmin*f1;
func f2min = Rmin*f2;
func f3min = Rmin*f3;
//savemesh(Th,"Th.mesh");
meshS Th3sph=movemesh23(Th,transfo=[f1min,f2min,f3min],orientation=-1);
//savemesh(Th3sph,"sphereR1.mesh");
//medit("sphereR1",wait=1,Th3sph);
plot(Th3sph);
~~~

Now the same operation is done to construct a bigger sphere

~~~freefem
real Rmax  = 2.;
func f1max = Rmax*f1;
func f2max = Rmax*f2;
func f3max = Rmax*f3;
meshS Th3sph2=movemesh23(Th,transfo=[f1max,f2max,f3max],orientation=1);
cout << "=====================" << endl;
//savemesh(Th3sph2,"sphereR2.mesh");
cout << "addition" << endl;
meshS Th3=Th3sph+Th3sph2;
//savemesh(Th3,"sphereAdd.mesh");
//medit("sphereSurfaceAdd",wait=1,Th3);
plot(Th3);
~~~

Now that we have a sphere inside another, tegen is called

~~~freefem
real[int] domain2 = [1.5,0.,0.,145,0.001,0.5,0.,0.,18,0.01];
cout << "==============================" << endl;
cout << " tetgen call without hole " << endl;
cout << "==============================" << endl;
mesh3 Th3fin=tetg(Th3,switch="paAAYYCCV",nbofregions=2,regionlist=domain2);
cout << "=============================" << endl;
cout << "finish: tetgen call without hole" << endl;
cout << "=============================" << endl;
//savemesh(Th3fin,"spherewithtworegion.mesh");
//medit("spherewithtworegion",wait=1,Th3fin);
plot(Th3fin);
~~~

Finally the inner sphere is defined as a hole and tetgen is called again

~~~freefem
real[int] hole = [0.,0.,0.];
real[int] domain = [1.5,0.,0.,53,0.001];
cout << "=============================" << endl;
cout << "  tetgen call with hole   " << endl;
cout << "=============================" << endl;
mesh3 Th3finhole=tetg(Th3,switch="paAAYCCV",nbofholes=1,holelist=hole,nbofregions=1,regionlist=domain);
cout << "=============================" << endl;
cout << "finish: tetgen call with hole   " << endl;
cout << "=============================" << endl;
//savemesh(Th3finhole,"spherewithahole.mesh");
//medit("spherewithahole",wait=1,Th3finhole);
plot(Th3finhole);
~~~

| The final result       |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/TetgenholeRegion/solution.png