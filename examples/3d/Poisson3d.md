---
name: Poisson3d
category: math
layout: 3d
---

## Construct a periodic 3d mesh of a cube minus a cylinder

A volumic mesh of a sphere is built by mapping a surface mesh of a square and then call $\texttt{tetgentransfo}$.
~~~freefem
load "tetgen"
load "medit"

mesh Th=square(10,20,[x*pi-pi/2,2*y*pi]);  //  $]\frac{-pi}{2},frac{-pi}{2}[\times]0,2\pi[ $
//  a parametrization of a sphere 
func f1 =cos(x)*cos(y);
func f2 =cos(x)*sin(y);
func f3 = sin(x);
//  de  partiel derivative of the parametrization DF
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
Th=adaptmesh(Th,m11*vv,m21*vv,m22*vv,IsMetric=1,periodic=perio);
plot(Th,wait=1);

verbosity=2;
real[int] domaine =[0.,0.,0.,1,0.01];
mesh3 Th3=tetgtransfo(Th,transfo=[f1,f2,f3],nbofregions=1,regionlist=domaine);
//savemesh(Th3,"sphere.meshb");
//medit("sphere",Th3);
plot(Th3,cmm="sphere");
~~~
The folowing problem is solved
$$
-\Delta u=f,\quad u=u_e ~~on~\Gamma_0\cup\Gamma1
$$
and homogeneous Neumann conditon on other boundaries. A $P^2$ FEM is used
~~~freefem
fespace Vh(Th3,P23d);
func ue =   2*x*x + 3*y*y + 4*z*z+ 5*x*y+6*x*z+1;
func f= -18. ;
Vh uhe = ue; // bug ..
cout << " uhe min:  " << uhe[].min << " max:" << uhe[].max << endl;
cout << uhe(0.,0.,0.) << endl;

Vh u,v;

macro Grad3(u) [dx(u),dy(u),dz(u)]  // EOM

problem Lap3d(u,v,solver=CG)=int3d(Th3)(Grad3(v)' *Grad3(u))//'
 - int3d(Th3)(f*v) + on(0,1,u=ue);
Lap3d;
cout << " u min::   " << u[]. min << "  max: " << u[].max << endl;
real err= int3d(Th3)( square(u-ue) );
cout << int3d(Th3)(1.) << " = " << Th3.measure << endl;
Vh d= ue-u;
cout <<  " err = " << err <<  " diff l^\intfy = " << d[].linfty << endl;
plot(u);
~~~
To compare with the exact solution a 2d mesh is built and the functions are displayed in a cut corresponding to this 2d mesh.
~~~freefem
border cc(t=0,2*pi){x=cos(t);y=sin(t);label=1;}
mesh Th2=buildmesh(cc(50));
fespace Vh2(Th2,P2);

Vh2 u2=u,u2e=ue-u2;
plot(u2e,wait=1, value=1);

assert(err < 1e-9);
~~~

| The 3d solution |
|-----------------|
|![][_solution]   |

| The 2d error |
|--------------|
|![][_error]   |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Poisson3d/solution.png

[_error]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Poisson3d/error.png
