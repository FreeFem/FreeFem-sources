---
name: Periodic-3d
category: math
layout: 3d
---

## Construct a periodic 3d mesh of a cube minus a cylinder

We begin with a 2d mesh of a square with a hole
~~~freefem
load "medit"

searchMethod=1;
verbosity=1;
real a=1, d=0.5, h=0.5;
border b1(t=0.5,-0.5) {x=a*t; y=-a/2; label=1;};
border b2(t=0.5,-0.5) {x=a/2; y=a*t; label=2;};
border b3(t=0.5,-0.5) {x=a*t; y=a/2; label=3;};
border b4(t=0.5,-0.5) {x=-a/2; y=a*t; label=4;};
border i1(t=0,2*pi) {x=d/2*cos(t); y=-d/2*sin(t); label=7;};
int nnb=10, nni=20; 
mesh Th=buildmesh(b1(-nnb)+b3(nnb)+b2(-nnb)+b4(nnb)+i1(nni));//, fixedborder=true);
//Th=adaptmesh(Th,0.1,IsMetric=1,periodic=[[1,x],[3,x],[2,y],[4,y]]);

~~~
Now we renumber the vertices
~~~freefem
int nz=10;
{ // Between braces to declare all var local, hence save on memory.
int[int] old2new(0:Th.nv-1);
fespace Vh2(Th,P1);
Vh2 sorder=x+y; 
sort(sorder[],old2new);
int[int]  new2old=old2new^-1;   // inverse the permuation 
//for(int i=0;i< Th.nv;++i) // so by hand. 
//  new2old[old2new[i]]=i;
Th= change(Th,renumv=new2old);
sorder[]=0:Th.nv-1;
}
~~~
Now the 3d mesh is built using $\texttt{buildlayers}$.
~~~freefem
{
  fespace Vh2(Th,P1);
  Vh2 nu;
  nu[]=0:Th.nv-1;
  plot(nu,cmm="nu=",wait=1);
}
int[int] rup=[0,5], rlow=[0,6], rmid=[1,1,2,2,3,3,4,4,7,7], rtet=[0,41];
func zmin=0;
func zmax=h;
mesh3 Th3=buildlayers(Th, nz, zbound=[zmin,zmax],
    reftet=rtet,reffacemid=rmid, reffaceup=rup, reffacelow=rlow);
for(int i=1;i<=6;++i)
  cout << " int " << i << " :  " << int2d(Th3,i)(1.) << " " << int2d(Th3,i)(1./area) << endl;
//savemesh(Th3,"Th3.mesh");
plot(Th3,wait=1);
//medit("Th3",Th3);
~~~

|The 3d mesh |
|------------|
|![][_mesh]  |

The finite element space with periodic mesh can be define and a PDE can be solved
~~~freefem
fespace Vh(Th3,P2, periodic=[[1,x,z],[3,x,z],[2,y,z],[4,y,z],[5,x,y],[6,x,y]]);
macro grad3(u) [dx(u),dy(u),dz(u)]//
func u3e =(sin(2*pi*x+1)*sin(2*pi*y+2)*sin(4*pi*z+3));
Vh u3,v3,u3h=u3e;
real ccc= pi^2*(4+4+16);
func f3= ccc*u3e;
solve P3(u3,v3) = int3d(Th3)(grad3(u3)'*grad3(v3)//'
                  +1e-6*u3*v3) - int3d(Th3)(f3*v3) 
                  + on(7,u3=u3e);
//plot(u3,u3h,cmm=1);
plot(u3,cmm=1);
u3[]-=u3h[];
cout << " err3  = " << u3[].linfty << endl;
assert(u3[].linfty < 0.07);
~~~

|The solution            |
|------------------------|
|![][_solution]          |

[_mesh]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/periodic-3d/mesh.png

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/periodic-3d/solution.png
