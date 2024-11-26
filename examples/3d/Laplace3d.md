---
name: Laplace3d
category: Applied Math
folder: 3d
---

## Solve the Laplace Equations in a Cube with P2 Elements

First build the mesh
~~~freefem
verbosity=2;

int nn=20;
mesh Th2=square(nn,nn,region=0);
fespace Vh2(Th2,P2);
Vh2 ux,uz,p2;
int[int] rup=[0,2],  rdown=[0,1], rmid=[1,1,2,1,3,1,4,1];
real zmin=0,zmax=1;
mesh3 Th=buildlayers(Th2,nn,
  zbound=[zmin,zmax],
  labelmid=rmid, 
  reffaceup = rup,
  reffacelow = rdown);
~~~
This method is more complex that $\texttt{ mesh3 Th = cube(nn,nn,nn);}$.  Now the following problem is solved
$$
\int_{Th}\nabla u\cdot\nabla v = \int_{Th}f v, \quad \frac{\partial u}{\partial n}|_{\Gamma_2}= \vec n\cdot\nabla u_e,\quad u|_{\Gamma_1}=u_e
$$
and $f$ is chosen  to be $-\Delta u_e$ so that the solution is $u_e$
~~~freefem
fespace Vh(Th,P2);

func ue =   2*x*x + 3*y*y + 4*z*z + 5*x*y+6*x*z+1;
func uex=   4*x+  5*y+6*z;
func uey=   6*y + 5*x;
func uez=   8*z +6*x;
func f= -18. ;
Vh uhe = ue; //

cout << " uhe min:  " << uhe[].min << " max:" << uhe[].max << endl;

Vh u,v;

macro Grad3(u) [dx(u),dy(u),dz(u)]  // EOM


  problem Lap3d(u,v,solver=CG)  =
  int3d(Th)(Grad3(v)' *Grad3(u)) //') for emacs 
  + int2d(Th,2)(u*v)  
  - int3d(Th)(f*v) 
  - int2d(Th,2) ( ue*v + (uex*N.x +uey*N.y +uez*N.z)*v )
  + on(1,u=ue);
Lap3d;
~~~
Various quantities are printed including the effect of a quadrature change.
~~~freefem
cout << " u min::   " << u[]. min << "  max: " << u[].max << endl;
real err= int3d(Th)( square(u-ue) );
real aa1= int3d(Th,qfV=qfV1)(u) ;
real aa2= int3d(Th,qfV=qfV1lump)(u) ;

cout << " aa1 = " << aa1 << endl;
cout << " aa2 = " << aa2 << endl;
cout << int3d(Th)(1.) << " = " << Th.measure << endl;
cout << int3d(Th,qfV=qfV1)(1.) << " = " << Th.measure << endl;
cout << int3d(Th,qfV=qfV1lump)(1.) << " = " << Th.measure << endl;
Vh d= ue-u;
cout <<  " err = " << err <<  " diff l^\intfy = " << d[].linfty << endl;
real  aire2=int2d(Th,2)(1.); // bug correct in version 3.0-4 
cout << " aire2 = " << aire2 << endl;
func uuu= 2.+x;
cout << uuu(1,1,1) << endl;
assert( abs(aire2-1.) < 1e-6);
plot(u,wait=1);

assert( err < 1e-6);
~~~

|The solution            |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Laplace3d/solution.png