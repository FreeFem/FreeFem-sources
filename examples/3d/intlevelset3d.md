---
name: intlevelset3d
category: mesh
folder: 3d
---

## Manipulation of data and PDE on surfaces defined by constant values of a function.

A cube with extreme points [-2,-2,-2] and [2,2,2] is triangulated with $16^3$ points. All labels are set to 1.
~~~freefem
load "medit"
include "cube.idp"
real surfS1 = 4*pi;
real volS1 =surfS1/3.; 
int nn= 16; 
int[int]  Nxyz=[nn,nn,nn];
real [int,int]  Bxyz=[[-2.,2.],[-2.,2.],[-2.,2.]];
int [int,int]  Lxyz=[[1,1],[1,1],[1,1]];
mesh3 Th=Cube(Nxyz,Bxyz,Lxyz);
~~~
A  surface of equal values of the function $x^2+y^2+z^2-1$  is measured.
~~~freefem
int err=0;
real eps = 0.5;
func r = sqrt(x*x +y*y+z*z);

real lc ;
verbosity=3;
lc = int2d(Th,levelset=r-1.)(1.) ; 
cout << " area of the level set = " <<  lc  << " =  surfS1 " << surfS1 ;
cout << ", Ok = " << (abs(lc-surfS1) < eps) << endl; 
if( abs(lc-surfS1) > eps) err++;
~~~
A bilinear form is defined on the surface and a vector $vv$  of values on the hat functions.
~~~freefem
fespace Vh(Th,P1);
// test linear and bilinear ... 
varf vl(u,v) = int2d(Th,levelset=r-1.)(v) + int2d(Th,levelset=r-1.)(u*v);
real[int] vv=vl(0,Vh);

cout << " area of the level set (varf linear ) = " <<  (lc=vv.sum)  << "=  surfS1 " << surfS1 ;
cout  << ", Ok = " << (abs(lc-surfS1) < eps) << endl;
if( abs(lc-surfS1) > eps) err++; 
~~~
Now let us prepare the data to solve a PDE on the level sets.

Just for test a idea approximation of int of negative part of levelset  to we just change the measure of the element not the quadrature point .
~~~freefem
real[int]  one(Vh.ndof); 
one=1.;
matrix VV=vl(Vh,Vh); //  matrix with levelset
vv = VV*one;
cout << " area of the level set (varf bilinear same) = " <<  (lc=vv.sum)  << "=  surfS1 " << surfS1;
cout << ", Ok = " << (abs(lc-surfS1) < eps) << endl;
if( abs(lc-surfS1) > eps) err++;
~~~
Test the new stuff for level set ...
~~~freefem
{
    macro grad(u) [dx(u),dy(u),dz(u)] //
    Vh u,v;
    solve Pxx(u,v) = int3d(Th) ( grad(u)'*grad(v)*1e-8 ) + int3d(Th, levelset= 1-r) ( grad(u)'*grad(v) ) + on(1,u=0) + int3d(Th, levelset= 1-r) ( 1*v);
    plot(u,wait=1);   
    varf vxx(u,v) =  int3d(Th, levelset= 1-r) ( u*v ) + int3d(Th, levelset= 1-r) ( 1*v);
  matrix XX=vxx(Vh,Vh);
  real[int] xx=vxx(0,Vh);
  real vol1= int3d(Th, levelset= 1-r)(1.);
  cout << "   vol1 = " << vol1 << "  ~= " << Th.measure - volS1 << endl;
  err += (abs(vol1-(Th.measure - volS1)) > eps); 
  cout << " xx.sum = " << xx.sum << " == " << vol1 <<endl;
  err += (abs(vol1-xx.sum) > 1e-8); 
  
  real[int] yy(Vh.ndof); yy=1;
  xx= XX*yy;
  cout << " XX.sum = " << xx.sum << " == " << vol1 << endl;
  err += (abs(vol1-xx.sum) > 1e-8); 
}
~~~
test on diff mesh3  not yet implemented
~~~freefem
if(0)
{
mesh3 Th1=Cube(Nxyz,Bxyz,Lxyz);
mesh3 Th2=Cube(Nxyz,Bxyz,Lxyz);
fespace Vh1(Th1,P1);
fespace Vh2(Th2,P1);

varf vl(u,v) = int2d(Th,levelset=r-1.)(v) + int2d(Th,levelset=r-1.)(u*v);
real[int] vv=vl(0,Vh2);

cout << " area of the level set (varf linear diff    ) = " <<  (lc=vv.sum)  << "=  surfS1 " << surfS1 ;
cout  << ", Ok = " << (abs(lc-surfS1) < eps) << endl;
if( abs(lc-surfS1) > eps) err++; 
real[int]  one(Vh1.ndof); 
one=1.;

matrix VV=vl(Vh1,Vh2);
vv = VV*one;
cout << " area of the level set (varf bilinear diff ) = " <<  (lc=vv.sum)  << "=  surfS1 " << surfS1;
cout << ", Ok = " << (abs(lc-surfS1) < eps) << endl;; 
if( abs(lc-surfS1) > eps) err++;

}
cout << " Nb err " << err << endl;
assert(err==0);
~~~

| The surface mesh extracted |
|----------------------------|
|![][_solution]              |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/intlevelset3d/solution.png