---
name: LapDG3d
category: Applied Math
folder: 3d
---

## Solve the Laplace Equations in a Cube with Discontinous Galerkin Method of degree 2

based on paper by Riviere, Beatrice; Wheeler, Mary F.; Girault, Vivette
called "A priori error estimates for finite element methods based on discontinuous approximation spaces for elliptic problems. SIAM J. Numer. Anal. 39 (2001), no. 3, 902--931 (electronic).

solve $ -\Delta u = f$ in $\Omega$ with $u= g$ on $\Gamma$.
The domain is the unit cube
~~~freefem
int nn=5;
mesh3 Th = cube(nn,nn,nn);
~~~
Two FE spaces are tested: $P^2$ and $P^2$-discontinuous element.  Notice that the formulation is non-symmetric.
~~~freefem
macro grad(u) [dx(u),dy(u),dz(u)] //
macro dn(u) (N'*grad(u) ) //'  def the normal derivative
int[int] labs=labels(Th);
fespace Vh(Th,P2dc3d);     // Discontinous P2 finite element
fespace Xh(Th,P2);
real pena=1e2; //to add penalisation for Vh, not needed if Xh is used
func f=1;
func g=0;
Vh u,v;
Xh uu,vv;
problem A(u,v,solver=sparsesolver) = 
   int3d(Th)( grad(u)'*grad(v)) //'
 + intallfaces(Th)(//  loop on all  edges of all triangles
           ( jump(v)*mean(dn(u)) -  jump(u)*mean(dn(v)) 
          + pena*jump(u)*jump(v) ) / nElementonB 
)
- int3d(Th)(f*v) 
- int2d(Th)(g*dn(v)  + pena*g*v) 
;
problem A1(uu,vv,solver=sparsesolver) 
= 
 int3d(Th)(grad(uu)'*grad(vv))//'
  - int3d(Th)(f*vv) + on(labs,uu=g);
 
 A; // solve  DG
 A1; // solve continuous

 real err= sqrt(int3d(Th) ( sqr(u-uu)));
 cout << " err= "<< err<< " " << u[].linfty << " " << uu[].linfty << endl; 
plot(u,uu,cmm="Discontinue Galerkin",wait=1,value=1);
plot(u,cmm="Discontinue Galerkin",wait=1,value=1,fill=1);
assert(err< 0.01);
~~~

| The solution           |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/LapDG3d/solution.png
