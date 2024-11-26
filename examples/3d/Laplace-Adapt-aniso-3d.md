---
name: Laplace-Adapt-3d
category: Applied Math
folder: 3d
---

## Solve the Laplace Equations in a Cube with Discontinous Galerkin Method of degree 1 and anisotropic mesh refinement

This is a continuation of Laplace-Adapt-3d.md. The problem is the same
$$
-\Delta u = f,\texttt{ in } \Omega\quad u|_{\partial\Omega}=g
$$
where $\Omega$ is the unit cube minus a half unit cube, $f=1$, $g=1$.
The finite element space chosen if the discontinuous $P^1$.

The unit cube is constructed in a different way with $\texttt{buildlayers}$.
~~~freefem
{
load "tetgen"
load "mshmet"
load "medit"
load "mmg"

int nn  = 6;

int[int] lc=[1,2,2,1,1,2]; //  label numbering 

mesh3 Th3=cube(nn,nn,nn,label=lc);
Th3 = trunc(Th3,(x<0.5) | (y < 0.5) | (z < 0.5) ,label=1);

fespace Vh(Th3,P1);
fespace Mh(Th3,[P1,P1,P1,P1,P1,P1]);
Vh u,v,usol,h3;
Mh [m11,m21,m22,m31,m32,m33];
macro Grad(u) [dx(u),dy(u),dz(u)] // EOM

problem Poisson(u,v,solver=CG) = int3d(Th3)( Grad(u)'*Grad(v) )  // ') for emacs 
  -int3d(Th3)( 1*v ) + on(1,u=0);

real lerr=0.05;
verbosity=4;

for(int ii=0; ii<4; ii++) //  BUG trap  in interation 3 
{
  Poisson;
  plot(u,wait=1);
  h3=0;
  [m11,m21,m22,m31,m32,m33]=[0,0,0,0,0,0];
  cout <<" u min, max = " <<  u[].min << " "<< u[].max << endl;
  real cc=(u[].max-u[].min);// rescale coefficiant 
 
  real[int] met=mshmet(Th3,u,hmin=1e-3,hmax=0.2,err=lerr,aniso=1);
  m11[]=met;
  Th3=mmg3d(Th3,metric=m11[],hgrad=2.3);//("oo/Th3.o.mesh");
  
  lerr *= 0.6;// change the level of error
  cout << " Th3" << Th3.nv < " " << Th3.nt << endl;
   u=u;
  if(ii>2) medit("U-adap-iso-"+ii,Th3,u,wait=1);
}
cout <<"end Laplace  Adapt aniso 3d. edp " <<endl;
}
~~~

| The third refined mesh |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/3d/Laplace-Adapt-aniso-3d/solution.png