---
name: beam-3d
category: solid
layout: module
---

# Compute the deflection of a 3d elastic beam
The geometry is a beam of rectangular cross section.
The equations are
$$
\displaystyle{
	\begin{align*}&
  \int_{\Omega}[\mu\epsilon(u)\cdot\epsilon(v) +\lambda (\nabla\cdot u )(\nabla\cdot v)]  = \int_{\partial\Omega}{g u_3},
  \quad \forall v\in V(\Omega)
  \cr&
  \nabla\cdot u= \partial_x u_1+\partial_y u_2+\partial_z u_3,
  \cr&
  \epsilon(u)=[\partial_x u_1,\partial_y u_2,\partial_z u_3, \frac1{\sqrt{2}}(\partial_z u_2 + \partial_y u_3), \frac1{\sqrt{2}}(\partial_y u_1 + \partial_x u_2)].
  \end{align*}
}
$$
where $V(\Omega)=\{ v\in (H(\Omega)^3): v|_S=0\}$ and $S\subset\partial\Omega$ is the portion of the boundary where the beam is clamped.
$\mu$ and $\lambda$ are computed from the Young and Poisson constants $E,\sigma$,
~~~freefem
real E = 21.5e4;
real sigma = 0.29;
real mu = E/(2*(1+sigma));
real lambda = E*sigma/((1+sigma)*(1-2*sigma));
real gravity = -0.05;
~~~
Gravity is vertical and the value is chosen to obtain a good picture.  The equatins hold for small displacements anyway.
The mesh is built with the function $\texttt{Cube}$ with $20\times5\times5$ edge vertices. The size of the beam is specified by [min x, max x],[min y, max y],[min z, max z] in Bxyz. The label of the beam faces are in Lxyz: [xmin face, xmax face],[ymin face, ymax face],[zmin, zmax]. Hence the  face at xmin has label 1. It will be where the beam is clamped.
~~~freefem
include "cube.idp"
int[int]  Nxyz=[20,5,5];
real [int,int]  Bxyz=[[0.,5.],[0.,1.],[0.,1.]];
int [int,int]  Lxyz=[[1,2],[2,2],[2,2]];
mesh3 Th=Cube(Nxyz,Bxyz,Lxyz);
~~~
A vector finite element space with 3 $P^1$ components is defined and 2 functions with vector values $u,v$:
~~~freefem
fespace Vh(Th,[P1,P1,P1]);
Vh [u1,u2,u3], [v1,v2,v3];
~~~
Everything is in place to solve the problem
~~~freefem
real sqrt2=sqrt(2.);
macro epsilon(u1,u2,u3)  [dx(u1),dy(u2),dz(u3),
                  (dz(u2)+dy(u3))/sqrt2,(dz(u1)+dx(u3))/sqrt2,(dy(u1)+dx(u2))/sqrt2] // EOM
macro div(u1,u2,u3) ( dx(u1)+dy(u2)+dz(u3) ) // EOM
 // EOM means End Of Macro
solve Lame([u1,u2,u3],[v1,v2,v3])=
  int3d(Th)(  
	    lambda*div(u1,u2,u3)*div(v1,v2,v3)	
	    +2.*mu*( epsilon(u1,u2,u3)'*epsilon(v1,v2,v3) ) //')
	      )
  - int3d(Th) (gravity*v3)
  + on(1,u1=0,u2=0,u3=0)
  ;
~~~
To visualize the result, the best is to move the mesh by the displacement just computed.
~~~freefem
real dmax= u1[].max;
cout << " max deplacement = " << dmax << endl;
real coef= 0.1/dmax;
int[int] ref2=[1,0,2,0];
mesh3 Thm=movemesh3(Th,transfo=[x+u1*coef,y+u2*coef,z+u3*coef],label=ref2);
Thm=change(Thm,label=ref2);
plot(Th,Thm, wait=1,cmm="coef  amplification = "+coef );
~~~
| The deflection of the beam |
|----------------------------|
|![][_deflection]            |


[_deflection]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/beam-3d/beam3d.png

