---
name: Stokes_composite
category: fluid
layout: example
---

# The Stokes system for creeping flow solved with composite FE
The system is set in 2 dimensions:
$$
-\Delta u + \nabla p =\vec f, \quad
 \nabla\cdot u = 0   \hbox{ in }\Omega ,~~
 u _{|\Gamma}= \vec g
 $$
with $\Omega$ the  square $(0,2\pi)\times(0,2\pi)$, $\vec f=[0,-4\cos x\sin y]^T$ and $\vec g=[\sin x\cos y,-\cos x \sin y]^T$.

A possible variational formulation is
$$∀(v,q),~~∫_Ω ∇u:∇v−∫_Ωp\nabla\cdot v=∫_Ωf\cdot v,
−∫_Ω\nabla\cdot u q−∫_Ωϵpq=0.
$$
~~~freefem
macro grad(u) [dx(u),dy(u)]//
macro Grad(u1,u2) [ grad(u1), grad(u2)]//
macro div(u1,u2) (dx(u1)+dy(u2))//

// definition of the boundary condition 
func g1 = sin(x)*cos(y);
func g2 = -cos(x)*sin(y);

// defintion of the right hand side
func f1 = 0;
func f2 = -4*cos(x)*sin(y);

int nn = 30; // number of edge in each direction
mesh Th=square(nn,nn,[2*pi*x,2*pi*y],flags=3);
mesh ThP=Th;                    // Pressure mesh
mesh ThU=trunc(ThP,1,split=2);  // Velocity mesh

fespace Uh(ThU,[P1,P1]);
fespace Ph(ThP,P1);
~~~
While it is unecessary to use flags=3, it does the following:
it will produce a mesh where all quads are split with diagonal $x-y=$ constant and forbids of a 3 vertices on the boundary at corners.

Using different meshes for the velocity and pressure is also not necessary but it is given to illustrate composite elements as explained in the documentation here

https://doc.freefem.org/documentation/composite.html

The solver is standard.  The penalization of $pq$ is needed to avoid the singularity of the system because in Stokes' equations, pressure is defined up to a constant.
~~~freefem
Uh [u1,u2],[v1,v2];
Ph p,q;

solve Stokes ( <[u1,u2],[p]>, <[v1,v2],[q]>) = int2d(ThU)( (Grad(u1,u2):Grad(v1,v2)) )
+ int2d(ThU)( - div(u1,u2)*q - div(v1,v2)*p )
+ int2d(ThP)( -1e-10*p*q )
- int2d(ThU) ( [f1,f2]'*[v1,v2] )
+ on(1,2,3,4,u1=g1,u2=g2);

plot( u1, cmm="u1" );
plot( u2, cmm="u2" );
plot( p,  cmm="p" );
~~~
## Results

| Isovalue lines of the pressure |
| ------------------------------ |
| ![][_solution]                 |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/stokes_composite/solution.png
