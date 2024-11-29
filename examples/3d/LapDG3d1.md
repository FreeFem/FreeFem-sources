---
name: LapDG3d1
category: Applied Math
folder: 3d
---

## Solve the Laplace Equations in a Cube with Discontinous Galerkin Method of degree 1

$$
-\Delta u = f,\texttt{ in } \Omega\quad u|_{\partial\Omega}=g
$$
where $\Omega$ is the unit cube, $f=1$, $g=1$.
The finite element space chosen if the discontinuous $P^1$.
~~~freefem
macro grad(u) [dx(u),dy(u),dz(u)] //
macro dn(u) (N'*grad(u) ) //'  def the normal derivative

int nn=10; 
mesh3 Th= cube(nn,nn,nn); // unit square

int[int] labs=labels(Th);
fespace Vh(Th,P1dc);     // Discontinous P1 finite element
real pena=1e5; // a paramater to add penalisation 

func f=1;
func g=1;
Vh u,v;
solve A(u,v,solver=sparsesolver) = 
  int3d(Th)( grad(u)'*grad(v))//'
+ intallfaces(Th)( pena*jump(u)*jump(v) ) 
- int2d(Th,labs)(pena*g*v)
+ int2d(Th,labs)(pena*u*v)
- int3d(Th)( f*v)
;

plot(u,cmm="Discontinue Galerkin",wait=1,value=1,fill=1);
~~~

The solution. This display is obtained by typing F (to cancel "fill") and typing L.

| The solution           |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/LapDG3d1/solution.png
