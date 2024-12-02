---
name: 3d-Leman
category: fluid
layout: module
---

# 3d-Potential Flow
Aims at computing the flow in the Leman lake with water input from the Rhone on the east side (label 1) and water output at the western tip near Geneva (label 2). The velocity is the gradient of the potential $u$ which solves
$$
\left\{
\begin{align*}
	-\Delta u  &= 0 \hbox{ in $\Omega_3$}\\
	\frac{\partial u}{\partial n}|_{\partial\Omega_3} &= q
\end{align*}
\right.
$$
where $\Omega_3$ is the region of space occupied by water.
The solution exists only iff $\int_{\partial\Omega_3}q=0$ and uniqueness holds up to a constant.
### Variational form
Given a small positive parameter $\epsilon$, find $u\in H^1(\Omega_3)$ such that
$$
\displaystyle{
	\int_{\Omega_3}{\nabla u\cdot\nabla v +\epsilon u v}  = \int_{\partial\Omega_3}{q v},
  \quad \forall v\in H^1(\Omega_3)
}
$$
It can be shown that $u$ converge to the solution when $\epsilon\to0$ iff $\int_{\partial\Omega_3}{q}=0$.

The 3D mesh is built layer by layer from a 2d mesh of the lake surface.  So one must first reads a 2d mesh, construct a 3d mesh from the 2d mesh and the depth and solves the potential flow equation.


## Step 1: Build a 3d mesh from a 2d mesh and and a depth function

To display the results we will call Medit so we need to load the corresponding module.
~~~freefem
load "medit"
~~~
Verbosity adjusts the messages sent by FreeFem++.
~~~freefem
verbosity=3;
~~~
Let us read the 2d mesh contained in file ("lac-leman-v4.msh". The format is
- On the first line : the number of vertices, the number of triangles and the number of boundary edges
- Then the coordinates of the vertices and their region number
- Then the identification number (i.e; the position in the above list) of the 3 vertices of each triangle
- Finally the connectivity of the boundary edges: the 2 vertices and the label number of the edge
~~~freefem
mesh Th2("lac-leman-v4.msh");
fespace Vh2(Th2,P1);
Vh2 d;
~~~
The 2 lines above have defined a Finite Element Space on the 2d mesh with $P^1$ triangular elements.

To construct a 3d mesh we need the bathymetry of the lake, i.e. a function for the lake's depth at each point of the surface, $d(x,y),~(x,y)\in \Omega_2$.
This is done in a very approximate way as the result of a Laplace equation
$$
\left\{
\begin{align*}
	\Delta d  &= 1 \hbox{ in }\Omega_2
  \\
	d|_{\partial\Omega_2} &= 0
\end{align*}
\right.
$$
In variational form, this is
~~~freefem
{  Vh2 v; 
	macro Grad(u) [dx(u),dy(u)] //
	solve P(d,v)= int2d(Th2)(Grad(d)'*Grad(v))+int2d(Th2)(v)
	+on(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,d=-1);
	d = d*5/abs(d[].min);
	plot(d,wait=1);
}
~~~
The macro is here to make the formula easier to read.
There are many labels for $\partial\Omega$ each corresponding to a portion of the lake boundary betwen two rivers.
The code has been surrounded by braces so as to make all variables local. This is by no means compulsory but this way one can reuse the names of the variables.

#### Build the 3d mesh using "buildlayers"
The mesh will have nn=10 layers:
~~~freefem
int[int] rup=[0,200],  rdown=[0,100];
real maxdeep = d[].min;
int nn=10; // controls the mesh fineness
mesh3 Th=buildlayers(Th2,nn,
  coef= d/maxdeep,
  zbound=[d,0],
  reffaceup = rup,
  reffacelow = rdown);
~~~
The purpose of rup and rdown is to replace 0 by 200 and 0 by 100 in the labels for the faces of the elements which are on the bottom of the lake and on the surface. For more details see
https://doc.freefem.org/documentation/mesh-generation.html#the-command-buildlayers
The mesh is displayed with "medit" (type "q" to exit medit)
### Step2: Solve the problem
As usual we need to work with a finite element space on the 3d mesh and define some functions
~~~freefem
fespace Vh(Th,P1);
Vh p,q;
macro Grad(u) [dx(u),dy(u),dz(u)] //

~~~
In the following $\texttt{din}$ and $\texttt{dout}$ are the debits from borders 1 and borders 2, in and out of the lake.  Some rescaling is needed to satisfy the necessary condition for existence of a solution.
~~~freefem
real ain=int2d(Th,1)(1.);
real aout=int2d(Th,2)(1.);
cout << " area " << ain << " " << aout << endl;
real din=1./ain; // inflow
real dout=-1./aout; // ouflow

~~~
Now the problem is solved and the solution is displayed statically by plot and interactively by medit:
~~~freefem
solve P(p,q)= int3d(Th)(Grad(p)'*Grad(q)+1e-5*p*q)-int2d(Th,1)(q*din)+int2d(Th,2)(q*dout);
plot(p,wait=1,nbiso=30,value=1);
/* if(!NoUseOfWait) medit("potentiel",Th,p,wait=1);*/
~~~
## Results

| The depth function d:|
|----------------------|
|![][_depth]           |

| The potential function p: the flow is slow except near Geneva at the tip of the lake.|
|--------------------|
|![][_flow]          |


[_depth]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/3d-leman/bathymetry.png

[_flow]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/3d-leman/flow.png
