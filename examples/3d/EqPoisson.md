---
name: EqPoisson
category: math
folder: 3d
---

## Laplace equation in 3D on a object of revolution

The domain $\Omega$ is obtained by an axisymmetric extrusion of a 2d domain which has a hole
~~~freefem
load "medit"
func f=2*((0.1+(((x/3))*(x-1)*(x-1)/1+x/100))^(1/3.)-(0.1)^(1/3.));
real yf=f(1.2,0); 
border up(t=1.2,0.){ x=t;y=f;label=0;}
border axe2(t=0.2,1.15) { x=t;y=0;label=0;}
border hole(t=pi,0) { x= 0.15 + 0.05*cos(t);y= 0.05*sin(t); label=1;}
border axe1(t=0,0.1) { x=t;y=0;label=0;}
border queue(t=0,1) { x= 1.15 + 0.05*t; y = yf*t; label =0;}
int np= 100;
func bord= up(np)+axe1(np/10)+hole(np/10)+axe2(8*np/10)+ queue(np/10);
plot( bord); 
mesh Th2=buildmesh(bord);
plot(Th2,wait=1);
int[int] l23=[0,0,1,1]; 
mesh3 Th=buildlayers(Th2,coef= max(.15,y/max(f,0.05)), 50 ,zbound=[0,2*pi]
   ,transfo=[x,y*cos(z),y*sin(z)],facemerge=1,labelmid=l23);
~~~

The problem is to find $u\in V_1$ such that
$$
\int_\Omega \nabla u\cdot\nabla v = \int_\Omega v \quad \forall v\in V_0
$$
$$
V_w=\{v\in H^1(\Omega),~v|_S=w\}.
$$
where $S$ is the first border defined above, called "up"
~~~freefem
macro Grad(u) [dx(u),dy(u),dz(u)] //
fespace Vh(Th,P1);  Vh u,v;
solve Poisson(u,v) = int3d(Th)( Grad(u)'*Grad(v) ) - int3d(Th)( v) + on(1,u=1);
plot(u,wait=1,nbiso=20,value=1);
medit("u",Th,u,wait=1);
~~~

| The initial mesh and the deformed mesh |
|------------------------|
|![][_solution]          |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/3d/EqPoisson/solution.png