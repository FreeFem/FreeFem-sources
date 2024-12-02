---
name: stokes_periodic_composite
category: fluid
layout: example
---

# The Stokes system for creeping flow with periodic conditions

This is a continuation of the previous example, but now there are
- periodic boundary condition in x-direction
- Dirichlet boundary condition in y-direction
These are chosen from a manufacture solution
$$
u=\sin x\cos y,~v=-\cos x\sin y,~p=2\cos x\cos y
$$
~~~freefem
mesh Th=square(60,60,[2*pi*x,2*pi*y],flags=3);
mesh ThU=trunc(Th,1,split=2);
mesh ThP=Th;

fespace UhPerio(ThU,[P1],periodic=[[2,y],[4,y]]);
fespace Uh(ThU,[P1]);
fespace Ph(Th,P1);

fespace Vh=UhPerio*Uh*Ph; // definition of the composite FE space

cout << "ndof = " << Vh.ndof << endl;

UhPerio u1;
Uh u2;
Ph p;

func g1 = sin(x)*cos(y);
func g2 = -cos(x)*sin(y);

func f1 = 0;
func f2 = -4*cos(x)*sin(y);
func fp = 2*cos(x)*cos(y);

~~~
To illustrate a different method, here "solve" or "problem" is not used; the linear system is built with "varf"
~~~freefem
macro grad(u) [dx(u),dy(u)]//
macro Grad(u1,u2) [grad(u1), grad(u2)]//
macro div(u1,u2) (dx(u1)+dy(u2))//

varf Stokes ( [u1,u2,p], [v1,v2,q] )
= int2d(Th)( (Grad(u1,u2):Grad(v1,v2)) ) + int2d(Th)( - div(u1,u2)*q - div(v1,v2)*p ) + int2d(Th)( -1e-10*p*q )
+ int2d(Th) ( [f1,f2]'*[v1,v2] )
+ on(1,2,3,4,u2=g2) + on(1,3,u1=g1);

matrix A = Stokes(Vh,Vh);
real[int] b1 = Stokes(0,Vh);
real[int] sol = A^(-1)*b1;
~~~
To retrieve the velocity and the pressure from the compound vector sol, do
~~~freefem
[u1[],u2[],p[]]=sol;

plot( u1, cmm="u1" );
plot( u2, cmm="u2");
plot( p,  cmm="p" );
~~~
The rest of the script is devoted to the computation of the $L^2$ error
$$\|[u,v,p]-[u_e,v_e,p_e]\|.$$
~~~freefem
fespace VhU(ThU,P1);
fespace VhP(ThP,P1);

VhU ue1=g1;
VhU ue2=g2;
VhP pe=fp;

cout << " int2d(Th) (dx(fu1)+dy(fu2)) =" << int2d(Th)(dx(ue1)+dy(ue2)) << endl;

cout << "error L2 (u) = " << int2d(Th)( (u1-ue1)^2 + (u2-ue2)^2 ) << endl;
cout << "error L2 (p) = " << int2d(Th)( (p-pe)^2 ) << endl;
cout << "error L2 relative (u) = " << int2d(Th)( (u1-ue1)^2 + (u2-ue2)^2 ) / int2d(Th)( ue1^2+ue2^2 ) << endl;
cout << "error L2 relative (p) = " << int2d(Th)( (p-pe)^2 ) / int2d(Th)( pe^2 ) << endl;
~~~
