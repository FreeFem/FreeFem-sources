---
name: potential
category:fluid mechanics
layout: example
---
## Potential flow around an airfoil with thermal effects.

Potential flow is
$$
\nabla\times\vec u=0,~~\nabla\cdot \vec u=0~~,i.e.~~\vec u=\nabla\times\psi,~~\Delta\psi=0.
$$
The flow is tangent to the airfoil, so $~~\vec u\cdot n=0,~i.e.~~\psi=a$ constant.
At infinity the flow is given: $\phi=y\cos\alpha-x\sin\alpha$ where $\alpha$ is the angle of incidence.  
The airfoil is a NACA0012 given as a function $x\to \pm y(x)$.
~~~freefem
real S = 99;

border C(t=0, 2*pi){x=3*cos(t); y=3*sin(t);} // Label 1,2 
border Splus(t=0, 1){x=t-0.5; y=0.17735*sqrt(t) - 0.075597*t - 0.212836*(t^2) + 0.17363*(t^3) - 0.06254*(t^4); label=S;}
border Sminus(t=1, 0){x=t-0.5; y=-(0.17735*sqrt(t) - 0.075597*t - 0.212836*(t^2) + 0.17363*(t^3) - 0.06254*(t^4)); label=S;}
mesh Th = buildmesh(C(50) + Splus(70) + Sminus(70));
// Fespace
fespace Vh(Th, P2);
~~~
Joukowski tells us that the pressure $p=-\frac12|\vec u|^2$ is continuous at the trailing edege P. As $\vec u\cdot\vec n|_P=0$, it implies that the tangent velocity is continuous: the jump from the lower side to the upper side at P is zero: $[\vec u\cdot\vec s]|_P=0$.  The trick is to define
$$
\begin{align*}&
\Delta\psi^0=0,~~ \psi^0|_C=0,~~\psi^0|_\infty=y\cos\alpha-x\sin\alpha
\cr&
\Delta\psi^1=0,~~ \psi^1|_C=1,~~\psi^1|_\infty=0
\end{align*}
$$
and search by superposition for $\psi=\psi_0+\beta\psi_1$ and adjust $\beta$ to satisfy the Joukowsky condition.
See https://doc.freefem.org/models/static-problems.html#aerodynamics

Here we assume that $\beta$ is given and we concentral on thermal effects.

Let $\Omega$ be a the bounded open set of $R^2$ approximating $\infty$ by a circle and having the airfoil near the center of the circle.  Consider the variational formulations to apply the finite element method
$$
\int_\Omega\nabla\psi\nabla\hat\psi =0~~\forall \hat\psi\in  H^1_0(\Omega);~~\psi-\psi_\Gamma^i\in H^1_0(\Omega).
$$
~~~freefem
Vh psi, w;
real cost = cos(5.*pi/180.), sint=sin(5.*pi/180.);// incidence 5 degres
// Problem
solve potential(psi, w)
  = int2d(Th)(dx(psi)*dx(w)+dy(psi)*dy(w))
  + on(C, psi = cost*y-sint*x) 
  + on(S, psi=0);

plot(psi, wait=1);
~~~
For the temperature equation we work with a different mesh on because the temperature varies also inside the airfoil.
~~~freefem
border D(t=0, 2.){x=0.5+t*cost; y=+t*sint;}
mesh Sh = buildmesh(C(25) + Splus(-90) + Sminus(-90) + D(200));
int steel = Sh(0.5, 0).region, air = Sh(-1, 0).region;
// Change label to put BC on In flow 
fespace Wh(Sh, P1);
Wh  vv;

fespace W0(Sh, P0);
W0 k = 0.01*(region == air) + 0.1*(region == steel);
W0 u1 = dy(psi)*(region == air), u2 = -dx(psi)*(region == air);
Wh v = 120*(region == steel), vold;
// pul label 10 on inflow boundary to inforce the temperature.
Sh = change(Sh,flabel = (label == C &&  [u1,u2]'*N<0) ? 10 : label);
~~~
The time dependent heat equation is solved by an implicit Euler time scheme.
$$
\frac1{dt} v -\nabla\cdot(k\nabla v) + a(u_1\partial_x v + u_2\partial_y v) =\frac1{dt} v_{old}
$$
where $[u_1,u_2]^T=\nabla\times\psi$.
~~~freefem
int i;
real dt = 0.005, nbT = 50;
problem thermic(v, vv, init=i, solver=LU)
  = int2d(Sh)(
      v*vv/dt + k*(dx(v)*dx(vv) + dy(v)*dy(vv))
    + 10*(u1*dx(v) + u2*dy(v))*vv
  )
  - int2d(Sh)(vold*vv/dt)
  + on(10, v= 0);
  

for(i = 0; i < nbT; i++) {
    vold[]= v[];
    thermic;
    plot(v);
}
plot(v, wait=1,fill=1,value=1);
~~~
