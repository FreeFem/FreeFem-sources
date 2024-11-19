---
name: heatex
category: physics
layout: example
---
# The time independent  Heat Equation

The heat equation  with a discontinuous thermal diffusion is integrated:
$$
-\nabla\cdot(\kappa\nabla u)=0, ~~u|_{\Gamma_1}=u_\Gamma.
$$
In this example the domain $\Omega$ is a circle minus 2 rectangles on which boundary the temperature is imposed.
The thermal diffusion is $\kappa=1+4~{\bf 1}_{x\in(-2,-1)\times(-3,3)}$ where ${\bf 1}_x$ is the Heaviside function. Finally $u=20$ on the outer circle and 100 one the right rectangle. The left rectangle is there only to define $\kappa$. It corresponds to the following FreeFem script

~~~freefem
int C1 = 99, C2 = 98; // could be anything

// Mesh
border C0(t=0, 2*pi){x=5*cos(t); y=5*sin(t);}

border C11(t=0, 1){x=1+t;  y=3;      label=C1;}
border C12(t=0, 1){x=2;    y=3-6*t;  label=C1;}
border C13(t=0, 1){x=2-t;  y=-3;     label=C1;}
border C14(t=0, 1){x=1;    y=-3+6*t; label=C1;}

border C21(t=0, 1){x=-2+t; y=3;      label=C2;}
border C22(t=0, 1){x=-1;   y=3-6*t;  label=C2;}
border C23(t=0, 1){x=-1-t; y=-3;     label=C2;}
border C24(t=0, 1){x=-2;   y=-3+6*t; label=C2;}

mesh Th=buildmesh(C0(50)
  + C11(5) + C12(20) + C13(5) + C14(20)
  + C21(-5) + C22(-20) + C23(-5) + C24(-20));
plot(Th, wait=1, ps="heatexTh.ps");

// Fespace
fespace Vh(Th, P1);
Vh u, v;
Vh kappa = 1 + 4*(x<-1)*(x>-2)*(y<3)*(y>-3);

// Problem
solve a(u, v)
  = int2d(Th)(
    kappa*(dx(u)*dx(v) + dy(u)*dy(v))
  )
  + on(C0, u=20)
  +on(C1, u=100);

// Plot
plot(u, value=true, wait=1, fill=true);
~~~

| The solution   |
| -------------- |
| ![][_solution] |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/heatex/solution.png
