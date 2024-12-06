---
name: condensator
category: electrostatics
layout: example
---

# The Black-Scholes equation for financial derivatives

In the documentation,
https://doc.freefem.org/tutorials/heatExchanger.html#heat-exchanger

It is quite easy to compute the electrostatic potential $u$ between 2 plates, one at potential -1 and the other at potential 1.  The plates $C_1,C_3$ are inside a round container $C_2$.

The geometry is a bit painful to describe with the keyword border.
~~~freefem
int C1 = 1; // labels to identify the 2 plates C1,C3 and the containner C2
int C2 = 2; // could be any number
int C3 = 3;
// Mesh
border C0(t=0, 2*pi){x=5*cos(t); y=5*sin(t); label=C2;}

border C11(t=0, 1){ x=1+t;  y=3;      label=C1;}
border C12(t=0, 1){ x=2;    y=3-6*t;  label=C1;}
border C13(t=0, 1){ x=2-t;  y=-3;     label=C1;}
border C14(t=0, 1){ x=1;    y=-3+6*t; label=C1;}

border C21(t=0, 1){ x=-2+t; y=3;      label=C3;}
border C22(t=0, 1){ x=-1;   y=3-6*t;  label=C3;}
border C23(t=0, 1){ x=-1-t; y=-3;     label=C3;}
border C24(t=1, 0){ x=-2;   y=-3+6*t; label=C3;}

mesh Th=buildmesh(C0(50)
    +C11(5)+C12(20)+C13(5)+C14(20)
    +C21(5)+C22(20)+C23(5)+C24(-20));
plot(Th, wait=true);
~~~

![][_mesh]

The electrostatic potential equation in absence of volumic charges is
$$
\Delta u =1,\quad u|_{C_1}=-1,\quad u|_{C_2}=0,\quad u|_{C_3}=1,
$$
The variational formulation is: find $u\in H^1(\Omega)$ with $u|_{C_i}=i-2, i=1,2,3$,
$$
\int_\Omega\nabla u\cdot\nabla v=0\quad \forall v\in H^1_0(\Omega).
$$
~~~freefem
fespace Vh(Th, P1);
Vh u, v;

// Problem
solve a(u, v)
    = int2d(Th)(
          dx(u)*dx(v) + dy(u)*dy(v)
    )
    + on(C1, u=1) + on(C2, u=0) + on(C3, u=-1);

// Plot
plot(u, value=true, ps="condersor.ps");
~~~
Notice that the solution is stored in a file condensor.eps in postscript format next to the program file. The extension .ps means it is a postcript format (don't use .eps at this elvel). To insert it in a publication you may use the linux command
    ps2eps condersor.ps
    epstopdf condensor.eps

![][_solution]

[_mesh]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/condensor/mesh.png

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/condensor/solution.png
