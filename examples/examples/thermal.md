---
name: thermal
category: thermodynamics
layout: example
---

# Heat Equation
The domain is a disk with 2 rectangles inside. One is a hole the other one is a region with a different heat diffusion parameter.
$$
\kappa= 1 + 4{\bf 1}_{(x<-1)\cap(x>-2)\cap(y<3)\cap(y>-3)}.
$$
~~~freefem
border C0(t=0, 2*pi){x=5*cos(t); y=5*sin(t);}
int C1 = 99, C2 = 98; // could be anything
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
plot(Th, wait=1);

// Fespace
fespace Vh(Th, P1);
Vh u, v;
Vh kappa = 1 + 4*(x<-1)*(x>-2)*(y<3)*(y>-3);
~~~

| The mesh   |
| ---------- |
| ![][_mesh] |

The heat equation is
$$
-\nabla\cdot(\kappa\nabla u)=0,~~u|_{C_0}=20,~~u|_{C_1}=100.
$$
~~~freefem
solve a(u, v)
  = int2d(Th)(kappa*(dx(u)*dx(v) + dy(u)*dy(v)))
  + on(C0, u=20)
  + on(C1, u=100);
plot(u, value=true, wait=1, fill=true);
~~~

| The temperature |
| --------------- |
| ![][_solution]  |

[_mesh]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/thermal/mesh.png

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/thermal/solution.png
