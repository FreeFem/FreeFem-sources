---
name: schwarz
category: mathematics
layout: example
---

## Domain Decomposition Method: solution by Schwarz iterations
The Laplace equation is solved in a domain which is the union of a circle and a rectangle.
~~~freefem
real n = 4; // Mesh quality

// Mesh
border a(t=0, 1){x=t; y=0;}
border a1(t=1, 2){x=t; y=0;}
border b(t=0, 1){x=2; y=t;}
border c(t=2, 0){x=t; y=1;}
border d(t=1, 0){x=0; y=t;}
border e(t=0, pi/2){x=cos(t); y=sin(t);}
border e1(t=pi/2, 2*pi){x=cos(t); y=sin(t);}

//Omega1 (rectangle)
mesh th = buildmesh(a(5*n) + a1(5*n) + b(5*n) + c(10*n) + d(5*n));
fespace Vh(th, P1);
Vh v, u=0;

//Omega2 (circle)
mesh TH = buildmesh(e(5*n) + e1(25*n));
fespace VH(TH,P1);
VH V, U=0;
~~~
The Schwarz algorithm solves iteratively and alternatively on the circle and the rectangle. The trace on the rectangle of the solution on the circle becomes the boundary conditon for thee problem on the rectangle and ssimilarly for the circle.
~~~freefem
for (int i = 0; i < 4; i++) {
  plot(U, u, wait=1, cmm="Iteration "+i);
  // Solve on Omega2
  solve AA(U, V)
    = int2d(TH)(dx(U)*dx(V) + dy(U)*dy(V))
    - int2d(TH)(V)
    + on(e, U=u)
    + on(e1, U=0);

  // Solve on Omega1
  solve aa(u, v)
    = int2d(th)(dx(u)*dx(v) + dy(u)*dy(v))
    - int2d(th)(v)
    + on(a, d, u=U)
    + on(a1, b, c, u=0);
}

plot(U, u, wait=1, cmm="Final solution");
~~~

| The solution at iteration 1 |
| --------------------------- |
| ![][_solone]                |

| The solution at iteration 4 |
| --------------------------- |
| ![][_solfour]               |

[_solone]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/schwarz/solone.png

[_solfour]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/schwarz/solfour.png
