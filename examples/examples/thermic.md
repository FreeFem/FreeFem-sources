---
name: thermal
category: thermodynamics
layout: example
---

# The time dependent nonlinear heat equation

The time dependent heat equation  with a discontinuous thermal diffusion and nonlinear dissipation is integrated in a rectangle:
~~~freefem
mesh Th = square(30, 5, [6*x, y]);
~~~
$$
\partial_t u - \nabla\cdot(k\nabla u)=0, ~~u|_{\Gamma_1}=u_\Gamma,~~k\frac{\partial u}{\partial n}|_{\Gamma_2}+b(u)(u-u_e)=0,~~u_{t=0}=u_0.
$$
In this example $k=1.8~{\bf 1}_{y<0.5}+0.2$ where ${\bf 1}_x$ is the Heaviside function.  The time varies from 0 to 5. Finally $u_e=20$.  The function $b$ corresponds to the linearization of a $T^4$ law.
~~~freefem
func u0 = 10+90*x/6;
func k = 1.8*(y<0.5) + 0.2;
real ue = 25, alpha = 0.25, T = 5, dt = 0.1;
real rad = 1e-8, uek = ue + 273.;

// Fespace
fespace Vh(Th,P1);
Vh vold, w, v=u0-ue, b;

problem thermradia(v, w)
  = int2d(Th)(
      v*w/dt
    + k*(dx(v)*dx(w) + dy(v)*dy(w))
  )
  + int1d(Th, 1, 3)(b*v*w)
  - int2d(Th)(vold*w/dt)
  + on(2, 4, v=u0-ue);
~~~
An implicit Euler scheme is used to integrate the equation in time with time step 0.1. A loop is used for the nonlinearity.
~~~freefem
for(real t = 0; t < T; t+=dt) {
  vold = v;
  for (int m = 0; m < 5; m++) {
    b = alpha + rad*(v + 2*uek)*((v+uek)^2 + uek^2);
    thermradia;
  }
}
vold = v + ue;

// Plot
plot(vold);
~~~

| The temperature   |
| ----------------- |
| ![][_solution]    |

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/thermic/solution.png
