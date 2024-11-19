---
name: heat
category: thermodynamics
layout: example
---
## Time dependent heat equation  with a discontinuous thermal diffusion
The time dependent heat equation  with a discontinuous thermal diffusion is integrated:
$$
\partial_t u - \nabla\cdot(k\nabla u)=0, ~~u|_{\Gamma_1}=u_\Gamma,~~k\frac{\partial u}{\partial n}|_{\Gamma_2}+k_f(u-u_e)=0,~~u_{t=0}=0.
$$
In this example the domain $\Omega$ is a rectangle $(0,3)\times(0,1)$; $k=1.8~{\bf 1}_{y<0.5}+0.2$, $k_f=1$ where ${\bf 1}_x$ is the Heaviside function.  The time varies from 0 to 5. Finally $u_e=20$. It corresponds to the FreeFem script
~~~freefem
func k = 1.8*(y < 0.5) + 0.2;
real kf = 1, ue = 20, T = 5, dt = 0.1;

// Mesh
mesh Th = square(150, 50, [3*x, y]);

// Fespace
fespace Vh(Th, P1);
Vh u, uold, v, usave;
~~~
An Euler imlicit scheme is used for the time. After time discretisation, the variational formulation for $u=u^{n+1}$ is
$$
\int_\Omega\frac u{\delta t} \hat u+ k\nabla u\nabla\hat u + \int_{\Gamma_2}k_f u \hat u =\int_\Omega\frac{u^n}{\delta t}\hat u + \int_{\Gamma_2}k_f u_e\hat u,~~~u|_{\Gamma_1}=u_\Gamma,
~~\forall \hat u~\hbox{ s.t. }\hat u_h|_{\Gamma_2}=0.
$$
~~~freefem
int kk=0;
problem Heat(u, v, init=kk)
  = int2d(Th)(
      u*v/dt
    + k*(dx(u)*dx(v) + dy(u)*dy(v))
  )
  + int1d(Th, 1, 3)(kf*v*u)
  - int1d(Th, 1, 3)(kf*v*ue)
  - int2d(Th)(uold*v/dt)
  + on(2, 4, u=30);

bool withplot = 0;
real cpu1, cpu2;
real cpu = clock();

// Initialization
u = 0;

// Basic time loop
for (real t = 0; t < T; t += dt) {
  uold = u;
  Heat;
  kk++;
  if (withplot) plot(u);
}

cpu1 = clock()-cpu;
plot(u);
usave[] = u[];
~~~
$\texttt{cpu}$ gives a measure of the computing time. $\texttt{usave}$ will be used later to compare with iteration integration method.  Notice that $u[]$ addresses the array of values at the vertices of $u$.

Now the same problem is solved using an explicit construction of the linear system.  We begin with the matrix A:
~~~freefem
varf vA(u, v) = int2d(Th)( u*v/dt
    + k*(dx(u)*dx(v) + dy(u)*dy(v))
    )
  + int1d(Th, 1, 3)(kf*v*u)
  + on(2, 4, u=30);
 matrix A = vA(Vh, Vh, solver=sparsesolver);
~~~
The righthandside is also a matrix B times a vector plus another vector vL.
~~~freefem
varf vB(u, v) = int2d(Th)(u*v/dt) ;
varf vRHS(u, v) = int1d(Th, 1, 3)(kf*v*ue);
varf vL(u, v) = on(2, 4, u=30);

  cpu = clock();
  real[int] rhsbc = vL(0, Vh);
  real[int] rhs0 = vRHS(0, Vh);
  matrix B = vB(Vh, Vh);

  // Initialization
  u = 0;

  // Optimized time loop (speed of C language)
  for (real t = 0; t < T; t += dt) {
    real[int] b = B*u[];
    b += rhs0;
    b = rhsbc ? rhsbc : b;
    u[] = A^-1*b;
    if(withplot) plot(u);
  }
  cpu2 = clock() - cpu;

plot(u, cmm="u2");
cout << " cpu method 1 = " << cpu1 << " cpu method matrix = " << cpu2 << " ratio = " <<  cpu1/cpu2 << endl;
~~~
![](https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/Heat/solution.png)
