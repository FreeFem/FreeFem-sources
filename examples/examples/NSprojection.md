---
name: NSprojection
category: Fluid Mechanics
layout: example
---

## Time dependent incompressible Navier-Stokes Equations solved with Newton's method

The Navier-Stokes equations are solved for a flow over a backward facing step.
The time independent Navier-Stokes equations for an incompressible fluid are
$$
\frac{\partial u}{\partial t}-\nu\Delta u +u\cdot\nabla u -\nabla p =0,~~~\nabla\cdot u=0~~~in ~ \Omega
$$

The velocity is specified on the left vertical side of the pipe (inlet) and free on the vertical side on the right (outlet). There is a noslip condition on the lateral walls.
The geometry is as follows. Note the trick $t^{1.2}$ to refine the mesh near the corner of the step.
~~~freefem
verbosity=0;
int n = 1;

// Mesh
border a0(t=1, 0){x=-2; y=t; label=1;} // inlet
border a1(t=-2, 0){x=t; y=0; label=2;}
border a2(t=0, -0.5){x=0; y=t; label=2;}
border a3(t=0, 1){x=18*t^1.2; y=-0.5; label=2;}
border a4(t=-0.5, 1){x=18; y=t; label=3;} // outlet
border a5(t=1, 0){x=-2+20*t; y=1; label=4;}

mesh Th = buildmesh(a0(3*n) + a1(20*n) + a2(10*n) + a3(150*n) + a4(5*n) + a5(100*n));
plot(Th);

// Fespace
fespace Vh(Th, P1);
Vh w, u = 0, v = 0, p = 0, q = 0;

real nu = 0.0025; // Reynolds=200
real dt = 0.2;
real epsv = 1e-6, epsu = 1e-6, epsp = 1e-6;
~~~

| The mesh   |
| ---------- |
| ![][_mesh] |

The matrices dtMx and dtMy are used to project $[u,v]^T$ on the space of divergence free functions
~~~freefem
matrix dtM1x, dtM1y;

macro  BuildMat()
  { /* for memory managenemt */
    varf vM(unused, v) = int2d(Th)(v) ;
    varf vdx(u, v) = int2d(Th)(v*dx(u)*dt) ;
    varf vdy(u, v) = int2d(Th)(v*dy(u)*dt) ;

    real[int] Mlump = vM(0, Vh);
    real[int] one(Vh.ndof); one = 1;
    real[int] M1 = one ./ Mlump;
    matrix dM1 = M1;
    matrix Mdx = vdx(Vh, Vh);
    matrix Mdy = vdy(Vh, Vh);
    dtM1x = dM1*Mdx;
    dtM1y = dM1*Mdy;
  }// end of macro
BuildMat
~~~
In the follwing time loop, ${\bf u}=[u,v]^T$ are computed at every time steps by using the method of characteristics implemented in the operator $\texttt{convect}$:
$$
\partial_t{\bf u}+{\bf u}\cdot\nabla{\bf u}|_{x,t}\approx \frac1{\delta t}[{\bf u}(x,t)-{\bf u}(x-{\bf u}(x,t-\delta t)\delta t,t-\delta t)]
$$
~~~freefem
real err = 1, outflux = 1;
for(int n = 0; n < 200; n++) {
  Vh uold = u, vold = v, pold = p;

  solve pb4u(u, w, init=n, solver=CG, eps=epsu)
    = int2d(Th)(
      u*w/dt + nu*(dx(u)*dx(w) + dy(u)*dy(w))
    )
    - int2d(Th)(
      (convect([uold, vold], -dt, uold)/dt - dx(p))*w
    )
    + on(1, u=4*y*(1-y))
    + on(2, 4, u=0);
  plot(u);

  solve pb4v(v, w, init=n, solver=CG, eps=epsv)
    = int2d(Th)(
      v*w/dt + nu*(dx(v)*dx(w) + dy(v)*dy(w))
    )
    - int2d(Th)(
      (convect([uold, vold], -dt, vold)/dt - dy(p))*w
    )
    + on(1, 2, 3, 4, v=0);

  solve pb4p(q, w, solver=CG, init=n, eps=epsp)
    = int2d(Th)(dx(q)*dx(w) + dy(q)*dy(w))
    - int2d(Th)((dx(u) + dy(v))*w/dt)
    + on(3, q=0);

  // to have absolute epsilon in CG algorithm.
  epsv = -abs(epsv);
  epsu = -abs(epsu);
  epsp = -abs(epsp);
~~~
And then ${\bf u}$ is projected with $q$ by $M$ , i.e. ${\bf u}+M(\nabla q)\delta t$ is divergence free.
~~~freefem
p = pold-q;
u[] += dtM1x*q[];
v[] += dtM1y*q[];
~~~
For better precosion the mesh is adapted to the flow
~~~freefem
if(n%50 == 49) {
    Th = adaptmesh(Th, [u, v], q, err=0.06, nbvx=100000);
    plot(Th, wait=true);
    BuildMat // rebuild mat.
}
~~~
The stopping criteria to exit the loop:
~~~freefem
err = sqrt(int2d(Th)(square(u - uold) + square(v - vold))/Th.area);
  outflux = int1d(Th)([u, v]'*[N.x, N.y]) ;
  cout << " iter " << n << " Err L2 = " << err << " - Outflow = " << outflux << endl;
  if (err < 1e-3) break;
}
assert(abs(outflux) < 5e-3); // verification
plot(p, wait=1, ps="NSprojP.ps");
plot(u, wait=1, ps="NSprojU.ps");
~~~

| The u-component of the velocity |
| ------------------------------- |
| ![][_u]                         |

| The pressure |
| ------------ |
| ![][_p]      |

[_mesh]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/NSprojection/mesh.png

[_u]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/NSprojection/u.png

[_p]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/NSprojection/p.png
