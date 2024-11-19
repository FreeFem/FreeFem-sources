---
name: Lamé
category: Solid Mechanics
layout: example
---

## Compute the deflection of a 2d elastic beam
The geometry is a rectangle, clamped on the left vertical side and pulled down by itsown weight.
The equations are
$$
\displaystyle{
	\begin{align*}&
  \int_{\Omega}[2\mu~\epsilon(u)\cdot\epsilon(v) +\lambda (\nabla\cdot u )(\nabla\cdot v)]  = \int_{\partial\Omega}{f u_2},
  \quad \forall v\in V(\Omega)
  \cr&
  \nabla\cdot u= \partial_x u_1+\partial_y u_2,
  \cr&
  \epsilon(u)=[\partial_x u_1,\partial_y u_2, \frac1{\sqrt{2}}(\partial_y u_1 + \partial_x u_2)].
  \end{align*}
}
$$
where $V(\Omega)=\{ v\in (H(\Omega)^2): v|_S=0\}$ and $S\subset\partial\Omega$ is the portion of the boundary where the beam is clamped.
$\mu$ and $\lambda$ are computed from the Young and Poisson constants $E,\sigma$,
$$
\mu = \frac{E}{(2(1 + \nu))}, \quad
\lambda =\frac{E~\nu}{((1 + \nu)(1 - 2\nu))}
$$


~~~freefem
real E = 21e5, nu = 0.28;
real f = -1;
real mu = E/(2*(1 + nu));
real lambda = E*nu/((1 + nu)*(1 - 2*nu));

mesh Th = square(10, 10, [20*x, 2*y-1]);

fespace Vh(Th, P2);
Vh u, v, uu, vv;

real sqrt2 = sqrt(2.);
macro epsilon(u1, u2) [dx(u1), dy(u2), (dy(u1)+dx(u2))/sqrt2] // EOM
macro div(u, v) (dx(u) + dy(v)) // EOM

solve lame([u, v], [uu, vv])
  = int2d(Th)(  lambda*div(u,v)*div(uu,vv)
              + 2.*mu*(epsilon(u,v)'*epsilon(uu,vv)) )
  - int2d(Th)(f*vv) + on(4, u=0, v=0);

real coef=100;
plot([u, v], wait=1, ps="lamevect.ps", coef=coef);
~~~

| The displacement vectors |
| ------------------------ |
| ![][_solution]           |

A better way to display the result is to move the mesh by the displacement $[u,v]^T$.
~~~freefem
mesh th1 = movemesh(Th, [x+u*coef, y+v*coef]);
plot(th1, wait=1, ps="lamedeform.eps");

real dxmin = u[].min;
real dymin = v[].min;
cout << "   displacement  max x = " << dxmin << " y = " << dymin << endl;
cout << "   displacement (20,0) = " << u(20,0) << " " << v(20,0) << endl;
~~~

| The displaced beam |
| ------------------ |
| ![][_dispbeam]     |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/lame/solution.png

[_dispbeam]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/lame/dispbeam.png
