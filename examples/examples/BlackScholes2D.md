---
name: BlackScholes2D
category: finance
layout: example
---

# The Black-Scholes equation for financial derivatives

In the documentation,
https://doc.freefem.org/models/evolution-problems.html#d-black-scholes-equation-for-an-european-put-option

the following modeling equation for an European Put $u$ is to be integrated in  $( 0 , T ) \times R^+ \times R$:
$$
\left\{
\begin{align*}&
\partial_t u + \frac{(\sigma_1 x)^2}2\frac{\partial^2 u}{\partial x^2} + \frac{(\sigma_2 y)^2}2\frac{\partial^2 u}{ \partial y^ 2} + \sigma_1\sigma_2 x y \frac{\partial^2 u}{\partial x \partial y} + r x \frac{\partial u}{\partial x} + r y \frac{\partial u}{ \partial y} +r u = 0 ,
\cr&
u ( x , y , T ) = ( K - \max ( x , y ) )^+
\end{align*}
\right.
$$
The 2 underlying asset (e.g. Renault and Peugeot shares) are valued at $x$ and $y$.

The interest rate is $r$, the volatilities are $\sigma_1,\sigma_2$, the correlation is $\rho$. $K$ is the strike price.
For example,
~~~freefem
real sigma1 = 0.3;
real sigma2 = 0.3;
real rho = 0.3;
real r = 0.05;
real K = 40;
~~~
Notice that the problem is backward in time because final conditions are given and the signs of the differential operator is positive.  We will change $t$ to $T-t$ to work forward in tome.

As usual we must work on the variational formulation: find $u\in V:=  L^2(0,T,H^1(\Omega))$ such that
$$
\begin{aligned}&
\int_\Omega\left\{\left(\partial_t u + r u+ ( -x r + x\sigma_1^2 + \frac12 x\rho\sigma_1\sigma_2)\partial_x u +(-y r + y\sigma_2^2 + \frac12 y\rho\sigma_1\sigma_2)\partial_y u\right) v
\right.\cr&\left. + \frac{(\sigma_1 x)^2}2\partial_x u\partial_x v + \frac{(\sigma_2 y)^2}2{\partial_y u}{\partial_y v} + \frac{ρ}2\sigma_1\sigma_2 x y[ {\partial_x u}{\partial_y v}  +{\partial_y u}{\partial_x v}]\right\} = 0\quad \forall v\in V,
\cr&
u ( x , y , 0) = ( K − \max ( x , y ) )^+
\end{aligned}
$$
Boundary conditions  are not needed at $x=0$ and $y=0$ because the differential operator degenerates. At infinity there are implicit Neumann conditions.  In practice the problem must be localized, so $R^+\times R^+$ is replaced by $(0,L)\times(0,LL)$. Hence the mesh and the FEM space are  built as follows:
~~~freefem
int m = 30;
int L = 80;
int LL = 80;
int j = 100;
mesh th = square(m, m, [L*x, LL*y]);
// Fespace
fespace Vh(th, P1);
Vh u = max(K - max(x, y), 0.);
Vh xveloc, yveloc, v, uold;
~~~
An implicit in time Euler scheme is used, but in addition for stability and avoiding upwinding, we use the method of characteristics:
$$
\begin{aligned}&
\partial_t u + r u+ ( -x r + x\sigma_1^2 + \frac12 x\rho\sigma_1\sigma_2)\partial_x u +(-y r + y\sigma_2^2 + \frac12 y\rho\sigma_1\sigma_2)\partial_y u
\cr&
\approx
\frac{1}{\delta t}(u^{n+1}-u^n(x-U\delta t,y-V\delta t) ) \text{ where }
\cr&
U=-x r + x\sigma_1^2 + \frac12 x\rho\sigma_1\sigma_2,
\quad
V=-y r + y\sigma_2^2 + \frac12 y\rho\sigma_1\sigma_2.
\end{aligned}
$$
So the implementation is
~~~freefem
real dt = 0.01;
for (int n = 0; n*dt <= 1.0; n++) {
    xveloc = -x*r + x*sigma1^2 + x*rho*sigma1*sigma2/2;
    yveloc = -y*r + y*sigma2^2 + y*rho*sigma1*sigma2/2;
    // Update
    uold = u;
    // Solve
    solve eq1(u, v, init=j, solver=LU)
        = int2d(th)(
              u*v*(r + 1/dt)
            + dx(u)*dx(v)*(x*sigma1)^2/2
            + dy(u)*dy(v)*(y*sigma2)^2/2
            + dy(u)*dx(v)*rho*sigma1*sigma2*x*y/2
            + dx(u)*dy(v)*rho*sigma1*sigma2*x*y/2
        )
        + int2d(th)(
            - v*convect([xveloc, yveloc], dt, uold)/dt
        )
        + on(2,3,u=0)
        ;
}
~~~
Notice the parameters in solve to indicate that the factorized matrix can be reused.  With the default solver of FreeFem this is no longer much of an optimization.
Visualization is as usual
~~~freefem
plot(u, fill=true, wait=true, value=true);
~~~
One could add mesh adaptation by inserting above the line

// update

the following:
~~~freefem
   // Mesh adaptation
    j = j + 1;
    if (j > 20) {
        th = adaptmesh(th, u, verbosity=1, abserror=1, nbjacoby=2,
        err=0.001, nbvx=5000, omega=1.8, ratio=1.8, nbsmooth=3,
        splitpbedge=1, maxsubdiv=5, rescaling=1) ;
        j = 0;
        u = u;
        plot(th, wait=true);
    }
~~~

## Results
| The price of the Put |
| -------------------- |
| ![][_solution]       |

| Mesh adapted to the solution |
| ---------------------------- |
| ![][_adaptedmesh]            |

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/BlackScholes2D/solution.png
[_adaptedmesh]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/BlackScholes2D/adaptedmesh.png
