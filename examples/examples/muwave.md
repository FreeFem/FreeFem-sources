---
name: Muwave
category:Electromagnetism
layout: example
---

## Compute the stationary electromagnetic wave for a microwave oven

The domain is a square box with a rectangle inside to simulate the object (called meast below) which is heated by the wave.  The electromagnetic source is from the top boundary
~~~freefem
int n = 2;
real a = 20, b = 20, c = 15, d = 8, e = 2, l = 12, f = 2, g = 2;

// Mesh
border a0(t=0, 1) {x=a*t; y=0; label=1;}
border a1(t=1, 2) {x=a; y= b*(t-1); label=1;}
border a2(t=2, 3) {x=a*(3-t); y=b; label=1;}
border a3(t=3, 4) {x=0; y=b-(b-c)*(t-3); label=1;}
border a4(t=4, 5) {x=0; y=c-(c-d)*(t-4); label=2;}
border a5(t=5, 6) {x=0; y= d*(6-t); label=1;}

border b0(t=0, 1) {x=a-f+e*(t-1); y=g; label=3;}
border b1(t=1, 4) {x=a-f; y=g+l*(t-1)/3; label=3;}
border b2(t=4, 5) {x=a-f-e*(t-4); y=l+g; label=3;}
border b3(t=5, 8) {x=a-e-f; y=l+g-l*(t-5)/3; label=3;}

mesh Th = buildmesh(a0(10*n) + a1(10*n) + a2(10*n) + a3(10*n) + a4(10*n) + a5(10*n)
  + b0(5*n) + b1(10*n) + b2(5*n) + b3(10*n));
~~~

| The geometry   |
| -------------- |
| ![][_geometry] |

The wave is solution of a PDE in the complex domain for which the variational formulation is
$$
\begin{align*}
\int_\Omega(
      (1+R)v w
    -(1-0.5i)(\partial_x v\partial_x w+ \partial_y v \partial_y w )=0
   \\
     ~~~~~ \forall  w~:~w|_{\Gamma_1\cup\Gamma_2}=0,
    \\
    v|_{\Gamma_1}=0,~~~~
 v|_{\Gamma_2}=\sin(\pi\frac{y-c}{c-d})
    \end{align*}
$$
and where $R$ takes a different value in the meat. This is implemented as follows
~~~freefem
real meat = Th(a-f-e/2, g+l/2).region,
     air= Th(0.01, 0.01).region;
plot(Th, wait=1);

fespace Vh(Th, P1);
Vh R = (region-air)/(meat-air);
Vh<complex> v, w;

solve muwave(v, w)
  = int2d(Th)(
      v*w*(1+R)
    -(dx(v)*dx(w)+dy(v)*dy(w))*(1-0.5i)
  )
  + on(1, v=0)
  + on(2, v=sin(pi*(y-c)/(c-d)));

Vh vr = real(v), vi = imag(v);
plot(vr, wait=1, ps="rmuonde.ps", fill=true);
plot(vi, wait=1, ps="imuonde.ps", fill=true);
~~~
Notice that the linear system is complex and could be singular (resonance).

| Real part      |
| -------------- |
| ![][_realpart] |

| Imaginary part    |
| ----------------- |
| ![][_imaginepart] |

To compute the temperature in the meat we solve
$$
\begin{align*}
\int_\Omega(
     (\partial_x u\partial_x w+ \partial_y u \partial_y w )= \int_\Omega|v|^2w
   \\
     ~~~~~ \forall  w~:~w|_{\Gamma_1\cup\Gamma_2}=0,
    \\
    v|_{\Gamma_1\cup\Gamma_2}=0,
    \end{align*}
$$
~~~freefem
fespace Uh(Th,P1);
Uh u, uu, ff=1e5*(vr^2 + vi^2)*R;

solve temperature(u, uu)
  = int2d(Th)(
      dx(u) * dx(uu)
    + dy(u) * dy(uu)
  )
  - int2d(Th)(ff*uu)
  + on(1, 2, u=0);

plot(u, wait=1, ps="tempmuonde.ps", fill=true);
~~~

| The temperature   |
| ----------------- |
| ![][_heat] |

[_geometry]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/muwave/geometry.png

[_realpart]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/muwave/realpart.png

[_imaginepart]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/muwave/imaginepart.png

[_heat]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/muwave/heat.png
