---
name: optimcontrol
category:optimization using BFGS
layout: example
---
# Optimal control with 3 parameters solved by BFGS

Let $\Omega$ be a bounded open set of $R^2$. Let $u_d$ and $w$ be given.  Consider
$$
\begin{align*}&
J(z)=\min_{z\in Z}\int_\Omega|u-u_d|^2 \hbox{where $u$ is solution of}
\cr&
\int_\Omega\nu(z)\nabla u\nabla\hat u =0~~~\forall \hat u\in  H^1_0(\Omega);~~u-w\in H^1_0(\Omega)
\cr&
\end{align*}
$$
where $\nu(z)$ is some given function of $z$.
$\texttt{BFGS}$ is an optimization module which can be called to find a local minimum of $J$ but it requires the derivative of $J$ with respect to $z$.

In this example $\Omega$ is a disk of radius 5.  The optimization space $Z$ is the set of functions which are constant on 3 smaller disks, $D_0,D_1,D_2$ of radius 1 in $\Omega$ and 1 outside these smaller disk.  Thus the dimension of $Z$ is 3:
$$
Z=\{z: ~ z(x,y)=1+z_0{\bf 1}_{D_0}+z_0{\bf 1}_{D_1}+z_0{\bf 1}_{D_2}\}
$$
where ${\bf 1}_{D_i}$ is the characteristic function of $D_i$.
In the follwoing example, $w=x^3-y^3$.
~~~freefem
border aa(t=0, 2*pi) {x = 5*cos(t); y = 5*sin(t);}
border bb(t=0, 2*pi) {x = cos(t); y = sin(t);}
border cc(t=0, 2*pi) {x = -3+cos(t); y = sin(t);}
border dd(t=0, 2*pi) {x = cos(t); y = -3+sin(t);}
mesh th = buildmesh(aa(70) + bb(35) + cc(35) + dd(35));

// Fespace
fespace Vh(th, P1);
Vh Ib = ((x^2+y^2)<1.0001),
   Ic = (((x+3)^2+ y^2)<1.0001),
   Id = ((x^2+(y+3)^2)<1.0001),
   Ie = (((x-1)^2+ y^2)<=4),
   ud, u, uh, du;

// Problem
real[int] z(3);
problem A(u, uh)
  = int2d(th)(
    (1 + z[0]*Ib + z[1]*Ic + z[2]*Id)*(dx(u)*dx(uh) + dy(u)*dy(uh))
  )
  + on(aa, u=x^3-y^3);

// A test
z[0] = 2; z[1] = 3; z[2] = 4;
A;
plot(u, wait=1);
~~~
The test will also be our target
![][_test]
BFGS requires $J$ and $J'_z$ to be defined as C-functions
~~~freefem
ud = u;
ofstream f("J.txt");

func real J(real[int] & Z) {
    for (int i = 0;i < z.n; i++) z[i] = Z[i];
    A;
    real s = int2d(th)(Ie*(u-ud)^2);
    f << s << "   "; // so that every time J is called this is printed
    return s;
}

real[int] dz(3), dJdz(3);

problem B(du, uh)
  =int2d(th)(
    (1 + z[0]*Ib + z[1]*Ic + z[2]*Id)*(dx(du)*dx(uh) + dy(du)*dy(uh))
  )
  +int2d(th)(
    (dz[0]*Ib + dz[1]*Ic + dz[2]*Id)*(dx(u)*dx(uh) + dy(u)*dy(uh))
  )
  +on(aa, du=0);

func real[int] DJ(real[int] &Z) {
  for(int i = 0; i < z.n; i++) {
    for(int j = 0; j < dz.n; j++)
      dz[j]=0;
    dz[i] = 1;
    B;
    dJdz[i] = 2*int2d(th)(Ie*(u-ud)*du);
  }
  return dJdz;
}
~~~
We are now ready to call BFGS with initial value for the optimization parameter $Z[i]=1$, $i=0,1,2$.
~~~freefem
real[int] Z(3);
for(int j = 0; j < z.n; j++) Z[j] = 1;

BFGS(J, DJ, Z, eps=1.e-6, nbiter=15, nbiterline=20);
~~~
The results are printed and plotted:
~~~freefem
cout << "BFGS: J(z) = " << J(Z) << endl;
for(int j = 0; j < z.n; j++) cout << z[j] << endl;
plot(ud, value=1, ps="u.ps");
~~~

| The solution   |
| -------------- |
| ![][_solution] |

[_test]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/optimcontrol/test.png

[_solution]: https://raw.githubusercontent.com/phtournier/ffmdtest/refs/heads/main/figures/examples/optimcontrol/solution.png
