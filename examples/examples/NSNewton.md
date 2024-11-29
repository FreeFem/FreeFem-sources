---
name: NSNewton
category: Fluid Mechanics
layout: example
---

## Stationary incompressible Navier-Stokes Equation solved with Newton's method
The time independent Navier-Stokes equations for an incompressible fluid are
$$
-\nu\Delta u +u\cdot\nabla u -\nabla p =0,~~~\nabla\cdot u=0~~~in ~ \Omega
$$
We are interested by the flow around a cylinder in an infinite fluid (approximated by a rectangle with a round side on the left and a circle in the center).
Dirichlet boundary conditons are imposed on boundaries 2 (the cylinder)  and 1, the left  (inflow) and lateral sides (tangent flow)  of the rectangle.

Here $\Omega$ is a cylinder.
~~~freefem
verbosity=0; //To minimize messages at execution

real R = 5, L = 15;

// Mesh
border cc(t=0, 2*pi) {x=cos(t)/2; y=sin(t)/2; label=2;}
border ce(t=pi/2, 3*pi/2) {x=cos(t)*R; y=sin(t)*R; label=1;}
border beb(tt=0, 1) {real t=tt^1.2; x=t*L; y=-R; label=1;}
border beu(tt=1, 0) {real t=tt^1.2; x=t*L; y=R; label=1;}
border beo(t=-R, R) {x=L; y=t; label=0;}
border bei(t=-R/4, R/4) {x=L/2; y=t; label=0;}
mesh Th = buildmesh(cc(-40) + ce(20) + beb(15) + beu(15) + beo(8) + bei(10));
plot(Th);
~~~

| The mesh   |
| ---------- |
| ![][_mesh] |

To resolve the nonlinearity we use a linearization
$$
\begin{align*}
-\nu\Delta (u+du) +u\cdot\nabla(u+du)+du\cdot\nabla v-\nabla q =0,
\\
\nabla\cdot (u+du)=0~~~in ~ \Omega
\end{align*}
$$
This, then, is the Newton method:  given $u=u^n$ compute $du,dq$ then set $u^{n+1}=u+du$
In variational form: find $du_1,du_2,dq$,
$$\begin{align*}
\int_{\Omega}\left(\nu\nabla du:\nabla v + (u\cdot\nabla du + du\cdot\nabla u +\nabla dq )v+\epsilon dq\cdot r\right))=0,~~\forall v, r
\end{align*}
$$
For stability the iterations are started with a large $\nu$ and then $\nu$ is decreased up to its desired value.
~~~freefem
real nu = 1./50, nufinal = 1/200., cnu = 0.5;
macro Grad(u1, u2) [dx(u1), dy(u1), dx(u2), dy(u2)]//
macro UgradV(u1, u2, v1, v2) [[u1, u2]'*[dx(v1), dy(v1)], [u1, u2]'*[dx(v2), dy(v2)] ]//
macro div(u1, u2)  (dx(u1) + dy(u2))//

fespace Xh(Th, P2);
Xh u1, u2, v1, v2, du1, du2, u1p, u2p;
fespace Mh(Th, P1);
Mh p, q, dp, pp;

// Intial guess with B.C.
u1 = (x^2+y^2) > 2;
u2 = 0;

// numerical parameters
real eps = 1e-4;

func bb = [[-1, -2], [4, 2]]; // bounding box for the plot

// Loop on vicosity
while(1) {
	int n;
	real err = 0;
	// Newton Loop
	for (n = 0; n < 15; n++) {
		solve Oseen([du1, du2, dp], [v1, v2, q])
			= int2d(Th)(
					nu*(Grad(du1,du2)'*Grad(v1,v2))
				+ UgradV(du1,du2, u1, u2)'*[v1,v2]
				+ UgradV( u1, u2,du1,du2)'*[v1,v2]
				- div(du1,du2)*q - div(v1,v2)*dp
				- 1e-8*dp*q // stabilization term
			)
			- int2d(Th)(
				  nu*(Grad(u1,u2)'*Grad(v1,v2))
				+ UgradV(u1,u2, u1, u2)'*[v1,v2]
				- div(u1,u2)*q - div(v1,v2)*p
				- 1e-8*p*q
			)
			+ on(1,2,du1=0,du2=0)
			;

		u1[] -= du1[];
		u2[] -= du2[];
		p[] -= dp[];

		real Lu1 = u1[].linfty, Lu2 = u2[].linfty, Lp = p[].linfty;
		err = du1[].linfty/Lu1 + du2[].linfty/Lu2 + dp[].linfty/Lp;

		cout << n << " err = " << err << " " << eps << " rey  = " << 1./nu << endl;
		if(err < eps) break; // converge
		if( n > 3 && err > 10.) break; // Blowup ?
	}
	if(err < eps) {// if converge  decrease nu (more difficult)
		plot([u1, u2], p, wait=1, cmm=" rey = " + 1./nu, coef=0.3, bb=bb);
		if(nu == nufinal) break;
		if(n < 4) cnu = cnu^1.5; // fast converge => change faster
		nu = max(nufinal, nu*cnu); // new vicosity
		u1p = u1;
		u2p = u2;
		pp = p; // save correct solution
	}
	else {  // if blowup, increase nu (more simple)
		assert(cnu < 0.95); // final blowup
		nu = nu/cnu; //  get previous value of viscosity
		cnu = cnu^(1./1.5); // no conv. => change lower
		nu = nu* cnu;  // new vicosity
		cout << " restart nu = " << nu << " Rey = " << 1./nu << "  (cnu = " << cnu << " ) \n";
		// restore correct solution
		u1 = u1p;
		u2 = u2p;
		p = pp;
	}
}
cout << " CPU "<< clock()<< " s " << endl;
~~~
| The solution   |
| -------------- |
| ![][_solution] |

[_mesh]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/NSNewton/mesh.png

[_solution]: https://raw.githubusercontent.com/FreeFem/FreeFem-markdown-figures/main/examples/examples/NSNewton/solution.png
