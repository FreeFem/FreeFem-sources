/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/

/// Characteristics Galerkin
// Parameters
verbosity = 1;
real dt = 0.17;
real t = 0;

// Mesh
border C(t=0, 2*pi){x=cos(t); y=sin(t);}
mesh Th = buildmesh(C(100));

// Fespace
fespace Uh(Th, P1);
Uh cold, c = exp(-10*((x-0.3)^2 +(y-0.3)^2));
Uh u1 = y, u2 = -x;

// Loop
for (int m = 0; m < 2*pi/dt; m++) {
	t += dt;
	cold = c;
	c = convect([u1, u2], -dt, cold);
	plot(c, cmm=" t="+t + ", min=" + c[].min + ", max=" +  c[].max);
}

/// Now with Discontinuous Galerkin
// Parameters
real u, al=0.5;
dt = 0.05;

// Fespace
fespace Vh(Th, P1dc);
Vh w, ccold, v1 = y, v2 = -x, cc = exp(-10*((x-0.3)^2 +(y-0.3)^2));

// Macro
macro n()(N.x*v1+N.y*v2) //

// Problem
problem  Adual(cc, w, init=t)
	= int2d(Th)((cc/dt + (v1*dx(cc) + v2*dy(cc)))*w)
	+ intalledges(Th)((1-nTonEdge)*w*(al*abs(n) - n/2)*jump(cc))
//  - int1d(Th, C)((n(u)<0)*abs(n(u))*cc*w)	// unused because cc=0 on d(Omega)^-
	- int2d(Th)(ccold*w/dt);

// Loop
for (t = 0; t < 2*pi; t += dt) {
	ccold = cc;
	Adual;
	plot(cc, fill=1, cmm="t="+t + ", min=" + cc[].min + ", max=" +  cc[].max);
}

// Plot
real [int] viso=[-0.1, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1];
plot(c, wait=1, fill=1, value=1, ps="convectCG.eps", viso=viso);
plot(cc, wait=1, fill=1, value=1, ps="convectDG.eps", viso=viso);

/// Discontinuous Galerkin very much faster
// Problem
varf aadual(cc, w)
	= int2d(Th)((cc/dt + (v1*dx(cc) + v2*dy(cc)))*w)
	+ intalledges(Th)((1-nTonEdge)*w*(al*abs(n) - n/2)*jump(cc));

varf bbdual(ccold, w) = -int2d(Th)(ccold*w/dt);

matrix  AA = aadual(Vh, Vh, verb=1);
matrix BB = bbdual(Vh, Vh);

// Loop
Vh rhs = 0;
for (t = 0; t < 2*pi ; t += dt) {
	ccold = cc;
	rhs[] = BB* ccold[];
	cc[] = AA^-1*rhs[];
	plot(cc, fill=0, cmm="t="+t + ", min=" + cc[].min + ", max=" +  cc[].max);
}

// Plot
plot(cc, wait=1, fill=1, value=1, ps="convectDG.eps", viso=viso);
