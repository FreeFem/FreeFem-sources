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

// Parameters
real kc2 = 1; // try this value 19.4256
func g = y*(1-y);
real sigma = 20; // value of the shift
int nev = 2; // number of computed eigen value close to sigma

// Mesh
border a0(t=0, 1) { x=5; y=1+2*t;}
border a1(t=0, 1) { x=5-2*t; y=3;}
border a2(t=0, 1) { x=3-2*t; y=3-2*t;}
border a3(t=0, 1) { x=1-t; y=1;}
border a4(t=0, 1) { x=0; y=1-t;}
border a5(t=0, 1) { x=t; y=0;}
border a6(t=0, 1) { x=1+4*t; y=t;}

mesh Th=buildmesh(a0(20) + a1(20) + a2(20) + a3(20) + a4(20) + a5(20) + a6(20));

// Fespace
fespace Vh(Th, P1);
Vh u,v;
Vh u1, u2;

// Problem
solve sound(u, v)
  =int2d(Th)(u*v*kc2 - dx(u)*dx(v) - dy(u)*dy(v))
  - int1d(Th, a4)(g*v);

plot(u, wait=1, ps="sound0.eps");

// Eigen values
varf op(u1, u2)
  = int2d(Th)(dx(u1)*dx(u2) + dy(u1)*dy(u2) - sigma* u1*u2);

varf b([u1], [u2])
  = int2d(Th)(u1*u2); // no boundary condition

matrix OP = op(Vh, Vh, solver=Crout);  // crout solver because the matrix in not positive
matrix B = b(Vh, Vh, solver=CG, eps=1e-20);

real[int] ev(nev); // to store the nev eigenvalue
Vh[int] eV(nev); // to store the nev eigenvector

int k = EigenValue(OP, B, sym=true, sigma=sigma, value=ev, vector=eV,
  tol=1e-10, maxit=0, ncv=0);
cout << ev(0) << " 2 eigen values " << ev(1) << endl;
v = eV[0];
plot(v, wait=1, ps="eigen.eps");
