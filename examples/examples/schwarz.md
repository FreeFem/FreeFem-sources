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
real n = 4; // Mesh quality

// Mesh
border a(t=0, 1){x=t; y=0;}
border a1(t=1, 2){x=t; y=0;}
border b(t=0, 1){x=2; y=t;}
border c(t=2, 0){x=t; y=1;}
border d(t=1, 0){x=0; y=t;}
border e(t=0, pi/2){x=cos(t); y=sin(t);}
border e1(t=pi/2, 2*pi){x=cos(t); y=sin(t);}

//Omega1 (rectangle)
mesh th = buildmesh(a(5*n) + a1(5*n) + b(5*n) + c(10*n) + d(5*n));
fespace Vh(th, P1);
Vh v, u=0;

//Omega2 (circle)
mesh TH = buildmesh(e(5*n) + e1(25*n));
fespace VH(TH,P1);
VH V, U=0;

for (int i = 0; i < 4; i++) {
  plot(U, u, wait=1, cmm="Iteration "+i);
  // Solve on Omega2
  solve AA(U, V)
    = int2d(TH)(dx(U)*dx(V) + dy(U)*dy(V))
    - int2d(TH)(V)
    + on(e, U=u)
    + on(e1, U=0);

  // Solve on Omega1
  solve aa(u, v)
    = int2d(th)(dx(u)*dx(v) + dy(u)*dy(v))
    - int2d(th)(v)
    + on(a, d, u=U)
    + on(a1, b, c, u=0);
}

plot(U, u, wait=1, cmm="Final solution");
