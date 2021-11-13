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
// Corrected by F. Hecht may 2021 
// Parameters
real S = 99;

border C(t=0, 2*pi){x=3*cos(t); y=3*sin(t);} // Label 1,2 
border Splus(t=0, 1){x=t-0.5; y=0.17735*sqrt(t) - 0.075597*t - 0.212836*(t^2) + 0.17363*(t^3) - 0.06254*(t^4); label=S;}
border Sminus(t=1, 0){x=t-0.5; y=-(0.17735*sqrt(t) - 0.075597*t - 0.212836*(t^2) + 0.17363*(t^3) - 0.06254*(t^4)); label=S;}
mesh Th = buildmesh(C(50) + Splus(70) + Sminus(70));
// Fespace
fespace Vh(Th, P2);
Vh psi, w;
real cost = cos(5.*pi/180.), sint=sin(5.*pi/180.);// incidence 5 degres
// Problem
solve potential(psi, w)
  = int2d(Th)(dx(psi)*dx(w)+dy(psi)*dy(w))
  + on(C, psi = cost*y-sint*x) 
  + on(S, psi=0);

// Plot
plot(psi, wait=1);

/// Thermic
// Parameters
real dt = 0.005, nbT = 50;

// Mesh
border D(t=0, 2.){x=0.5+t*cost; y=+t*sint;}
mesh Sh = buildmesh(C(25) + Splus(-90) + Sminus(-90) + D(200));
int steel = Sh(0.5, 0).region, air = Sh(-1, 0).region;
// Change label to put BC on In flow 
// Fespace
fespace Wh(Sh, P1);
Wh  vv;

fespace W0(Sh, P0);
W0 k = 0.01*(region == air) + 0.1*(region == steel);
W0 u1 = dy(psi)*(region == air), u2 = -dx(psi)*(region == air);
Wh v = 120*(region == steel), vold;
// pul label 10 on inflow boundary to inforce the temperature.
Sh = change(Sh,flabel = (label == C &&  [u1,u2]'*N<0) ? 10 : label);
int i;
problem thermic(v, vv, init=i, solver=LU)
  = int2d(Sh)(
      v*vv/dt + k*(dx(v)*dx(vv) + dy(v)*dy(vv))
    + 10*(u1*dx(v) + u2*dy(v))*vv
  )
  - int2d(Sh)(vold*vv/dt)
  + on(10, v= 0);
  

for(i = 0; i < nbT; i++) {
    vold[]= v[];
    thermic;
    plot(v);
}
plot(v, wait=1,fill=1,value=1);
