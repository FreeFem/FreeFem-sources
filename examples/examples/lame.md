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
real E = 21e5, nu = 0.28;
real f = -1;

// Mesh
mesh Th = square(10, 10, [20*x, 2*y-1]);

// Fespace
fespace Vh(Th, P2);
Vh u, v, uu, vv;

// Macro
real sqrt2 = sqrt(2.);
macro epsilon(u1, u2) [dx(u1), dy(u2), (dy(u1)+dx(u2))/sqrt2] // EOM
// Remark: the 1/sqrt2 term in (dy(u1)+dx(u2)) is here
// to obtain a 1/2 when we use epsilon(u1, u2)'*epsilon(v1, v2)
macro div(u, v) (dx(u) + dy(v)) // EOM

// Problem
real mu = E/(2*(1 + nu));
real lambda = E*nu/((1 + nu)*(1 - 2*nu));

solve lame([u, v], [uu, vv])
  = int2d(Th)(
      lambda*div(u,v)*div(uu,vv)
    + 2.*mu*(epsilon(u,v)'*epsilon(uu,vv))
  )
  - int2d(Th)(f*vv)
  + on(4, u=0, v=0);

// Plot
real coef = 100;
plot([u, v], wait=1, ps="lamevect.eps", coef=coef);

// Move mesh
mesh th1 = movemesh(Th, [x+u*coef, y+v*coef]);
plot(th1, wait=1, ps="lamedeform.eps");

// Display
real dxmin = u[].min;
real dymin = v[].min;
cout << " - dep.  max x = " << dxmin << " y = " << dymin << endl;
cout << "   dep. (20,0) = " << u(20,0) << " " << v(20,0) << endl;
