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
func u0 = 10+90*x/6;
func k = 1.8*(y<0.5) + 0.2;
real ue = 25, alpha = 0.25, T = 5, dt = 0.1;
real rad = 1e-8, uek = ue + 273.;

// Mesh
mesh Th = square(30, 5, [6*x, y]);

// Fespace
fespace Vh(Th,P1);
Vh vold, w, v=u0-ue, b;

// Problem
problem thermradia(v, w)
  = int2d(Th)(
      v*w/dt
    + k*(dx(v)*dx(w) + dy(v)*dy(w))
  )
  + int1d(Th, 1, 3)(b*v*w)
  - int2d(Th)(vold*w/dt)
  + on(2, 4, v=u0-ue);

for(real t = 0; t < T; t+=dt) {
  vold = v;
  for (int m = 0; m < 5; m++) {
    b = alpha + rad*(v + 2*uek)*((v+uek)^2 + uek^2);
    thermradia;
  }
}
vold = v + ue;

// Plot
plot(vold);
