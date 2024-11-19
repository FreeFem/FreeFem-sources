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
int n = 3; // mesh quality

// Mesh
mesh Th = square(10*n, 10*n);

// Fespace
fespace Uh(Th, P1b);
Uh u, v;
Uh uu, vv;
fespace Ph(Th, P1);
Ph p, pp;

// Problem
solve stokes([u, v, p], [uu, vv, pp])
  = int2d(Th)(
      dx(u)*dx(uu) + dy(u)*dy(uu)
    + dx(v)*dx(vv) + dy(v)*dy(vv)
    + dx(p)*uu + dy(p)*vv
    + pp*(dx(u) + dy(v))
    -1e-10*p*pp
  )
  + on(1, 2, 4, u=0, v=0)
  + on(3, u=1, v=0)
  ;

// Plot
plot([u,v],p,wait=1);
