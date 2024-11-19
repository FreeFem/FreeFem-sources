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

// Mesh
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
z[0] = 2; z[1] = 3; z[2] = 4;

// Solve
A; ud = u;

ofstream f("J.txt");
func real J(real[int] & Z) {
    for (int i = 0;i < z.n; i++) z[i] = Z[i];
    A;
    real s = int2d(th)(Ie*(u-ud)^2);
    f << s << "   ";
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

real[int] Z(3);
for(int j = 0; j < z.n; j++) Z[j] = 1;

BFGS(J, DJ, Z, eps=1.e-6, nbiter=15, nbiterline=20);
cout << "BFGS: J(z) = " << J(Z) << endl;
for(int j = 0; j < z.n; j++) cout << z[j] << endl;
plot(ud, value=1, ps="u.eps");
