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
real theta = 4.*pi/3.;
real a = 2., b = 1.; // the length of the semimajor axis and semiminor axis
func z = x;

// Mesh
border Gamma1(t=0, theta)    {x=a*cos(t); y=b*sin(t);}
border Gamma2(t=theta, 2*pi) {x=a*cos(t); y=b*sin(t);}
mesh Th = buildmesh(Gamma1(100) + Gamma2(50));

// Fespace
fespace Vh(Th,P2); // P2 conforming triangular finite element space
Vh phi, w, f = 1;

// Problem (resolution of laplace equation)
solve Laplace(phi, w)
  =int2d(Th)(
    dx(phi)*dx(w) + dy(phi)*dy(w)
  )
  - int2d(Th)(f*w)
  + on(Gamma1, phi=z);

// Plot
plot(phi, wait=true, ps="membrane.eps"); //Plot Th and phi
plot(Th, wait=true, ps="membraneTh.eps"); //Plot Th

// Export to gnupot
{
  ofstream ff("graph.txt");
  for (int i = 0; i < Th.nt; i++) {
    for (int j = 0; j < 3; j++)
      ff << Th[i][j].x  << "    " << Th[i][j].y << "  " << phi[][Vh(i,j)] << endl;
  ff << Th[i][0].x << "    " << Th[i][0].y << "  " << phi[][Vh(i,0)] << endl << endl << endl;
  }
}

// Save the mesh
savemesh(Th, "Th.msh");
