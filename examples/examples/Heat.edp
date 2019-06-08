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
bool withplot = 0;
real cpu1, cpu2;
func k = 1.8*(y < 0.5) + 0.2;
real kf = 1, ue = 20, T = 5, dt = 0.1;

// Mesh
mesh Th = square(150, 50, [3*x, y]);

// Fespace
fespace Vh(Th, P1);
Vh u, uold, v, usave;

// Problem
int kk=0;
problem Heat(u, v, init=kk)
  = int2d(Th)(
      u*v/dt
    + k*(dx(u)*dx(v) + dy(u)*dy(v))
  )
  + int1d(Th, 1, 3)(kf*v*u)
  - int1d(Th, 1, 3)(kf*v*ue)
  - int2d(Th)(uold*v/dt)
  + on(2, 4, u=30);

real cpu = clock();

// Initialization
u = 0;

// Basic time loop
for (real t = 0; t < T; t += dt) {
  uold = u;
  Heat;
  kk++;
  if (withplot) plot(u);
}

cpu1 = clock()-cpu;

plot(u);
usave[] = u[];

/// Optimized code
// Problem definition
varf vA(u, v)
  = int2d(Th)(
      u*v/dt
    + k*(dx(u)*dx(v) + dy(u)*dy(v))
    )
  + int1d(Th, 1, 3)(kf*v*u)
  + on(2, 4, u=30);
varf vB(u, v) = int2d(Th)(u*v/dt) ;
varf vRHS(u, v) = int1d(Th, 1, 3)(kf*v*ue);
varf vL(u, v) = on(2, 4, u=30);

{
  cpu = clock();
  real[int] rhsbc = vL(0, Vh);
  real[int] rhs0 = vRHS(0, Vh);
  matrix A = vA(Vh, Vh, solver=sparsesolver);
  matrix B = vB(Vh, Vh);

  // Initialization
  u = 0;

  // Optimized time loop (speed of C language)
  for (real t = 0; t < T; t += dt) {
    real[int] b = B*u[];
    b += rhs0;
    b = rhsbc ? rhsbc : b;
    u[] = A^-1*b;
    if(withplot) plot(u);
  }
  cpu2 = clock() - cpu;
}

plot(u, cmm="u2");
cout << " cpu method 1 = " << cpu1 << " cpu method matrix = " << cpu2 << " ratio = " <<  cpu1/cpu2 << endl;
