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
int m = 30;
int L = 80;
int LL = 80;
int j = 100;
real sigmax = 0.3;
real sigmay = 0.3;
real rho = 0.3;
real r = 0.05;
real K = 40;
real dt = 0.01;

// Mesh
mesh th = square(m, m, [L*x, LL*y]);

// Fespace
fespace Vh(th, P1);
Vh u = max(K - max(x, y), 0.);
Vh xveloc, yveloc, v, uold;

// Time loop
for (int n = 0; n*dt <= 1.0; n++) {
    // Mesh adaptation
    if (j > 20) {
        th = adaptmesh(th, u, verbosity=1, abserror=1, nbjacoby=2,
            err=0.001, nbvx=5000, omega=1.8, ratio=1.8, nbsmooth=3,
            splitpbedge=1, maxsubdiv=5, rescaling=1) ;
        j = 0;
        xveloc = -x*r + x*sigmax^2 + x*rho*sigmax*sigmay/2;
        yveloc = -y*r + y*sigmay^2 + y*rho*sigmax*sigmay/2;
        u = u;
    }

    // Update
    uold = u;

    // Solve
    solve eq1(u, v, init=j, solver=LU)
        = int2d(th)(
              u*v*(r + 1/dt)
            + dx(u)*dx(v)*(x*sigmax)^2/2
            + dy(u)*dy(v)*(y*sigmay)^2/2
            + dy(u)*dx(v)*rho*sigmax*sigmay*x*y/2
            + dx(u)*dy(v)*rho*sigmax*sigmay*x*y/2
        )
        + int2d(th)(
            - v*convect([xveloc, yveloc], dt, uold)/dt
        )
        + on(2,3,u=0)
        ;

    j = j + 1;
}

// Plot
plot(u, wait=true, value=true);
plot(th, wait=true);
