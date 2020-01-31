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
// Author: F. Hecht
// (jul. 2014)

verbosity=0;

// Parameters
int n = 1;
real nu = 0.0025; // Reynolds=200
real dt = 0.2;
real epsv = 1e-6, epsu = 1e-6, epsp = 1e-6;// Eps CG ..

// Mesh
border a0(t=1, 0){x=-2; y=t; label=1;} // inlet
border a1(t=-2, 0){x=t; y=0; label=2;}
border a2(t=0, -0.5){x=0; y=t; label=2;}
border a3(t=0, 1){x=18*t^1.2; y=-0.5; label=2;}
border a4(t=-0.5, 1){x=18; y=t; label=3;} // outlet
border a5(t=1, 0){x=-2+20*t; y=1; label=4;}

mesh Th = buildmesh(a0(3*n) + a1(20*n) + a2(10*n) + a3(150*n) + a4(5*n) + a5(100*n));
plot(Th);

// Fespace
fespace Vh(Th, P1);
Vh w, u = 0, v = 0, p = 0, q = 0;

// Definitions of Matrix dtMx and dtMy
matrix dtM1x, dtM1y;

macro  BuildMat()
  { /* for memory managenemt */
    varf vM(unused, v) = int2d(Th)(v) ;
    varf vdx(u, v) = int2d(Th)(v*dx(u)*dt) ;
    varf vdy(u, v) = int2d(Th)(v*dy(u)*dt) ;

    real[int] Mlump = vM(0, Vh);
    real[int] one(Vh.ndof); one = 1;
    real[int] M1 = one ./ Mlump;
    matrix dM1 = M1;
    matrix Mdx = vdx(Vh, Vh);
    matrix Mdy = vdy(Vh, Vh);
    dtM1x = dM1*Mdx;
    dtM1y = dM1*Mdy;
  }//

BuildMat

real err = 1, outflux = 1;
for(int n = 0; n < 200; n++) {
  Vh uold = u, vold = v, pold = p;

  solve pb4u(u, w, init=n, solver=CG, eps=epsu)
    = int2d(Th)(
      u*w/dt + nu*(dx(u)*dx(w) + dy(u)*dy(w))
    )
    - int2d(Th)(
      (convect([uold, vold], -dt, uold)/dt - dx(p))*w
    )
    + on(1, u=4*y*(1-y))
    + on(2, 4, u=0);
  plot(u);

  solve pb4v(v, w, init=n, solver=CG, eps=epsv)
    = int2d(Th)(
      v*w/dt + nu*(dx(v)*dx(w) + dy(v)*dy(w))
    )
    - int2d(Th)(
      (convect([uold, vold], -dt, vold)/dt - dy(p))*w
    )
    + on(1, 2, 3, 4, v=0);

  solve pb4p(q, w, solver=CG, init=n, eps=epsp)
    = int2d(Th)(dx(q)*dx(w) + dy(q)*dy(w))
    - int2d(Th)((dx(u) + dy(v))*w/dt)
    + on(3, q=0);

  // to have absolute epsilon in CG algorithm.
  epsv = -abs(epsv);
  epsu = -abs(epsu);
  epsp = -abs(epsp);

  p = pold-q;
  u[] += dtM1x*q[];
  v[] += dtM1y*q[];

  if(n%50 == 49) {
    Th = adaptmesh(Th, [u, v], q, err=0.06, nbvx=100000);
    plot(Th, wait=true);
    BuildMat // rebuild mat.
  }

  err = sqrt(int2d(Th)(square(u - uold) + square(v - vold))/Th.area);
  outflux = int1d(Th)([u, v]'*[N.x, N.y]) ;
  cout << " iter " << n << " Err L2 = " << err << " - Outflow = " << outflux << endl;
  if (err < 1e-3) break;
}
assert(abs(outflux) < 5e-3); // verification
plot(p, wait=1, ps="NSprojP.eps");
plot(u, wait=1, ps="NSprojU.eps");
