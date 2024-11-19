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
int C1 = 99, C2 = 98; // could be anything

// Mesh
border C0(t=0, 2*pi){x=5*cos(t); y=5*sin(t);}

border C11(t=0, 1){x=1+t;  y=3;      label=C1;}
border C12(t=0, 1){x=2;    y=3-6*t;  label=C1;}
border C13(t=0, 1){x=2-t;  y=-3;     label=C1;}
border C14(t=0, 1){x=1;    y=-3+6*t; label=C1;}

border C21(t=0, 1){x=-2+t; y=3;      label=C2;}
border C22(t=0, 1){x=-1;   y=3-6*t;  label=C2;}
border C23(t=0, 1){x=-1-t; y=-3;     label=C2;}
border C24(t=0, 1){x=-2;   y=-3+6*t; label=C2;}

mesh Th=buildmesh(C0(50)
  + C11(5) + C12(20) + C13(5) + C14(20)
  + C21(-5) + C22(-20) + C23(-5) + C24(-20));
plot(Th, wait=1, ps="heatexTh.ps");

// Fespace
fespace Vh(Th, P1);
Vh u, v;
Vh kappa = 1 + 4*(x<-1)*(x>-2)*(y<3)*(y>-3);

// Problem
solve a(u, v)
  = int2d(Th)(
    kappa*(dx(u)*dx(v) + dy(u)*dy(v))
  )
  + on(C0, u=20)
  +on(C1, u=100);

// Plot
plot(u, value=true, wait=1, fill=true);
