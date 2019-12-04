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
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Jean-Marie Mirebeau
// E-MAIL  : jean-marie.mirebeau@math.u-psud.fr

#ifndef EXAMPLE_METRICS_H
#define EXAMPLE_METRICS_H

#include "math.h"
#include "RZ.h"

// ********** functions ***********

/*
 * Some riemannian metrics on the unit square or cube, for the purpose of testing algorithms.
 */

template< int whichMetric >
const sym2 ExampleMetric(const R2 &P);
template< int whichMetric >
const sym3 ExampleMetric3D(const R3 &P);

/********************** 2D *****************/

template<>
const sym2 ExampleMetric< 0 >(const R2 &P) {
  return sym2(1, 0, 1);
}    // identity

template<>
const sym2 ExampleMetric< 1 >(const R2 &P) {    // A piecewise constant anisotropic metric
  const double scal = fabs(P.x - 1 / 2.) < 1 / 6. ? 4 : 1;

  return sym2(scal, -scal, 4 * scal);
}

template<>
const sym2 ExampleMetric< 2 >(const R2 &P) {    // circle, regularity Graded.
  const double delta = 0.03;                    // paramètre
  const R2 Q = P - R2(0.5, 0.5);
  const double r = Q.norm( );
  const double h = max(fabs(r - 1 / 2.), delta);

  return sym2(1 / (h * h), 1 / h, Q);
}

template<>
const sym2 ExampleMetric< 3 >(const R2 &P) {    // circle, regularity QuasiAcute.
  const double delta = 0.4;                     // paramètre
  const R2 Q = P - R2(0.5, 0.5);
  const double r = Q.norm( );
  const double h = max(fabs(r - 1 / 2.), delta);
  const double k = max(fabs(r - 1 / 2.), delta * delta);

  return sym2(1 / (k * k), 1 / (h * h), Q);
}

template<>
const sym2 ExampleMetric< 4 >(const R2 &P) {
  return sym2(10, 0, 1);
}    // diagonal

template<>
const sym2 ExampleMetric< 5 >(
  const R2 &P) {    // High anisotropy along the spiral r=k(theta+2 mu Pi), mu in {0,1,2}.
  const double pi = 4 * atan(1);
  const double width = 0.006;
  const double k = 0.4 / (6 * pi);
  const double mu = 100.;
  const R2 Q = P - R2(0.5, 0.5);
  const double r = Q.norm( );
  double theta =
    Q.x == -r ? pi : 2 * atan(Q.y / (r + Q.x));    // theta = theta >= 0 ? theta : theta+pi;

  if (fabs(r - k * theta) <= width)
    theta = theta + 0 * pi;
  else if (fabs(r - k * (theta + 2 * pi)) <= width)
    theta = theta + 2 * pi;
  else if (fabs(r - k * (theta + 4 * pi)) <= width)
    theta = theta + 4 * pi;
  else if (fabs(r - k * (theta + 6 * pi)) <= width && theta <= 0)
    theta = theta + 6 * pi;
  else
    return sym2(1, 0, 1);    // {metric[0]=1; metric[1]=0; metric[2]=1; break;}

  double c = cos(theta) - theta * sin(theta),
         s = sin(theta) + theta * cos(theta);    // tangente à la spirale
  double cOld = c;
  c = -s;
  s = cOld;

  return sym2(1, 1 / (mu * mu), R2(c, s));
}

template<>
const sym2 ExampleMetric< 6 >(const R2 &P) {    // high but constant anisotropy
  const double mu = 30., t = 0.3;
  R2 Q(cos(t), sin(t));

  return sym2(1, 1 / (mu * mu), Q);
}

template<>
const sym2 ExampleMetric< 7 >(const R2 &P) {
  const double s = 0.1 + (P - R2(0.1, 0.2)).norm( );
  return sym2( ) / square(s);
}

template<>
const sym2 ExampleMetric< 8 >(const R2 &P) {
  const double s = 0.1 + (P - R2(0.1, 0.2)).norm( );
  return sym2(100, 1, R2(1 / 2., sqrt(3.) / 2)) / square(s);
}

template<>
const sym2 ExampleMetric< 9 >(const R2 &P) {
  const double s = 0.1 + fabs(P.x);
  return sym2(100, 0, 1) / square(s);
}

/************************ 3D *************************/

template<>
const sym3 ExampleMetric3D< 0 >(const R3 &P) {
  return sym3(1, 1, 1, 0, 0, 0);
}

template<>
const sym3 ExampleMetric3D< 1 >(const R3 &P) {
  return sym3(1, 10, 100, 0, 0, 0);
}

template<>
const sym3 ExampleMetric3D< 2 >(const R3 &P) {
  return sym3(1, 10, R3(0.1, -0.2, 0.4));
}

template<>
const sym3 ExampleMetric3D< 3 >(const R3 &P) {    // tire bouchon...
  const double r0 = 0.33;
  const double theta0 = 4 * M_PI;
  const double delta0 = 0.06;
  const double mu = 1 / 8.;
  const R3 Q(P.x - 0.5, P.y - 0.5, P.z - 0.5);
  const double r = sqrt(Q.x * Q.x + Q.y * Q.y);

  if (fabs(r - r0) > delta0) {
    return sym3( );
  }

  if (square(Q.x - r * cos(theta0 * Q.z)) + square(Q.y - r * sin(theta0 * Q.z)) >
      square(r * delta0)) {
    return sym3( );
  }

  return sym3(mu * mu, 1, R3(-r0 * theta0 * sin(theta0 * Q.z), r0 * theta0 * cos(theta0 * Q.z), 1));
}

#endif
