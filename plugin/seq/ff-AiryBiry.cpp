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
// AUTHORS : ...
// E-MAIL  : ...

#include <iostream>
#include <math.h>
using namespace std;

double airy(double x, long df) {
  double f, y, a, b, s;
  int p;
  double u = .258819403792807, v = .355028053887817;

  if (x <= 1.7 && x >= -6.9) {
    y = x * x * x / 9.;
    if (df) {
      a = 2. / 3.;
      b = -a;
      v *= x * x / 2.;
      u = -u;
    } else {
      b = 1. / 3.;
      a = -b;
      u *= -x;
    }

    for (p = 1, f = u + v;; ++p) {
      a += 1.;
      b += 1.;
      v *= y / (p * a);
      u *= y / (p * b);
      s = u + v;
      f += s;
      if (fabs(s) < 1.e-14) {
        break;
      }
    }
  } else {
    s = 1. / sqrt(v = 3.14159265358979);
    y = fabs(x);
    if (df) {
      s *= pow(y, .25);
    } else {
      s /= pow(y, .25);
    }

    y *= 2. * sqrt(y) / 3.;
    if (x > 0.) {
      a = 12. / pow(y, .333);
      p = a * a;
      if (df) {
        a = -7. / 36.;
      } else {
        a = 5. / 36.;
      }

      b = 2. * (p + y);
      f = 1.;
      u = x = 0.;
      s *= exp(-y) / 2.;

      for (; p > 0; --p, b -= 2.) {
        y = (b * f - (p + 1) * x) / (p - 1 + a / p);
        x = f;
        f = y;
        u += f;
      }

      if (df) {
        f *= -s / u;
      } else {
        f *= s / u;
      }
    } else {
      x = y - v / 4.;
      y *= 2.;
      b = .5;
      f = s;
      v = 0.;
      if (df) {
        a = 2. / 3.;
      } else {
        a = 1. / 3.;
      }

      for (p = 1; (u = fabs(s)) > 1.e-14; ++p, b += 1.) {
        s *= (a + b) * (a - b) / (p * y);
        if (fabs(s) >= u) {
          break;
        }

        if (!(p & 1)) {
          s = -s;
          f += s;
        } else {
          v += s;
        }
      }

      if (df) {
        f = f * sin(x) + v * cos(x);
      } else {
        f = f * cos(x) - v * sin(x);
      }
    }
  }

  return f;
}

double biry(double x, long df) {
  double f, y, a, b, s;
  int p;
  double u = .258819403792807, v = .355028053887817;

  if (x <= 7.6 && x >= -6.9) {
    y = x * x * x / 9.;
    if (df) {
      b = -(a = 2. / 3.);
      u *= (f = sqrt(3.));
      v *= f * x * x / 2.;
    } else {
      a = -(b = 1. / 3.);
      v *= (f = sqrt(3.));
      u *= f * x;
    }

    for (p = 1, f = u + v;; ++p) {
      v *= y / (p * (a += 1.));
      u *= y / (p * (b += 1.));
      f += (s = u + v);
      if (fabs(s) < 1.e-14 * (1. + fabs(f))) {
        break;
      }
    }
  } else {
    s = 1. / sqrt(v = 3.14159265358979);
    y = fabs(x);
    if (df) {
      s *= pow(y, .25);
    } else {
      s /= pow(y, .25);
    }

    y *= 2. * sqrt(y) / 3.;
    b = .5;
    if (df) {
      a = 2. / 3.;
    } else {
      a = 1. / 3.;
    }

    if (x > 0.) {
      s *= exp(y);
      f = s;
      y *= -2.;

      for (p = 1; (u = fabs(s)) > 1.e-14; ++p, b += 1.) {
        s *= (a + b) * (a - b) / (p * y);
        if (fabs(s) >= u) {
          break;
        }

        f += s;
      }
    } else {
      x = y - v / 4.;
      y *= 2.;
      f = s;
      v = 0.;

      for (p = 1; (u = fabs(s)) > 1.e-14; ++p, b += 1.) {
        s *= (a + b) * (a - b) / (p * y);
        if (fabs(s) >= u) {
          break;
        }

        if (!(p & 1)) {
          s = -s;
          f += s;
        } else {
          v += s;
        }
      }

      if (df) {
        f = f * cos(x) - v * sin(x);
      } else {
        f = -(f * sin(x) + v * cos(x));
      }
    }
  }

  return f;
}

#include "ff++.hpp"
static void InitFF( ) {
  Global.Add("airy", "(", new OneOperator2< double, double, long >(airy));
  Global.Add("biry", "(", new OneOperator2< double, double, long >(biry));
}

LOADFUNC(InitFF)
