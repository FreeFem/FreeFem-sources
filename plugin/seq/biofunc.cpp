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

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

#include <ff++.hpp>
#include <AFunction_ext.hpp>
#include <cstdlib>
using namespace Fem2D;
typedef double R;
double fmonod(const double &xx, const double &k, const double &kb) {
  double x = max(xx, 0.);
  return (k * x) / (kb + x);
}
// on array ..
KN<double> * fmonod( KN<double> * const&  pf , KN<double> *const &  px, const double &k, const double &kb)
{
    KN<double> &f = *pf, &xx=*px;
    long n = f.N();
    ffassert(n == xx.N());
    for( int i=0; i<n; ++i)
    {
    double x = max(xx[i], 0.);
    f[i]=  (k * x) / (kb + x);
    }
    return pf;
}
KN<double> * fmonod( KN<double> * const&  pf , KN<double> *const &  px, KN<double> *const &pk, const double &kb)
{
    KN<double> &f = *pf, &xx=*px , & k= *pk;
    long n = f.N();
    ffassert(n == xx.N() && n == k.N());
    for( int i=0; i<n; ++i)
    {
        double x = max(xx[i], 0.);
        f[i]=  (k[i] * x) / (kb + x);
    }
    return pf;
}
KN<double> * dfmonod( KN<double> * const&  pf , KN<double> *const &  px, const double &k, const double &kb)
{
    KN<double> &f = *pf, &xx=*px;
    long n = f.N();
    ffassert(n == xx.N());
    for( int i=0; i<n; ++i)
    {
        double x = max(xx[i], 0.), a = kb + x;
        f[i]=  k / a - k * x / (a * a);
    }
    return pf;
}
KN<double> * dfmonod( KN<double> * const&  pf , KN<double> *const &  px,KN<double> *const &pk, const double &kb)
{
    KN<double> &f = *pf, &xx=*px , & k= *pk;
    long n = f.N();
    ffassert(n == xx.N() && n == k.N());    for( int i=0; i<n; ++i)
    {
        double x = max(xx[i], 0.), a = kb + x;
        f[i]=  k[i] / a - k[i] * x / (a * a);
    }
    return pf;
}
double dfmonod(const double &xx, const double &k, const double &kb) {
  double x = max(xx, 0.), a = kb + x;
  return k / a - k * x / (a * a);
}
double fmonod(const double &xx, const double &kb) {
  double x = max(xx, 0.);
  return (x) / (kb + x);
}
double dfmonod(const double &xx, const double &kb) {
  double x = max(xx, 0.), a = kb + x;
  return 1. / a - x / (a * a);
}

static void init( ) {
  Global.Add("fmonod", "(", new OneOperator4_< KN<double> *, KN<double> *, KN<double> *, R, R >(fmonod));
  Global.Add("fmonod", "(", new OneOperator4_< KN<double> *, KN<double> *, KN<double> *, KN<double> *, R >(fmonod));
  Global.Add("dfmonod", "(", new OneOperator4_< KN<double> *, KN<double> *, KN<double> *, R, R >(dfmonod));
  Global.Add("dfmonod", "(", new OneOperator4_< KN<double> *, KN<double> *, KN<double> *, KN<double> *, R >(dfmonod));
  Global.Add("fmonod", "(", new OneOperator3_< R, R, R, R >(fmonod));
  Global.Add("dfmonod", "(", new OneOperator3_< R, R, R, R >(dfmonod));
  Global.Add("fmonod", "(", new OneOperator2_< R, R, R >(fmonod));
  Global.Add("dfmonod", "(", new OneOperator2_< R, R, R >(dfmonod));
}

LOADFUNC(init);
