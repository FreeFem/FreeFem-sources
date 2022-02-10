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

// Example C++ function "myfunction", dynamically loaded into "ff-c++ dfft.cpp "

/* clang-format off */
//ff-c++-LIBRARY-dep: fftw3
//ff-c++-cpp-dep:
/* clang-format on */

#include "ff++.hpp"
#include "AFunction_ext.hpp"
#include <fftw3.h>

template< class Complex >
class DFFT_1d2dor3d {
 public:
  Complex *x;
  int n, m, k;
  int sign;
  DFFT_1d2dor3d(KN< Complex > *xx, long signn, long nn = 1, long kk = 1)
    : x(*xx), n(nn), m(xx->N( ) / (nn * kk)), k(kk), sign(signn) {
    cout << xx << " " << signn << " " << nn << " " << xx->N( ) << " n: " << n << " m:" << m
         << " k:  " << k << endl;
    ffassert(n > 0 && (n * m * k == xx->N( )));
  }

  DFFT_1d2dor3d(KNM< Complex > *xx, long signn) : x(*xx), n(xx->M( )), m(xx->N( )), sign(signn) {}
};

DFFT_1d2dor3d< Complex > dfft(KN< Complex > *const &x, const long &sign) {
  return DFFT_1d2dor3d< Complex >(x, sign);
}

DFFT_1d2dor3d< Complex > dfft(KN< Complex > *const &x, const long &nn, const long &sign) {
  return DFFT_1d2dor3d< Complex >(x, sign, nn);
}

DFFT_1d2dor3d< Complex > dfft(KN< Complex > *const &x, const long &nn, const long &kk,
                              const long &sign) {
  return DFFT_1d2dor3d< Complex >(x, sign, nn, kk);
}

DFFT_1d2dor3d< Complex > dfft(KNM< Complex > *const &x, const long &sign) {
  return DFFT_1d2dor3d< Complex >(x, sign);
}

bool ff_execute(fftw_plan *p) {
  if (*p) {
    fftw_execute(*p);
  }

  return 0;
}

bool ff_delete(fftw_plan *p) {
  if (*p) {
    fftw_destroy_plan(*p);
  }

  *p = 0;
  return 0;
}

KN< Complex > *dfft_eq(KN< Complex > *const &x, const DFFT_1d2dor3d< Complex > &d) {
  ffassert(x->N( ) == d.n * d.m * d.k);
  Complex *px = *x;
  fftw_plan p;
  if (d.k == 1) {
    if (d.n > 1) {
      p = fftw_plan_dft_2d(d.n, d.m, reinterpret_cast< fftw_complex * >(d.x),
                           reinterpret_cast< fftw_complex * >(px), d.sign, FFTW_ESTIMATE);
    } else {
      p = fftw_plan_dft_1d(d.m, reinterpret_cast< fftw_complex * >(d.x),
                           reinterpret_cast< fftw_complex * >(px), d.sign, FFTW_ESTIMATE);
    }
  } else {
    if (d.n > 1) {
      p = fftw_plan_dft_3d(d.n, d.m, d.k, reinterpret_cast< fftw_complex * >(d.x),
                           reinterpret_cast< fftw_complex * >(px), d.sign, FFTW_ESTIMATE);
    } else {
      p = fftw_plan_dft_2d(d.m, d.k, reinterpret_cast< fftw_complex * >(d.x),
                           reinterpret_cast< fftw_complex * >(px), d.sign, FFTW_ESTIMATE);
    }
  }

  fftw_execute(p);
  fftw_destroy_plan(p);
  return x;
}

KN< double > *dfft_eq(KN< double > *const &x, const DFFT_1d2dor3d< double > &d) {
  ffassert(0);
  return x;
}

/*  class Init { public:
 * Init();
 * };*/
// bofbof ..
struct fftw_plan_s {};

// ...

template<>
inline AnyType DeletePtr< fftw_plan * >(Stack, const AnyType &x) {
  fftw_plan *a = PGetAny< fftw_plan >(x);

  if (*a) {
    fftw_destroy_plan(*a);
  }

  *a = 0;
  return Nothing;
};

fftw_plan *plan__eq(fftw_plan *a, fftw_plan b) {
  if (*a) {
    fftw_destroy_plan(*a);
  }

  *a = b;
  return a;
}

fftw_plan *plan_set(fftw_plan *a, fftw_plan b) {
  *a = b;
  return a;
}

fftw_plan plan_dfft(KN< Complex > *const &x, KN< Complex > *const &y, const long &sign) {
  return fftw_plan_dft_1d(x->N( ), reinterpret_cast< fftw_complex * >(&x[0]),
                          reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
}

fftw_plan plan_dfft(KNM< Complex > *const &x, KNM< Complex > *const &y, const long &sign) {
  long m = x->N( ), n = x->M( );

  fftw_plan_dft_2d(n, m, reinterpret_cast< fftw_complex * >(&x[0]),
                   reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
  return 0;
}

fftw_plan plan_dfft(KN< Complex > *const &x, KN< Complex > *const &y, const long &n,
                    const long &sign) {
  long nn = n, mm = y->N( ) / nn;

  ffassert(mm * nn == y->N( ) && x->N( ) == y->N( ));

  return fftw_plan_dft_2d(nn, mm, reinterpret_cast< fftw_complex * >(&x[0]),
                          reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
}

fftw_plan plan_dfft(KN< Complex > *const &x, KN< Complex > *const &y, const long &n, const long &k,
                    const long &sign) {
  int nn = n, mm = y->N( ) / (k * n), kk = k;

  ffassert(y->N( ) == nn * mm * kk);
  if (nn > 1) {
    return fftw_plan_dft_3d(nn, mm, kk, reinterpret_cast< fftw_complex * >(&x[0]),
                            reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
  } else {
    return fftw_plan_dft_2d(nn, mm, reinterpret_cast< fftw_complex * >(&x[0]),
                            reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
  }
}

class Mapkk : public E_F0mps {
 public:
  typedef Complex R;
  typedef KN_< R > Result;
  ;
  static basicAC_F0::name_and_type *name_param;
  static const int n_name_param = 0;
  Expression expv, expm, exp;
  Expression nargs[n_name_param];

  Mapkk(const basicAC_F0 &args) : expv(0), expm(0), exp(0) {
    args.SetNameParam(n_name_param, name_param, nargs);
    expv = to< KN< R > * >(args[0]);    // a the expression to get the mesh
    expm = to< long >(args[1]);
    exp = to< R >(args[2]);    // a the expression to get the mesh
  }

  ~Mapkk( ) {}

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< KN< R > * >( ), atype< long >( ), atype< R >( ));
  }

  static E_F0 *f(const basicAC_F0 &args) { return new Mapkk(args); }

  AnyType operator( )(Stack s) const;
};

basicAC_F0::name_and_type *Mapkk::name_param = 0;
AnyType Mapkk::operator( )(Stack s) const {
  // correct July 2015 ... not tested before..
  MeshPoint *mp(MeshPointStack(s)), mps = *mp;

  KN< R > *pv = GetAny< KN< R > * >((*expv)(s));
  KN< R > v(*pv);

  long nn = v.N( );
  long m = GetAny< long >((*expm)(s));
  if (verbosity > 10) {
    cout << "  map: expm " << expm << " m = " << m << endl;
  }

  long n = nn / m;
  double ki = 1. / n;
  double kj = 1. / m;
  double ki0 = 0., kj0 = 0;
  if (verbosity > 10) {
    cout << " map: " << n << " " << m << " " << nn << " == " << n * m << endl;
  }

  ffassert(m * n == nn);
  long n2 = (n + 1) / 2, m2 = (m + 1) / 2;

  for (long j = 0, kk = 0; j < m; ++j) {
    for (long k = 0, i = 0; i < n; ++i) {
      R2 P(i * ki + ki0, j * kj + kj0);
      mp->set(P.x, P.y);
      v[kk++] = GetAny< R >((*exp)(s));
    }
  }

  *mp = mps;
  return 0L;
}

static void Load_Init( ) {
  typedef DFFT_1d2dor3d< Complex > DFFT_C;
  typedef DFFT_1d2dor3d< double > DFFT_R;

  cout << " load: init dfft " << endl;
  Dcl_Type< DFFT_C >( );
  Dcl_Type< DFFT_R >( );

  Dcl_Type< fftw_plan * >(::InitializePtr< fftw_plan * >, ::DeletePtr< fftw_plan * >);
  Dcl_Type< fftw_plan >( );
  zzzfff->Add("fftwplan", atype< fftw_plan * >( ));

  TheOperators->Add("=", new OneOperator2< fftw_plan *, fftw_plan *, fftw_plan >(plan__eq));
  TheOperators->Add("<-", new OneOperator2< fftw_plan *, fftw_plan *, fftw_plan >(plan_set));

  Global.Add("plandfft", "(",
             new OneOperator3_< fftw_plan, KN< Complex > *, KN< Complex > *, long >(plan_dfft));
  Global.Add(
    "plandfft", "(",
    new OneOperator4_< fftw_plan, KN< Complex > *, KN< Complex > *, long, long >(plan_dfft));
  Global.Add(
    "plandfft", "(",
    new OneOperator5_< fftw_plan, KN< Complex > *, KN< Complex > *, long, long, long >(plan_dfft));
  Global.Add("plandfft", "(",
             new OneOperator3_< fftw_plan, KNM< Complex > *, KNM< Complex > *, long >(plan_dfft));

  Global.Add("execute", "(", new OneOperator1< bool, fftw_plan * >(ff_execute));
  Global.Add("delete", "(", new OneOperator1< bool, fftw_plan * >(ff_delete));

  Global.Add("dfft", "(", new OneOperator2_< DFFT_C, KN< Complex > *, long >(dfft));
  Global.Add("dfft", "(", new OneOperator3_< DFFT_C, KN< Complex > *, long, long >(dfft));
  Global.Add("dfft", "(", new OneOperator4_< DFFT_C, KN< Complex > *, long, long, long >(dfft));
  Global.Add("dfft", "(", new OneOperator2_< DFFT_C, KNM< Complex > *, long >(dfft));
  Global.Add("map", "(", new OneOperatorCode< Mapkk >( ));
  TheOperators->Add("=", new OneOperator2_< KN< Complex > *, KN< Complex > *, DFFT_C >(dfft_eq));
}

LOADFUNC(Load_Init)
