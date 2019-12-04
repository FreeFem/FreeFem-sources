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

// Example C++ function "myfunction", dynamically loaded into "load.edp"

/* clang-format off */
//ff-c++-LIBRARY-dep: gsl
//ff-c++-cpp-dep:
/* clang-format on */

#include <ff++.hpp>
#include <AFunction_ext.hpp>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_clausen.h>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_dawson.h>
#include <gsl/gsl_sf_debye.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_gegenbauer.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_mathieu.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_transport.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_poly.h>

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_spline.h>

#include <gsl/gsl_multifit.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "ff_gsl_awk.hpp"
#include "gsl/gsl_version.h"
struct GSLInterpolation {
  gsl_interp_accel *acc;
  gsl_spline *spline;
  double *xy;
  size_t n;
  const gsl_interp_type *splinetype;
  void init( ) {
    spline = 0;
    acc = 0;
    n = 0, xy = 0;
  }

  GSLInterpolation( ) : acc(0), spline(0), n(0), xy(0) {}

  void init(const KN_< double > &x, const KN_< double > &f, bool INIT = false, long cas = 0) {
    typedef const gsl_interp_type *cgsl_interpp;
#if (GSL_MAJOR_VERSION < 2)
    static const cgsl_interpp gsl_interp_steffen = gsl_interp_cspline;
#endif

    static cgsl_interpp interp[] = {gsl_interp_cspline,       gsl_interp_akima,
                                    gsl_interp_steffen,       gsl_interp_linear,
                                    gsl_interp_polynomial,    gsl_interp_cspline_periodic,
                                    gsl_interp_akima_periodic};
    if (INIT) {
      destroy( );
    }

    ffassert(x.N( ) == f.N( ));
    n = x.N( );
    splinetype = interp[cas];
    xy = new double[n * 2];

    for (long k = 0; k < n; ++k) {
      xy[k] = x[k];
      xy[k + n] = f[k];
    }

    spline = gsl_spline_alloc(splinetype, n);
    gsl_spline_init(spline, xy, xy + n, n);
  }

  void init(const KNM_< double > &kxy, bool INIT = false, long cas = 0) {
    init(kxy(0, ':'), kxy(1, ':'), INIT, cas);
  }

  void init(GSLInterpolation *g, bool INIT = false) {
    if (INIT) {
      destroy( );
    }

    n = g->n;

    xy = new double[n * 2];
    splinetype = g->splinetype;

    for (int i = 0; i < 2 * n; ++i) {
      xy[i] = g->xy[i];
    }

    spline = gsl_spline_alloc(splinetype, n);
    gsl_spline_init(spline, xy, xy + n, n);
  }

  double bb(double xi) const { return max(min(xi, xy[n - 1]), xy[0]); }

  double eval(double xi) { return gsl_spline_eval(spline, bb(xi), acc); }

  double deval(double xi) { return gsl_spline_eval_deriv(spline, bb(xi), acc); }

  double ddeval(double xi) { return gsl_spline_eval_deriv2(spline, bb(xi), acc); }

  void destroy( ) {
    if (spline) {
      gsl_spline_free(spline);
    }

    if (acc) {
      gsl_interp_accel_free(acc);
    }

    if (xy) {
      delete[] xy;
    }

    init( );
  }

  ~GSLInterpolation( ) { destroy( ); }

 private:    // no copy operator
  GSLInterpolation(const GSLInterpolation &);
  void operator=(const GSLInterpolation &);
};
struct dGSLInterpolation {
  GSLInterpolation *p;
  dGSLInterpolation(GSLInterpolation *pp) : p(pp) {}
};
struct ddGSLInterpolation {
  GSLInterpolation *p;
  ddGSLInterpolation(GSLInterpolation *pp) : p(pp) {}
};

dGSLInterpolation dGSLInterpolationedef(GSLInterpolation *p) { return dGSLInterpolation(p); }

ddGSLInterpolation ddGSLInterpolationedef(GSLInterpolation *p) { return ddGSLInterpolation(p); }

GSLInterpolation *init_GSLInterpolation(GSLInterpolation *const &gi, KNM_< double > const &a) {
  gi->init(a, true);
  return gi;
}

GSLInterpolation *init_GSLInterpolation(GSLInterpolation *const &gi, KN_< double > const &a,
                                        KN_< double > const &b) {
  gi->init(a, b, true);
  return gi;
}

GSLInterpolation *init_GSLInterpolation(GSLInterpolation *const &gi, long const &cas,
                                        KNM_< double > const &a) {
  gi->init(a, true, cas);
  return gi;
}

GSLInterpolation *init_GSLInterpolation(GSLInterpolation *const &gi, long const &cas,
                                        KN_< double > const &a, KN_< double > const &b) {
  gi->init(a, b, true, cas);
  return gi;
}

GSLInterpolation *set_GSLInterpolation(GSLInterpolation *const &gi, KNM_< double > const &a) {
  gi->init(a, false);
  return gi;
}

GSLInterpolation *set_GSLInterpolation(GSLInterpolation *const &gi, GSLInterpolation *const &gj) {
  gi->init(gj, false);
  return gi;
}

GSLInterpolation *set_GSLInterpolation(GSLInterpolation *const &gi, KN_< double > const &a,
                                       KN_< double > const &b) {
  gi->init(a, b, false);
  return gi;
}

GSLInterpolation *set_GSLInterpolation(GSLInterpolation *const &gi, long const &cas,
                                       KNM_< double > const &a) {
  gi->init(a, false, cas);
  return gi;
}

GSLInterpolation *set_GSLInterpolation(GSLInterpolation *const &gi, long const &cas,
                                       KN_< double > const &a, KN_< double > const &b) {
  gi->init(a, b, false, cas);
  return gi;
}

double GSLInterpolationeval(GSLInterpolation *gi, double x) { return gi->eval(x); }

double dGSLInterpolationeval(dGSLInterpolation gi, double x) { return gi.p->deval(x); }

double ddGSLInterpolationeval(ddGSLInterpolation gi, double x) { return gi.p->ddeval(x); }

long gslpolysolvequadratic(KN_< double > a, KN_< double > x) {
  ffassert(a.N( ) > 2 && x.N( ) > 1);
  return gsl_poly_solve_quadratic(a[2], a[1], a[0], &(x[0]), &(x[1]));
}

long gslpolysolvecubic(KN_< double > a, KN_< double > x) {
  ffassert(a.N( ) > 2 && x.N( ) > 2);
  return gsl_poly_solve_cubic(a[2], a[1], a[0], &(x[0]), &(x[1]), &(x[2]));
}

long gslpolycomplexsolve(KN_< double > a, KN_< Complex > x) {
  int n = a.N( );

  ffassert(n - 1 <= x.N( ));
  KN< double > z(n * 2);
  gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(n);
  int ok = gsl_poly_complex_solve(&a[0], n, w, &z[0]);
  gsl_poly_complex_workspace_free(w);

  for (long i = 0; i < n - 1; i++) {
    x[i] = Complex(z[2 * i], z[2 * i + 1]);
  }

  return ok;
}

// Ramdom part..
AnyType init_GSLInterpolation(Stack, const AnyType &x) {
  GSLInterpolation *p = PGetAny< GSLInterpolation >(x);

  p->init( );
  return x;
}

AnyType delete_GSLInterpolation(Stack, const AnyType &x) {
  GSLInterpolation *p = PGetAny< GSLInterpolation >(x);

  p->destroy( );
  return x;
}

AnyType init_gsl_rng(Stack, const AnyType &x) {
  gsl_rng **pp = PGetAny< gsl_rng * >(x);

  *pp = gsl_rng_alloc(gsl_rng_default);
  return x;
};
AnyType delete_gsl_rng(Stack, const AnyType &x) {
  gsl_rng **pp = PGetAny< gsl_rng * >(x);

  if (*pp) {
    gsl_rng_free(*pp);
  }

  *pp = 0;
  return Nothing;
};

gsl_rng **init_gsl_rng_type(gsl_rng **pp, const gsl_rng_type *g) {
  *pp = gsl_rng_alloc(g);
  return pp;
}

gsl_rng **set_gsl_rng_type(gsl_rng **pp, const gsl_rng_type *g) {
  if (*pp) {
    gsl_rng_free(*pp);
  }

  *pp = gsl_rng_alloc(g);
  return pp;
}

gsl_rng **set_gsl_cpy(gsl_rng **pp, gsl_rng **gg) {
  if (*pp) {
    gsl_rng_free(*pp);
  }

  *pp = gsl_rng_clone(*gg);
  return pp;
}

double gslrnguniform(gsl_rng **pr) { return gsl_rng_uniform(*pr); }

double gslrnguniformpos(gsl_rng **pr) { return gsl_rng_uniform_pos(*pr); }

long gsl_rng_get(gsl_rng **pr) { return gsl_rng_get(*pr); }

long gsl_rng_min(gsl_rng **pr) { return gsl_rng_min(*pr); }

long gsl_rng_max(gsl_rng **pr) { return gsl_rng_max(*pr); }

long gsl_rng_set(gsl_rng **pr, long s) {
  gsl_rng_set(*pr, s);
  return 0;
}

string *gsl_name(Stack s, const gsl_rng_type *const &pr) {
  return Add2StackOfPtr2Free(s, new string((*pr).name));
}

long ngslrng = 0;
long gslabort = 1;
static const gsl_rng_type **gsl_rngpp;
const gsl_rng_type *gslrngtype(long i) {
  ffassert(i >= 0 && i < ngslrng);
  return gsl_rngpp[i];
}

extern "C" {
void ffhandler(const char *reason, const char *file, int line, int gsl_errno);
}
void ffhandler(const char *reason, const char *file, int line, int gsl_errno) {
  cerr << "\n GSL Error = " << reason << " in " << file << " at " << line << " err= " << gsl_errno
       << endl;
  if (gslabort) {
    ExecError("Gsl errorhandler");
  }
}

using namespace Fem2D;
static void Load_Init( ) {
  Global.Add("gslpolysolvequadratic", "(",
             new OneOperator2< long, KN_< double >, KN_< double > >(gslpolysolvequadratic));
  Global.Add("gslpolysolvecubic", "(",
             new OneOperator2< long, KN_< double >, KN_< double > >(gslpolysolvecubic));
  Global.Add("gslpolycomplexsolve", "(",
             new OneOperator2< long, KN_< double >, KN_< Complex > >(gslpolycomplexsolve));
  /* spline gsl and June 2013 */
  /*
   * Dcl_Type<gsl_bspline_workspace**>(::InitializePtr<gsl_bspline_workspace
   * **>,::DeletePtr<gsl_bspline_workspace **>);
   * zzzfff->Add("gslbspline",atype<gsl_bspline_workspace ** >());
   * TheOperators->Add("<-",
   *                new OneOperator3_<gsl_bspline_workspace **,gsl_bspline_workspace **,KNM_<double>
   * >(pBuilQFd<R1>),
   */
  // a faire ... interface randon of gsl ...
  gsl_rng_env_setup( );
  gsl_rngpp = gsl_rng_types_setup( );

  for (long i = 0; gsl_rngpp[i];) {
    ngslrng = ++i;
  }

  Dcl_Type< gsl_rng ** >(init_gsl_rng, delete_gsl_rng);
  Dcl_Type< GSLInterpolation * >(init_GSLInterpolation, delete_GSLInterpolation);
  Dcl_Type< dGSLInterpolation >( );
  Dcl_Type< ddGSLInterpolation >( );

  Dcl_Type< const gsl_rng_type * >( );    // gsl_rng_type
  Global.New("ngslrng", CConstant< long >(ngslrng));

  zzzfff->Add("gslrng", atype< gsl_rng ** >( ));
  zzzfff->Add("gslspline", atype< GSLInterpolation * >( ));

  TheOperators->Add(
    "<-", new OneOperator2< gsl_rng **, gsl_rng **, const gsl_rng_type * >(init_gsl_rng_type));
  TheOperators->Add("<-",
                    new OneOperator2_< GSLInterpolation *, GSLInterpolation *, KNM_< double > >(
                      init_GSLInterpolation));
  TheOperators->Add("=",
                    new OneOperator2_< GSLInterpolation *, GSLInterpolation *, KNM_< double > >(
                      set_GSLInterpolation));
  TheOperators->Add("=",
                    new OneOperator2_< GSLInterpolation *, GSLInterpolation *, GSLInterpolation * >(
                      set_GSLInterpolation));

  TheOperators->Add(
    "<-", new OneOperator3_< GSLInterpolation *, GSLInterpolation *, KN_< double >, KN_< double > >(
            init_GSLInterpolation));
  TheOperators->Add(
    "<-", new OneOperator3_< GSLInterpolation *, GSLInterpolation *, long, KNM_< double > >(
            init_GSLInterpolation));
  TheOperators->Add(
    "<-",
    new OneOperator4_< GSLInterpolation *, GSLInterpolation *, long, KN_< double >, KN_< double > >(
      init_GSLInterpolation));
  Add< GSLInterpolation * >(
    "(", "", new OneOperator2< double, GSLInterpolation *, double >(GSLInterpolationeval));
  Add< GSLInterpolation * >(
    "d", ".", new OneOperator1< dGSLInterpolation, GSLInterpolation * >(dGSLInterpolationedef));
  Add< GSLInterpolation * >(
    "dd", ".", new OneOperator1< ddGSLInterpolation, GSLInterpolation * >(ddGSLInterpolationedef));
  Add< dGSLInterpolation >(
    "(", "", new OneOperator2< double, dGSLInterpolation, double >(dGSLInterpolationeval));
  Add< ddGSLInterpolation >(
    "(", "", new OneOperator2< double, ddGSLInterpolation, double >(ddGSLInterpolationeval));

  TheOperators->Add(
    "=", new OneOperator2< gsl_rng **, gsl_rng **, const gsl_rng_type * >(set_gsl_rng_type));
  TheOperators->Add("=", new OneOperator2< gsl_rng **, gsl_rng **, gsl_rng ** >(set_gsl_cpy));

  Global.Add("gslrnguniform", "(", new OneOperator1< double, gsl_rng ** >(gslrnguniform));
  Global.Add("gslrnguniformpos", "(", new OneOperator1< double, gsl_rng ** >(gslrnguniformpos));

  Global.Add("gslname", "(", new OneOperator1s_< string *, const gsl_rng_type * >(gsl_name));
  Global.Add("gslrngget", "(", new OneOperator1< long, gsl_rng ** >(gsl_rng_get));
  Global.Add("gslrngmin", "(", new OneOperator1< long, gsl_rng ** >(gsl_rng_min));
  Global.Add("gslrngmax", "(", new OneOperator1< long, gsl_rng ** >(gsl_rng_max));
  Global.Add("gslrngset", "(", new OneOperator2< long, gsl_rng **, long >(gsl_rng_set));
  Global.Add("gslrngtype", "(", new OneOperator1< const gsl_rng_type *, long >(gslrngtype));
  init_gsl_sf( );
  gslabort = 1;
  Global.New("gslabortonerror", CConstant< long * >(&gslabort));
  // static  cgsl_interpp interp[] =
  // {gsl_interp_cspline,gsl_interp_akima,gsl_interp_steffen,gsl_interp_linear,gsl_interp_polynomial,gsl_interp_cspline_periodic,gsl_interp_akima_periodic};

  // type of spline of gsl ?????
  Global.New("gslinterpcspline", CConstant< long >(0));
  Global.New("gslinterpakima", CConstant< long >(1));
  Global.New("gslinterpsteffen", CConstant< long >(2));
  Global.New("gslinterplinear", CConstant< long >(3));
  Global.New("gslinterppolynomial", CConstant< long >(4));
  Global.New("gslinterpcsplineperiodic", CConstant< long >(5));
  Global.New("gslinterpakimaperiodic", CConstant< long >(6));

  gsl_set_error_handler(ffhandler);
}

LOADFUNC(Load_Init)
