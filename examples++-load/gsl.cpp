// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$
//ff-c++-LIBRARY-dep:   gsl
//ff-c++-cpp-dep:  
#include <ff++.hpp>
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
#include <gsl/gsl_multifit.h>

#include "ff_gsl_awk.hpp"

long gslpolysolvequadratic( KN_<double> a,  KN_<double> x)
{
  ffassert(a.N()>2 && x.N()>1);
  return gsl_poly_solve_quadratic (a[2],a[1],a[0],&(x[0]),&(x[1]));
}
long gslpolysolvecubic( KN_<double> a,  KN_<double> x)
{
  ffassert(a.N()>2 && x.N()>2);
  return gsl_poly_solve_cubic (a[2],a[1],a[0],&(x[0]),&(x[1]),&(x[2]));
  }

long gslpolycomplexsolve( KN_<double> a,  KN_<Complex> x)
{
  int n = a.N();
  ffassert( n-1 <=  x.N()); 
  KN<double> z(n*2); 
  gsl_poly_complex_workspace * w= gsl_poly_complex_workspace_alloc (n);
  int ok=gsl_poly_complex_solve (&a[0], n, w, &z[0]);
  gsl_poly_complex_workspace_free (w);
  for (long i = 0; i < n-1; i++)
    x[i] = Complex(z[2*i], z[2*i+1]);
  return ok; 
}

class Init { public:
  Init();
};
LOADINIT(Init);
Init::Init(){
  init_gsl_sf() ;
  Global.Add("gslpolysolvequadratic","(",new OneOperator2<long,KN_<double>,KN_<double> >( gslpolysolvequadratic));
  Global.Add("gslpolysolvecubic","(",new OneOperator2<long,KN_<double>,KN_<double> >(gslpolysolvecubic));
  Global.Add("gslpolycomplexsolve","(",new OneOperator2<long,KN_<double>,KN_<Complex> >( gslpolycomplexsolve));
/* spline gsl and June 2013 */
    /*
    Dcl_Type<gsl_bspline_workspace**>(::InitializePtr<gsl_bspline_workspace **>,::DeletePtr<gsl_bspline_workspace **>);
    zzzfff->Add("gslbspline",atype<gsl_bspline_workspace ** >());
    TheOperators->Add("<-",
                      new OneOperator3_<gsl_bspline_workspace **,gsl_bspline_workspace **,KNM_<double>  >(pBuilQFd<R1>),
*/
}
