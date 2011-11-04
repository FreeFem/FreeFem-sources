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

/*
// wrapper gsl (all int in feefem++ are long)
double gsl_sf_bessel_Jn (long n, double x) {return gsl_sf_bessel_Jn((int) n,x); }
double gsl_sf_bessel_jn (long n, double x) {return gsl_sf_bessel_jn((int) n,x); }
double gsl_sf_bessel_yl (long n, double x) {return gsl_sf_bessel_yl((int) n,x); }
double gsl_sf_bessel_kl_scaled (long n, double x) {return gsl_sf_bessel_kl_scaled((int) n,x); }
double gsl_sf_bessel_Yn (long n, double x) {return gsl_sf_bessel_Yn((int) n,x); }
double gsl_sf_bessel_Kn (long n, double x) {return gsl_sf_bessel_Kn((int) n,x); }
double gsl_sf_bessel_In (long n, double x) {return gsl_sf_bessel_In((int) n,x); }
double gsl_sf_bessel_zero_J0 (long s) {return gsl_sf_bessel_zero_J0( (unsigned int) s); }
double gsl_sf_bessel_zero_J1 (long s) {return gsl_sf_bessel_zero_J1( (unsigned int) s); }
double gsl_sf_bessel_zero_Jnu (double nu,long s) {return gsl_sf_bessel_zero_Jnu(nu, (unsigned int) s); }
*/
#include "ff_gsl_awk.hpp"


class Init { public:
  Init();
};
LOADINIT(Init);
Init::Init(){
  init_gsl_sf() ; 
  /*
  Global.Add("gslSfBesselJ0","(",new OneOperator1<double>(gsl_sf_bessel_J0));
  Global.Add("gslSfBesselJ1","(",new OneOperator1<double>(gsl_sf_bessel_J1));
  Global.Add("gslSfBesselJn","(",new OneOperator2<double,long,double>(gsl_sf_bessel_Jn));
  Global.Add("gslSfBesselJnu","(",new OneOperator2<double,double,double>(gsl_sf_bessel_Jnu));

  Global.Add("gslSfBesselj0","(",new OneOperator1<double>(gsl_sf_bessel_j0));
  Global.Add("gslSfBesselj1","(",new OneOperator1<double>(gsl_sf_bessel_j1));
  Global.Add("gslSfBesselj2","(",new OneOperator1<double>(gsl_sf_bessel_j2));
  Global.Add("gslSfBesseljn","(",new OneOperator2<double,long,double>(gsl_sf_bessel_jn));


  Global.Add("gslSfBessely0","(",new OneOperator1<double>(gsl_sf_bessel_y0));
  Global.Add("gslSfBessely1","(",new OneOperator1<double>(gsl_sf_bessel_y1));
  Global.Add("gslSfBessely2","(",new OneOperator1<double>(gsl_sf_bessel_y2));
  Global.Add("gslSfBesselyl","(",new OneOperator2<double,long,double>(gsl_sf_bessel_yl));


  Global.Add("gslSfBesselk0Scaled","(",new OneOperator1<double>(gsl_sf_bessel_k0_scaled));
  Global.Add("gslSfBesselk1Scaled","(",new OneOperator1<double>(gsl_sf_bessel_k1_scaled));
  Global.Add("gslSfBesselk2Scaled","(",new OneOperator1<double>(gsl_sf_bessel_k2_scaled));
  Global.Add("gslSfBesselklScaled","(",new OneOperator2<double,long,double>(gsl_sf_bessel_kl_scaled));

  Global.Add("gslSfBesselY0","(",new OneOperator1<double>(gsl_sf_bessel_Y0));
  Global.Add("gslSfBesselY1","(",new OneOperator1<double>(gsl_sf_bessel_Y1));
  Global.Add("gslSfBesselYn","(",new OneOperator2<double,long,double>(gsl_sf_bessel_Yn));
  Global.Add("gslSfBesselYnu","(",new OneOperator2<double,double,double>(gsl_sf_bessel_Ynu));

  Global.Add("gslSfBesselK0","(",new OneOperator1<double>(gsl_sf_bessel_K0));
  Global.Add("gslSfBesselK1","(",new OneOperator1<double>(gsl_sf_bessel_K1));
  Global.Add("gslSfBesselKn","(",new OneOperator2<double,long,double>(gsl_sf_bessel_Kn));
  Global.Add("gslSfBesselKnu","(",new OneOperator2<double,double,double>(gsl_sf_bessel_Knu));

  Global.Add("gslSfBesselI0","(",new OneOperator1<double>(gsl_sf_bessel_I0));
  Global.Add("gslSfBesselI1","(",new OneOperator1<double>(gsl_sf_bessel_I1));
  Global.Add("gslSfBesselIn","(",new OneOperator2<double,long,double>(gsl_sf_bessel_In));
  Global.Add("gslSfBesselInu","(",new OneOperator2<double,double,double>(gsl_sf_bessel_Inu));

  Global.Add("gslSfBesselZeroJ0","(",new OneOperator1<double,long>(gsl_sf_bessel_zero_J0));
  Global.Add("gslSfBesselZeroJ1","(",new OneOperator1<double,long>(gsl_sf_bessel_zero_J1));
  Global.Add("gslSfBesselZeroJnu","(",new OneOperator2<double,double,long>(gsl_sf_bessel_zero_Jnu));
  */
}
