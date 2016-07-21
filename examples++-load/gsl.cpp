// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$
//ff-c++-LIBRARY-dep:   gsl
//ff-c++-cpp-dep:  
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

struct  GSLInterpolation  {
    
    gsl_interp_accel *acc;
    gsl_spline * spline ;
    double *xy;
    size_t n;
    const gsl_interp_type *splinetype;
    void init() {spline=0;acc=0;n=0,xy=0;}
    GSLInterpolation() :acc(0),spline(0),n(0),xy(0) {}
    void init(const KN_<double> & x,const KN_<double> & f,bool INIT=false,long cas=0)
    {
        typedef const gsl_interp_type * cgsl_interpp ;
        static  cgsl_interpp interp[] = {gsl_interp_cspline,gsl_interp_akima,gsl_interp_steffen,gsl_interp_linear,gsl_interp_polynomial,gsl_interp_cspline_periodic,gsl_interp_akima_periodic};
        if(INIT) destroy();
        ffassert(x.N()==f.N());
        n = x.N();
        splinetype=interp[cas];
        xy= new double[n*2];
        
        for(long k=0;k<n;++k)
        {
            xy[k]=x[k];
            xy[k+n] = f[k];
        }
        spline = gsl_spline_alloc (splinetype, n);
        gsl_spline_init (spline, xy,xy+n,n);
    }
    void init(const KNM_<double> & kxy,bool INIT=false,long cas=0)
    {
        init(kxy(0,':'),kxy(1,':'),INIT,cas);
    }
    void init(GSLInterpolation * g,bool INIT=false)
    {
        if(INIT) destroy();
        n=g->n;
        
        xy= new double[n*2];
        splinetype=g->splinetype;
        for(int i=0; i< 2*n; ++i)
            xy[i]= g->xy[i];
        spline = gsl_spline_alloc (splinetype, n);
        gsl_spline_init (spline, xy,xy+n,n);

    }

    double eval(double xi) { return gsl_spline_eval (spline, xi, acc);}
    double deval(double xi) { return gsl_spline_eval_deriv (spline, xi, acc);}
    double ddeval(double xi) { return gsl_spline_eval_deriv2 (spline, xi, acc);}
    
    void destroy()
    {
        if(spline) gsl_spline_free (spline);
        if(acc) gsl_interp_accel_free (acc);
        if(xy) delete [] xy;
        init();
    }

    ~GSLInterpolation() {destroy();}
private: // no copy operator
    GSLInterpolation(const GSLInterpolation  &);
    void operator=(const GSLInterpolation  &);
};
struct dGSLInterpolation {
    GSLInterpolation * p;
    dGSLInterpolation(GSLInterpolation *pp) : p(pp){}
};
struct ddGSLInterpolation { GSLInterpolation * p;
    ddGSLInterpolation(GSLInterpolation *pp) : p(pp){}
};

dGSLInterpolation dGSLInterpolationedef(GSLInterpolation *p)
{return dGSLInterpolation(p);}
ddGSLInterpolation ddGSLInterpolationedef(GSLInterpolation *p)
{return ddGSLInterpolation(p);}

GSLInterpolation * init_GSLInterpolation(GSLInterpolation * const & gi, KNM_<double>  const & a)
{
    gi->init(a,true);
    return gi;
}
GSLInterpolation * init_GSLInterpolation(GSLInterpolation *  const & gi, KN_<double>  const &  a,KN_<double>  const &  b)
{
    gi->init(a,b,true);
    return gi;
}
GSLInterpolation * init_GSLInterpolation(GSLInterpolation * const & gi,long const & cas , KNM_<double>  const & a)
{
    gi->init(a,true,cas);
    return gi;
}
GSLInterpolation * init_GSLInterpolation(GSLInterpolation *  const & gi,long const & cas ,  KN_<double>  const &  a,KN_<double>  const &  b)
{
    gi->init(a,b,true,cas );
    return gi;
}

GSLInterpolation * set_GSLInterpolation(GSLInterpolation * const & gi, KNM_<double>  const & a)
{
    gi->init(a,false);
    return gi;
}
GSLInterpolation * set_GSLInterpolation(GSLInterpolation * const & gi,GSLInterpolation * const & gj)
{
    gi->init(gj,false);
    return gi;
}

GSLInterpolation * set_GSLInterpolation(GSLInterpolation *  const & gi, KN_<double>  const &  a,KN_<double>  const &  b)
{
    gi->init(a,b,false);
    return gi;
}
GSLInterpolation * set_GSLInterpolation(GSLInterpolation * const & gi,long const & cas , KNM_<double>  const & a)
{
    gi->init(a,false,cas);
    return gi;
}
GSLInterpolation * set_GSLInterpolation(GSLInterpolation *  const & gi,long const & cas ,  KN_<double>  const &  a,KN_<double>  const &  b)
{
    gi->init(a,b,false,cas );
    return gi;
}


double GSLInterpolationeval(GSLInterpolation *  gi,double x)
{
    return gi->eval(x);
}
double dGSLInterpolationeval(dGSLInterpolation  gi,double x)
{
    return gi.p->deval(x);
}
double ddGSLInterpolationeval(ddGSLInterpolation  gi,double x)
{
    return gi.p->ddeval(x);
}

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

//  Ramdom part..
AnyType init_GSLInterpolation(Stack,const AnyType &x){
    GSLInterpolation *p = PGetAny< GSLInterpolation>(x);
    p->init();
    return x;
}

AnyType delete_GSLInterpolation(Stack,const AnyType &x){
    GSLInterpolation *p = PGetAny< GSLInterpolation>(x);
    p->destroy();
    return x;
}

AnyType  init_gsl_rng(Stack,const AnyType &x){
  gsl_rng ** pp = PGetAny< gsl_rng *>(x);
    *pp = gsl_rng_alloc(gsl_rng_default) ;
    return  x;
};
AnyType delete_gsl_rng(Stack,const AnyType &x)
{
 gsl_rng ** pp = PGetAny< gsl_rng *>(x);
 if(*pp) gsl_rng_free(*pp);
 *pp=0;
 return  Nothing;
};

gsl_rng ** init_gsl_rng_type( gsl_rng ** pp, const gsl_rng_type * g)
{
    *pp = gsl_rng_alloc(g) ;
    return pp;
}

gsl_rng ** set_gsl_rng_type( gsl_rng ** pp, const gsl_rng_type * g)
{
    if(*pp) gsl_rng_free(*pp);
    *pp = gsl_rng_alloc(g) ;
    return pp;
}
gsl_rng ** set_gsl_cpy( gsl_rng ** pp, gsl_rng ** gg)
{
    if(*pp) gsl_rng_free(*pp);
    *pp =  gsl_rng_clone(*gg);
    return pp;
}

double gslrnguniform( gsl_rng ** pr) { return gsl_rng_uniform(*pr);}
double gslrnguniformpos( gsl_rng ** pr) { return gsl_rng_uniform_pos(*pr);}
long gsl_rng_get(gsl_rng ** pr){ return gsl_rng_get(*pr);}
long gsl_rng_min(gsl_rng ** pr){ return gsl_rng_min(*pr);}
long gsl_rng_max(gsl_rng ** pr){ return gsl_rng_max(*pr);}
long gsl_rng_set(gsl_rng ** pr, long s){  gsl_rng_set(*pr,s);return 0; }
string * gsl_name(Stack s,const gsl_rng_type * const & pr)
  {return Add2StackOfPtr2Free(s,new string((*pr).name));}

long  ngslrng =0;
long  gslabort =1;
static const gsl_rng_type ** gsl_rngpp; 

const gsl_rng_type * gslrngtype(long i)
{ 
  ffassert(i >=0 && i < ngslrng);
  return gsl_rngpp[i]; 
}

extern "C" {
    void ffhandler (const char * reason,
                  const char * file,
                  int line,
                    int gsl_errno);

}
void ffhandler (const char * reason,
                const char * file,
                int line,
                int gsl_errno)
{
    cerr << "\n GSL Error = " << reason << " in " <<file << " at " << line << " err= " <<gsl_errno << endl;
    if(gslabort) ExecError("Gsl errorhandler");
}

/*  class Init { public:
  Init();
};
$1 */
using  namespace Fem2D ;
static void Load_Init(){
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
    // a faire ... interface randon of gsl ... 
  gsl_rng_env_setup();
  gsl_rngpp =gsl_rng_types_setup(); 
  for(long  i=0; gsl_rngpp[i]; )
    ngslrng=++i; 
   Dcl_Type< gsl_rng **  > (init_gsl_rng,delete_gsl_rng);
   Dcl_Type< GSLInterpolation *  > (init_GSLInterpolation,delete_GSLInterpolation);
   Dcl_Type< dGSLInterpolation   > ();
   Dcl_Type< ddGSLInterpolation   > ();
    
   Dcl_Type< const gsl_rng_type * > (); //gsl_rng_type
   Global.New("ngslrng",CConstant<long>(ngslrng)); 
//
// all gsl random generator .???? .
/*
Global.New("gslrngborosh13",CConstant<const gsl_rng_type *>(gsl_rng_borosh13));
Global.New("gslrngcoveyou",CConstant<const gsl_rng_type *>(gsl_rng_coveyou));
Global.New("gslrngcmrg",CConstant<const gsl_rng_type *>(gsl_rng_cmrg));
Global.New("gslrngfishman18",CConstant<const gsl_rng_type *>(gsl_rng_fishman18));
Global.New("gslrngfishman20",CConstant<const gsl_rng_type *>(gsl_rng_fishman20));
Global.New("gslrngfishman2x",CConstant<const gsl_rng_type *>(gsl_rng_fishman2x));
Global.New("gslrnggfsr4",CConstant<const gsl_rng_type *>(gsl_rng_gfsr4));
Global.New("gslrngknuthran",CConstant<const gsl_rng_type *>(gsl_rng_knuthran));
Global.New("gslrngknuthran2",CConstant<const gsl_rng_type *>(gsl_rng_knuthran2));
Global.New("gslrngknuthran2002",CConstant<const gsl_rng_type *>(gsl_rng_knuthran2002));
Global.New("gslrnglecuyer21",CConstant<const gsl_rng_type *>(gsl_rng_lecuyer21));
Global.New("gslrngminstd",CConstant<const gsl_rng_type *>(gsl_rng_minstd));
Global.New("gslrngmrg",CConstant<const gsl_rng_type *>(gsl_rng_mrg));
Global.New("gslrngmt19937",CConstant<const gsl_rng_type *>(gsl_rng_mt19937));
Global.New("gslrngmt199371999",CConstant<const gsl_rng_type *>(gsl_rng_mt19937_1999));
Global.New("gslrngmt199371998",CConstant<const gsl_rng_type *>(gsl_rng_mt19937_1998));
Global.New("gslrngr250",CConstant<const gsl_rng_type *>(gsl_rng_r250));
Global.New("gslrngran0",CConstant<const gsl_rng_type *>(gsl_rng_ran0));
Global.New("gslrngran1",CConstant<const gsl_rng_type *>(gsl_rng_ran1));
Global.New("gslrngran2",CConstant<const gsl_rng_type *>(gsl_rng_ran2));
Global.New("gslrngran3",CConstant<const gsl_rng_type *>(gsl_rng_ran3));
Global.New("gslrngrand",CConstant<const gsl_rng_type *>(gsl_rng_rand));
Global.New("gslrngrand48",CConstant<const gsl_rng_type *>(gsl_rng_rand48));
Global.New("gslrngrandom128bsd",CConstant<const gsl_rng_type *>(gsl_rng_random128_bsd));
Global.New("gslrngrandom128glibc2",CConstant<const gsl_rng_type *>(gsl_rng_random128_glibc2));
Global.New("gslrngrandom128libc5",CConstant<const gsl_rng_type *>(gsl_rng_random128_libc5));
Global.New("gslrngrandom256bsd",CConstant<const gsl_rng_type *>(gsl_rng_random256_bsd));
Global.New("gslrngrandom256glibc2",CConstant<const gsl_rng_type *>(gsl_rng_random256_glibc2));
Global.New("gslrngrandom256libc5",CConstant<const gsl_rng_type *>(gsl_rng_random256_libc5));
Global.New("gslrngrandom32bsd",CConstant<const gsl_rng_type *>(gsl_rng_random32_bsd));
Global.New("gslrngrandom32glibc2",CConstant<const gsl_rng_type *>(gsl_rng_random32_glibc2));
Global.New("gslrngrandom32libc5",CConstant<const gsl_rng_type *>(gsl_rng_random32_libc5));
Global.New("gslrngrandom64bsd",CConstant<const gsl_rng_type *>(gsl_rng_random64_bsd));
Global.New("gslrngrandom64glibc2",CConstant<const gsl_rng_type *>(gsl_rng_random64_glibc2));
Global.New("gslrngrandom64libc5",CConstant<const gsl_rng_type *>(gsl_rng_random64_libc5));
Global.New("gslrngrandom8bsd",CConstant<const gsl_rng_type *>(gsl_rng_random8_bsd));
Global.New("gslrngrandom8glibc2",CConstant<const gsl_rng_type *>(gsl_rng_random8_glibc2));
Global.New("gslrngrandom8libc5",CConstant<const gsl_rng_type *>(gsl_rng_random8_libc5));
Global.New("gslrngrandombsd",CConstant<const gsl_rng_type *>(gsl_rng_random_bsd));
Global.New("gslrngrandomglibc2",CConstant<const gsl_rng_type *>(gsl_rng_random_glibc2));
Global.New("gslrngrandomlibc5",CConstant<const gsl_rng_type *>(gsl_rng_random_libc5));
Global.New("gslrngrandu",CConstant<const gsl_rng_type *>(gsl_rng_randu));
Global.New("gslrngranf",CConstant<const gsl_rng_type *>(gsl_rng_ranf));
Global.New("gslrngranlux",CConstant<const gsl_rng_type *>(gsl_rng_ranlux));
Global.New("gslrngranlux389",CConstant<const gsl_rng_type *>(gsl_rng_ranlux389));
Global.New("gslrngranlxd1",CConstant<const gsl_rng_type *>(gsl_rng_ranlxd1));
Global.New("gslrngranlxd2",CConstant<const gsl_rng_type *>(gsl_rng_ranlxd2));
Global.New("gslrngranlxs0",CConstant<const gsl_rng_type *>(gsl_rng_ranlxs0));
Global.New("gslrngranlxs1",CConstant<const gsl_rng_type *>(gsl_rng_ranlxs1));
Global.New("gslrngranlxs2",CConstant<const gsl_rng_type *>(gsl_rng_ranlxs2));
Global.New("gslrngranmar",CConstant<const gsl_rng_type *>(gsl_rng_ranmar));
Global.New("gslrngslatec",CConstant<const gsl_rng_type *>(gsl_rng_slatec));
Global.New("gslrngtaus",CConstant<const gsl_rng_type *>(gsl_rng_taus));
Global.New("gslrngtaus2",CConstant<const gsl_rng_type *>(gsl_rng_taus2));
Global.New("gslrngtaus113",CConstant<const gsl_rng_type *>(gsl_rng_taus113));
Global.New("gslrngtransputer",CConstant<const gsl_rng_type *>(gsl_rng_transputer));
Global.New("gslrngtt800",CConstant<const gsl_rng_type *>(gsl_rng_tt800));
Global.New("gslrnguni",CConstant<const gsl_rng_type *>(gsl_rng_uni));
Global.New("gslrnguni32",CConstant<const gsl_rng_type *>(gsl_rng_uni32));
Global.New("gslrngvax",CConstant<const gsl_rng_type *>(gsl_rng_vax));
Global.New("gslrngwaterman14",CConstant<const gsl_rng_type *>(gsl_rng_waterman14));
Global.New("gslrngzuf",CConstant<const gsl_rng_type *>(gsl_rng_zuf));
Global.New("gslrngdefault",CConstant<const gsl_rng_type *>(gsl_rng_default));
*/    
    
zzzfff->Add("gslrng",atype<gsl_rng ** >());
zzzfff->Add("gslspline",atype<GSLInterpolation * >());
    
TheOperators->Add("<-",new OneOperator2<gsl_rng  **,gsl_rng  **, const gsl_rng_type *  >(init_gsl_rng_type));
TheOperators->Add("<-",new OneOperator2_<GSLInterpolation *,GSLInterpolation *,KNM_<double> >(init_GSLInterpolation));
TheOperators->Add("=",new OneOperator2_<GSLInterpolation *,GSLInterpolation *,KNM_<double> >(set_GSLInterpolation));
TheOperators->Add("=",new OneOperator2_<GSLInterpolation *,GSLInterpolation *,GSLInterpolation *>(set_GSLInterpolation));
    
TheOperators->Add("<-",new OneOperator3_<GSLInterpolation *,GSLInterpolation *,KN_<double> ,KN_<double> >(init_GSLInterpolation));
TheOperators->Add("<-",new OneOperator3_<GSLInterpolation *,GSLInterpolation *,long, KNM_<double> >(init_GSLInterpolation));
TheOperators->Add("<-",new OneOperator4_<GSLInterpolation *,GSLInterpolation *,long,KN_<double> ,KN_<double> >(init_GSLInterpolation));
Add<GSLInterpolation * >("(","",new OneOperator2<double ,GSLInterpolation *,double >(GSLInterpolationeval));
Add<GSLInterpolation * >("d",".",new OneOperator1<dGSLInterpolation  ,GSLInterpolation*>(dGSLInterpolationedef));
Add<GSLInterpolation * >("dd",".",new OneOperator1<ddGSLInterpolation  ,GSLInterpolation*>(ddGSLInterpolationedef));
Add<dGSLInterpolation  >("(","",new OneOperator2<double ,dGSLInterpolation ,double >(dGSLInterpolationeval));
Add<ddGSLInterpolation  >("(","",new OneOperator2<double ,ddGSLInterpolation ,double >(ddGSLInterpolationeval));

TheOperators->Add("=",new OneOperator2<gsl_rng  **,gsl_rng  **, const gsl_rng_type *  >(set_gsl_rng_type));
TheOperators->Add("=",new OneOperator2<gsl_rng  **,gsl_rng  **, gsl_rng  **   >(set_gsl_cpy));
//map_type[typeid(gsl_rng *).name()]->AddCast(   new E_F1_funcT<gsl_rng *,gsl_rng **>(UnRef<gsl_rng*>) );
//map_type[typeid(gsl_rng *).name()]->AddCast(   new E_F1_funcT<gsl_rng *,gsl_rng **>(UnRef<gsl_rng*>) );

Global.Add("gslrnguniform","(",new OneOperator1<double,gsl_rng **>( gslrnguniform));
Global.Add("gslrnguniformpos","(",new OneOperator1<double,gsl_rng **>( gslrnguniformpos));
    
Global.Add("gslname","(",new OneOperator1s_<string * ,const gsl_rng_type *>( gsl_name));
Global.Add("gslrngget","(",new OneOperator1<long   ,gsl_rng **>( gsl_rng_get));
Global.Add("gslrngmin","(",new OneOperator1<long   ,gsl_rng **>( gsl_rng_min));
Global.Add("gslrngmax","(",new OneOperator1<long   ,gsl_rng **>( gsl_rng_max));
Global.Add("gslrngset","(",new OneOperator2<long   ,gsl_rng **, long>(gsl_rng_set));
 Global.Add("gslrngtype","(",new OneOperator1<const gsl_rng_type * ,long>(gslrngtype));     
  init_gsl_sf() ;
 gslabort=1;
 Global.New("gslabortonerror",CConstant<long*>(&gslabort));
    //          static  cgsl_interpp interp[] = {gsl_interp_cspline,gsl_interp_akima,gsl_interp_steffen,gsl_interp_linear,gsl_interp_polynomial,gsl_interp_cspline_periodic,gsl_interp_akima_periodic};

 // type of spline of gsl ?????
 Global.New("gslinterpcspline",CConstant<long>(0));
 Global.New("gslinterpakima",CConstant<long>(1));
 Global.New("gslinterpsteffen",CConstant<long>(2));
 Global.New("gslinterplinear",CConstant<long>(3));
 Global.New("gslinterppolynomial",CConstant<long>(4));
 Global.New("gslinterpcsplineperiodic",CConstant<long>(5));
 Global.New("gslinterpakimaperiodic",CConstant<long>(6));
    
    
    gsl_set_error_handler(ffhandler);
}
LOADFUNC(Load_Init)
