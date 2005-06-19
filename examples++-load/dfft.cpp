// Example C++ function "myfunction", dynamically loaded into "load.edp"
// ---------------------------------------------------------------------
// $Id$

#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;  
#include "rgraph.hpp"
#include "RNM.hpp"
#include <fftw3.h>
#include <cmath>

//using namespace Fem2D;


class DFFT_1dor2d 
{  
public: 
  Complex * x;
  int n, m;
  int sign;   
  DFFT_1dor2d(KN<Complex> * xx,long signn,long nn=1) : x(*xx),n(nn),m(xx->N()/nn),sign(signn) {
    cout << xx << " " << signn << " " << nn << " " << xx->N() << " n: " << n << " m:" << m << endl;
    ffassert(n>0 && (n*m==xx->N()));
    
  } 


  DFFT_1dor2d(KNM<Complex> * xx,long signn) : x(*xx),n(xx->M()),m(xx->N()),sign(signn) {
  }
  
};

DFFT_1dor2d dfft(KN<Complex> * const  & x,const long &sign)
{
  return DFFT_1dor2d(x,sign);
}

DFFT_1dor2d dfft(KN<Complex> *const  &x,const long &nn,const long &sign)
{
  return DFFT_1dor2d(x,sign,nn);
}

DFFT_1dor2d dfft(KNM<Complex> * const &  x,const long &sign)
{
  return DFFT_1dor2d(x,sign);
}

KN<Complex> * dfft_eq(  KN<Complex> * const   &x,const DFFT_1dor2d & d)
{
  ffassert(x->N()==d.n*d.m);
  Complex *px =  *x;
  fftw_plan p; 
  cout << " dfft " << px << " = " << d.x << " n = " << d.n << " " << d.m << " sign = " << d.sign << endl; 
  if ( d.n > 1)
    p = fftw_plan_dft_2d(d.n,d.m,reinterpret_cast<fftw_complex*>(d.x),reinterpret_cast<fftw_complex*> (px),d.sign,FFTW_ESTIMATE);
  else
    p = fftw_plan_dft_1d(d.n,reinterpret_cast<fftw_complex*>(d.x),reinterpret_cast<fftw_complex*> (px),d.sign,FFTW_ESTIMATE);
  cout << " ---" ;
  fftw_execute(p);
  cout << " ---" ;
  fftw_destroy_plan(p);
  cout << " ---" ;
  return  x;
}

class Init { public:
  Init();
};
Init init;
Init::Init(){
  cout << " lood: init dfft " << endl;
  Dcl_Type<DFFT_1dor2d >();
  Global.Add("dfft","(", new OneOperator2_<DFFT_1dor2d,KN<Complex>*,long >(dfft ));
  Global.Add("dfft","(", new OneOperator3_<DFFT_1dor2d,KN<Complex>*,long,long >(dfft ));
  Global.Add("dfft","(", new OneOperator2_<DFFT_1dor2d,KNM<Complex>*,long >(dfft ));
  TheOperators->Add("=", new OneOperator2_<KN<Complex>*,KN<Complex>*,DFFT_1dor2d>(dfft_eq));
  // TheOperators->Add("=", new OneOperator2_<KNM<Complex>*,KNM<Complex>*,DFFT_1dor2d>(dfft_eq));

}
