#include "ff++.hpp"
#include "RNM.hpp"
using namespace std;
#include <vecLib/clapack.h>
class Init { public:
  Init();
};

bool lapack_inv(KNM<double>* A)
{
  long n=A->N();
  long m=A->M();
  double *a=&(*A)(0,0);
  long info;
  long lda=n;
  KN<long> ipiv(n);
  long lw=10*n;
  KN<double> w(lw);
  ffassert(n==m);
  dgetrf_(&n,&n,a,&lda,ipiv,&info);
  if(info) return info;
  dgetri_(&n,a,&lda,ipiv,w,&lw,&info);
    return info;
}
static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  Global.Add("inv","(",new  OneOperator1<bool,KNM<double>*>(lapack_inv));  
}

