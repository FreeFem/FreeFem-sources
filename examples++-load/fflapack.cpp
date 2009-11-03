//ff-c++-LIBRARY-dep:   lapack
//ff-c++-LIBRARY-dep:   blas
#include "ff++.hpp"
#include "RNM.hpp"
using namespace std;

#ifdef __LP64__
  typedef int intblas;
#else
  typedef long intblas;
#endif

#include "clapack.h"
class Init { public:
  Init();
};

bool lapack_inv(KNM<double>* A)
{
  intblas n=A->N();
  intblas m=A->M();
  double *a=&(*A)(0,0);
  intblas info;
  intblas lda=n;
  KN<intblas> ipiv(n);
  intblas  lw=10*n;
  KN<double> w(lw);
  ffassert(n==m);
  dgetrf_(&n,&n,a,&lda,ipiv,&info);
  if(info) return info;
  dgetri_(&n,a,&lda,ipiv,w,&lw,&info);
    return info;
}
long lapack_dgeev(KNM<double> *const &A,KN<Complex> *const &vp,KNM<Complex> *const &vectp)
{
  intblas nvp =0,zero=0;
  intblas n= A->N();
  ffassert(A->M()==n);
  ffassert(vectp->M()>=n);
  ffassert(vectp->N()>=n);
  ffassert(vp->N()>=n);
  KN<double> wr(n),wi(n),vr(n*n),vl(n*n);
  KNM<double> mat(*A);
  intblas info,lw=-1;  
  KN<double> w(1);
  char N='N',V='V';
  dgeev_(&N,&V,&n,mat,&n,wr,wi,vl,&n,vr,&n,w,&lw,&info);
  lw=w[0];
  cout << lw << endl;
  w.resize(lw);
  dgeev_(&N,&V,&n,mat,&n,wr,wi,vl,&n,vr,&n,w,&lw,&info);
  if(info)
    cout << " info =  " << info << endl;
  if(!info)
    {
      int k=0;
      for(int i=0;i<n;++i)
	{
	  (*vp)[i]=Complex(wr[i],wi[i]);
	  if(verbosity>2)
	    cout << "   dgeev: vp "<< i << " : "  << (*vp)[i] << endl;
	  if( wi[i] == 0)
	    for(int j=0;j<n;++j)
	      (*vectp)(j,i)=vr[k++];
	  else if (  wi[i] >  0)
	    {
	      int ki= k+n;
	      for(int j=0;j<n;++j)
		(*vectp)(j,i)=Complex(vr[k++],vr[ki++]);	      
	    }
	  else 
	    {
	      int kr= k-n;
	      for(int j=0;j<n;++j)
		(*vectp)(j,i)=Complex(vr[kr++],-vr[k++]);	      
	    }
	  if(verbosity>5)
	    cout << "   dgeev :   " << (*vectp)(':',i) <<endl;
	}
    }
  else
    {
      (*vp)=Complex();
      (*vectp)=Complex();
    }
  return nvp;
}
static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  Global.Add("inv","(",new  OneOperator1<bool,KNM<double>*>(lapack_inv));  
  Global.Add("dgeev","(",new  OneOperator3_<long,KNM<double>*,KN<Complex>*,KNM<Complex>*>(lapack_dgeev));  
}

