//ff-c++-LIBRARY-dep:   lapack
//ff-c++-LIBRARY-dep:   blas
#include "ff++.hpp"
#include "RNM.hpp"
#include "AFunction_ext.hpp" // Extension of "AFunction.hpp" to deal with more than 3 parameters function

using namespace std;


#ifdef __LP64__
  typedef int intblas;
  typedef int integer;
#else
  typedef long intblas;
  typedef long integer;
#endif


typedef integer  logical;
typedef float   LAPACK_real;
typedef double   doublereal;
typedef logical  (* L_fp)();
typedef integer      ftnlen;

typedef complex<float> LAPACK_complex;
typedef complex<double> doublecomplex;
typedef void VOID; 
#define complex LAPACK_complex 
#define real LAPACK_real 

#include "clapack.h"
#undef real
#undef complex 
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
  // cout << lw << endl;
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
      nvp=0;
      (*vp)=Complex();
      (*vectp)=Complex();
    }
  return nvp;
}

long lapack_zgeev(KNM<Complex> *const &A,KN<Complex> *const &vp,KNM<Complex> *const &vectp)
{
  intblas nvp =0,zero=0;
  intblas n= A->N();
  ffassert(A->M()==n);
  ffassert(vectp->M()>=n);
  ffassert(vectp->N()>=n);
  ffassert(vp->N()>=n);
  KN<Complex> w(n),vr(n*n),vl(n*n);
  KNM<Complex> mat(*A);
  intblas info,lw=n*(n+1)*10;
  KN<Complex> wk(lw);
  KN<double> rwk(2*n);

  char N='N',V='V';
  // lw=1;// to get opt size value 
  zgeev_(&N,&V,&n, mat,&n, w, vl,&n, vr,&n,wk,&lw,rwk,&info);
  //  cout << lw << " " << wk[0] << " " << info <<   endl;
  /* lw=wk[0].real();
  w.resize(lw);
  zgeev_(&N,&V,&n, mat,&n, w, vl,&n, vr,&n,wk,&lw,rwk,&info);
  */
  if(info)
    cout << " info =  " << info << endl;
  if(!info)
    {
      int k=0;
      for(int i=0;i<n;++i)
        {
          (*vp)[i]=w[i];
          if(verbosity>2)
            cout << "   zgeev: vp "<< i << " : "  << (*vp)[i] << endl;
	  for(int j=0;j<n;++j)
              (*vectp)(j,i)=vr[k++];
          if(verbosity>5)
            cout << "   zgeev :   " << (*vectp)(':',i) <<endl;
        }
    }
  else
    {
      nvp=0;
      (*vp)=Complex();
      (*vectp)=Complex();
    }
  return nvp;
}

// VL, 10/02/2010
long lapack_dggev(KNM<double> *const &A,KNM<double> *const &B,KN<Complex> *const &vpa,KN<double> *const &vpb,KNM<Complex> *const &vectp)
{
    intblas nvp =0,zero=0;
    intblas n= A->N();
    ffassert(A->M()==n);
    ffassert(B->M()==n);
    ffassert(B->N()==n);
    ffassert(vectp->M()>=n);
    ffassert(vectp->N()>=n);
    ffassert(vpa->N()>=n);
    ffassert(vpb->N()>=n);
    
    KN<double> war(n),wai(n),wb(n),vr(n*n),vl(n*n);
    KNM<double> matA(*A);
    KNM<double> matB(*B);
    intblas info,lw=-1;  
    KN<double> w(1);
    //char N='N',V='V'; VL: do not compute eigenvectors (if yes, switch with following line)
    char VL='N',VR='N';
    
    dggev_(&VL,&VR,&n,matA,&n,matB,&n,war,wai,wb,vl,&n,vr,&n,w,&lw,&info);
    lw=w[0];
    // cout << lw << endl;
    w.resize(lw);
    dggev_(&VL,&VR,&n,matA,&n,matB,&n,war,wai,wb,vl,&n,vr,&n,w,&lw,&info);
    if(info)
	cout << " info =  " << info << endl;
    if(!info)
      {
	int k=0;
	for(int i=0;i<n;++i)
	  {
	    (*vpa)[i]=Complex(war[i],wai[i]);
	    (*vpb)[i]=wb[i];
	    if(verbosity>2)
		cout << "   dggev: vp "<< i << " : "  << (*vpa)[i] << " ; " << (*vpb)[i] << endl;
	    if( wai[i] == 0)
		for(int j=0;j<n;++j)
		    (*vectp)(j,i)=vr[k++];
	    else if (  wai[i] >  0)
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
		cout << "   dggev :   " << (*vectp)(':',i) <<endl;
	  }
      }
    else
      {
	nvp=0;
	(*vpa)=Complex();
	(*vectp)=Complex();
      }
    return nvp;
}

template<class T>
class Inverse{ public:
  T  t;
  Inverse( T  v)
   : t(v) {}
  template<class TT> Inverse( TT  v) : t(v) {}  
  template<class TT> Inverse( TT * v) : t(*v) {}  
  operator const T & () const {return t;}
};


template<class K>
class OneBinaryOperatorRNM_inv : public OneOperator { public:  
    OneBinaryOperatorRNM_inv() 
      : OneOperator( atype< Inverse< KNM<K>* > >(),atype<KNM<K> *>(),atype<long>()) {}
  E_F0 * code(const basicAC_F0 & args) const 
  { Expression p=args[1];
    if ( ! p->EvaluableWithOutStack() ) 
      { 
	bool bb=p->EvaluableWithOutStack();
	cout << bb << " " <<  * p <<  endl;
	CompileError(" A^p, The p must be a constant == -1, sorry");}
    long pv = GetAny<long>((*p)(0));
    if (pv !=-1)   
      { char buf[100];
	sprintf(buf," A^%ld, The pow must be  == -1, sorry",pv);
	CompileError(buf);}     
    return  new E_F_F0<Inverse< KNM<K>* > ,KNM<K> *>(Build<Inverse< KNM<K>* > ,KNM<K> *>,t[0]->CastTo(args[0])); 
  }
};
/* 
class Init { public:
  Init();
};
*/

KNM<R>* Solve(KNM<R>* a,Inverse<KNM<R >*> b) 
{
  /*
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
   *  N       (input) INTEGER
   *          The number of linear equations, i.e., the order of the
   *          matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the N-by-N coefficient matrix A.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U; the unit diagonal elements of L are not stored.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (output) INTEGER array, dimension (N)
   *          The pivot indices that define the permutation matrix P;
   *          row i of the matrix was interchanged with row IPIV(i).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the N-by-NRHS matrix of right hand side matrix B.
   *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
   *                has been completed, but the factor U is exactly
   *                singular, so the solution could not be computed.
   *
   */
  typedef double R;
  integer info;
  KNM<R> B(*b);
  integer  n= B.N();
  KN<integer> p(n);
  ffassert(B.M()==n);
  a->resize(n,n);
  *a=0.;
  for(int i=0;i<n;++i)
    (*a)(i,i)=(R) 1;;

  dgesv_(&n,&n,B,&n,p,*a,&n,&info);
  if(info) cerr << " error:  dgesv_ "<< info << endl;
  return a;
}


KNM<Complex>* SolveC(KNM<Complex>* a,Inverse<KNM<Complex >*> b) 
{
  /*
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
   *  N       (input) INTEGER
   *          The number of linear equations, i.e., the order of the
   *          matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the N-by-N coefficient matrix A.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U; the unit diagonal elements of L are not stored.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (output) INTEGER array, dimension (N)
   *          The pivot indices that define the permutation matrix P;
   *          row i of the matrix was interchanged with row IPIV(i).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the N-by-NRHS matrix of right hand side matrix B.
   *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
   *                has been completed, but the factor U is exactly
   *                singular, so the solution could not be computed.
   *
   */
  typedef Complex R;
  integer info;
  KNM<R> B(*b);
  integer   n= B.N();
  KN<integer> p(n);
  ffassert(B.M()==n);
  a->resize(n,n);
  *a=0.;
  for(int i=0;i<n;++i)
    (*a)(i,i)=(R) 1;;

  zgesv_(&n,&n,(R*) B,&n,p, (R*) *a,&n,&info);
  if(info) cerr << " error:  zgesv_ "<< info << endl;
  return a;
}
/*
static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  // avec de matrice plein 
  //  B = A ^-1 * C;
  //  A = A ^-1;
  Dcl_Type< Inverse<KNM<double >* > > ();
  Dcl_Type< Inverse<KNM<Complex >* > > ();
  
  TheOperators->Add("^", new OneBinaryOperatorRNM_inv<double>());
  TheOperators->Add("^", new OneBinaryOperatorRNM_inv<Complex>());
  TheOperators->Add("=", new OneOperator2<KNM<double>*,KNM<double>*,Inverse<KNM<double >*> >( Solve) );
  TheOperators->Add("=", new OneOperator2<KNM<Complex>*,KNM<Complex>*,Inverse<KNM<Complex >*> >( SolveC) );
  

 
}

*/
static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 

  if( map_type.find(typeid(Inverse<KNM<double >* >).name() ) == map_type.end() )
    {
      if(verbosity) 
	cout << " Add lapack interface ..." ;
      Dcl_Type< Inverse<KNM<double >* > > ();
      Dcl_Type< Inverse<KNM<Complex >* > > ();
      
      TheOperators->Add("^", new OneBinaryOperatorRNM_inv<double>());
      TheOperators->Add("^", new OneBinaryOperatorRNM_inv<Complex>());
      TheOperators->Add("=", new OneOperator2<KNM<double>*,KNM<double>*,Inverse<KNM<double >*> >( Solve) );
      TheOperators->Add("=", new OneOperator2<KNM<Complex>*,KNM<Complex>*,Inverse<KNM<Complex >*> >( SolveC) );
      
      Global.Add("inv","(",new  OneOperator1<bool,KNM<double>*>(lapack_inv));  
      Global.Add("dgeev","(",new  OneOperator3_<long,KNM<double>*,KN<Complex>*,KNM<Complex>*>(lapack_dgeev));  
      Global.Add("zgeev","(",new  OneOperator3_<long,KNM<Complex>*,KN<Complex>*,KNM<Complex>*>(lapack_zgeev));  
      Global.Add("dggev","(",new  OneOperator5_<long,KNM<double>*,KNM<double>*,KN<Complex>*,KN<double>*,KNM<Complex>*>(lapack_dggev));  

    }
 else
   cout << "( load: lapack <=> fflapack , skeep ) ";
}

