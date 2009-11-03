//ff-c++-LIBRARY-dep:   lapack

#include "ff++.hpp"
//ff-c++-LIBRARY-dep:   lapack
//ff-c++-LIBRARY-dep:   blas

extern "C" {
typedef int integer;
typedef Complex doublecomplex ;

/* Subroutine */ int zgesv_(integer *n, integer *nrhs, doublecomplex *a, 
			    integer *lda, integer *ipiv, doublecomplex *b, integer *ldb, integer *
			    info);

/* Subroutine */ int dgesv_(integer *n, integer *nrhs, double *a, integer *
			    lda, integer *ipiv, double *b, integer *ldb, integer *info);

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
 
class Init { public:
  Init();
};


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

