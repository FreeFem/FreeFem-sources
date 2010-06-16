//   for automatic  compilation with ff-c++
//ff-c++-LIBRARY-dep:  umfpack amd blas 
//ff-c++-cpp-dep: 
//  
//  file to add UMFPACK solver with dynamic load.
#include  <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"


#include "MatriceCreuse_tpl.hpp"

 
#ifdef HAVE_LIBUMFPACK
extern "C" {
#ifdef HAVE_UMFPACK_H
#include <umfpack.h>
#else
#ifdef HAVE_UMFPACK_UMFPACK_H
#include <umfpack/umfpack.h>
#else
#ifdef HAVE_BIG_UMFPACK_UMFPACK_H
#include <UMFPACK/umfpack.h>
#else
#ifdef HAVE_UFSPARSE_UMFPACK_H
#include <ufsparse/umfpack.h>
#else
#ifdef HAVE_SUITESPARSE_UMFPACK_H
#include <suitesparse/umfpack.h>
#else

  // Defaults to a local version of the UMFPACK headers
#include "../../download/include/umfpack.h"

#endif // HAVE_SUITESPARSE_UMFPACK_H
#endif // HAVE_UFSPARSE_UMFPACK_H
#endif // HAVE_BIG_UMFPACK_UMFPACK_H
#endif // HAVE_UMFPACK_UMFPACK_H
#endif // HAVE_UMFPACK_H
}
#endif
template<class R>
class SolveUMFPACK64 :   public MatriceMorse<R>::VirtualSolver  {
  double eps;
  mutable double  epsr;
  double tgv;
  void *Symbolic, *Numeric ;
  int umfpackstrategy;
  double tol_pivot_sym,tol_pivot; //Add 31 oct 2005
public:
  SolveUMFPACK64(const MatriceMorse<R> &A,int strategy,double ttgv, double epsilon=1e-6,
	       double pivot=-1.,double pivot_sym=-1.  ) : 
    eps(epsilon),epsr(0),
    tgv(ttgv),
    Symbolic(0),Numeric(0)  ,
    umfpackstrategy(strategy),
    tol_pivot_sym(pivot_sym),tol_pivot(pivot)
  { 
    
    int status;
    throwassert( !A.sym() && Numeric == 0 && Symbolic==0 );
    int n=A.n;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    
    for(int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0;
    for(int i=0;i<UMFPACK_INFO;i++) Info[i]=0;
    
    umfpack_dl_defaults (Control) ;
    Control[UMFPACK_PRL]=1;
   // Control[UMFPACK_PIVOT_TOLERANCE]=1E-10;
    
    if(verbosity>4) Control[UMFPACK_PRL]=2;
    if(tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=pivot_sym;
    if(tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=pivot;
    if(umfpackstrategy>=0)   Control[UMFPACK_STRATEGY]=umfpackstrategy;
    if(verbosity>3) { 
      cout << "  UMFPACK (long) real  Solver Control :" ;
      cout << "\n\t SYM_PIVOT_TOLERANCE "<< Control[UMFPACK_SYM_PIVOT_TOLERANCE];
      cout << "\n\t PIVOT_TOLERANCE     "<< Control[UMFPACK_PIVOT_TOLERANCE];
      cout << "\n\t PRL                 "<< Control[UMFPACK_PRL];
      cout << "\n";      
    }
    //  convert   array in long ...
    KN<long> Alg(n+1),Acl(A.nbcoef);
    for(int i=0;i<=n;++i)
      Alg[i]=A.lg[i];

    for(int i=0;i<A.nbcoef;++i)
      Acl[i]=A.cl[i];

    status = umfpack_dl_symbolic (n, n, Alg, Acl, A.a, &Symbolic,Control,Info) ;
    if (status !=  0)
    {
      (void) umfpack_dl_report_matrix (n, n, Alg, Acl, A.a, 1, Control) ;

	umfpack_dl_report_info (Control, Info) ;
	umfpack_dl_report_status (Control, status) ;
	cerr << "umfpack_dl_symbolic failed" << endl;
	ExecError("umfpack_dl_symbolic failed");
	//ffassert(0);
    }

    status = umfpack_dl_numeric (Alg, Acl, A.a, Symbolic, &Numeric,Control,Info) ;
    if (status !=  0)
    {
	umfpack_dl_report_info (Control, Info) ;
	umfpack_dl_report_status (Control, status) ;
	cerr << "umfpack_dl_numeric failed" << endl;
	ExecError("umfpack_dl_numeric failed");
	ffassert(0);
    }

    if (Symbolic) umfpack_dl_free_symbolic (&Symbolic),Symbolic=0; 
    if(verbosity>3)
    cout << "  -- umfpack_dl_build LU " << n <<  endl;
    if(verbosity>5)     (void)  umfpack_dl_report_info(Control,Info);

  }
  void Solver(const MatriceMorse<R> &A,KN_<R> &x,const KN_<R> &b) const  {
    ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    for(int i=0;i<UMFPACK_CONTROL;i++) Control[i]=0;
    for(int i=0;i<UMFPACK_INFO;i++) Info[i]=0;
    int n= b.N(); 
     ffassert(A.ChecknbLine( n) && n == x.N() && A.ChecknbColumn(n) );
    
     umfpack_dl_defaults (Control) ;
    //  convert   array in long ...
    KN<long> Alg(n+1),Acl(A.nbcoef);
    for(int i=0;i<=n;++i)
      Alg[i]=A.lg[i];

    for(int i=0;i<A.nbcoef;++i)
      Acl[i]=A.cl[i];

     // change UMFPACK_At to UMFPACK_Aat in complex 
    int status = umfpack_dl_solve (UMFPACK_Aat, Alg, Acl, A.a, x, b, Numeric,Control,Info) ;
    if (status != 0)
    {
	umfpack_dl_report_info (Control, Info) ;
	umfpack_dl_report_status (Control, status) ;
	cerr << "umfpack_dl_solve failed" << endl;
	ExecError("umfpack_dl_solve failed");
	
	ffassert(0);
    }
     if(verbosity>2)
    cout << " -- umfpack_dl_solve " << endl;
    if(verbosity>3)
    cout << "   b min max " << b.min() << " " <<b.max() << endl;
    if(verbosity>3)     (void)  umfpack_dl_report_info(Control,Info);
     if(verbosity>1) cout << "   x min max " << x.min() << " " <<x.max() << endl;
  }

  ~SolveUMFPACK64() { 
   if(verbosity>3)
    cout << "~SolveUMFPACK 64:" << Symbolic << " N:" << Numeric <<endl;
    if (Symbolic)   umfpack_dl_free_symbolic  (&Symbolic),Symbolic=0; 
    if (Numeric)    umfpack_dl_free_numeric (&Numeric),Numeric=0;
  }
  void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<R> &) (*this) * x; 
  }
     
}; 


template<>
class SolveUMFPACK64<Complex> :   public MatriceMorse<Complex>::VirtualSolver  {
  double eps;
  mutable double  epsr;
  int umfpackstrategy;
  double tgv;
  void *Symbolic, *Numeric ;
  double *ar,*ai;


    double tol_pivot_sym,tol_pivot; //Add 31 oct 2005

public:
  SolveUMFPACK64(const MatriceMorse<Complex> &A,int strategy,double ttgv, double epsilon=1e-6,
     double pivot=-1.,double pivot_sym=-1.
) : 
    eps(epsilon),epsr(0),umfpackstrategy(strategy),tgv(ttgv),
    Symbolic(0),Numeric(0),
    ar(0),ai(0),
    tol_pivot_sym(pivot_sym), 
    tol_pivot(pivot)
   { 
    int status;
    throwassert( !A.sym());
    int n=A.n;
    //  copy the coef of the matrice ---
     ar= new double[A.nbcoef];
     ai= new double[A.nbcoef];
     ffassert(ar && ai);
     C2RR(A.nbcoef,A.a,ar,ai);
        
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
    umfpack_zl_defaults (Control) ;
    Control[UMFPACK_PRL]=1;
    if(verbosity>4) Control[UMFPACK_PRL]=2;
   //    Control[UMFPACK_SYM_PIVOT_TOLERANCE]=1E-10;
  //  Control[UMFPACK_PIVOT_TOLERANCE]=1E-10;
    if(tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=pivot_sym;
    if(tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=pivot;
    if(umfpackstrategy>=0) Control[UMFPACK_STRATEGY]=umfpackstrategy;
    if(verbosity>3) { 
      cout << "  UMFPACK(64) complex Solver Control :" ;
      cout << "\n\t SYM_PIVOT_TOLERANCE "<< Control[UMFPACK_SYM_PIVOT_TOLERANCE];
      cout << "\n\t PIVOT_TOLERANCE     "<< Control[UMFPACK_PIVOT_TOLERANCE];
      cout << "\n\t PRL                 "<< Control[UMFPACK_PRL];
      cout << "\n";      
    }

    //  convert   array in long ...
    KN<long> Alg(n+1),Acl(A.nbcoef);
    for(int i=0;i<=n;++i)
      Alg[i]=A.lg[i];

    for(int i=0;i<A.nbcoef;++i)
      Acl[i]=A.cl[i];


    status = umfpack_zl_symbolic (n, n, Alg, Acl, ar,ai, &Symbolic,Control,Info) ;
    if (status < 0)
    {
      (void) umfpack_zl_report_matrix (n, n, Alg, Acl, ar,ai, 1, Control) ;

	umfpack_zl_report_info (Control, Info) ;
	umfpack_zl_report_status (Control, status) ;
	cerr << "umfpack_zl_symbolic failed" << endl;
	ExecError("umfpack_zl_symbolic failed");
	ffassert(0);
	exit(2);
    }

    status = umfpack_zl_numeric (Alg, Acl, ar,ai, Symbolic, &Numeric,Control,Info) ;
    if (status < 0)
    {
	umfpack_zl_report_info (Control, Info) ;
	umfpack_zl_report_status (Control, status) ;
	cerr << "umfpack_zl_numeric failed" << endl;
	ExecError("umfpack_zl_numeric failed");
	ffassert(0);
	exit(2);
    }

    if (Symbolic) umfpack_zl_free_symbolic (&Symbolic),Symbolic=0; 
    if(verbosity>3)
    cout << "umfpack_zl_build LU " << n <<  endl;
    if(verbosity>5)     (void)  umfpack_zl_report_info(Control,Info);

  }
  void Solver(const MatriceMorse<Complex> &A,KN_<Complex> &x,const KN_<Complex> &b) const  {
        ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
     umfpack_zl_defaults (Control) ;
     int n = b.N();
     ffassert(A.ChecknbLine( n) && n == x.N() && A.ChecknbColumn(n) );
     KN<double> xr(n),xi(n),br(n),bi(n);
     C2RR(n,b,br,bi);
     // change UMFPACK_At to UMFPACK_Aat in complex  oct 2005

    //  convert   array in long ...
    KN<long> Alg(n+1),Acl(A.nbcoef);
    for(int i=0;i<=n;++i)
      Alg[i]=A.lg[i];

    for(int i=0;i<A.nbcoef;++i)
      Acl[i]=A.cl[i];

 
    int status = umfpack_zl_solve (UMFPACK_Aat, Alg, Acl, ar,ai, xr, xi, br,bi, Numeric,Control,Info) ;
    if (status < 0)
    {
	umfpack_zl_report_info (Control, Info) ;
	umfpack_zl_report_status (Control, status) ;
	cerr << "umfpack_zl_solve failed" << endl;
	ExecError("umfpack_zl_numeric failed");
	ffassert(0);
	exit(2);
    }
    RR2C(n,xr,xi,x);
    if(verbosity>1)
    {
     cout << "  -- umfpack_zl_solve " << endl;
     if(verbosity>3)     (void)  umfpack_zl_report_info(Control,Info);
    
      cout << "   b min max " << b.min() << " " <<b.max() << endl;
      cout << "   x min max " << x.min() << " " <<x.max() << endl;
    }
  }

  ~SolveUMFPACK64() { 
    if(verbosity>5)
    cout << "~SolveUMFPACK64 " << endl;
    if (Symbolic)   umfpack_zl_free_symbolic  (&Symbolic),Symbolic=0; 
    if (Numeric)    umfpack_zl_free_numeric (&Numeric),Numeric=0;
    delete [] ar;
    delete [] ai;   
  }
  void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<Complex> &) (*this) * x; 
  }
     

}; 

inline MatriceMorse<double>::VirtualSolver *
BuildSolverIUMFPack64(DCL_ARG_SPARSE_SOLVER(double,A))
{
    cout << " BuildSolverUMFPack64<double>" << endl;
    return new SolveUMFPACK64<double>(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym);
}

inline MatriceMorse<Complex>::VirtualSolver *
BuildSolverIUMFPack64(DCL_ARG_SPARSE_SOLVER(Complex,A))
{
    cout << " BuildSolverUMFPack64<Complex>" << endl;
    return new SolveUMFPACK64<Complex>(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym);
}


//  the 2 default sparse solver double and complex
DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; ;
DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetDefault()
{
    if(verbosity>1)
	cout << " SetDefault sparse to default" << endl;
    DefSparseSolver<double>::solver =SparseMatSolver_R;
    DefSparseSolver<Complex>::solver =SparseMatSolver_C;
    TypeSolveMat::defaultvalue =TypeSolveMat::SparseSolver;
}

bool SetUMFPACK64()
{
    if(verbosity>1)
	cout << " SetDefault sparse solver to IUMFPack64" << endl;
    DefSparseSolver<double>::solver  =BuildSolverIUMFPack64;
    DefSparseSolver<Complex>::solver =BuildSolverIUMFPack64;    
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
}

class Init { public:
    Init();
};
Init init;
Init::Init(){    
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  if(verbosity>1)
    cout << "\n Add: UMFPACK64:  defaultsolver defaultsolverUMFPACK64" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  
  DefSparseSolver<double>::solver =BuildSolverIUMFPack64;
  DefSparseSolver<Complex>::solver =BuildSolverIUMFPack64;
  if(! Global.Find("defaultsolver").NotNull() )
    {    cout << "\n add defaultsolver (64)" << endl;
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  }
  if(! Global.Find("defaulttoUMFPACK64").NotNull() )
    Global.Add("defaulttoUMFPACK64","(",new OneOperator0<bool>(SetUMFPACK64));  
}


