//  file to add UMFPACK solver with dynamic load.
//ff-c++-LIBRARY-dep:  amd  umfpack blas 
//ff-c++-cpp-dep: 

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
class SolveUMFPACK :   public MatriceMorse<R>::VirtualSolver  {
  double eps;
  mutable double  epsr;
  double tgv;
  void *Symbolic, *Numeric ;
  int umfpackstrategy;
  double tol_pivot_sym,tol_pivot; //Add 31 oct 2005
public:
  SolveUMFPACK(const MatriceMorse<R> &A,int strategy,double ttgv, double epsilon=1e-6,
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
    
    umfpack_di_defaults (Control) ;
    Control[UMFPACK_PRL]=1;
   // Control[UMFPACK_PIVOT_TOLERANCE]=1E-10;
    
    if(verbosity>4) Control[UMFPACK_PRL]=2;
    if(tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=pivot_sym;
    if(tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=pivot;
    if(umfpackstrategy>=0)   Control[UMFPACK_STRATEGY]=umfpackstrategy;
    if(verbosity>3) { 
      cout << "  UMFPACK real  Solver Control :" ;
      cout << "\n\t SYM_PIVOT_TOLERANCE "<< Control[UMFPACK_SYM_PIVOT_TOLERANCE];
      cout << "\n\t PIVOT_TOLERANCE     "<< Control[UMFPACK_PIVOT_TOLERANCE];
      cout << "\n\t PRL                 "<< Control[UMFPACK_PRL];
      cout << "\n";      
    }
    
    status = umfpack_di_symbolic (n, n, A.lg, A.cl, A.a, &Symbolic,Control,Info) ;
    if (status !=  0)
    {
      (void) umfpack_di_report_matrix (n, n, A.lg, A.cl, A.a, 1, Control) ;

	umfpack_di_report_info (Control, Info) ;
	umfpack_di_report_status (Control, status) ;
	cerr << "umfpack_di_symbolic failed" << endl;
	ExecError("umfpack_di_symbolic failed");
	//ffassert(0);
    }

    status = umfpack_di_numeric (A.lg, A.cl, A.a, Symbolic, &Numeric,Control,Info) ;
    if (status !=  0)
    {
	umfpack_di_report_info (Control, Info) ;
	umfpack_di_report_status (Control, status) ;
	cerr << "umfpack_di_numeric failed" << endl;
	ExecError("umfpack_di_numeric failed");
	ffassert(0);
    }

    if (Symbolic) umfpack_di_free_symbolic (&Symbolic),Symbolic=0; 
    if(verbosity>3)
    cout << "  -- umfpack_di_build LU " << n <<  endl;
    if(verbosity>5)     (void)  umfpack_di_report_info(Control,Info);

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
    
     umfpack_di_defaults (Control) ;
     // change UMFPACK_At to UMFPACK_Aat in complex 
     int status = umfpack_di_solve (UMFPACK_Aat, A.lg, A.cl, A.a,KN_2Ptr<R> (x), KN_2Ptr<R>(b), Numeric,Control,Info) ;
    if (status != 0)
    {
	umfpack_di_report_info (Control, Info) ;
	umfpack_di_report_status (Control, status) ;
	cerr << "umfpack_di_solve failed" << endl;
	ExecError("umfpack_di_solve failed");
	
	ffassert(0);
    }
     if(verbosity>2)
    cout << " -- umfpack_di_solve " << endl;
    if(verbosity>3)
    cout << "   b min max " << b.min() << " " <<b.max() << endl;
    if(verbosity>3)     (void)  umfpack_di_report_info(Control,Info);
     if(verbosity>1) cout << "   x min max " << x.min() << " " <<x.max() << endl;
  }

  ~SolveUMFPACK() { 
   if(verbosity>3)
    cout << "~SolveUMFPACK S:" << Symbolic << " N:" << Numeric <<endl;
    if (Symbolic)   umfpack_di_free_symbolic  (&Symbolic),Symbolic=0; 
    if (Numeric)    umfpack_di_free_numeric (&Numeric),Numeric=0;
  }
  void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<R> &) (*this) * x; 
  }
     
}; 


template<>
class SolveUMFPACK<Complex> :   public MatriceMorse<Complex>::VirtualSolver  {
  double eps;
  mutable double  epsr;
  int umfpackstrategy;
  double tgv;
  void *Symbolic, *Numeric ;
  double *ar,*ai;


    double tol_pivot_sym,tol_pivot; //Add 31 oct 2005

public:
  SolveUMFPACK(const MatriceMorse<Complex> &A,int strategy,double ttgv, double epsilon=1e-6,
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
    umfpack_zi_defaults (Control) ;
    Control[UMFPACK_PRL]=1;
    if(verbosity>4) Control[UMFPACK_PRL]=2;
   //    Control[UMFPACK_SYM_PIVOT_TOLERANCE]=1E-10;
  //  Control[UMFPACK_PIVOT_TOLERANCE]=1E-10;
    if(tol_pivot_sym>0) Control[UMFPACK_SYM_PIVOT_TOLERANCE]=pivot_sym;
    if(tol_pivot>0) Control[UMFPACK_PIVOT_TOLERANCE]=pivot;
    if(umfpackstrategy>=0) Control[UMFPACK_STRATEGY]=umfpackstrategy;
    if(verbosity>3) { 
      cout << "  UMFPACK complex Solver Control :" ;
      cout << "\n\t SYM_PIVOT_TOLERANCE "<< Control[UMFPACK_SYM_PIVOT_TOLERANCE];
      cout << "\n\t PIVOT_TOLERANCE     "<< Control[UMFPACK_PIVOT_TOLERANCE];
      cout << "\n\t PRL                 "<< Control[UMFPACK_PRL];
      cout << "\n";      
    }
    status = umfpack_zi_symbolic (n, n, A.lg, A.cl, ar,ai, &Symbolic,Control,Info) ;
    if (status < 0)
    {
      (void) umfpack_zi_report_matrix (n, n, A.lg, A.cl, ar,ai, 1, Control) ;

	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	cerr << "umfpack_zi_symbolic failed" << endl;
	ExecError("umfpack_zi_symbolic failed");
	ffassert(0);
	exit(2);
    }

    status = umfpack_zi_numeric (A.lg, A.cl, ar,ai, Symbolic, &Numeric,Control,Info) ;
    if (status < 0)
    {
	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	cerr << "umfpack_zi_numeric failed" << endl;
	ExecError("umfpack_zi_numeric failed");
	ffassert(0);
	exit(2);
    }

    if (Symbolic) umfpack_zi_free_symbolic (&Symbolic),Symbolic=0; 
    if(verbosity>3)
    cout << "umfpack_zi_build LU " << n <<  endl;
    if(verbosity>5)     (void)  umfpack_zi_report_info(Control,Info);

  }
  void Solver(const MatriceMorse<Complex> &A,KN_<Complex> &x,const KN_<Complex> &b) const  {
        ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    // cout << " epsr = " << epsr << endl;
    double Control[UMFPACK_CONTROL];
    double Info[UMFPACK_INFO];
     umfpack_zi_defaults (Control) ;
     int n = b.N();
     ffassert(A.ChecknbLine( n) && n == x.N() && A.ChecknbColumn(n) );
     KN<double> xr(n),xi(n),br(n),bi(n);
     C2RR(n,b,br,bi);
     // change UMFPACK_At to UMFPACK_Aat in complex  oct 2005 
    int status = umfpack_zi_solve (UMFPACK_Aat, A.lg, A.cl, ar,ai, xr, xi, br,bi, Numeric,Control,Info) ;
    if (status < 0)
    {
	umfpack_zi_report_info (Control, Info) ;
	umfpack_zi_report_status (Control, status) ;
	cerr << "umfpack_zi_solve failed" << endl;
	ExecError("umfpack_zi_numeric failed");
	ffassert(0);
	exit(2);
    }
    RR2C(n,xr,xi,x);
    if(verbosity>1)
    {
     cout << "  -- umfpack_zi_solve " << endl;
     if(verbosity>3)     (void)  umfpack_zi_report_info(Control,Info);
    
      cout << "   b min max " << b.min() << " " <<b.max() << endl;
      cout << "   x min max " << x.min() << " " <<x.max() << endl;
    }
  }

  ~SolveUMFPACK() { 
    if(verbosity>5)
    cout << "~SolveUMFPACK " << endl;
    if (Symbolic)   umfpack_zi_free_symbolic  (&Symbolic),Symbolic=0; 
    if (Numeric)    umfpack_zi_free_numeric (&Numeric),Numeric=0;
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
BuildSolverIUMFPack(DCL_ARG_SPARSE_SOLVER(double,A))
{
  if( verbosity>9)
    cout << " BuildSolverUMFPack<double>" << endl;
    return new SolveUMFPACK<double>(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym);
}

inline MatriceMorse<Complex>::VirtualSolver *
BuildSolverIUMFPack(DCL_ARG_SPARSE_SOLVER(Complex,A))
{
  if( verbosity>9)
    cout << " BuildSolverUMFPack<Complex>" << endl;
    return new SolveUMFPACK<Complex>(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym);
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

bool SetUMFPACK()
{
    if(verbosity>1)
	cout << " SetDefault sparse solver to IUMFPack" << endl;
    DefSparseSolver<double>::solver  =BuildSolverIUMFPack;
    DefSparseSolver<Complex>::solver =BuildSolverIUMFPack;    
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
}

class Init { public:
    Init();
};
LOADINIT(Init);
Init::Init(){    
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  if(verbosity>1)
    cout << "\n Add: UMFPACK:  defaultsolver defaultsolverUMFPACK" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  
  DefSparseSolver<double>::solver =BuildSolverIUMFPack;
  DefSparseSolver<Complex>::solver =BuildSolverIUMFPack;
  if(! Global.Find("defaultsolver").NotNull() )
    {    cout << "\n add defaultsolver" << endl;
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  }
  if(! Global.Find("defaulttoUMFPACK").NotNull() )
    Global.Add("defaulttoUMFPACK","(",new OneOperator0<bool>(SetUMFPACK));  
}


