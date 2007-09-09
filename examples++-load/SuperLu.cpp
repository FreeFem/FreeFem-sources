#include  <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"

//#include "lex.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "slu_ddefs.h"
#include "slu_zdefs.h"

template <class R> struct SuperLUDriver
{
    
};

template <> struct SuperLUDriver<double>
{
    /* Driver routines */
    static  Dtype_t R_SLU_T() { return SLU_D;} 
    static void
    gssv(superlu_options_t * p1, SuperMatrix * p2, int * p3, int * p4, SuperMatrix * p5,
	  SuperMatrix * p6, SuperMatrix * p7 , SuperLUStat_t * p8, int * p9)
    { dgssv( p1,p2,p3,p4,p5,p6,p7,p8,p9); }
    
    
    static void
	gssvx(superlu_options_t * p1, SuperMatrix * p2, int * p3, int * p4, int * p5,
	       char * p6, double * p7, double * p8, SuperMatrix * p9, SuperMatrix * p10,
	       void * p11, int p12, SuperMatrix * p13, SuperMatrix * p14,
	       double * p15, double * p16, double * p17, double * p18,
	       mem_usage_t * p19, SuperLUStat_t * p20, int * p21)
    { dgssvx( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,  p11,p12,p13,p14,p15,p16,p17,p18,p19,p20, p21); }
    
    
    
    /* Supernodal LU factor related */
    static void
	Create_CompCol_Matrix(SuperMatrix * p1, int p2 , int p3, int p4, double * p5,
			       int * p6, int * p7, Stype_t p8, Dtype_t p9 , Mtype_t p10)
    {
	    dCreate_CompCol_Matrix( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10);
    }
    
    
    static void
	Create_CompRow_Matrix(SuperMatrix * p1, int p2, int p3, int p4, double * p5,
			       int * p6, int * p7, Stype_t p8, Dtype_t p9, Mtype_t p10)
    {
	dCreate_CompRow_Matrix( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10);
    }
    
    
    static void
	Create_Dense_Matrix(SuperMatrix * p1, int p2, int p3, double * p4, int p5,
			     Stype_t p6, Dtype_t p7, Mtype_t p8)
    {
	dCreate_Dense_Matrix( p1,p2,p3,p4,p5,p6,p7,p8);
    }
    
    
    static void
	Create_SuperNode_Matrix(SuperMatrix * p1, int p2, int p3, int p4, double * p5, 
				 int * p6, int * p7, int * p8, int * p9, int * p10,
				 Stype_t p11, Dtype_t p12, Mtype_t p13)
    {
	    dCreate_SuperNode_Matrix( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,  p11,p12,p13);
    }
    
    
};



template <> struct SuperLUDriver<Complex>
{
    /* Driver routines */
        static  Dtype_t R_SLU_T() { return SLU_Z;} 
    static doublecomplex *dc(Complex *p)  { return (doublecomplex *) (void *) p;}
    
    static void
    gssv(superlu_options_t * p1, SuperMatrix * p2, int * p3, int * p4, SuperMatrix * p5,
	 SuperMatrix * p6, SuperMatrix * p7 , SuperLUStat_t * p8, int * p9)
    { zgssv( p1,p2,p3,p4,p5,p6,p7,p8,p9); }
    
    
    static void
    gssvx(superlu_options_t * p1, SuperMatrix * p2, int * p3, int * p4, int * p5,
	  char * p6, double * p7, double * p8, SuperMatrix * p9, SuperMatrix * p10,
	  void * p11, int p12, SuperMatrix * p13, SuperMatrix * p14,
	  double * p15, double * p16, double * p17, double * p18,
	  mem_usage_t * p19, SuperLUStat_t * p20, int * p21)
    { zgssvx( p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,  p11,p12,p13,p14,p15,p16,p17,p18,p19,p20, p21); }
    
    
    
    /* Supernodal LU factor related */
    static void
    Create_CompCol_Matrix(SuperMatrix * p1, int p2 , int p3, int p4, Complex * p5,
			  int * p6, int * p7, Stype_t p8, Dtype_t p9 , Mtype_t p10)
    {
	zCreate_CompCol_Matrix( p1,p2,p3,p4,dc(p5),p6,p7,p8,p9,p10);
    }
    
    
    static void
    Create_CompRow_Matrix(SuperMatrix * p1, int p2, int p3, int p4, Complex * p5,
			  int * p6, int * p7, Stype_t p8, Dtype_t p9, Mtype_t p10)
    {
	zCreate_CompRow_Matrix( p1,p2,p3,p4,dc(p5),p6,p7,p8,p9,p10);
    }
    
    
    static void
    Create_Dense_Matrix(SuperMatrix * p1, int p2, int p3, Complex * p4, int p5,
			Stype_t p6, Dtype_t p7, Mtype_t p8)
    {
	zCreate_Dense_Matrix( p1,p2,p3,dc(p4),p5,p6,p7,p8);
    }
    
    
    static void
    Create_SuperNode_Matrix(SuperMatrix * p1, int p2, int p3, int p4, Complex * p5, 
			    int * p6, int * p7, int * p8, int * p9, int * p10,
			    Stype_t p11, Dtype_t p12, Mtype_t p13)
    {
	zCreate_SuperNode_Matrix( p1,p2,p3,p4,dc(p5),p6,p7,p8,p9,p10,  p11,p12,p13);
    }
    
    
};


template<class R>
class SolveSuperLU :   public MatriceMorse<R>::VirtualSolver, public SuperLUDriver<R>   {
  double eps;
  mutable double  epsr;
  double tgv;
   double tol_pivot_sym,tol_pivot; //Add 31 oct 2005
   
   
   mutable char           equed[1];
   yes_no_t       equil;
   mutable SuperMatrix    A, L, U;
   NCformat       *Astore;
   NCformat       *Ustore;
   SCformat       *Lstore;
   R         *a;
   int            *asub, *xa;
   int            *perm_c; /* column permutation vector */
   int            *perm_r; /* row permutations from partial pivoting */
   int            *etree;
   R         *rhsb, *rhsx, *xact;
   double         *RR, *CC;
   int m, n, nnz;
   
   mutable superlu_options_t options;
   mutable mem_usage_t    mem_usage;
   
public:
  SolveSuperLU(const MatriceMorse<R> &AA,int strategy,double ttgv, double epsilon=1e-6,
	       double pivot=-1.,double pivot_sym=-1.  ) : 
    eps(epsilon),epsr(0),
    tgv(ttgv),
    etree(0),perm_r(0)  ,perm_c(0),
    RR(0), CC(0),

    tol_pivot_sym(pivot_sym),tol_pivot(pivot)
  { 
     SuperMatrix    B, X;
     SuperLUStat_t stat;
     void           *work=0;
     int            info, lwork=0, nrhs=0;
     int            i;
     double         ferr, berr;
     double          rpg, rcond;
    
	
     A.Store=0;
     B.Store=0;
     X.Store=0;
     L.Store=0;
     U.Store=0;
	
    int status;
     n=AA.n;
     m=AA.m;
     nnz=AA.nbcoef;
     a=AA.a;
     asub=AA.cl;
     xa=AA.lg;
	 
    /* Defaults */
    lwork = 0;
    nrhs  = 0;
    
    /* Set the default values for options argument:
	options.Fact = DOFACT;
    options.Equil = YES;
    options.ColPerm = COLAMD;
    options.DiagPivotThresh = 1.0;
    options.Trans = NOTRANS;
    options.IterRefine = NOREFINE;
    options.SymmetricMode = NO;
    options.PivotGrowth = NO;
    options.ConditionNumber = NO;
    options.PrintStat = YES;
    */
    set_default_options(&options);
    options.Trans = TRANS;

    Dtype_t R_SLU = SuperLUDriver<R>::R_SLU_T(); 
    Create_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, R_SLU, SLU_GE);
    if(verbosity>4)
    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, nnz);
    
    Create_Dense_Matrix(&B, m, 0, (R*) 0, m, SLU_DN, R_SLU, SLU_GE);
    Create_Dense_Matrix(&X, m, 0, (R*) 0, m, SLU_DN, R_SLU, SLU_GE);
    
    if ( !(etree  = new int[n]) ) ABORT("Malloc fails for etree[].");
    if ( !(perm_r = new int[n]) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = new int[n]) ) ABORT("Malloc fails for perm_c[].");
    
    if ( !(RR = new double[n]) ) 
        ABORT("SUPERLU_MALLOC fails for R[].");
    if ( !(CC = new double[m]) )
        ABORT("SUPERLU_MALLOC fails for C[].");
    ferr=0;
    berr=0;
    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    /* ONLY PERFORM THE LU DECOMPOSITION */
    B.ncol = 0;  /* Indicate not to solve the system */
    SuperLUDriver<R>::gssvx(&options, &A, perm_c, perm_r, etree, equed, RR, CC,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, &ferr, &berr,
           &mem_usage, &stat, &info);
    if(verbosity>2)
    printf("LU factorization: dgssvx() returns info %d\n", info);
    if(verbosity>3)
    {
    if ( info == 0 || info == n+1 ) {
	
	if ( options.PivotGrowth ) printf("Recip. pivot growth = %e\n", rpg);
	if ( options.ConditionNumber )
	    printf("Recip. condition number = %e\n", rcond);
        Lstore = (SCformat *) L.Store;
        Ustore = (NCformat *) U.Store;
	printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
    	printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
    	printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
	printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	       mem_usage.for_lu/1e6, mem_usage.total_needed/1e6,
	       mem_usage.expansions);
	fflush(stdout);
	
    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
    }
    }
    if ( verbosity>5 ) StatPrint(&stat);
    StatFree(&stat);
    if( B.Store)  Destroy_SuperMatrix_Store(&B);
    if( X.Store)  Destroy_SuperMatrix_Store(&X);
    options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */


  }
  void Solver(const MatriceMorse<R> &AA,KN_<R> &x,const KN_<R> &b) const  {
    SuperMatrix    B, X;
    SuperLUStat_t stat;
    void           *work=0;
    int            info=0, lwork=0, nrhs=1;
    int            i;
    double         ferr, berr;
    double         rpg, rcond;
    
    B.Store=0;
    X.Store=0;
    ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    Dtype_t R_SLU = SuperLUDriver<R>::R_SLU_T(); 

    Create_Dense_Matrix(&B, m, 1, b, m, SLU_DN, R_SLU, SLU_GE);
    Create_Dense_Matrix(&X, m, 1, x, m, SLU_DN, R_SLU, SLU_GE);
    
    B.ncol = nrhs;  /* Set the number of right-hand side */
    
    /* Initialize the statistics variables. */
    StatInit(&stat);
    
    SuperLUDriver<R>::gssvx(&options, &A, perm_c, perm_r, etree, equed, RR, CC,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, &ferr, &berr,
           &mem_usage, &stat, &info);
    if(verbosity>2)
    printf("Triangular solve: dgssvx() returns info %d\n", info);
    
    if(verbosity>3)
    {
    if ( info == 0 || info == n+1 ) {
	
        /* This is how you could access the solution matrix. */
        R *sol = (R*) ((DNformat*) X.Store)->nzval; 
	
	if ( options.IterRefine ) {
            printf("Iterative Refinement:\n");
	    printf("%8s%8s%16s%16s\n", "rhs", "Steps", "FERR", "BERR");
	    printf("%8d%8d%16e%16e\n", i+1, stat.RefineSteps, ferr, berr);
	}
	fflush(stdout);
    } else if ( info > 0 && lwork == -1 ) {
        printf("** Estimated memory: %d bytes\n", info - n);
    }
    }
    if(verbosity>1) cout << "   x min max " << x.min() << " " <<x.max() << endl;
    if( B.Store)  Destroy_SuperMatrix_Store(&B);
    if( X.Store)  Destroy_SuperMatrix_Store(&X);
  }

  ~SolveSuperLU() { 
   if(verbosity>3)
       cout << "~SolveSuperLU S:" << endl;
      if (etree)    delete[] etree; 
      if (perm_r)   delete[] perm_r; 
      if (perm_c)   delete[] perm_c; 
      if (RR)   delete[] RR; 
      if (CC)   delete[] CC; 
      if( A.Store)  Destroy_SuperMatrix_Store(&A);
      if( L.Store)  Destroy_SuperNode_Matrix(&L);
      if( U.Store)  Destroy_CompCol_Matrix(&U);
      
  }
  void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<R> &) (*this) * x; 
  }
     
}; 




 MatriceMorse<double>::VirtualSolver *
BuildSolverSuperLU(const MatriceMorse<double> *A,int strategy,double tgv, double eps, double tol_pivot,double tol_pivot_sym,
		   int NbSpace,int itmax ,const  void * precon, void * stack)
{
    if(verbosity>9)
    cout << " BuildSolverSuperLU<double>" << endl;
    return new SolveSuperLU<double>(*A,strategy,tgv,eps,tol_pivot,tol_pivot_sym);
}


 MatriceMorse<Complex>::VirtualSolver *
BuildSolverSuperLU(const MatriceMorse<Complex> *A,int strategy,double tgv, double eps, double tol_pivot,double tol_pivot_sym,
		   int NbSpace,int itmax ,const  void * precon, void * stack)
{
    if(verbosity>9)
    cout << " BuildSolverSuperLU<Complex>" << endl;
    return new SolveSuperLU<Complex>(*A,strategy,tgv,eps,tol_pivot,tol_pivot_sym);
}


class Init { public:
    Init();
};

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

bool SetSuperLU()
{
    if(verbosity>1)
	cout << " SetDefault sparse solver to SuperLU" << endl;
    DefSparseSolver<double>::solver  =BuildSolverSuperLU;
    DefSparseSolver<Complex>::solver =BuildSolverSuperLU;    
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
}



Init init;
Init::Init()
{ 
  
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  
  if(verbosity>1)
    cout << "\n Add: SuperLU,  defaultsolver defaultsolverSuperLU" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  DefSparseSolver<double>::solver =BuildSolverSuperLU;
  DefSparseSolver<Complex>::solver =BuildSolverSuperLU;
  if(! Global.Find("defaultsolver").NotNull() )
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  Global.Add("defaulttoSuperLU","(",new OneOperator0<bool>(SetSuperLU));
}

