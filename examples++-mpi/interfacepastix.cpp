/*
  Interface entre freefem++ et pastix
*/
//ff-c++-LIBRARY-dep: double_pastix   blas parmetis metis scotch  mpi fc
//ff-c++-cpp-dep: 
#include <mpi.h>
#include  <iostream>
using namespace std;
   
#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"

//#include "lex.hpp"
#include "MatriceCreuse_tpl.hpp"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <time.h>
#include <sys/time.h>




// #include <ctype.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <unistd.h>
// //#include <pthread.h>
// #include <string.h>
// #include <time.h>
// #include <sys/time.h>
// #include "mpi.h"
   
// #include <assert.h>
// #include "pastix.h"
// #include "cscd_utils.h"
// #include "read_matrix.h"

#include <assert.h>
extern "C"{
#include "pastix.h"
}
//#include "cscd_utils.h"
//#include "read_matrix.h"

#undef memFree_null
#define memFree_null(x) {if (x ==NULL) {fprintf(stdout,"%s:%d freeing NULL\n",__FILE__,__LINE__);} free(x); x=NULL;}

#define STR_SIZE 256

static pastix_int_t * pastixint(int * ii){ return (pastix_int_t*) (void *) ii;} 
static pastix_float_t * pastixfloat(double * ii){ return (pastix_float_t*) (void *) ii;} 

typedef struct pastix_param {
  pastix_data_t          *pastix_data; /*Pointer used by PaStiX to keep information alive between calls */
  MPI_Comm                comm;        /* Communicator used by PaStiX                                    */
  pastix_int_t            Ncol;        /* Size of the Matrix                                             */
  pastix_int_t           *ia;          /* Index of first element of each column in ja and avals          */  
  pastix_int_t           *ja;          /* Rows of the unknows of the matrix                              */
  pastix_float_t         *avals;       /* Values of the matrix                                           */
  pastix_int_t           *perm;        /* Permutation used for re-numbering of the unknowns              */
  pastix_int_t           *invp;        /* Inverse permutation                                            */
  pastix_float_t         *rhs;         /* Right hand side                                                */
  pastix_int_t           *iparm;       /* Integer parameters                                             */
  double                 *dparm;       /* Floating parameters                                            */
} pastix_param_t;



void
Morse_to_CSC(int m, int n, int nnz, 
	     double *a, int *colind, int  *rowptr,
	     pastix_float_t **at, pastix_int_t **rowind, 
	     pastix_int_t **colptr)
{
    register int i, j, col, relpos;
    pastix_int_t *marker;

    /* Allocate storage for another copy of the matrix. */
    *at     = (pastix_float_t *) malloc(sizeof(pastix_float_t)*nnz);
    *rowind = (pastix_int_t *) malloc(sizeof(pastix_int_t)*nnz);
    *colptr = (pastix_int_t *) malloc(sizeof(pastix_int_t)*(n+1));
    marker  = (pastix_int_t *) malloc(sizeof(pastix_int_t)*n);
    
    for (i = 0; i < n; ++i)
      marker[i] = 0;
    /* Get counts of each column of A, and set up column pointers */
    for (i = 0; i < m; ++i)
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) ++marker[colind[j]];
    (*colptr)[0] = 0;
    for (j = 0; j < n; ++j) {
	(*colptr)[j+1] = (*colptr)[j] + marker[j];
	marker[j] = (*colptr)[j];
    }

    /* Transfer the matrix into the compressed column storage. */
    for (i = 0; i < m; ++i) {
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) {
	    col = colind[j];
	    relpos = marker[col];
	    (*rowind)[relpos] = i;
	    (*at)[relpos] = a[j];
	    ++marker[col];
	}
    }

    free(marker);
}

static const int MAX_CHAR_PER_LINE=256;
void read_datafile_pastixff(const string &datafile, pastix_int_t *iparmtab, double *dparmtab){
 
  FILE*   m_File;
  int     i = 0;
  char    szbuff[MAX_CHAR_PER_LINE];
  char*   token;

  char filename[datafile.size()+1];  
  strcpy( filename, datafile.c_str()); 

  m_File = fopen(filename,"rt");

  if(!m_File)
    {
      printf("error in reading filename %s\n",&filename);
    }

  fgets(szbuff,MAX_CHAR_PER_LINE,m_File);
  token = strtok(szbuff," /#!\t\n");
  
  if( !(strcmp(token,"iparm") == 0) ){
    printf("freefem++: error in reading iparm parameter for pastix (see strcuture of ffpastix_iparm_dparm.txt) \n");
    exit(1);
  }
  else{
    printf("freefem++: reading iparm parameter for pastix \n");    
  }
  while(!feof(m_File) && i < 64)
    {   
      fgets(szbuff,MAX_CHAR_PER_LINE,m_File);
      token = strtok(szbuff," /#!\t\n");
      iparmtab[i] = (pastix_int_t)atol(token);
      i++;
    }

  i=0;
  fgets(szbuff,MAX_CHAR_PER_LINE,m_File);
  token = strtok(szbuff," /#!\t\n");  
  if( !(strcmp(token,"dparm") == 0) ){
    printf("freefem++: error in reading dparm parameter for pastix (see strcuture of ffpastix_iparm_dparm.txt) \n");
    exit(1);
  }
  else{
    printf("freefem++: reading dparm parameter for pastix \n");    
  }
  while(!feof(m_File) && i < 64)
    {   
      fgets(szbuff,MAX_CHAR_PER_LINE,m_File);
      token = strtok(szbuff," /#!\t\n");
      dparmtab[i] = atof(token);
      i++;
    }
 
  fclose(m_File);

#ifdef OOC
/*   if (iparmtab[IPARM_OOC_THREAD] > 1) */
    iparmtab[IPARM_OOC_THREAD] = 1;
#endif
  /* On empeche le 2d avec NUMA_ALLOC */
#ifdef NUMA_ALLOC
  if (iparmtab[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
      errorPrint("2D not available with NUMA allocation\n");
      exit(-1);
    }
#endif
}

// ATTENTION :: pastix_float_t  
//      peut être soit un complex ou un reel cela depend de la maniere dont on a compiler pastix

// CAS DOUBLE SEULEMENT 


class dSolvepastixmpi :   public MatriceMorse<double>::VirtualSolver   {
  
  double eps;
  mutable double  epsr;
  double tgv;
  double tol_pivot_sym,tol_pivot; //Add 31 oct 2005
  

  int paraoption;
  int myid, mpi_size;
  int Nrow;
  int mpi_flag;
  int init_raff;
  int thrd_flag;
  int SYM;
  
  string data_option;
  
  mutable pastix_int_t    iparm[64];
  mutable double          dparm[64];
  mutable pastix_int_t    Ncol;
  mutable pastix_int_t   *ia;
  mutable pastix_int_t   *ja;
  mutable pastix_float_t *avals;
  mutable pastix_int_t   *loc2glob;
  //char           *Type    = NULL;
  //char           *RhsType = NULL;
  mutable pastix_float_t *rhs;
  mutable pastix_int_t   *perm;
  mutable pastix_int_t   *invp;
  mutable pastix_data_t  *pastix_data;
  

public:

  dSolvepastixmpi(const MatriceMorse<double> &AA, string datafile, KN<long> &param_int, KN<double> &param_double, 
		  KN<long> &pperm_r, KN<long> &pperm_c) : 
    data_option(datafile) 
  { 
    //int m;
    //int ierr;
    struct timeval  tv1, tv2;
   
    ia    = NULL;
    ja    = NULL;
    avals   = NULL;
    loc2glob = NULL;
    rhs     = NULL;
    pastix_data = NULL;
    
    // matrix assembled on host
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    printf("- Rang MPI : %d\n", myid);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // SYMETRIQUE
    mpi_flag  = 0;
    thrd_flag = 0;
        
    Ncol = AA.m;
    Nrow = AA.n;
  
    // Avant : on ecrit la transposée
    /*
      ia = (pastix_int_t *) malloc( (Ncol+1)*sizeof(pastix_int_t));
      for(int ii=0; ii < Ncol+1; ii++){
      ia[ii] = AA.lg[ii]+1;
      }
      assert( ia[Ncol]-1 == AA.nbcoef );
      
      ja = (pastix_int_t *) malloc((ia[Ncol]-1)*sizeof(pastix_int_t));
      for(int ii=0; ii < ia[Ncol]-1; ii++)
      ja[ii] = AA.cl[ii]+1;
      
      
      if( sizeof(pastix_float_t) == sizeof(double) ){
      avals = (pastix_float_t *) malloc( (ia[Ncol]-1)*sizeof(pastix_float_t));
      for(int ii=0; ii < ia[Ncol]-1; ii++)
      avals[ii] = AA.a[ii];
      }
    */    
    // AA.cl : indices des colonnes
    // AA.lg : pointeurs des lignes
    Morse_to_CSC( AA.n , AA.m, AA.nbcoef, AA.a, AA.cl, AA.lg, &avals, &ja, &ia);
    // ia : pointeurs des colonnes
    // ja : indices des lignes
    
    cout << "AA.n= "<< AA.n << " AA.m=" <<  AA.m << " AA.nbcoef=" << AA.nbcoef << endl;
    
    for(int ii=0; ii < Ncol+1; ii++){
      ia[ii] = ia[ii]+1;
    }
    assert( ia[Ncol]-1 == AA.nbcoef );
    for(int ii=0; ii < ia[Ncol]-1; ii++){
      ja[ii] = ja[ii]+1; 
    }
       
    perm = (pastix_int_t *) malloc(Ncol*sizeof(pastix_int_t));
    invp = (pastix_int_t *) malloc(Ncol*sizeof(pastix_int_t));
    
    rhs = (pastix_float_t *) malloc(Ncol*sizeof(pastix_float_t));
    
    // reading permutation given by the user
    if(pperm_r) 
      for(int ii=0; ii < Ncol; ii++)
	perm[ii] = pperm_r[ii];
    if(pperm_c)  
      for(int ii=0; ii < Ncol; ii++)
	invp[ii] = pperm_c[ii];
   

    // CAS DE LA  MATRICE NON DISTRIBUER
    pastix_int_t init_raff;
    fprintf(stdout,"-- INIT PARAMETERS --\n");
    
    

    // reading iparm from array    
    if(!data_option.empty()) read_datafile_pastixff(data_option,iparm,dparm);
    else if(param_int || param_double){
      if( param_int ) 
      {
	cout << "read param_int" << endl;
	assert(param_int.N() == 64);
	for(int ii=0; ii<64; ii++) 
	  iparm[ii] = param_int[ii];
	iparm[IPARM_MODIFY_PARAMETER] = API_YES;
      }
      if( param_double ) 
      {
	cout << "read param_double" << endl;
	assert(param_double.N() == 64);
	for(int ii=0; ii<64; ii++) 
	  dparm[ii] = param_double[ii];
      }
    }  
    else{
      iparm[IPARM_MODIFY_PARAMETER] = API_NO;
      cout << "initialize parameter" << endl;
    }
  
    iparm[IPARM_START_TASK] = API_TASK_INIT;
    iparm[IPARM_END_TASK]   = API_TASK_INIT;
    iparm[IPARM_SYM] = API_SYM_NO; // Matrix is considered nonsymetric
    
    pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm); 
    fprintf(stdout,"-- FIN INIT PARAMETERS --\n");
    init_raff = iparm[IPARM_ITERMAX];
    
    fflush(stdout);
    /* Passage en mode verbose */
    
    iparm[IPARM_RHS_MAKING] = API_RHS_B;
    if( !param_int && data_option.empty() ){
      iparm[IPARM_MATRIX_VERIFICATION] = API_YES;
      iparm[IPARM_REFINEMENT] = API_RAF_GMRES;
      iparm[IPARM_INCOMPLETE] = API_NO;
    }

    if( !param_double && data_option.empty()){
      dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
      dparm[DPARM_EPSILON_MAGN_CTRL] = 1e-32;
    }

    SYM = AA.symetrique; 
    cout << "SYM = "<< SYM << endl;
    // SYMETRIQUE
    if( SYM == 1 ){
      iparm[IPARM_SYM] = API_SYM_YES;
      //iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    }
    if( SYM == 0 ){
      iparm[IPARM_SYM] = API_SYM_NO;
      //iparm[IPARM_FACTORIZATION] = API_FACT_LU;
    }
    
    /* Scotch */
    fprintf(stdout,"-- Scotch --\n");
    fflush(stdout);
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK]   = API_TASK_ORDERING; 
    pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    
    /* Fax */
    fprintf(stdout,"-- Fax --\n");
    iparm[IPARM_START_TASK] = API_TASK_SYMBFACT;
    iparm[IPARM_END_TASK]   = API_TASK_SYMBFACT;
    pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    
    /* Blend */
    fprintf(stdout,"-- Blend --\n");
    iparm[IPARM_START_TASK] = API_TASK_ANALYSE;
    iparm[IPARM_END_TASK]   = API_TASK_ANALYSE;
    pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    
    if( SYM == 1 ){
      //iparm[IPARM_SYM] = API_SYM_YES;
      iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    }
    if( SYM == 0 ){
      //iparm[IPARM_SYM] = API_SYM_NO;
      iparm[IPARM_FACTORIZATION] = API_FACT_LU;
    }

    /* Factorisation */
    iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
    iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
    gettimeofday(&tv1, NULL);
    fprintf(stdout,"-- SOPALIN --\n");
    pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    gettimeofday(&tv2, NULL);
    fprintf(stdout,"Time to call factorization : %ld usec\n", 
	    (long)((tv2.tv_sec  - tv1.tv_sec ) * 1000000 + 
		   tv2.tv_usec - tv1.tv_usec));
    
  
    
    for(int ii=0; ii < Ncol+1; ii++)
      ia[ii] = ia[ii]-1;
    for(int ii=0; ii < ia[Ncol]-1; ii++)
      ja[ii] = ja[ii]-1;
  }
  void Solver(const MatriceMorse<double> &AA,KN_<double> &x,const KN_<double> &b) const  {
  
    struct timeval  tv1, tv2;
     // index for pastix
    for(int ii=0; ii < Ncol+1; ii++)
      ia[ii] = ia[ii]+1;
    for(int ii=0; ii < ia[Ncol]-1; ii++)
      ja[ii] = ja[ii]+1;

    // give value of the second member
    for(int ii=0; ii < Ncol; ii++)
      rhs[ii] = b[ii];
    
    //fprintf(stdout,"SOLVE STEP %ld (in FACTORIZE STEP %ld)\n",(long)ii,(long)jj);
    
    /* updo */
    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    iparm[IPARM_END_TASK]   = API_TASK_SOLVE;
    iparm[IPARM_RHS_MAKING] = API_RHS_B;
    gettimeofday(&tv1, NULL);
    pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    gettimeofday(&tv2, NULL);
    fprintf(stdout,"Time to call updown : %ld usec\n", 
	    (long)((tv2.tv_sec  - tv1.tv_sec ) * 1000000 + 
		   tv2.tv_usec - tv1.tv_usec));    
    
    //for(int jj=0; jj < Ncol; jj++)
    //  cout << "rhs["<< jj << "]=" << rhs[jj] << endl;
    
    
    //fprintf(stdout,"RAFF STEP %ld (in FACTORIZE STEP %ld)\n",(long)ii,(long)jj);
    /* raff */

    
    iparm[IPARM_START_TASK] = API_TASK_REFINE;
    iparm[IPARM_END_TASK]   = API_TASK_REFINE;
    //iparm[IPARM_RHS_MAKING] = API_RHS_B;
    iparm[IPARM_ITERMAX]    = init_raff;
    gettimeofday(&tv1, NULL);
    pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    gettimeofday(&tv2, NULL);
    fprintf(stdout,"Time to call refinement : %ld usec\n", 
	    (long)((tv2.tv_sec  - tv1.tv_sec ) * 1000000 + 
		   tv2.tv_usec - tv1.tv_usec));
    
    for(int ii=0; ii < Ncol; ii++)
      x[ii] = rhs[ii];
   
    // index for freefem
    for(int ii=0; ii < Ncol+1; ii++)
      ia[ii] = ia[ii]-1;
    for(int ii=0; ii < ia[Ncol]-1; ii++)
      ja[ii] = ja[ii]-1;
    
  }

  ~dSolvepastixmpi(){
    /* mem free */
    iparm[IPARM_START_TASK] = API_TASK_CLEAN;
    iparm[IPARM_END_TASK]   = API_TASK_CLEAN;
    pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    
    if( sizeof(pastix_int_t) != sizeof(int) )
      {
	memFree_null(ia);
	memFree_null(ja);
      }
    /* Free mem no longer necessary */
    memFree_null(rhs);
   
  }
  
  void addMatMul(const KN_<double> & x, KN_<double> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<double> &) (*this) * x; 
  }

};

MatriceMorse<double>::VirtualSolver *
BuildSolverpastixmpi(DCL_ARG_SPARSE_SOLVER(double,A))
{
    if(verbosity>9)
    cout << " BuildSolverpastixmpi<double>" << endl;
    return new dSolvepastixmpi(*A,ds.data_filename, 
			       ds.lparams, ds.dparams, ds.perm_r, ds.perm_c);
}


class Init { public:
    Init();
};

//  the 2 default sparse solver double and complex
DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; ;
//DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetDefault()
{
    if(verbosity)
	cout << " SetDefault sparse to default" << endl;
    DefSparseSolver<double>::solver =SparseMatSolver_R;
    //DefSparseSolver<Complex>::solver =SparseMatSolver_C;
    TypeSolveMat::defaultvalue =TypeSolveMat::SparseSolver;
}

bool Setpastixmpi()
{
    if(verbosity)
	cout << " SetDefault sparse solver to pastixmpi" << endl;
    DefSparseSolver<double>::solver  =BuildSolverpastixmpi;
    //DefSparseSolver<Complex>::solver =BuildSolverpastixmpi;    
    TypeSolveMat::defaultvalue  = TypeSolveMatdefaultvalue;
}



LOADINIT(Init);
Init::Init()
{ 
  
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  //SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  
  if(verbosity>1)
    cout << "\n Add: pastixmpi,  defaultsolver defaultsolverpastixmpi" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  DefSparseSolver<double>::solver =BuildSolverpastixmpi;
  //DefSparseSolver<Complex>::solver =BuildSolverpastixmpi;
  if(! Global.Find("defaultsolver").NotNull() )
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  Global.Add("defaulttopastixmpi","(",new OneOperator0<bool>(Setpastixmpi));
}
