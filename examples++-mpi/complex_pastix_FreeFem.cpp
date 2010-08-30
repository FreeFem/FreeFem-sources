// ORIG-DATE: 02/2009
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
//

/* 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */
//ff-c++-LIBRARY-dep: complex_pastix   blas parmetis metis scotch mpi fc
//ff-c++-cpp-dep: 
/*
  Interface entre freefem++ et pastix
*/
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


#include <mpi.h>
#include "pastix_long_complex.h"

#undef memFree_null
#define memFree_null(x) {if (x ==NULL) {fprintf(stdout,"%s:%d freeing NULL\n",__FILE__,__LINE__);} free(x); x=NULL;}

#define STR_SIZE 256

static pastix_int_t * pastixint(int * ii){ return (pastix_int_t*) (void *) ii;} 
static pastix_float_t * pastixfloat(Complex * ii){ return (pastix_float_t*) (void *) ii;} 

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
	     Complex *a, int *colind, int  *rowptr,
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
//void read_datafile_pastixff(const string &datafile, pastix_int_t *iparmtab, double *dparmtab){
void read_datafile_pastixff(const string &datafile, int &mpi_flag, pastix_int_t *iparmtab, double *dparmtab){
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

  
  if( !(strcmp(token,"matrix") == 0) ){
    printf("freefem++: error in reading matrix parameter for pastix (see strcuture of ffpastix_iparm_dparm.txt) \n");
    exit(1);
  }
  else{
    printf("freefem++: reading matrix parameter for pastix \n");    
  }

  fgets(szbuff,MAX_CHAR_PER_LINE,m_File);
  token = strtok(szbuff," /#!\t\n");
  
  if(strcmp(token,"assembled") == 0)
    mpi_flag = 0;
  else if(strcmp(token,"distributedglobal") == 0) 
    mpi_flag = 1;
  else if(strcmp(token,"distributed") == 0) 
    mpi_flag = 2;
  else{
    printf("value of parameter matrix is not correct %s \n", token );
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

// CAS COMPLEX SEULEMENT 


class zSolvepastixmpi :   public MatriceMorse<Complex>::VirtualSolver   {
  
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

  zSolvepastixmpi(const MatriceMorse<Complex> &AA,int strategy,double ttgv, double epsilon,
		  double pivot,double pivot_sym, string datafile, KN<long> &param_int, KN<double> &param_double, 
		  KN<long> &pperm_r, KN<long> &pperm_c) : 
    eps(epsilon),epsr(0),
    tgv(ttgv),tol_pivot_sym(pivot_sym),tol_pivot(pivot),
    data_option(datafile) 
  { 
    //KN_<long> param_int(pparam_int);
    //KN_<double> param_double(pparam_double);

    //int m;
    //int ierr;
    struct timeval  tv1, tv2;
    int nnz;
    // time variables
    long int starttime,finishtime;
    long int timeused;
    if(verbosity) starttime = clock();

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
    // MPI_flag need to unselect for non distributed matrix
    mpi_flag  = 0;
    thrd_flag = 0;

    // ######################  
    //pastix_int_t init_raff;
    fprintf(stdout,"-- INIT PARAMETERS --\n");
    
    // reading iparm from array    
    if(!data_option.empty()){
      read_datafile_pastixff(data_option,mpi_flag,iparm,dparm);
      if(mpi_flag != 0) 
	cerr << "ERROR :: GLOBAT INPUT MATRIX FOR ALL PROCS  matrix=assembled" << endl;
    }
    else if( !(param_int==NULL) || !(param_double==NULL)){
	if( ! (param_int==NULL) ) 
      {
	cout << "internal param_int" << endl;
	assert(param_int.N() == 64);
	for(int ii=0; ii<64; ii++) 
	  iparm[ii] = param_int[ii];
	iparm[IPARM_MODIFY_PARAMETER] = API_YES;
      }
	if( !(param_double==NULL) ) 
      {
	cout << "internal param_double" << endl;
	assert(param_double.N() == 64);
	for(int ii=0; ii<64; ii++) 
	  dparm[ii] = param_double[ii];
      }
    }  
    else{
      iparm[IPARM_MODIFY_PARAMETER] = API_NO;
      cout << "initialize default parameter" << endl;
    }
    
    //################################
    if(myid == 0){
      Ncol = AA.m;
      Nrow = AA.n;
      nnz  = AA.nbcoef;
      // Avant : on ecrit la transposée
      
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
      MPI_Bcast( &Ncol,   1,    MPI_INT,   0, MPI_COMM_WORLD );
      MPI_Bcast( &Nrow,   1,    MPI_INT,   0, MPI_COMM_WORLD );
      MPI_Bcast( &nnz,    1,    MPI_INT,   0, MPI_COMM_WORLD );

      MPI_Bcast( avals, nnz,    MPI_PASTIX_FLOAT, 0, MPI_COMM_WORLD );
      MPI_Bcast(    ia, Ncol+1, MPI_PASTIX_INT,   0, MPI_COMM_WORLD );
      MPI_Bcast(    ja, nnz,    MPI_PASTIX_INT,   0, MPI_COMM_WORLD );
    }
    else{
      MPI_Bcast( &Ncol, 1,        MPI_INT,  0, MPI_COMM_WORLD );
      MPI_Bcast( &Nrow, 1,        MPI_INT,  0, MPI_COMM_WORLD );
      MPI_Bcast( &nnz,  1,        MPI_INT,  0, MPI_COMM_WORLD );
      
      avals = (pastix_float_t *) malloc( nnz*sizeof(pastix_float_t) );
      ia = (pastix_int_t *) malloc( (Ncol+1)*sizeof(pastix_int_t) );
      ja = (pastix_int_t *) malloc( nnz*sizeof(pastix_int_t) );

      MPI_Bcast( avals, nnz,  MPI_PASTIX_FLOAT,   0, MPI_COMM_WORLD );
      MPI_Bcast(    ia, Ncol+1, MPI_PASTIX_INT,   0, MPI_COMM_WORLD );
      MPI_Bcast(    ja, nnz,    MPI_PASTIX_INT,   0, MPI_COMM_WORLD );
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
   
  
    iparm[IPARM_START_TASK] = API_TASK_INIT;
    iparm[IPARM_END_TASK]   = API_TASK_INIT;
    iparm[IPARM_SYM] = API_SYM_NO; // Matrix is considered nonsymetric    
    if(mpi_flag == 0)
      pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm); 
    else
      cerr << "error :: mpi_flag = 0 for calling pastix" << endl; 
    fprintf(stdout,"-- FIN INIT PARAMETERS --\n");
    init_raff = iparm[IPARM_ITERMAX];
    cout << "init_raff=" << init_raff << endl;
    fflush(stdout);
    /* Passage en mode verbose */
    
    iparm[IPARM_RHS_MAKING] = API_RHS_B;
    if( (param_int==NULL) && data_option.empty() ){
      iparm[IPARM_MATRIX_VERIFICATION] = API_YES;
      iparm[IPARM_REFINEMENT] = API_RAF_GMRES;
      iparm[IPARM_INCOMPLETE] = API_NO;
    }

    if( (param_double==NULL) && data_option.empty()){
      dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
      dparm[DPARM_EPSILON_MAGN_CTRL] = 1e-32;
    }

  
 //    cscd_checksym(Ncol, ia, ja, loc2glob, MPI_COMM_WORLD);
    
//     if (iparm[IPARM_SYM]==API_SYM_YES)
//       {
// 	/* Symetric problem */
// 	/* Build non oriented graph */
// 	/* build non symmetric csc from symmetric csc */
// 	/*maillage global*/
// 	INT *tmpia;
// 	INT *tmpja;
// 	INT  tmpn;
	
// 	cscd_symgraph_int(*n2,   *col2,  *row2 , NULL,
// 			  &tmpn, &tmpia, &tmpja, NULL,
// 			  *loc2glob2, pastix_comm, API_YES);
	
// 	memFree_null(*col2);
// 	*col2 = tmpia;
// 	memFree_null(*row2);
// 	*row2 = tmpja;
// 	*n2   = tmpn;
//       }
    

    SYM = AA.symetrique; 
    cout << "SYM = "<< SYM << endl;
    // SYMETRIQUE
    if( SYM == 1 ){
      iparm[IPARM_SYM] = API_SYM_YES;
      iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    }
    if( SYM == 0 ){
      iparm[IPARM_SYM] = API_SYM_NO;
      iparm[IPARM_FACTORIZATION] = API_FACT_LU;
    }
    
    /* Scotch */
    fprintf(stdout,"-- Scotch --\n");
    fflush(stdout);
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK]   = API_TASK_ORDERING; 
    if(mpi_flag == 0)
      pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    else
      cerr << "error :: mpi_flag = 0 for calling pastix" << endl;  
    iparm[IPARM_SYM] = API_SYM_NO;
    /* Fax */
    fprintf(stdout,"-- Fax --\n");
    iparm[IPARM_START_TASK] = API_TASK_SYMBFACT;
    iparm[IPARM_END_TASK]   = API_TASK_SYMBFACT;
    if(mpi_flag == 0)
      pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    else
      cerr << "error :: mpi_flag = 0 for calling pastix" << endl; 
    /* Blend */
    fprintf(stdout,"-- Blend --\n");
    iparm[IPARM_START_TASK] = API_TASK_ANALYSE;
    iparm[IPARM_END_TASK]   = API_TASK_ANALYSE;
    if( SYM == 1 ){
      iparm[IPARM_SYM] = API_SYM_YES;
      iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    }
    if( SYM == 0 ){
      iparm[IPARM_SYM] = API_SYM_NO;
      iparm[IPARM_FACTORIZATION] = API_FACT_LU;
    }
    if(mpi_flag == 0)
      pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    else
      cerr << "error :: mpi_flag = 0 for calling pastix" << endl; 
   
    /* Factorisation */
    iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
    iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
    gettimeofday(&tv1, NULL);
    fprintf(stdout,"-- SOPALIN --\n");
    if(mpi_flag == 0)
      pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    else
       cerr << "error :: mpi_flag = 0 for calling pastix" << endl; 
    gettimeofday(&tv2, NULL);
    fprintf(stdout,"Time to call factorization : %ld usec\n", 
	    (long)((tv2.tv_sec  - tv1.tv_sec ) * 1000000 + 
		   tv2.tv_usec - tv1.tv_usec));
    
    for(int ii=0; ii < ia[Ncol]-1; ii++)
      ja[ii] = ja[ii]-1;
    
    for(int ii=0; ii < Ncol+1; ii++)
      ia[ii] = ia[ii]-1;
    
    //for(int ii=0; ii < ia[Ncol]-1; ii++)
    //  ja[ii] = ja[ii]-1;

    if(myid==0){
      finishtime = clock();
      timeused= (finishtime-starttime)/(1000 );
      printf("=====================================================\n");
      cout << " pastix : time factorization  :: " << timeused << " ms" <<endl;
      printf("=====================================================\n");
    }

    
  }
  void Solver(const MatriceMorse<Complex> &AA,KN_<Complex> &x,const KN_<Complex> &b) const  {
  
    struct timeval  tv1, tv2;
    // time variables
    long int starttime,finishtime;
    long int timeused;
    if(verbosity) starttime = clock();
    
    // index for pastix    
    for(int ii=0; ii < Ncol+1; ii++)
      ia[ii] = ia[ii]+1;
    for(int ii=0; ii < ia[Ncol]-1; ii++)
      ja[ii] = ja[ii]+1;
    
    
    // give value of the second member
    for(int ii=0; ii < Ncol; ii++){
      rhs[ii] = b[ii];  
    }
    
  


    //fprintf(stdout,"SOLVE STEP %ld (in FACTORIZE STEP %ld)\n",(long)ii,(long)jj);
    
    /* updo */
    iparm[IPARM_START_TASK] = API_TASK_SOLVE;
    iparm[IPARM_END_TASK]   = API_TASK_SOLVE;
    iparm[IPARM_RHS_MAKING] = API_RHS_B;
    gettimeofday(&tv1, NULL);
    if(mpi_flag == 0)
      pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    else
      cerr << "error :: mpi_flag = 0 for calling pastix" << endl; 
    gettimeofday(&tv2, NULL);
    fprintf(stdout,"Time to call updown : %ld usec\n", 
	    (long)((tv2.tv_sec  - tv1.tv_sec ) * 1000000 + 
		   tv2.tv_usec - tv1.tv_usec));    
    
    if(verbosity > 1)
      for(int jj=0; jj < Ncol; jj++)
	cout << "rhs["<< jj << "]=" << rhs[jj] << endl;
    
    
    //fprintf(stdout,"RAFF STEP %ld (in FACTORIZE STEP %ld)\n",(long)ii,(long)jj);
    /* raff */
    
    
    iparm[IPARM_START_TASK] = API_TASK_REFINE;
    iparm[IPARM_END_TASK]   = API_TASK_REFINE;
    iparm[IPARM_RHS_MAKING] = API_RHS_B;
    iparm[IPARM_ITERMAX]    = init_raff;
    gettimeofday(&tv1, NULL);
    if(mpi_flag == 0)
      pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
    else
      cerr << "error :: mpi_flag = 0 for calling pastix" << endl; 
    gettimeofday(&tv2, NULL);
    fprintf(stdout,"Time to call refinement : %ld usec\n", 
	    (long)((tv2.tv_sec  - tv1.tv_sec ) * 1000000 + 
		   tv2.tv_usec - tv1.tv_usec));

    
    for(int ii=0; ii < Ncol; ii++)
      x[ii] = rhs[ii];
       
    // index for freefem
    for(int ii=0; ii < ia[Ncol]-1; ii++)
      ja[ii] = ja[ii]-1;
    for(int ii=0; ii < Ncol+1; ii++)
      ia[ii] = ia[ii]-1;
    
    //for(int ii=0; ii < ia[Ncol]-1; ii++)
    //  ja[ii] = ja[ii]-1;

    if(myid==0){
      finishtime = clock();
      timeused= (finishtime-starttime)/(1000 );
      printf("=====================================================\n");
      cout << " pastix : time solve  :: " << timeused << " ms" <<endl;
      printf("=====================================================\n");
    }
    
    
  }

  ~zSolvepastixmpi(){
    /* mem free */
    iparm[IPARM_START_TASK] = API_TASK_CLEAN;
    iparm[IPARM_END_TASK]   = API_TASK_CLEAN;
    
    pastix(&pastix_data, MPI_COMM_WORLD, Ncol,ia,ja,avals,perm,invp,rhs,1,iparm,dparm);
   
    memFree_null(ia);
    memFree_null(ja);
    
    /* Free mem no longer necessary */
    memFree_null(perm);
    memFree_null(rhs);
    
  }
  
  void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<Complex> &) (*this) * x; 
  }

};

MatriceMorse<Complex>::VirtualSolver *
BuildSolverpastix_complex_mpi(DCL_ARG_SPARSE_SOLVER(Complex,A))
{
    if(verbosity>9)
    cout << " BuildSolverpastix_complex_mpi<complex>" << endl;
    return new zSolvepastixmpi(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym, ds.data_filename, 
			       ds.lparams, ds.dparams, ds.perm_r, ds.perm_c);
}


class Init { public:
    Init();
};

//  the 2 default sparse solver double and complex
//DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; ;
DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetDefault()
{
    if(verbosity)
	cout << " SetDefault sparse to default" << endl;
    //DefSparseSolver<double>::solver =SparseMatSolver_R;
    DefSparseSolver<Complex>::solver =SparseMatSolver_C;
    TypeSolveMat::defaultvalue =TypeSolveMat::SparseSolver;
}

bool Setpastixmpi()
{
    if(verbosity)
	cout << " SetDefault sparse solver to pastixmpi" << endl;
    //DefSparseSolver<double>::solver  =BuildSolverpastix_complex_mpi;
    DefSparseSolver<Complex>::solver =BuildSolverpastix_complex_mpi;    
    TypeSolveMat::defaultvalue  = TypeSolveMatdefaultvalue;
}



Init init;
Init::Init()
{ 
  
  //SparseMatSolver_R= DefSparseSolver<d>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  
  if(verbosity>1)
    cout << "\n Add: pastix,  defaultsolver defaultsolverpastix" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  //DefSparseSolver<double>::solver =BuildSolverpastix_complex_mpi;
  DefSparseSolver<Complex>::solver =BuildSolverpastix_complex_mpi;
  if(! Global.Find("defaultsolver").NotNull() )
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  Global.Add("complexdefaulttopastix","(",new OneOperator0<bool>(Setpastixmpi));
}
