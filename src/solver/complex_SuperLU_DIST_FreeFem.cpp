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

/*
  Interface freefem++ et SuperLU_DIST_2.3 

  /bin/sh ff-mpic++ zSuperLU_DIST.cpp -I/Users/morice/librairie/SuperLU_DIST_2.3/SRC/ -L/Users/morice/librairie/openmpi/lib/ -lmpi -lopen-pal -lopen-rte -L/Users/morice/librairie/PATCHVECLIB/ -lwrapperdotblas -framework veclib -L/Users/morice/librairie/ParMetis-3.1/ -lparmetis -lmetis -L/Users/morice/librairie/SuperLU_DIST_2.3/lib/ -lsuperlu_dist_2.3

*/
#include  <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"

//#include "lex.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "mpi.h"
#include "superlu_zdefs.h"
#include "ffsuperludistoption.hpp"


template <class R> struct SuperLUmpiDISTDriver
{
    
};
template <> struct SuperLUmpiDISTDriver<Complex>
{
  /* Driver routines */
   /* Driver routines */
  static  Dtype_t R_SLU_T() { return SLU_Z;} 
  static doublecomplex *dc(Complex *p)  { return (doublecomplex *) (void *) p;}
  static doublecomplex **dc(Complex **p)  { return (doublecomplex **) (void *) p;}

  // Remplacement doublecomplex par Complex
  static void
  pgssvx(superlu_options_t *p1, SuperMatrix *p2, 
	  ScalePermstruct_t *p3,
	  Complex *p4, int p5, int p6, gridinfo_t *p7,
	  LUstruct_t *p8, SOLVEstruct_t *p9, double *p10,
	  SuperLUStat_t *p11, int *p12)
  {
     pzgssvx( p1,p2,p3,
	      dc(p4),p5,p6,p7,p8,p9,p10,p11,p12 );
  }
    
  static void
  pgssvx_ABglobal(superlu_options_t *p1, SuperMatrix *p2, 
		  ScalePermstruct_t *p3, 
		  Complex *p4, int p5, int p6, gridinfo_t *p7,
		  LUstruct_t *p8, double *p9,
		  SuperLUStat_t *p10, int *p11)
   
  {
    pzgssvx_ABglobal( p1,p2,p3,
		      dc(p4),p5,p6,p7,p8,p9,p10,p11 );
  }


  static void
  Create_CompCol_Matrix_dist(SuperMatrix *p1, int_t p2, int_t p3, int_t p4, 
			      Complex *p5, int_t *p6, int_t *p7,
			      Stype_t p8, Dtype_t p9, Mtype_t p10)
  {
    zCreate_CompCol_Matrix_dist( p1,p2,p3,
				 p4,dc(p5),p6,p7,p8,p9,p10 );
  }

  static  void
  Create_CompRowLoc_Matrix_dist(SuperMatrix *p1, int_t p2, int_t p3,
				 int_t p4, int_t p5, int_t p6,
				 Complex *p7, int_t *p8, int_t *p9,
				 Stype_t p10, Dtype_t p11, Mtype_t p12)
  {
    zCreate_CompRowLoc_Matrix_dist( p1,p2,p3,
				    p4,p5,p6,dc(p7),p8,p9,p10,p11,p12 );
  }

  static  void
  CompRow_to_CompCol_dist(int_t p1, int_t p2, int_t p3, 
			   Complex *p4, int_t *p5, int_t *p6,
			   Complex **p7, int_t **p8, int_t **p9)
  {
    zCompRow_to_CompCol_dist( p1,p2,p3,
			      dc(p4),p5,p6,dc(p7),p8,p9 );
  }


  static void
  Create_Dense_Matrix_dist(SuperMatrix *p1, int_t p2, int_t p3, Complex *p4,
			  int_t p5, Stype_t p6, Dtype_t p7,Mtype_t p8)
  {
     zCreate_Dense_Matrix_dist( p1,p2,p3,
				dc(p4),p5,p6,p7,p8  );
  }



  static void
  Create_SuperNode_Matrix_dist(SuperMatrix *p1, int_t p2, int_t p3, int_t p4, 
			       Complex *p5, int_t *p6,
			       int_t *p7, int_t *p8,
			       int_t *p9, int_t *p10,
			       Stype_t p11, Dtype_t p12, Mtype_t p13)
  {
    zCreate_SuperNode_Matrix_dist( p1,p2,p3,
				   p4,dc(p5),p6,p7,p8,p9,p10,p11,p12,p13 );
  }
  
  static void
  Print_CompRowLoc_Matrix_dist(SuperMatrix *p1)
  {
    zPrint_CompRowLoc_Matrix_dist(p1);
  }


};

template<class R>
class ZSolveSuperLUmpi :   public MatriceMorse<R>::VirtualSolver, public SuperLUmpiDISTDriver<R>   {
  
  double eps;
  mutable double  epsr;
  double tgv;
  double tol_pivot_sym,tol_pivot; //Add 31 oct 2005
   
   
  //mutable char           equed[1];
  //yes_no_t       equil;
  mutable SuperMatrix    A;
  NCformat       *Astore;
  //NCformat       *Ustore;
  //SCformat       *Lstore;

  mutable superlu_options_t options;
  mutable mem_usage_t    mem_usage;
  mutable ScalePermstruct_t ScalePermstruct;
  mutable LUstruct_t        LUstruct;
  mutable SOLVEstruct_t     SOLVEstruct;
  mutable gridinfo_t        grid;

  string string_option;
  string data_option;
  R             *a;
  int           *asub, *xa;
  int_t m, n, nnz;
  // rajout pour //
  int_t nprow,npcol;  /* process rows and process columns*/


  int matrixdist; // type of distributed matrix
  MPI_Comm commworld ;

  static const int assembled =0;
  static const int distributedglobal =1;
  static const int distributed =2;

public:
  ZSolveSuperLUmpi(const MatriceMorse<R> &AA,int strategy,double ttgv, double epsilon=1e-6,
		  double pivot=-1.,double pivot_sym=-1., string datafile,
		   string param_char, KN<long> &pperm_r, KN<long> &pperm_c, void * ccommworld=0 ) : 
    eps(epsilon),epsr(0),
    tgv(ttgv),string_option(param_char),data_option(datafile),
    tol_pivot_sym(pivot_sym),tol_pivot(pivot)
  { 
    commworld = ccommworld ? *static_cast<MPI_Comm*>( ccommworld) : MPI_COMM_WORLD;    
    
    R*      B;
    //R*      X;
    SuperLUStat_t stat;
    int            info, ldb, nrhs=0;
    int            i;
    double*        berr;
    
    int iam;

    // Add for distributed matrix
    int_t         m_loc, m_loc_fst, fst_row, nnz_loc, fst_nnz;
    R             *aloc;
    int           *asubloc, *xaloc;
    // End Add for distributed matrix
    
    // time variables
    long int starttime,finishtime;
    long int timeused;
    if(verbosity) starttime = clock();

    A.Store=0;  
    /* Defaults */
    nrhs  = 0;

    /* lecture de nprow and npcol */
    // Cas max deux procs
    nprow = 1;
    npcol = 1;
    matrixdist=0;
    /* set the default options */
    set_default_options_dist(&options);
    DiagScale_t optionDiagScale;
    
    //if(verbosity > 10) print_options_dist(&options);

    if(!string_option.empty()) read_nprow_npcol_freefem( &string_option, &nprow, &npcol, &matrixdist);
    if(!string_option.empty()) read_options_freefem(&string_option,&options,&optionDiagScale);

    if(!data_option.empty()) read_options_superlu_datafile(&data_option,&options,&nprow, &npcol, &matrixdist,&optionDiagScale);

    //if(verbosity > 10) print_options_dist(&options);
    
     /* ------------------------------------------------------------
	 INITIALIZE THE SUPERLU PROCESS GRID. 
	 ------------------------------------------------------------*/
    cout << "Complex superlu_gridinit " << commworld << " "<< ccommworld << " : "  << nprow << "X" <<  npcol <<endl;
    superlu_gridinit(commworld , nprow, npcol, &grid);
    cout << " --\n";
    /* Bail out if I do not belong in the grid. */
    iam = grid.iam;
    if ( iam >= nprow * npcol ){
      //superlu_gridexit(&grid);
      printf("this process is not used in superlu %d \n",iam);
    }
    else
      {
	// matrix to procs and vectors
	if( matrixdist == assembled ){
	  
	  if(!iam){
	    cout <<  "iam=" << iam << endl;
	    printf("\tProcess grid\t%d X %d\n", grid.nprow, grid.npcol);
	    
	    /* create the matrix for superlu_dist */
	    n=AA.n;
	    m=AA.m;
	    nnz=AA.nbcoef;
	  
	    assert( AA.lg[n] == nnz );	   
	    printf("\tDimension\t%dx%d\t # nonzeros %d\n", m, n, nnz);
	    
	    /* transform Row to Col */
	    // cela coute cher comme fonction //
	    // dallocateA_dist(n, nnz, &a, &asub, &xa);
	    // dCompRow_to_CompCol_dist(m,n,nnz,arow,asubrow,xarow,&a,&asub,&xa);
	    
	    CompRow_to_CompCol_dist(m,n,nnz,AA.a,AA.cl,AA.lg,&a,&asub,&xa);
	  
	    /* Broadcast matrix A to the other PEs. */
	    MPI_Bcast( &m,   1,   mpi_int_t,  0, grid.comm );
	    MPI_Bcast( &n,   1,   mpi_int_t,  0, grid.comm );
	    MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid.comm );

	    MPI_Bcast( a,    nnz, SuperLU_MPI_DOUBLE_COMPLEX, 0, grid.comm );
	    MPI_Bcast( asub, nnz, mpi_int_t,  0, grid.comm );
	    MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid.comm );	    
	    
	  }
	  else{
	    /*
	      printf("\tProcess grid\t%d X %d\n", grid.nprow, grid.npcol);
	      /* Receive matrix A from PE 0. */
	    MPI_Bcast( &m,   1,   mpi_int_t,  0, grid.comm );
	    MPI_Bcast( &n,   1,   mpi_int_t,  0, grid.comm );
	    MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid.comm );
	    
	    /* Allocate storage for compressed column representation. */
	    zallocateA_dist(n, nnz, dc(&a), &asub, &xa);
	    
	    MPI_Bcast( a, nnz, SuperLU_MPI_DOUBLE_COMPLEX, 0, grid.comm );
	    MPI_Bcast( asub, nnz, mpi_int_t,  0, grid.comm );
	    MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid.comm );
	    
	  }
	  
	  Dtype_t R_SLU = SuperLUmpiDISTDriver<R>::R_SLU_T(); 
	  
	  cout << "Debut: Create_CompCol_Matrix_dist" <<endl;
	  Create_CompCol_Matrix_dist(&A, m, n, nnz, a, asub, xa, SLU_NC, R_SLU, SLU_GE);      
	  cout << "Fin: Create_CompCol_Matrix_dist" <<endl;
	  /* creation of pseudo solution + second member */
	  
	  if ( !(B = new Complex[m] ) ){
	    //if ( !(B = doublecomplexMalloc_dist(m )) ){
	    printf("probleme d allocation\n");
	    exit(1);
	  }
	  
	  if(verbosity)
	    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, nnz);

	 
// 	  /* set the default options */
// 	  set_default_options_dist(&options);
// 	  DiagScale_t optionDiagScale;
// 	  if(!string_option.empty()) read_options_freefem(&string_option,&options,&optionDiagScale);
	  
	  /* Initialize ScalePermstruct and LUstruct. */
	  ScalePermstructInit(m, n, &ScalePermstruct);
	  if(!(pperm_r==NULL)  || !(pperm_c==NULL) ) ScalePermstruct.DiagScale=optionDiagScale;
	  if( !(pperm_r==NULL) ) 
	    for(int ii=0; ii<m; ii++) ScalePermstruct.perm_r[ii] = pperm_r[ii];
	  if( !(pperm_c==NULL) )
	    for(int ii=0; ii<n; ii++) ScalePermstruct.perm_c[ii]= pperm_c[ii];
	  
	  if( ScalePermstruct.DiagScale != NOEQUIL ){
	    printf("FreeFem++ doesn't support change of the original matrix"); 
	    exit(1);
	  }
	  LUstructInit(m, n, &LUstruct);
	  
	  /* Initialize the statistics variables. */
	  PStatInit(&stat);
	  
	  ldb = m;
	  nrhs=1;
	  if ( !(berr = doubleMalloc_dist(nrhs )) ){
	    printf("probleme d allocation\n");
	    exit(1);
	  }
	  berr[0]=0.;
    	
	  if(verbosity)
	    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, nnz);
	
	  /* INIT LU struct*/
	  
	  /* ONLY PERFORM THE LU DECOMPOSITION */
	  //B.ncol = 0;  /* Indicate not to solve the system */
	  
	  nrhs=0;
	  SuperLUmpiDISTDriver<R>::pgssvx_ABglobal(&options, &A,  &ScalePermstruct, B, ldb, nrhs, &grid,
					       &LUstruct, berr, &stat, &info);
	
	  if(verbosity)
	    printf("LU factorization: pdgssvx()/p returns info %d\n", info);
	  
	  if ( verbosity) PStatPrint(&options,&stat,&grid);
	  PStatFree(&stat);
	  	 
	}
	//##########################################################
	//
	//       matrix distributed with matrix global given
	//
	//##########################################################
	else if( matrixdist == distributedglobal) {
	   if(!iam){

	     printf("\tProcess grid\t%d X %d\n", grid.nprow, grid.npcol);
	
	     /* create the matrix for superlu_dist */
	     n=AA.n;
	     m=AA.m;
	     nnz=AA.nbcoef;
	     a=AA.a;
	     asub=AA.cl;
	     xa=AA.lg;
	     
	     xa[n] = nnz;
	     printf("\tDimension\t%dx%d\t # nonzeros %d\n", m, n, nnz);
	     
	     /* Broadcast matrix A to the other PEs. */
	     MPI_Bcast( &m,   1,   mpi_int_t,  0, grid.comm );
	     MPI_Bcast( &n,   1,   mpi_int_t,  0, grid.comm );
	     MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid.comm );
	     
	     MPI_Bcast( AA.a,    nnz, SuperLU_MPI_DOUBLE_COMPLEX, 0, grid.comm );
	     MPI_Bcast( AA.cl, nnz, mpi_int_t,  0, grid.comm );
	     MPI_Bcast( AA.lg,   n+1, mpi_int_t,  0, grid.comm );
	     
	     
	   }
	   else
	     {
	     
	     printf("\tProcess grid\t%d X %d\n", grid.nprow, grid.npcol);
	     /* Receive matrix A from PE 0. */
	     MPI_Bcast( &m,   1,   mpi_int_t,  0, grid.comm );
	     MPI_Bcast( &n,   1,   mpi_int_t,  0, grid.comm );
	     MPI_Bcast( &nnz, 1,   mpi_int_t,  0, grid.comm );
	     
	     /* Allocate storage for compressed column representation. */
	     zallocateA_dist(n, nnz, dc(&a), &asub, &xa);
	     
	     MPI_Bcast( a, nnz, SuperLU_MPI_DOUBLE_COMPLEX, 0, grid.comm );
	     MPI_Bcast( asub, nnz, mpi_int_t,  0, grid.comm );
	     MPI_Bcast( xa,   n+1, mpi_int_t,  0, grid.comm );

	   }
	   
	   /* Compute the number of rows to be distributed to local process */
	   m_loc = m / (grid.nprow * grid.npcol); 
	   m_loc_fst = m_loc;
	   /* When m / procs is not an integer */
	   if ((m_loc * grid.nprow * grid.npcol) != m) {
	     /*m_loc = m_loc+1;
	       m_loc_fst = m_loc;*/
	     if (iam == (grid.nprow * grid.npcol - 1)) /* last proc. gets all*/
	       m_loc = m - m_loc * (grid.nprow * grid.npcol - 1);
	   }
	   
	   fst_row = iam * m_loc_fst;
	   
	   nnz_loc = xa[fst_row+m_loc]-xa[fst_row];
	   zallocateA_dist(m_loc, nnz_loc, dc(&aloc), &asubloc, &xaloc);
	   
	   //xaloc = (int_t*) intMalloc_dist(m_loc+1);
	   for(int ii=0; ii < m_loc; ii++){
	     xaloc[ii] = xa[fst_row+ii]-xa[fst_row];	
	   }
	   
	   xaloc[m_loc]=nnz_loc;
	   
	   fst_nnz = xa[fst_row];
	   //aloc    = new R[nnz_loc];
	   //aloc    = (Complex*) doubleMalloc_dist(nnz_loc);
	   //asubloc = (int_t*)  intMalloc_dist(nnz_loc);
	   
	   for(int ii=0; ii < nnz_loc; ii++){
	     aloc[ii] = a[fst_nnz+ii];
	     asubloc[ii] = asub[fst_nnz+ii];
	   }
	   
	   if( iam ){
	     SUPERLU_FREE( a );
	     SUPERLU_FREE( asub );
	     SUPERLU_FREE( xa );
	   }
	   
	   Dtype_t R_SLU = SuperLUmpiDISTDriver<R>::R_SLU_T(); 
	   
	   cout << "Debut: Create_CompRowCol_Matrix_dist" <<endl;
	   Create_CompRowLoc_Matrix_dist(&A, m, n, nnz_loc, m_loc, fst_row, aloc, asubloc, xaloc, SLU_NR_loc, R_SLU, SLU_GE);
	   
	   cout << "Fin: Create_CompRowCol_Matrix_dist" <<endl;
	   /* creation of pseudo solution + second member */
	   
	   
	   if ( !(B = new Complex[m_loc]) ){
	     printf("probleme d allocation\n");
	     exit(1);
	   }
	   
	   for(int ii=0; ii < m_loc; ii++){
	     B[ii] = 1.; //BB[fst_row+ii];
	   }
     
	   if(verbosity)
	     printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, nnz);
	   
	   /* set the default options */
	   set_default_options_dist(&options);
	   DiagScale_t optionDiagScale;
	   if(!string_option.empty()) read_options_freefem(&string_option,&options,&optionDiagScale);
	   	   
	   m=A.nrow;
	   n=A.ncol;
	   printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, nnz);
	   /* Initialize ScalePermstruct and LUstruct. */
	   ScalePermstructInit(m, n, &ScalePermstruct);
	   if(pperm_r  ||  pperm_c ) ScalePermstruct.DiagScale=optionDiagScale;
	   if(pperm_r) 
	     for(int ii=0; ii<m; ii++) ScalePermstruct.perm_r[ii] = pperm_r[fst_row+ii];
	   if(pperm_c) 
	     for(int ii=0; ii<n; ii++) ScalePermstruct.perm_c[ii] = pperm_c[ii];
	   
	   LUstructInit(m, n, &LUstruct);
	   


	   /* Initialize the statistics variables. */
	   PStatInit(&stat);
	   
	   ldb = m_loc;
	   //ldx = m_loc;
	   
	   nrhs=1;
	   if ( !(berr = doubleMalloc_dist(nrhs )) ){
	     printf("probleme d allocation\n");
	     exit(1);
	   }
	   berr[0]=0.;
	   
	   /* ONLY PERFORM THE LU DECOMPOSITION */
    
	   nrhs=0;
	   SuperLUmpiDISTDriver<R>::pgssvx(&options, &A,  &ScalePermstruct, B, ldb, nrhs, &grid,
					   &LUstruct, &SOLVEstruct, berr, &stat, &info);
	   
	   if(verbosity)
	     printf("LU factorization: pdgssvx()/p returns info %d\n", info);
	   
	   if ( verbosity) PStatPrint(&options,&stat,&grid);
	   PStatFree(&stat);
	}
	else if( matrixdist == distributed) {
	  printf("in construction\n");
	  exit(1);
	}
	else{
	  printf("matrix choice for SuperLU_DIST is assembled, distributedglobal and distributed \n");
	  exit(1);
	}
	
	delete [] B;
	options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
	nrhs=1;
	SUPERLU_FREE(berr);  	

	if(iam==0){
	  finishtime = clock();
	  timeused= (finishtime-starttime)/(1000 );
	  printf("=====================================================\n");
	  cout << "SuperLU_DIST : time factorisation :: " << timeused << " ms" <<endl;
	  printf("=====================================================\n");
	}
	 
      }
  }

  void Solver(const MatriceMorse<R> &AA,KN_<R> &x,const KN_<R> &b) const  {
    R*        B;
    SuperLUStat_t  stat;
    int            iam;
    int            info=0, ldb=m, nrhs=1;
    int            i;
    double*        berr;
    double         ferr; 
    double         rpg, rcond;
      
    int_t    m_loc,m_loc_fst,fst_row;
    
    // time variables
    long int starttime,finishtime;
    long int timeused;
    if(verbosity) starttime = clock();

    if(n != m) exit(1);

    ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;

    Dtype_t R_SLU = SuperLUmpiDISTDriver<R>::R_SLU_T(); 
    nrhs= 1;
    
     
    /* Initialize the statistics variables. */
    PStatInit(&stat);

    iam = grid.iam;
    if( iam < nprow*npcol){
      /* cas matrix assembled */ 
      if( matrixdist == assembled ){
	
	if( !(B = new R[m*nrhs] ) ){
	  printf("probleme d allocation\n");
	  exit(1);
	}
	
	for(int ii=0; ii<n; ii++){
	  B[ii]=b[ii];
	}
	
	if ( !(berr = doubleMalloc_dist(nrhs )) ){
	  printf("probleme d allocation\n");
	  exit(1);
	}
	berr[0]=0.;
	
	options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */   
	ldb = m;
	//nrhs= 1;
	SuperLUmpiDISTDriver<R>::pgssvx_ABglobal (&options, &A, &ScalePermstruct, B, ldb, nrhs, &grid,
						  &LUstruct, berr, &stat, &info );
	
	if(verbosity)
	  printf("Triangular solve: dgssvx() returns info %d\n", info);
	
	if(verbosity) PStatPrint(&options, &stat, &grid);   
	
	for(int ii=0; ii<n; ii++){
	  x[ii] = B[ii]; 
	}
	
	if(verbosity) cout << "   x min max " << x.min() << " " <<x.max() << endl;
	
      }
      else if( matrixdist == distributedglobal) {
	R*    xtemp;
	iam = grid.iam;
	/* Compute the number of rows to be distributed to local process */
	m_loc = m / (grid.nprow * grid.npcol); 
	m_loc_fst = m_loc;
	/* When m / procs is not an integer */
	if ((m_loc * grid.nprow * grid.npcol) != m) {
	  /*m_loc = m_loc+1;
	    m_loc_fst = m_loc;*/
	  if (iam == (grid.nprow * grid.npcol - 1)) /* last proc. gets all*/
	    m_loc = m - m_loc * (grid.nprow * grid.npcol - 1);
	}
	
	fst_row = iam * m_loc_fst;
	
	if ( !(B = new R[m_loc] ) ){
	  printf("probleme d allocation\n");
	  exit(1);
	}
	
	//printf("initilisation B:");
	for(int ii=0; ii<m_loc; ++ii){
	  B[ii] = b[ii+fst_row];
	  //printf("  B[%d]= %f  ",ii,B[ii]);
	}
	//printf(" :: fin \n");
	//fflush(stdout);
	
      
	if ( !(berr = doubleMalloc_dist(nrhs )) ){
	  printf("probleme d allocation\n");
	  exit(1);
	}
	berr[0]=0.;
	
	options.Fact = FACTORED; /* Indicate the factored form of A is supplied. */
	//options.Equil = YES;
	//options.Trans = TRANS;
	
	
	ldb = m;
	SuperLUmpiDISTDriver<R>::pgssvx(&options, &A, &ScalePermstruct, B, ldb, nrhs, &grid,
					&LUstruct, &SOLVEstruct, berr, &stat, &info );
	
	if(verbosity)
	  printf("Triangular solve: dgssvx() returns info %d\n", info);
	
	if ( !(xtemp = new R[AA.n]) ){
	  printf("probleme d allocation de xtemp\n");
	  exit(1);
	}
	
      
	int disp[nprow*npcol];
	MPI_Allgather(&fst_row, 1, MPI_INT, disp, 1, MPI_INT, grid.comm);
	
	int recv[nprow*npcol];
	MPI_Allgather(&m_loc, 1, MPI_INT, recv, 1, MPI_INT, grid.comm);
	
	MPI_Allgatherv(B, m_loc, SuperLU_MPI_DOUBLE_COMPLEX, xtemp, recv, disp, SuperLU_MPI_DOUBLE_COMPLEX, grid.comm);
	
	for(int ii= 0; ii< AA.n ; ii++)
	  x[ii] = xtemp[ii];
	
	if(verbosity) cout << "   x min max " << x.min() << " " <<x.max() << endl;
	
	delete [] xtemp;
      
      }
      else if( matrixdist == distributed) {
	printf("in construction\n");
	exit(1);
      }
      else{
	printf("matrix choice for SuperLU_DIST is assembled, distributedglobal and distributed \n");
	exit(1);
      }
      

      delete [] B;
      SUPERLU_FREE( berr );
      
      PStatFree(&stat);

      if(iam==0){
	finishtime = clock();
	timeused= (finishtime-starttime)/(1000 );
	printf("=====================================================\n");
	cout << " SuperLU_DIST : time solve  :: " << timeused << " ms" <<endl;
	printf("=====================================================\n");
      }

    }
  }
    
  ~ZSolveSuperLUmpi() { 
    int iam;
    iam = grid.iam;
    if(iam < nprow*npcol){
      if(verbosity)
	cout << "~SolveSuperLUmpi Complex:" << endl;
      
      if( matrixdist == assembled) {
	//if( A.Store)  Destroy_CompCol_Matrix_dist(&A);
	//if( L.Store && U.Store )  {
	Destroy_LU(n, &grid, &LUstruct);
	ScalePermstructFree(&ScalePermstruct);
	LUstructFree(&LUstruct);
	//}
	if ( options.SolveInitialized ) {
	  zSolveFinalize(&options, &SOLVEstruct);
	}
      }
      else if( matrixdist == distributedglobal) {
	if( A.Store)  Destroy_CompRowLoc_Matrix_dist(&A);
	
	Destroy_LU(n, &grid, &LUstruct);
	ScalePermstructFree(&ScalePermstruct);
	LUstructFree(&LUstruct);
	
	if ( options.SolveInitialized ) {
	  zSolveFinalize(&options, &SOLVEstruct);
	}
      }
      else if( matrixdist == distributed) {
	printf("in construction\n");
	exit(1);
      }
      else{
	printf("matrix choice for SuperLU_DIST is assembled, distributedglobal and distributed \n");
	exit(1);
      }
    }
    printf("superlu_gridexit(&grid), %d\n",iam);
    superlu_gridexit(&grid); 
   
  }
  void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<R> &) (*this) * x; 
  }
     
}; 




MatriceMorse<Complex>::VirtualSolver *
BuildSolverSuperLUmpi(DCL_ARG_SPARSE_SOLVER(Complex,A))
{
  if(verbosity>9)
      cout << " BuildSolverSuperLUmpi<double>" << endl;
  return new ZSolveSuperLUmpi<Complex>(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym,
				       ds.data_filename, ds.sparams, ds.perm_r, ds.perm_c);
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

bool SetSuperLUmpi()
{
    if(verbosity)
	cout << " SetDefault sparse solver to SuperLUmpi double" << endl;
    //DefSparseSolver<double>::solver  =BuildSolverSuperLUmpi;
    DefSparseSolver<Complex>::solver =BuildSolverSuperLUmpi;    
    TypeSolveMat::defaultvalue  = TypeSolveMatdefaultvalue;
}



Init init;
Init::Init()
{ 
  
  //SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  
  if(verbosity>1)
    cout << "\n Add: Complex SuperLU_DIST,  defaultsolver defaultsolverSuperLUdist" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  //DefSparseSolver<double>::solver =BuildSolverSuperLUmpi;
  DefSparseSolver<Complex>::solver =BuildSolverSuperLUmpi;
  if(! Global.Find("defaultsolver").NotNull() )
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  Global.Add("complexdefaulttoSuperLUdist","(",new OneOperator0<bool>(SetSuperLUmpi));
}

