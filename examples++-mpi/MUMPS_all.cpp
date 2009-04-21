#include  <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"

//#include "lex.hpp"
#include "MatriceCreuse_tpl.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"
#include "dmumps_c.h"
#include "zmumps_c.h"

// read options for MUMPS in freefem++
int s_(char* str, const char* cmp[])
{
  int i = 0;
  while( cmp[i] != 0){
    if( strcmp(str, cmp[i]) == 0){
      //cout << *str << " return" << i << endl;
      return i+1 ;
    }
    i++;
  }
  //cout << *str << " return 0" << endl;
  return 0;
}

void read_options_freefem(string *string_option, int *SYM, int *PAR){
  
  static const char* comp[] = {"SYM", "PAR", 0};

  char data[string_option->size()+1];  
  strcpy( data, string_option->c_str()); 
  cout << "data=" << data << endl;
  char * tictac;
  tictac = strtok(data," =,\t\n");
  cout << "tictac=" << data << endl;

  while(tictac != NULL){
    
    int id_option = s_(tictac, comp);
    tictac = strtok(NULL," =,\t\n");
    int val_options;

    switch (id_option)
      { 
      case 1 : // SYM
	*SYM = atoi(tictac);
	// strtol ???
	if(*SYM != 0 && *SYM !=1) 
	  cout << "SYM must be equal to 1 or 0 for MUMPS" << endl;
	
	break;
      case 2:  // PAR
	*PAR = atoi(tictac);
	if(*PAR != 0 && *PAR !=1) 
	  cout << "PAR must be equal to 1 or 0" << endl;
	//strtol ???
	break;
      case 0: // Equivalent of case default
	break;
      }  
    tictac = strtok(NULL," =,\t\n");
  }
  
}

class dSolveMUMPSmpi :   public MatriceMorse<double>::VirtualSolver   {
  
  double eps;
  mutable double  epsr;
  double tgv;
  double tol_pivot_sym,tol_pivot; //Add 31 oct 2005

  double            *a;
  int       *irn, *jcn;
  int          n, m, nz; 
 
  // parameter MUMPS
  KN<int>        iparam;
  KN<int>        perm_r; /* row permutations from partial pivoting */
  KN<int>        perm_c;
  KN<double>     scale_r;
  KN<double>     scale_c;
  string string_option;
  string data_option;
  int SYM;
  int PAR;
  int myid;

  // distribuer
  int nz_loc;
  int *jcn_loc, *irn_loc;
  double *a_loc;
  

  long int starttime,finishtime;
  long int timeused;

  static const int JOB_INIT=-1;
  static const int JOB_END=-2;
  static const int USE_COMM_WORLD= -987654;

  // variable reel
  mutable DMUMPS_STRUC_C id;
  
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
public:
  dSolveMUMPSmpi(const MatriceMorse<double> &AA,int strategy,double ttgv, double epsilon=1e-6,
		 double pivot=-1.,double pivot_sym=-1., string param_string, string datafile, KN<long> &param_int, 
		 KN<long> &pperm_r, KN<long> &pperm_c, KN<double> &pscale_r,KN<double> &pscale_c
		) : 
    eps(epsilon),epsr(0),
    tgv(ttgv), string_option(param_string), data_option(datafile), perm_r(pperm_r), perm_c(pperm_c), 
    tol_pivot_sym(pivot_sym),tol_pivot(pivot), scale_r(pscale_r), scale_c(pscale_c)
  { 
    if(verbosity) starttime = clock();
    int dataint[40];
    int ierr;
  
    n    = AA.n;
    m    = AA.m; 
    nz   = AA.nbcoef;
    

    if(param_int) 
      assert( param_int.N() == 40);
    if(pperm_r)
      assert( perm_r.N() == n);
    if(pscale_r) 
      assert( scale_r.N() == n);
    if(pscale_c) 
      assert( scale_c.N() == m);

    if( n != m )
      cerr << "only square matrix are supported by MUMPS" << endl;
   
     /* ------------------------------------------------------------
       INITIALIZE THE MUMPS PROCESS GRID. 
       ------------------------------------------------------------*/
   
    //ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // initialisation par defaut
    SYM=0; PAR=1;
    
    if(!string_option.empty()) 
      {
	if(myid==0) cout << "read string option" <<endl;
	read_options_freefem(&string_option,&SYM,&PAR);
      }
    
    if(!data_option.empty())
      {	
	//if(myid==0) cout << "lecture data option file" <<endl;
	char * retfile= new char[data_option.size()+1];
	strcpy(retfile, (&data_option)->c_str());
	if(myid==0) printf("lecture data option file %s\n",retfile);
	FILE *pFile=fopen(retfile,"rt");
	
	int     i_data=0;
	char    data[256];
	char   *tictac;
	
	fgets(data,256,pFile);
	tictac = strtok(data," /!#\t\n");
	SYM = atoi(tictac);
	
	fgets(data,256,pFile);
	tictac = strtok(data," /!#\t\n");
	PAR = atoi(tictac);
	
	while( !feof(pFile) && i_data < 40){
	  fgets(data,256,pFile);
	  tictac = strtok(data," /!#\t\n");
	  dataint[i_data] = (int)atol(tictac);
	  i_data++;
	}  
	
	fclose(pFile);
	delete [] retfile;
      }
    
  
    /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
    id.job=JOB_INIT; 
    id.par=PAR; 
    id.sym=SYM;
    id.comm_fortran=USE_COMM_WORLD;

    dmumps_c(&id);

    /* set parameter of mumps */
    if(param_int){
      if(!data_option.empty()){ 
	printf("read option before with the file %s\n",&data_option);
	exit(1);
      }
      for(int ii=0; ii<40; ii++)
	id.ICNTL(ii+1) = param_int[ii];
    }
    else 
      if(!data_option.empty()){
	for(int ii=0; ii<40; ii++)
	  id.ICNTL(ii+1) = dataint[ii];
      }
      else{
      // parameter by default
      id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
    }

    // uniquement donner au host 
    if(myid==0){
      if(perm_r && id.ICNTL(7)==1){
	for(int ii=0; ii<n; ii++) id.perm_in[ii] = pperm_r[ii];
      }
      if( scale_r && scale_c && id.ICNTL(8)==-1 ){
	// param_double[0::n-1] :: row  
	// param_double[n::n+m-1] :: column 
	for(int ii=0; ii<n; ii++) id.rowsca[ii] = scale_r[ii]; 
	for(int ii=0; ii<m; ii++) id.colsca[ii] = scale_c[ii];
      }
    }

    irn = NULL;
    jcn = NULL;

    /* valeur par defaut mis dans le dur */      
    //id.ICNTL(5)  = 0;  
    //id.ICNTL(18) = 0;
    
    // Distribution depend de la valeur de ICNTL(18)
    if( id.ICNTL(5) != 0 ){
      printf("we consider only assembled format \n");
      printf("forced assembled format ==> id.ICNTL(5)  = 0;  \n");
      id.ICNTL(5)  = 0;  
      exit(1);
    }

    if (myid == 0) {
      id.n = n; id.nz =nz; 
      //id.irn=irn; id.jcn=jcn;
      //id.a = a; //id.rhs = rhs;
    }

    if( id.ICNTL(18) == 0)
      {
	// CASE:: NON DISTRIBUTED MATRIX
	a=AA.a;
	// ATTENTION 
	// AA.cl :: indice des colonnes (exacte) et AA.lg :: indice des lignes 
	// index of row and colummn by 1
	jcn = AA.cl;
	for(int ii=0; ii<nz; ii++)
	  jcn[ii] = jcn[ii]+1;
	
	if( !(irn = (int*) malloc(sizeof(int)*nz)) ){
	  printf("problem allocation jcn ");
	  exit(1);
	}
	
	assert(AA.lg[n] == nz);
	for(int ii=0; ii< n; ii++)
	  for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
	    irn[ii1] = ii+1;  
		
	if (myid == 0) {
	  id.irn=irn; id.jcn=jcn;
	  id.a = a; //id.rhs = rhs;
	}

	/* Call the MUMPS package. */
	// Analyse + Factorisation 
	id.job=4;
	dmumps_c(&id);
	
      }


    if( id.ICNTL(18) == 1 || id.ICNTL(18) == 2 )
      {

	cout <<"id.ICNTL(18) = 1 || id.ICNTL(18) == 2 "<< endl;;
	// ATTENTION 
	// AA.cl :: indice des colonnes (exacte) et AA.lg :: indice des lignes 
	// index of row and column by 1
	jcn = AA.cl;
	for(int ii=0; ii<nz; ii++)
	  jcn[ii] = jcn[ii]+1;
		
	if( !(irn = (int*) malloc(sizeof(int)*nz)) ){
	  printf("problem allocation irn ");
	  exit(1);
	}
	
	assert(AA.lg[n] == nz);
	for(int ii=0; ii< n; ii++)
	  for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
	    irn[ii1] = ii+1;  


	if (myid == 0) {
	  id.irn=irn; id.jcn=jcn;
	}
	
	/* Call the MUMPS package. */
	// Analyse   
	id.job=1;
	dmumps_c(&id);

	if(id.ICNTL(18) == 1 ){
  
	  if(myid==0){

	    MPI_Bcast( id.mapping, nz, MPI_INT,  0, MPI_COMM_WORLD );

	    nz_loc=0;
	    for(int ii=0;ii<nz; ii++){
	      if( id.mapping[ii] == myid) nz_loc++;
	    }
	   
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc   = (double*) malloc(sizeof(double)*nz_loc);

	    int jj=0;
	    for(int ii=0;ii<nz; ii++)
	      if( id.mapping[ii] == myid){
		irn_loc[jj] = irn[ ii ];
		jcn_loc[jj] = jcn[ ii ];
		a_loc[jj] = AA.a[ ii ];
		jj++;
	      }
	    assert(jj==nz_loc);
	    
	    if(PAR==1){
	      id.nz_loc = nz_loc;
	      id.irn_loc = irn_loc;
	      id.jcn_loc = jcn_loc;
	      id.a_loc = a_loc;
	    }
	    
	  }
	  else{
	    int *mapping;
	    mapping = (int*) malloc(sizeof(int)*nz);
	    MPI_Bcast( mapping, nz, MPI_INT,  0, MPI_COMM_WORLD );
	    nz_loc=0;

	    for(int ii=0;ii<nz; ii++)
	      if( mapping[ii] == myid) nz_loc++;
	    
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc    = (double*) malloc(sizeof(double)* nz_loc);
	    
	    int jj=0.;
	    for(int ii=0;ii<nz; ii++)
	      if( mapping[ii] == myid){
		irn_loc[jj] = irn[ ii ];
		jcn_loc[jj] = jcn[ ii ];
		a_loc[jj] = AA.a[ ii ];
		jj++;
	      }
	    assert(jj==nz_loc);

	    free(mapping);

	    id.nz_loc = nz_loc;
	    id.irn_loc = irn_loc;
	    id.jcn_loc = jcn_loc;
	    id.a_loc = a_loc;



	  }
	 
	  /* Call the MUMPS package. */
	  // Factorisation   
	  id.job=2;
	  dmumps_c(&id);	  
	  
	}
      

	if(id.ICNTL(18) == 2 ){

	  if(PAR == 0){ 
	    printf("id.ICNTL(18)==2 pas prevus \n");
	    exit(1);
	    if(myid !=0) {
	      int commSize;	    
	      ierr=MPI_Comm_size(MPI_COMM_WORLD,&commSize);
	      commSize=commSize-1;
	      int myidpar=myid-1;
	      int m_loc_fst = AA.m/commSize;
	      int m_loc;
	      if( myidpar == commSize-1 && ( m_loc_fst*commSize != AA.m ) )  
		m_loc = m-m_loc_fst*( commSize-1 );
	      else
		m_loc = m_loc_fst;
	      
	      int fst_row= myidpar*m_loc_fst;
	      nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
	      
	      // allocation des tableaux
	      irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	      jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	      a_loc    = (double*) malloc(sizeof(double)*nz_loc);
	      
	      int fst_nnz = AA.lg[fst_row];
	      for(int ii=0; ii < nz_loc; ii++){
		a_loc[ii] = AA.a[fst_nnz+ii];
		jcn_loc[ii] = AA.cl[fst_nnz+ii]; // jcn=AA.cl a ete augmenter de 1 avant => pas ajouter 1
	      }
	      
	      for(int ii=fst_row; ii< fst_row+m_loc; ii++){
		for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
		  irn_loc[ii1-fst_nnz] = ii+1;  
	      }
	      
	      id.nz_loc = nz_loc;
	      id.irn_loc = irn_loc;
	      id.jcn_loc = jcn_loc;
	      id.a_loc = a_loc;
	    }
	  }
	  if(PAR == 1){

	    cout << "value of m=" << m << endl;
	    cout << "value of AA.m=" << AA.m << endl;
	    int commSize;	    
	    ierr=MPI_Comm_size(MPI_COMM_WORLD,&commSize);
	    int m_loc_fst = AA.m/commSize;
	    int m_loc;
	    if( myid == commSize-1 && ( m_loc_fst*commSize != AA.m ) )  
	      m_loc = m-m_loc_fst*( commSize-1 );
	    else
	      m_loc = m_loc_fst;
	    
	    int fst_row= myid*m_loc_fst;
	    nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
	    
	    // allocation des tableaux
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc    = (double*) malloc(sizeof(double)*nz_loc);
	    
	    int fst_nnz = AA.lg[fst_row];
	    for(int ii=0; ii < nz_loc; ii++){
	      a_loc[ii] = AA.a[fst_nnz+ii];
	      jcn_loc[ii] = AA.cl[fst_nnz+ii];  // jcn=AA.cl a ete augmenter de 1 avant => pas ajouter 1
	    }
	    
	    for(int ii=fst_row; ii< fst_row+m_loc; ii++){
	      for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
		irn_loc[ii1-fst_nnz] = ii+1;  
	    }
	  
	    id.nz_loc = nz_loc;
	    id.irn_loc = irn_loc;
	    id.jcn_loc = jcn_loc;
	    id.a_loc = a_loc;
	    
// 	    cout << "value of nz_loc=" << nz_loc << endl; 
// 	    for(int ii=0; ii < nz_loc; ii++){
// 	     cout << "a_loc[" << ii << "] =" <<  a_loc[ii] << " ";
// 	     cout << "jcn_loc[" << ii << "] =" <<  jcn_loc[ii] << " ";
// 	     cout << "irn_loc[" << ii << "] =" <<  irn_loc[ii] << endl;
// 	    }

	  }

	  /* Call the MUMPS package. */
	  // Factorisation   
	  id.job=2;
	  dmumps_c(&id);
	}
      }
	
    if( id.ICNTL(18) == 3 )
      {	
	// indices et colonnes de la matrice
	
// 	//  Cas Matrice parallele ::
// 	//  ========================
//	
// 	// Cas stockage Morse parallele
// 	m_loc = AA.m_loc;       // Nombre de lignes prise en compte
// 	nz_loc = AA.nbcoef_loc;   // Nombre de coefficients non nulles 
// 	// indice des colonnes
// 	jcn_loc = AA.cl_loc;       // indices des colonnes dans la matrice locale
//	
// 	if( !(irn_loc = (int*) malloc(sizeof(int)*nz_loc)) ){
// 	  printf("problem allocation jcn ");
// 	  exit(1);
// 	}
// 	assert(AA.lg_loc[nrow_loc] == nz_loc);
// 	for(int ii=0; ii< nrow_loc; ii++)
// 	  for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
// 	    irn_loc[ii1] = ii+1;
//
// 	a_loc=AA.a_loc;
//
	// Pas de matrice parallele ==> utilisation astuce de SuperLU
	// Matrice :: distribution bloc continue de lignes :: voir SuperLU 
	// Attention :: performance ???
	

	cout <<"id.ICNTL(18) = 3,  PAR="<< PAR << endl;

	 if(PAR == 0){ 
	    printf("id.ICNTL(18)==2 pas prevus \n");
	    exit(1);
	    if(myid !=0) {
	      int commSize;	    
	      ierr=MPI_Comm_size(MPI_COMM_WORLD,&commSize);
	      commSize=commSize-1;
	      int myidpar=myid-1;
	      int m_loc_fst;
	      m_loc_fst= AA.m/commSize;
	      int m_loc;
	      if( myidpar == commSize-1 && ( m_loc_fst*commSize != AA.m ) )  
		m_loc = m-m_loc_fst*( commSize-1 );
	      else
		m_loc = m_loc_fst;
	      
	      int fst_row;
	      fst_row= myidpar*m_loc_fst;
	      nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
	      
	      // allocation des tableaux
	      irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	      jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	      a_loc    = (double*) malloc(sizeof(double)* nz_loc);
	      
	      int fst_nnz;
	      fst_nnz = AA.lg[fst_row];
	      for(int ii=0; ii < nz_loc; ii++){
		a_loc[ii] = AA.a[fst_nnz+ii];
		jcn_loc[ii] = AA.cl[fst_nnz+ii]+1;
	      }
	      
	      for(int ii=fst_row; ii< fst_row+m_loc; ii++){
		for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
		  irn_loc[ii1-fst_nnz] = ii+1;  
	      }
	      
	      id.nz_loc = nz_loc;
	      id.irn_loc = irn_loc;
	      id.jcn_loc = jcn_loc;
	      id.a_loc = a_loc;
	    }
	 }
	 if(PAR ==1) {
	   if(verbosity > 10){ 
	     cout << "avant repar" << endl;
	     cout << "value of m=" << m << endl;
	     cout << "value of AA.m=" << AA.m << endl;
	   }
	   int commSize;
	   ierr=MPI_Comm_size(MPI_COMM_WORLD,&commSize);
	   int m_loc_fst;
	   m_loc_fst= AA.m/commSize;
	   int m_loc;
	   if( myid == commSize-1 && ( m_loc_fst*commSize != AA.m ) )  
	     m_loc = m-m_loc_fst*( commSize-1 );
	   else
	     m_loc = m_loc_fst;
	   
	   int fst_row;
	   fst_row = myid*m_loc_fst;
	   nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
	    
	   // allocation des tableaux
	   irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	   jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	   a_loc    = (double*) malloc(sizeof(double)* nz_loc);
	   
	   int fst_nnz;
	   fst_nnz= AA.lg[fst_row];
	   for(int ii=0; ii < nz_loc; ii++){
	     a_loc[ii] = AA.a[fst_nnz+ii];
	     jcn_loc[ii] = AA.cl[fst_nnz+ii]+1; 
	   }
	   
	   for(int ii=fst_row; ii< fst_row+m_loc; ii++){
	     for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
	       irn_loc[ii1-fst_nnz] = ii+1;  
	   }
	   

// 	   for(int ii=0; ii < nz_loc; ii++){
// 	     cout << "a_loc[" << ii << "] =" <<  a_loc[ii] << " ";
// 	     cout << "jcn_loc[" << ii << "] =" <<  jcn_loc[ii] << " ";
// 	     cout << "irn_loc[" << ii << "] =" <<  irn_loc[ii] << endl;
// 	   }
	   

	   id.nz_loc = nz_loc;
	   id.irn_loc = irn_loc;
	   id.jcn_loc = jcn_loc;
	   id.a_loc = a_loc;
	 }
	/* Call the MUMPS package. */
	// Analyse + Factorisation
	 
	 id.job=1;
	 dmumps_c(&id);

	 id.job=2;
	 dmumps_c(&id);
      }
    
    
    

    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if( jcn != NULL )
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]-1;

    if( irn != NULL && id.ICNTL(18) >0 ){
      free(irn); 
      irn=NULL;
    }
    if(verbosity){
      finishtime = clock();
      timeused= (finishtime-starttime)/(1000 );
      
      cout << "Factorisation with MUMPS (rank "<< myid << ") :: " << timeused << " ms" <<endl;
    }
    if(verbosity>10){
      cout << "apres repar" << endl;
      cout << "value of m=" << m << endl;
      cout << "value of AA.m=" << AA.m << endl;
    }
    
      FILE *pfile=fopen("ffmumps_fileparam2.txt","w");
      fprintf(  pfile," %d         /* SYM */\n", SYM);
      fprintf(  pfile," %d         /* PAR */\n", PAR);
      for(int ii=0; ii<40; ii++)
	fprintf(pfile," %d         /* ICNTL(%d) */\n",id.ICNTL(ii+1),ii+1); 
      fclose(pfile);

  }
  void Solver(const MatriceMorse<double> &AA,KN_<double> &x,const KN_<double> &b) const  {
	
    double *rhs;
    int job;
    
    ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
  
    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if(jcn != NULL)
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]+1;
   

    if ( !(rhs = (double*) malloc(sizeof(double)*AA.m) ) ){
      printf("Pb allocate rhs in MUMPS\n");
      exit(1);
    }
   
    for(int ii=0; ii<AA.m; ++ii){
      rhs[ii] = b[ii];
    }

    if( myid == 0 )
      id.rhs=rhs;

    /* solve linear problem */

    id.job=3;
    dmumps_c(&id);

    if( myid==0 ){
      x=id.rhs;
      MPI_Bcast( x, AA.n, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    }
    else
      MPI_Bcast( x, AA.n, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    
    // deallocation de rhs
    free(rhs);

    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if(jcn != NULL)
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]-1;
    
    if(verbosity) cout << "   x min max " << x.min() << " " <<x.max() << endl;
  }

  ~dSolveMUMPSmpi() { 
    if(verbosity)
      cout << "~SolveMUMPS S:" << endl;
    
     id.job=JOB_END; 
     dmumps_c(&id); /* Terminate instance */

     if( irn != NULL){
      free(irn); 
      irn=NULL;
     }
     /*
       free(jcn_loc);
       free(irn_loc);
       free(a_loc);
     */
  }
  void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<R> &) (*this) * x; 
  }
     
}; 


static mumps_double_complex *mumps_dc(Complex *p)  { return (mumps_double_complex *) (void *) p;}
static Complex *inv_mumps_dc(mumps_double_complex *p)  { return (Complex *) (void *) p;}

class zSolveMUMPSmpi :   public MatriceMorse<Complex>::VirtualSolver   {
  
  double eps;
  mutable double  epsr;
  double tgv;
  double tol_pivot_sym,tol_pivot; //Add 31 oct 2005

  Complex           *a;
  int       *irn, *jcn;
  int          n, m, nz; 
   
  // parameter MUMPS
  KN<int>        perm_r; /* row permutations from partial pivoting */
  KN<int>        perm_c;
  KN<double>     scale_r;
  KN<double>     scale_c;
  string  string_option;
  string data_option;
  int SYM;
  int PAR;
  int myid;

  // distribuer
  int nz_loc;
  int *jcn_loc, *irn_loc;
  Complex *a_loc;


  long int starttime,finishtime;
  long int timeused;
  
  static const int JOB_INIT=-1;
  static const int JOB_END=-2;
  static const int USE_COMM_WORLD= -987654;

  // variable complex 
  mutable ZMUMPS_STRUC_C id;
  
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
public:
  zSolveMUMPSmpi(const MatriceMorse<Complex> &AA,int strategy,double ttgv, double epsilon=1e-6,
		 double pivot=-1.,double pivot_sym=-1., string param_string, string datafile, KN<long> &param_int, 
		 KN<long> &pperm_r, KN_<long> &pperm_c, KN<double> &pscale_r,KN<double> &pscale_c) : 
    eps(epsilon),epsr(0),
    tgv(ttgv), string_option(param_string), data_option(datafile), perm_r(pperm_r), perm_c(pperm_c), 
    tol_pivot_sym(pivot_sym),tol_pivot(pivot), scale_r(pscale_r), scale_c(pscale_c)
  { 
    if(verbosity) starttime = clock();
    int dataint[40];
    int ierr;

    n    = AA.n;
    m    = AA.m; 
    nz   = AA.nbcoef;

    if(param_int) 
      assert( param_int.N() == 40);
    if(pperm_r)
      assert( perm_r.N() == n);
    if(pscale_r) 
      assert( scale_r.N() == n+m );
    if(pscale_c) 
       assert( scale_c.N() == n+m );
    if( n != m )
      cerr << "only square matrix are supported by MUMPS" << endl;
   
      /* ------------------------------------------------------------
       INITIALIZE THE MUMPS PROCESS GRID. 
       ------------------------------------------------------------*/
   
    //ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // initialisation par defaut
 
     SYM=0; PAR=1;
     if(!string_option.empty()) 
      {
	
	read_options_freefem(&string_option,&SYM,&PAR);
      }
    else 
      if(!data_option.empty())
	{
	  
	  char * retfile= new char[data_option.size()+1];
	  strcpy(retfile, (&data_option)->c_str());
	  printf("read data from file %s\n", retfile);
	  FILE *pFile=fopen(retfile,"rt");
	  
	  int     i_data=0;
	  char    data[256];
	  char   *tictac;
	  
	  fgets(data,256,pFile);
	  tictac = strtok(data," /!#\t\n");
	  SYM = atoi(tictac);
	  
	  fgets(data,256,pFile);
	  tictac = strtok(data," /!#\t\n");
	  PAR = atoi(tictac);
	  
	  while( !feof(pFile) && i_data < 40){
	    fgets(data,256,pFile);
	    tictac = strtok(data," /!#\t\n");
	    dataint[i_data] = (int)atol(tictac);
	    i_data++;
	  }  
	  
	  fclose(pFile);
	  delete [] retfile;
	}
    
   

    /* Initialize a MUMPS instance. Use MPI_COMM_WORLD */
    id.job=JOB_INIT; 
    id.par=PAR; 
    id.sym=SYM;
    id.comm_fortran=USE_COMM_WORLD;

    zmumps_c(&id);

     /* set parameter of mumps */
    if(param_int){
      if(!data_option.empty()){ 
	printf("read option before with the file %s\n",&data_option);
	exit(1);
      }
      for(int ii=0; ii<40; ii++)
	id.ICNTL(ii+1) = param_int[ii];
    }    
    else 
      if(!data_option.empty()){
	for(int ii=0; ii<40; ii++)
	  id.ICNTL(ii+1) = dataint[ii];
      }
      else{
	// parameter by default
	id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
      }

    // uniquement donner au host 
    if(myid==0){
      if(perm_r && id.ICNTL(7)==1){
	for(int ii=0; ii<n; ii++) id.perm_in[ii] = pperm_r[ii];
      }
      if( scale_r && scale_c && id.ICNTL(8)==-1 ){
	// param_double[0::n-1] :: row  
	// param_double[n::n+m-1] :: column 
	for(int ii=0; ii<n; ii++) id.rowsca[ii] = scale_r[ii]; 
	for(int ii=0; ii<m; ii++) id.colsca[ii] = scale_c[ii];
      }
    }

    
    /* valeur par defaut mis dans le dur */
   
    irn = NULL;
    jcn = NULL;


    if( id.ICNTL(5) != 0 ){
      printf("we consider only assembled format \n");
      exit(1);
    }
   
    
    /* Define the problem on the host */
    if (myid == 0) {
      id.n = n; id.nz =nz; 
      //id.irn=irn; id.jcn=jcn;
      //id.a = mumps_dc(a); //id.rhs = rhs;
    }

    if( id.ICNTL(18) == 0)
      {
	// CASE:: NON DISTRIBUTED MATRIX
	a=AA.a;
	// ATTENTION 
	// AA.cl :: indice des colonnes (exacte) et AA.lg :: indice des lignes 
	// index of row and colummn by 1
	jcn = AA.cl;
	for(int ii=0; ii<nz; ii++)
	  jcn[ii] = jcn[ii]+1;
	
	if( !(irn = (int*) malloc(sizeof(int)*nz)) ){
	  printf("problem allocation jcn ");
	  exit(1);
	}
	
	assert(AA.lg[n] == nz);
	for(int ii=0; ii< n; ii++)
	  for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
	    irn[ii1] = ii+1;  
		
	if (myid == 0) {
	  id.irn=irn; id.jcn=jcn;
	  id.a = mumps_dc(a); //id.rhs = rhs;
	}

	/* Call the MUMPS package. */
	// Analyse + Factorisation 
	id.job=4;
	zmumps_c(&id);
	
      }

    
    if( id.ICNTL(18) == 1 || id.ICNTL(18) == 2 )
      {

	cout <<"id.ICNTL(18) = 1 || id.ICNTL(18) == 2 "<< endl;;
	// ATTENTION 
	// AA.cl :: indice des colonnes (exacte) et AA.lg :: indice des lignes 
	// index of row and column by 1
	jcn = AA.cl;
	for(int ii=0; ii<nz; ii++)
	  jcn[ii] = jcn[ii]+1;
		
	if( !(irn = (int*) malloc(sizeof(int)*nz)) ){
	  printf("problem allocation irn ");
	  exit(1);
	}
	
	assert(AA.lg[n] == nz);
	for(int ii=0; ii< n; ii++)
	  for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
	    irn[ii1] = ii+1;  


	if (myid == 0) {
	  id.irn=irn; id.jcn=jcn;
	}
	
	/* Call the MUMPS package. */
	// Analyse   
	id.job=1;
	zmumps_c(&id);

	if(id.ICNTL(18) == 1 ){
	  //printf("id.ICNTL(18)==1 pas prevus a construire \n");
	  //exit(1);

	  
	  if(myid==0){
// 	    for(int ii=0;ii<nz; ii++){
// 	      cout << "proc0::  id.mapping["<< ii << "]=" << id.mapping[ii] << endl;
// 	    }
	    MPI_Bcast( id.mapping, nz, MPI_INT,  0, MPI_COMM_WORLD );

	    nz_loc=0;
	    for(int ii=0;ii<nz; ii++){
	      if( id.mapping[ii] == myid) nz_loc++;
	    }
	   
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc   = (Complex*) malloc(sizeof(Complex)*nz_loc);

	    int jj=0;
	    for(int ii=0;ii<nz; ii++)
	      if( id.mapping[ii] == myid){
		irn_loc[jj] = irn[ ii ];
		jcn_loc[jj] = jcn[ ii ];
		a_loc[jj] = AA.a[ ii ];
		jj++;
	      }
	    assert(jj==nz_loc);
	    
	    if(PAR==1){
	      id.nz_loc = nz_loc;
	      id.irn_loc = irn_loc;
	      id.jcn_loc = jcn_loc;
	      id.a_loc = mumps_dc(a_loc);
	    }
  
	  }
	  else{
	    int *mapping;
	    mapping = (int*) malloc(sizeof(int)*nz);
	    MPI_Bcast( mapping, nz, MPI_INT,  0, MPI_COMM_WORLD );
	    nz_loc=0;
// 	    for(int ii=0;ii<nz; ii++){
// 	      cout << "proc 1 id.mapping["<< ii << "]=" << mapping[ii] << endl;
// 	    }
	    for(int ii=0;ii<nz; ii++)
	      if( mapping[ii] == myid) nz_loc++;
	    
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc    = (Complex*) malloc(sizeof(Complex)* nz_loc);
	    
	    int jj=0.;
	    for(int ii=0;ii<nz; ii++)
	      if( mapping[ii] == myid){
		irn_loc[jj] = irn[ ii ];
		jcn_loc[jj] = jcn[ ii ];
		a_loc[jj] = AA.a[ ii ];
		jj++;
	      }
	    assert(jj==nz_loc);

	    free(mapping);

	    id.nz_loc = nz_loc;
	    id.irn_loc = irn_loc;
	    id.jcn_loc = jcn_loc;
	    id.a_loc = mumps_dc(a_loc);



	  }
	 
	  /* Call the MUMPS package. */
	  // Factorisation   
	  id.job=2;
	  zmumps_c(&id);	  
	  
	}
      

	if(id.ICNTL(18) == 2 ){

	  if(PAR == 0){ 
	    printf("id.ICNTL(18)==2 pas prevus \n");
	    exit(1);
	    if(myid !=0) {
	      int commSize;	    
	      ierr=MPI_Comm_size(MPI_COMM_WORLD,&commSize);
	      commSize=commSize-1;
	      int myidpar=myid-1;
	      int m_loc_fst = AA.m/commSize;
	      int m_loc;
	      if( myidpar == commSize-1 && ( m_loc_fst*commSize != AA.m ) )  
		m_loc = m-m_loc_fst*( commSize-1 );
	      else
		m_loc = m_loc_fst;
	      
	      int fst_row= myidpar*m_loc_fst;
	      nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
	      
	      // allocation des tableaux
	      irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	      jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	      a_loc    = (Complex*) malloc(sizeof(Complex)*nz_loc);
	      
	      int fst_nnz = AA.lg[fst_row];
	      for(int ii=0; ii < nz_loc; ii++){
		a_loc[ii] = AA.a[fst_nnz+ii];
		jcn_loc[ii] = AA.cl[fst_nnz+ii]; // jcn=AA.cl a ete augmenter de 1 avant => pas ajouter 1
	      }
	      
	      for(int ii=fst_row; ii< fst_row+m_loc; ii++){
		for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
		  irn_loc[ii1-fst_nnz] = ii+1;  
	      }
	      
	      id.nz_loc = nz_loc;
	      id.irn_loc = irn_loc;
	      id.jcn_loc = jcn_loc;
	      id.a_loc = mumps_dc(a_loc);
	    }
	  }
	  if(PAR == 1){

	    cout << "value of m=" << m << endl;
	    cout << "value of AA.m=" << AA.m << endl;
	    int commSize;	    
	    ierr=MPI_Comm_size(MPI_COMM_WORLD,&commSize);
	    int m_loc_fst = AA.m/commSize;
	    int m_loc;
	    if( myid == commSize-1 && ( m_loc_fst*commSize != AA.m ) )  
	      m_loc = m-m_loc_fst*( commSize-1 );
	    else
	      m_loc = m_loc_fst;
	    
	    int fst_row= myid*m_loc_fst;
	    nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
	    
	    // allocation des tableaux
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc    = (Complex*) malloc(sizeof(Complex)*nz_loc);
	    
	    int fst_nnz = AA.lg[fst_row];
	    for(int ii=0; ii < nz_loc; ii++){
	      a_loc[ii] = AA.a[fst_nnz+ii];
	      jcn_loc[ii] = AA.cl[fst_nnz+ii];  // jcn=AA.cl a ete augmenter de 1 avant => pas ajouter 1
	    }
	    
	    for(int ii=fst_row; ii< fst_row+m_loc; ii++){
	      for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
		irn_loc[ii1-fst_nnz] = ii+1;  
	    }
	  
	    id.nz_loc = nz_loc;
	    id.irn_loc = irn_loc;
	    id.jcn_loc = jcn_loc;
	    id.a_loc = mumps_dc(a_loc);

	  }

	  /* Call the MUMPS package. */
	  // Factorisation   
	  id.job=2;
	  zmumps_c(&id);
	}
      }
	
    if( id.ICNTL(18) == 3 )
      {	
	// indices et colonnes de la matrice
	
// 	//  Cas Matrice parallele ::
// 	//  ========================
//	
// 	// Cas stockage Morse parallele
// 	m_loc = AA.m_loc;       // Nombre de lignes prise en compte
// 	nz_loc = AA.nbcoef_loc;   // Nombre de coefficients non nulles 
// 	// indice des colonnes
// 	jcn_loc = AA.cl_loc;       // indices des colonnes dans la matrice locale
//	
// 	if( !(irn_loc = (int*) malloc(sizeof(int)*nz_loc)) ){
// 	  printf("problem allocation jcn ");
// 	  exit(1);
// 	}
// 	assert(AA.lg_loc[nrow_loc] == nz_loc);
// 	for(int ii=0; ii< nrow_loc; ii++)
// 	  for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
// 	    irn_loc[ii1] = ii+1;
//
// 	a_loc=AA.a_loc;
//
	// Pas de matrice parallele ==> utilisation astuce de SuperLU
	// Matrice :: distribution bloc continue de lignes :: voir SuperLU 
	// Attention :: performance ???
	

	cout <<"id.ICNTL(18) = 3,  PAR="<< PAR << endl;

	 if(PAR == 0){ 
	    if(myid !=0) {
	      int commSize;	    
	      ierr=MPI_Comm_size(MPI_COMM_WORLD,&commSize);
	      commSize=commSize-1;
	      int myidpar=myid-1;
	      int m_loc_fst;
	      m_loc_fst= AA.m/commSize;
	      int m_loc;
	      if( myidpar == commSize-1 && ( m_loc_fst*commSize != AA.m ) )  
		m_loc = m-m_loc_fst*( commSize-1 );
	      else
		m_loc = m_loc_fst;
	      
	      int fst_row;
	      fst_row= myidpar*m_loc_fst;
	      nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
	      
	      // allocation des tableaux
	      irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	      jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	      a_loc    = (Complex*) malloc(sizeof(Complex)* nz_loc);
	      
	      int fst_nnz;
	      fst_nnz = AA.lg[fst_row];
	      for(int ii=0; ii < nz_loc; ii++){
		a_loc[ii] = AA.a[fst_nnz+ii];
		jcn_loc[ii] = AA.cl[fst_nnz+ii]+1;
	      }
	      
	      for(int ii=fst_row; ii< fst_row+m_loc; ii++){
		for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
		  irn_loc[ii1-fst_nnz] = ii+1;  
	      }
	      
	      id.nz_loc = nz_loc;
	      id.irn_loc = irn_loc;
	      id.jcn_loc = jcn_loc;
	      id.a_loc = mumps_dc(a_loc);
	    }
	 }
	 if(PAR ==1) {
	   if(verbosity > 10){ 
	     cout << "avant repar" << endl;
	     cout << "value of m=" << m << endl;
	     cout << "value of AA.m=" << AA.m << endl;
	   }
	   int commSize;
	   ierr=MPI_Comm_size(MPI_COMM_WORLD,&commSize);
	   int m_loc_fst;
	   m_loc_fst= AA.m/commSize;
	   int m_loc;
	   if( myid == commSize-1 && ( m_loc_fst*commSize != AA.m ) )  
	     m_loc = m-m_loc_fst*( commSize-1 );
	   else
	     m_loc = m_loc_fst;
	   
	   int fst_row;
	   fst_row = myid*m_loc_fst;
	   nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
	    
	   // allocation des tableaux
	   irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	   jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	   a_loc    = (Complex*) malloc(sizeof(Complex)*nz_loc);
	   
	   int fst_nnz;
	   fst_nnz= AA.lg[fst_row];
	   for(int ii=0; ii < nz_loc; ii++){
	     a_loc[ii] = AA.a[fst_nnz+ii];
	     jcn_loc[ii] = AA.cl[fst_nnz+ii]+1; 
	   }
	   
	   for(int ii=fst_row; ii< fst_row+m_loc; ii++){
	     for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
	       irn_loc[ii1-fst_nnz] = ii+1;  
	   }
	   
	   id.nz_loc = nz_loc;
	   id.irn_loc = irn_loc;
	   id.jcn_loc = jcn_loc;
	   id.a_loc = mumps_dc(a_loc);
	 }
	/* Call the MUMPS package. */
	// Analyse + Factorisation
	 
	 id.job=1;
	 zmumps_c(&id);

	 id.job=2;
	 zmumps_c(&id);
      }
    
    
    

    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if( jcn != NULL )
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]-1;

    if( irn != NULL ) free(irn); 

    if(verbosity){
      finishtime = clock();
      timeused= (finishtime-starttime)/(1000 );
      
      cout << "Factorisation with MUMPS (rank "<< myid << ") :: " << timeused << " ms" <<endl;
    }

        
  
  }

  void Solver(const MatriceMorse<Complex> &AA,KN_<Complex> &x,const KN_<Complex> &b) const  {
    //*******************************************************************//
    //    depend pas de la forme de la matrice: distribuer ou assembler
    Complex *rhs;
    int job;
    
    ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
  
    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if(jcn != NULL)
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]+1;
    
    if ( !(rhs = (Complex*) malloc(sizeof(Complex)*AA.m) ) ){
      printf("Pb allocate rhs in MUMPS\n");
      exit(1);
    }
   
    for(int ii=0; ii<AA.m; ++ii){
      rhs[ii] = b[ii];
    }

    if( myid == 0 )
      id.rhs=mumps_dc(rhs);

    /* solve linear problem */
    id.job=3;
    zmumps_c(&id);

   
    if( myid==0 ){
      x=inv_mumps_dc(id.rhs); 
      MPI_Bcast( x,  AA.n, MPI_DOUBLE_COMPLEX,  0, MPI_COMM_WORLD );
    }
    else
      MPI_Bcast( x,  AA.n, MPI_DOUBLE_COMPLEX,  0, MPI_COMM_WORLD );
    
  
    
    // deallocation de rhs
    free(rhs);

    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if(jcn != NULL)
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]-1;
    
   
  }

  ~zSolveMUMPSmpi() { 
    //*******************************************************************//
    //    depend pas de la forme de la matrice: distribuer ou assembler
    if(verbosity)
      cout << "~SolveMUMPS Z:" << endl;
    
     id.job=JOB_END; 
     zmumps_c(&id); /* Terminate instance */
     
  }
  void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<Complex> &) (*this) * x; 
  }
     
}; 


MatriceMorse<double>::VirtualSolver *
BuildSolverMUMPSmpi(DCL_ARG_SPARSE_SOLVER(double,A))
{
    if(verbosity>9)
      cout << " BuildSolverMUMPSmpi<double>" << endl;
    return new dSolveMUMPSmpi(*A,ds.strategy, ds.tgv, ds.epsilon, ds.tol_pivot, ds.tol_pivot_sym, ds.sparams, ds.data_filename,
			      ds.lparams, ds.perm_r, ds.perm_c, ds.scale_r, ds.scale_c);
}




MatriceMorse<Complex>::VirtualSolver *
BuildSolverMUMPSmpi(DCL_ARG_SPARSE_SOLVER(Complex,A))
{
    if(verbosity>9)
      cout << " BuildSolverMUMPSmpi<Complex>" << endl;
    return new zSolveMUMPSmpi(*A,ds.strategy, ds.tgv, ds.epsilon, ds.tol_pivot, ds.tol_pivot_sym, ds.sparams, ds.data_filename,  
			      ds.lparams, ds.perm_r, ds.perm_c, ds.scale_r, ds.scale_c);
}


class Init { public:
    Init();
};

//  the 2 default sparse solver double and complex
DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; 
DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetDefault()
{
    if(verbosity)
	cout << " SetDefault sparse to default" << endl;
    DefSparseSolver<double>::solver =SparseMatSolver_R;
    DefSparseSolver<Complex>::solver =SparseMatSolver_C;
    TypeSolveMat::defaultvalue =TypeSolveMat::SparseSolver;
}

bool SetMUMPSmpi()
{
    if(verbosity)
	cout << " SetDefault sparse solver to MUMPSmpi" << endl;
    DefSparseSolver<double>::solver  =BuildSolverMUMPSmpi;
    DefSparseSolver<Complex>::solver =BuildSolverMUMPSmpi;    
    TypeSolveMat::defaultvalue  = TypeSolveMatdefaultvalue;
}



Init init;
Init::Init()
{ 
  
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  
  if(verbosity>1)
    cout << "\n Add: MUMPSmpi,  defaultsolver defaultsolverMUMPSmpi" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  DefSparseSolver<double>::solver =BuildSolverMUMPSmpi;
  DefSparseSolver<Complex>::solver =BuildSolverMUMPSmpi;
  if(! Global.Find("defaultsolver").NotNull() )
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  Global.Add("defaulttoMUMPSmpi","(",new OneOperator0<bool>(SetMUMPSmpi));
}

