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
  mutable MPI_Comm comm;
  double            *a;
  int       *irn, *jcn;
  int          n, m, nz; 
 
  // parameter MUMPS
  
  KN_<long>        perm_r; /* row permutations from partial pivoting */
  KN_<long>        perm_c;
  KN_<double>     scale_r;
  KN_<double>     scale_c;
  string string_option;
  string data_option;
  int SYM;
  int PAR;
  int myid;

  // distribuer
  int nz_loc;
  int *jcn_loc, *irn_loc;
  double *a_loc;
  


  static const int JOB_INIT=-1;
  static const int JOB_END=-2;
  static const int USE_COMM_WORLD= -987654;

  // variable reel
  mutable DMUMPS_STRUC_C id;
  
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define CNTL(I) cntl[(I)-1] /* macro s.t. indices match documentation */
#define RINFOG(I) rinfog[(I)-1] /* macro s.t. indices match documentation */
#define INFOG(I) infog[(I)-1] /* macro s.t. indices match documentation */
public:
  dSolveMUMPSmpi(const MatriceMorse<double> &AA,int strategy,double ttgv, double epsilon,
		 double pivot,double pivot_sym, string param_string, string datafile, KN<long> &param_int, 
		 KN<double> &param_double, KN<long> &pperm_r, KN<long> &pperm_c, KN<double> &pscale_r,KN<double> &pscale_c, MPI_Comm  * mpicommw
		) : 
    eps(epsilon),epsr(0),
    tgv(ttgv), string_option(param_string), data_option(datafile), perm_r(pperm_r), perm_c(pperm_c), 
    tol_pivot_sym(pivot_sym),tol_pivot(pivot), scale_r(pscale_r), scale_c(pscale_c)
  { 
    long int starttime,finishtime;
    long int timeused;

    if(verbosity) starttime = clock();
    int dataint[40];
    double datadouble[15];
    int ierr;

     if(mpicommw==0){
	comm=MPI_COMM_WORLD;
	}
	else
	comm= *mpicommw;


    /* ------------------------------------------------------------
       INITIALIZE THE MUMPS PROCESS GRID. 
       ------------------------------------------------------------*/
    ierr = MPI_Comm_rank(comm, &myid);
    
    if( myid ==0){
      n    = AA.n;
      m    = AA.m; 
      nz   = AA.nbcoef;
    
      MPI_Bcast(  &n, 1, MPI_INT,  0, comm );
      MPI_Bcast(  &m, 1, MPI_INT,  0, comm );
      MPI_Bcast( &nz, 1, MPI_INT,  0, comm );
    }
    else{
      MPI_Bcast(  &n, 1, MPI_INT,  0, comm );
      MPI_Bcast(  &m, 1, MPI_INT,  0, comm );
      MPI_Bcast( &nz, 1, MPI_INT,  0, comm );
    }
    
    if( !(param_int==NULL) ) 
      assert( param_int.N() == 42);
    if( !(param_double==NULL) ) 
      assert( param_double.N() == 15);
    if(perm_r)
      assert( perm_r.N() == n);
    if(perm_c)
      assert( perm_c.N() == m);
    if(scale_r) 
      assert( scale_r.N() == n);
    if(scale_c) 
      assert( scale_c.N() == m);

    if( n != m )
      cerr << "only square matrix are supported by MUMPS" << endl;
   
  
    // initialisation par defaut
    SYM=0; PAR=1;
    
    /*
      if(!string_option.empty()) 
      {
      if(myid==0){
      cout << "read string option" <<endl;
      read_options_freefem(&string_option,&SYM,&PAR);
      
      MPI_Bcast(  &SYM, 1, MPI_INT,  0, comm );
      MPI_Bcast(  &PAR, 1, MPI_INT,  0, comm );
      }
      else{
      MPI_Bcast(  &SYM, 1, MPI_INT,  0, comm );
      MPI_Bcast(  &PAR, 1, MPI_INT,  0, comm );
      }
      }
    */
    if( !(param_int==NULL) ){
      SYM = param_int[0];
      PAR = param_int[1];
      cout << "param :: myid =" << myid << endl;
    }
    else if( !data_option.empty() )
      {	
	cout << "myid =" << myid << endl;
	if(myid==0){
	  
	  char * retfile= new char[data_option.size()+1];
	  strcpy(retfile, (&data_option)->c_str());
	  printf("read data option file %s\n",retfile);
	  FILE *pFile=fopen(retfile,"rt");
	
	  int     i_data=0;
	  int     d_data=0;
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
	  assert(i_data == 40);
	  while( !feof(pFile) && d_data < 15){
	    fgets(data,256,pFile);
	    tictac = strtok(data," /!#\t\n");
	    datadouble[d_data] = (double)atof(tictac);
	    d_data++;
	  }  
	  assert(d_data == 15);
	  /*
	    for(int ii=0; ii< 40; ii++){
	    cout << "double int["<< ii <<"] ="<< dataint[ii] << endl;
	    }  
	    for(int ii=0; ii< 15; ii++){
	    cout << "double data["<< ii <<"] ="<< datadouble[ii] << endl;
	    }  
	  */
	  MPI_Bcast(  &SYM, 1, MPI_INT,  0, comm );
	  MPI_Bcast(  &PAR, 1, MPI_INT,  0, comm );

	  cout << "myid =" << myid << " init parameter :: PAR & SYM " << PAR << " " << SYM << endl; 

	  MPI_Bcast(  dataint, 40, MPI_INT,  0, comm );
	  MPI_Bcast(  datadouble, 15, MPI_DOUBLE,  0, comm );


	  fclose(pFile);
	  delete [] retfile;
	}
	else{
	  

	  MPI_Bcast(  &SYM, 1, MPI_INT,  0, comm );
	  MPI_Bcast(  &PAR, 1, MPI_INT,  0, comm );
	  
	  cout << "myid =" << myid << "  init parameter :: PAR & SYM " << PAR << " " << SYM << endl; 


	  MPI_Bcast(  dataint, 40, MPI_INT,  0, comm );
	  MPI_Bcast(  datadouble, 15, MPI_DOUBLE,  0, comm );
	}
      }
    
  
    /* Initialize a MUMPS instance. Use comm */
    id.job=JOB_INIT; 
    id.par=PAR; 
    id.sym=SYM;
    id.comm_fortran= (MUMPS_INT) MPI_Comm_c2f( comm );
    //id.comm_fortran= (F_INT) comm;

    if(verbosity) cout << "init parameter :: PAR & SYM " << PAR << " " << SYM << endl; 

    dmumps_c(&id);

    if(verbosity) cout << "fin init parameter" << endl; 

    /* set parameter of mumps */
    if( !(param_int == NULL) || !(param_double == NULL) ){
      if(!data_option.empty()){ 
	printf("MUMPS ERROR:  parameters are given on the file %s and in the array lparams and dparams => double definition of parameters.",&data_option);
	exit(1);
      }
      
      if( !(param_int == NULL) ){
	cout << "internal parameter" << endl;
	for(int ii=0; ii<40; ii++)	  
	  id.ICNTL(ii+1) = param_int[ii+2];
      }
      else{

	// parameter by default
	id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
      }
      if( !(param_double == NULL) ){
	cout << "internal parameter" << endl;
	for(int ii=0; ii<15; ii++)
	  id.CNTL(ii+1) = param_double[ii];
      }
    }    
    else 
      if(!data_option.empty()){
	for(int ii=0; ii<40; ii++)
	  id.ICNTL(ii+1) = dataint[ii];
	for(int ii=0; ii<15; ii++)
	  id.CNTL(ii+1) = datadouble[ii];
      }
      else{
	// parameter by default
	cout << "default parameter" << endl;
	id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
    }

    // uniquement donner au host 
    if(myid==0){
      if( !(perm_r==NULL) && id.ICNTL(7)==1){
	for(int ii=0; ii<n; ii++) id.perm_in[ii] = pperm_r[ii];
      }
      // a decommenter
      //if( !(perm_c==NULL) && id.ICNTL(6)==1){
      //for(int ii=0; ii<m; ii++) id.perm_in[ii] = pperm_c[ii];
      //}
      if( !(scale_r==NULL) && !(scale_c==NULL) && id.ICNTL(8)==-1 ){
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
	if (myid == 0) { // nouveau

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
		
	  // if (myid == 0) {   // ancien
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

	if(verbosity > 1) cout <<"id.ICNTL(18) = 1 || id.ICNTL(18) == 2 "<< endl;
	// ATTENTION 
	// AA.cl :: indice des colonnes (exacte) et AA.lg :: indice des lignes 
	// index of row and column by 1
	
	if (myid == 0) { // new host process 
	
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
	  
	  
	  //if (myid == 0) {  // enlever new host
	  id.irn=irn; id.jcn=jcn;
	}
	
	/* Call the MUMPS package. */
	// Analyse   
	id.job=1;
	dmumps_c(&id);

	if(id.ICNTL(18) == 1 ){
	  
	   if( PAR == 0 ){ 
	    int *nz_loc_procs;
	    int *fst_nz_procs;
	    int *irn_g;
	    int *jcn_g;
	    double *a_g;
	    int commSize;
	    
	    MPI_Comm_size(comm,&commSize);
	    cout << commSize << "commSize" << "nz =" << nz << endl;
	    if(myid==0){
	      // allocation des differents tableaux
	      nz_loc_procs = (int*) malloc ( commSize*sizeof(int) );
	      for(int ii=0; ii<commSize; ii++){
		nz_loc_procs[ii]=0;
	      }
	      for(int ii=0; ii<nz; ii++){
		nz_loc_procs[ id.mapping[ii] ]++;
	      }
	      assert(nz_loc_procs[0] == 0);
	      nz_loc_procs[0] = 2;

	      fst_nz_procs = (int*) malloc ( commSize*sizeof(int) );	      
	      fst_nz_procs[0] = 0;
	      for(int ii=1; ii<commSize; ii++){
		fst_nz_procs[ii] = fst_nz_procs[ii-1]+nz_loc_procs[ii-1];
	      }
	      
	      irn_g = (int*) malloc( sizeof(int)*(nz+2) );
	      jcn_g = (int*) malloc( sizeof(int)*(nz+2) );
	      a_g   = (double*) malloc( sizeof(double)*(nz+2) );
	      
	      int *index_p;
	      index_p = (int*) malloc ( commSize*sizeof(int) );
	      for(int ii=0; ii<commSize; ii++)
		index_p[ii] =0;
	      
	      irn_g[ 0 ] = 1;
	      jcn_g[ 0 ] = 1;
	      a_g  [ 0 ] = 1.;
	      
	      irn_g[ 1 ] = 1;
	      jcn_g[ 1 ] = 1;
	      a_g  [ 1 ] = 1.;
	      
	      for(int ii=0;ii<nz; ii++){	      
		int jj1 = id.mapping[ii];
		int jj2 = fst_nz_procs[jj1] + index_p[jj1];
		assert(jj2 > 1);
		irn_g[ jj2 ] =  irn[ ii ];
		jcn_g[ jj2 ] =  jcn[ ii ];
		a_g  [ jj2 ] = AA.a[ ii ];
		cout << "jj2= " << jj2 << endl;
		assert( jj2 < nz+2);
		index_p[jj1]++;
	      }
	      free(index_p);
	      
	    }
	    
	    MPI_Scatter( nz_loc_procs, 1, MPI_INT, &nz_loc, 1, MPI_INT, 0, comm);
	    
	    // allocation des tableaux locaux
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc    = (double*) malloc(sizeof(double)*nz_loc);
	    
	    MPI_Scatterv(   a_g, nz_loc_procs, fst_nz_procs, MPI_DOUBLE, a_loc, nz_loc, MPI_DOUBLE, 0, comm);
	    MPI_Scatterv( jcn_g, nz_loc_procs, fst_nz_procs, MPI_INT, jcn_loc, nz_loc, MPI_INT, 0, comm);
	    MPI_Scatterv( irn_g, nz_loc_procs, fst_nz_procs, MPI_INT, irn_loc, nz_loc, MPI_INT, 0, comm);
	    cout << "myid=" << myid <<" nz_loc=" << nz_loc << endl;
	    if( myid == 1){
	      cout << "nz_loc=" << nz_loc << endl;
	      for(int ii=0;ii<nz_loc; ii++){
		cout << "a_loc[ii]" << a_loc[ii] << endl;
	
	      }
	    }
	    
	    if( myid > 0){
	      id.nz_loc = nz_loc;
	      id.irn_loc = irn_loc;
	      id.jcn_loc = jcn_loc;
	      id.a_loc = a_loc;
	    }
	    
	    if( myid == 0){
	      //free( irn_loc );
	      //free( jcn_loc );
	      //free( a_loc );
	      free( nz_loc_procs );
	      free( fst_nz_procs );
	      free( irn_g );
	      free( jcn_g );
	      free( a_g );
	    }
	   }
	
	  
	  if( PAR == 1 ){ 
	    int *nz_loc_procs;
	    int *fst_nz_procs;
	    int *irn_g;
	    int *jcn_g;
	    double *a_g;
	    int commSize;
	    
	    MPI_Comm_size(comm,&commSize);
	    
	    if(myid==0){
	      // allocation des differents tableaux
	      nz_loc_procs = (int*) malloc ( commSize*sizeof(int) );
	      for(int ii=0; ii<commSize; ii++){
		nz_loc_procs[ii]=0;
	      }
	      for(int ii=0; ii<nz; ii++){
		nz_loc_procs[ id.mapping[ii] ]++;
	      }
	      
	      fst_nz_procs = (int*) malloc ( commSize*sizeof(int) );
	      
	      fst_nz_procs[0] = 0;
	      for(int ii=1; ii<commSize; ii++){
		fst_nz_procs[ii] = fst_nz_procs[ii-1]+nz_loc_procs[ii-1];
	      }
	      
	      irn_g = (int*) malloc(sizeof(int)*nz);
	      jcn_g = (int*) malloc(sizeof(int)*nz);
	      a_g   = (double*) malloc(sizeof(double)*nz);
	      
	      int *index_p;
	      index_p = (int*) malloc ( commSize*sizeof(int) );
	      for(int ii=0; ii<commSize; ii++)
		index_p[ii] =0;
	      
	      for(int ii=0;ii<nz; ii++){	      
		int jj1 = id.mapping[ii];
		int jj2 = fst_nz_procs[jj1] + index_p[jj1];

		irn_g[ jj2 ] =  irn[ ii ];
		jcn_g[ jj2 ] =  jcn[ ii ];
		a_g  [ jj2 ] = AA.a[ ii ];
		index_p[jj1]++;
	      }
	      free(index_p);
	      
	    }
	    
	    MPI_Scatter( nz_loc_procs, 1, MPI_INT, &nz_loc, 1, MPI_INT, 0, comm);
	    
	    // allocation des tableaux locaux
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc    = (double*) malloc(sizeof(double)*nz_loc);
	    
	    MPI_Scatterv(   a_g, nz_loc_procs, fst_nz_procs, MPI_DOUBLE, a_loc, nz_loc, MPI_DOUBLE, 0, comm);
	    MPI_Scatterv( jcn_g, nz_loc_procs, fst_nz_procs, MPI_INT, jcn_loc, nz_loc, MPI_INT, 0, comm);
	    MPI_Scatterv( irn_g, nz_loc_procs, fst_nz_procs, MPI_INT, irn_loc, nz_loc, MPI_INT, 0, comm);
	    
	    id.nz_loc = nz_loc;
	    id.irn_loc = irn_loc;
	    id.jcn_loc = jcn_loc;
	    id.a_loc = a_loc;
	    
	    
	    if( myid == 0){
	      free( nz_loc_procs );
	      free( fst_nz_procs );
	      free( irn_g );
	      free( jcn_g );
	      free( a_g );
	    }
	  }
	  // version all procs
// 	    if(myid==0){
//	    
// 	    MPI_Bcast( id.mapping, nz, MPI_INT,  0, comm );
//    
// 	    nz_loc=0;
// 	    for(int ii=0;ii<nz; ii++){
// 	      if( id.mapping[ii] == myid) nz_loc++;
// 	    }
//   
// 	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	    a_loc   = (double*) malloc(sizeof(double)*nz_loc);
//
// 	    int jj=0;
// 	    for(int ii=0;ii<nz; ii++)
// 	      if( id.mapping[ii] == myid){
// 		irn_loc[jj] = irn[ ii ];
// 		jcn_loc[jj] = jcn[ ii ];
// 		a_loc[jj] = AA.a[ ii ];
// 		jj++;
// 	      }
// 	    assert(jj==nz_loc);
//	    
// 	    if(PAR==1){
// 	      id.nz_loc = nz_loc;
// 	      id.irn_loc = irn_loc;
// 	      id.jcn_loc = jcn_loc;
// 	      id.a_loc = a_loc;
// 	    }
//	    
// 	  }
// 	  else{
// 	    int *mapping;
// 	    mapping = (int*) malloc(sizeof(int)*nz);
// 	    MPI_Bcast( mapping, nz, MPI_INT,  0, comm );
// 	    nz_loc=0;
//
// 	    for(int ii=0;ii<nz; ii++)
// 	      if( mapping[ii] == myid) nz_loc++;
//	    
// 	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	    a_loc    = (double*) malloc(sizeof(double)* nz_loc);
//	    
// 	    //==============================
// 	    //    besoin de irn, jcn, AA.a
// 	    //==============================
//
// 	    int jj=0.;
// 	    for(int ii=0;ii<nz; ii++)
// 	      if( mapping[ii] == myid){
// 		irn_loc[jj] = irn[ ii ];  
// 		jcn_loc[jj] = jcn[ ii ];
// 		a_loc[jj] = AA.a[ ii ];
// 		jj++;
// 	      }
// 	    assert(jj==nz_loc);
//
// 	    free(mapping);
//
// 	    id.nz_loc = nz_loc;
// 	    id.irn_loc = irn_loc;
// 	    id.jcn_loc = jcn_loc;
// 	    id.a_loc = a_loc;
// 	  }
	 
	  /* Call the MUMPS package. */
	  // Factorisation   
	  id.job=2;
	  dmumps_c(&id);	  
	  
	}
      

	if(id.ICNTL(18) == 2 ){
	  printf("id.ICNTL(18)==2 not avaible yet \n");
	  exit(1);


	  if(PAR == 0){ 
	    printf("id.ICNTL(18)==2 with PAR=0 not available yet \n");
	    exit(1);
	    if(myid !=0) {
	      
	      //==============================
	      //    besoin de irn, jcn, AA.a
	      //==============================
	      
	      int commSize;	    
	      ierr=MPI_Comm_size(comm,&commSize);
	      commSize=commSize-1;
	      int myidpar=myid-1;
	      int m_loc_fst = m/commSize;
	      int m_loc;
	      if( myidpar == commSize-1 && ( m_loc_fst*commSize != m ) )  
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

	    //==============================
	    //    besoin de irn, jcn, AA.a
	    //==============================
	      
	    int commSize;	    
	    ierr=MPI_Comm_size(comm,&commSize);
	    int m_loc_fst = m/commSize;
	    int m_loc;
	    if( myid == commSize-1 && ( m_loc_fst*commSize != m ) )  
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
	
// 	    //==================
// 	    // pour un proc :
// 	    // a:     AA.a[fstow]   , ..., AA.a[fstrow+nz_loc-1] 
// 	    // jcn:   AA.cl[fstrow] , ..., AA.cl[fstrow+nz_loc-1]
// 	    // irn:   AA.lg[fstrow] , ..., AA.lg[fstrow+m_loc] 
// 	    //==================
// 	    // apres reception : 
// 	    // irn_reel:  
// 	    //  int jj=0, 
// 	    //  do ii=0,m_loc-1
// 	    //      do ii1=irn_donnee[ii],irn_donnee[ii+1]-1
// 	    //         irn_reel[jj] = fst_row+ii+1;
// 	    //         jj++
// 	    //      end do
// 	    //  end do
// 	    //=================

	cout <<"id.ICNTL(18) = 3,  PAR="<< PAR << endl;

	if(PAR == 0){
	  	  
	    
// 	  if(myid != 0) {
// 	    int commSize;	    
// 	    ierr=MPI_Comm_size(comm,&commSize);
// 	    commSize=commSize-1;
// 	    int myidpar=myid-1;
// 	    int m_loc_fst;
// 	    m_loc_fst= m/commSize;
// 	    int m_loc;
// 	    if( myidpar == commSize-1 && ( m_loc_fst*commSize != m ) )  
// 	      m_loc = m-m_loc_fst*( commSize-1 );
// 	    else
// 	      m_loc = m_loc_fst;
	    
// 	    int fst_row;
// 	    fst_row= myidpar*m_loc_fst;
// 	    nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
	    
// 	    // allocation des tableaux
// 	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	    a_loc    = (double*) malloc(sizeof(double)* nz_loc);
	    
// 	    int fst_nnz;
// 	    fst_nnz = AA.lg[fst_row];
// 	    for(int ii=0; ii < nz_loc; ii++){
// 	      a_loc[ii] = AA.a[fst_nnz+ii];
// 	      jcn_loc[ii] = AA.cl[fst_nnz+ii]+1;
// 	    }
	    
// 	    for(int ii=fst_row; ii< fst_row+m_loc; ii++){
// 	      for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
// 		irn_loc[ii1-fst_nnz] = ii+1;  
// 	    }
	 
// 	    id.nz_loc = nz_loc;
// 	    id.irn_loc = irn_loc;
// 	    id.jcn_loc = jcn_loc;
// 	    id.a_loc = a_loc;
// 	  }


	  // definition de variables
	  int commSize;
	  int m_loc_fst;	 
	  int m_loc;
	  int fst_row;
	  

	  int *nz_loc_procs;
	  int *fst_nz_procs;
	  int *m_loc_procs;
	  int *fst_row_procs;
	 

	  double *tab_a;
	  int *tab_cl;
	  int *tab_lg;
	  int *tab_lg_loc;

	  MPI_Comm_size(comm,&commSize);
	  
	  if( myid !=0){
	    int commSizemm;
	    int myidpar=myid-1;

	    commSizemm = commSize-1;
	    m_loc_fst= m/commSizemm;
	  
	    if( myidpar == commSizemm-1 && ( m_loc_fst*commSizemm != m ) )  
	      m_loc = m-m_loc_fst*( commSizemm-1 );
	    else
	      m_loc = m_loc_fst;
	  
	    if(verbosity > 5){
	      fst_row = myidpar*m_loc_fst;
	      cout << "   myid = " << myid << endl;
	      cout <<"   m_loc = " << m_loc << endl;
	      cout <<" fst_row = " << fst_row << endl;
	    }

	  }
	  if( myid ==0){

	    int commSizemm;
	    commSizemm = commSize-1;
	    m_loc_fst= m/commSizemm;

	    fst_row_procs = (int* ) malloc( commSize*sizeof(int) );
	    m_loc_procs = (int* ) malloc( commSize*sizeof(int) );
	    fst_nz_procs = (int* ) malloc( commSize*sizeof(int) );
	    nz_loc_procs = (int* ) malloc ( commSize*sizeof(int) );
	    
	    
	    fst_row_procs [0] = 0;
	    m_loc_procs   [0] = 0;

	    for( int ii= 1; ii<commSize; ii++){
	      fst_row_procs [ii] = (ii-1)*m_loc_fst;
	      m_loc_procs [ii] = m_loc_fst;
	    }
	    
	    if( m_loc_fst*(commSize-1) != m ) 
	      m_loc_procs [commSize-1] = m-m_loc_fst*( (commSize-1)-1 );


	    nz_loc_procs [0] = 0;
	    fst_nz_procs [0] = 0;
	    
	    for( int ii= 1; ii<commSize; ii++){
	      nz_loc_procs [ii] = AA.lg[fst_row_procs[ii]+m_loc_procs[ii] ] - AA.lg[fst_row_procs[ii]];
	      fst_nz_procs [ii] = AA.lg[fst_row_procs[ii]];
	    }

	   
	    /*
	      tab_a= (int* ) malloc( nz*sizeof(double) );
	      tab_cl = (int* ) malloc( nz*sizeof(int) );
	      tab_lg = (int* ) malloc ( n*sizeof(int) );
	    */
	    tab_a  = AA.a;
	    tab_cl = AA.cl;
	    tab_lg = AA.lg;
	  }

	  MPI_Scatter( nz_loc_procs, 1, MPI_INT, &nz_loc, 1, MPI_INT, 0, comm);
	  MPI_Scatter( m_loc_procs,  1, MPI_INT, &m_loc, 1, MPI_INT, 0, comm);
	  MPI_Scatter( fst_row_procs,  1, MPI_INT, &fst_row, 1, MPI_INT, 0, comm);

	  if(verbosity > 5){
	    cout << "after scatter " << myid << endl;
	    cout << "   myid = " << myid << endl;
	    cout <<"   m_loc = " << m_loc << endl;
	    cout <<" fst_row = " << fst_row << endl;
	  }
	  // allocation des tableaux locaux
	  irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	  jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	  a_loc    = (double*) malloc(sizeof(double)*nz_loc);
	  tab_lg_loc = (int*) malloc(sizeof(int)*(m_loc) );
	  
	

	  MPI_Scatterv(  tab_a, nz_loc_procs, fst_nz_procs, MPI_DOUBLE, a_loc, nz_loc, MPI_DOUBLE, 0, comm);
	  MPI_Scatterv( tab_cl, nz_loc_procs, fst_nz_procs, MPI_INT, jcn_loc, nz_loc, MPI_INT, 0, comm);
	  MPI_Scatterv( tab_lg,  m_loc_procs, fst_row_procs, MPI_INT, tab_lg_loc, m_loc, MPI_INT, 0, comm);
	  
	
	  int jj=0;
	  for(int  ii=0; ii<m_loc-1; ii++)
	    for(int ii1= tab_lg_loc[ii]; ii1 < tab_lg_loc[ii+1]; ii1++){
	      irn_loc[jj] = fst_row+ii+1;
	      jj++;
	    }
	  
	  for(int ii1= tab_lg_loc[m_loc-1]; ii1 < tab_lg_loc[0]+nz_loc; ii1++){
	    irn_loc[jj] = fst_row+(m_loc-1)+1;
	    jj++;
	  }
	  
	  for(int ii=0; ii < nz_loc; ii++){	    
	    jcn_loc[ii] = jcn_loc[ii]+1; 
	  }

	  assert( jj == nz_loc );
	  
	  free( tab_lg_loc );
	  
	  id.nz_loc = nz_loc;
	  id.irn_loc = irn_loc;
	  id.jcn_loc = jcn_loc;
	  id.a_loc = a_loc;
	 
	  if( myid == 0 ){
	    free( fst_row_procs );
	    free( m_loc_procs );
	    free( fst_nz_procs );
	    free( nz_loc_procs );
	  }

	}
	if(PAR ==1) {
	  /*
	    int commSize;
	    ierr=MPI_Comm_size(comm,&commSize);
	    int m_loc_fst;
	    m_loc_fst= m/commSize;
	    int m_loc;
	    if( myid == commSize-1 && ( m_loc_fst*commSize != m ) )  
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
	    for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++)
	    irn_loc[ii1-fst_nnz] = ii+1;  
	    }
	    
	  */

	  // definition de variables
	  int commSize;
	  int m_loc_fst;	 
	  int m_loc;
	  int fst_row;
	  

	  int *nz_loc_procs;
	  int *fst_nz_procs;
	  int *m_loc_procs;
	  int *fst_row_procs;
	 

	  double *tab_a;
	  int *tab_cl;
	  int *tab_lg;
	  int *tab_lg_loc;

	  MPI_Comm_size(comm,&commSize);
	  m_loc_fst= m/commSize;
	  
	  if( myid == commSize-1 && ( m_loc_fst*commSize != m ) )  
	    m_loc = m-m_loc_fst*( commSize-1 );
	  else
	    m_loc = m_loc_fst;	  
	  
	  fst_row = myid*m_loc_fst;
	  
	  if( myid ==0){
	    fst_row_procs = (int* ) malloc( commSize*sizeof(int) );
	    m_loc_procs = (int* ) malloc( commSize*sizeof(int) );
	    fst_nz_procs = (int* ) malloc( commSize*sizeof(int) );
	    nz_loc_procs = (int* ) malloc ( commSize*sizeof(int) );

	    for( int ii= 0; ii<commSize; ii++){
	      fst_row_procs [ii] = ii*m_loc_fst;
	      m_loc_procs [ii] = m_loc_fst;
	    }
	    
	    if( m_loc_fst*commSize != m ) 
	      m_loc_procs [commSize-1] = m-m_loc_fst*( commSize-1 );
	    
	    for( int ii= 0; ii<commSize; ii++){
	      nz_loc_procs [ii] = AA.lg[fst_row_procs[ii]+m_loc_procs[ii] ] - AA.lg[fst_row_procs[ii]];
	      fst_nz_procs [ii] = AA.lg[fst_row_procs[ii]];
	    }

	   
	    /*
	      tab_a= (int* ) malloc( nz*sizeof(double) );
	      tab_cl = (int* ) malloc( nz*sizeof(int) );
	      tab_lg = (int* ) malloc ( n*sizeof(int) );
	    */
	    tab_a  = AA.a;
	    tab_cl = AA.cl;
	    tab_lg = AA.lg;
	  }

	  MPI_Scatter( nz_loc_procs, 1, MPI_INT, &nz_loc, 1, MPI_INT, 0, comm);
	  cout << "nz_loc("<<myid<<")="<< nz_loc << endl;
	  cout << "m_loc("<<myid<<")="<< m_loc << endl;
	  // allocation des tableaux locaux
	  irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	  jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	  a_loc    = (double*) malloc(sizeof(double)*nz_loc);
	  tab_lg_loc = (int*) malloc(sizeof(int)*(m_loc) );
	
	  MPI_Scatterv(  tab_a, nz_loc_procs, fst_nz_procs, MPI_DOUBLE, a_loc, nz_loc, MPI_DOUBLE, 0, comm);
	  MPI_Scatterv( tab_cl, nz_loc_procs, fst_nz_procs, MPI_INT, jcn_loc, nz_loc, MPI_INT, 0, comm);
	  MPI_Scatterv( tab_lg,  m_loc_procs, fst_row_procs, MPI_INT, tab_lg_loc, m_loc, MPI_INT, 0, comm);
	
	  int jj=0;
	  for(int  ii=0; ii<m_loc-1; ii++)
	    for(int ii1= tab_lg_loc[ii]; ii1 < tab_lg_loc[ii+1]; ii1++){
	      irn_loc[jj] = fst_row+ii+1;
	      jj++;
	    }
	  
	  for(int ii1= tab_lg_loc[m_loc-1]; ii1 < tab_lg_loc[0]+nz_loc; ii1++){
	    irn_loc[jj] = fst_row+(m_loc-1)+1;
	    jj++;
	  }
	  
	  for(int ii=0; ii < nz_loc; ii++){	    
	    jcn_loc[ii] = jcn_loc[ii]+1; 
	  }

	  assert( jj == nz_loc );
	  
	  free( tab_lg_loc );
	  
	  id.nz_loc = nz_loc;
	  id.irn_loc = irn_loc;
	  id.jcn_loc = jcn_loc;
	  id.a_loc = a_loc;
	 
	  if( myid == 0 ){
	    free( fst_row_procs );
	    free( m_loc_procs );
	    free( fst_nz_procs );
	    free( nz_loc_procs );
	  }
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


    
    if( verbosity > 1){
      /* information given by mumps*/
      int Rinfo=20;
      int Sinfo=40;
      // in Freefem++ we give only global information
      if(myid == 0){
	printf("Global Output Information of MUMPS: RINFOG and INFOG \n");
	printf("=============  After Factorisation ==================\n");
	for(int ii=0; ii< Rinfo; ii++) 
	  printf( "RINFOG[%d]= %f \n", ii, id.RINFOG(ii+1) );
	printf("=====================================================\n");
	for(int ii=0; ii< Sinfo; ii++) 
	  printf( "INFOG[%d]= %f \n", ii, id.INFOG(ii+1) );
	printf("=====================================================\n");
      }
    }
    
    if( verbosity )
      if(myid==0){
	finishtime = clock();
	timeused= (finishtime-starttime)/(1000 );
	printf("=====================================================\n");
	cout << "MUMPS : time factorisation  :: " << timeused << " ms" <<endl;
	printf("=====================================================\n");
      }


  }
  void Solver(const MatriceMorse<double> &AA,KN_<double> &x,const KN_<double> &b) const  {
    long int starttime,finishtime;
    long int timeused;
    /////////////////////////////
    double *rhs;
    int job;

    if(verbosity) starttime = clock();
    
    ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
  
    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if(jcn != NULL)
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]+1;
   

    if ( !(rhs = (double*) malloc(sizeof(double)*m) ) ){
      printf("Pb allocate rhs in MUMPS\n");
      exit(1);
    }
   
    for(int ii=0; ii<m; ++ii){
      rhs[ii] = b[ii];
    }

    if( myid == 0 )
      id.rhs=rhs;

    /* solve linear problem */

    id.job=3;
    dmumps_c(&id);

    if( myid==0 ){
      x=id.rhs;
      MPI_Bcast( x, n, MPI_DOUBLE, 0, comm );
    }
    else
      MPI_Bcast( x, n, MPI_DOUBLE, 0, comm );
    
    // deallocation de rhs
    free(rhs);

    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if(jcn != NULL)
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]-1;
    
    if(verbosity) cout << "   x min max " << x.min() << " " <<x.max() << endl;

    


    if( verbosity >1){
      /* information given by mumps*/
      int Rinfo=20;
      int Sinfo=40;
      // in Freefem++ we give only global information
      if(myid == 0){
	printf("Global Output Information of MUMPS: RINFOG and INFOG \n");
	printf("=============  After Solving       ==================\n");
	for(int ii=0; ii< Rinfo; ii++) 
	  printf( "RINFOG[%d]= %f \n", ii, id.RINFOG(ii+1) );
	printf("=====================================================\n");
	for(int ii=0; ii< Sinfo; ii++) 
	  printf( "INFOG[%d]= %f \n", ii, id.INFOG(ii+1) );
	printf("=====================================================\n");
      }
    }

    if( verbosity ){
      if(myid==0){
	finishtime = clock();
	timeused= (finishtime-starttime)/(1000 );
	printf("=====================================================\n");
	cout << "MUMPS : time solve step  :: " << timeused << " ms" <<endl;
	printf("=====================================================\n");
      }
    }
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
  mutable MPI_Comm comm;

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


  static const int JOB_INIT=-1;
  static const int JOB_END=-2;
  static const int USE_COMM_WORLD= -987654;

  // variable complex 
  mutable ZMUMPS_STRUC_C id;

  /* variable d'informations */ 
  
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define CNTL(I) cntl[(I)-1] /* macro s.t. indices match documentation */
#define RINFOG(I) rinfog[(I)-1] /* macro s.t. indices match documentation */
#define INFOG(I) infog[(I)-1] /* macro s.t. indices match documentation */
public:
  zSolveMUMPSmpi(const MatriceMorse<Complex> &AA,int strategy,double ttgv, double epsilon,
		 double pivot,double pivot_sym, string param_string, string datafile, KN<long> &param_int, 
		 KN<double> &param_double, KN<long> &pperm_r, KN_<long> &pperm_c, KN<double> &pscale_r,KN<double> &pscale_c, MPI_Comm  * mpicommw) : 
    eps(epsilon),epsr(0),
    tgv(ttgv), string_option(param_string), data_option(datafile), perm_r(pperm_r), perm_c(pperm_c), 
    tol_pivot_sym(pivot_sym),tol_pivot(pivot), scale_r(pscale_r), scale_c(pscale_c)
  { 
    long int starttime,finishtime;
    long int timeused;
    if(verbosity) starttime = clock();
    int dataint[40];
    double datadouble[15];
    int ierr;

    if(mpicommw==0){
      comm=MPI_COMM_WORLD;
    }
    else
      comm= *mpicommw;
    
    /* ------------------------------------------------------------
       INITIALIZE THE MUMPS PROCESS GRID. 
       ------------------------------------------------------------*/
    ierr = MPI_Comm_rank(comm, &myid);
    if( myid ==0){
      n    = AA.n;
      m    = AA.m; 
      nz   = AA.nbcoef;
    
      MPI_Bcast(  &n, 1, MPI_INT,  0, comm );
      MPI_Bcast(  &m, 1, MPI_INT,  0, comm );
      MPI_Bcast( &nz, 1, MPI_INT,  0, comm );
    }
    else{
      MPI_Bcast(  &n, 1, MPI_INT,  0, comm );
      MPI_Bcast(  &m, 1, MPI_INT,  0, comm );
      MPI_Bcast( &nz, 1, MPI_INT,  0, comm );
    }
   

    if( !(param_int==NULL) ) 
      assert( param_int.N() == 42);
    if( !(param_double==NULL) ) 
      assert( param_double.N() == 15);
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
   
    // initialisation par defaut
 
    SYM=0; PAR=1;
    
    //      if(!string_option.empty()) 
    //       {	
    // 	read_options_freefem(&string_option,&SYM,&PAR);
    //       }
    if( !(param_int==NULL) ){
      SYM = param_int[0];
      PAR = param_int[1];
    }
    else 
      if(!data_option.empty())
	{
	  if(myid==0){
	    char * retfile= new char[data_option.size()+1];
	    strcpy(retfile, (&data_option)->c_str());
	    printf("read data from file %s\n", retfile);
	    FILE *pFile=fopen(retfile,"rt");
	    
	    int     i_data=0;
	    int     d_data=0.;
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
	    assert(i_data == 40);
	    while( !feof(pFile) && d_data < 15){
	      fgets(data,256,pFile);
	      tictac = strtok(data," /!#\t\n");
	      datadouble[d_data] = (double)atof(tictac);
	      d_data++;
	    }  
	    assert(d_data == 15);
	    fclose(pFile);
	    delete [] retfile;
	    
	    MPI_Bcast(  &SYM, 1, MPI_INT,  0, comm );
	    MPI_Bcast(  &PAR, 1, MPI_INT,  0, comm );
	    
	    MPI_Bcast(  dataint, 40, MPI_INT,  0, comm );
	    MPI_Bcast(  datadouble, 15, MPI_DOUBLE,  0, comm );
	    
	    fclose(pFile);
	    delete [] retfile;
	  }
	  else{
	    
	    MPI_Bcast(  &SYM, 1, MPI_INT,  0, comm );
	    MPI_Bcast(  &PAR, 1, MPI_INT,  0, comm );
	    
	    MPI_Bcast(  dataint, 40, MPI_INT,  0, comm );
	    MPI_Bcast(  datadouble, 15, MPI_DOUBLE,  0, comm );
	  }
	}
    
    /* Initialize a MUMPS instance. Use comm */
    id.job=JOB_INIT; 
    id.par=PAR; 
    id.sym=SYM;
    id.comm_fortran=(MUMPS_INT) MPI_Comm_c2f(comm); 

    zmumps_c(&id);

     /* set parameter of mumps */
    if( !(param_int==NULL) || !(param_double==NULL) ){
      if(!data_option.empty()){ 
	printf("read option before with the file %s\n",&data_option);
	exit(1);
      }
      
      if( !(param_int==NULL) ){ 
	for(int ii=0; ii<40; ii++)
	  id.ICNTL(ii+1) = param_int[ii+2];
      }
      else
	// parameter by default
	id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;

      if( !(param_double==NULL) )
	for(int ii=0; ii<15; ii++)
	  id.CNTL(ii+1) = param_double[ii];
      
    }    
    else 
      if(!data_option.empty()){
	for(int ii=0; ii<40; ii++)
	  id.ICNTL(ii+1) = dataint[ii];
	for(int ii=0; ii<15; ii++)
	  id.CNTL(ii+1) = datadouble[ii];
      }
      else{
	// parameter by default
	id.ICNTL(1)=-1; id.ICNTL(2)=-1; id.ICNTL(3)=-1; id.ICNTL(4)=0;
      }

    // uniquement donner au host 
    if(myid==0){
      if( !(perm_r==NULL) && id.ICNTL(7)==1){
	for(int ii=0; ii<n; ii++) id.perm_in[ii] = pperm_r[ii];
      }
      // a decommenter
      //if( !(perm_c==NULL) && id.ICNTL(6)==1){
      //for(int ii=0; ii<m; ii++) id.perm_in[ii] = pperm_c[ii];
      //}
      if( !(scale_r==NULL) && !(scale_c==NULL) && id.ICNTL(8)==-1 ){
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

	if(myid == 0) { // uniquement sur le proc 0
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


	  //if (myid == 0) { // changement uniquement sur le proc 0
	  id.irn=irn; id.jcn=jcn;
	}
	
	/* Call the MUMPS package. */
	// Analyse   
	id.job=1;
	zmumps_c(&id);

	if(id.ICNTL(18) == 1 ){

	   if( PAR == 0 ){ 
	    int *nz_loc_procs;
	    int *fst_nz_procs;
	    int *irn_g;
	    int *jcn_g;
	    Complex *a_g;
	    int commSize;
	    
	    MPI_Comm_size(comm,&commSize);
	    cout << commSize << "commSize" << "nz =" << nz << endl;
	    if(myid==0){
	      // allocation des differents tableaux
	      nz_loc_procs = (int*) malloc ( commSize*sizeof(int) );
	      for(int ii=0; ii<commSize; ii++){
		nz_loc_procs[ii]=0;
	      }
	      for(int ii=0; ii<nz; ii++){
		nz_loc_procs[ id.mapping[ii] ]++;
	      }
	      assert(nz_loc_procs[0] == 0);
	      nz_loc_procs[0] = 2;

	      fst_nz_procs = (int*) malloc ( commSize*sizeof(int) );	      
	      fst_nz_procs[0] = 0;
	      for(int ii=1; ii<commSize; ii++){
		fst_nz_procs[ii] = fst_nz_procs[ii-1]+nz_loc_procs[ii-1];
	      }
	      
	      irn_g = (int*) malloc( sizeof(int)*(nz+2) );
	      jcn_g = (int*) malloc( sizeof(int)*(nz+2) );
	      a_g   = (Complex*) malloc( 2*sizeof(double)*(nz+2) );
	      
	      int *index_p;
	      index_p = (int*) malloc ( commSize*sizeof(int) );
	      for(int ii=0; ii<commSize; ii++)
		index_p[ii] =0;
	      
	      irn_g[ 0 ] = 1;
	      jcn_g[ 0 ] = 1;
	      a_g  [ 0 ] = 1.;
	      
	      irn_g[ 1 ] = 1;
	      jcn_g[ 1 ] = 1;
	      a_g  [ 1 ] = 1.;
	      
	      for(int ii=0;ii<nz; ii++){	      
		int jj1 = id.mapping[ii];
		int jj2 = fst_nz_procs[jj1] + index_p[jj1];
		assert(jj2 > 1);
		irn_g[ jj2 ] =  irn[ ii ];
		jcn_g[ jj2 ] =  jcn[ ii ];
		a_g  [ jj2 ] = AA.a[ ii ];
		cout << "jj2= " << jj2 << endl;
		assert( jj2 < nz+2);
		index_p[jj1]++;
	      }
	      free(index_p);
	      
	    }
	    
	    MPI_Scatter( nz_loc_procs, 1, MPI_INT, &nz_loc, 1, MPI_INT, 0, comm);
	    
	    // allocation des tableaux locaux
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc    = (Complex*) malloc(2*sizeof(double)*nz_loc);
	    
	    MPI_Scatterv(   a_g, nz_loc_procs, fst_nz_procs, MPI_DOUBLE_COMPLEX, a_loc, nz_loc, MPI_DOUBLE_COMPLEX, 0, comm);
	    MPI_Scatterv( jcn_g, nz_loc_procs, fst_nz_procs, MPI_INT, jcn_loc, nz_loc, MPI_INT, 0, comm);
	    MPI_Scatterv( irn_g, nz_loc_procs, fst_nz_procs, MPI_INT, irn_loc, nz_loc, MPI_INT, 0, comm);
	    	    
	    if( myid > 0){
	      id.nz_loc = nz_loc;
	      id.irn_loc = irn_loc;
	      id.jcn_loc = jcn_loc;
	      id.a_loc = mumps_dc(a_loc);
	    }
	    
	    if( myid == 0){
	      //free( irn_loc );
	      //free( jcn_loc );
	      //free( a_loc );
	      free( nz_loc_procs );
	      free( fst_nz_procs );
	      free( irn_g );
	      free( jcn_g );
	      free( a_g );
	    }
	   }
	
	  
	  if( PAR == 1 ){ 
	    int *nz_loc_procs;
	    int *fst_nz_procs;
	    int *irn_g;
	    int *jcn_g;
	    Complex *a_g;
	    int commSize;
	    
	    MPI_Comm_size(comm,&commSize);
	    
	    if(myid==0){
	      // allocation des differents tableaux
	      nz_loc_procs = (int*) malloc ( commSize*sizeof(int) );
	      for(int ii=0; ii<commSize; ii++){
		nz_loc_procs[ii]=0;
	      }
	      for(int ii=0; ii<nz; ii++){
		nz_loc_procs[ id.mapping[ii] ]++;
	      }
	      
	      fst_nz_procs = (int*) malloc ( commSize*sizeof(int) );
	      
	      fst_nz_procs[0] = 0;
	      for(int ii=1; ii<commSize; ii++){
		fst_nz_procs[ii] = fst_nz_procs[ii-1]+nz_loc_procs[ii-1];
	      }
	      
	      irn_g = (int*) malloc(sizeof(int)*nz);
	      jcn_g = (int*) malloc(sizeof(int)*nz);
	      a_g   = (Complex*) malloc(2*sizeof(double)*nz);
	      
	      int *index_p;
	      index_p = (int*) malloc ( commSize*sizeof(int) );
	      for(int ii=0; ii<commSize; ii++)
		index_p[ii] =0;
	      
	      for(int ii=0;ii<nz; ii++){	      
		int jj1 = id.mapping[ii];
		int jj2 = fst_nz_procs[jj1] + index_p[jj1];

		irn_g[ jj2 ] =  irn[ ii ];
		jcn_g[ jj2 ] =  jcn[ ii ];
		a_g  [ jj2 ] = AA.a[ ii ];
		index_p[jj1]++;
	      }
	      free(index_p);
	      
	    }
	    
	    MPI_Scatter( nz_loc_procs, 1, MPI_INT, &nz_loc, 1, MPI_INT, 0, comm);
	    
	    // allocation des tableaux locaux
	    irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	    a_loc   = (Complex*) malloc(2*sizeof(double)*nz_loc);
	    
	    MPI_Scatterv(   a_g, nz_loc_procs, fst_nz_procs, MPI_DOUBLE_COMPLEX, a_loc, nz_loc, MPI_DOUBLE_COMPLEX, 0, comm);
	    MPI_Scatterv( jcn_g, nz_loc_procs, fst_nz_procs, MPI_INT, jcn_loc, nz_loc, MPI_INT, 0, comm);
	    MPI_Scatterv( irn_g, nz_loc_procs, fst_nz_procs, MPI_INT, irn_loc, nz_loc, MPI_INT, 0, comm);
	    
	    id.nz_loc = nz_loc;
	    id.irn_loc = irn_loc;
	    id.jcn_loc = jcn_loc;
	    id.a_loc = mumps_dc(a_loc);
	    
	    
	    if( myid == 0){
	      free( nz_loc_procs );
	      free( fst_nz_procs );
	      free( irn_g );
	      free( jcn_g );
	      free( a_g );
	    }
	  }




	  //printf("id.ICNTL(18)==1 pas prevus a construire \n");
	  //exit(1);

	  /* // version matrice sur tous les processeurs
	  if(myid==0){

	    MPI_Bcast( id.mapping, nz, MPI_INT,  0, comm );

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
	    MPI_Bcast( mapping, nz, MPI_INT,  0, comm );
	    nz_loc=0;

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
	  */

	  /* Call the MUMPS package. */
	  // Factorisation   
	  id.job=2;
	  zmumps_c(&id);	  
	  
	}
      

	if(id.ICNTL(18) == 2 ){
	   printf("id.ICNTL(18)==2 not yet available \n");
	   exit(1);

	  if(PAR == 0){ 
	    printf("id.ICNTL(18)==2 pas prevus \n");
	    exit(1);
	    if(myid !=0) {
	      int commSize;	    
	      ierr=MPI_Comm_size(comm,&commSize);
	      commSize=commSize-1;
	      int myidpar=myid-1;
	      int m_loc_fst = m/commSize;
	      int m_loc;
	      if( myidpar == commSize-1 && ( m_loc_fst*commSize != m ) )  
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

	    int commSize;	    
	    ierr=MPI_Comm_size(comm,&commSize);
	    int m_loc_fst = m/commSize;
	    int m_loc;
	    if( myid == commSize-1 && ( m_loc_fst*commSize != m ) )  
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
// 	    if(myid !=0) {
// 	      int commSize;	    
// 	      ierr=MPI_Comm_size(comm,&commSize);
// 	      commSize=commSize-1;
// 	      int myidpar=myid-1;
// 	      int m_loc_fst;
// 	      m_loc_fst= m/commSize;
// 	      int m_loc;
// 	      if( myidpar == commSize-1 && ( m_loc_fst*commSize != m ) )  
// 		m_loc = m-m_loc_fst*( commSize-1 );
// 	      else
// 		m_loc = m_loc_fst;
//	      
// 	      int fst_row;
// 	      fst_row= myidpar*m_loc_fst;
// 	      nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
//	      
// 	      // allocation des tableaux
// 	      irn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	      jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	      a_loc    = (Complex*) malloc(sizeof(Complex)* nz_loc);
//	      
// 	      int fst_nnz;
// 	      fst_nnz = AA.lg[fst_row];
// 	      for(int ii=0; ii < nz_loc; ii++){
// 		a_loc[ii] = AA.a[fst_nnz+ii];
// 		jcn_loc[ii] = AA.cl[fst_nnz+ii]+1;
// 	      }	      
// 	      for(int ii=fst_row; ii< fst_row+m_loc; ii++){
// 		for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
// 		  irn_loc[ii1-fst_nnz] = ii+1;  
// 	      }	      
// 	      id.nz_loc = nz_loc;
// 	      id.irn_loc = irn_loc;
// 	      id.jcn_loc = jcn_loc;
// 	      id.a_loc = mumps_dc(a_loc);
// 	    }

	    // definition de variables
	  int commSize;
	  int m_loc_fst;	 
	  int m_loc;
	  int fst_row;
	  

	  int *nz_loc_procs;
	  int *fst_nz_procs;
	  int *m_loc_procs;
	  int *fst_row_procs;
	 

	  Complex *tab_a;
	  int *tab_cl;
	  int *tab_lg;
	  int *tab_lg_loc;

	  MPI_Comm_size(comm,&commSize);
	  
	  if( myid !=0){
	    int commSizemm;
	    int myidpar=myid-1;

	    commSizemm = commSize-1;
	    m_loc_fst= m/commSizemm;
	  
	    if( myidpar == commSizemm-1 && ( m_loc_fst*commSizemm != m ) )  
	      m_loc = m-m_loc_fst*( commSizemm-1 );
	    else
	      m_loc = m_loc_fst;
	  
	    if(verbosity > 5){
	      fst_row = myidpar*m_loc_fst;
	      cout << "   myid = " << myid << endl;
	      cout <<"   m_loc = " << m_loc << endl;
	      cout <<" fst_row = " << fst_row << endl;
	    }

	  }
	  if( myid ==0){

	    int commSizemm;
	    commSizemm = commSize-1;
	    m_loc_fst= m/commSizemm;

	    fst_row_procs = (int* ) malloc( commSize*sizeof(int) );
	    m_loc_procs = (int* ) malloc( commSize*sizeof(int) );
	    fst_nz_procs = (int* ) malloc( commSize*sizeof(int) );
	    nz_loc_procs = (int* ) malloc ( commSize*sizeof(int) );
	    
	    
	    fst_row_procs [0] = 0;
	    m_loc_procs   [0] = 0;

	    for( int ii= 1; ii<commSize; ii++){
	      fst_row_procs [ii] = (ii-1)*m_loc_fst;
	      m_loc_procs [ii] = m_loc_fst;
	    }
	    
	    if( m_loc_fst*(commSize-1) != m ) 
	      m_loc_procs [commSize-1] = m-m_loc_fst*( (commSize-1)-1 );


	    nz_loc_procs [0] = 0;
	    fst_nz_procs [0] = 0;
	    
	    for( int ii= 1; ii<commSize; ii++){
	      nz_loc_procs [ii] = AA.lg[fst_row_procs[ii]+m_loc_procs[ii] ] - AA.lg[fst_row_procs[ii]];
	      fst_nz_procs [ii] = AA.lg[fst_row_procs[ii]];
	    }

	   
	    /*
	      tab_a= (int* ) malloc( nz*sizeof(double) );
	      tab_cl = (int* ) malloc( nz*sizeof(int) );
	      tab_lg = (int* ) malloc ( n*sizeof(int) );
	    */
	    tab_a  = AA.a;
	    tab_cl = AA.cl;
	    tab_lg = AA.lg;
	  }

	  MPI_Scatter( nz_loc_procs, 1, MPI_INT, &nz_loc, 1, MPI_INT, 0, comm);
	  MPI_Scatter( m_loc_procs,  1, MPI_INT, &m_loc, 1, MPI_INT, 0, comm);
	  MPI_Scatter( fst_row_procs,  1, MPI_INT, &fst_row, 1, MPI_INT, 0, comm);

	  if(verbosity > 5){
	    cout << "after scatter " << myid << endl;
	    cout << "   myid = " << myid << endl;
	    cout <<"   m_loc = " << m_loc << endl;
	    cout <<" fst_row = " << fst_row << endl;
	  }
	  // allocation des tableaux locaux
	  irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	  jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	  a_loc    = (Complex*) malloc(2*sizeof(double)*nz_loc);
	  tab_lg_loc = (int*) malloc(sizeof(int)*(m_loc) );
	  
	

	  MPI_Scatterv(  tab_a, nz_loc_procs, fst_nz_procs, MPI_DOUBLE_COMPLEX, a_loc, nz_loc, MPI_DOUBLE_COMPLEX, 0, comm);
	  MPI_Scatterv( tab_cl, nz_loc_procs, fst_nz_procs, MPI_INT, jcn_loc, nz_loc, MPI_INT, 0, comm);
	  MPI_Scatterv( tab_lg,  m_loc_procs, fst_row_procs, MPI_INT, tab_lg_loc, m_loc, MPI_INT, 0, comm);
	  
	
	  int jj=0;
	  for(int  ii=0; ii<m_loc-1; ii++)
	    for(int ii1= tab_lg_loc[ii]; ii1 < tab_lg_loc[ii+1]; ii1++){
	      irn_loc[jj] = fst_row+ii+1;
	      jj++;
	    }
	  
	  for(int ii1= tab_lg_loc[m_loc-1]; ii1 < tab_lg_loc[0]+nz_loc; ii1++){
	    irn_loc[jj] = fst_row+(m_loc-1)+1;
	    jj++;
	  }
	  
	  for(int ii=0; ii < nz_loc; ii++){	    
	    jcn_loc[ii] = jcn_loc[ii]+1; 
	  }

	  assert( jj == nz_loc );
	  
	  free( tab_lg_loc );
	  
	  id.nz_loc = nz_loc;
	  id.irn_loc = irn_loc;
	  id.jcn_loc = jcn_loc;
	  id.a_loc = mumps_dc(a_loc);
	 
	  if( myid == 0 ){
	    free( fst_row_procs );
	    free( m_loc_procs );
	    free( fst_nz_procs );
	    free( nz_loc_procs );
	  }
	   

	 }
	 if(PAR ==1) {
	  
	   
// 	   int commSize;
// 	   ierr=MPI_Comm_size(comm,&commSize);
// 	   int m_loc_fst;
// 	   m_loc_fst= m/commSize;
// 	   int m_loc;
// 	   if( myid == commSize-1 && ( m_loc_fst*commSize != m ) )  
// 	     m_loc = m-m_loc_fst*( commSize-1 );
// 	   else
// 	     m_loc = m_loc_fst;
//	   
// 	   int fst_row;
// 	   fst_row = myid*m_loc_fst;
// 	   nz_loc = AA.lg[fst_row+m_loc]-AA.lg[fst_row];
//	    
// 	   allocation des tableaux
// 	   irn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	   jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
// 	   a_loc    = (Complex*) malloc(sizeof(Complex)*nz_loc);
//	   
// 	   int fst_nnz;
// 	   fst_nnz= AA.lg[fst_row];
// 	   for(int ii=0; ii < nz_loc; ii++){
// 	     a_loc[ii] = AA.a[fst_nnz+ii];
// 	     jcn_loc[ii] = AA.cl[fst_nnz+ii]+1; 
// 	   }
//   
// 	   for(int ii=fst_row; ii< fst_row+m_loc; ii++){
// 	     for(int ii1=AA.lg[ii]; ii1 < AA.lg[ii+1]; ii1++ )
// 	       irn_loc[ii1-fst_nnz] = ii+1;  
// 	   }
//	   
// 	   id.nz_loc = nz_loc;
// 	   id.irn_loc = irn_loc;
// 	   id.jcn_loc = jcn_loc;
// 	   id.a_loc = mumps_dc(a_loc);

	   // definition de variables
	   int commSize;
	   int m_loc_fst;	 
	   int m_loc;
	   int fst_row;
	   

	   int *nz_loc_procs;
	   int *fst_nz_procs;
	   int *m_loc_procs;
	   int *fst_row_procs;
	   
	   
	   Complex *tab_a;
	   int *tab_cl;
	   int *tab_lg;
	   int *tab_lg_loc;
	   
	   MPI_Comm_size(comm,&commSize);
	   m_loc_fst= m/commSize;
	   
	   if( myid == commSize-1 && ( m_loc_fst*commSize != m ) )  
	     m_loc = m-m_loc_fst*( commSize-1 );
	   else
	     m_loc = m_loc_fst;	  
	   
	   fst_row = myid*m_loc_fst;
	   
	   if( myid ==0){
	     fst_row_procs = (int* ) malloc( commSize*sizeof(int) );
	     m_loc_procs = (int* ) malloc( commSize*sizeof(int) );
	     fst_nz_procs = (int* ) malloc( commSize*sizeof(int) );
	     nz_loc_procs = (int* ) malloc ( commSize*sizeof(int) );

	     for( int ii= 0; ii<commSize; ii++){
	       fst_row_procs [ii] = ii*m_loc_fst;
	       m_loc_procs [ii] = m_loc_fst;
	     }
	     
	     if( m_loc_fst*commSize != m ) 
	       m_loc_procs [commSize-1] = m-m_loc_fst*( commSize-1 );
	     
	     for( int ii= 0; ii<commSize; ii++){
	       nz_loc_procs [ii] = AA.lg[fst_row_procs[ii]+m_loc_procs[ii] ] - AA.lg[fst_row_procs[ii]];
	       fst_nz_procs [ii] = AA.lg[fst_row_procs[ii]];
	     }
	     
	     
	     /*
	       tab_a= (int* ) malloc( nz*sizeof(double) );
	       tab_cl = (int* ) malloc( nz*sizeof(int) );
	       tab_lg = (int* ) malloc ( n*sizeof(int) );
	     */
	     tab_a  = AA.a;
	     tab_cl = AA.cl;
	     tab_lg = AA.lg;
	   }
	   
	   MPI_Scatter( nz_loc_procs, 1, MPI_INT, &nz_loc, 1, MPI_INT, 0, comm);
	  
	   // allocation des tableaux locaux
	   irn_loc = (int*) malloc(sizeof(int)*nz_loc);
	   jcn_loc = (int*) malloc(sizeof(int)*nz_loc);
	   a_loc    = (Complex*) malloc(2*sizeof(double)*nz_loc);
	   tab_lg_loc = (int*) malloc(sizeof(int)*(m_loc) );
	   
	   MPI_Scatterv(  tab_a, nz_loc_procs, fst_nz_procs, MPI_DOUBLE_COMPLEX, a_loc, nz_loc, MPI_DOUBLE_COMPLEX, 0, comm);
	   MPI_Scatterv( tab_cl, nz_loc_procs, fst_nz_procs, MPI_INT, jcn_loc, nz_loc, MPI_INT, 0, comm);
	   MPI_Scatterv( tab_lg,  m_loc_procs, fst_row_procs, MPI_INT, tab_lg_loc, m_loc, MPI_INT, 0, comm);
	
	   int jj=0;
	   for(int  ii=0; ii<m_loc-1; ii++)
	     for(int ii1= tab_lg_loc[ii]; ii1 < tab_lg_loc[ii+1]; ii1++){
	       irn_loc[jj] = fst_row+ii+1;
	      jj++;
	     }
	   
	   for(int ii1= tab_lg_loc[m_loc-1]; ii1 < tab_lg_loc[0]+nz_loc; ii1++){
	     irn_loc[jj] = fst_row+(m_loc-1)+1;
	     jj++;
	   }
	   
	   for(int ii=0; ii < nz_loc; ii++){	    
	     jcn_loc[ii] = jcn_loc[ii]+1; 
	   }
	   
	   assert( jj == nz_loc );
	   
	   free( tab_lg_loc );
	   
	   id.nz_loc  = nz_loc;
	   id.irn_loc = irn_loc;
	   id.jcn_loc = jcn_loc;
	   id.a_loc = mumps_dc(a_loc);
	 
	   if( myid == 0 ){
	     free( fst_row_procs );
	     free( m_loc_procs );
	     free( fst_nz_procs );
	     free( nz_loc_procs );
	   }
	   

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
 

    if( verbosity >1){
      /* information given by mumps*/
      int Rinfo=20;
      int Sinfo=40;
      // in Freefem++ we give only global information
      if(myid == 0){
	printf("Global Output Information of MUMPS: RINFOG and INFOG \n");
	printf("=============  After Factorisation ==================\n");
	for(int ii=0; ii< Rinfo; ii++) 
	  printf( "RINFOG[%d]= %f \n", ii, id.RINFOG(ii+1) );
	printf("=====================================================\n");
	for(int ii=0; ii< Sinfo; ii++) 
	  printf( "INFOG[%d]= %f \n", ii, id.INFOG(ii+1) );
	printf("=====================================================\n");
      }
    }
    if( verbosity){
      if(myid==0){
	finishtime = clock();
	timeused= (finishtime-starttime)/(1000 );
	printf("=====================================================\n");
	cout << "MUMPS : time factorisation :: " << timeused << " ms" <<endl;
	printf("=====================================================\n");
      }

    }
        
  
  }

  void Solver(const MatriceMorse<Complex> &AA,KN_<Complex> &x,const KN_<Complex> &b) const  {
    long int starttime,finishtime;
    long int timeused;
    //*******************************************************************//
    //    depend pas de la forme de la matrice: distribuer ou assembler
    Complex *rhs;
    int job;
    
    if(verbosity) starttime = clock();

    ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
  
    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if(jcn != NULL)
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]+1;
    
    if ( !(rhs = (Complex*) malloc(sizeof(Complex)*m) ) ){
      printf("Pb allocate rhs in MUMPS\n");
      exit(1);
    }
   
    for(int ii=0; ii<m; ++ii){
      rhs[ii] = b[ii];
    }

    if( myid == 0 )
      id.rhs=mumps_dc(rhs);

    /* solve linear problem */
    id.job=3;
    zmumps_c(&id);

   
    if( myid==0 ){
      x=inv_mumps_dc(id.rhs); 
      MPI_Bcast( x,  n, MPI_DOUBLE_COMPLEX,  0, comm );
    }
    else
      MPI_Bcast( x,  n, MPI_DOUBLE_COMPLEX,  0, comm );
    
  
    
    // deallocation de rhs
    free(rhs);

    // indices des colonnes commence par 1 avec mumps 
    //  et 0 dans freefem ==> renumerotation
    if(jcn != NULL)
      for(int ii=0; ii<nz; ii++)
	jcn[ii] = jcn[ii]-1;


    if( verbosity > 1){
      /* information given by mumps*/
      int Rinfo=20;
      int Sinfo=40;
      // in Freefem++ we give only global information
      if(myid == 0){
	printf("Global Output Information of MUMPS: RINFOG and INFOG \n");
	printf("=============  After Solving       ==================\n");
	for(int ii=0; ii< Rinfo; ii++) 
	  printf( "RINFOG[%d]= %f \n", ii, id.RINFOG(ii+1) );
	printf("=====================================================\n");
	for(int ii=0; ii< Sinfo; ii++) 
	  printf( "INFOG[%d]= %f \n", ii, id.INFOG(ii+1) );
	printf("=====================================================\n");
      }
    }

    if(verbosity)
      if(myid==0){
	finishtime = clock();
	timeused= (finishtime-starttime)/(1000 );
	printf("=====================================================\n");
	cout << " MUMPS : time solve  :: " << timeused << " ms" <<endl;
	printf("=====================================================\n");
      }
    
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
      cout << " BuildSolverMUMPS<double>" << endl;
    return new dSolveMUMPSmpi(*A,ds.strategy, ds.tgv, ds.epsilon, ds.tol_pivot, ds.tol_pivot_sym, ds.sparams, ds.data_filename,
			      ds.lparams, ds.dparams, ds.perm_r, ds.perm_c, ds.scale_r, ds.scale_c,(MPI_Comm *)ds.commworld);
}




MatriceMorse<Complex>::VirtualSolver *
BuildSolverMUMPSmpi(DCL_ARG_SPARSE_SOLVER(Complex,A))
{
    if(verbosity>9)
      cout << " BuildSolverMUMPS<Complex>" << endl;
    return new zSolveMUMPSmpi(*A,ds.strategy, ds.tgv, ds.epsilon, ds.tol_pivot, ds.tol_pivot_sym, ds.sparams, ds.data_filename,  
			      ds.lparams, ds.dparams, ds.perm_r, ds.perm_c, ds.scale_r, ds.scale_c,(MPI_Comm *)ds.commworld);
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
    cout << "\n Add: MUMPS ,  defaultsolver defaultsolverMUMPS " << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  DefSparseSolver<double>::solver =BuildSolverMUMPSmpi;
  DefSparseSolver<Complex>::solver =BuildSolverMUMPSmpi;
  if(! Global.Find("defaultsolver").NotNull() )
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  Global.Add("defaulttoMUMPS","(",new OneOperator0<bool>(SetMUMPSmpi));
}

