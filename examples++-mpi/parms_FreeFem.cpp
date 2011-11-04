// ORIG-DATE: 04/2009
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : INRIA FUTUR
// AUTHOR   : Guy Atenekeng
// E-MAIL   : Guy_Antoine_Atenekeng_Kahou@lri.fr
//
//ff-c++-LIBRARY-dep: metis parms  blas mpifc 
//ff-c++-cpp-dep: 

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

#include <mpi.h>
#include  <iostream>
using namespace std;


#define MCW MPI_COMM_WORLD


#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"

//#include "lex.hpp"
#include "MatriceCreuse_tpl.hpp"


/* Explain here what foobar does */
#ifdef __cplusplus
extern "C" {
#endif
#include "psparslib.h"
#include "generaldefs.h"

#ifdef __cplusplus
}
#endif

extern "C" {
  #include "metis.h"
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define BUFLEN 100
#define NORHS 0

#define BUFLEN 100

//This functions come from pARMS package and consist to scale matrix 







class parm_param {
public:
	mutable PrePar prepar;
	mutable IterPar ipar;
	mutable int sol_type;
	mutable int  iov,scale, unsymm,method;
	mutable int solver;
	mutable int VERBOSE;
	

	public: parm_param(){
		PARMS_malloc(prepar,1,_PrePar);
		PARMS_malloc(ipar,1,_IterPar);  
		/*----------------------------------------------------------------------  *
		 * This function  sets some default  parameters just to get  started with *
		 * pARMS.  There are  two arrays  which define  the parameters  in pARMS. *
		 * Because ther\ e are so many methods, the number of parameters may seem *
		 * overwhelming.. However\ , not all  are used. For example when add_ilut *
		 * is  used  only  two parameters  are  r\  equired  which are  the  drop *
		 * tolerance and the  max fill per row.. This  function s\ ets everything *
		 * so that all the methods can be safely tested. YS -                     *
		 *                                                                        *
		 *------------------------------------------------------------------------*/
		  int i; 
		/*------------------------------------------------------------------------*/
		  for (i=1; i<18; i++)
		    ipar->ipar[i] = 0;
		/* parameters to be set for the solver -- not all of these are used --  
		     depending on the method selected -                                   */
		/* parameters associated with various accelerators (inner/outer)          */
		  ipar->ipar[0] = 3;       /* nlev in arms when used                      */
		  ipar->ipar[1] = 20;      /* block size in arms when used                */
		  ipar->ipar[5] = 30;      /* Outer iteration Krylov subs. dimension      */
		  ipar->ipar[6] = 200;     /* Max. outer iterations                       */
		  ipar->ipar[3] = 3;       /* Inner iteration Krylov subs. dimension      */
		  ipar->ipar[4] = 3;       /* Max. inner iterations  */
		  
		  ipar->pgfpar[0] = 0.01;  /* tolerance for inner iteration when used     */
		  ipar->pgfpar[1] = 1.e-10;/* tolerance for outer iteration               */
		/* preconditioning parameters                                             */
		  prepar->tolind = 0.1;    /* tolerance used for diag dom. filtration alg.*/
		  prepar->lfil[0] = 20;    /* lfil0 (ilut, iluk, and arms)                */
		  prepar->lfil[4] = 20;    /* lfil for Schur complement const.            */ 
		  prepar->lfil[5] = 20;    /* lfil for Schur complement const.            */ 
		  prepar->droptol[0]=.001; /* droptol0 (ilut, iluk, and arms)             */
		  prepar->droptol[4]=.001; /* droptol for Schur complement const.         */ 
		  prepar->droptol[5]=.001; /* droptol for Schur complement const.         */ 
		  prepar->mc = 1;          /* multicoloring or not in ILU when used       */
		 //Integer parameter
		 iov=0; scale=1; unsymm=1; method= 2; VERBOSE=0; /*no statistic infos has print*/
		 /*VERBOSE=1 Only informations on convergence will be print*/
		 /*VERBOSE=2 Only informations on time will be print*/
		  solver=0;
		 
		  memcpy(prepar->ipar, ipar->ipar, sizeof(int)*18);

		  for (i=1; i<4; i++) {
		    prepar->droptol[i] = prepar->droptol[0]; 
		    prepar->lfil[i] = prepar->lfil[0];
		  }
		  prepar->lfil[6] = prepar->lfil[5];
		  ipar->in_iters = 0;
		 
		/*-------------------- done  */                       
	}
	public: parm_param(const KN<long> &param_int, const KN<double> &param_double)
	{
		 int i; 
		PARMS_malloc(prepar,1,_PrePar);
		PARMS_malloc(ipar,1,_IterPar);  
		for (i=1; i<18; i++)
		 ipar->ipar[i] = 0;
		
		if(param_int.N()>0){if((param_int[0]<0)||(param_int[0]>2)) 
		{printf("%s","WRONG SOLVER INDEX WE SET DEFAULT ONE \n"); solver=0; } else solver=param_int[0];}else solver=0;

		if(param_int.N()>1){if((param_int[1]<0)||(param_int[1]>13)) 
		{printf("%s","WRONG INDEX FOR PRECONDITIONER, WE SET DEFAULT ONE \n"); method=2; } 
		else if((param_int[1]>=0)||(param_int[1]<=13)) method=param_int[1];}else method=2;
		
		if(param_int.N()>2){if(param_int[2]>0) ipar->ipar[5]=param_int[2]; else ipar->ipar[5]=30;}
		else ipar->ipar[5]=30; /* Outer iteration Krylov subs. dimension      */


		if(param_int.N()>3){if(param_int[3]>0) ipar->ipar[6]=param_int[3]; else ipar->ipar[6]=1000;}
		else ipar->ipar[6]=1000; /* Max. outer iterations    */


		if((method==3)||(method==7)||(method==11)) {if(param_int.N()>4) {if(param_int[4]>0) ipar->ipar[0]=param_int[4]; 
		else ipar->ipar[0]=3;}	else ipar->ipar[0]=3;  } /* nlev in arms when used                      */

		if(param_int.N()>5){ if(param_int[5]>0) ipar->ipar[3]=param_int[5]; else ipar->ipar[3]=3;}
		else ipar->ipar[3]=3; /* Inner iteration Krylov subs. dimension */

		if(param_int.N()>6){if(param_int[6]>0) ipar->ipar[4]=param_int[6]; else ipar->ipar[4]=3;}
		else ipar->ipar[4]=3;  /* Max. inner iterations  */

		if(param_int.N()>7){if(param_int[7]>=0) unsymm=param_int[7]; else unsymm=0;} else unsymm=0;

		if(param_int.N()>8) {if(param_int[8]>=0) iov=param_int[8]; else iov=0;} else iov=0;
		if(param_int.N()>9) {if(param_int[9]>0) scale=param_int[9]; else scale=1;}else scale=1;

		if(param_int.N()>10){if(param_int[10]>0) ipar->ipar[1]=param_int[10]; else ipar->ipar[1]=20; }else ipar->ipar[1]=20;
		 
			

		if(param_int.N()>11){if(param_int[11]>0) prepar->lfil[0]=param_int[11]; else prepar->lfil[0]=20;}
		else prepar->lfil[0]=20; /* lfil0(ilut, iluk, and arms)   */

		if(param_int.N()>12){if(param_int[12]>0) prepar->lfil[4]=param_int[12]; else prepar->lfil[4]=20;}
		else prepar->lfil[4]=20; /* lfil for Schur complement const.   */

		if(param_int.N()>13){if(param_int[13]>0) prepar->lfil[5]=param_int[13]; else prepar->lfil[13]=20;}
		else prepar->lfil[5]=20; /* lfil for Schur complement const.   */		

		if(param_int.N()>14){if(param_int[14]>0) prepar->mc=param_int[14]; else prepar->mc=1; }
		else prepar->mc=1; 

		if(param_int.N()>15){if(param_int[15]>0) ipar->in_iters=param_int[15]; else ipar->in_iters=0; }
		else ipar->in_iters=0;
		
		if(param_int.N()>16){if(param_int[16]>0) VERBOSE=param_int[16]; else VERBOSE=0; }else VERBOSE=0;


			


		if(param_double.N()>0){if(param_double[0]>0) ipar->pgfpar[1]=param_double[0]; else ipar->pgfpar[1]=1.e-08;}
		ipar->pgfpar[1]=1.e-08; /* tolerance for outer iteration               */

		if(param_double.N()>1){if(param_double[1]>0) ipar->pgfpar[0]=param_double[1]; else ipar->pgfpar[0]=0.01;}
		else ipar->pgfpar[0]=0.01;  /* tolerance for inner iteration when used     */

		if(param_double.N()>2) {if(param_double[2]>0) prepar->tolind = param_double[2]; else prepar->tolind=0.1;}
		else prepar->tolind=0.1; /* tolerance used for diag dom. filtration alg.*/

		if(param_double.N()>3){if(param_double[3]>0) prepar->droptol[0]=param_double[3]; else prepar->droptol[0]=.001;}
		else prepar->droptol[0]=.001;	/* droptol0 (ilut, iluk, and arms) */
		if(param_double.N()>4) {if(param_double[4]>0) prepar->droptol[4]=param_double[4]; else prepar->droptol[4]=.001;}
		else prepar->droptol[4]=.001; /* droptol for Schur complement const.         */ 

		if(param_double.N()>5){if(param_double[5]>0) prepar->droptol[5]=param_double[5]; else prepar->droptol[5]=.001;}
		else prepar->droptol[5]=.001; /* droptol for Schur complement const.         */ 

		

		memcpy(prepar->ipar, ipar->ipar, sizeof(int)*18);
		  for (i=1; i<4; i++) {
		    prepar->droptol[i] = prepar->droptol[0]; 
		    prepar->lfil[i] = prepar->lfil[0];
		  }
		  prepar->lfil[6] = prepar->lfil[5];
		 
		 
	}

  public: parm_param(string datafile,DistMatrix dm)
	{
	  FILE *fp;
	  char buf[BUFLEN], meth_buf[BUFLEN];
	  int num;
	  PARMS_malloc(prepar,1,_PrePar);
  	  PARMS_malloc(ipar,1,_IterPar);  

          parm_param();

	  char * filename=new char[datafile.length()+1]; 
                strcpy(filename,datafile.c_str()); 
 	 /* read parameters for preconditioner and iteration from file  'filename' */
  /*  ---- start modification by MS
   */
 	 if( (fp = fopen(filename, "r")) == NULL ){
          fprintf(stderr, "Cannot open file inputs\n");
    	  PARMS_Final();exit(1);
    	}
	  for (num=0; num< 18; num++)
	    ipar->ipar[num] = 0; 

  	num = 1;

  while(fgets(buf, BUFLEN, fp) != NULL) {
    switch(num) {
     case 1:                             /* solver */
      sscanf(buf, "%d", &solver);
	if((solver!=0)&&(solver!=1)&&(solver!=2))
	{printf("%s","WRONG SOLVER INDEX, WE SET DEFAULT ONE \n"); solver=0;}
      break;
    case 2:                             /* preconditionner */
      sscanf(buf,"%s",meth_buf); 
      method = assignprecon(meth_buf, dm);
	if((method<0)||(method>13)) 
	{printf("%s","WRONG INDEX FOR PRECONDITIONER, WE SET DEFAULT ONE \n"); method=2; }		
      break;
    case 3:                             /* im (Outer Krylov subs dim) */
      sscanf(buf,"%d",&ipar->ipar[5]);
      break;	
    case 4:                             /* im (Outer Krylov subs dim) */
      sscanf(buf,"%d",&ipar->ipar[6]);
      break;	
    case 5:                             /* outer tol */
      sscanf(buf,"%lf",&ipar->pgfpar[1]);
      break;	
    case 6:                             /* inner tol */
      sscanf(buf,"%lf",&ipar->pgfpar[0]);
      break;
    case 7:                             /* unsym */
      sscanf(buf, "%d", &unsymm);
      break;
    case 8:                             /* inim (inned Krylov subs dim) */
      sscanf(buf, "%d", &ipar->ipar[3]);
      break;
    case 9:                             /* Max. inner iterations     */ 
      sscanf(buf,"%d",&ipar->ipar[4]);
      break;
    case 10:                             /* nlev   */
      sscanf(buf, "%d", &iov);
      break;
    case 11:                             /*  scale  */ 
      sscanf(buf,"%d",&scale);
      break;
    case 12:                             /* For printing result */ 
      sscanf(buf, "%d", &VERBOSE);
      break;
    case 13:                             /* lfil0 (ilut, iluk, and arms) */
      sscanf(buf, "%d", &prepar->lfil[0]);
      break;
    case 14:                             /* lfil 4 (schur construction) */
      sscanf(buf, "%d", &prepar->lfil[4]);
      break;
    case 15:                             /* lfil 5 (ILUT for Schur) */
      sscanf(buf, "%d", &prepar->lfil[5]);
      break;
    case 16:                             /* bsize   */
      sscanf(buf,"%d",&ipar->ipar[1]);
      break;
    case 17:                              /* tolerance used for diag dom. filtration alg.*/  
      sscanf(buf,"%lf",&prepar->tolind);
      break;
    case 18:                             /*droptol (0) -- simliar to lfil0*/
      sscanf(buf, "%lf", &prepar->droptol[0]);
      break;
    case 19:                             /*droptol (4) -- simliar to lfil4*/
      sscanf(buf, "%lf", &prepar->droptol[4]);
      break;
    case 20:                             /*droptol (5) -- simliar to lfil5*/
      sscanf(buf, "%lf", &prepar->droptol[5]);
      break;
    case 21:                             /* multicoloring or not */ 
      sscanf(buf, "%d", &prepar->mc);
      break;
    
    default:
      break;
    }
    num++;
  }
  fclose(fp);
  
  memcpy(prepar->ipar, ipar->ipar, sizeof(int)*18);
  
  prepar->droptol[1] = prepar->droptol[2] = prepar->droptol[3] =
    prepar->droptol[0]; 
  prepar->droptol[6] = prepar->droptol[5];
  
  prepar->lfil[1] = prepar->lfil[2] = prepar->lfil[3] =
    prepar->lfil[0];
  prepar->lfil[6] = prepar->lfil[5];
  ipar->in_iters = 0; free(filename);

	}
	~parm_param(){
		free(prepar); free(ipar); 
	
	}
};


int assignprecon( char *precon_str, DistMatrix dm) 
{
/*------------------------------------------------------------------ 
     create preconditioner handler  
   * precon    --  preconditioning handler
   * add_ilu0  --  additive schwarz preconditioner with ilu0 as local
                   preconditioner
   * add_ilut  --  additive schwarz preconditioner with ilut as local
                   preconditioner
   * add_iluk  --  additive schwarz preconditioner with iluk as local
                   preconditioner
   * add_arms  --  additive schwarz preconditioner with arms as local
                   preconditioner
   * lsch_ilu0 --  schur complement preconditioner with ilu0 as local
                   preconditioner		  
   * lsch_ilut --  schur complement preconditioner with ilut as local
                   preconditioner		  
   * lsch_iluk --  schur complement preconditioner with iluk as local
                   preconditioner		  
   * lsch_arms --  schur complement preconditioner with arms as local
                   preconditioner
   * sch_gilu0 --  parallel ilu0 preconditioner
   * sch_sgs   --  Schur-Symmetric GS preconditioner
   *------------------------------------------------------------------*/
  /* set the method label for preconditioning and parameter for inner
     iteration -- actual labels (int) are defined in SRC/ARMS/data.h */ 

  int ierr, method;

  ierr = 0; method=0;

  if(!strcmp(precon_str, "add_ilu0")) method = add_ilu0;
  else if(!strcmp(precon_str, "add_ilut")) method = add_ilut;
  else if(!strcmp(precon_str, "add_iluk")) method = add_iluk;
  else if(!strcmp(precon_str, "add_arms")) method = add_arms;
  else if(!strcmp(precon_str, "lsch_ilu0")) method = lsch_ilu0;
  else if(!strcmp(precon_str, "lsch_ilut")) method = lsch_ilut;
  else if(!strcmp(precon_str, "lsch_iluk")) method = lsch_iluk;
  else if(!strcmp(precon_str, "lsch_arms")) method = lsch_arms;
  else if(!strcmp(precon_str, "rsch_ilu0")) method = rsch_ilu0;
  else if(!strcmp(precon_str, "rsch_ilut")) method = rsch_ilut;
  else if(!strcmp(precon_str, "rsch_iluk")) method = rsch_iluk;
  else if(!strcmp(precon_str, "rsch_arms")) method = rsch_arms;
  else if(!strcmp(precon_str, "sch_sgs")) method = sch_sgs;
  else if(!strcmp(precon_str, "sch_gilu0")) method = sch_gilu0;
  else ierr = 1; 
  char pcrM[40];
  strcpy(pcrM,"invalid choice for the preconditioner \n");
 
  if (ErrHand(ierr, dm, pcrM)) 
      exit(1) ; 
  return method;cout << "Cette resolution semble prendre du temps" << endl;

  }



#define minint(a, b)       ( (a) < (b) ? (a) : (b) )
#define maxint(a, b)       ( (a) > (b) ? (a) : (b) )
#define max(a, b)       ( (a) > (b) ? (a) : (b) )
#define min(a, b)   ( a < b ? (a) : (b) )
#define fabsmin(a, b)   ( fabs(a) < fabs(b) ? (a) : (b) )
#define fabsmax(a, b)   ( fabs(a) > fabs(b) ? (a) : (b) )



class dSolvePARMS :   public MatriceMorse<double>::VirtualSolver   {

	char 					mat_domain[BUFLEN];
    	int *riord;
	mutable int                     n, nnz, nloc;

	int 					iov, scale, unsymm, method,solver,VERBOSE;
	string data_option;
	int                     rk, size;
	mutable double 					*rhs1, res;
	//double 					t1, t2, t3, t4;
	double 					*u, *v;
	int 					job,i1,i2;
	mutable DistMatrix 				dm;  		//distributed matrix object
	mutable PrePar 					prepar;  
	mutable PreCon 					precon;  	//structure for preconditioner
	mutable IterPar 				ipar;		//structure for iteration
	mutable  int * maptmp1, *mapptr1;
	double eps, tol_pivot_sym,tgv, tol_pivot,epsr;
	mutable double t1,t2,t3,t4;
	mutable char * meth[14];
	mutable int *maptmp, *mapptr,*iwork1,*iwork;
	mutable double *scaletmpr, *scaletmpc;
	mutable MPI_Comm comm;
	char *p;
//Constructor construct the distribute matrix and also preconditionner
public:
  dSolvePARMS(const MatriceMorse<double> &AA,string datafile,  KN<long> &param_int,  KN<double> &param_double,  MPI_Comm  * mpicommw)
	{
	int n,i,job,nnz,ierr,tmp;
	int *ptr2,*id_rows2;
	double *vals2;
	double *AAv;
	int *p, *pr;	
	

	int 	 j,    node;
	/* Initialize PARMSARSLIB environment */
	 if(mpicommw==NULL){
	comm=MPI_COMM_WORLD;
	}
	else
	comm= *mpicommw;
	// comm=MPI_COMM_WORLD;
	 ierr = MPI_Comm_rank(comm, &rk);
	 ierr = MPI_Comm_size(comm, &size);
	 parm_param * pp;
	/*Differents preconditionners use*/
	 meth[0]=(char *)malloc(sizeof(char)*9); strcpy(meth[0],"add_ilu0"); meth[1]=(char *)malloc(sizeof(char)*9); strcpy(meth[1],"add_ilut"); 
	 meth[2]=(char *)malloc(sizeof(char)*9); strcpy(meth[2],"add_iluk"); meth[3]=(char *)malloc(sizeof(char)*9); strcpy(meth[3],"add_arms");
	 meth[4]=(char *)malloc(sizeof(char)*9); strcpy(meth[4],"lsch_ilu0"); meth[5]=(char *)malloc(sizeof(char)*9); strcpy(meth[5],"lsch_ilut");  
	 meth[6]=(char *)malloc(sizeof(char)*9); strcpy(meth[6],"lsch_iluk"); meth[7]=(char *)malloc(sizeof(char)*9); strcpy(meth[7],"lsch_arms"); 	
	 meth[8]=(char *)malloc(sizeof(char)*9); strcpy(meth[8],"rsch_ilu0"); meth[9]=(char *)malloc(sizeof(char)*9); strcpy(meth[9],"rsch_ilut");
         meth[10]=(char *)malloc(sizeof(char)*9); strcpy(meth[10],"rsch_iluk"); meth[11]=(char *)malloc(sizeof(char)*9); strcpy(meth[11],"rsch_arms");
	 meth[12]=(char *)malloc(sizeof(char)*9); strcpy(meth[12],"sch_gilu0"); meth[13]=(char *)malloc(sizeof(char)*9); strcpy(meth[13],"sch_sgs");
	/*storage format of the matrix*/
         char pcrM[4];
	strcpy(pcrM,"csr");
	
	/*- Create Distributed Matrix dm in CSR format */
		
	CreateMat(&dm, pcrM);
	
	dm->comm->mpi_comm=comm;dm->comm->myproc=rk; dm->comm->npro=size; 

	/*------ Create PrePar /iterPar pointer */
	
		
	/*---- parameters for preconditioning and iteration*/
	if((datafile.empty())&&(param_double==NULL)&&(param_int==NULL)){
		if(dm->comm->myproc==0)
		printf("%s","We are going to set default parameters because user did not specify any one  \n \n ");
		
		parm_param * pp= new parm_param(); 
		 
		iov=pp->iov; scale=pp->scale; unsymm=pp->unsymm; 
		
		method= assignprecon(meth[pp->method], dm);
		
		prepar=pp->prepar; ipar=pp->ipar; VERBOSE=pp->VERBOSE; solver=pp->solver;
	
	}
	if(((param_double!=NULL)||(param_int!=NULL))&&(datafile.empty()))
	{
		if(dm->comm->myproc==0)
  		  printf("%s","User have set parameter inside vector of parameter  \n");	

	        parm_param * pp= new parm_param(param_int, param_double); 
		iov=pp->iov; scale=pp->scale; unsymm=pp->unsymm; method= assignprecon(meth[pp->method], dm); 
		prepar=pp->prepar; ipar=pp->ipar; VERBOSE=pp->VERBOSE; solver=pp->solver;
              
	}
	if((!datafile.empty())&&((param_double==NULL)&&(param_int==NULL)))
	{
		if(dm->comm->myproc==0)
		printf("%s","User have set parameter inside file of parameter  \n");
		parm_param * pp= new parm_param(datafile, dm); 
		iov=pp->iov; scale=pp->scale; unsymm=pp->unsymm; method=pp->method; 
		prepar=pp->prepar; ipar=pp->ipar; VERBOSE=pp->VERBOSE; solver=pp->solver;
	}
	if(((solver==1)||(solver==2))&&((method!=0)||(method!=4))) 
	{
		if(dm->comm->myproc==0) 
		printf("%s%s%s", "WE NOT GARANTI THE PRECONDITIONNER WILL WORK BECAUSE ACCELARATOR ", meth[method] ," NO NEED INNER ITERATION \n");
		//MPI_Finalize();
	}
	if((dm->comm->myproc==0)&&(VERBOSE>=0)){
		printf("###########################################\n");
		printf("######### CALLING PARMS PACKAGE ###########\n");
		if(solver==0)
		printf("########### SOLVER : FGMRES #######\n");
		if(solver==1)
		printf("########### SOLVER : BiCGStab #######\n");
		if(solver==2)
		printf("########### SOLVER : DGMRES #######\n");
		printf("%s%s%s","########### PRECONDITIONNER : ",  meth[method], "\n" );
		printf("###########################################\n");
	}

	/*----from zero-based to 1-based before calling pARMS routine----*/
	n=AA.n; nnz=AA.nbcoef;
	PARMS_malloc(pr,n+1,int) ;
	PARMS_malloc(p,nnz,int) ;
	PARMS_malloc(AAv,nnz,double) ;

	for(i=0;i<nnz;i++){
	AAv[i]=AA.a[i];
	p[i]=AA.cl[i];
	if(i<=n) pr[i]=AA.lg[i];
	}
        
	for(j=0; j<pr[n]; j++){
		p[j]+=1;
	}	
	for(i=0; i<=n; i++){
		pr[i]+=1;
	}
//	pr[0]=1;

	double * vals1;
	int * id_rows1, *ptr1;
	/*-------- transpose the matrix to CSC format--*/
	PARMS_malloc(vals1,nnz,double) ;
	PARMS_malloc(id_rows1,nnz,int) ;
	PARMS_malloc(ptr1,n+1,int) ;
	job=1;
	i=1;
	
	csrcsc(&n, &job, &i, AAv, p, pr, vals1, id_rows1, ptr1);
	
 /*----------- compute C=A+B where B=transpose(A) -------*/	
	
	 if(unsymm) {
	
	PARMS_malloc(iwork,n,int);
	PARMS_malloc(vals2, 2*nnz, double);
	PARMS_malloc(id_rows2, 2*nnz, int);
	PARMS_malloc(ptr2, n+1, int);
	
	for(i=0; i<nnz; i++) AAv[i]=0.0;
	job=1;
	i1=2*nnz; 

	aplb(&n, &n, &job, vals1, id_rows1, ptr1, AAv, p, pr, 
		  vals2, id_rows2, ptr2, &i1, iwork, &ierr);

	

	if(ierr) printf("error: in aplb, ierr=%d\n",ierr);
	nnz=ptr2[n]-1;
	if (!rk)
        printf("Matrix pattern has been symmetrized; added %d nonzeros.\n",
           ptr2[n]-pr[n]);

	memcpy(pr, ptr2, (n+1)*sizeof(int));
	p=(int *)realloc(p, nnz*sizeof(int));
	AAv=(double *)realloc(AAv, nnz*sizeof(double));
	memcpy(p, id_rows2, nnz*sizeof(int));
	memcpy(AAv, vals2, nnz*sizeof(double));
	free(vals1); free(vals2); free(id_rows1); free(id_rows2); free(ptr1); free(ptr2);
	
	}

	 else {
/*-------- simply overwrite A with B (transposed matrix)  */
          nnz = ptr1[n]-1;
      memcpy(pr, ptr1, (n+1)*sizeof(int));
      memcpy(p, id_rows1, nnz*sizeof(int));
      memcpy(AAv,  vals1, nnz*sizeof(double));
    }
	

	/*------ scale the matrix */
	  if(scale) {
	    job = 1;
	    tmp = 2; /*-- compute 2-norm */
	    PARMS_malloc(scaletmpr,n,double) ; PARMS_malloc(scaletmpc,n,double) ;
	    roscal(&n,&job,&tmp,AAv,p,pr,scaletmpr,AAv,p,pr,&ierr);
	    if (ierr) fprintf(stderr, "Error: in roscal, ierr = %d\n", ierr);
	
	    coscal(&n,&job,&tmp,AAv,p,pr,scaletmpc,AAv,p,pr,&ierr);
	    if (ierr) fprintf(stderr, "Error: in coscal, ierr = %d\n", ierr);
	  } /*--- end of branch on scaling */

	/* -------Matrix partitioning : Use Metis ----------*/
	
	PARMS_malloc(riord, n, int);
	int wgtflag=0, numflag=1, volume;
	
	int option[5];
	option[0]=0;
	if(size>1){
		METIS_PartGraphVKway(&n, pr, p, NULL, NULL, &wgtflag, &numflag,&size, option, &volume, riord);
	}
	else if(size==1){
		for (i=0; i<n; i++) riord[i]=1;
	}
	
	PARMS_malloc(iwork, n, int);
	PARMS_malloc(maptmp, n, int);
	PARMS_malloc(mapptr, size+1, int);
	PARMS_malloc(iwork1, size+1, int);	
	for(i=0; i<size; i++){
		iwork[i]=0;
	}
	for(j=0; j<n; j++){
		iwork[riord[j]-1]++;
		
	}
	numflag=1;
	for(i=0; i<size; i++){
		mapptr[i]=numflag;
		numflag+=iwork[i];
	}
	mapptr[size]=numflag;
	
	for (i=0; i<size; i++){
		iwork[i]=mapptr[i]-1;
	}
	for(i=0; i<n; i++){
		maptmp[iwork[riord[i]-1]]=i+1;
		iwork[riord[i]-1]++;
	}
  	
	
	if(iov == 0) {
/*------ in non-overlapping case, simply copy the node-to-processor info */
 	   maptmp1= maptmp ;
   	   mapptr1 = mapptr;
  	}

	if(iov != 0) {
/*----- expand sub-domain if overlapping is desired */
	int maxmp;
   	 maxmp = 15*n;
	
    	 PARMS_malloc(maptmp1,maxmp,int) ;
    	 PARMS_malloc(mapptr1,size+1,int) ;
         expnddom(&n, &size, p, pr, maptmp, mapptr, maptmp1, &maxmp,mapptr1, iwork);   	
	free(maptmp); free(mapptr);
  }
	
	
/*---map from global node lables to local ---*/
	getmap(dm, maptmp1, mapptr1, &n);	
/*---prepare for extraction of local submatrix  --*/
	int tmp0=1;
	tmp =1;
	i1=mapptr1[rk];
	i2=mapptr1[rk+1]-1;
	nloc=i2-i1+1;
	nnz=0;

	for(i=i1; i <=i2; i++){
		node=maptmp1[i-1];
		nnz += pr[node]-pr[node-1];
		
	}
	PARMS_malloc(vals1, nnz, double);
	PARMS_malloc(id_rows1, nnz, int);
	PARMS_malloc(ptr1, nloc+1, int);
	/*---- extract the submatrix to be owned by rk */
	

	DPERM1(&i1,&i2,AAv,p,pr,vals1,id_rows1,ptr1,maptmp1,&tmp0,&tmp);
	
	
	

	 free(p); p=NULL; free(pr); pr=NULL; free(AAv); AAv=NULL;

/* --- CSR to Distributed matrix structure ---*/
	
	 CopyCsrToDm(dm, vals1, id_rows1, ptr1);
	 free(vals1); free(id_rows1); free(ptr1);
	 vals1=NULL; id_rows1=NULL; ptr1=NULL;
	
/*--- create boundary information */
 	 bdry(dm);

/*--- set up the local data strucutre for the sparse matrix */
	 setup(dm);
	
	
	
/*check consistency of input parameters for 'rsch*' and 'lsch*' preconditioners.
	 -- at least one iteration on inner iterative solver should be performed for
	 -- them, e.g., inner convergence tolerance; inner number of iterations and
	 -- the inner subspace size should be nonzero.
*/
	
	
	 if (!strncmp(meth[method], "lsch", 4) || !strncmp(meth[method], "rsch", 4) ){
		 ierr = 0;
		 if (ipar->pgfpar[0] == 0.0){
			 if (rk == 0)
				 fprintf(stderr, "Error: Tolerance for inner solver\n");
			 ierr = 1;
		 }
		 if (ipar->ipar[3] == 0){
			 if (rk == 0)
				 fprintf(stderr, "Error: Krylov subspace size for inner solver\n");
			 ierr = 1;
		 }
		 if (ipar->ipar[4] == 0){
			 if (rk == 0)
				 fprintf(stderr, "Error: Maximum number of inner iterations\n");
			 ierr = 1;
		 }
		 if (ierr == 1){
			 if (rk == 0)
				 fprintf(stderr,"should be nonzero to invoke the Schur Complement iteration\n ");
			 PARMS_Final();
			 exit(1);
		 }
	 }
	
	
	 /*----- create preconditioner */
	if(VERBOSE==3){
	
	  MPI_Barrier(dm->comm->mpi_comm);
	 t1 = dwalltime();}
	
	 ierr = CreatePrec(dm,&precon,method,prepar,ipar);
	
	 
	 if(VERBOSE==3) {
	  MPI_Barrier(dm->comm->mpi_comm);
	 t2 = dwalltime();
	 t2 = fabs(t2-t1);
	 double tmax=0.0;
          MPI_Reduce(&t2, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0,dm->comm->mpi_comm);
   	  MPI_Bcast(&tmax, 1, MPI_INT,0,dm->comm->mpi_comm);
	t2=tmax;
	 }
	 
	 /*----- check for errors in preconditioning construction, quit if any */
	 if (ierr != 0)
	 fprintf(stderr, "Fatal Error (ierr = %d) in preconditioner construction on processor %d\nExiting...\n", ierr, rk);
	 ierr = ierr*ierr;
	 int tmp_ierr=0;

	 MPI_Allreduce(&ierr, &tmp_ierr, 1, MPI_INT, MPI_SUM, dm->comm->mpi_comm);
	 
	 if ( tmp_ierr > 0 ){ 
		 /* delete distributed local matrix */
		 DeleteMat(&dm);
		 /* delete distributed local vector */	 
		 free(prepar);
		 free(ipar);
		 free(pp); 
	         
		 /* exit PARMS and MPI environment */
		 PARMS_Final();
		 exit(1);
	 }	
	
}
	void Solver(const MatriceMorse<double> &AA,KN_<double> &x,const KN_<double> &b) const  {
	/* Vec structure
	 rhs -- right hand side,
	 sol -- solution,   
	 xsol -- exact solution */
	
         Vec 	sol, rhs;
	 int i,comt,node;
	 int * poloc;
	 double res1;
	 double 	dnnz_pre, dgprec_nnz;
	 nnz=AA.nbcoef;	

	 CreateVec(&rhs);
	 CreateVec(&sol);
	 n=AA.n;
	  double * rhsb=(double *)malloc(sizeof(double)*n);
	 for(i=0;i<n;i++) rhsb[i]=b[i];
	 /* Copy communication structure to Vec */
	 CopyComm(dm,rhs);
	 CopyComm(dm,sol);
	 
	 PARMS_malloc(rhs->vec, nloc, double); PARMS_malloc(poloc, nloc, int);
	comt=0;
	 if(scale)
	     { for(i = 0; i < n; i++) 
		rhsb[i] *= scaletmpr[i];
		
	}
	

	for(i=i1;i<=i2;i++)
	{
		  node = maptmp1[i-1];poloc[comt]=node-1;
		  rhs->vec[comt]=rhsb[node-1];comt++;
		  
	}
/*----Iteration count-------------------------------------*/
	ipar->in_iters=0;
	ipar->iters=0;
	
/*--- Permute RHS according to result of setup routine ---*/
	 setuprhs(rhs);
	   
/*------- populate the initial guess with values */
	 VecSetVal(sol,0.0); 
/*------- calculate the norm of the residual on input */
	 res1 = ResiNorm2(dm, sol, rhs);
	 

	 if(VERBOSE==3){
	MPI_Barrier(dm->comm->mpi_comm);
	 t3 = dwalltime(); }
	double dgnnz;
	if(VERBOSE==3){
	dnnz_pre = precon->nnz_pre;
        dgnnz = (double)nnz;
	dgprec_nnz=0;
	 MPI_Reduce(&dgprec_nnz, &dnnz_pre, 1, MPI_DOUBLE, MPI_SUM, 0, dm->comm->mpi_comm); 
         }


	 if(rk == 0) 	{
			printf("Total NNZ(PreCon) / (NNZ*npro) = %6.2f\n",dgprec_nnz/dgnnz);
			
			
		}
		

	if(solver==0) fgmresd(dm,precon,ipar,rhs,sol);
	if(solver==1) dgmresd(dm,precon,ipar,rhs,sol);
	if(solver==2) bcgstabd(dm,precon,ipar,rhs,sol);
	
	
	 res = ResiNorm2(dm, sol, rhs);
	

	

	if(VERBOSE==3){	
	MPI_Barrier(dm->comm->mpi_comm);
	 t4 = dwalltime(); 
	 t4 = fabs(t4-t3); 
	  double tmax=0;
         MPI_Reduce(&t4, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0,dm->comm->mpi_comm);
   	 MPI_Bcast(&tmax, 1, MPI_DOUBLE,0,dm->comm->mpi_comm);
	t4=tmax;
	 }
	 
	 
	 /*----find the norm of the solution error */
	 i=1;
	 t3 = -1.0;
	 int j;
	  for (i=0; i<nloc; i++){
       		j=sol->node[i]-1; /* get the global node corresponding to node i*/
		poloc[i]=j;
    	 }
	
	 if (scale){
/*--------- apply permutations to scaletmp */
   	for (i=0; i<nloc; i++){
       		j=sol->node[i]-1; /* get the global node corresponding to node i*/
	        sol->vec[i] = sol->vec[i]*scaletmpc[j];
	/*---- find the residual of the computed solution */  
	}
  }	
	 /*-----  compute the relative error of computed solution */

	 if((dm->comm->myproc == 0)&&(VERBOSE==3)) {
		 fprintf(stdout,"################   SOLVER STATISTICS     ####################\n");
		 fprintf(stdout, "OUTER ITERATIONS COUNT IS %d\n", ipar->iters);
		 fprintf(stdout, "INNER ITERATION COUNT IS %d\n", ipar->in_iters);
		 fprintf(stdout, "THE TOTAL TIME IS %16.8f\n", t2+t4);
		 fprintf(stdout, "THE TIME FOR CREATING PRECONDITIONER IS %16.8f\n", t2);
		 fprintf(stdout, "THE TIME FOR SOLVING PROCESS is %16.8f\n", t4);
		 fprintf(stdout, "The 2-NORM OF THE RELATIVE RESIDUAL IS %16.8g\n", res/res1);
	 }
	 if((dm->comm->myproc == 0)&&(VERBOSE==2)) {
		 fprintf(stdout, "THE TOTAL TIME IS %16.8f\n", t2+t4);
		 fprintf(stdout, "THE TIME FOR CREATING PRECONDITIONER IS %16.8f\n", t2);
		 fprintf(stdout, "THE TIME FOR SOLVING PROCESS is %16.8f\n", t4);
	 }
	 if((dm->comm->myproc == 0)&&(VERBOSE==1)) {
		 fprintf(stdout, "OUTER ITERATIONS COUNT IS %d\n", ipar->iters);
		 fprintf(stdout, "INNER ITERATION COUNT IS %d\n", ipar->in_iters);
		 fprintf(stdout, "The 2-NORM OF THE RELATIVE RESIDUAL IS %16.8g\n", res/res1);
		 
		
	 }

  	  double *  xx= (double *)malloc(sizeof(double)*n);
	 
	 
	   comt=0;
	    
	    
	   for(i=0;i<dm->comm->npro;i++) mapptr[i]--;
	   int *displs, *perm; 
	   PARMS_malloc(displs, nloc, int);PARMS_malloc(perm, n, int);
	 
	
	   MPI_Gatherv(&(sol->vec[0]), nloc, MPI_DOUBLE, &(xx[0]), iwork, mapptr ,  MPI_DOUBLE, 0,comm );
	   MPI_Gatherv(&(poloc[0]), nloc, MPI_INT, &(perm[0]), iwork, mapptr ,  MPI_INT, 0,comm );
	   MPI_Bcast(perm,AA.n,MPI_INT,0, comm);
	   
	   int *invp=(int *)malloc(sizeof(int)*n);
	   for(i=0;i<n;i++) invp[perm[i]]=i;

	   if(dm->comm->myproc==0){for(i=0;i<n;i++) x[i]=xx[invp[i]];  } 
	   MPI_Bcast(x,AA.n,MPI_DOUBLE,0, comm);
	   for(i=0;i<dm->comm->npro;i++) mapptr[i]++;
	  /*Delete use vectors*/
	   DeleteVec(&sol); DeleteVec(&rhs); free(xx);
	   //This should be in Destructor
	   free(perm); free(invp);
	  	
	
	}
	~dSolvePARMS()
	{
	 if(VERBOSE==3){
	  cout << "~SolvePARMS:" << endl;
	 free(mapptr);
	 DeletePrec(precon); free(scaletmpc); free(scaletmpr);
	 /*---- Delete distributed local matrix */
	 DeleteMat(&dm);
	 /*---- Delete distributed local vector */
	free(prepar); free(ipar);
	/*Delete matrix and right hand side*/
	
	}
	 
	  //PARMS_Final();

         }
     void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
     {  
      ffassert(x.N()==Ax.N());
      Ax +=  (const MatriceMorse<Complex> &) (*this) * x; 
     }
};


inline MatriceMorse<double>::VirtualSolver *
BuilddSolvePARMS(DCL_ARG_SPARSE_SOLVER(double,A))
{
	 if(verbosity>9)
      cout << " BuildSolverMUMPSmpi<double>" << endl;

	

    return new dSolvePARMS(*A,ds.data_filename, ds.lparams, ds.dparams,(MPI_Comm *)ds.commworld);
}

class Init { public:
    Init();
};

//  the 2 default sparse solver double and complex
DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; 

// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetDefault()
{
    if(verbosity>1)
    cout << " SetDefault sparse to default" << endl;
    DefSparseSolver<double>::solver =SparseMatSolver_R;
    TypeSolveMat::defaultvalue =TypeSolveMat::SparseSolver;
    return TRUE;
}

bool SetdSolvePARMS()
{
    if(verbosity>1)
	cout << " SetDefault sparse solver to PARMS" << endl;
      DefSparseSolver<double>::solver  =BuilddSolvePARMS;
     TypeSolveMat::defaultvalue  = TypeSolveMatdefaultvalue;
     return TRUE;
}
LOADINIT(Init);
Init::Init()
{ 
  
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  //SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  
  if(verbosity>1)
    cout << "\n Add: pARMSmpi,  defaultsolver defaultsolverpARMSmpi" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  DefSparseSolver<double>::solver =BuilddSolvePARMS;
  //DefSparseSolver<Complex>::solver =BuildSolverMUMPSmpi;
  if(! Global.Find("defaultsolver").NotNull() )
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  Global.Add("defaulttopARMSmpi","(",new OneOperator0<bool>(SetdSolvePARMS));
}


