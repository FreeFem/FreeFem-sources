// ORIG-DATE: 02/2009
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : INRIA FUTUR 
// AUTHOR   : Guy Atenekeng
// E-MAIL   : Guy_Antoine_Atenekeng_Kahou@lri.fr
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

// #ifdef __cplusplus
// extern "C" {
// #include "metis.h"
// #endif
// #ifdef __cplusplus
// }
// #endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#ifdef MPI_WTIME_IS_GLOBAL
#undef MPI_WTIME_IS_GLOBAL
#endif
#define MPI_WTIME_IS_GLOBAL 1
#define STATS

#include "HYPRE.h"
#include "HYPRE_krylov.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#define CHECK_ZERO
#define MCW MPI_COMM_WORLD
#define BUFLEN 100
#define SCALE 0


#ifdef SUN 
/*
 * 	It uses the system call gethrtime(3C), which is accurate to 
 *	nanoseconds. 
*/
#include <sys/time.h>
 
double dwalltime() {
			 return ( (double)gethrtime() / 1e9 );
				 }

#else

#ifndef NO_TIMER
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <sys/time.h>
#endif

#ifndef CLK_TCK
#define CLK_TCK 60
#endif

				 double dwalltime()
		 {
#ifdef NO_TIMER
			 /* no sys/times.h on WIN32 */
			 double tmp;
			 tmp = 0.0;
#else
			 struct tms use;
			 double tmp;
			 times(&use);
			 tmp = use.tms_utime;
			 tmp += use.tms_stime;
#endif
			 return (double)(tmp) / CLK_TCK;
		 }

#endif



#ifdef __cplusplus
extern "C" {
#include "metis.h"
#endif
#ifdef __cplusplus
}
#endif










typedef struct sparse_mat_loc
{
        int *ptr;  //index of the beginning of rows in id_cols and vals
        int *rows; //index of non empty rows
        int *id_cols;
        double *vals;
	int ilower; // lower index limit
	int iupper; //upper index limit
        int n_loc ; /*number of rows*/
        int *ncols; /*number of columns in each row*/
}sparse_mat_loc;


/* 
  **Function to distribute a sparse matrix as blocks rows on several processes**
  
  A: (input) sparse matrix.
    (matrix A as input is available on each process
  type :(input) 0=CSR format, any other value=CSC
  size: (input) size of the communicator
  rk: (input) rank of the process in the communicator
  A_loc: (output) sparse matrix in CSR reduced on local process
  */
  
int dist_matrix(int n, int *ptr, int* id_rows, double *vals, int type, int size, int rk,int * mapptr,int *maptmp, int *iwork1 ,sparse_mat_loc *A_loc)
{
    int     i, j, ilower, iupper;
	int     n_loc, n_loc_fst, nnz_loc, relpos;
	int     *marker;
	
	/* Compute the number of rows to distribute to local process */
    n_loc =iwork1[rk];
   n_loc_fst=n_loc;
    nnz_loc=0;
   int i1,i2,node,ncols;	
    //save the lower index (ilower) and upper (iupper) for each process
    (A_loc)->ilower=mapptr[rk];
    ilower=(A_loc)->ilower;
    (A_loc)->iupper=mapptr[rk+1]-1; 
    iupper=(A_loc)->iupper;
    (A_loc)->n_loc=n_loc;
	
        if( !((A_loc)->ptr=(int *)malloc((n_loc+1)*sizeof(int))) ) {printf("%s","Malloc fails for ptr \n"); exit(1);}
	if( !((A_loc)->rows=(int *)malloc((n_loc)*sizeof(int))) ) {printf("%s","Malloc fails for rows \n");exit(1);}
	if( !((A_loc)->ncols=(int *)malloc((n_loc)*sizeof(int))) ) {printf("%s","Malloc fails for ncols \n");exit(1);}
	
    //Change global Input matrix (A) to local (A_loc) on each process
	if(type==0){  //Matrix A is in CSR format
		
		//Gets local nnz 
		 i1=mapptr[rk];
	         i2=mapptr[rk+1]-1;
		 (A_loc)->ilower=i1; (A_loc)->iupper=i2;
		 for(i=i1;i<=i2;i++){
			node=maptmp[i]; nnz_loc+=ptr[node+1]-ptr[node];
		 }
		//Allocate memory for local matrix
		if( !((A_loc)->id_cols=(int *)malloc(nnz_loc*sizeof(int))) ) {printf("%s","Malloc fails for id_cols \n");exit(1);}
		if( !((A_loc)->vals=(double *)malloc(nnz_loc*sizeof(double))) ) {printf("%s","Malloc fails for vals"); exit(1);}
		
		//Transfer the corresponding values from global to local
		relpos=0; 
		//int ncols; //count number of elements in each row
		for(i=i1;i<=i2;i++){
			(A_loc)->rows[i-i1]=maptmp[i];
			(A_loc)->ptr[i-i1]=relpos;
			ncols=relpos;
			for(j=ptr[maptmp[i]];j<ptr[maptmp[i]+1];j++){
				(A_loc)->id_cols[relpos]=id_rows[j];
				(A_loc)->vals[relpos]=vals[j];
				relpos++;
			}
			(A_loc)->ncols[i-i1]=relpos-ncols;
		}
	}
	else{ //matrix A is in CSC format
		marker= (int *)calloc(n, sizeof(int)); //count number of elements in each row
		for(i=0; i<n; i++)
			for(j=ptr[i]; j<ptr[i+1]; j++)
				marker[id_rows[j]]++;
		
		(A_loc)->ptr[0]=0; //set up the beginning of each row
		for(i=0; i<n_loc; i++){
			(A_loc)->ptr[i+1] = (A_loc)->ptr[i] + marker[i+ilower];
				(A_loc)->id_cols[relpos]=id_rows[j];
				(A_loc)->vals[relpos]=vals[j];
				relpos++;
			}
			(A_loc)->ncols[i-ilower]=relpos-ncols;
		}
    return 0;
}




class hypreParam {
 //Solveur and preconditionner
  public: 
   char  solver[BUFLEN];
   char  precon[BUFLEN];
 //BoomerAMG parameter
   int                  amg_coarsentype ;    /* Falgout coarsening */
   int                  amg_relaxtype;      /* hybrid Gauss-Seidel or SOR */ 
   int			amg_interptype;		/* default*/
    int 			amg_maxlevels;
   int                  amg_numsweeps;        /*default*/
   double               amg_strongthreshold;/*suitable for 3D Laplace Operator*/
   double		amg_truncfactor; 
   int                  amg_prntlevel;        /* print setup info */
   int			amg_tol;			//BoomerAMG Tolerance
   int			amg_maxiter;
   
   //More complex smoothers (Schwarz methods, Pilut, Parasails, Euclid)
   int					smooth_type;
   int 					smooth_numlevels; 
   int 					smooth_numsweeps;
   double				pilut_droptol;
   double 				pilut_maxnz;
   int					schwarz_overlap;
   int 					schwarz_variant;
   int					schwarz_domaintype;
   
   //parasails parameter
    int                  sai_max_levels ;
   double               sai_threshold ;

   double               sai_filter ;
   int                  sai_sym ;
   int                  sai_log ;
   int                  VERBOSE;
   
   //Euclid Parameter
   char                 euclid_param[13]; /*Edit Euclid parameter in this file*/
   
   //Solver parameter (used for GMRES , PCG or BiCGStab)
    double               solv_tol ;
   int                  solv_maxiter;
   int                  solv_kdim;
   int                  solv_log;
   int                  solv_prntlevel;
   int                     precond_id, solver_id;
 public : hypreParam(const KN<long> &param_int, const KN<double> &param_double)
	{


		amg_coarsentype =6;    /* Falgout coarsening */
   		amg_relaxtype =3;      /* hybrid Gauss-Seidel or SOR */ 
	 	amg_interptype=0;		/* default*/
   		amg_maxlevels=25;
   		amg_numsweeps=1;        /*default*/
		amg_strongthreshold = 0.25;/*suitable for 3D Laplace Operator*/
   		amg_truncfactor=0.3; 
		amg_prntlevel =1;        /* print setup info */
		amg_tol=0.0;			//BoomerAMG Tolerance
		amg_maxiter=1;
   //More complex smoothers (Schwarz methods, Pilut, Parasails, Euclid)
		smooth_type=9;
		smooth_numlevels=3; 
		smooth_numsweeps=1;
		pilut_droptol=1.0e-4;
 		pilut_maxnz=100;
		schwarz_overlap=0;
		schwarz_variant=2;
		schwarz_domaintype=2;
   //parasails parameter
  		sai_max_levels = 1;
		sai_threshold = 0.1;
		sai_filter = 0.05;
	        sai_sym =0;
	        sai_log = 1;
		strcpy(euclid_param,"euclid_param");		
   //Solver parameter (used for GMRES or BiCGStab)
   		solv_tol = 1.0e-11;
		solv_maxiter = 1000;
		solv_kdim =40;

		solv_log = 0;
		solv_prntlevel = 2;
		precond_id=0;//BOOMER AMG 	
		solver_id=1; //GMRES as solver
	        VERBOSE=0;


		
		
		if(param_int.N()>0) {if((param_int[0]>=0)&&(param_int[0]<=2)) solver_id=param_int[0]; else solver_id=1;} //BOOMER AMG }
		if(param_int.N()>1) {if((param_int[1]>=0)&&(param_int[1]<=2)) solver_id=param_int[1]; else precond_id=0;}//GMRES as solver
		if(param_int.N()>2) {if(param_int[2]>0) solv_maxiter = param_int[2];else solv_maxiter=1000;}
		if(param_int.N()>3) {if(param_int[3]>0) solv_kdim =param_int[3];}
		if(param_double.N()>4) {if(param_double[4]>0) solv_tol = param_double[4];}
		if(param_int.N()>4) {if(param_int[4]>=0) solv_log = param_int[4];}
		if(param_int.N()>5) {if(param_int[5]>=0) solv_prntlevel = param_int[5];}
		if(param_int.N()>6) {if(param_int[6]>=0) amg_coarsentype =param_int[6];}    /* Falgout coarsening */
		if(param_int.N()>7) {if(param_int[7]>=0) amg_relaxtype =param_int[7];}
		if(param_int.N()>8) {if(param_int[8]>=0) amg_interptype =param_int[8];}
		if(param_int.N()>9) {if(param_int[9]>=0) amg_maxlevels =param_int[9];}
		if(param_int.N()>10) {if(param_int[10]>=0) amg_numsweeps =param_int[10];}

		if(param_double.N()>5) {if(param_double[5]>0) amg_strongthreshold = param_double[5];}
		if(param_double.N()>6) {if(param_double[6]>0) amg_truncfactor = param_double[6];}
		if(param_int.N()>11) {if(param_int[11]>=0) amg_prntlevel= param_int[11];}
		if(param_int.N()>12)  {if(param_int[12]>0) amg_maxiter=param_int[12];}
		if(param_double.N()>7) {if(param_double[7]>0) amg_tol=param_double[7];}
		
		
		if(param_int.N()>13) {if(param_int[13]>0) smooth_type=param_int[13];}
		if(param_int.N()>14) {if(param_int[14]>0) smooth_numlevels=param_int[14];}
		if(param_int.N()>15) {if(param_int[15]>0) smooth_numsweeps=param_int[15];}
		if(param_double.N()>8) {if(param_double[8]>0) pilut_droptol=param_double[8];}
		if(param_int.N()>16) {if(param_int[16]>0) pilut_maxnz=param_int[16];}
		
		if(param_int.N()>17) {if(param_int[17]>0) schwarz_overlap=param_int[17];}
		if(param_int.N()>18) {if(param_int[18]>0) schwarz_variant=param_int[18];}
		if(param_int.N()>19) {if(param_int[19]>0) schwarz_domaintype=param_int[19];}
   //parasails parameter
		
  		if(param_int.N()>20) {if(param_int[20]>0) sai_max_levels = param_int[20];}
		if(param_double.N()>9) {if(param_double[9]>0) sai_threshold = param_double[9];}
		if(param_double.N()>10) {if(param_double[10]>0) sai_filter = param_double[10];}
	        if(param_int.N()>21) {if(param_int[21]>0) sai_sym =param_int[21];}
	        if(param_int.N()>22) {if(param_int[22]>0) sai_log = param_int[22];}
		if(param_int.N()>23) {if(param_int[23]>0) VERBOSE = param_int[23];}
		//strcpy(euclid_param,"euclid_param");		
                
	}	
  public : hypreParam()
	{
		amg_coarsentype =6;    /* Falgout coarsening */
   		amg_relaxtype =3;      /* hybrid Gauss-Seidel or SOR */ 
	 	amg_interptype=0;		/* default*/
   		amg_maxlevels=25;
   		amg_numsweeps=1;        /*default*/
		amg_strongthreshold = 0.25;/*suitable for 3D Laplace Operator*/
   		amg_truncfactor=0.3; 
		amg_prntlevel =1;        /* print setup info */
		amg_tol=0.0;			//BoomerAMG Tolerance
		amg_maxiter=1;
   //More complex smoothers (Schwarz methods, Pilut, Parasails, Euclid)
		smooth_type=9;
		smooth_numlevels=3; 
		smooth_numsweeps=1;
		pilut_droptol=1.0e-4;
 		pilut_maxnz=100;
		schwarz_overlap=0;
		schwarz_variant=2;
		schwarz_domaintype=2;
   //parasails parameter
  		sai_max_levels = 1;
		sai_threshold = 0.1;
		sai_filter = 0.05;
	        sai_sym =0;
	        sai_log = 1;
		strcpy(euclid_param,"euclid_param");		
   //Solver parameter (used for GMRES or BiCGStab)
   		solv_tol = 1.0e-11;
		solv_maxiter = 1000;
		solv_kdim =40;
  //  int                  solv_stopcrit = 1; //only for BiCGSTAB
		solv_log = 0;
		solv_prntlevel = 2;
		precond_id=0;//BOOMER AMG 	
		solver_id=2; //GMRES as solver
		
	}
  public : hypreParam(char * fileparameter, MPI_Comm comm)
	{
		FILE *f;
        	char buf[BUFLEN];
		int num;
		int rk, size;
		MPI_Comm_rank(comm,&rk);
		MPI_Comm_size(comm, &size);

		
		amg_coarsentype =6;    /* Falgout coarsening */
   		amg_relaxtype =3;      /* hybrid Gauss-Seidel or SOR */ 
	 	amg_interptype=0;		/* default*/
   		amg_maxlevels=25;
   		amg_numsweeps=1;        /*default*/
		amg_strongthreshold = 0.25;/*suitable for 3D Laplace Operator*/
   		amg_truncfactor=0.3; 
		amg_prntlevel =1;        /* print setup info */
		amg_tol=0.0;			//BoomerAMG Tolerance
		amg_maxiter=1;
   //More complex smoothers (Schwarz methods, Pilut, Parasails, Euclid)
		smooth_type=9;
		smooth_numlevels=3; 
		smooth_numsweeps=1;
		pilut_droptol=1.0e-4;
 		pilut_maxnz=100;
		schwarz_overlap=0;
		schwarz_variant=2;
		schwarz_domaintype=2;
   //parasails parameter
  		sai_max_levels = 1;
		sai_threshold = 0.1;
		sai_filter = 0.05;
	        sai_sym =0;
	        sai_log = 1;
		strcpy(euclid_param,"euclid_param");		
   //Solver parameter (used for GMRES or BiCGStab)
   		solv_tol = 1.0e-11;
		solv_maxiter = 1000;
		solv_kdim =40;
  //  int                  solv_stopcrit = 1; //only for BiCGSTAB
		solv_log = 0;
		solv_prntlevel = 0;
		precond_id=0;//BOOMER AMG 	
		solver_id=1; //GMRES as solver
		VERBOSE=0;
		



		if(fileparameter==NULL)
		{
			if(rk==0) printf("%s","Set default parameter because you not precise the parameter file \n \n");
			solver_id=1; //GMRES as solver
			precond_id=0;//BOOMER AMG 
		
		}
		else if( (f=fopen(fileparameter,"r") )==NULL)
			{
			if(rk==0)printf("%s","Set default parameter because your parameter file not exist \n \n");
			solver_id=1; //GMRES as solver
			precond_id=0;//BOOMER AMG 	
			}	
			else
			{
			 if(rk==0) printf("%s%s%s","Read parameter from file ", fileparameter, "\n \n");
			 num =0;
			 while(fgets(buf, BUFLEN, f) != NULL) 
			 {
				switch(num) {
				case 0:
				sscanf(buf, "%s", solver); 
				break;
				case 1:
 				sscanf(buf, "%d", &solver_id);
			        break;
				case 2:
 				sscanf(buf, "%lf", &solv_tol);
			        break;
				case 3:
 				sscanf(buf, "%d", &solv_maxiter);
			        break;
				case 4:
 				sscanf(buf, "%d", &solv_prntlevel);
			        break;
				case 5:
 				sscanf(buf, "%d", &solv_log);
				if(solver_id!=1) {fgets(buf, BUFLEN, f);num++;}
				break;				
				case 6:
				sscanf(buf, "%d", &solv_kdim);
			        break;
				case 7:
				sscanf(buf, "%s", precon); 
				break;
				case 8:
				sscanf(buf, "%d", &precond_id); 
				 	
				if(precond_id==2) //The parameter of preconditionner is inside file
				{
					fclose(f);
					exit(1);
				}
				break;
				case 9:
				if(precond_id==0)
				sscanf(buf, "%d", &amg_coarsentype); 
				if(precond_id==1)
				sscanf(buf, "%lf", &sai_threshold); 
				break;
				case 10:
				if(precond_id==0)
				sscanf(buf, "%d", &amg_prntlevel); 
				if(precond_id==1)
				sscanf(buf, "%d", &sai_max_levels); 
				break;
				case 11:
				if(precond_id==0)
				sscanf(buf, "%d", &amg_interptype); 
				if(precond_id==1)
				sscanf(buf, "%lf", &sai_filter); 
				break;
				case 12:
				if(precond_id==0)
				sscanf(buf, "%d", &amg_maxlevels); 
				if(precond_id==1)
				case 17:
				sscanf(buf, "%d", &amg_prntlevel); 
				break;
				case 18:
				sscanf(buf, "%d", &amg_tol); 
				break;
				case 19:
				sscanf(buf, "%d", &amg_maxiter); 
				break;				
				default:
				break;
				}	
         		 num++;	
			 }
			 if(fclose(f)==EOF) printf("%s","Error while closing the file \n"); else printf("%s","File is well close \n");			
			}
	}
};

class HypreSolver :   public MatriceMorse<double>::VirtualSolver   {
   mutable HYPRE_IJMatrix       ij_A;
   mutable HYPRE_IJVector       ij_B;
   mutable HYPRE_IJVector       ij_X;
   mutable HYPRE_ParCSRMatrix   par_A;
   mutable HYPRE_ParVector      par_B;
   mutable HYPRE_ParVector      par_X;
   string data_option;
   mutable HYPRE_Solver         solver;
   mutable HYPRE_Solver         precond;
   mutable hypreParam           *param;
   int rk,size;
    int                  jlower, jupper;
   int                  ilower, iupper;
   mutable    int num_iter;
   mutable   double 		final_res_norm;
   mutable   double 		tgv,eps,tol_pivot,tol_pivot_sym,epsr;
   sparse_mat_loc	A_loc;
   mutable int *iwork, *maptmp, *mapptr, *iwork1,*riord;
   mutable int n_loc,n,VERBOSE;
   public:
  HypreSolver(const MatriceMorse<double> &AA,string datafile, const KN<long> &param_int, const KN<double> &param_double)  	
  { 
  int i,j;
    /****INITIALIZE MPI ENVIRONMENT***/
    MPI_Comm_rank(MCW,&rk);
    MPI_Comm_size(MCW, &size);
    n=AA.n;
   
	if((data_option.empty())&&((param_int!=NULL)&&(param_double!=NULL)))
	param=new hypreParam(param_int,param_double);
	if((data_option.empty())&&((param_int==NULL)&&(param_double==NULL)))
	param= new hypreParam();
	if((!data_option.empty())&&((param_int==NULL)&&(param_double==NULL)))
	{
	char *p;
	p=new char[data_option.length()+1]; strcpy(p,data_option.c_str()); 
	param= new hypreParam(p,MCW);
	}
	 /*###################################################
                        USING HYPRE
         ####################################################*/
	if((rk==0)&&(VERBOSE>=0)){
	   printf("###########################################\n");
	   printf("######### CALLING HYPRE PACKAGE ###########\n");
	   
	   switch(param->precond_id){
		   case 0: 
			   printf("####### PRECONDITIONER : BoomerAMG #######\n"); 
			   break;
		   case 1 :
			   printf("####### PRECONDITIONER : ParaSails #######\n");
			   break;
		   case 2 :
			   printf("####### PRECONDITIONER : EUCLID   #######\n");
			   break;
		   default: printf("##### Wrong Precond_id\n");break;
	   }
	   switch(param->solver_id){
		   case 0: 
			   printf("########### SOLVER : BiCGStab #######\n"); 
			   break;
			   
		   case 1 :
			   printf("########### SOLVER : GMRES #######\n");
			   break;
		  case 2  : 
			  printf("########### SOLVER : PCG #######\n");
			  break;
		   default: printf("##### Wrong Solver_id\n"); break;
	   }
	   printf("###########################################\n");
	}
	/*************************************************************
	Distribute input matrix into local structures
	*************************************************************/
		
		int type=0; //0=CSR; 1=CSC
		int wgtflag=0, numflag=0, volume;
        	int option[5];  option[0]=0;
		riord=(int *)malloc(sizeof(int)*AA.n);
		if(riord==NULL){if(rk==0) printf("%s","In partition,  memory allocation failed \n");exit(1);}
		if(size>1)
		METIS_PartGraphKway(&(n), AA.lg, AA.cl, NULL, NULL, &wgtflag, &numflag,&size, option, &volume, riord);
		else if(size==1){
                for (i=0; i<n; i++) riord[i]=1;
        }
		iwork=(int *)malloc(sizeof(int)*AA.n);
	        maptmp=(int *)malloc(sizeof(int)*AA.n);
		mapptr=(int *)malloc(sizeof(int)*(size+1));
		iwork1=(int *)malloc(sizeof(int)*(size+1));
	        for(i=0; i<size; i++){
	        iwork[i]=0;iwork1[i]=0; }
	        for(j=0; j<n; j++){
                iwork[riord[j]]++;
	        iwork1[riord[j]]++;
       		 }
		 numflag=0;
		for(i=0; i<size; i++){
                mapptr[i]=numflag;
                numflag+=iwork[i];
	        }
		 mapptr[size]=numflag;

	        for (i=0; i<size; i++){
                iwork[i]=mapptr[i];
       		 }
	        for(i=0; i<n; i++){
                maptmp[iwork[riord[i]]]=i;
                iwork[riord[i]]++;
       		 }
                                                                    
		dist_matrix(AA.n,AA.lg,AA.cl,AA.a, type, size, rk,mapptr,maptmp,iwork1, &A_loc); 	

		
		/***Distribute vector***/
		n_loc=A_loc.n_loc;
	/**** Preparing Matrix, X and RHS *****/
   	ilower=A_loc.ilower; jlower=ilower;
   	iupper=A_loc.iupper; jupper=iupper;
   
	HYPRE_IJMatrixCreate(MCW, ilower, iupper, jlower, jupper, &ij_A);
	HYPRE_IJMatrixSetObjectType(ij_A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(ij_A);
	
	
   
   /* Set matrix coefficients */
  	HYPRE_IJMatrixSetValues(ij_A, n_loc, A_loc.ncols, A_loc.rows, A_loc.id_cols, A_loc.vals);
	
	HYPRE_IJMatrixAssemble(ij_A);
	
   /*Extract HYPRE Object data structures suitable for the solvers*/
   	HYPRE_IJMatrixGetObject(ij_A, (void **) &par_A);
	
  /*Free local matrix in A_loc, B and X*/
        free(A_loc.rows); free(A_loc.id_cols); free(A_loc.vals); free(A_loc.ncols); 
        double t3,t4;
	
	/*********************************************
                Create preconditioner 
   ********************************************/
	if(VERBOSE>=1){
         MPI_Barrier(MCW);
         t3 = dwalltime();

        }

	

   switch(param->precond_id){
       case 0 : //Set up BoomerAMG Preconditioner
           HYPRE_BoomerAMGCreate(&precond);
	    
           HYPRE_BoomerAMGSetCoarsenType(precond, param->amg_coarsentype); 
           HYPRE_BoomerAMGSetPrintLevel(precond, param->amg_prntlevel);
           HYPRE_BoomerAMGSetRelaxType(precond, param->amg_relaxtype); 
           HYPRE_BoomerAMGSetInterpType(precond, param->amg_interptype);
	   HYPRE_BoomerAMGSetNumSweeps(precond, param->amg_numsweeps); 
           HYPRE_BoomerAMGSetStrongThreshold(precond, param->amg_strongthreshold);
	   HYPRE_BoomerAMGSetTol(precond, param->amg_tol);
	   HYPRE_BoomerAMGSetMaxIter(precond, param->amg_maxiter);
	   HYPRE_BoomerAMGSetMaxLevels(precond, param->amg_maxlevels);
	   
	   break;
       case 1 : //Set up Parasails preconditioner
           HYPRE_ParaSailsCreate(MCW, &precond);
           HYPRE_ParaSailsSetParams(precond, param->sai_threshold, param->sai_max_levels); 
           HYPRE_ParaSailsSetFilter(precond, param->sai_filter);
           HYPRE_ParaSailsSetSym(precond, param->sai_sym);
           HYPRE_ParaSailsSetLogging(precond, param->sai_log);
           break;

       case 2 : //Set up Euclid Preconditioner
           HYPRE_EuclidCreate(MCW, &precond);
           HYPRE_EuclidSetParamsFromFile(precond, param->euclid_param);
           break;
       default : 
           if(rk==0) printf("Wrong precond_id\n");
   }

	 if(VERBOSE>=1){  
	 MPI_Barrier(MCW);
         t4 = dwalltime();
	 t4=t4-t3;
	 MPI_Reduce(&t4, &t3, 1, MPI_DOUBLE, MPI_MAX, 0,MCW);
   	 MPI_Bcast(&t3, 1, MPI_DOUBLE,0,MCW);
         if((rk==0)&&(VERBOSE>=1)) printf("%s%16.8g%s","TIME FOR CREATING PRECONDITIONNER",t3,"\n");
	}
   /********************************************
                  Setup and Use solver
   *******************************************/
   	if(VERBOSE>=1){
         MPI_Barrier(MCW);
         t3 = dwalltime();
         
        }

	

   switch(param->solver_id){
       case 0 : //BiCGStab solver
            HYPRE_ParCSRBiCGSTABCreate(MCW, &solver);
            //set various parameters
            HYPRE_ParCSRBiCGSTABSetTol(solver, param->solv_tol);
            HYPRE_ParCSRBiCGSTABSetMaxIter(solver, param->solv_maxiter);
            //HYPRE_ParCSRBiCGSTABSetStopCrit(solver, solv_stopcrit);
            HYPRE_ParCSRBiCGSTABSetPrintLevel(solver, param->solv_prntlevel);
            HYPRE_ParCSRBiCGSTABSetLogging(solver, param->solv_log);
            
            //Set Preconditioner
            switch (param->precond_id){
                case 0 : //set BoomerAMG as preconditioner
                    HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup, precond);
                    break;
                case 1 : //set ParaSails as preconditioner
                    HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_ParaSailsSolve, (HYPRE_PtrToParSolverFcn) HYPRE_ParaSailsSetup, precond);
                    break;
                case 2 : //Set Euclid as preconditioner
                    HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_EuclidSolve, (HYPRE_PtrToParSolverFcn) HYPRE_EuclidSetup, precond);
                    break;
            }
           
            break;
           
       case 1 :  //GMRES Solver

	    HYPRE_ParCSRGMRESCreate(MCW, &solver);
	   
            //Set various parameters
            HYPRE_ParCSRGMRESSetMaxIter(solver, param->solv_maxiter );
            HYPRE_ParCSRGMRESSetKDim(solver,param->solv_kdim);
            HYPRE_ParCSRGMRESSetTol(solver,param->solv_tol );
            HYPRE_ParCSRGMRESSetPrintLevel(solver,param->solv_prntlevel);
            HYPRE_ParCSRGMRESSetLogging(solver, param->solv_log);      
	  
            //Set Preconditioner
            switch (param->precond_id){
                case 0 : //set BoomerAMG as preconditioner
                    HYPRE_ParCSRGMRESSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup, precond);
		    
                    break;
                case 1 : //set ParaSails as preconditioner
                    HYPRE_ParCSRGMRESSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_ParaSailsSolve, (HYPRE_PtrToParSolverFcn) HYPRE_ParaSailsSetup, precond);
                    break;
                case 2 : //Set Euclid as preconditioner
                    HYPRE_ParCSRGMRESSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_EuclidSolve, (HYPRE_PtrToParSolverFcn) HYPRE_EuclidSetup, precond);
                    break;
            }
            
            break;
	case 2:	  //PCG Solver
	
	 HYPRE_ParCSRPCGCreate(MCW, &solver);
	 HYPRE_ParCSRPCGSetMaxIter(solver, param->solv_maxiter); /* max iterations */
	 HYPRE_ParCSRPCGSetTol(solver, param->solv_tol); /* conv. tolerance */
	 HYPRE_ParCSRPCGSetPrintLevel(solver, param->solv_prntlevel); /* prints out the iteration info */
	 HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */
	//Set Preconditioner
            switch (param->precond_id){
                case 0 : //set BoomerAMG as preconditioner
                    HYPRE_ParCSRPCGSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSolve, (HYPRE_PtrToParSolverFcn) HYPRE_BoomerAMGSetup, precond);
	
		    
                    break;
                case 1 : //set ParaSails as preconditioner
                    HYPRE_ParCSRPCGSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_ParaSailsSolve, (HYPRE_PtrToParSolverFcn) HYPRE_ParaSailsSetup, precond);
                    break;
                case 2 : //Set Euclid as preconditioner
                    HYPRE_ParCSRPCGSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_EuclidSolve, (HYPRE_PtrToParSolverFcn) HYPRE_EuclidSetup, precond);
                    break;
            }
	break;
        default :
        if(rk==0) printf("Wrong solver_id\n");
	break;
	}
	
	if(verbosity==0){
         MPI_Barrier(MCW);
         t4 = dwalltime();
	 t4=t4-t3;
	 MPI_Reduce(&t4, &t3, 1, MPI_DOUBLE, MPI_MAX, 0,MCW);
   	 MPI_Bcast(&t3, 1, MPI_DOUBLE,0,MCW);
         if((rk==0)&&(VERBOSE>=1)) printf("%s%16.8g%s","TIME FOR SETTING SOLVER",t3,"\n");


        }


}
	
	void Solver(const MatriceMorse<double> &AA,KN_<double> &x,const KN_<double> &b) const  {
	 ffassert ( &x[0] != &b[0]);
	 epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
	 int i, i1,i2,node;
	 int *row,*row1;
	 double * b_loc,*X_loc,*rhs,*xx;
	 n=AA.n;
	// n_loc =n/size; n_loc_fst=n_loc;n_loc1=n/size; n_loc2=n-size*n_loc1+n_loc1;
	 rhs= (double *)malloc(sizeof(double)*n);
	 xx= (double *)malloc(sizeof(double)*n);
	 x= (double *)malloc(sizeof(double)*n);
	 i1=mapptr[rk];
        i2=mapptr[rk+1]-1;
	double t2,t1;

	 if(VERBOSE>=2) {
         MPI_Barrier(MCW);
         t1 = dwalltime();
       
         }

	for(i=0;i<n;i++) rhs[i]=b[i];
   	
	b_loc=(double *)malloc(n_loc*sizeof(double));
	X_loc=(double *)malloc(n_loc*sizeof(double));
	row=(int *)malloc(n_loc*sizeof(int));
	row1=(int *)malloc(n_loc*sizeof(int));
	for(i=i1; i<=i2; i++){
	 node = maptmp[i];
	b_loc[i-i1]=rhs[node];
	X_loc[i-i1]=0.0;  //Initial Guest
	row[i-i1]=node;  //used to get results later
	row1[i-i1]=A_loc.ilower+i-i1;
	}
     
	HYPRE_IJVectorCreate(MCW, jlower, jupper, &ij_B);
	HYPRE_IJVectorCreate(MCW, jlower, jupper, &ij_X);
	
	 HYPRE_IJVectorSetObjectType(ij_B, HYPRE_PARCSR);
	 HYPRE_IJVectorSetObjectType(ij_X, HYPRE_PARCSR);
	 HYPRE_IJVectorInitialize(ij_B);
	 HYPRE_IJVectorInitialize(ij_X);
	 
	HYPRE_IJVectorSetValues(ij_B, n_loc, row, b_loc); 
	
   	HYPRE_IJVectorSetValues(ij_X, n_loc, row, X_loc);
	  
	HYPRE_IJVectorAssemble(ij_B);
	HYPRE_IJVectorAssemble(ij_X);
       
	HYPRE_IJVectorGetObject(ij_B, (void **) &par_B); 
        HYPRE_IJVectorGetObject(ij_X, (void **) &par_X);
	
	switch (param->solver_id){
                case 0 : //BICGSTAB solver
          	  HYPRE_ParCSRBiCGSTABSetup (solver, par_A, par_B, par_X);
	          HYPRE_ParCSRBiCGSTABSolve (solver, par_A, par_B, par_X);
	          HYPRE_ParCSRBiCGSTABGetNumIterations(solver, &num_iter);
                  HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
		
		  HYPRE_ParCSRBiCGSTABDestroy(solver);
                    break;
                case 1 : //GMRES Solver

           	 HYPRE_ParCSRGMRESSetup(solver, par_A, par_B, par_X);
           	 HYPRE_ParCSRGMRESSolve(solver, par_A, par_B, par_X);
            	HYPRE_ParCSRGMRESGetNumIterations(solver, &num_iter);
            	HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);
		
		HYPRE_ParCSRGMRESDestroy(solver);
                    break;
                case 2 : //PCG
			 
			 HYPRE_ParCSRPCGSetup(solver, par_A, par_B, par_X);
			 HYPRE_ParCSRPCGSolve(solver, par_A, par_B, par_X);
			 
			 HYPRE_PCGGetNumIterations(solver, &num_iter);
     			 HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
			HYPRE_ParCSRPCGDestroy(solver);
			
                  break;
		  default :
		  break;
		}

	   int *perm=(int *)malloc(sizeof(int)*n);

	    HYPRE_IJVectorGetValues(ij_X, n_loc, row1, X_loc);
	    MPI_Gatherv(X_loc,n_loc, MPI_DOUBLE, xx,iwork1,mapptr,MPI_DOUBLE,0,MCW);
	     MPI_Gatherv(row1,n_loc, MPI_INT, perm,iwork1,mapptr,MPI_INT,0,MCW);
	    MPI_Bcast(xx,n,MPI_DOUBLE,0, MCW); MPI_Bcast(perm,n,MPI_INT,0, MCW);
	    for(i=0;i<n;i++) x[i]=xx[perm[i]];

	    if(verbosity==0) {
	     MPI_Barrier(MCW);
             t2 = dwalltime();
	     t2=t2-t1;
	     MPI_Reduce(&t2, &t1, 1, MPI_DOUBLE, MPI_MAX, 0,MCW);
	     MPI_Bcast(&t1, 1, MPI_DOUBLE,0,MCW);
             
	     if((rk==0)&&(VERBOSE>=2)){printf("%s%18.6g%s","TIME FOR SOLVING ", fabs(t1) ,"\n \n");
		       printf("%s%18.6g%s","RELATIVE RESIDU  ", final_res_norm ,"\n \n");	
			printf("%s%d%s","NUMBER OF ITERATION ", num_iter ," \n \n");
		      }
         }



	    if(X_loc!=NULL) free(X_loc);  if(rhs!=NULL) free(rhs); if(xx!=NULL) free(xx);
	    if(b_loc!=NULL) free(b_loc);  if(row!=NULL) free(row);
            }
	  ~HypreSolver() { 
    		if(verbosity){
	        cout << "HypreSolver():" << endl;
		
  		HYPRE_IJMatrixDestroy(ij_A);
		HYPRE_IJVectorDestroy(ij_B);
		HYPRE_IJVectorDestroy(ij_X);
	   
		  }
	}
 void addMatMul(const KN_<double> & x, KN_<double> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<double> &) (*this) * x; 
  }
};

inline MatriceMorse<double>::VirtualSolver *
BuildHypreSolver(DCL_ARG_SPARSE_SOLVER(double,A))
{
    if(verbosity>9)
      cout << " BuildSolverHypre>" << endl;
    return new HypreSolver(*A,ds.data_filename, ds.lparams, ds.dparams);
}


class Init { public:
    Init();
};

//  the 2 default sparse solver double and complex
DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; 
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
        return 1;
}

bool BuildHypreSolver()
{
    if(verbosity)
	cout << " SetDefault sparse solver to Hyprempi" << endl;
    DefSparseSolver<double>::solver  =BuildHypreSolver;
    TypeSolveMat::defaultvalue  = TypeSolveMatdefaultvalue;
    return 1;
}
Init init;
Init::Init()
{ 
  
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  
  
  if(verbosity>1)
    cout << "\n Add: Hyprempi,  defaultsolver defaultsolverHyprempi" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  DefSparseSolver<double>::solver =BuildHypreSolver;
  if(! Global.Find("defaultsolver").NotNull() )
    Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  Global.Add("defaulttoHyprempi","(",new OneOperator0<bool>(BuildHypreSolver));
}




