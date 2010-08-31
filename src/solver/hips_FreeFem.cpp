// ORIG-DATE: 04/2009
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : INRIA FUTUR
// AUTHOR   : Guy Atenekeng
// E-MAIL   : Guy_Antoine_Atenekeng_Kahou@lri.fr
//
//ff-c++-LIBRARY-dep: metis hips  blas 
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

#include  <iostream>
using namespace std;

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"
#include "MatriceCreuse_tpl.hpp"
#include "mpi.h"

extern "C" {
  #include "hips.h"  
 // #include "io.h"
  #include "metis.h"
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFLEN 200
#define MCW MPI_COMM_WORLD

//roscal(&n,&job,&tmp,AAv,p,pr,scaletmpr,AAv,p,pr,&ierr);

int roscal(int n, int job,int nrm, double *AAv, int *p, int *pr, double * scaletmpr , int *ierr)
{
/*---------------------------------------------------------------------
|
| This routine scales each row of mata so that the norm is 1.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SparRow form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> row j is a zero row
|--------------------------------------------------------------------*/
/*   local variables    */
   int i, k;
   double  scal;
   
   for (i=0; i<n; i++) {
      scal = 0.0;
    // kr = &(AAv[pr[i]]);
      if (nrm == 0) {
	 for (k=pr[i]; k<pr[i+1]; k++)
	    if (fabs(AAv[k]) > scal) scal = fabs(AAv[k]);
      }
      else if (nrm == 1) {
         for (k=pr[i]; k<pr[i+1]; k++)
            scal += fabs(AAv[k]);
      }
      else {  /* nrm = 2 */
         for (k=pr[i]; k<(pr[i+1]); k++)
            scal += AAv[k]*AAv[k];
      }
      if (nrm == 2) scal = sqrt(scal);
      if (scal == 0.0) {
	*ierr=i;
	 return i+1;
      }
      else 
	 scal = 1.0 / scal;
      scaletmpr[i] = scal;
      for (k=pr[i]; k<(pr[i+1]); k++)
	 AAv[k] = AAv[k] * scal;
     
   }
   *ierr=0;
   return 0;
}
/*---------------end of roscalC-----------------------------------------
----------------------------------------------------------------------*/
int coscal(int n, int job,int nrm, double *AAv, int *p, int *pr, double * scaletmpc , int * ierr)
{
/*---------------------------------------------------------------------
|
| This routine scales each column of mata so that the norm is 1.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in SparRow form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> column j is a zero column
|--------------------------------------------------------------------*/
/*   local variables    */
   int i, j, k;
   double *kr;
   int *ki;
   for (i=0; i<n; i++)
      scaletmpc[i] = 0.0;
/*---------------------------------------
|   compute the norm of each column
|--------------------------------------*/
   for (i=0; i<n; i++) {
      kr = &(AAv[pr[i]]);
      ki = &(pr[i]);
      if (nrm == 0) {
	 for (k=pr[i]; k<pr[i+1]; k++) {
	    j = pr[i];
	    if (fabs(AAv[k]) > scaletmpc[p[k]]) scaletmpc[p[k]] = fabs(AAv[k]);
	 }
      }
      else if (nrm == 1) {
         for (k=pr[i]; k<pr[i+1]; k++)
            scaletmpc[p[k]] += fabs(AAv[k]);
      }
      else {
         for (k=pr[i]; k<pr[i+1]; k++)
            scaletmpc[p[k]] += fabs(AAv[k])*fabs(AAv[k]);
      }
   }
   if (nrm == 2) {
      for (i=0; i<n; i++)
	 scaletmpc[i] = sqrt(scaletmpc[i]);
   }
/*---------------------------------------
|   invert
|--------------------------------------*/
   for (i=0; i<n; i++) {
      if (scaletmpc[i] == 0.0)
	{
	 *ierr=i+1;
	 return i+1;
	}
      else 
	 scaletmpc[i] = 1.0 / scaletmpc[i];
   }
/*---------------------------------------
|   C = A * D
|--------------------------------------*/
   for (i=0; i<n; i++) {
    
      for (k=pr[i]; k<pr[i+1]; k++)
	AAv[k]=AAv[k]*scaletmpc[p[k]];
	
   }
   *ierr=0;
   return 0;
}
/*---------------end of coscalC-----------------------------------------
----------------------------------------------------------------------*/






void parm_param(string datafile,KN<long> param_int,KN<double> param_double)
	{
	char buf[BUFLEN];
	int num,in_val;
	double val;
	  FILE *fp;
	  char * filename=new char[datafile.length()+1]; 
          strcpy(filename,datafile.c_str()); 
	  int i;
	  for(i=0;i<16;i++) param_int[i]=-1; for(i=0;i<9;i++) param_double[i]=-1.0;
 	 /* read parameters for preconditioner and iteration from file  'filename' */
  	/*  ---- start modification by MS   */
 	 if( (fp = fopen(filename, "r")) == NULL ){
          fprintf(stderr, "Cannot open file inputs\n");
    	  exit(1);
    	}
  	num = 0;

  while(fgets(buf, BUFLEN, fp) != NULL) {
	if(num<=15) {sscanf(buf, "%d", &in_val); param_int[num]=in_val;}
	else 
	{sscanf(buf, "%lf", &val); param_double[num]=val;}
	num++;
	}
  fclose(fp);
}




class HipsSolver :   public MatriceMorse<double>::VirtualSolver   {
	double eps;
	mutable double  epsr;
	double tgv;	
	double tol_pivot_sym,tol_pivot; //Add 31 oct 2005
	string data_option;
	MPI_Comm  comm;
	mutable INTS id, idnbr, i, j;
	mutable INTS *unknownlist;
	mutable double *x;
	mutable INTS   ln;
	mutable INTS ierr;
	mutable INTS n,nnz;
	mutable double * a;
	mutable INTS *ia, *ja;
	mutable int *pp;
	int loc_size,pbegin, pend;
	INTS domsize, nproc,proc_id;
	mutable int sym_pattern, sym_matrix;
	KN<long> param_int;
	KN<double> param_double;
	mutable double *scaletmpr, *scaletmpc;
	
	mutable int *mapptr,*maptmp,*iwork,*riord,*iwork1,scale;
	mutable int *pr, *p;
	mutable double * AAv;
	
 public:
  HipsSolver(const MatriceMorse<double> &AA,string datafile, const KN<long> &param_int1, const KN<double> &param_double1,  MPI_Comm  * mpicommw  ) : data_option(datafile) 
  { 
	int argc,sym,symm;
	argc=1; 

	

	param_int=param_int1;
	param_double=param_double1;	
  	
	
        if(mpicommw==0){
	comm=MPI_COMM_WORLD;
	}
	else
	comm= *mpicommw;

	  
      
  	MPI_Comm_rank(comm, &proc_id);
  	MPI_Comm_size(comm, &nproc);
	
	id = 0; /** id of the linear system **/
	
	 
  /***************************************/
  /* Initialize HIPS for one problem     */
  /***************************************/
	idnbr = 1; /* total */
	ierr = HIPS_Initialize(idnbr);
	HIPS_ExitOnError(ierr);
	id = 0; /** id of the linear system **/

	int ic;
	
	if((!data_option.empty())&&((param_int==NULL)||(param_double==NULL))) 
		parm_param(datafile,param_int,param_double);
	else{
		if((param_double.N()>0)||(param_int.N()>0)){
		if(proc_id==0)      cout << "    WE SET PARAMETER FROM VECTORS   " << endl;}
		else if(proc_id==0) cout << "    WE SET DEFAULT PARAMETERS       " << endl;
	}
	
	 if(param_int.N()>0) {
	if((param_int[0]==HIPS_ITERATIVE)||(param_int[0]==HIPS_HYBRID)) // Using strategy. Input is ok
	{ HIPS_SetDefaultOptions(id,param_int[0]); ic=param_int[0];}
	else {	HIPS_SetDefaultOptions(id,HIPS_ITERATIVE );ic=HIPS_ITERATIVE; } // Strategy is not they existing one

	}else {HIPS_SetDefaultOptions(id,HIPS_ITERATIVE);ic=HIPS_ITERATIVE;} //Default strategy is ITERATIVE
	
	if(param_double.N()>0) {if(param_double[0]>0) HIPS_SetOptionREAL(id,HIPS_PREC,param_double[0]);
	else HIPS_SetOptionREAL(id,HIPS_PREC,1e-09);} else  HIPS_SetOptionREAL(id,HIPS_PREC,1e-09);	
	
	if((param_int.N()>1)&&(ic==HIPS_ITERATIVE)) {if((param_int[1]==1)||(param_int[1]==0)) HIPS_SetOptionINT(id,HIPS_KRYLOV_METHOD,param_int[1]); 
	else HIPS_SetOptionINT(id,HIPS_KRYLOV_METHOD,0);}else {if(ic==HIPS_ITERATIVE) HIPS_SetOptionINT(id,HIPS_KRYLOV_METHOD,0);}
	
	if(param_int.N()>2) { if(param_int[2]>0) HIPS_SetOptionINT(id,HIPS_ITMAX,param_int[2]); else  HIPS_SetOptionINT(id,HIPS_ITMAX,1000);}
	else  HIPS_SetOptionINT(id,HIPS_ITMAX,1000);

	if(param_int.N()>3) { if(param_int[3]>0) HIPS_SetOptionINT(id,HIPS_KRYLOV_RESTART,param_int[3]); 
	else  HIPS_SetOptionINT(id,HIPS_KRYLOV_RESTART,40);} else HIPS_SetOptionINT(id,HIPS_KRYLOV_RESTART,40);

	if(param_int.N()>4) {if(param_int[4]>0) {HIPS_SetOptionINT(id,HIPS_SYMMETRIC,param_int[4]); symm=param_int[4];}}
	else {HIPS_SetOptionINT(id,HIPS_SYMMETRIC,1);symm=1;}

	if(param_int.N()>5) {if((param_int[5]==0)||(param_int[5]==1)) {HIPS_SetOptionINT(id,HIPS_GRAPH_SYM,param_int[5]);sym=param_int[5];}
	else {HIPS_SetOptionINT(id,HIPS_GRAPH_SYM,1);sym=1;}} else  {HIPS_SetOptionINT(id,HIPS_GRAPH_SYM,1);sym=1;}	

	if(param_int.N()>6) { if(param_int[6]>0) HIPS_SetOptionINT(id,HIPS_PARTITION_TYPE,param_int[6]); 
	else  HIPS_SetOptionINT(id,HIPS_PARTITION_TYPE,0);}else  HIPS_SetOptionINT(id,HIPS_PARTITION_TYPE,0);			
		
	if(param_int.N()>7) {if(param_int[7]>0) HIPS_SetOptionINT(id,HIPS_LOCALLY,param_int[7]); 
	else  HIPS_SetOptionINT(id,HIPS_LOCALLY,2);}else HIPS_SetOptionINT(id,HIPS_LOCALLY,2);

	if(param_int.N()>8) {if(param_int[8]>0) HIPS_SetOptionINT(id,HIPS_FORTRAN_NUMBERING,param_int[8]);
	else HIPS_SetOptionINT(id,HIPS_FORTRAN_NUMBERING,0);} else  HIPS_SetOptionINT(id,HIPS_FORTRAN_NUMBERING,0);	
	
	if(param_int.N()>9) {if(param_int[9]>0) HIPS_SetOptionINT(id,HIPS_SCALE,param_int[9]);
	else { HIPS_SetOptionINT(id,HIPS_SCALE,1); scale=1;}} else  {HIPS_SetOptionINT(id,HIPS_SCALE,1); scale=1;}

	if(param_int.N()>10) {if(param_int[10]>0) HIPS_SetOptionINT(id,HIPS_REORDER,param_int[10]);
	else HIPS_SetOptionINT(id,HIPS_REORDER,1);} else  HIPS_SetOptionINT(id,HIPS_REORDER,1);	

	if(param_int.N()>11) {if(param_int[11]>0) HIPS_SetOptionINT(id,HIPS_DOF,param_int[11]);
	else HIPS_SetOptionINT(id,HIPS_DOF,1);} else  HIPS_SetOptionINT(id,HIPS_DOF,1);

	if(param_int.N()>12) {if(param_int[12]>0) HIPS_SetOptionINT(id,HIPS_SCALENBR,param_int[12]);
	else HIPS_SetOptionINT(id,HIPS_SCALENBR,2);} else  HIPS_SetOptionINT(id,HIPS_SCALENBR,2);
		
	if(param_int.N()>13) {if(param_int[13]>0) HIPS_SetOptionINT(id,HIPS_VERBOSE,param_int[13]);
	else HIPS_SetOptionINT(id,HIPS_VERBOSE,5);} else  HIPS_SetOptionINT(id,HIPS_VERBOSE,5);		

	if(param_int.N()>14) {if(param_int[14]>0) HIPS_SetOptionINT(id,HIPS_DOMSIZE,param_int[14]);
	else HIPS_SetOptionINT(id,HIPS_DOMSIZE,2);} else  HIPS_SetOptionINT(id,HIPS_DOMSIZE,2);	

	// if(param_int.N()>15) {if(param_int[15]>0) HIPS_SetOptionINT(id,HIPS_SCHUR_METHOD,param_int[15]);
        //else HIPS_SetOptionINT(id,HIPS_SCHUR_METHOD,2);} else  HIPS_SetOptionINT(id,HIPS_SCHUR_METHOD,2);

      //   if(param_int.N()>16) {if(param_int[16]>0) HIPS_SetOptionINT(id,HIPS_ITMAX_SCHUR,param_int[16]);
       // else HIPS_SetOptionINT(id,HIPS_ITMAX_SCHUR,2);} else  HIPS_SetOptionINT(id,HIPS_ITMAX_SCHUR,2);


			
	
	if(param_double.N()>1){if(param_double[1]>0.0) HIPS_SetOptionREAL(id,HIPS_DROPTOL0,param_double[1]);
	else  HIPS_SetOptionREAL(id,HIPS_DROPTOL0,0.005);} else   HIPS_SetOptionREAL(id,HIPS_DROPTOL0,0.005);

	if(param_double.N()>2) {if(param_double[2]>0.0) HIPS_SetOptionREAL(id,HIPS_DROPTOL1,param_double[2]);
	else  HIPS_SetOptionREAL(id,HIPS_DROPTOL1,0.005);} else  HIPS_SetOptionREAL(id,HIPS_DROPTOL1,0.005);

	if(param_double.N()>3){if(param_double[3]>0.0) HIPS_SetOptionREAL(id,HIPS_DROPTOLE,param_double[3]);
	else  HIPS_SetOptionREAL(id,HIPS_DROPTOLE,0.005);} HIPS_SetOptionREAL(id,HIPS_DROPTOLE,0.005);

	if(param_double.N()>4){if(param_double[4]>0.0) HIPS_SetOptionREAL(id,HIPS_AMALG,param_double[4]);
	else  HIPS_SetOptionREAL(id,HIPS_AMALG,0.005);} HIPS_SetOptionREAL(id,HIPS_AMALG,0.005);
	
	/*if(param_double.N()>5){if(param_double[5]>0.0) HIPS_SetOptionREAL(id,HIPS_DROPTOLSCHUR ,param_double[5]);
        else  HIPS_SetOptionREAL(id,HIPS_DROPTOLSCHUR ,0.005);} HIPS_SetOptionREAL(id,HIPS_DROPTOLSCHUR ,0.005);*/


	

    

     

        

	
		
		
	HIPS_SetCommunicator(id,comm);

	n=AA.n; nnz=AA.nbcoef;

	
	
	int ierr;
	
	pr= new int[n+1];
	p=  new int[nnz];
	AAv=new double[nnz];
	

	for(i=0;i<nnz;i++){
	AAv[i]=AA.a[i];
	p[i]=AA.cl[i];
	if(i<=n) pr[i]=AA.lg[i];
	}

	
	int job, tmp;
	 if(scale) {
	  	  job = 1;
	   	  tmp = 2; /*-- compute 2-norm */
		 scaletmpr=new double[n]; 		 scaletmpc=new double[n]; 
	   
	    roscal(n,job,tmp,AAv,p,pr,scaletmpr,&ierr);
	    if (ierr) fprintf(stderr, "Error: in roscal, ierr = %d\n", ierr);
	/*------- scale the RHS according to row scaling coefficients */
		
	    coscal(n,job,tmp,AAv,p,pr,scaletmpc,&ierr);
	    if (ierr) fprintf(stderr, "Error: in coscal, ierr = %d\n", ierr);
	  
	  } /*--- end of branch on scaling */


	int wgtflag=0, numflag=0, volume;
	
	riord=(int *)malloc(sizeof(int)*n);
	if(riord==NULL) {if(nproc==0)printf("%s","Memory allocation failed in partition stage \n"); exit(1);}
	int option[5];	option[0]=0;
	if(nproc>1){
		METIS_PartGraphKway(&n, AA.lg, AA.cl, NULL, NULL, &wgtflag, &numflag,&nproc, option, &volume, riord);
	}
	else if(nproc==1){
		for (i=0; i<n; i++) riord[i]=0;
	}
	iwork=(int *)malloc(sizeof(int)*n);
	maptmp=(int *)malloc(sizeof(int)*n);
	mapptr=(int *)malloc(sizeof(int)*(nproc+1));
	iwork1=(int *)malloc(sizeof(int)*(nproc+1));
	
       if(nproc==1){iwork[0]=0;iwork1[0]=0;}
	
	for(i=0; i<nproc; i++){
		iwork[i]=0;iwork1[i]=0;
	}

	for(j=0; j<n; j++){
		iwork[riord[j]]++;
		iwork1[riord[j]]++;
	}
	numflag=0;
	for(i=0; i<nproc; i++){
		mapptr[i]=numflag;
		numflag+=iwork[i];
	}
	
	mapptr[nproc]=numflag;
	
	for (i=0; i<nproc; i++){
		iwork[i]=mapptr[i];
	}
	if(nproc==0) iwork[0]=mapptr[0];
	
	for(i=0; i<n; i++){
		maptmp[iwork[riord[i]]]=i;
		iwork[riord[i]]++;
	}
	int nnzz;
	nnzz=0;
	for(i=0;i<n;i++) if(riord[i]==proc_id){ nnzz=(AA.lg[i+1]-AA.lg[i])+nnzz;}
	ierr = HIPS_GraphBegin(id, n, nnzz);
	HIPS_ExitOnError(ierr);
	
       

	for(i=0;i<n;i++)
	{
		if(riord[i]==proc_id){
		for(j=AA.lg[i];j<AA.lg[i+1];j++)
		{
			ierr = HIPS_GraphEdge(id, i, AA.cl[j]);
			HIPS_ExitOnError(ierr);	
		}
		}
		
	}
	
	ierr = HIPS_GraphEnd(id);
  	HIPS_ExitOnError(ierr);
	if(proc_id==0)
	{
		ierr = HIPS_SetPartition(id, nproc, mapptr, maptmp);
	        HIPS_ExitOnError(ierr);
	}
	 ierr = HIPS_AssemblyBegin(id, nnzz, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_FOOL,symm);
 	 HIPS_ExitOnError(ierr);
 
	 for(i=0;i<n;i++)
	{
		if(riord[i]==proc_id){
		for(j=AA.lg[i];j<AA.lg[i+1];j++)
		{					
			ierr = HIPS_AssemblySetValue(id, i, AA.cl[j], AAv[j]);
			HIPS_ExitOnError(ierr);
		}
		}
		
	}
    	ierr = HIPS_AssemblyEnd(id);
  	HIPS_ExitOnError(ierr);
	if(pr!=NULL)	free(pr); if(p!=NULL)	free(p); if(AAv!=NULL) free(AAv);
  }
	void Solver(const MatriceMorse<double> &AA,KN_<double> &x,const KN_<double> &b) const  {
	/***************************************************/
  	/*                                                 */
  	/*          ENTER THE RIGHT-HAND-SIDE              */
  	/*                                                 */
  	/***************************************************/
	int i,nloc;
	nloc=0; 
	int nnsize;
	MPI_Comm_size(comm,&nnsize);
	

 	COEF * rhsloc = (COEF *)malloc(sizeof(COEF)*iwork1[proc_id]);
 	COEF * xx = (COEF *)malloc(sizeof(COEF)*iwork1[proc_id]);
	int * unknownlist = (INTS *)malloc(sizeof(INTS)*iwork1[proc_id]);
	COEF * xz = (COEF *)malloc(sizeof(COEF)*n);
  	
	nloc=0;
	if(scale){
	for(i=0;i<n;i++)
	{
		if(riord[i]==proc_id){
		rhsloc[nloc]=b[i]*scaletmpr[i]; unknownlist[nloc]=i; nloc++; 
		}
	}
	}
	for(i=0;i<iwork1[proc_id];i++) xx[i]=0.0;
	 ierr = HIPS_SetRHS(id, nloc, unknownlist, rhsloc, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_FOOL);
  	HIPS_ExitOnError(ierr);
	
	 /****************************************************/
  	/* Get the local solution                           */
  	/****************************************************/ 
	
  	ierr = HIPS_GetSolution(id, nloc, unknownlist, xx, HIPS_ASSEMBLY_FOOL);
	
  	HIPS_ExitOnError(ierr);
	
	int *perm, *invp;
	perm=(int *)malloc(sizeof(int)*n); invp=(int *)malloc(sizeof(int)*n);



	MPI_Gatherv(xx,iwork1[proc_id], MPI_DOUBLE, xz,iwork1,mapptr,MPI_DOUBLE,0,comm);
        MPI_Gatherv(unknownlist,iwork1[proc_id], MPI_INT, perm,iwork1,mapptr,MPI_INT,0,comm);
        MPI_Bcast(xz,n,MPI_DOUBLE,0, comm); MPI_Bcast(perm,n,MPI_INT,0, comm);
	
	for(i=0;i<n;i++) invp[perm[i]]=i;
	if(scale){for(i=0;i<n;i++) {x[i]=xz[invp[i]]; x[i]=x[i]*scaletmpc[i];}}


	/**************************************************/
	/* Free HIPS internal structure for problem "id"  */
	/*************************************************/

	
             
		free(xz); free(perm); free(invp);free(rhsloc); free(unknownlist);  free(xx);
		
		HIPS_SetOptionINT(id,HIPS_DISABLE_PRECOND,0);HIPS_ExitOnError(ierr);
 		ierr = HIPS_Clean(id);
  		HIPS_ExitOnError(ierr);	
	
	}
	~HipsSolver()
	{
		if(verbosity>3)
       		cout << "~Hips_Solver S:" << endl;
		free(iwork1); free(mapptr); 
		free(iwork); free(maptmp);
		
		HIPS_SetOptionINT(id,HIPS_DISABLE_PRECOND,0);HIPS_ExitOnError(ierr);
 		ierr = HIPS_Clean(id);
  		HIPS_ExitOnError(ierr);	
		
	}	
   void addMatMul(const KN_<double> & x, KN_<double> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<double> &) (*this) * x; 
  }
};

//BuildSolverIUMFPack(DCL_ARG_SPARSE_SOLVER(double,A))
inline MatriceMorse<double>::VirtualSolver *
BuildSolverHipsSolvermpi(DCL_ARG_SPARSE_SOLVER(double,A))
{
    if(verbosity>9)
    cout << " BuildSolverSuperLU<double>" << endl;
    return new HipsSolver(*A,ds.data_filename, ds.lparams, ds.dparams,(MPI_Comm *)ds.commworld);
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
    if(verbosity>1)
	cout << " SetDefault sparse to default" << endl;
    DefSparseSolver<double>::solver =SparseMatSolver_R;
   //DefSparseSolver<Complex>::solver =SparseMatSolver_C;
    TypeSolveMat::defaultvalue =TypeSolveMat::SparseSolver;
	return 1;
}

bool SetHipsSolver()
{
    if(verbosity>1)
	cout << " SetDefault sparse solver to Hips" << endl;
    DefSparseSolver<double>::solver  =BuildSolverHipsSolvermpi;
    //DefSparseSolver<Complex>::solver =BuildSolverHipsSolvermpi;    
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
	return 1;
}
Init init;
Init::Init()
{   
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  //SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  
  if(verbosity>1)
    cout << "\n Add: Hips,  defaultsolver defaultsolverHips" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver;
  DefSparseSolver<double>::solver  =BuildSolverHipsSolvermpi;
//  DefSparseSolver<Complex>::solver =BuildSolverHipsSolver;
  if(! Global.Find("defaultsolver").NotNull() )
  Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
  Global.Add("defaulttoHips","(",new OneOperator0<bool>(SetHipsSolver));
}




