// ORIG-DATE: 04/2009
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
#include "MatriceCreuse_tpl.hpp"
#include "mpi.h"

extern "C" {
  #include "hips.h"  
  #include "io.h"
  #include "metis.h"
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFLEN 200
#define MCW MPI_COMM_WORLD

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
	mutable MPI_Comm comm;
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
	
	mutable int *mapptr,*maptmp,*iwork,*riord,*iwork1;
	
 public:
  HipsSolver(const MatriceMorse<double> &AA,string datafile, const KN<long> &param_int1, const KN<double> &param_double1) : data_option(datafile) 
  { 
	int argc,sym,symm;
	argc=1; 

	

	param_int=param_int1;
	param_double=param_double1;	
  	
	comm=MPI_COMM_WORLD;
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

	//char * filename=new char[data_option.length()+1]; 
         //strcpy(filename,data_option.c_str()); 
	
	
	int ic;


	if(!data_option.empty()) 
	parm_param(datafile,param_int,param_double);
	

	 if(param_int.N()>0) {
	if((param_int[0]==HIPS_ITERATIVE)||(param_int[0]==HIPS_ILUT)||(param_int[0]==HIPS_HYBRID)||(param_int[0]==HIPS_BLOCK)|| (param_int[0]==HIPS_DIRECT)) 
	{ HIPS_SetDefaultOptions(id,param_int[0]); ic=param_int[0];}
	else {
		HIPS_SetDefaultOptions(id,HIPS_ITERATIVE );ic=HIPS_ITERATIVE; 
		if(proc_id==0) printf("%s", "WARNING:  THE METHOD INDEX DO NOT EXIST. WE SET HIPS_ITERATIVE AS METHOD \n ");  
	 }
	}else {HIPS_SetDefaultOptions(id,HIPS_ITERATIVE);ic=HIPS_ITERATIVE;}

	if((ic==HIPS_ITERATIVE)||(ic==HIPS_ILUT))  {
	
	    if(param_double.N()==0){HIPS_SetOptionREAL(id,HIPS_DROPTOL0,0.005); HIPS_SetOptionREAL(id,HIPS_DROPTOL1,0.005); 
	    HIPS_SetOptionREAL(id,HIPS_DROPTOLE,0.005); HIPS_SetOptionREAL(id,HIPS_AMALG,0.005);
	    //set default double parameter
	  }

	    if(param_double.N()>1){if(param_double[1]!=-1.0) HIPS_SetOptionREAL(id,HIPS_DROPTOL0,param_double[1]);
	    else  HIPS_SetOptionREAL(id,HIPS_DROPTOL0,0.005);} else   HIPS_SetOptionREAL(id,HIPS_DROPTOL0,0.005);

	    if(param_double.N()>2) {if(param_double[2]!=-1.0) HIPS_SetOptionREAL(id,HIPS_DROPTOL1,param_double[2]);
	    else  HIPS_SetOptionREAL(id,HIPS_DROPTOL1,0.005);} else  HIPS_SetOptionREAL(id,HIPS_DROPTOL1,0.005);

	    if(param_double.N()>4){if(param_double[4]!=-1.0) HIPS_SetOptionREAL(id,HIPS_DROPTOLE,param_double[4]);
	    else  HIPS_SetOptionREAL(id,HIPS_DROPTOLE,0.005);} HIPS_SetOptionREAL(id,HIPS_DROPTOLE,0.005);	

	    if(param_double.N()>5){if(param_double[5]!=-1.0) HIPS_SetOptionREAL(id,HIPS_AMALG,param_double[5]);
	    else  HIPS_SetOptionREAL(id,HIPS_AMALG,0.005);} HIPS_SetOptionREAL(id,HIPS_AMALG,0.005);	

	   
	}
		
		
	else 
	{
		if(param_double.N()>1){if(param_double[1]!=-1.0) HIPS_SetOptionREAL(id,HIPS_DROPTOL0,param_double[1]);
		else  HIPS_SetOptionREAL(id,HIPS_DROPTOL0,0.0);} else   HIPS_SetOptionREAL(id,HIPS_DROPTOL0,0.05);

		if(param_double.N()>2) {if(param_double[2]!=-1.0) HIPS_SetOptionREAL(id,HIPS_DROPTOL1,param_double[2]);
		else  HIPS_SetOptionREAL(id,HIPS_DROPTOL1,0.0);} else  HIPS_SetOptionREAL(id,HIPS_DROPTOL1,0.05);

		if(param_double.N()>4){if(param_double[4]!=-1.0) HIPS_SetOptionREAL(id,HIPS_DROPTOLE,param_double[4]);
		else  HIPS_SetOptionREAL(id,HIPS_DROPTOLE,0.0);} HIPS_SetOptionREAL(id,HIPS_DROPTOLE,0.05);

		if(param_double.N()>5){if(param_double[5]!=-1.0) HIPS_SetOptionREAL(id,HIPS_AMALG,param_double[5]);
		else  HIPS_SetOptionREAL(id,HIPS_AMALG,0.0);} HIPS_SetOptionREAL(id,HIPS_AMALG,0.05);	

		
	}	
		if(param_double.N()>0) {if(param_double[0]>0) HIPS_SetOptionREAL(id,HIPS_PREC,param_double[0]);
		else HIPS_SetOptionREAL(id,HIPS_PREC,1e-08);} else  HIPS_SetOptionREAL(id,HIPS_PREC,1e-08);		

		if(param_int.N()>9) {if((param_int[9]==1)||(param_int[9]==0)) {HIPS_SetOptionINT(id,HIPS_GRAPH_SYM,param_int[9]);sym=param_int[9];}
		else {HIPS_SetOptionINT(id,HIPS_GRAPH_SYM,0);sym=0;}} else  {HIPS_SetOptionINT(id,HIPS_GRAPH_SYM,0);sym=0;}
		
		if(param_int.N()>1) {if((param_int[1]==1)||(param_int[1]==0)) {HIPS_SetOptionINT(id,HIPS_SYMMETRIC,param_int[1]); symm=0;}
		else {HIPS_SetOptionINT(id,HIPS_SYMMETRIC,0);symm=0;}} else  {HIPS_SetOptionINT(id,HIPS_SYMMETRIC,0);symm=0;}


		if(param_int.N()>5) {if((param_int[5]>0)&&(param_int[1]==1)) HIPS_SetOptionINT(id,HIPS_KRYLOV_METHOD,1); 
		else if((param_int[5]>0)&&(param_int[1]==0))  HIPS_SetOptionINT(id,HIPS_KRYLOV_METHOD,0);}else HIPS_SetOptionINT(id,HIPS_KRYLOV_METHOD,0);
	
		
		if(param_int.N()>2) {if(param_int[2]>0) HIPS_SetOptionINT(id,HIPS_FORTRAN_NUMBERING,param_int[2]);
		else HIPS_SetOptionINT(id,HIPS_FORTRAN_NUMBERING,0);} else  HIPS_SetOptionINT(id,HIPS_FORTRAN_NUMBERING,0);


		if(param_int.N()>3) { if(param_int[3]>0) HIPS_SetOptionINT(id,HIPS_ITMAX,param_int[3]); else  HIPS_SetOptionINT(id,HIPS_ITMAX,500);}
		else  HIPS_SetOptionINT(id,HIPS_ITMAX,500);

		if(param_int.N()>4) { if(param_int[4]>0) HIPS_SetOptionINT(id,HIPS_KRYLOV_RESTART,param_int[4]); 
		else  HIPS_SetOptionINT(id,HIPS_KRYLOV_RESTART,30);} else HIPS_SetOptionINT(id,HIPS_KRYLOV_RESTART,30);
				
		if(param_int.N()>6) { if(param_int[6]>0) HIPS_SetOptionINT(id,HIPS_PARTITION_TYPE,param_int[6]); 
		else  HIPS_SetOptionINT(id,HIPS_PARTITION_TYPE,0);}else  HIPS_SetOptionINT(id,HIPS_PARTITION_TYPE,0);		

		
		if(param_int.N()>7) {if(param_int[7]>0) HIPS_SetOptionINT(id,HIPS_LOCALLY,param_int[7]); 
		else  HIPS_SetOptionINT(id,HIPS_LOCALLY,0);}else HIPS_SetOptionINT(id,HIPS_LOCALLY,0);

		if(param_int.N()>8) {if(param_int[8]>0) HIPS_SetOptionINT(id,HIPS_VERBOSE,param_int[8]);
		else HIPS_SetOptionINT(id,HIPS_VERBOSE,0);} else  HIPS_SetOptionINT(id,HIPS_VERBOSE,5);

		
		if(param_int.N()>10) {if(param_int[10]>0) HIPS_SetOptionINT(id,HIPS_REORDER,param_int[10]);
		else HIPS_SetOptionINT(id,HIPS_REORDER,0);} else  HIPS_SetOptionINT(id,HIPS_REORDER,0);
		
		if(param_int.N()>11) {if(param_int[11]>0) HIPS_SetOptionINT(id,HIPS_SCALENBR,param_int[11]);
		else HIPS_SetOptionINT(id,HIPS_SCALENBR,1);} else  HIPS_SetOptionINT(id,HIPS_SCALENBR,1);

		if(param_int.N()>12) {if(param_int[12]>0) HIPS_SetOptionINT(id,HIPS_DOF,param_int[12]);
		else HIPS_SetOptionINT(id,HIPS_DOF,1);} else  HIPS_SetOptionINT(id,HIPS_DOF,1);


	n=AA.n; nnz=AA.nbcoef;
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
			ierr = HIPS_AssemblySetValue(id, i, AA.cl[j], AA.a[j]);
			HIPS_ExitOnError(ierr);
		}
		}
		
	}
    	ierr = HIPS_AssemblyEnd(id);
  	HIPS_ExitOnError(ierr);
  }
	void Solver(const MatriceMorse<double> &AA,KN_<double> &x,const KN_<double> &b) const  {
	/***************************************************/
  	/*                                                 */
  	/*          ENTER THE RIGHT-HAND-SIDE              */
  	/*                                                 */
  	/***************************************************/
	int i,nloc;
	nloc=0; 

 	COEF * rhsloc = (COEF *)malloc(sizeof(COEF)*iwork1[proc_id]);
 	COEF * xx = (COEF *)malloc(sizeof(COEF)*iwork1[proc_id]);
	int * unknownlist = (INTS *)malloc(sizeof(INTS)*iwork1[proc_id]);
	COEF * xz = (COEF *)malloc(sizeof(COEF)*n);
  	x = (COEF *)malloc(sizeof(COEF)*n);
	nloc=0;

	for(i=0;i<n;i++)
	{
		if(riord[i]==proc_id){
		rhsloc[nloc]=b[i]; unknownlist[nloc]=i; nloc++; 
		}
		
	}
		

	 ierr = HIPS_SetRHS(id, nloc, unknownlist, rhsloc, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_FOOL);
  	HIPS_ExitOnError(ierr);
	
	 /****************************************************/
  	/* Get the local solution                           */
  	/****************************************************/ 
	
  	ierr = HIPS_GetSolution(id, nloc, unknownlist, xx, HIPS_ASSEMBLY_FOOL);
  	HIPS_ExitOnError(ierr);
	
	int *perm, *invp;
	perm=(int *)malloc(sizeof(int)*n); invp=(int *)malloc(sizeof(int)*n);



	MPI_Gatherv(xx,iwork1[proc_id], MPI_DOUBLE, xz,iwork1,mapptr,MPI_DOUBLE,0,MCW);
        MPI_Gatherv(unknownlist,iwork1[proc_id], MPI_INT, perm,iwork1,mapptr,MPI_INT,0,MCW);
        MPI_Bcast(xz,n,MPI_DOUBLE,0, MCW); MPI_Bcast(perm,n,MPI_INT,0, MCW);
	
	for(i=0;i<n;i++) invp[perm[i]]=i;
	for(i=0;i<n;i++) x[i]=xz[invp[i]];


	/**************************************************/
	/* Free HIPS internal structure for problem "id"  */
	/*************************************************/
	free(xz); free(perm); free(invp);free(rhsloc); free(unknownlist); free(iwork1); free(mapptr); free(xx);
	free(iwork); free(maptmp);
	HIPS_SetOptionINT(id,HIPS_DISABLE_PRECOND,0);HIPS_ExitOnError(ierr);
 	ierr = HIPS_Clean(id);
  	HIPS_ExitOnError(ierr);


	}
	~HipsSolver()
	{
		if(verbosity>3)
       		cout << "~Hips_Solver S:" << endl;
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
    return new HipsSolver(*A,ds.data_filename, ds.lparams, ds.dparams);
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




