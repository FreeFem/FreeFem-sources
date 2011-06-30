// ORIG-DATE: 04/2009
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : INRIA FUTUR
// AUTHOR   : Guy Atenekeng
// E-MAIL   : Guy_Antoine_Atenekeng_Kahou@lri.fr
//
//ff-c++-LIBRARY-dep:  hips metis  blas  mpi
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

#include "rgraph.hpp"
#include "error.hpp"
#include "AFunction.hpp"
#include "MatriceCreuse_tpl.hpp"

#ifndef  MPI_SUCCESS
#define  MPI_SUCCESS
#endif

extern "C" {
#include "hips.h"  
#include "metis.h"
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFLEN 200
#define MCW MPI_COMM_WORLD


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
    mutable INTS id,  i, j;
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

   
private:
  static const int MaxIds=100;
  static int Ids[MaxIds];
  static int  GetId() {
    static bool Initialized=false;
    if(!Initialized)
      {
	Initialized=true;
	if(verbosity)
	cout << " Hips HIPS_Initialize " << MaxIds <<endl;
	INTS ierr = HIPS_Initialize(MaxIds);
	HIPS_ExitOnError(ierr);
	for(int i=0;i<MaxIds;++i)
	  {
	    Ids[i]=-1; // ununsed 
	  }
      }
    INTS id =-1;
    for(int i=0;i<MaxIds;++i)
      if( Ids[i] <0)  
	{
	  Ids[i]=i;
	  if(verbosity) cout << "   find HipsSoler :  id = " << i << "/" <<   MaxIds << endl;
	  return i;
	}
    cerr<< " All id of Hips are busy " << MaxIds << " try to store less matrix or change MaxIds (FH.) in " << endl; 
    ffassert(0); 
    return -1; 
  }
public:
  static void Def_iopt(long * param_intd)
  {
    param_intd[0]= HIPS_ITERATIVE  ; //  HIPS_STRA
    param_intd[1]= 0  ; //  HIPS_KRYLOV_METHOD
    param_intd[2]= 1000  ; //  HIPS_ITMAX
    param_intd[3]= 40  ; //  HIPS_KRYLOV_RESTART
    param_intd[4]= 1  ; //  HIPS_SYMMETRIC
    param_intd[5]= 1  ; //  HIPS_GRAPH_SYM
    param_intd[6]= 0  ; //  HIPS_PARTITION_TYPE
    param_intd[7]= 2  ; //  HIPS_LOCALLY
    param_intd[8]= 0  ; //  HIPS_FORTRAN_NUMBERING
    param_intd[9]= 1  ; //  HIPS_SCALE
    param_intd[10]= 1  ; //  HIPS_REORDER
    param_intd[11]= 1  ; //  HIPS_DOF
    param_intd[12]= 2  ; //  HIPS_SCALENBR
    param_intd[13]= verbosity  ; //  HIPS_VERBOSE
    param_intd[14]= 2  ; //  HIPS_DOMSIZE
    param_intd[15]= 2  ; //  HIPS_SCHUR_METHOD
    param_intd[16]= 2  ; //  HIPS_ITMAX_SCHUR
  }

  static void Def_dopt(double *d)
  {
    d[0]= 1e-09  ; //  HIPS_PREC
    d[1]= 0.005  ; //  HIPS_DROPTOL0
    d[2]= 0.005  ; //  HIPS_DROPTOL1
    d[3]= 0.005  ; //  HIPS_DROPTOLE
    d[4]= 0.005  ; //  HIPS_AMALG
    d[5]= 0.005  ; //  HIPS_DROPSCHUR
  }

public:
  
  
  HipsSolver(const MatriceMorse<double> &AA,string datafile, const KN<long> &param_int1, const KN<double> &param_double1,  MPI_Comm  * mpicommw  ) 
    : data_option(datafile) ,param_int(17), param_double(6),id(GetId())
  {
    
    if(mpicommw==0)
	comm=MPI_COMM_WORLD;
    else 
	comm= *mpicommw;
    
    MPI_Comm_rank(comm, &proc_id);
    MPI_Comm_size(comm, &nproc);
    if(proc_id==0  || verbosity>2)
	cout << "  Hips Comm " << proc_id<< " / "<< nproc << endl;
    
    Def_iopt(param_int);
    Def_dopt(param_double);

    static int dopt_wrapper[6] = {
      HIPS_PREC ,
      HIPS_DROPTOL0 ,
      HIPS_DROPTOL1 ,
      HIPS_DROPTOLE ,
      HIPS_AMALG ,
      HIPS_DROPSCHUR 
 };

    static int iopt_wrapper[17] = {
      -1,// DEf STRATEGIC not in H
      HIPS_KRYLOV_METHOD ,
      HIPS_ITMAX,
      HIPS_KRYLOV_RESTART,
      HIPS_SYMMETRIC,
      HIPS_GRAPH_SYM,
      HIPS_PARTITION_TYPE, 
      HIPS_LOCALLY,
      HIPS_FORTRAN_NUMBERING,
      HIPS_SCALE,
      HIPS_REORDER,
      HIPS_DOF,
      HIPS_SCALENBR, 
      HIPS_VERBOSE,
      HIPS_DOMSIZE,
      HIPS_SCHUR_METHOD, 
      HIPS_ITMAX_SCHUR
    }     ; 
    
    int ic,sym=AA.symetrique ,symm=AA.symetrique; 

 
    if(!data_option.empty())
      parm_param(datafile,param_int,param_double);
    else
      {
	for(int i=0;i< min(param_int.N(),param_int1.N()); ++i) 
	  if(param_int1[i]>-1)  param_int[i]=param_int1[i];	    
	for(int i=0;i< min(param_double.N(),param_double1.N()); ++i) 
	  if(param_double1[i]>-0.9999)  param_double[i]=param_double1[i];	    
      }
    // force param  value ... 
    param_int[0]=max(min(param_int[0],2L),0L);
    param_int[5]= sym;
    param_int[4]= symm;
    
    ic = param_int[0];
    scale=param_int[9];
    
    if(verbosity>3 && proc_id==0  )
      {
	cout << " Hips INT  opts " << param_int << endl;
	cout << " Hips REAL  opts " << param_double << endl;
      }

    HIPS_SetDefaultOptions(id, param_int[0] );
    
    for(int i=1;i<param_int.N();++i)     // the fist value is teh STRATEGIE KING not aoption
      HIPS_SetOptionINT(id,iopt_wrapper[i],param_int[i] );
    for(int i=0;i<param_double.N();++i)    
      HIPS_SetOptionREAL(id,dopt_wrapper[i],param_double[i] );
 
    if(!data_option.empty()) 
      parm_param(datafile,param_int,param_double);
    
    HIPS_SetCommunicator(id,comm);
    
    n=AA.n; nnz=AA.nbcoef;
    
    int ierr;
    /*
    pr= new int[n+1];
    p=  new int[nnz];
    AAv=new double[nnz];
    
    
    for(int i=0;i<nnz;i++)
      {
	AAv[i]=AA.a[i];
	p[i]=AA.cl[i];
	if(i<=n) pr[i]=AA.lg[i];
      }
    */
    
    int job, tmp;
    if(scale) {
      job = 1;
      tmp = 2; /*-- compute 2-norm */
      scaletmpr=new double[n];
      scaletmpc=new double[n]; 
      
      roscal(n,job,tmp,AA.a,AA.cl,AA.lg,scaletmpr,&ierr);
      if (ierr) fprintf(stderr, "Error: in roscal, ierr = %d\n", ierr);
      /*------- scale the RHS according to row scaling coefficients */
      
      coscal(n,job,tmp,AA.a,AA.cl,AA.lg,scaletmpc,&ierr);
      if (ierr) fprintf(stderr, "Error: in coscal, ierr = %d\n", ierr);
      
    } /*--- end of branch on scaling */
    
    
    int wgtflag=0, numflag=0, volume;
    
    riord= new int[n]; //(int *)malloc(sizeof(int)*n);
    if(riord==NULL) {
      if(nproc==0)
	printf("%s","Memory allocation failed in partition stage \n"); 
      exit(1);
    }
    int option[5];	option[0]=0;
    if(nproc>1){
      METIS_PartGraphKway(&n, AA.lg, AA.cl, NULL, NULL, &wgtflag, &numflag,&nproc, option, &volume, riord);
    }
    else if(nproc==1){
      for (int i=0; i<n; i++) 
	riord[i]=0;
    }

    iwork= new int[nproc+1];// (int *)malloc(sizeof(int)*n);
    maptmp= new int [n];//(int *)malloc(sizeof(int)*n);
    mapptr= new int [nproc+1];//(int *)malloc(sizeof(int)*(nproc+1));
    iwork1= new int[nproc+1];//(int *)malloc(sizeof(int)*(nproc+1));
    
    for(int i=0; i<=nproc; i++)
      iwork[i]=iwork1[i]=0;

    for(int j=0; j<n; j++)
      {
	iwork[riord[j]]++;
	iwork1[riord[j]]++;
      }
    numflag=0;
    for(int i=0; i<nproc; i++)
      {
	mapptr[i]=numflag;
	numflag+=iwork[i];
      }
    
    mapptr[nproc]=numflag;
    
    for (int i=0; i<nproc; i++){
      iwork[i]=mapptr[i];
    }
    if(nproc==0) 
      iwork[0]=mapptr[0];
    
    for(int i=0; i<n; i++){
      maptmp[iwork[riord[i]]]=i;
      iwork[riord[i]]++;
    }
    int nnzz;
    nnzz=0;
      for(int i=0;i<n;i++) 
	if(riord[i]==proc_id)
	  {
	    nnzz+=(AA.lg[i+1]-AA.lg[i]);
	  }
      ierr = HIPS_GraphBegin(id, n, nnzz);
      HIPS_ExitOnError(ierr);
      if(verbosity > 5)
	  cout << "   Hips : proc " << proc_id << " / nzz = " << nnzz << endl;
      
      
      for(int i=0;i<n;i++)
	{
	  if(riord[i]==proc_id){
	    for(int j=AA.lg[i];j<AA.lg[i+1];j++)
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
      //symm=1;
      if( nproc ==1)
	{
	  ierr = HIPS_MatrixGlobalCSR(id, n, AA.lg, AA.cl, AA.a, 0, HIPS_ASSEMBLY_OVW, sym_matrix);
	  HIPS_ExitOnError(ierr);
	  
	}
      else 
	{	  
	  ierr = HIPS_AssemblyBegin(id, nnzz, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_FOOL,symm);
	  HIPS_ExitOnError(ierr);
	    int kkk;
	  for(int i=0;i<n;i++)
	    {
	      if(riord[i]==proc_id){
		  for(int k=AA.lg[i];k<AA.lg[i+1];k++)
		    {			
			kkk++;
			if(verbosity >100) cout << "       " << proc_id << " a( " << i << ", " <<AA.cl[k] << ")= " << AA.a[k] << endl;
			ierr = HIPS_AssemblySetValue(id, i, AA.cl[k], AA.a[k]);
			HIPS_ExitOnError(ierr);
		    }
	      }
	      
	    }
	    ffassert(kkk);
	  ierr = HIPS_AssemblyEnd(id);
	  
	  HIPS_ExitOnError(ierr);
	}
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
    
    
    COEF * rhsloc = new COEF[iwork1[proc_id]] ;//(COEF *)malloc(sizeof(COEF)*iwork1[proc_id]);
    COEF * xx =  new COEF[iwork1[proc_id]] ;// (COEF *)malloc(sizeof(COEF)*iwork1[proc_id]);
    INTS * unknownlist =  new INTS [iwork1[proc_id]] ;//(INTS *)malloc(sizeof(INTS)*iwork1[proc_id]);
    COEF * xz = new COEF[n]; // (COEF *)malloc(sizeof(COEF)*n);
    
    nloc=0;
    // if(scale)
      {
	for(i=0;i<n;i++)
	  {
	    if(riord[i]==proc_id){
	      if(scale) rhsloc[nloc]=b[i]*scaletmpr[i]; 
	      unknownlist[nloc++]=i;
	    }
	  }
      }
    for(i=0;i<iwork1[proc_id];i++) 
      xx[i]=0.0;
    ierr = HIPS_SetRHS(id, nloc, unknownlist, rhsloc, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_FOOL);
    HIPS_ExitOnError(ierr);
    
    /****************************************************/
    /* Get the local solution                           */
    /****************************************************/ 
    
    ierr = HIPS_GetSolution(id, nloc, unknownlist, xx, HIPS_ASSEMBLY_FOOL);
    
    HIPS_ExitOnError(ierr);
    
    int *perm = new int[n], *invp= new int[n];
    
    MPI_Gatherv(xx,iwork1[proc_id], MPI_DOUBLE, xz,iwork1,mapptr,MPI_DOUBLE,0,comm);
    MPI_Gatherv(unknownlist,iwork1[proc_id], MPI_INT, perm,iwork1,mapptr,MPI_INT,0,comm);
    MPI_Bcast(xz,n,MPI_DOUBLE,0, comm);
    MPI_Bcast(perm,n,MPI_INT,0, comm);
    
    for(int i=0;i<n;i++) 
      invp[perm[i]]=i;
    if(scale)
      {
	for(int i=0;i<n;i++) 
	  {
	    x[i]=xz[invp[i]];
	    x[i]=x[i]*scaletmpc[i];
	  }
      }
    
    
    /**************************************************/
    /* Free HIPS internal structure for problem "id"  */
    /*************************************************/
      /*
      INTL nnzp;
     ierr= HIPS_GetInfoINT(id,HIPS_INFO_NNZ_PEAK,&nnzp);
    if(verbosity>1 && ierr==HIPS_SUCCESS ) 
	cout << "  Hips Peak nnz  = " << nnzp  << "   ( id = " << id << ")" << endl; ;
    */
    delete [] xz;
    delete [] perm;
    delete [] invp;
    delete [] rhsloc; 
    delete [] unknownlist;
    delete [] xx;
    
    
  }
							 
  ~HipsSolver()
  {
    assert(id>=0);
    if( (verbosity>3 && proc_id==0 ) ||(verbosity>9) )
	cout << "   ~Hips_Solver S:" << id << endl;
  //  HIPS_SetOptionINT(id,HIPS_DISABLE_PRECOND,0);
  //  HIPS_ExitOnError(ierr);
    ierr = HIPS_Clean(id);
    HIPS_ExitOnError(ierr);	

    delete [] iwork1;
    delete [] mapptr;

    delete [] iwork;
    delete [] maptmp;
    
    if(id>0 && id< MaxIds) Ids[id]=-2;
    id=-2; 
 	
  }
  
  
  void addMatMul(const KN_<double> & x, KN_<double> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<double> &) (*this) * x; 
  }
private:// no copy
    HipsSolver(const HipsSolver &);  
    HipsSolver & operator=(const HipsSolver &); 
};  // CLASS HipsSolver
  
int HipsSolver::Ids[HipsSolver::MaxIds];




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
bool HipsDefaults(KN<long>* piop,KN<double> * pdop)
{
  if(piop)
    {
      piop->resize(17);
      HipsSolver::Def_iopt(*piop);
    }
  if(pdop)
    {
      pdop->resize(6);
      HipsSolver::Def_dopt(*pdop);      
    }
  
  return true;
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
  Global.Add("HipsDefaults","(",new OneOperator2<bool,KN<long>*,KN<double> *>(HipsDefaults));
  
}
		      
		      
		      
		      
