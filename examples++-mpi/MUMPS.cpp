// SUMMARY  :   simpler MUMPS interface
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : P. Jolivet 
// E-MAIL   : Pierre Jolivet <pierre.jolivet@ljll.math.upmc.fr>
//
//ff-c++-LIBRARY-dep:  mumps parmetis ptscotch scotch  scalapack blas  mpifc  fc mpi  pthread 
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

 */

#include <mpi.h>
#include "rgraph.hpp"
#include "AFunction.hpp"

// FFCS - 23/4/13 - instanciate some global symbols which are not found by default in MS MPI Fortran libraries
#ifdef _WIN32
__declspec(dllexport) int toto;
MPI_Fint* _imp__MPI_F_STATUS_IGNORE;
MPI_Fint* _imp__MPI_F_STATUSES_IGNORE;
#endif

#include "MatriceCreuse.hpp"

#include "dmatrix.hpp"
#include <dmumps_c.h>
#include <zmumps_c.h>
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
//#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
//#define INFOG(I) infog[(I)-1] /* macro s.t. indices match documentation */
//#define INFO(I) info[(I)-1]   /* macro s.t. indices match documentation */

template<typename RR> struct MUMPS_STRUC_TRAIT {typedef void MUMPS;  typedef void R; };
template<> struct MUMPS_STRUC_TRAIT<double>  {typedef DMUMPS_STRUC_C MUMPS; typedef double R;};
template<> struct MUMPS_STRUC_TRAIT<Complex>  {typedef ZMUMPS_STRUC_C MUMPS; typedef ZMUMPS_COMPLEX R;};
void mumps_c(DMUMPS_STRUC_C *id) { dmumps_c(id);}  
void mumps_c(ZMUMPS_STRUC_C *id) { zmumps_c(id);}  

template<class T> struct MPI_TYPE {static const MPI_Datatype  TYPE(){return MPI_BYTE;}};;
template<> struct MPI_TYPE<long>      {static const MPI_Datatype  TYPE(){return MPI_LONG;}};
template<> struct MPI_TYPE<int>      {static const MPI_Datatype TYPE(){return MPI_INT;}};
template<> struct MPI_TYPE<double>    {static const MPI_Datatype TYPE(){return MPI_DOUBLE;}};
template<> struct MPI_TYPE<char>    {static const MPI_Datatype TYPE(){return MPI_BYTE;}};
template<> struct MPI_TYPE<Complex>    {static const MPI_Datatype TYPE(){return MPI_DOUBLE_COMPLEX;}};


static std::string analysis[] = {"AMD", "", "AMF", "SCOTCH", "PORD", "METIS", "QAMD", "automatic sequential", "automatic parallel", "PT-SCOTCH", "ParMetis"};

template<class R>
class SolverMumps : public MatriceMorse<R>::VirtualSolver {



    private:
         mutable typename MUMPS_STRUC_TRAIT<R>::MUMPS * _id; 
  //mutable DMUMPS_STRUC_C* _id;
        mutable unsigned char   _strategy;
    bool distributed;
    MPI_Comm comm;
    int mpirank;
    KN<double> * rinfog;
    KN<long> * infog;

    public:
    
    int & ICNTL(int i) const  { return _id->icntl[i-1];}
    double & CNTL(int i) const  { return _id->cntl[i-1];}
    int & INFO(int i) const { return _id->info[i-1];}
    double & RINFO(int i) const { return _id->rinfo[i-1];}
    int & INFOG(int i) const { return _id->infog[i-1];}
    double & RINFOG(int i) const { return _id->rinfog[i-1];}
    void SetVerb(int i=verbosity) const {
        if( verbosity < 5)
        {
            ICNTL(1)=6;
            ICNTL(3)=0;
            ICNTL(4)=0;
        }
        else
        {
            ICNTL(1)=6;
            ICNTL(2)=0;
            ICNTL(3)=6;
            ICNTL(4)=0;
            if(verbosity < 10) ICNTL(4)=1;
            else if(verbosity < 15) ICNTL(4)=1;
            else if(verbosity < 20) ICNTL(4)=2;
            else if(verbosity < 25) ICNTL(4)=3;
            else ICNTL(4)=4;
        }
        //cout << ICNTL(1) << " " << ICNTL(2) << " "<< ICNTL(3) << " "<< ICNTL(4) << endl; 
    }
  
    
        typedef typename  MUMPS_STRUC_TRAIT<R>::R MR; 
  SolverMumps(const MatriceMorse<R> &A, KN<long> &param_int, KN<double> &param_R, MPI_Comm* pcomm,int strategy=3,int matrank=0,KN<double> *rinfogg=0  ,KN<long> *infogg=0)
    : comm( pcomm ? *pcomm :MPI_COMM_WORLD ),
      distributed(matrank<0),
      rinfog(rinfogg),infog(infogg)
  {
    
    MPI_Comm_rank(comm, &mpirank);
    int master = mpirank==matrank;  
    _id = new typename MUMPS_STRUC_TRAIT<R>::MUMPS ;
    _id->job = JOB_INIT;
    _id->par = 1;
    
    _id->comm_fortran = MPI_Comm_c2f(comm);
    _id->sym = A.symetrique;
    _strategy = strategy;
    mumps_c(_id);
    int* I = 0;
    int* J = 0;
    R * C = 0;
    long nnz=0;
    if( distributed || (mpirank == matrank) )
        {
            	
	if(_id->sym == 0) 
	  {
	    nnz = A.nbcoef;
	    I = new int[A.nbcoef];
	    CSR2COO<'C', 'U'>(A.n, A.lg, I);
	    C =  A.a;
	    J = A.cl;
	    for(unsigned int i = 0; i < A.nbcoef; ++i)
	      ++J[i];
	  }
	else 
	  {
	    if(A.symetrique) 
	      {
		nnz = A.nbcoef;
		I = new int[A.nbcoef];
		J = new int[A.nbcoef];
		C = new R[A.nbcoef];
		for(unsigned int i = 0; i < A.n; ++i)
		  C[i] = A.a[A.lg[i + 1] - 1];
		std::generate(I, I + A.n, step(0, 1));
		CSR2COO<'C', 'L'>(A.n, A.lg, I + A.n);
		std::copy(I, I + A.n, J);
		for(unsigned int i = 1; i < A.n; ++i) {
		  for(unsigned int j = A.lg[i]; j < A.lg[i + 1] - 1; ++j) {
		    J[A.n + j - i] = A.cl[j] + 1;
		    C[A.n + j - i] = A.a[j];
		  }
		}
	      }
	    else 
	      {
		nnz = A.n + (A.nbcoef - A.n) / 2;
		I = new int[A.n + (A.nbcoef - A.n) / 2];
		J = new int[A.n + (A.nbcoef - A.n) / 2];
		C = new R[A.n + (A.nbcoef - A.n) / 2];
		trimCSR<false, 'F',R>(A.n, I + A.n, A.lg, J + A.n, A.cl, C + A.n, A.a);
		for(unsigned int i = 0; i < A.n - 1; ++i)
		  C[i] = A.a[A.lg[i + 1] - (I[i + 1 + A.n] - I[i + A.n]) - 1];
		C[A.n - 1] = A.a[A.nbcoef - 1];
		std::generate(I, I + A.n, step(0, 1));
		CSR2COO<'F', 'U'>(A.n - 1, I + A.n, I + A.n);
		std::copy(I, I + A.n, J);
	      }
	  }
        _id->n = A.n;
        if(!distributed)
	  {
	    _id->nz=nnz;
	    _id->a =reinterpret_cast<MR *>( C);
	    _id->irn = I;
	    _id->jcn = J;
	  }
        else
	  {
	    
            _id->nz_loc=nnz;
            _id->a_loc =reinterpret_cast<MR *>( C);
            _id->irn_loc = I;
            _id->jcn_loc = J;
	  }
	
	}
    else 
      { // no matrix ...
	_id->nz=0;
	_id->a =0; 
	_id->irn = 0;;
	_id->jcn = 0;
      }
    _id->nrhs = 1;
    ICNTL(1) = 0;
    ICNTL(2) = 0;
    ICNTL(3) = verbosity > 1 ? 6 : 0;
    ICNTL(4) = 0; // verbose level
    ICNTL(5) = 0;                                                          // assembled format
    if(_strategy > 0 && _strategy < 9 && _strategy != 2) 
      {
	ICNTL(28) = 1;             // 1: sequential analysis
	ICNTL(7)  = _strategy - 1; //     0: AMD
      }     
    //     1:
    //     2: AMF
    //     3: SCOTCH
    //     4: PORD
    //     5: METIS
    //     6: QAMD
    //     7: automatic
    else
      {
	ICNTL(28) = 1;
	ICNTL(7)  = 7;
      }
    if(_strategy > 8 && _strategy < 12) 
      {
	ICNTL(28) = 2;              // 2: parallel analysis
	ICNTL(29) = _strategy - 9;  //     0: automatic
      }                                   //     1: PT-STOCH
    //     2: ParMetis
    ICNTL(9)  = 1;
    ICNTL(11) = 0;                 // verbose level
    ICNTL(18) = distributed ? 3: 0;        // centralized matrix input if !distributed
    ICNTL(20) = 0;                 // dense RHS
    ICNTL(14) = 30;                // percentage increase in the estimated working space
    _id->job = 4;
    mumps_c(_id);
    if(INFOG(1) != 0)
      std::cout << "BUG MUMPS, INFOG(1) = " << INFOG(1) << " distributed: " << distributed << " master " << matrank << std::endl;
    if(I) {
      if(_id->sym == 0) {
	for(unsigned int i = 0; i < A.nbcoef; ++i)
	  --J[i];
      }
      else {
	delete [] C;
	delete [] J;
      }
      delete [] I;
    }
      if( rinfog)
      {
          // copy rinfog
          if(rinfog->N() <40) rinfog->resize(40);
          for(int i=0; i<40;++i)
              (*rinfog)[i]=RINFOG(i+1);
      }
      if( infog)
      {
          // copy ginfo
          if(infog->N() <40) infog->resize(40);
          for(int i=0; i<40;++i)
              (*infog)[i]=INFOG(i+1);
      }
     
  };
  
  void Solver(const MatriceMorse<R> &A, KN_<R> &x, const KN_<R> &b) const
    {
     ICNTL(20) = 0; // dense RHS
     ICNTL(21) = 0; // centralized dense solution 
     if(distributed)
     {
       MPI_Reduce( (void *) (R*) b,(void *) (R*) x  , x.N() , MPI_TYPE<R>::TYPE(),MPI_SUM,0,comm);
     }
     else if(mpirank==0)  x = b;
    ICNTL(3) = verbosity > 1 ? 6 : 0;
    
    _id->rhs = reinterpret_cast<MR*>((R*) x);
    _id->job = 3;
    mumps_c(_id);
   if(distributed)
        {
          MPI_Bcast(reinterpret_cast<void*> ( (R*) x), x.N(), MPI_TYPE<R>::TYPE(), 0,comm);
        }
        if( rinfog)
        {
            // copy rinfog
            if(rinfog->N() <40) rinfog->resize(40);
            for(int i=0; i<40;++i)
                (*rinfog)[i]=RINFOG(i+1);
        }
        if( infog)
        {
            // copy ginfo
            if(infog->N() <40) infog->resize(40);
            for(int i=0; i<40;++i)
                (*infog)[i]=INFOG(i+1);
        }
        
    
  };
  
  ~SolverMumps() {
    _id->job = JOB_END;
    mumps_c(_id);
    if(_id)
      delete _id;
  };
};


template<class  R>
typename MatriceMorse<R>::VirtualSolver* buildSolver(DCL_ARG_SPARSE_SOLVER(R, A)) 
{// gestion de la star
 
  MPI_Comm cw = MPI_COMM_WORLD, * pcw= (MPI_Comm*) ds.commworld;
  if(!pcw) pcw = & cw;
  int mpirank ;
  MPI_Comm_rank(*pcw, &mpirank);
  int strategy = ds.strategy,matrank=0;
  if(Data_Sparse_Solver_version()>0 ) matrank=ds.master;
  if( !strategy && ds.lparams.N() > 0) strategy = ds.lparams[0]; 
  if(ds.lparams.N()>1) matrank = ds.lparams[1];// <0 => distri mat ..
  if( ! strategy) strategy=3;
  int mat = mpirank == matrank || matrank < 0; 
  if(A)
    return new SolverMumps<R>(*A, ds.lparams, ds.dparams, pcw,strategy,matrank,ds.rinfo,ds.info);
  else 
    ffassert(0); 
  return new SolverMumps<R>(*A, ds.lparams, ds.dparams, pcw,strategy,matrank,ds.rinfo,ds.info);
}

//  the 2 default sparse solver double and complex


DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; ;
DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
DefSparseSolverSym<double>::SparseMatSolver SparseMatSolverSym_R ; ;
DefSparseSolverSym<Complex>::SparseMatSolver SparseMatSolverSym_C;
// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetMUMPS()
{
    if(verbosity>1)
      cout << " SetDefault sparse solver to MUMPS (mpi) " << endl;
    DefSparseSolver<double>::solver  = buildSolver<double>;
    DefSparseSolver<Complex>::solver = buildSolver<Complex>;
    DefSparseSolverSym<double>::solver  = buildSolver<double>;
    DefSparseSolverSym<Complex>::solver = buildSolver<Complex>; 
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
    return  true;
}


void initMUMPS()
{
  if(verbosity && mpirank==0)  cout << "\n MUMPS (mpi) "<< endl;
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  SparseMatSolverSym_R= DefSparseSolverSym<double>::solver;
  SparseMatSolverSym_C= DefSparseSolverSym<Complex>::solver;
  
  if(verbosity>1)
    cout << "\n Add: MUMPS(mpi):  defaultsolver defaultsolverMUMPS_" << endl;
  DefSparseSolver<double>::solver  = buildSolver;
  DefSparseSolver<Complex>::solver = buildSolver;
  DefSparseSolverSym<double>::solver  = buildSolver;
  DefSparseSolverSym<Complex>::solver = buildSolver; 
  TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
  if(! Global.Find("defaulttoMUMPS").NotNull() )
    Global.Add("defaulttoMUMPS","(",new OneOperator0<bool>(SetMUMPS));  
}


LOADFUNC(initMUMPS);
