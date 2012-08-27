// SUMMARY  :   add interface with partionning library scotch 
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : F. Hecht
// E-MAIL   : F. Hecht <hecht@ljll.math.upmc.fr>
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

 */

//   for automatic  compilation with ff-c++
//ff-c++-LIBRARY-dep:  mumps-seq blas  libseq  fc  pthread
//ff-c++-cpp-dep: 
//
// F. Hecht  december 2011
// ----------------------------
//  file to add MUMPS sequentiel interface for sparce linear solver with dynamic load.

#include  <iostream>
using namespace std;


#include "ff++.hpp"

#include "mpi.h"
#include "dmumps_c.h"
#include "zmumps_c.h"

const int  JOB_INIT = -1;
const int  JOB_END =-2;
const int JOB_ANA_FAC = 4;
const int JOB_SOLVE = 3;
const int USE_COMM_WORLD = -987654; 

template<typename RR> struct MUMPS_STRUC_TRAIT {typedef void MUMPS;  typedef void R; };
template<> struct MUMPS_STRUC_TRAIT<double>  {typedef DMUMPS_STRUC_C MUMPS; typedef double R;};
template<> struct MUMPS_STRUC_TRAIT<Complex>  {typedef ZMUMPS_STRUC_C MUMPS; typedef ZMUMPS_COMPLEX R;};
void mumps_c(DMUMPS_STRUC_C *id) { dmumps_c(id);}  
void mumps_c(ZMUMPS_STRUC_C *id) { zmumps_c(id);}  


template<typename  R>
class SolveMUMPS_seq :   public MatriceMorse<R>::VirtualSolver 
{
public:
  // typedef double R;
  double eps;
  mutable double  epsr;
  double tgv;
  typedef typename  MUMPS_STRUC_TRAIT<R>::R MR; 
  mutable typename MUMPS_STRUC_TRAIT<R>::MUMPS id; 

  int & ICNTL(int i) const  { return id.icntl[i-1];}
  R & CNTL(int i) const  { return id.cntl[i-1];}
  int & INFO(int i) const { return id.info[i-1];}
  R & RINFO(int i) const { return id.rinfo[i-1];}
  int & INFOG(int i) const { return id.infog[i-1];}
  R & RINFOG(int i) const { return id.rinfog[i-1];}
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
  void Check(const char * msg="mumps-seq") const 
  {
    if( INFO(1) !=0) 
      {
	cout << " Erreur Mumps number " << INFO(1) << endl;
	cout << " Fatal Erreur  " << msg << endl;
	Clean();
	id.job=JOB_END;
	mumps_c(&id); /* Terminate instance */
	int ierr = MPI_Finalize();
	ErrorExec(msg,INFO(1) ); 
      }
  }
  SolveMUMPS_seq(const MatriceMorse<R> &A,int strategy,double ttgv, double epsilon=1e-6,
	       double pivot=-1.,double pivot_sym=-1.  ) : 
    eps(epsilon),epsr(0),
    tgv(ttgv)
  { 
    int myid=0;
    int ierr=0;
    int argc=0;
    char ** argv = 0;;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    id.job=JOB_INIT;
    id.par=1;
    id.sym=A.sym() ;
    id.comm_fortran=USE_COMM_WORLD;
    mumps_c(&id);
    SetVerb();

    ICNTL(5)=0; // input matrix type 
    ICNTL(7)=7; // NUMBERING ...

    ICNTL(9)=1; // 1: A x = b, !1 : tA x = b 
    ICNTL(18) = 0; 
    id.nrhs=0; // 
    int  n = A.n;
    int nz = A.nbcoef;
    ffassert( A.n == A.m) ;

    int * irn = new int[nz];
    int * jcn = new int[nz];
    R  * a= new  R [nz];
    for(int i = 0;i<n;++i)
      for(int k = A.lg[i]; k<A.lg[i+1];++k)
      {
	irn[k]=i+1;
	jcn[k]=A.cl[k]+1;
	a[k]= A.a[k];
      }
   
    id.n = n; 
    id.nz =nz; 
    id.irn=irn; 
    id.jcn=jcn;
    id.a = (MR *) (void  *)a; 
    id.rhs = 0;



    id.job=JOB_ANA_FAC; // performs the analysis. and performs the factorization. 
    mumps_c(&id);
    Check("MUMPS-seq analayse and Factorize");
    if(verbosity>3)
      cout << "  -- MUMPS LU   n=  " << n << ", peak Mem: " << INFOG(22) << " Mb" <<  endl;

  }
  void Solver(const MatriceMorse<R> &A,KN_<R> &x,const KN_<R> &b) const  {
     ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    ffassert(A.ChecknbLine( id.n) && id.n == x.N() && A.ChecknbColumn(id.n) );
    //  convert   array in long ...
    if(verbosity>2)
      cout << " -- MUMPS solve,  peak Mem : "<< INFOG(22) << " Mb,   n = " 
	   << id.n << " sym =" << id.sym <<  endl; 
    id.nrhs = 1;
    x = b;
    id.rhs= (MR *) (void  *) (R*) x;
    SetVerb();
    id.job=JOB_SOLVE; // performs the analysis. and performs the factorization. 
    mumps_c(&id);
    Check("MUMPS-seq Solve");
    
    if(verbosity>3)
      cout << "   b min max " << b.min() << " " <<b.max() << endl;
    if(verbosity>1) cout << "   x min max " << x.min() << " " <<x.max() << endl;
  }
  void Clean() const 
  {
    if(verbosity>10)
      cout << "~SolveMUMPS_seq:" << this << endl;
    delete [] id.irn;
    delete [] id.jcn;
    delete [] id.a;
    SetVerb();
  }

  ~SolveMUMPS_seq() 
  { 
    Clean();
    id.job=JOB_END;
    mumps_c(&id); /* Terminate instance */
    int ierr = MPI_Finalize();  
  }

  void addMatMul(const KN_<R> & x, KN_<R> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<R> &) (*this) * x; 
  }
     
}; 

/*
template<>
class SolveMUMPS_seq<Complex> :   public MatriceMorse<Complex>::VirtualSolver  {
  double eps;
  mutable double  epsr;
  double tgv;



  
public:
  SolveMUMPS_seq(const MatriceMorse<Complex> &A,int strategy,double ttgv, double epsilon=1e-6,
     double pivot=-1.,double pivot_sym=-1.
) : 
    eps(epsilon),epsr(0),tgv(ttgv),
    ar(0),ai(0),

   { 
    int status;
    throwassert( !A.sym());

  }
  void Solver(const MatriceMorse<Complex> &A,KN_<Complex> &x,const KN_<Complex> &b) const  {
        ffassert ( &x[0] != &b[0]);
    epsr = (eps < 0) ? (epsr >0 ? -epsr : -eps ) : eps ;
    if(verbosity>1)
    {
      cout << "  -- MUMPS _solve,  peak Mem : " <<  -1 << "Mbytes " << endl;
      cout << "   b min max " << b.min() << " " <<b.max() << endl;
      cout << "   x min max " << x.min() << " " <<x.max() << endl;
    }
  }

  ~SolveMUMPS_seq() { 
    if(verbosity>5)
    cout << "~SolveMUMPS_seq " << endl;
  }
  void addMatMul(const KN_<Complex> & x, KN_<Complex> & Ax) const 
  {  
    ffassert(x.N()==Ax.N());
    Ax +=  (const MatriceMorse<Complex> &) (*this) * x; 
  }
     

}; 
*/

inline MatriceMorse<double>::VirtualSolver *
BuildSolverIMUMPSseq(DCL_ARG_SPARSE_SOLVER(double,A))
{
  if(verbosity>3)
    cout << " BuildSolverMUMPSseq<double>" << endl;
    return new SolveMUMPS_seq<double>(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym);
}

inline MatriceMorse<Complex>::VirtualSolver *
BuildSolverIMUMPSseq(DCL_ARG_SPARSE_SOLVER(Complex,A))
{
  if(verbosity>3)
    cout << " BuildSolverMUMPSseq<Complex>" << endl;
    return new SolveMUMPS_seq<Complex>(*A,ds.strategy,ds.tgv,ds.epsilon,ds.tol_pivot,ds.tol_pivot_sym);
}


//  the 2 default sparse solver double and complex
DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R ; ;
DefSparseSolver<Complex>::SparseMatSolver SparseMatSolver_C;
// the default probleme solver 
TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue=TypeSolveMat::defaultvalue;

bool SetDefault()
{
    if(verbosity>1)
	cout << " SetDefault sparse to default" << endl;
    DefSparseSolver<double>::solver =SparseMatSolver_R;
    DefSparseSolver<Complex>::solver =SparseMatSolver_C;
    TypeSolveMat::defaultvalue =TypeSolveMat::SparseSolver;
}

bool SetMUMPS_seq()
{
    if(verbosity>1)
	cout << " SetDefault sparse solver to MUMPSseq" << endl;
    DefSparseSolver<double>::solver  =BuildSolverIMUMPSseq;
    DefSparseSolver<Complex>::solver =BuildSolverIMUMPSseq;    
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
}


void init22()
{    
  SparseMatSolver_R= DefSparseSolver<double>::solver;
  SparseMatSolver_C= DefSparseSolver<Complex>::solver;
  if(verbosity>1)
    cout << "\n Add: MUMPS_seq:  defaultsolver defaultsolverMUMPS_seq" << endl;
  TypeSolveMat::defaultvalue=TypeSolveMat::SparseSolver; 
  DefSparseSolver<double>::solver =BuildSolverIMUMPSseq;
  DefSparseSolver<Complex>::solver =BuildSolverIMUMPSseq;
  if(! Global.Find("defaultsolver").NotNull() )
    { 
      if(verbosity)
	cout << "\n add defaultsolver (64)" << endl;
      Global.Add("defaultsolver","(",new OneOperator0<bool>(SetDefault));
    }
  if(! Global.Find("defaulttoMUMPSseq").NotNull() )
    Global.Add("defaulttoMUMPSseq","(",new OneOperator0<bool>(SetMUMPS_seq));  
}


LOADFUNC(init22);
