// SUMMARY  :   PARDISO interface
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : P. Jolivet 
// E-MAIL   : Pierre Jolivet <pierre.jolivet@ljll.math.upmc.fr>
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
//ff-c++-LIBRARY-dep:  mkl
//#include <mpi.h>
#include <mkl_pardiso.h>
#include <mkl_spblas.h>
#include <mkl_types.h>
#if 0
#include <omp.h>
#else

extern "C" {
extern int    omp_get_max_threads  (void);
extern int     omp_get_num_threads  (void);
extern void    omp_set_num_threads (int);
}
#endif
#include "rgraph.hpp"
#include "AFunction.hpp"

#include "MatriceCreuse.hpp"
#include "dmatrix.hpp"

template<typename RR> struct  PARDISO_STRUC_TRAIT {  typedef void R;static  const int unSYM=0;static const int SYM=0;};
template<> struct PARDISO_STRUC_TRAIT<double>  { typedef double R;static const int unSYM=11;static const int SYM=2;};
template<> struct PARDISO_STRUC_TRAIT<Complex>  { typedef Complex R;static const int unSYM=13;static const int SYM=-4;};

void mkl_csrcsc (MKL_INT * job, MKL_INT * n, Complex *Acsr, MKL_INT * AJ0, MKL_INT * AI0, Complex *Acsc, MKL_INT * AJ1, MKL_INT * AI1, MKL_INT * info)
{ mkl_zcsrcsc (job,n,reinterpret_cast<MKL_Complex16*>(Acsr),AJ0,AI0,reinterpret_cast<MKL_Complex16*>(Acsc),AJ1,AI1,info);}

void mkl_csrcsc (MKL_INT * job, MKL_INT * n, double *Acsr, MKL_INT * AJ0, MKL_INT * AI0, double *Acsc, MKL_INT * AJ1, MKL_INT * AI1, MKL_INT * info)
{ mkl_dcsrcsc (job,n,Acsr,AJ0,AI0,Acsc,AJ1,AI1,info);}


template<class R>
class SolverPardiso : public MatriceMorse<R>::VirtualSolver {
    private:
        mutable void* _pt[64];
        mutable MKL_INT _mtype;
        mutable MKL_INT _iparm[64];
        mutable MatriceMorse<R>* ptA;
        MKL_INT* _I;
        MKL_INT* _J;
    
    public:
    typedef typename  PARDISO_STRUC_TRAIT<R>::R MR;

        SolverPardiso(const MatriceMorse<R> &A, KN<long> &param_int, KN<double> &param_double) {
            ptA = (MatriceMorse<R>*)(&A);
            MKL_INT phase, error, msglvl;
            MR ddum;
            R* _C;
            if(A.symetrique)
                _mtype =  PARDISO_STRUC_TRAIT<R>::SYM;
            else {
                if(param_int)
                    _mtype = param_int[0];
                else
                    _mtype = PARDISO_STRUC_TRAIT<R>::unSYM;
            }
            for (unsigned short i = 0; i < 64; ++i) {
                _iparm[i] = 0;
                _pt[i] = NULL;
            }
            _iparm[0] = 1;      /* No solver default */
            _iparm[1] = 3;      /* Fill-in reordering from METIS */
            _iparm[2] = 1;
            _iparm[3] = 0;      /* No iterative-direct algorithm */
            _iparm[4] = 0;      /* No user fill-in reducing permutation */
            _iparm[5] = 0;      /* Write solution into rhs */
            _iparm[6] = 0;      /* Not in use */
            _iparm[7] = 0;      /* Max numbers of iterative refinement steps */
            _iparm[8] = 0;      /* Not in use */
            _iparm[9] = 13;     /* Perturb the pivot elements with 1E-13 */
            _iparm[10] = 1;     /* Use nonsymmetric permutation and scaling MPS */
            _iparm[11] = 0;     /* Not in use */
            _iparm[12] = 0;     /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try _iparm[12] = 1 in case of inappropriate accuracy */
            _iparm[13] = 0;     /* Output: Number of perturbed pivots */
            _iparm[14] = 0;     /* Not in use */
            _iparm[15] = 0;     /* Not in use */
            _iparm[16] = 0;     /* Not in use */
            _iparm[17] = -1;    /* Output: Number of nonzeros in the factor LU */
            _iparm[18] = -1;    /* Output: Mflops for LU factorization */
            _iparm[19] = 0;     /* Output: Numbers of CG Iterations */
            _iparm[34] = 1;
            msglvl = 0;         /* Print statistical information in file */
            if(verbosity > 1)
                msglvl = 1;
            error = 0;          /* Initialize error flag */
            phase = 12;
            MKL_INT one = 1;
            MKL_INT n = A.n;
            if(_mtype != 2) {
                _I = A.lg;
                _J = A.cl;
                _C = A.a;
            }
            else {
                if(A.symetrique) {
                    int job[6] = { 0, 0, 0, 0, 0, 1 };
                    _I = new MKL_INT[n + 1];
                    _J = new MKL_INT[A.nbcoef];
                    _C = new R[A.nbcoef];
                    mkl_csrcsc(job, &n, reinterpret_cast<MR*>( A.a), A.cl, A.lg, reinterpret_cast<MR*>(_C), _J, _I, &error);
                }
                else {
                    _I = new MKL_INT[n + 1];
                    _J = new MKL_INT[n + (A.nbcoef - n) / 2];
                    _C = new R[n + (A.nbcoef - n) / 2];
                    trimCSR<true, 'C',R>(n, _I, A.lg, _J, A.cl, _C, A.a);
                }
            }
            PARDISO(_pt, &one, &one, &_mtype, &phase,
                    &n, reinterpret_cast<MR*>(_C), _I, _J, &one, &one, _iparm, &msglvl, &ddum, &ddum, &error);
            if(_C != A.a)
                delete [] _C;
        };

        void Solver(const MatriceMorse<R> &A, KN_<R> &x, const KN_<R> &b) const  {
            MKL_INT one    = 1;
            MKL_INT msglvl = 0;
            if(verbosity > 1)
                msglvl = 1;
            MKL_INT error  = 0;
            MKL_INT phase = 33;
            MR ddum;
            MKL_INT n = A.n;
            PARDISO(_pt, &one, &one, &_mtype, &phase, &n, &ddum, _I, _J, &one, &one, _iparm, &msglvl, reinterpret_cast<MR*>((R*)b), reinterpret_cast<MR*>((R*)x), &error);
        };

        ~SolverPardiso() {
            MKL_INT phase = -1;
            MKL_INT one = 1;
            MKL_INT msglvl = 0;
            MKL_INT error;
            MR ddum;
            MKL_INT idum;
            MKL_INT n = ptA->n;
            PARDISO(_pt, &one, &one, &_mtype, &phase, &n, &ddum, &idum, &idum, &one, &one, _iparm, &msglvl, &ddum, &ddum, &error);
            if(_mtype == 2) {
                if(_I)
                    delete [] _I;
                if(_J)
                    delete [] _J;
            }
        };
};

template<class R>
typename MatriceMorse<R>::VirtualSolver* buildSolver(DCL_ARG_SPARSE_SOLVER(R, A)) {
    return new SolverPardiso<R>(*A, ds.lparams, ds.dparams);
}
extern  TypeSolveMat::TSolveMat  TypeSolveMatdefaultvalue;
bool SetPARDISO()
{
    if(verbosity>1)
        cout << " SetDefault sparse solver to MUMPS" << endl;
    DefSparseSolver<double>::solver  = buildSolver<double>;
    DefSparseSolver<Complex>::solver = buildSolver<Complex>;
    DefSparseSolverSym<double>::solver  = buildSolver<double>;
    DefSparseSolverSym<Complex>::solver = buildSolver<Complex>;
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
    return 0;
}
/*  class Init {
    public:
        Init();
};



$1 */




static long ffompgetnumthreads(){return omp_get_num_threads();}
static long ffompgetmaxthreads(){return omp_get_max_threads();}

static long ffompsetnumthreads(long n){omp_set_num_threads(n); return n;}

static void Load_Init() {
  //}static void initPARDISO()
  //{
    
    if(verbosity>1)
        cout << "\n Add: PARDISO:  defaultsolver defaultsolverPARDISO" << endl;
    TypeSolveMat::defaultvalue = TypeSolveMat::SparseSolver;
    DefSparseSolver<double>::solver = buildSolver<double>;
    DefSparseSolver<Complex>::solver = buildSolver<Complex>;

    DefSparseSolver<double>::solver  = buildSolver;
    DefSparseSolver<Complex>::solver = buildSolver;
    DefSparseSolverSym<double>::solver  = buildSolver;
    DefSparseSolverSym<Complex>::solver = buildSolver;
    TypeSolveMat::defaultvalue =TypeSolveMatdefaultvalue;
     if(! Global.Find("defaulttoPARDISO").NotNull() )
        Global.Add("defaulttoPARDISO","(",new OneOperator0<bool>(SetPARDISO));
    if(! Global.Find("ompsetnumthreads").NotNull() )
    Global.Add("ompsetnumthreads","(",new OneOperator1<long,long>(ffompsetnumthreads));
    if(! Global.Find("ompgetnumthreads").NotNull() )
    Global.Add("ompgetnumthreads","(",new OneOperator0<long>(ffompgetnumthreads));
    if(! Global.Find("ompgetmaxthreads").NotNull() )
        Global.Add("ompgetmaxthreads","(",new OneOperator0<long>(ffompgetmaxthreads));

}


//LOADFUNC(initPARDISO);
LOADFUNC(Load_Init)
