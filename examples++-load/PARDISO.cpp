/****************************************************************************/
/* This file is part of FreeFem++.                                          */
/*                                                                          */
/* FreeFem++ is free software: you can redistribute it and/or modify        */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFem++ is distributed in the hope that it will be useful,             */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFem++. If not, see <http://www.gnu.org/licenses/>.        */
/****************************************************************************/
// SUMMARY : PARDISO interface
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Pierre Jolivet
// E-MAIL  : pierre.joliver@enseeiht.fr

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep: mkl
//ff-c++-cpp-dep:
// *INDENT-ON* //

// #include <mpi.h>
#include <mkl_pardiso.h>
#include <mkl_spblas.h>
#include <mkl_types.h>
#if 0
#include <omp.h>
#else

extern "C" {
    extern int omp_get_max_threads (void);
    extern int omp_get_num_threads (void);
    extern void omp_set_num_threads (int);
}
#endif
#include "rgraph.hpp"
#include "AFunction.hpp"

#include "MatriceCreuse.hpp"
#include "dmatrix.hpp"

template<typename RR> struct  PARDISO_STRUC_TRAIT {typedef void R; static const int unSYM = 0; static const int SYM = 0;};
template<> struct PARDISO_STRUC_TRAIT<double>  {typedef double R; static const int unSYM = 11; static const int SYM = 2;};
template<> struct PARDISO_STRUC_TRAIT<Complex>  {typedef Complex R; static const int unSYM = 13; static const int SYM = -4;};

void mkl_csrcsc (MKL_INT *job, MKL_INT *n, Complex *Acsr, MKL_INT *AJ0, MKL_INT *AI0, Complex *Acsc, MKL_INT *AJ1, MKL_INT *AI1, MKL_INT *info)
{mkl_zcsrcsc(job, n, reinterpret_cast<MKL_Complex16 *>(Acsr), AJ0, AI0, reinterpret_cast<MKL_Complex16 *>(Acsc), AJ1, AI1, info);}

void mkl_csrcsc (MKL_INT *job, MKL_INT *n, double *Acsr, MKL_INT *AJ0, MKL_INT *AI0, double *Acsc, MKL_INT *AJ1, MKL_INT *AI1, MKL_INT *info)
{mkl_dcsrcsc(job, n, Acsr, AJ0, AI0, Acsc, AJ1, AI1, info);}

template<class R>
class SolverPardiso: public  VirtualSolver<int,R> {
private:
    typedef HashMatrix<int,R> HMat;
    mutable void *_pt[64];
    mutable MKL_INT mtype;
    mutable MKL_INT _iparm[64];
    mutable HMat *ptA;
    MKL_INT *_I;
    MKL_INT *_J;
    R *_C; // Coef
    long verb;
    long cn,cs;
    MKL_INT pmtype;
        MKL_INT phase, error, msglvl, maxfct, mnum;
public:
     static const int orTypeSol = 1&2&8&16;// Do all
    static const MKL_INT pmtype_unset= -1000000000;
    typedef typename PARDISO_STRUC_TRAIT<R>::R MR;
    
    SolverPardiso (HMat  &AH, const Data_Sparse_Solver & ds,Stack stack )
    : ptA(&AH),_I(0),_J(0),_C(0),verb(ds.verb),cn(0),cs(0),pmtype(pmtype_unset)
    {
        if( ds.lparams.N()>1) pmtype=ds.lparams[0]; // bof bof ...
        msglvl=0;
        error=0;
        phase =0;
    //(const MatriceMorse<R> &A, KN<long> &param_int, KN<double> &param_double) {
    }

void SetState()
{
    if(verb>2 || verbosity> 9) std::cout << " SetState "<< ptA-> re_do_numerics << " " << ptA-> re_do_symbolic <<std::endl;
    if( ptA->GetReDoNumerics() ) cn++;
    if( ptA->GetReDoSymbolic()) cs++;
    this->CheckState(ptA->n,cs,cn);
}
void fac_init()
{
    MR ddum;
    maxfct = 1; /* Maximum number of numerical factorizations. */
    mnum = 1; /* Which factorization to use. */
    if (ptA->half) {
        mtype = PARDISO_STRUC_TRAIT<R>::SYM;
    } else {
        if (pmtype != pmtype_unset )
            mtype = pmtype;// Debil .... FH
        else
            mtype = PARDISO_STRUC_TRAIT<R>::unSYM;
    }
    
    for (unsigned short i = 0; i < 64; ++i) {
        _iparm[i] = 0;
        _pt[i] = NULL;
    }
    
    _iparm[0] = 1;    /* No solver default */
    _iparm[1] = 3;    /* Fill-in reordering from METIS */
    _iparm[2] = 1;
    _iparm[3] = 0;    /* No iterative-direct algorithm */
    _iparm[4] = 0;    /* No user fill-in reducing permutation */
    _iparm[5] = 0;    /* Write solution into rhs */
    _iparm[6] = 0;    /* Not in use */
    _iparm[7] = 0;    /* Max numbers of iterative refinement steps */
    _iparm[8] = 0;    /* Not in use */
    _iparm[9] = 13;    /* Perturb the pivot elements with 1E-13 */
    _iparm[10] = 1;    /* Use nonsymmetric permutation and scaling MPS */
    _iparm[11] = 0;    /* Not in use */
    _iparm[12] = 0;    /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try _iparm[12] = 1 in case of inappropriate accuracy */
    _iparm[13] = 0;    /* Output: Number of perturbed pivots */
    _iparm[14] = 0;    /* Not in use */
    _iparm[15] = 0;    /* Not in use */
    _iparm[16] = 0;    /* Not in use */
    _iparm[17] = -1;/* Output: Number of nonzeros in the factor LU */
    _iparm[18] = -1;/* Output: Mflops for LU factorization */
    _iparm[19] = 0;    /* Output: Numbers of CG Iterations */
    _iparm[34] = 1;
    msglvl = 0;    /* Print statistical information in file */
    if (verbosity > 1) {
        msglvl = 1;
    }
    
}
void fac_symbolic()
{// phase 11

    

    error = 0;    /* Initialize error flag */
    MKL_INT one = 1, nrhs=0;
    MKL_INT n = ptA->n;
    
    if (mtype != 2) {
        cout << " Pardiso :CSR  "<< endl;
        ptA->CSR();
        _I = ptA->p;
        _J = ptA->i;
        _C = ptA->aij;
    } else {
        // WARNING PARDISO USE just UPPER TRAINGular PART in CSR format
        
        if (ptA->half) {
             cout << " Pardiso :do2Triangular "<< endl;
            ptA->do2Triangular(false);// triangular Upper ...
            ptA->CSR();
            _I = ptA->p;
            _J = ptA->i;
            _C = ptA->aij;

        } else {
            cout << " Pardiso : Copie mat ????? "<< endl;
            ptA->CSR();
            _I = new MKL_INT[n + 1];
            _J = new MKL_INT[n + (ptA->nnz - n) / 2];
            _C = new R[n + (ptA->nnz - n) / 2];
            trimCSR<true, 'C', R>(n, _I, ptA->p, _J, ptA->j, _C, ptA->aij);
        }
    }
    if(verb>2 || verbosity> 9) cout << "fac_symbolic PARDISO  nnz U " << " nnz= "  << n<< " " << ptA->nnz << " "
        <<   endl;
    MR ddum;
    MR * _CC =  reinterpret_cast<MR *>((R *)_C);
    phase = 11;
    PARDISO(_pt,  &maxfct, &mnum, &mtype, &phase,
            &n, _CC, _I, _J, &one, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error);
    printf("\nReordering completed ... ");
    printf("\nNumber of nonzeros in factors = %d", _iparm[17]);
    printf("\nNumber of factorization MFLOPS = %d", _iparm[18]);


    cout << " PARDISO" << error << endl;
}
void fac_numeric()
{
    MKL_INT nrhs = 0;
    MKL_INT idum;
    error = 0;
    phase = 22;
    MR ddum;
    MKL_INT n = ptA->n;
    if(verb>2 || verbosity> 9) cout << "fac_NUMERIC PARDISO  R: nnz U " << " nnz= "  << ptA->nnz << endl;
    MR * _CC =  reinterpret_cast<MR *>((R *)_C);
    PARDISO(_pt, &maxfct, &mnum, &mtype, &phase,
            &n, _CC , _I, _J,  &idum, &nrhs, _iparm, &msglvl, &ddum, &ddum, &error);
    /*
     phase = 22;
     PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
     &n, a, ia, ja, &idum, &nrhs,
     iparm, &msglvl, &ddum, &ddum, &error);

     */
  

}
void dosolver(R *x,R*b,int N,int trans)
{
    MKL_INT nrhs = N;
     msglvl = 0;
    MKL_INT one = 1;
    MKL_INT idum;
    
    if (verbosity > 1) {
        msglvl = 1;
    }
    
     error = 0;
     phase = 33;
    MR ddum;
    MKL_INT n = ptA->n;
       if(verb>2 || verbosity> 9) cout << "dosolver PARDISO  R: nnz U " << " nnz= "  << ptA->nnz << " "<< N << endl;
    MR * _CC =  reinterpret_cast<MR *>((R *)_C);
    PARDISO(_pt, &maxfct, &mnum, &mtype, &phase, &n,_CC, _I, _J, &idum, &nrhs, _iparm, &msglvl, reinterpret_cast<MR *>((R *)b), reinterpret_cast<MR *>((R *)x), &error);
}

~SolverPardiso () {
    MKL_INT phase = -1;
    MKL_INT one = 1;
    MKL_INT msglvl = 0;
    MKL_INT error;
    MR ddum;
    MKL_INT idum;
    MKL_INT n = ptA->n;
    
    PARDISO(_pt,&maxfct, &mnum, &mtype, &phase, &n, &ddum, &idum, &idum, &one, &one, _iparm, &msglvl, &ddum, &ddum, &error);
    if (mtype == 2) {
        if (_I)
            delete [] _I;
        
        
        if (_J)
            delete [] _J;
        
        if (_C != ptA->aij)
            delete [] _C;

    }
}
};




static long ffompgetnumthreads () {return omp_get_num_threads();}

static long ffompgetmaxthreads () {return omp_get_max_threads();}

static long ffompsetnumthreads (long n) {omp_set_num_threads(n); return n;}

static void Load_Init ()
{
    // }static void initPARDISO()
    // {
    addsolver<SolverPardiso<double>>("PARDISO",50,3);
    addsolver<SolverPardiso<Complex>>("PARDISO",50,3);

    
    
    if (!Global.Find("ompsetnumthreads").NotNull()) {
        Global.Add("ompsetnumthreads", "(", new OneOperator1<long, long>(ffompsetnumthreads));
    }
    
    if (!Global.Find("ompgetnumthreads").NotNull()) {
        Global.Add("ompgetnumthreads", "(", new OneOperator0<long>(ffompgetnumthreads));
    }
    
    if (!Global.Find("ompgetmaxthreads").NotNull()) {
        Global.Add("ompgetmaxthreads", "(", new OneOperator0<long>(ffompgetmaxthreads));
    }
}

// LOADFUNC(initPARDISO);
LOADFUNC(Load_Init)
