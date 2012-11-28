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

#include <mpi.h>
#include <mkl_pardiso.h>
#include <mkl_spblas.h>
#include <mkl_types.h>
#include <omp.h>

#include "rgraph.hpp"
#include "AFunction.hpp"

#include "MatriceCreuse.hpp"
#include "dmatrix.hpp"

class SolverPardiso : public MatriceMorse<double>::VirtualSolver {
    private:
        mutable void* _pt[64];
        mutable MKL_INT _mtype;
        mutable MKL_INT _iparm[64];
        mutable MatriceMorse<double>* ptA;
        MKL_INT* _I;
        MKL_INT* _J;

    public:
        SolverPardiso(const MatriceMorse<double> &A, KN<long> &param_int, KN<double> &param_double, MPI_Comm* comm) {
            ptA = (MatriceMorse<double>*)(&A);
            MKL_INT phase, error, msglvl;
            double ddum;
            double* _C;
            if(A.symetrique)
                _mtype = 2;
            else {
                if(param_int)
                    _mtype = param_int[0];
                else
                    _mtype = 11;
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
                    _C = new double[A.nbcoef];
                    mkl_dcsrcsc(job, &n, A.a, A.cl, A.lg, _C, _J, _I, &error);
                }
                else {
                    _I = new MKL_INT[n + 1];
                    _J = new MKL_INT[n + (A.nbcoef - n) / 2];
                    _C = new double[n + (A.nbcoef - n) / 2];
                    trimCSR<true, 'C'>(n, _I, A.lg, _J, A.cl, _C, A.a);
                }
            }
            PARDISO(_pt, &one, &one, &_mtype, &phase,
                    &n, _C, _I, _J, &one, &one, _iparm, &msglvl, &ddum, &ddum, &error);
            if(_C != A.a)
                delete [] _C;
        };

        void Solver(const MatriceMorse<double> &A, KN_<double> &x, const KN_<double> &b) const  {
            MKL_INT one    = 1;
            MKL_INT msglvl = 0;
            if(verbosity > 1)
                msglvl = 1;
            MKL_INT error  = 0;
            MKL_INT phase = 33;
            double ddum;
            MKL_INT n = A.n;
            PARDISO(_pt, &one, &one, &_mtype, &phase, &n, &ddum, _I, _J, &one, &one, _iparm, &msglvl, b, x, &error);
        };

        ~SolverPardiso() {
            MKL_INT phase = -1;
            MKL_INT one = 1;
            MKL_INT msglvl = 0;
            MKL_INT error;
            double ddum;
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


MatriceMorse<double>::VirtualSolver* buildSolver(DCL_ARG_SPARSE_SOLVER(double, A)) {
    return new SolverPardiso(*A, ds.lparams, ds.dparams, (MPI_Comm*)ds.commworld);
}

class Init {
    public:
        Init();
};

DefSparseSolver<double>::SparseMatSolver SparseMatSolver_R;

LOADINIT(Init);

Init::Init() {
    TypeSolveMat::defaultvalue = TypeSolveMat::SparseSolver;
    DefSparseSolver<double>::solver = buildSolver;
}
