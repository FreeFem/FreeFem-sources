// SUMMARY  :   simpler MUMPS interface
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : P. Jolivet 
// E-MAIL   : Pierre Jolivet <pierre.jolivet@ljll.math.upmc.fr>
//
//ff-c++-LIBRARY-dep:  mumps parmetis ptscotch  scalapack blacs blas  mpifc  fc mpi  pthread 
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

#include "MatriceCreuse.hpp"

#include "dmatrix.hpp"
#include <dmumps_c.h>
#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654
#define ICNTL(I) icntl[(I)-1] /* macro s.t. indices match documentation */
#define INFOG(I) infog[(I)-1] /* macro s.t. indices match documentation */
#define INFO(I) info[(I)-1]   /* macro s.t. indices match documentation */

static std::string analysis[] = {"AMD", "", "AMF", "SCOTCH", "PORD", "METIS", "QAMD", "automatic sequential", "automatic parallel", "PT-SCOTCH", "ParMetis"};

class SolverMumps : public MatriceMorse<double>::VirtualSolver {
    private:
        mutable DMUMPS_STRUC_C* _id;
        mutable unsigned char   _strategy;

    public:
        SolverMumps(const MatriceMorse<double> &A, KN<long> &param_int, KN<double> &param_double, MPI_Comm* comm) {
            _id = new DMUMPS_STRUC_C;
            _id->job = JOB_INIT; _id->par = 1;
            if(comm)
                _id->comm_fortran = MPI_Comm_c2f(*comm);
            else
                _id->comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
            _id->sym = A.symetrique;
            _strategy = param_int.n > 0 ? param_int[0] : 3;
            dmumps_c(_id);
            int* I = NULL;
            int* J = NULL;
            double* C = NULL;
            if((param_int.n > 1 && mpirank == param_int[1]) || (param_int.n < 2 && mpirank == 0)) {
                _id->n = A.n;

                if(_id->sym == 0) {
                    _id->nz = A.nbcoef;
                    I = new int[A.nbcoef];
                    CSR2COO<'C', 'U'>(A.n, A.lg, I);
                    C = A.a;
                    J = A.cl;
                    for(unsigned int i = 0; i < A.nbcoef; ++i)
                        ++J[i];
                    _id->a = C;
                }
                else {
                    if(A.symetrique) {
                        _id->nz = A.nbcoef;
                        I = new int[A.nbcoef];
                        J = new int[A.nbcoef];
                        C = new double[A.nbcoef];
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
                        _id->a = C;
                    }
                    else {
                        _id->nz = A.n + (A.nbcoef - A.n) / 2;
                        I = new int[A.n + (A.nbcoef - A.n) / 2];
                        J = new int[A.n + (A.nbcoef - A.n) / 2];
                        C = new double[A.n + (A.nbcoef - A.n) / 2];
                        trimCSR<false, 'F'>(A.n, I + A.n, A.lg, J + A.n, A.cl, C + A.n, A.a);
                        for(unsigned int i = 0; i < A.n - 1; ++i)
                            C[i] = A.a[A.lg[i + 1] - (I[i + 1 + A.n] - I[i + A.n]) - 1];
                        C[A.n - 1] = A.a[A.nbcoef - 1];
                        std::generate(I, I + A.n, step(0, 1));
                        CSR2COO<'F', 'U'>(A.n - 1, I + A.n, I + A.n);
                        std::copy(I, I + A.n, J);
                        _id->a = C;
                    }
                }
                _id->irn = I;
                _id->jcn = J;
            }
            _id->nrhs = 1;
            _id->ICNTL(1) = 0; _id->ICNTL(2) = 0; _id->ICNTL(3) = verbosity > 1 ? 6 : 0; _id->ICNTL(4) = 0; // verbose level
            _id->ICNTL(5) = 0;                                                          // assembled format
            if(_strategy > 0 && _strategy < 9 && _strategy != 2) {
                _id->ICNTL(28) = 1;             // 1: sequential analysis
                _id->ICNTL(7)  = _strategy - 1; //     0: AMD
            }                                   //     1:
                                                //     2: AMF
                                                //     3: SCOTCH
                                                //     4: PORD
                                                //     5: METIS
                                                //     6: QAMD
                                                //     7: automatic
            else {
                _id->ICNTL(28) = 1;
                _id->ICNTL(7)  = 7;
            }
	    if(_strategy > 8 && _strategy < 12) {
	      _id->ICNTL(28) = 2;              // 2: parallel analysis
	      _id->ICNTL(29) = _strategy - 9;  //     0: automatic
            }                                   //     1: PT-STOCH
                                                //     2: ParMetis
            _id->ICNTL(9)  = 1;
            _id->ICNTL(11) = 0;                 // verbose level
            _id->ICNTL(18) = 0;                 // centralized matrix input
            _id->ICNTL(20) = 0;                 // dense RHS
            _id->ICNTL(14) = 30;                // percentage increase in the estimated working space
            _id->job = 4;
            dmumps_c(_id);
            if(_id->INFOG(1) != 0)
                std::cout << "BUG MUMPS, INFOG(1) = " << _id->INFOG(1) << std::endl;
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
        };

        void Solver(const MatriceMorse<double> &A, KN_<double> &x, const KN_<double> &b) const  {
            _id->ICNTL(21) = 0; _id->ICNTL(3) = verbosity > 1 ? 6 : 0;
            x = b;
            _id->rhs = x;
            _id->job = 3;
            dmumps_c(_id);
        };

        ~SolverMumps() {
            _id->job = JOB_END;
            dmumps_c(_id);
            if(_id)
                delete _id;
        };
};


MatriceMorse<double>::VirtualSolver* buildSolver(DCL_ARG_SPARSE_SOLVER(double, A)) {
    if(A)
        return new SolverMumps(*A, ds.lparams, ds.dparams, (MPI_Comm*)ds.commworld);
    else {
        MatriceMorse<double> empty;
        return new SolverMumps(empty, ds.lparams, ds.dparams, (MPI_Comm*)ds.commworld);
    }
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
