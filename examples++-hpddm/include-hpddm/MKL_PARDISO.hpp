/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2012-10-07

   Copyright (C) 2011-2014 Université de Grenoble
                 2015      Eidgenössische Technische Hochschule Zürich
                 2016-     Centre National de la Recherche Scientifique

   HPDDM is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   HPDDM is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with HPDDM.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _HPDDM_MKL_PARDISO_
#define _HPDDM_MKL_PARDISO_

#ifdef DMKL_PARDISO
#include <mkl_cluster_sparse_solver.h>
#endif
#ifdef MKL_PARDISOSUB
#include <mkl_pardiso.h>
#endif

namespace HPDDM {
template<class K>
struct prds {
    static constexpr int SPD = !Wrapper<K>::is_complex ? 2 : 4;
    static constexpr int SYM = !Wrapper<K>::is_complex ? -2 : 6;
    static constexpr int UNS = !Wrapper<K>::is_complex ? 1 : 3;
};

#ifdef DMKL_PARDISO
#define COARSEOPERATOR HPDDM::MklPardiso
/* Class: MKL Pardiso
 *
 *  A class inheriting from <DMatrix> to use <MKL Pardiso>.
 *
 * Template Parameter:
 *    K              - Scalar type. */
template<class K>
class MklPardiso : public DMatrix {
    private:
        /* Variable: pt
         *  Internal data pointer. */
        void*      _pt[64];
        /* Variable: a
         *  Array of data. */
        K*              _C;
        /* Variable: I
         *  Array of row pointers. */
        int*            _I;
        /* Variable: J
         *  Array of column indices. */
        int*            _J;
        /* Variable: w
         *  Workspace array. */
        K*              _w;
        /* Variable: mtype
         *  Matrix type. */
        int         _mtype;
        /* Variable: iparm
         *  Array of parameters. */
        int     _iparm[64];
        /* Variable: comm
         *  MPI communicator. */
        int          _comm;
    protected:
        /* Variable: numbering
         *  0-based indexing. */
        static constexpr char _numbering = 'C';
    public:
        MklPardiso() : _pt(), _C(), _I(), _J(), _w(), _comm(-1) { }
        ~MklPardiso() {
            delete [] _w;
            int phase = -1;
            int error;
            K ddum;
            int idum;
            if(_comm != -1)
                CLUSTER_SPARSE_SOLVER(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), &_mtype, &phase, &(DMatrix::_n), &ddum, &idum, &idum, const_cast<int*>(&i__1), const_cast<int*>(&i__1), _iparm, const_cast<int*>(&i__0), &ddum, &ddum, const_cast<int*>(&_comm), &error);
            delete [] _I;
            if(DMatrix::_communicator != MPI_COMM_NULL && DMatrix::_n == _iparm[41] - _iparm[40] + 1 && _mtype != prds<K>::SPD)
                delete [] _C;
        }
        /* Function: numfact
         *
         *  Initializes <MKL Pardiso::pt> and <MKL Pardiso::iparm>, and factorizes the supplied matrix.
         *
         * Template Parameter:
         *    S              - 'S'ymmetric or 'G'eneral factorization.
         *
         * Parameters:
         *    ncol           - Number of local rows.
         *    I              - Array of row pointers.
         *    loc2glob       - Lower and upper bounds of the local domain.
         *    J              - Array of column indices.
         *    C              - Array of data. */
        template<char S>
        void numfact(unsigned int ncol, int* I, int* loc2glob, int* J, K* C) {
            _I = I;
            _J = J;
            _C = C;
            Option& opt = *Option::get();
            if(S == 'S')
                _mtype = opt.val<unsigned short>("master_not_spd", 0) ? prds<K>::SYM : prds<K>::SPD;
            else
                _mtype = prds<K>::UNS;
            int phase, error;
            K ddum;
            std::fill_n(_iparm, 64, 0);
            _iparm[0] = 1;
#ifdef _OPENMP
            _iparm[1] = omp_get_num_threads() > 1 ? 3 : 2;
#else
            _iparm[1] = 2;
#endif
            _iparm[5] = 1;
            _iparm[9] = 13;
            _iparm[10] = 1;
            _iparm[27] = std::is_same<double, underlying_type<K>>::value ? 0 : 1;
            _iparm[34] = (_numbering == 'C');
            _iparm[39] = 2;
            _iparm[40] = loc2glob[0];
            _iparm[41] = loc2glob[1];
            delete [] loc2glob;
            phase = 12;
            CLUSTER_SPARSE_SOLVER(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), &_mtype, &phase, &(DMatrix::_n), C, _I, _J, const_cast<int*>(&i__1), const_cast<int*>(&i__1), _iparm, opt.val<int>("verbosity") < 2 ? const_cast<int*>(&i__0) : const_cast<int*>(&i__1), &ddum, &ddum, const_cast<int*>(&_comm), &error);
            _w = new K[_iparm[41] - _iparm[40] + 1];
        }
        /* Function: solve
         *
         *  Solves the system in-place.
         *
         * Template Parameter:
         *    D              - Distribution of right-hand sides and solution vectors.
         *
         * Parameters:
         *    rhs            - Input right-hand sides, solution vectors are stored in-place.
         *    n              - Number of right-hand sides. */
        template<DMatrix::Distribution D>
        void solve(K* rhs, const unsigned short& n) {
            int error;
            int phase = 33;
            int nrhs = n;
            K ddum;
            CLUSTER_SPARSE_SOLVER(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), &_mtype, &phase, &(DMatrix::_n), _C, _I, _J, const_cast<int*>(&i__1), &nrhs, _iparm, const_cast<int*>(&i__0), rhs, _w, const_cast<int*>(&_comm), &error);
        }
        void initialize() {
            if(DMatrix::_communicator != MPI_COMM_NULL)
                _comm = MPI_Comm_c2f(DMatrix::_communicator);
            DMatrix::initialize("PARDISO", { DISTRIBUTED_SOL_AND_RHS });
        }
};
#endif // DMKL_PARDISO

#ifdef MKL_PARDISOSUB
#define SUBDOMAIN HPDDM::MklPardisoSub
template<class K>
class MklPardisoSub {
    private:
        mutable void*  _pt[64];
        K*                  _C;
        int*                _I;
        int*                _J;
        K*                  _w;
        int             _mtype;
        mutable int _iparm[64];
        int                 _n;
        int           _partial;
    public:
        MklPardisoSub() : _pt(), _C(), _I(), _J(), _w(), _partial() { }
        MklPardisoSub(const MklPardisoSub&) = delete;
        ~MklPardisoSub() {
            delete [] _w;
            _w = nullptr;
            int phase = -1;
            int error;
            int idum;
            K ddum;
            _n = 1;
            PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), &_mtype, &phase, &_n, &ddum, &idum, &idum, const_cast<int*>(&i__1), const_cast<int*>(&i__1), _iparm, const_cast<int*>(&i__0), &ddum, &ddum, &error);
            if(_mtype == prds<K>::SPD || _mtype == prds<K>::SYM) {
                delete [] _I;
                delete [] _J;
            }
            if(_mtype == prds<K>::SYM)
                delete [] _C;
        }
        static constexpr char _numbering = 'F';
        template<char N = HPDDM_NUMBERING>
        void numfact(MatrixCSR<K>* const& A, bool detection = false, K* const& schur = nullptr) {
            static_assert(N == 'C' || N == 'F', "Unknown numbering");
            int* perm = nullptr;
            int phase, error;
            K ddum;
            if(!_w) {
                _n = A->_n;
                std::fill_n(_iparm, 64, 0);
                _iparm[0] = 1;
                Option& opt = *Option::get();
                for(unsigned short i : { 1, 7, 9, 10, 12, 20, 23, 26 }) {
                    int val = opt.val<int>("mkl_pardiso_iparm_" + to_string(i + 1));
                    if(val != std::numeric_limits<int>::lowest())
                        _iparm[i] = val;

                }
                _iparm[27] = std::is_same<double, underlying_type<K>>::value ? 0 : 1;
                _iparm[34] = (N == 'C');
                phase = 12;
                if(A->_sym) {
                    _I = new int[_n + 1];
                    _J = new int[A->_nnz];
                    _C = new K[A->_nnz];
                }
                else
                    _mtype = prds<K>::UNS;
                if(schur != nullptr) {
                    _iparm[35] = 2;
                    perm = new int[_n];
                    _partial = static_cast<int>(std::real(schur[1]));
                    std::fill_n(perm, _partial, 0);
                    std::fill(perm + _partial, perm + _n, 1);
                }
                _w = new K[_n];
            }
            else {
                if(_mtype == prds<K>::SPD)
                    _C = new K[A->_nnz];
                phase = 22;
            }
            if(A->_sym) {
                _mtype = Wrapper<K>::is_complex || detection ? prds<K>::SYM : prds<K>::SPD;
                Wrapper<K>::template csrcsc<N, N>(&_n, A->_a, A->_ja, A->_ia, _C, _J, _I);
            }
            else {
                _I = A->_ia;
                _J = A->_ja;
                _C = A->_a;
            }
            PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), &_mtype, &phase,
                    const_cast<int*>(&_n), _C, _I, _J, perm, const_cast<int*>(&i__1), _iparm, const_cast<int*>(&i__0), &ddum, schur, &error);
            delete [] perm;
            if(_mtype == prds<K>::SPD)
                delete [] _C;
        }
        void solve(K* x) const {
            int error;
            _iparm[5] = 1;
            if(!_partial) {
                int phase = 33;
                PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), const_cast<int*>(&_mtype), &phase, const_cast<int*>(&_n), _C, _I, _J, const_cast<int*>(&i__1), const_cast<int*>(&i__1), _iparm, const_cast<int*>(&i__0), x, const_cast<K*>(_w), &error);
            }
            else {
                int phase = 331;
                PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), const_cast<int*>(&_mtype), &phase, const_cast<int*>(&_n), _C, _I, _J, const_cast<int*>(&i__1), const_cast<int*>(&i__1), _iparm, const_cast<int*>(&i__0), x, const_cast<K*>(_w), &error);
                std::fill(x + _partial, x + _n, K());
                phase = 333;
                PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), const_cast<int*>(&_mtype), &phase, const_cast<int*>(&_n), _C, _I, _J, const_cast<int*>(&i__1), const_cast<int*>(&i__1), _iparm, const_cast<int*>(&i__0), x, const_cast<K*>(_w), &error);
            }
        }
        void solve(const K* const b, K* const x) const {
            int error;
            if(!_partial) {
                _iparm[5] = 0;
                int phase = 33;
                PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), const_cast<int*>(&_mtype), &phase, const_cast<int*>(&_n), _C, _I, _J, const_cast<int*>(&i__1), const_cast<int*>(&i__1), _iparm, const_cast<int*>(&i__0), const_cast<K*>(b), x, &error);
            }
            else {
                _iparm[5] = 1;
                int phase = 331;
                std::copy_n(b, _partial, x);
                PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), const_cast<int*>(&_mtype), &phase, const_cast<int*>(&_n), _C, _I, _J, const_cast<int*>(&i__1), const_cast<int*>(&i__1), _iparm, const_cast<int*>(&i__0), x, const_cast<K*>(_w), &error);
                std::fill(x + _partial, x + _n, K());
                phase = 333;
                PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), const_cast<int*>(&_mtype), &phase, const_cast<int*>(&_n), _C, _I, _J, const_cast<int*>(&i__1), const_cast<int*>(&i__1), _iparm, const_cast<int*>(&i__0), x, const_cast<K*>(_w), &error);
            }
        }
        void solve(K* const x, const unsigned short& n) const {
            int error;
            int phase = 33;
            int nrhs = n;
            _iparm[5] = 1;
            K* w = new K[_n * n];
            PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), const_cast<int*>(&_mtype), &phase, const_cast<int*>(&_n), _C, _I, _J, const_cast<int*>(&i__1), &nrhs, _iparm, const_cast<int*>(&i__0), x, w, &error);
            delete [] w;
        }
        void solve(const K* const b, K* const x, const unsigned short& n) const {
            int error;
            int phase = 33;
            int nrhs = n;
            _iparm[5] = 0;
            PARDISO(_pt, const_cast<int*>(&i__1), const_cast<int*>(&i__1), const_cast<int*>(&_mtype), &phase, const_cast<int*>(&_n), _C, _I, _J, const_cast<int*>(&i__1), &nrhs, _iparm, const_cast<int*>(&i__0), const_cast<K*>(b), x, &error);
        }
};
#endif // MKL_PARDISOSUB
} // HPDDM
#endif // _HPDDM_MKL_PARDISO_
