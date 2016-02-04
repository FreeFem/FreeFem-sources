/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2012-10-04

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

#ifndef _HPDDM_PASTIX_
#define _HPDDM_PASTIX_

extern "C" {
#include <pastix.h>
#include <cscd_utils.h>
}

#define HPDDM_GENERATE_PASTIX(C, T)                                                                          \
template<>                                                                                                   \
struct pstx<T> {                                                                                             \
    static void dist(pastix_data_t** pastix_data, MPI_Comm pastix_comm, pastix_int_t n, pastix_int_t* colptr,\
                     pastix_int_t* row, T* avals, pastix_int_t* loc2glob, pastix_int_t* perm,                \
                     pastix_int_t* invp, T* b, pastix_int_t rhs, pastix_int_t* iparm, double* dparm) {       \
        C ## _dpastix(pastix_data, pastix_comm,                                                              \
                      n, colptr, row, avals, loc2glob, perm, invp, b, rhs, iparm, dparm);                    \
    }                                                                                                        \
    static void seq(pastix_data_t** pastix_data, MPI_Comm pastix_comm, pastix_int_t n, pastix_int_t* colptr, \
                    pastix_int_t* row, T* avals, pastix_int_t* perm, pastix_int_t* invp, T* b,               \
                    pastix_int_t rhs, pastix_int_t* iparm, double* dparm) {                                  \
        C ## _pastix(pastix_data, pastix_comm, n, colptr, row, avals, perm, invp, b, rhs, iparm, dparm);     \
    }                                                                                                        \
    static pastix_int_t cscd_redispatch(pastix_int_t n, pastix_int_t* ia, pastix_int_t* ja, T* a, T* rhs,    \
                                        pastix_int_t nrhs, pastix_int_t* l2g, pastix_int_t dn,               \
                                        pastix_int_t** dia, pastix_int_t** dja, T** da, T** drhs,            \
                                        pastix_int_t* dl2g, MPI_Comm comm, pastix_int_t dof) {               \
        return C ## _cscd_redispatch(n, ia, ja, a, rhs, nrhs, l2g, dn, dia, dja, da, drhs, dl2g, comm, dof); \
    }                                                                                                        \
    static void initParam(pastix_int_t* iparm, double* dparm) {                                              \
        C ## _pastix_initParam(iparm, dparm);                                                                \
    }                                                                                                        \
    static pastix_int_t getLocalNodeNbr(pastix_data_t** pastix_data) {                                       \
        return C ## _pastix_getLocalNodeNbr(pastix_data);                                                    \
    }                                                                                                        \
    static pastix_int_t getLocalNodeLst(pastix_data_t** pastix_data, pastix_int_t* nodelst) {                \
        return C ## _pastix_getLocalNodeLst(pastix_data, nodelst);                                           \
    }                                                                                                        \
    static pastix_int_t setSchurUnknownList(pastix_data_t* pastix_data, pastix_int_t n, pastix_int_t* list) {\
        return C ## _pastix_setSchurUnknownList(pastix_data, n, list);                                       \
    }                                                                                                        \
    static pastix_int_t setSchurArray(pastix_data_t* pastix_data, T* array) {                                \
        return C ## _pastix_setSchurArray(pastix_data, array);                                               \
    }                                                                                                        \
};

namespace HPDDM {
template<class>
struct pstx { };
HPDDM_GENERATE_PASTIX(s, float)
HPDDM_GENERATE_PASTIX(d, double)
#ifdef PASTIX_HAS_COMPLEX
HPDDM_GENERATE_PASTIX(c, std::complex<float>)
HPDDM_GENERATE_PASTIX(z, std::complex<double>)
#endif

#ifdef DPASTIX
#define COARSEOPERATOR HPDDM::Pastix
/* Class: Pastix
 *
 *  A class inheriting from <DMatrix> to use <Pastix>.
 *
 * Template Parameter:
 *    K              - Scalar type. */
template<class K>
class Pastix : public DMatrix {
    private:
        /* Variable: data
         *  Internal data pointer. */
        pastix_data_t*      _data;
        /* Variable: values2
         *  Array of data. */
        K*               _values2;
        /* Variable: dparm
         *  Array of double-precision floating-point parameters. */
        double*            _dparm;
        /* Variable: ncol2
         *  Number of local rows. */
        pastix_int_t       _ncol2;
        /* Variable: colptr2
         *  Array of row pointers. */
        pastix_int_t*    _colptr2;
        /* Variable: rows2
         *  Array of column indices. */
        pastix_int_t*      _rows2;
        /* Variable: loc2glob2
         *  Local to global numbering. */
        pastix_int_t*  _loc2glob2;
        /* Variable: iparm
         *  Array of integer parameters. */
        pastix_int_t*      _iparm;
    protected:
        /* Variable: numbering
         *  1-based indexing. */
        static constexpr char _numbering = 'F';
    public:
        Pastix() : _data(), _values2(), _dparm(), _colptr2(), _rows2(), _loc2glob2(), _iparm() { }
        ~Pastix() {
            free(_rows2);
            free(_values2);
            delete [] _loc2glob2;
            free(_colptr2);
            if(_iparm) {
                _iparm[IPARM_START_TASK]          = API_TASK_CLEAN;
                _iparm[IPARM_END_TASK]            = API_TASK_CLEAN;

                pstx<K>::dist(&_data, DMatrix::_communicator,
                              0, NULL, NULL, NULL, NULL,
                              NULL, NULL, NULL, 1, _iparm, _dparm);
                delete [] _iparm;
                delete [] _dparm;
            }
        }
        /* Function: numfact
         *
         *  Initializes <Pastix::iparm> and <Pastix::dparm>, and factorizes the supplied matrix.
         *
         * Template Parameter:
         *    S              - 'S'ymmetric or 'G'eneral factorization.
         *
         * Parameters:
         *    ncol           - Number of local rows.
         *    I              - Array of row pointers.
         *    loc2glob       - Local to global numbering.
         *    J              - Array of column indices.
         *    C              - Array of data. */
        template<char S>
        void numfact(unsigned int ncol, int* I, int* loc2glob, int* J, K* C) {
            _iparm = new pastix_int_t[IPARM_SIZE];
            _dparm = new double[DPARM_SIZE];

            pstx<K>::initParam(_iparm, _dparm);
            Option& opt = *Option::get();
            int val = opt.val<int>("verbosity");
            if(val < 2)
                _iparm[IPARM_VERBOSE]         = API_VERBOSE_NOT;
            else
                _iparm[IPARM_VERBOSE]         = val - 1;
            _iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
            _iparm[IPARM_START_TASK]          = API_TASK_INIT;
            _iparm[IPARM_END_TASK]            = API_TASK_INIT;
            if(S == 'S') {
                _iparm[IPARM_SYM]             = API_SYM_YES;
                _iparm[IPARM_FACTORIZATION]   = opt.val<unsigned short>("master_not_spd", 0) ? API_FACT_LDLT : API_FACT_LLT;
            }
            else {
                _iparm[IPARM_SYM]             = API_SYM_NO;
                _iparm[IPARM_FACTORIZATION]   = API_FACT_LU;
                if(Wrapper<K>::is_complex)
                    _iparm[IPARM_TRANSPOSE_SOLVE] = API_YES;
            }
            _iparm[IPARM_RHSD_CHECK]          = API_NO;
            pastix_int_t* perm = new pastix_int_t[ncol];
            pstx<K>::dist(&_data, DMatrix::_communicator,
                          ncol, I, J, NULL, loc2glob,
                          perm, NULL, NULL, 1, _iparm, _dparm);

            _iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
            _iparm[IPARM_END_TASK]            = API_TASK_ANALYSE;

            pstx<K>::dist(&_data, DMatrix::_communicator,
                          ncol, I, J, NULL, loc2glob,
                          perm, NULL, NULL, 1, _iparm, _dparm);
            delete [] perm;

            _iparm[IPARM_VERBOSE]             = API_VERBOSE_NOT;

            _ncol2 = pstx<K>::getLocalNodeNbr(&_data);

            _loc2glob2 = new pastix_int_t[_ncol2];
            pstx<K>::getLocalNodeLst(&_data, _loc2glob2);

            pstx<K>::cscd_redispatch(ncol, I, J, C, NULL, 0, loc2glob,
                                     _ncol2, &_colptr2, &_rows2, &_values2, NULL, _loc2glob2,
                                     DMatrix::_communicator, 1);

            _iparm[IPARM_START_TASK]          = API_TASK_NUMFACT;
            _iparm[IPARM_END_TASK]            = API_TASK_NUMFACT;

            pstx<K>::dist(&_data, DMatrix::_communicator,
                          _ncol2, _colptr2, _rows2, _values2, _loc2glob2,
                          NULL, NULL, NULL, 1, _iparm, _dparm);

            _iparm[IPARM_CSCD_CORRECT]        = API_YES;
            delete [] I;
            delete [] loc2glob;
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
            K* rhs2 = new K[n * _ncol2];
            if(!DMatrix::_mapOwn && !DMatrix::_mapRecv) {
                int nloc = DMatrix::_ldistribution[DMatrix::_rank];
                DMatrix::initializeMap<1>(_ncol2, _loc2glob2, rhs2, rhs);
                DMatrix::_ldistribution = new int[1];
                *DMatrix::_ldistribution = nloc;
            }
            else
                DMatrix::redistribute<1>(rhs2, rhs);
            for(unsigned short nu = 1; nu < n; ++nu)
                DMatrix::redistribute<1>(rhs2 + nu * _ncol2, rhs + nu * *DMatrix::_ldistribution);

            _iparm[IPARM_START_TASK] = API_TASK_SOLVE;
            _iparm[IPARM_END_TASK]   = API_TASK_SOLVE;
            pstx<K>::dist(&_data, DMatrix::_communicator,
                          _ncol2, _colptr2, _rows2, _values2, _loc2glob2,
                          NULL, NULL, rhs2, 1, _iparm, _dparm);

            for(unsigned short nu = 0; nu < n; ++nu)
                DMatrix::redistribute<2>(rhs + nu * *DMatrix::_ldistribution, rhs2 + nu * _ncol2);
            delete [] rhs2;
        }
        void initialize() {
            DMatrix::initialize("PaStiX", { DISTRIBUTED_SOL_AND_RHS });
        }
};
#endif // DPASTIX

#ifdef PASTIXSUB
#define SUBDOMAIN HPDDM::PastixSub
template<class K>
class PastixSub {
    private:
        pastix_data_t*    _data;
        K*              _values;
        double*          _dparm;
        pastix_int_t      _ncol;
        pastix_int_t*   _colptr;
        pastix_int_t*     _rows;
        pastix_int_t*    _iparm;
    public:
        PastixSub() : _data(), _values(), _dparm(), _colptr(), _rows(), _iparm() { }
        PastixSub(const PastixSub&) = delete;
        ~PastixSub() {
            if(_iparm) {
                if(_iparm[IPARM_SYM] == API_SYM_YES || _iparm[IPARM_SYM] == API_SYM_HER) {
                    delete [] _rows;
                    delete [] _colptr;
                    delete [] _values;
                }
                _iparm[IPARM_START_TASK]          = API_TASK_CLEAN;
                _iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
                pstx<K>::seq(&_data, MPI_COMM_SELF,
                             0, NULL, NULL, NULL,
                             NULL, NULL, NULL, 1, _iparm, _dparm);
                delete [] _iparm;
                _iparm = nullptr;
                delete [] _dparm;
            }
        }
        static constexpr char _numbering = 'F';
        template<char N = HPDDM_NUMBERING>
        void numfact(MatrixCSR<K>* const& A, bool detection = false, K* const& schur = nullptr) {
            static_assert(N == 'C' || N == 'F', "Unknown numbering");
            if(!_iparm) {
                _iparm = new pastix_int_t[IPARM_SIZE];
                _dparm = new double[DPARM_SIZE];
                _ncol = A->_n;
                pstx<K>::initParam(_iparm, _dparm);
                _iparm[IPARM_VERBOSE]             = API_VERBOSE_NOT;
                _iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
                _iparm[IPARM_START_TASK]          = API_TASK_INIT;
                _iparm[IPARM_END_TASK]            = API_TASK_INIT;
                _iparm[IPARM_SCHUR]               = schur ? API_YES : API_NO;
                _iparm[IPARM_RHSD_CHECK]          = API_NO;
                _dparm[DPARM_EPSILON_MAGN_CTRL]   = -1.0 / HPDDM_PEN;
                if(A->_sym) {
                    _values = new K[A->_nnz];
                    _colptr = new int[_ncol + 1];
                    _rows = new int[A->_nnz];
                    _iparm[IPARM_SYM]             = API_SYM_YES;
                }
                else  {
                    _iparm[IPARM_SYM]             = API_SYM_NO;
                    _iparm[IPARM_FACTORIZATION]   = API_FACT_LU;
                }
            }
            if(A->_sym) {
                _iparm[IPARM_FACTORIZATION]       = Wrapper<K>::is_complex || detection ? API_FACT_LDLT : API_FACT_LLT;
                Wrapper<K>::template csrcsc<N, 'F'>(&_ncol, A->_a, A->_ja, A->_ia, _values, _rows, _colptr);
            }
            else {
                _values = A->_a;
                _colptr = A->_ia;
                _rows = A->_ja;
                if(N == 'C') {
                    std::for_each(_colptr, _colptr + _ncol + 1, [](int& i) { ++i; });
                    std::for_each(_rows, _rows + A->_nnz, [](int& i) { ++i; });
                }
            }
            pastix_int_t* perm = new pastix_int_t[2 * _ncol];
            pastix_int_t* iperm = perm + _ncol;
            int* listvar = nullptr;
            if(_iparm[IPARM_START_TASK] == API_TASK_INIT) {
                pstx<K>::seq(&_data, MPI_COMM_SELF,
                             _ncol, _colptr, _rows, NULL,
                             NULL, NULL, NULL, 1, _iparm, _dparm);
                if(schur != nullptr) {
                    listvar = new int[static_cast<int>(std::real(schur[0]))];
                    std::iota(listvar, listvar + static_cast<int>(std::real(schur[0])), static_cast<int>(std::real(schur[1])));
                    pstx<K>::setSchurUnknownList(_data, static_cast<int>(std::real(schur[0])), listvar);
                    pstx<K>::setSchurArray(_data, schur);
                }
                _iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
                _iparm[IPARM_END_TASK]            = API_TASK_NUMFACT;
            }
            else {
                _iparm[IPARM_START_TASK]          = API_TASK_NUMFACT;
                _iparm[IPARM_END_TASK]            = API_TASK_NUMFACT;
            }
            pstx<K>::seq(&_data, MPI_COMM_SELF,
                         _ncol, _colptr, _rows, _values,
                         perm, iperm, NULL, 1, _iparm, _dparm);
            delete [] listvar;
            delete [] perm;
            if(N == 'C' && _iparm[IPARM_SYM] == API_SYM_NO) {
                std::for_each(_colptr, _colptr + _ncol + 1, [](int& i) { --i; });
                std::for_each(_rows, _rows + A->_nnz, [](int& i) { --i; });
            }
        }
        void solve(K* const x, const unsigned short& n = 1) const {
            _iparm[IPARM_START_TASK] = API_TASK_SOLVE;
            _iparm[IPARM_END_TASK]   = API_TASK_SOLVE;
            pstx<K>::seq(const_cast<pastix_data_t**>(&_data), MPI_COMM_SELF,
                         // _ncol, _colptr, _rows, _values,
                         _ncol, NULL, NULL, NULL,
                         NULL, NULL, x, n, _iparm, _dparm);
        }
        void solve(const K* const b, K* const x, const unsigned short& n = 1) const {
            std::copy_n(b, n * _ncol, x);
            solve(x, n);
        }
};
#endif // PASTIXSUB
} // HPDDM
#endif // _HPDDM_PASTIX_
