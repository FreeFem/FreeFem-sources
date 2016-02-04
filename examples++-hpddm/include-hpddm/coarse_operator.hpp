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

#ifndef _HPDDM_COARSE_OPERATOR_
#define _HPDDM_COARSE_OPERATOR_

#if defined(DPASTIX) || defined(DMKL_PARDISO) || defined(DSUITESPARSE) || defined(DHYPRE)
# define HPDDM_CSR_CO
#endif
#if defined(DPASTIX) || defined(DMKL_PARDISO) || defined(DHYPRE)
# define HPDDM_LOC2GLOB
#endif
#if defined(DMKL_PARDISO) || defined(DHYPRE)
# define HPDDM_CONTIGUOUS
#endif

namespace HPDDM {
/* Class: Coarse operator
 *
 *  A class for handling coarse corrections.
 *
 * Template Parameters:
 *    Solver         - Solver used for the factorization of the coarse operator.
 *    S              - 'S'ymmetric or 'G'eneral coarse operator.
 *    K              - Scalar type. */
template<template<class> class Solver, char S, class K>
class CoarseOperator : public Solver<K> {
    private:
        /* Variable: gatherComm
         *  Communicator used for assembling right-hand sides. */
        MPI_Comm               _gatherComm;
        /* Variable: scatterComm
         *  Communicator used for distributing solution vectors. */
        MPI_Comm              _scatterComm;
        /* Variable: rankWorld
         *  Rank of the current subdomain in the global communicator supplied as an argument of <Coarse operator::constructionCommunicator>. */
        int                     _rankWorld;
        /* Variable: sizeWorld
         *  Size of <Subdomain::communicator>. */
        int                     _sizeWorld;
        int                     _sizeSplit;
        /* Variable: local
         *  Local number of coarse degrees of freedom (usually set to <Eigensolver::nu> after a call to <Eigensolver::selectNu>). */
        int                         _local;
        /* Variable: sizeRHS
         *  Local size of right-hand sides and solution vectors. */
        unsigned int              _sizeRHS;
        bool                       _offset;
        /* Function: constructionCommunicator
         *  Builds both <Coarse operator::scatterComm> and <DMatrix::communicator>. */
        template<bool exclude>
        void constructionCommunicator(const MPI_Comm&);
        /* Function: constructionCollective
         *
         *  Builds the buffers <DMatrix::gatherCounts>, <DMatrix::displs>, <DMatrix::gatherSplitCounts>, and <DMatrix::displsSplit> for all collective communications involving coarse corrections.
         *
         * Template Parameters:
         *    U              - True if the distribution of the coarse operator is uniform, false otherwise.
         *    D              - <DMatrix::Distribution> of right-hand sides and solution vectors.
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise. */
        template<bool U, typename Solver<K>::Distribution D, bool excluded>
        void constructionCollective(const unsigned short* = nullptr, unsigned short = 0, const unsigned short* = nullptr);
        /* Function: constructionMap
         *
         *  Builds the maps <DMatrix::ldistribution> and <DMatrix::idistribution> necessary for sending and receiving distributed right-hand sides or solution vectors.
         *
         * Template Parameters:
         *    T              - Coarse operator distribution topology.
         *    U              - True if the distribution of the coarse operator is uniform, false otherwise.
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise. */
        template<char T, bool U, bool excluded>
        void constructionMap(unsigned short, const unsigned short* = nullptr);
        /* Function: constructionMatrix
         *
         *  Builds and factorizes the coarse operator.
         *
         * Template Parameters:
         *    T              - Coarse operator distribution topology.
         *    U              - True if the distribution of the coarse operator is uniform, false otherwise.
         *    excluded       - True if the master processes are excluded from the domain decomposition, false otherwise.
         *    Operator       - Operator used in the definition of the Galerkin matrix. */
        template<char T, unsigned short U, unsigned short excluded, class Operator>
        std::pair<MPI_Request, const K*>* constructionMatrix(Operator&);
        /* Function: constructionCommunicatorCollective
         *
         *  Builds both communicators <Coarse operator::gatherComm> and <DMatrix::scatterComm> needed for coarse corrections.
         *
         * Template Parameter:
         *    countMasters   - True if the master processes must be taken into consideration, false otherwise. */
        template<bool countMasters>
        void constructionCommunicatorCollective(const unsigned short* const pt, unsigned short size, MPI_Comm& in, MPI_Comm* const out = nullptr) {
            unsigned short sizeComm = std::count_if(pt, pt + size, [](const unsigned short& nu) { return nu != 0; });
            if(sizeComm != size && in != MPI_COMM_NULL) {
                MPI_Group oldComm, newComm;
                MPI_Comm_group(in, &oldComm);
                if(*pt == 0)
                    ++sizeComm;
                int* array = new int[sizeComm];
                array[0] = 0;
                for(unsigned short i = 1, j = 1, k = 0; j < sizeComm; ++i) {
                    if(pt[i] != 0)
                        array[j++] = i - k;
                    else if(countMasters && Solver<K>::_ldistribution[k + 1] == i)
                        ++k;
                }
                MPI_Group_incl(oldComm, sizeComm, array, &newComm);
                if(out)
                    MPI_Comm_create(in, newComm, out);
                else {
                    MPI_Comm tmp;
                    MPI_Comm_create(in, newComm, &tmp);
                    MPI_Comm_free(&in);
                    if(tmp != MPI_COMM_NULL) {
                        MPI_Comm_dup(tmp, &in);
                        MPI_Comm_free(&tmp);
                    }
                    else
                        in = MPI_COMM_NULL;
                }
                delete [] array;
            }
            else if(out)
                MPI_Comm_dup(in, out);
        }
        /* Function: transfer
         *
         *  Transfers vectors from the fine grid to the coarse grid, and vice versa.
         *
         * Template Parameter:
         *    T              - True if fine to coarse, false otherwise.
         *
         * Parameters:
         *    counts         - Array of integers <DMatrix::gatherSplitCounts> or <DMatrix::gatherCounts> used for MPI collectives.
         *    n              - Number of vectors or size of the communicator <Coarse operator::gatherComm> for MPI collectives.
         *    m              - Size of the communicator <Coarse operator::gatherComm> for MPI collectives or number of vectors.
         *    ab             - Array to transfer. */
        template<bool T>
        void transfer(int* const counts, const int n, const int m, K* const ab) const {
            if(!T) {
                std::for_each(counts, counts + 2 * n, [&](int& i) { i *= m; });
                MPI_Gatherv(MPI_IN_PLACE, 0, Wrapper<K>::mpi_type(), ab, counts, counts + n, Wrapper<K>::mpi_type(), 0, _gatherComm);
            }
            if(n != 1 && m != 1) {
                int size = T ? m : n;
                K* ba = new K[counts[size - 1] + counts[2 * size - 1]];
                if(!T) {
                    for(int i = 0; i < size; ++i)
                        for(int j = 0; j < m; ++j)
                            std::copy_n(ab + counts[size + i] + j * (counts[i] / m), counts[i] / m, ba + counts[size + i] / m + j * ((counts[size - 1] + counts[2 * size - 1]) / m));
                }
                else {
                    for(int j = 0; j < n; ++j)
                        for(int i = 0; i < size; ++i)
                            std::copy_n(ab + counts[size + i] / n + j * ((counts[size - 1] + counts[2 * size - 1]) / n), counts[i] / n, ba + counts[size + i] + j * (counts[i] / n));
                }
                std::copy_n(ba, counts[size - 1] + counts[2 * size - 1], ab);
                delete [] ba;
            }
            if(T) {
                MPI_Scatterv(ab, counts, counts + m, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, Wrapper<K>::mpi_type(), 0, _gatherComm);
                std::for_each(counts, counts + 2 * m, [&](int& i) { i /= n; });
            }
        }
    public:
        CoarseOperator() : _gatherComm(MPI_COMM_NULL), _scatterComm(MPI_COMM_NULL), _rankWorld(), _sizeWorld(), _sizeSplit(), _local(), _sizeRHS(), _offset(false) {
            static_assert(S == 'S' || S == 'G', "Unknown symmetry");
            static_assert(!Wrapper<K>::is_complex || S != 'S', "Symmetric complex coarse operators are not supported");
        }
        ~CoarseOperator() {
            if(_gatherComm != _scatterComm && _gatherComm != MPI_COMM_NULL)
                MPI_Comm_free(&_gatherComm);
            if(_scatterComm != MPI_COMM_NULL)
                MPI_Comm_free(&_scatterComm);
        }
        /* Function: construction
         *  Wrapper function to call all needed subroutines. */
        template<unsigned short U, unsigned short excluded, class Operator>
        std::pair<MPI_Request, const K*>* construction(Operator&&, const MPI_Comm&);
        /* Function: callSolver
         *
         *  Solves a coarse system.
         *
         * Parameter:
         *    rhs            - Input right-hand side, solution vector is stored in-place. */
        template<bool = false>
        void callSolver(K* const, const unsigned short& = 1);
#if HPDDM_ICOLLECTIVE
        template<bool = false>
        void IcallSolver(K* const, MPI_Request*);
#endif
        /* Function: getRank
         *  Simple accessor that returns <Coarse operator::rankWorld>. */
        int getRank() const { return _rankWorld; }
        /* Function: getLocal
         *  Returns the value of <Coarse operator::local>. */
        int getLocal() const { return _local; }
        /* Function: getAddrLocal
         *  Returns the address of <Coarse operator::local>. */
        const int* getAddrLocal() const { return &_local; }
        /* Function: setLocal
         *  Sets the value of <Coarse operator::local>. */
        void setLocal(int l) { _local = l; }
        /* Function: getSizeRHS
         *  Returns the value of <Coarse operator::sizeRHS>. */
        unsigned int getSizeRHS() const { return _sizeRHS; }
};
} // HPDDM
#endif // _HPDDM_COARSE_OPERATOR_
