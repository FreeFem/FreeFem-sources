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
        void constructionCollective(const unsigned short* = nullptr, unsigned short p = 0, const unsigned short* = nullptr);
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
        void callSolver(K* const, const int& = 0);
#if HPDDM_ICOLLECTIVE
        template<bool = false>
        void IcallSolver(K* const, MPI_Request*, const int& = 0);
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
        /* Function: reallocateRHS
         *
         *  Reallocates the array for storing right-hand sides and solution vectors.
         *
         * Parameters:
         *    rhs            - Reference to the pointer to reallocate.
         *    n              - Additional space needed, see also <Coarse operator::sizeRHS>. */
        void reallocateRHS(K*& rhs, const unsigned short& n) const {
            if(rhs)
                delete [] rhs;
            if(Solver<K>::_communicator != MPI_COMM_NULL)
                rhs = new K[_sizeRHS + _sizeSplit * n];
            else
                rhs = new K[_sizeRHS + n];
        }
};
} // HPDDM
#endif // _HPDDM_COARSE_OPERATOR_
