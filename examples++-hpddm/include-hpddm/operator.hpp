/*
   This file is part of HPDDM.

   Author(s): Pierre Jolivet <pierre.jolivet@enseeiht.fr>
        Date: 2013-03-12

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

#ifndef _HPDDM_OPERATOR_
#define _HPDDM_OPERATOR_

#include <queue>

namespace HPDDM {
template<bool> class Members { };
template<> class Members<true> {
    protected:
        std::unordered_map<unsigned short, unsigned int> _offsets;
        std::vector<std::vector<unsigned short>>     _vecSparsity;
        const unsigned short                                _rank;
        unsigned short                               _consolidate;
        Members(unsigned short r) : _rank(r), _consolidate() { }
};
template<char P, class Preconditioner, class K>
class OperatorBase : protected Members<P != 's'> {
    private:
        template<class T>
        class has_LDR {
            private:
                typedef char one;
                typedef one (&two)[2];
                template<class C> static one test(decltype(&C::getLDR));
                template<class C> static two test(...);
            public:
                static constexpr bool value = (sizeof(test<T>(0)) == sizeof(one));
        };
        template<class Q = Preconditioner> typename std::enable_if<has_LDR<typename std::remove_reference<Q>::type>::value, bool>::type
        offsetDeflation() {
            const unsigned int offset = *_p.getLDR() - _n;
            if(_deflation && offset)
                std::for_each(_deflation, _deflation + _local, [&](K*& v) { v -= offset; });
            return true;
        }
        template<class Q = Preconditioner> typename std::enable_if<!has_LDR<typename std::remove_reference<Q>::type>::value, bool>::type
        offsetDeflation() { return false; }
    protected:
        const Preconditioner&                 _p;
        K** const                     _deflation;
        const vectorNeighbor&               _map;
        std::vector<unsigned short>    _sparsity;
        const int                             _n;
        const int                         _local;
        unsigned short                   _signed;
        unsigned short             _connectivity;
        template<char Q = P, typename std::enable_if<Q == 's'>::type* = nullptr>
        OperatorBase(const Preconditioner& p, const unsigned short& c) : _p(p), _deflation(p.getVectors()), _map(p.getMap()), _n(p.getDof()), _local(p.getLocal()), _connectivity(c) {
            static_assert(Q == P, "Wrong sparsity pattern");
            _sparsity.reserve(_map.size());
            for(const pairNeighbor& n : _map)
                _sparsity.emplace_back(n.first);
        }
        template<char Q = P, typename std::enable_if<Q != 's'>::type* = nullptr>
        OperatorBase(const Preconditioner& p, const unsigned short& c) : Members<true>(p.getRank()), _p(p), _deflation(p.getVectors()), _map(p.getMap()), _n(p.getDof()), _local(p.getLocal()), _signed(_p.getSigned()), _connectivity(c) {
            const unsigned int offset = *_p.getLDR() - _n;
            if(_deflation && offset)
                std::for_each(_deflation, _deflation + _local, [&](K*& v) { v += offset; });
            static_assert(Q == P, "Wrong sparsity pattern");
            if(!_map.empty()) {
                unsigned short** recvSparsity = new unsigned short*[_map.size() + 1];
                *recvSparsity = new unsigned short[(_connectivity + 1) * _map.size()];
                unsigned short* sendSparsity = *recvSparsity + _connectivity * _map.size();
                MPI_Request* rq = _p.getRq();
                for(unsigned short i = 0; i < _map.size(); ++i) {
                    sendSparsity[i] = _map[i].first;
                    recvSparsity[i] = *recvSparsity + _connectivity * i;
                    MPI_Irecv(recvSparsity[i], _connectivity, MPI_UNSIGNED_SHORT, _map[i].first, 4, _p.getCommunicator(), rq + i);
                }
                for(unsigned short i = 0; i < _map.size(); ++i)
                    MPI_Isend(sendSparsity, _map.size(), MPI_UNSIGNED_SHORT, _map[i].first, 4, _p.getCommunicator(), rq + _map.size() + i);
                Members<true>::_vecSparsity.resize(_map.size());
                for(unsigned short i = 0; i < _map.size(); ++i) {
                    int index, count;
                    MPI_Status status;
                    MPI_Waitany(_map.size(), rq, &index, &status);
                    MPI_Get_count(&status, MPI_UNSIGNED_SHORT, &count);
                    Members<true>::_vecSparsity[index].assign(recvSparsity[index], recvSparsity[index] + count);
                }
                MPI_Waitall(_map.size(), rq + _map.size(), MPI_STATUSES_IGNORE);

                delete [] *recvSparsity;
                delete [] recvSparsity;

                _sparsity.reserve(_map.size());
                if(P == 'c') {
                    std::vector<unsigned short> neighbors;
                    neighbors.reserve(_map.size());
                    std::for_each(_map.cbegin(), _map.cend(), [&](const pairNeighbor& n) { neighbors.emplace_back(n.first); });
                    typedef std::pair<std::vector<unsigned short>::const_iterator, std::vector<unsigned short>::const_iterator> pairIt;
                    auto comp = [](const pairIt& lhs, const pairIt& rhs) { return *lhs.first > *rhs.first; };
                    std::priority_queue<pairIt, std::vector<pairIt>, decltype(comp)> pq(comp);
                    pq.push({ neighbors.cbegin(), neighbors.cend() });
                    for(const std::vector<unsigned short>& v : Members<true>::_vecSparsity)
                        pq.push({ v.cbegin(), v.cend() });
                    while(!pq.empty()) {
                        pairIt p = pq.top();
                        pq.pop();
                        if(*p.first != Members<true>::_rank && (_sparsity.empty() || (*p.first != _sparsity.back())))
                            _sparsity.emplace_back(*p.first);
                        if(++p.first != p.second)
                            pq.push(p);
                    }
                }
                else {
                    for(const pairNeighbor& n : _map)
                        _sparsity.emplace_back(n.first);
                    for(std::vector<unsigned short>& v : Members<true>::_vecSparsity) {
                        unsigned short i = 0, j = 0, k = 0;
                        while(i < v.size() && j < _sparsity.size()) {
                            if(v[i] == Members<true>::_rank) {
                                v[k++] = Members<true>::_rank;
                                ++i;
                            }
                            else if(v[i] < _sparsity[j])
                                ++i;
                            else if(v[i] > _sparsity[j])
                                ++j;
                            else {
                                v[k++] = _sparsity[j++];
                                ++i;
                            }
                        }
                        v.resize(k);
                    }
                }
            }
        }
        ~OperatorBase() { offsetDeflation(); }
        template<char S, bool U, class T>
        void initialize(T& in, const unsigned short* info, T const& out, MPI_Request* const& rqRecv, unsigned short*& infoNeighbor) {
            static_assert(P == 'c' || P == 'f', "Unsupported constructor with such a sparsity pattern");
            if(!U) {
                if(P == 'c') {
                    infoNeighbor = new unsigned short[_map.size()];
                    std::vector<unsigned short>::const_iterator begin = _sparsity.cbegin();
                    for(unsigned short i = 0; i < _map.size(); ++i) {
                        std::vector<unsigned short>::const_iterator idx = std::lower_bound(begin, _sparsity.cend(), _map[i].first);
                        infoNeighbor[i] = info[std::distance(_sparsity.cbegin(), idx)];
                        begin = idx + 1;
                    }
                }
                else
                    infoNeighbor = const_cast<unsigned short*>(info);
            }
            std::vector<unsigned int> displs;
            displs.reserve(2 * _map.size());
            if(S != 'S') {
                if(!U) {
                    unsigned short size = std::accumulate(infoNeighbor, infoNeighbor + _map.size(), _local);
                    if(!_map.empty())
                        displs.emplace_back(size * _map[0].second.size());
                    for(unsigned short i = 1; i < _map.size(); ++i)
                        displs.emplace_back(displs.back() + size * _map[i].second.size());
                    for(unsigned short i = 0; i < _map.size(); ++i) {
                        size = infoNeighbor[i];
                        std::vector<unsigned short>::const_iterator begin = _sparsity.cbegin();
                        for(const unsigned short& rank : Members<true>::_vecSparsity[i]) {
                            if(rank == Members<true>::_rank)
                                size += _local;
                            else {
                                std::vector<unsigned short>::const_iterator idx = std::lower_bound(begin, _sparsity.cend(), rank);
                                size += info[std::distance(_sparsity.cbegin(), idx)];
                                begin = idx + 1;
                            }
                        }
                        if(_local)
                            displs.emplace_back(displs.back() + size * _map[i].second.size());
                        else
                            rqRecv[i] = MPI_REQUEST_NULL;
                    }
                }
                else {
                    if(!_map.empty())
                        displs.emplace_back(_local * (_map.size() + 1) * _map[0].second.size());
                    for(unsigned short i = 1; i < _map.size(); ++i)
                        displs.emplace_back(displs.back() + _local * (_map.size() + 1) * _map[i].second.size());
                    for(unsigned short i = 0; i < _map.size(); ++i)
                        displs.emplace_back(displs.back() + _local * (Members<true>::_vecSparsity[i].size() + 1) * _map[i].second.size());
                }
            }
            else {
                if(!U) {
                    unsigned short size = std::accumulate(infoNeighbor, infoNeighbor + _map.size(), 0);
                    if(!_map.empty()) {
                        displs.emplace_back((size + _local * (0 < _signed)) * _map[0].second.size());
                        size -= infoNeighbor[0];
                    }
                    for(unsigned short i = 1; i < _map.size(); ++i) {
                        displs.emplace_back(displs.back() + (size + _local * (i < _signed)) * _map[i].second.size());
                        size -= infoNeighbor[i];
                    }
                    for(unsigned short i = 0; i < _map.size(); ++i) {
                        size = infoNeighbor[i] * !(i < _signed) + _local;
                        std::vector<unsigned short>::const_iterator end = _sparsity.cend();
                        for(std::vector<unsigned short>::const_reverse_iterator rit = Members<true>::_vecSparsity[i].rbegin(); *rit > Members<true>::_rank; ++rit) {
                            std::vector<unsigned short>::const_iterator idx = std::lower_bound(_sparsity.cbegin(), end, *rit);
                            size += info[std::distance(_sparsity.cbegin(), idx)];
                            end = idx - 1;
                        }
                        if(_local)
                            displs.emplace_back(displs.back() + size * _map[i].second.size());
                        else
                            rqRecv[i] = MPI_REQUEST_NULL;
                    }
                }
                else {
                    if(!_map.empty())
                        displs.emplace_back(_local * (_map.size() + (0 < _signed)) * _map[0].second.size());
                    for(unsigned short i = 1; i < _map.size(); ++i)
                        displs.emplace_back(displs.back() + _local * (_map.size() + (i < _signed) - i) * _map[i].second.size());
                    for(unsigned short i = 0; i < _map.size(); ++i) {
                        unsigned short size = std::distance(std::lower_bound(Members<true>::_vecSparsity[i].cbegin(), Members<true>::_vecSparsity[i].cend(), Members<true>::_rank), Members<true>::_vecSparsity[i].cend()) + !(i < _signed);
                        displs.emplace_back(displs.back() + _local * size * _map[i].second.size());
                    }
                }
            }
            if(!displs.empty()) {
                *in = new K[displs.back()];
                for(unsigned short i = 1; i < _map.size(); ++i)
                    in[i] = *in + displs[i - 1];
                if(U == 1 || _local)
                    for(unsigned short i = 0; i < _map.size(); ++i) {
                        if(displs[i + _map.size()] != displs[i - 1 + _map.size()]) {
                            out[i] = *in + displs[i - 1 + _map.size()];
                            MPI_Irecv(out[i], displs[i + _map.size()] - displs[i - 1 + _map.size()], Wrapper<K>::mpi_type(), _map[i].first, 2, _p.getCommunicator(), rqRecv + i);
                        }
                        else
                            out[i] = nullptr;
                    }
            }
            else
                *in = nullptr;
        }
        template<char S, char N, bool U, char Q = P, typename std::enable_if<Q != 's'>::type* = nullptr>
        void assembleOperator(int* I, int* J, int coefficients, unsigned int offsetI, unsigned int* offsetJ, unsigned short* const& infoNeighbor) {
            if(Members<true>::_consolidate == _map.size()) {
                unsigned short between = std::distance(_sparsity.cbegin(), std::lower_bound(_sparsity.cbegin(), _sparsity.cend(), _p.getRank()));
                unsigned int offset = 0;
                if(S != 'S')
                    for(unsigned short k = 0; k < between; ++k)
                        for(unsigned short i = 0; i < _local; ++i) {
                            unsigned int l = offset + coefficients * i;
                            for(unsigned short j = 0; j < (U ? _local : infoNeighbor[k]); ++j) {
#ifndef HPDDM_CSR_CO
                                I[l + j] = offsetI + i;
#endif
                                J[l + j] = (U ? _sparsity[k] * _local + (N == 'F') : offsetJ[k]) + j;
                            }
                            offset += U ? _local : infoNeighbor[k];
                        }
                else
                    coefficients += _local - 1;
                for(unsigned short i = 0; i < _local; ++i) {
                    unsigned int l = offset + coefficients * i - (S == 'S') * ((i * (i - 1)) / 2);
                    for(unsigned short j = (S == 'S') * i; j < _local; ++j) {
#ifndef HPDDM_CSR_CO
                        I[l + j] = offsetI + i;
#endif
                        J[l + j] = offsetI + j;
                    }
                }
                offset += _local;
                for(unsigned short k = between; k < _sparsity.size(); ++k) {
                    for(unsigned short i = 0; i < _local; ++i) {
                        unsigned int l = offset + coefficients * i - (S == 'S') * ((i * (i - 1)) / 2);
                        for(unsigned short j = 0; j < (U ? _local : infoNeighbor[k]); ++j) {
#ifndef HPDDM_CSR_CO
                            I[l + j] = offsetI + i;
#endif
                            J[l + j] = (U ? _sparsity[k] * _local + (N == 'F') : offsetJ[k - (S == 'S') * between]) + j;
                        }
                    }
                    offset += U ? _local : infoNeighbor[k];
                }
            }
        }
    public:
        static constexpr char _pattern = P != 's' ? 'c' : 's';
        void adjustConnectivity(const MPI_Comm& comm) {
            if(P == 'c') {
#if 0
                _connectivity *= _connectivity - 1;
#else
                _connectivity = _sparsity.size();
                MPI_Allreduce(MPI_IN_PLACE, &_connectivity, 1, MPI_UNSIGNED_SHORT, MPI_MAX, comm);
#endif
            }
        }
        const std::vector<unsigned short>& getPattern() const { return _sparsity; }
        unsigned short getConnectivity() const { return _connectivity; }
        template<char Q = P, typename std::enable_if<Q != 's'>::type* = nullptr>
        void initialize(unsigned int, K*&, unsigned short) { }
};

#if HPDDM_SCHWARZ
template<class Preconditioner, class K>
class MatrixMultiplication : public OperatorBase<'s', Preconditioner, K> {
    private:
        typedef OperatorBase<'s', Preconditioner, K> super;
        const MatrixCSR<K>* const                       _A;
        MatrixCSR<K>*                                   _C;
        const underlying_type<K>* const                 _D;
        K*                                           _work;
        template<char S, bool U>
        void applyFromNeighbor(const K* in, unsigned short index, K*& work, unsigned short* infoNeighbor) {
            int m = U ? super::_local : *infoNeighbor;
            std::fill_n(work, m * super::_n, K());
            for(unsigned short i = 0; i < m; ++i)
                Wrapper<K>::sctr(super::_map[index].second.size(), in + i * super::_map[index].second.size(), super::_map[index].second.data(), work + i * super::_n);
            Wrapper<K>::diag(super::_n, _D, work, _work, m);
            Blas<K>::gemm(&(Wrapper<K>::transc), "N", &(super::_local), &m, &(super::_n), &(Wrapper<K>::d__1), *super::_deflation, &(super::_n), _work, &(super::_n), &(Wrapper<K>::d__0), work, &(super::_local));
        }
    public:
        template<template<class> class Solver, char S, class T> friend class CoarseOperator;
        MatrixMultiplication(const Preconditioner& p, const unsigned short c) : super(p, std::move(c)), _A(p.getMatrix()), _C(), _D(p.getScaling()) { }
        void initialize(unsigned int k, K*& work, unsigned short s) {
            if(_A->_sym) {
                std::vector<std::vector<std::pair<unsigned int, K>>> v(_A->_n);
                unsigned int nnz = std::floor((_A->_nnz + _A->_n - 1) / _A->_n) * 2;
                std::for_each(v.begin(), v.end(), [&](std::vector<std::pair<unsigned int, K>>& r) { r.reserve(nnz); });
                nnz = 0;
                for(unsigned int i = 0; i < _A->_n; ++i) {
                    const underlying_type<K> scal = _D[i];
                    unsigned int j = _A->_ia[i] - (HPDDM_NUMBERING == 'F');
                    while(j < _A->_ia[i + 1] - (HPDDM_NUMBERING == 'F' ? 2 : 1)) {
                        if(_D[_A->_ja[j] - (HPDDM_NUMBERING == 'F')] > HPDDM_EPS) {
                            v[i].emplace_back(_A->_ja[j], _A->_a[j] * _D[_A->_ja[j] - (HPDDM_NUMBERING == 'F')]);
                            ++nnz;
                        }
                        if(scal > HPDDM_EPS) {
                            v[_A->_ja[j] - (HPDDM_NUMBERING == 'F')].emplace_back(i + (HPDDM_NUMBERING == 'F'), _A->_a[j] * scal);
                            ++nnz;
                        }
                        ++j;
                    }
                    if(i != _A->_ja[j] - (HPDDM_NUMBERING == 'F')) {
                        if(_D[_A->_ja[j] - (HPDDM_NUMBERING == 'F')] > HPDDM_EPS) {
                            v[i].emplace_back(_A->_ja[j], _A->_a[j] * _D[_A->_ja[j] - (HPDDM_NUMBERING == 'F')]);
                            ++nnz;
                        }
                    }
                    if(scal > HPDDM_EPS) {
                        v[_A->_ja[j] - (HPDDM_NUMBERING == 'F')].emplace_back(i + (HPDDM_NUMBERING == 'F'), _A->_a[j] * scal);
                        ++nnz;
                    }
                }
                _C = new MatrixCSR<K>(_A->_n, _A->_n, nnz, false);
                _C->_ia[0] = (Wrapper<K>::I == 'F');
                nnz = 0;
#pragma omp parallel for schedule(static, HPDDM_GRANULARITY)
                for(unsigned int i = 0; i < _A->_n; ++i)
                    std::sort(v[i].begin(), v[i].end(), [](const std::pair<unsigned int, K>& lhs, const std::pair<unsigned int, K>& rhs) { return lhs.first < rhs.first; });
                for(unsigned int i = 0; i < _A->_n; ++i) {
                    for(const std::pair<unsigned int, K>& p : v[i]) {
                        _C->_ja[nnz] = p.first + (Wrapper<K>::I == 'F' && HPDDM_NUMBERING != Wrapper<K>::I);
                        _C->_a[nnz++] = p.second;
                    }
                    _C->_ia[i + 1] = nnz + (Wrapper<K>::I == 'F');
                }
            }
            else {
                _C = new MatrixCSR<K>(_A->_n, _A->_n, _A->_nnz, false);
                _C->_ia[0] = (Wrapper<K>::I == 'F');
                unsigned int nnz = 0;
                for(unsigned int i = 0; i < _A->_n; ++i) {
                    for(unsigned int j = _A->_ia[i] - (HPDDM_NUMBERING == 'F'); j < _A->_ia[i + 1] - (HPDDM_NUMBERING == 'F'); ++j)
                        if(_D[_A->_ja[j] - (HPDDM_NUMBERING == 'F')] > HPDDM_EPS) {
                            _C->_ja[nnz] = _A->_ja[j] + (Wrapper<K>::I == 'F' && HPDDM_NUMBERING != Wrapper<K>::I);
                            _C->_a[nnz++] = _A->_a[j] * _D[_A->_ja[j] - (HPDDM_NUMBERING == 'F')];
                        }
                    _C->_ia[i + 1] = nnz + (Wrapper<K>::I == 'F');
                }
                _C->_nnz = nnz;
            }
            work = new K[2 * k];
            _work = work + k;
            super::_signed = s;
        }
        template<char S, bool U, class T>
        void applyToNeighbor(T& in, K*& work, MPI_Request*& rq, const unsigned short* info, T = nullptr, MPI_Request* = nullptr) {
            Wrapper<K>::template csrmm<Wrapper<K>::I>(false, &(super::_n), &(super::_local), _C->_a, _C->_ia, _C->_ja, *super::_deflation, _work);
            delete _C;
            for(unsigned short i = 0; i < super::_signed; ++i) {
                if(U || info[i]) {
                    for(unsigned short j = 0; j < super::_local; ++j)
                        Wrapper<K>::gthr(super::_map[i].second.size(), _work + j * super::_n, in[i] + j * super::_map[i].second.size(), super::_map[i].second.data());
                    MPI_Isend(in[i], super::_map[i].second.size() * super::_local, Wrapper<K>::mpi_type(), super::_map[i].first, 2, super::_p.getCommunicator(), rq++);
                }
            }
            Wrapper<K>::diag(super::_n, _D, _work, work, super::_local);
        }
        template<char S, bool U>
        void assembleForMaster(K* C, const K* in, const int& coefficients, unsigned short index, K* arrayC, unsigned short* const& infoNeighbor = nullptr) {
            applyFromNeighbor<S, U>(in, index, arrayC, infoNeighbor);
            for(unsigned short j = 0; j < (U ? super::_local : *infoNeighbor); ++j) {
                K* pt = C + j;
                for(unsigned short i = 0; i < super::_local; pt += coefficients - (S == 'S') * i++)
                    *pt = arrayC[j * super::_local + i];
            }
        }
        template<char S, char N, bool U>
        void applyFromNeighborMaster(const K* in, unsigned short index, int* I, int* J, K* C, int coefficients, unsigned int offsetI, unsigned int* offsetJ, K* arrayC, unsigned short* const& infoNeighbor = nullptr) {
            applyFromNeighbor<S, U>(in, index, arrayC, infoNeighbor);
            unsigned int offset = U ? super::_map[index].first * super::_local + (N == 'F') : *offsetJ;
            for(unsigned short i = 0; i < super::_local; ++i) {
                unsigned int l = coefficients * i - (S == 'S') * (i * (i - 1)) / 2;
                for(unsigned short j = 0; j < (U ? super::_local : *infoNeighbor); ++j) {
#ifndef HPDDM_CSR_CO
                    I[l + j] = offsetI + i;
#endif
                    J[l + j] = offset + j;
                    C[l + j] = arrayC[j * super::_local + i];
                }
            }
        }
};
#endif // HPDDM_SCHWARZ

#if HPDDM_FETI
template<class Preconditioner, FetiPrcndtnr Q, class K>
class FetiProjection : public OperatorBase<Q == FetiPrcndtnr::SUPERLUMPED ? 'f' : 'c', Preconditioner, K> {
    private:
        typedef OperatorBase<Q == FetiPrcndtnr::SUPERLUMPED ? 'f' : 'c', Preconditioner, K> super;
        template<char S, bool U>
        void applyFromNeighbor(const K* in, unsigned short index, K*& work, unsigned short* info) {
            std::vector<unsigned short>::const_iterator middle = std::lower_bound(super::_vecSparsity[index].cbegin(), super::_vecSparsity[index].cend(), super::_rank);
            unsigned int accumulate = 0;
            if(!(index < super::_signed)) {
                for(unsigned short k = 0; k < (U ? super::_local : info[std::distance(super::_sparsity.cbegin(), std::lower_bound(super::_sparsity.cbegin(), super::_sparsity.cend(), super::_map[index].first))]); ++k)
                    for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                        work[super::_offsets[super::_map[index].first] + super::_map[index].second[j] + k * super::_n] += in[k * super::_map[index].second.size() + j];
                accumulate += (U ? super::_local : info[std::distance(super::_sparsity.cbegin(), std::lower_bound(super::_sparsity.cbegin(), super::_sparsity.cend(), super::_map[index].first))]) * super::_map[index].second.size();
            }
            else if(S != 'S') {
                for(unsigned short k = 0; k < (U ? super::_local : info[std::distance(super::_sparsity.cbegin(), std::lower_bound(super::_sparsity.cbegin(), super::_sparsity.cend(), super::_map[index].first))]); ++k)
                    for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                        work[super::_offsets[super::_map[index].first] + super::_map[index].second[j] + k * super::_n] -= in[k * super::_map[index].second.size() + j];
                accumulate += (U ? super::_local : info[std::distance(super::_sparsity.cbegin(), std::lower_bound(super::_sparsity.cbegin(), super::_sparsity.cend(), super::_map[index].first))]) * super::_map[index].second.size();
            }
            std::vector<unsigned short>::const_iterator begin = super::_sparsity.cbegin();
            if(S != 'S')
                for(std::vector<unsigned short>::const_iterator it = super::_vecSparsity[index].cbegin(); it != middle; ++it) {
                    if(!U) {
                        std::vector<unsigned short>::const_iterator idx = std::lower_bound(begin, super::_sparsity.cend(), *it);
                        if(*it > super::_map[index].first || super::_signed > index)
                            for(unsigned short k = 0; k < info[std::distance(super::_sparsity.cbegin(), idx)]; ++k)
                                for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                    work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] -= in[accumulate + k * super::_map[index].second.size() + j];
                        else
                            for(unsigned short k = 0; k < info[std::distance(super::_sparsity.cbegin(), idx)]; ++k)
                                for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                    work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
                        accumulate += info[std::distance(super::_sparsity.cbegin(), idx)] * super::_map[index].second.size();
                        begin = idx + 1;
                    }
                    else {
                        if(*it > super::_map[index].first || super::_signed > index)
                            for(unsigned short k = 0; k < super::_local; ++k)
                                for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                    work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] -= in[accumulate + k * super::_map[index].second.size() + j];
                        else
                            for(unsigned short k = 0; k < super::_local; ++k)
                                for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                    work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
                        accumulate += super::_local * super::_map[index].second.size();
                    }
                }
            if(index < super::_signed)
                for(unsigned short k = 0; k < super::_local; ++k)
                    for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                        work[super::_offsets[super::_rank] + super::_map[index].second[j] + k * super::_n] -= in[accumulate + k * super::_map[index].second.size() + j];
            else
                for(unsigned short k = 0; k < super::_local; ++k)
                    for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                        work[super::_offsets[super::_rank] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
            accumulate += super::_local * super::_map[index].second.size();
            for(std::vector<unsigned short>::const_iterator it = middle + 1; it < super::_vecSparsity[index].cend(); ++it) {
                if(!U) {
                    std::vector<unsigned short>::const_iterator idx = std::lower_bound(begin, super::_sparsity.cend(), *it);
                    if(*it > super::_map[index].first && super::_signed > index)
                        for(unsigned short k = 0; k < info[std::distance(super::_sparsity.cbegin(), idx)]; ++k)
                            for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] -= in[accumulate + k * super::_map[index].second.size() + j];
                    else
                        for(unsigned short k = 0; k < info[std::distance(super::_sparsity.cbegin(), idx)]; ++k)
                            for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
                    accumulate += info[std::distance(super::_sparsity.cbegin(), idx)] * super::_map[index].second.size();
                    begin = idx + 1;
                }
                else {
                    if(*it > super::_map[index].first && super::_signed > index)
                        for(unsigned short k = 0; k < super::_local; ++k)
                            for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] -= in[accumulate + k * super::_map[index].second.size() + j];
                    else
                        for(unsigned short k = 0; k < super::_local; ++k)
                            for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
                    accumulate += super::_local * super::_map[index].second.size();
                }
            }
        }
    public:
        template<template<class> class Solver, char S, class T> friend class CoarseOperator;
        FetiProjection(const Preconditioner& p, const unsigned short c) : super(p, std::move(c)) { }
        template<char S, bool U, class T>
        void applyToNeighbor(T& in, K*& work, MPI_Request*& rq, const unsigned short* info, T const& out = nullptr, MPI_Request* const& rqRecv = nullptr) {
            unsigned short* infoNeighbor;
            super::template initialize<S, U>(in, info, out, rqRecv, infoNeighbor);
            MPI_Request* rqMult = new MPI_Request[2 * super::_map.size()];
            unsigned int* offset = new unsigned int[super::_map.size() + 2];
            offset[0] = 0;
            offset[1] = super::_local;
            for(unsigned short i = 2; i < super::_map.size() + 2; ++i)
                offset[i] = offset[i - 1] + (U ? super::_local : infoNeighbor[i - 2]);
            const int nbMult = super::_p.getMult();
            K* mult = new K[offset[super::_map.size() + 1] * nbMult];
            unsigned short* displs = new unsigned short[super::_map.size() + 1];
            displs[0] = 0;
            for(unsigned short i = 0; i < super::_map.size(); ++i) {
                MPI_Irecv(mult + offset[i + 1] * nbMult + displs[i] * (U ? super::_local : infoNeighbor[i]), super::_map[i].second.size() * (U ? super::_local : infoNeighbor[i]), Wrapper<K>::mpi_type(), super::_map[i].first, 11, super::_p.getCommunicator(), rqMult + i);
                displs[i + 1] = displs[i] + super::_map[i].second.size();
            }

            K* tmp = new K[offset[super::_map.size() + 1] * super::_n]();
            const underlying_type<K>* const* const m = super::_p.getScaling();
            for(unsigned short i = 0; i < super::_signed; ++i) {
                for(unsigned short k = 0; k < super::_local; ++k)
                    for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                        tmp[super::_map[i].second[j] + k * super::_n] -= m[i][j] * (mult[displs[i] * super::_local + j + k * super::_map[i].second.size()] = - super::_deflation[k][super::_map[i].second[j]]);
                MPI_Isend(mult + displs[i] * super::_local, super::_map[i].second.size() * super::_local, Wrapper<K>::mpi_type(), super::_map[i].first, 11, super::_p.getCommunicator(), rqMult + super::_map.size() + i);
            }
            for(unsigned short i = super::_signed; i < super::_map.size(); ++i) {
                for(unsigned short k = 0; k < super::_local; ++k)
                    for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                        tmp[super::_map[i].second[j] + k * super::_n] += m[i][j] * (mult[displs[i] * super::_local + j + k * super::_map[i].second.size()] =   super::_deflation[k][super::_map[i].second[j]]);
                MPI_Isend(mult + displs[i] * super::_local, super::_map[i].second.size() * super::_local, Wrapper<K>::mpi_type(), super::_map[i].first, 11, super::_p.getCommunicator(), rqMult + super::_map.size() + i);
            }

            for(unsigned short i = 0; i < super::_map.size(); ++i) {
                int index;
                MPI_Waitany(super::_map.size(), rqMult, &index, MPI_STATUS_IGNORE);
                if(index < super::_signed)
                    for(unsigned short k = 0; k < (U ? super::_local : infoNeighbor[index]); ++k)
                        for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                            tmp[super::_map[index].second[j] + (offset[index + 1] + k) * super::_n] = - m[index][j] * mult[offset[index + 1] * nbMult + displs[index] * (U ? super::_local : infoNeighbor[index]) + j + k * super::_map[index].second.size()];
                else
                    for(unsigned short k = 0; k < (U ? super::_local : infoNeighbor[index]); ++k)
                        for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                            tmp[super::_map[index].second[j] + (offset[index + 1] + k) * super::_n] =   m[index][j] * mult[offset[index + 1] * nbMult + displs[index] * (U ? super::_local : infoNeighbor[index]) + j + k * super::_map[index].second.size()];
            }

            delete [] displs;

            if(offset[super::_map.size() + 1])
                super::_p.applyLocalPreconditioner(tmp, offset[super::_map.size() + 1]);

            MPI_Waitall(super::_map.size(), rqMult + super::_map.size(), MPI_STATUSES_IGNORE);
            delete [] rqMult;
            delete [] mult;

            unsigned int accumulate = 0;
            unsigned short stop = std::distance(super::_sparsity.cbegin(), std::upper_bound(super::_sparsity.cbegin(), super::_sparsity.cend(), super::_rank));
            if(S != 'S') {
                super::_offsets.reserve(super::_sparsity.size() + 1);
                for(unsigned short i = 0; i < stop; ++i) {
                    super::_offsets.emplace(super::_sparsity[i], accumulate);
                    accumulate += super::_n * (U ? super::_local : info[i]);
                }
            }
            else
                super::_offsets.reserve(super::_sparsity.size() + 1 - stop);
            super::_offsets.emplace(super::_rank, accumulate);
            accumulate += super::_n * super::_local;
            for(unsigned short i = stop; i < super::_sparsity.size(); ++i) {
                super::_offsets.emplace(super::_sparsity[i], accumulate);
                accumulate += super::_n * (U ? super::_local : info[i]);
            }

            work = new K[accumulate]();

            for(unsigned short i = 0; i < super::_signed; ++i) {
                accumulate = super::_local;
                for(unsigned short k = 0; k < super::_local; ++k)
                    for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                        work[super::_offsets[super::_rank] + super::_map[i].second[j] + k * super::_n] -= (in[i][k * super::_map[i].second.size() + j] = - m[i][j] * tmp[super::_map[i].second[j] + k * super::_n]);
                for(unsigned short l = (S != 'S' ? 0 : i); l < super::_map.size(); ++l) {
                    if(Q == FetiPrcndtnr::SUPERLUMPED && l != i && !std::binary_search(super::_vecSparsity[i].cbegin(), super::_vecSparsity[i].cend(), super::_map[l].first)) {
                        if(S != 'S' || !(l < super::_signed))
                            for(unsigned short k = 0; k < (U ? super::_local : infoNeighbor[l]); ++k)
                                for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                                    work[super::_offsets[super::_map[l].first] + super::_map[i].second[j] + k * super::_n] -= - m[i][j] * tmp[super::_map[i].second[j] + (offset[l + 1] + k) * super::_n];
                        continue;
                    }
                    for(unsigned short k = 0; k < (U ? super::_local : infoNeighbor[l]); ++k)
                        for(unsigned int j = 0; j < super::_map[i].second.size(); ++j) {
                            if(S != 'S' || !(l < super::_signed))
                                work[super::_offsets[super::_map[l].first] + super::_map[i].second[j] + k * super::_n] -= (in[i][(accumulate + k) * super::_map[i].second.size() + j] = - m[i][j] * tmp[super::_map[i].second[j] + (offset[l + 1] + k) * super::_n]);
                            else
                                in[i][(accumulate + k) * super::_map[i].second.size() + j] = - m[i][j] * tmp[super::_map[i].second[j] + (offset[l + 1] + k) * super::_n];
                        }
                    accumulate += U ? super::_local : infoNeighbor[l];
                }
                if(U || infoNeighbor[i])
                    MPI_Isend(in[i], super::_map[i].second.size() * accumulate, Wrapper<K>::mpi_type(), super::_map[i].first, 2, super::_p.getCommunicator(), rq++);
            }
            for(unsigned short i = super::_signed; i < super::_map.size(); ++i) {
                if(S != 'S') {
                    accumulate = super::_local;
                    for(unsigned short k = 0; k < super::_local; ++k)
                        for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                            work[super::_offsets[super::_rank] + super::_map[i].second[j] + k * super::_n] += (in[i][k * super::_map[i].second.size() + j] =   m[i][j] * tmp[super::_map[i].second[j] + k * super::_n]);
                }
                else {
                    accumulate = 0;
                    for(unsigned short k = 0; k < super::_local; ++k)
                        for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                            work[super::_offsets[super::_rank] + super::_map[i].second[j] + k * super::_n] += m[i][j] * tmp[super::_map[i].second[j] + k * super::_n];
                }
                for(unsigned short l = S != 'S' ? 0 : super::_signed; l < super::_map.size(); ++l) {
                    if(Q == FetiPrcndtnr::SUPERLUMPED && l != i && !std::binary_search(super::_vecSparsity[i].cbegin(), super::_vecSparsity[i].cend(), super::_map[l].first)) {
                        if(S != 'S' || !(l < i))
                            for(unsigned short k = 0; k < (U ? super::_local : infoNeighbor[l]); ++k)
                                for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                                    work[super::_offsets[super::_map[l].first] + super::_map[i].second[j] + k * super::_n] +=   m[i][j] * tmp[super::_map[i].second[j] + (offset[l + 1] + k) * super::_n];
                        continue;
                    }
                    for(unsigned short k = 0; k < (U ? super::_local : infoNeighbor[l]); ++k)
                        for(unsigned int j = 0; j < super::_map[i].second.size(); ++j) {
                            if(S != 'S' || !(l < i))
                                work[super::_offsets[super::_map[l].first] + super::_map[i].second[j] + k * super::_n] += (in[i][(accumulate + k) * super::_map[i].second.size() + j] =   m[i][j] * tmp[super::_map[i].second[j] + (offset[l + 1] + k) * super::_n]);
                            else
                                work[super::_offsets[super::_map[l].first] + super::_map[i].second[j] + k * super::_n] += m[i][j] * tmp[super::_map[i].second[j] + (offset[l + 1] + k) * super::_n];
                        }
                    if(S != 'S' || !(l < i))
                        accumulate += U ? super::_local : infoNeighbor[l];
                }
                if(U || infoNeighbor[i])
                    MPI_Isend(in[i], super::_map[i].second.size() * accumulate, Wrapper<K>::mpi_type(), super::_map[i].first, 2, super::_p.getCommunicator(), rq++);
            }
            delete [] tmp;
            delete [] offset;
            if(!U && Q != FetiPrcndtnr::SUPERLUMPED)
                delete [] infoNeighbor;
        }
        template<char S, bool U>
        void assembleForMaster(K* C, const K* in, const int& coefficients, unsigned short index, K* arrayC, unsigned short* const& infoNeighbor = nullptr) {
            applyFromNeighbor<S, U>(in, index, arrayC, infoNeighbor);
            if(++super::_consolidate == super::_map.size()) {
                if(S != 'S')
                    Blas<K>::gemm(&(Wrapper<K>::transc), "N", &coefficients, &(super::_local), &(super::_n), &(Wrapper<K>::d__1), arrayC, &(super::_n), *super::_deflation, super::_p.getLDR(), &(Wrapper<K>::d__0), C, &coefficients);
                else
                    for(unsigned short j = 0; j < super::_local; ++j) {
                        int local = coefficients + super::_local - j;
                        Blas<K>::gemv(&(Wrapper<K>::transc), &(super::_n), &local, &(Wrapper<K>::d__1), arrayC + super::_n * j, &(super::_n), super::_deflation[j], &i__1, &(Wrapper<K>::d__0), C - (j * (j - 1)) / 2 + j * (coefficients + super::_local), &i__1);
                    }
            }
        }
        template<char S, char N, bool U>
        void applyFromNeighborMaster(const K* in, unsigned short index, int* I, int* J, K* C, int coefficients, unsigned int offsetI, unsigned int* offsetJ, K* arrayC, unsigned short* const& infoNeighbor = nullptr) {
            assembleForMaster<S, U>(C, in, coefficients, index, arrayC, infoNeighbor);
            super::template assembleOperator<S, N, U>(I, J, coefficients, offsetI, offsetJ, infoNeighbor);
        }
};
#endif // HPDDM_FETI

#if HPDDM_BDD
template<class Preconditioner, class K>
class BddProjection : public OperatorBase<'c', Preconditioner, K> {
    private:
        typedef OperatorBase<'c', Preconditioner, K> super;
        template<char S, bool U>
        void applyFromNeighbor(const K* in, unsigned short index, K*& work, unsigned short* info) {
            std::vector<unsigned short>::const_iterator middle = std::lower_bound(super::_vecSparsity[index].cbegin(), super::_vecSparsity[index].cend(), super::_rank);
            unsigned int accumulate = 0;
            if(S != 'S' || !(index < super::_signed)) {
                for(unsigned short k = 0; k < (U ? super::_local : info[std::distance(super::_sparsity.cbegin(), std::lower_bound(super::_sparsity.cbegin(), super::_sparsity.cend(), super::_map[index].first))]); ++k)
                    for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                        work[super::_offsets[super::_map[index].first] + super::_map[index].second[j] + k * super::_n] += in[k * super::_map[index].second.size() + j];
                accumulate += (U ? super::_local : info[std::distance(super::_sparsity.cbegin(), std::lower_bound(super::_sparsity.cbegin(), super::_sparsity.cend(), super::_map[index].first))]) * super::_map[index].second.size();
            }
            std::vector<unsigned short>::const_iterator begin = super::_sparsity.cbegin();
            if(S != 'S')
                for(std::vector<unsigned short>::const_iterator it = super::_vecSparsity[index].cbegin(); it != middle; ++it) {
                    if(!U) {
                        std::vector<unsigned short>::const_iterator idx = std::lower_bound(begin, super::_sparsity.cend(), *it);
                        for(unsigned short k = 0; k < info[std::distance(super::_sparsity.cbegin(), idx)]; ++k)
                            for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
                        accumulate += info[std::distance(super::_sparsity.cbegin(), idx)] * super::_map[index].second.size();
                        begin = idx + 1;
                    }
                    else {
                        for(unsigned short k = 0; k < super::_local; ++k)
                            for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                                work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
                        accumulate += super::_local * super::_map[index].second.size();
                    }
                }
            for(unsigned short k = 0; k < super::_local; ++k) {
                for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                    work[super::_offsets[super::_rank] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
            }
            accumulate += super::_local * super::_map[index].second.size();
            for(std::vector<unsigned short>::const_iterator it = middle + 1; it < super::_vecSparsity[index].cend(); ++it) {
                if(!U) {
                    std::vector<unsigned short>::const_iterator idx = std::lower_bound(begin, super::_sparsity.cend(), *it);
                    for(unsigned short k = 0; k < info[std::distance(super::_sparsity.cbegin(), idx)]; ++k)
                        for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                            work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
                    accumulate += info[std::distance(super::_sparsity.cbegin(), idx)] * super::_map[index].second.size();
                    begin = idx + 1;
                }
                else {
                    for(unsigned short k = 0; k < super::_local; ++k)
                        for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                            work[super::_offsets[*it] + super::_map[index].second[j] + k * super::_n] += in[accumulate + k * super::_map[index].second.size() + j];
                    accumulate += super::_local * super::_map[index].second.size();
                }
            }
        }
    public:
        template<template<class> class Solver, char S, class T> friend class CoarseOperator;
        BddProjection(const Preconditioner& p, const unsigned short c) : super(p, std::move(c)) { }
        template<char S, bool U, class T>
        void applyToNeighbor(T& in, K*& work, MPI_Request*& rq, const unsigned short* info, T const& out = nullptr, MPI_Request* const& rqRecv = nullptr) {
            unsigned short* infoNeighbor;
            super::template initialize<S, U>(in, info, out, rqRecv, infoNeighbor);
            MPI_Request* rqMult = new MPI_Request[2 * super::_map.size()];
            unsigned int* offset = new unsigned int[super::_map.size() + 2];
            offset[0] = 0;
            offset[1] = super::_local;
            for(unsigned short i = 2; i < super::_map.size() + 2; ++i)
                offset[i] = offset[i - 1] + (U ? super::_local : infoNeighbor[i - 2]);
            const int nbMult = super::_p.getMult();
            K* mult = new K[offset[super::_map.size() + 1] * nbMult];
            unsigned short* displs = new unsigned short[super::_map.size() + 1];
            displs[0] = 0;
            for(unsigned short i = 0; i < super::_map.size(); ++i) {
                MPI_Irecv(mult + offset[i + 1] * nbMult + displs[i] * (U ? super::_local : infoNeighbor[i]), super::_map[i].second.size() * (U ? super::_local : infoNeighbor[i]), Wrapper<K>::mpi_type(), super::_map[i].first, 11, super::_p.getCommunicator(), rqMult + i);
                displs[i + 1] = displs[i] + super::_map[i].second.size();
            }

            K* tmp = new K[offset[super::_map.size() + 1] * super::_n]();
            const underlying_type<K>* const m = super::_p.getScaling();
            for(unsigned short i = 0; i < super::_map.size(); ++i) {
                for(unsigned short k = 0; k < super::_local; ++k)
                    for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                        tmp[super::_map[i].second[j] + k * super::_n] = (mult[displs[i] * super::_local + j + k * super::_map[i].second.size()] = m[super::_map[i].second[j]] * super::_deflation[k][super::_map[i].second[j]]);
                MPI_Isend(mult + displs[i] * super::_local, super::_map[i].second.size() * super::_local, Wrapper<K>::mpi_type(), super::_map[i].first, 11, super::_p.getCommunicator(), rqMult + super::_map.size() + i);
            }

            for(unsigned short i = 0; i < super::_map.size(); ++i) {
                int index;
                MPI_Waitany(super::_map.size(), rqMult, &index, MPI_STATUS_IGNORE);
                for(unsigned short k = 0; k < (U ? super::_local : infoNeighbor[index]); ++k)
                    for(unsigned int j = 0; j < super::_map[index].second.size(); ++j)
                        tmp[super::_map[index].second[j] + (offset[index + 1] + k) * super::_n] = mult[offset[index + 1] * nbMult + displs[index] * (U ? super::_local : infoNeighbor[index]) + j + k * super::_map[index].second.size()];
            }

            delete [] displs;

            if(offset[super::_map.size() + 1])
                super::_p.applyLocalSchurComplement(tmp, offset[super::_map.size() + 1]);

            MPI_Waitall(super::_map.size(), rqMult + super::_map.size(), MPI_STATUSES_IGNORE);
            delete [] rqMult;
            delete [] mult;

            unsigned int accumulate = 0;
            unsigned short stop = std::distance(super::_sparsity.cbegin(), std::upper_bound(super::_sparsity.cbegin(), super::_sparsity.cend(), super::_rank));
            if(S != 'S') {
                super::_offsets.reserve(super::_sparsity.size() + 1);
                for(unsigned short i = 0; i < stop; ++i) {
                    super::_offsets.emplace(super::_sparsity[i], accumulate);
                    accumulate += super::_n * (U ? super::_local : info[i]);
                }
            }
            else
                super::_offsets.reserve(super::_sparsity.size() + 1 - stop);
            super::_offsets.emplace(super::_rank, accumulate);
            accumulate += super::_n * super::_local;
            for(unsigned short i = stop; i < super::_sparsity.size(); ++i) {
                super::_offsets.emplace(super::_sparsity[i], accumulate);
                accumulate += super::_n * (U ? super::_local : info[i]);
            }

            work = new K[accumulate]();

            for(unsigned short i = 0; i < super::_map.size(); ++i) {
                if(i < super::_signed || S != 'S') {
                    accumulate = super::_local;
                    for(unsigned short k = 0; k < super::_local; ++k)
                        for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                            work[super::_offsets[super::_rank] + super::_map[i].second[j] + k * super::_n] = in[i][k * super::_map[i].second.size() + j] = tmp[super::_map[i].second[j] + k * super::_n];
                }
                else {
                    accumulate = 0;
                    for(unsigned short k = 0; k < super::_local; ++k)
                        for(unsigned int j = 0; j < super::_map[i].second.size(); ++j)
                            work[super::_offsets[super::_rank] + super::_map[i].second[j] + k * super::_n] = tmp[super::_map[i].second[j] + k * super::_n];
                }
                for(unsigned short l = S != 'S' ? 0 : std::min(i, super::_signed); l < super::_map.size(); ++l) {
                    for(unsigned short k = 0; k < (U ? super::_local : infoNeighbor[l]); ++k)
                        for(unsigned int j = 0; j < super::_map[i].second.size(); ++j) {
                            if(S != 'S' || !(l < std::max(i, super::_signed)))
                                work[super::_offsets[super::_map[l].first] + super::_map[i].second[j] + k * super::_n] = in[i][(accumulate + k) * super::_map[i].second.size() + j] = tmp[super::_map[i].second[j] + (offset[l + 1] + k) * super::_n];
                            else {
                                if(i < super::_signed)
                                    in[i][(accumulate + k) * super::_map[i].second.size() + j] = tmp[super::_map[i].second[j] + (offset[l + 1] + k) * super::_n];
                                else
                                    work[super::_offsets[super::_map[l].first] + super::_map[i].second[j] + k * super::_n] = tmp[super::_map[i].second[j] + (offset[l + 1] + k) * super::_n];
                            }
                        }
                    if(S != 'S' || !(l < i) || i < super::_signed)
                        accumulate += U ? super::_local : infoNeighbor[l];
                }
                if(U || infoNeighbor[i])
                    MPI_Isend(in[i], super::_map[i].second.size() * accumulate, Wrapper<K>::mpi_type(), super::_map[i].first, 2, super::_p.getCommunicator(), rq++);
            }
            delete [] tmp;
            delete [] offset;
            if(!U)
                delete [] infoNeighbor;
        }
        template<char S, bool U>
        void assembleForMaster(K* C, const K* in, const int& coefficients, unsigned short index, K* arrayC, unsigned short* const& infoNeighbor = nullptr) {
            applyFromNeighbor<S, U>(in, index, arrayC, infoNeighbor);
            if(++super::_consolidate == super::_map.size()) {
                const underlying_type<K>* const m = super::_p.getScaling();
                for(unsigned short j = 0; j < coefficients + (S == 'S') * super::_local; ++j)
                    Wrapper<K>::diag(super::_n, m, arrayC + j * super::_n);
                if(S != 'S')
                    Blas<K>::gemm(&(Wrapper<K>::transc), "T", &coefficients, &(super::_local), &(super::_n), &(Wrapper<K>::d__1), arrayC, &(super::_n), *super::_deflation, super::_p.getLDR(), &(Wrapper<K>::d__0), C, &coefficients);
                else
                    for(unsigned short j = 0; j < super::_local; ++j) {
                        int local = coefficients + super::_local - j;
                        Blas<K>::gemv(&(Wrapper<K>::transc), &(super::_n), &local, &(Wrapper<K>::d__1), arrayC + super::_n * j, &(super::_n), super::_deflation[j], &i__1, &(Wrapper<K>::d__0), C - (j * (j - 1)) / 2 + j * (coefficients + super::_local), &i__1);
                    }
            }
        }
        template<char S, char N, bool U>
        void applyFromNeighborMaster(const K* in, unsigned short index, int* I, int* J, K* C, int coefficients, unsigned int offsetI, unsigned int* offsetJ, K* arrayC, unsigned short* const& infoNeighbor = nullptr) {
            assembleForMaster<S, U>(C, in, coefficients, index, arrayC, infoNeighbor);
            super::template assembleOperator<S, N, U>(I, J, coefficients, offsetI, offsetJ, infoNeighbor);
        }
};
#endif // HPDDM_BDD
} // HPDDM
#endif // _HPDDM_OPERATOR_
