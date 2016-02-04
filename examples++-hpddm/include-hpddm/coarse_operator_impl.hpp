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

#ifndef _HPDDM_COARSE_OPERATOR_IMPL_
#define _HPDDM_COARSE_OPERATOR_IMPL_

#include "coarse_operator.hpp"

namespace HPDDM {
template<template<class> class Solver, char S, class K>
template<bool exclude>
inline void CoarseOperator<Solver, S, K>::constructionCommunicator(const MPI_Comm& comm) {
    MPI_Comm_size(comm, &_sizeWorld);
    MPI_Comm_rank(comm, &_rankWorld);
    Option& opt = *Option::get();
    unsigned short p = opt["master_p"];
#ifndef DSUITESPARSE
    if(p > _sizeWorld / 2 && _sizeWorld > 1) {
        p = opt["master_p"] = _sizeWorld / 2;
        if(_rankWorld == 0)
            std::cout << "WARNING -- the number of master processes was set to a value greater than MPI_Comm_size / 2, the value has been reset to " << p << std::endl;
    }
    else if(p < 1)
#endif
        p = opt["master_p"] = 1;
    if(p == 1) {
        MPI_Comm_dup(comm, &_scatterComm);
        _gatherComm = _scatterComm;
        if(_rankWorld != 0)
            Solver<K>::_communicator = MPI_COMM_NULL;
        else
            Solver<K>::_communicator = MPI_COMM_SELF;
        Solver<K>::_rank = 0;
        Solver<K>::_ldistribution = new int[1]();
    }
    else {
        MPI_Group master, split;
        MPI_Group world;
        MPI_Comm_group(comm, &world);
        int* ps;
        unsigned int tmp;
        Solver<K>::_ldistribution = new int[p];
        char T = opt["master_topology"];
        if(T == 2) {
            // Here, it is assumed that all subdomains have the same number of coarse degrees of freedom as the rank 0 ! (only true when the distribution is uniform)
            float area = _sizeWorld *_sizeWorld / (2.0 * p);
            *Solver<K>::_ldistribution = 0;
            for(unsigned short i = 1; i < p; ++i)
                Solver<K>::_ldistribution[i] = static_cast<int>(_sizeWorld - std::sqrt(std::max(_sizeWorld * _sizeWorld - 2 * _sizeWorld * Solver<K>::_ldistribution[i - 1] - 2 * area + Solver<K>::_ldistribution[i - 1] * Solver<K>::_ldistribution[i - 1], 1.0f)) + 0.5);
            int* idx = std::upper_bound(Solver<K>::_ldistribution, Solver<K>::_ldistribution + p, _rankWorld);
            unsigned short i = idx - Solver<K>::_ldistribution;
            tmp = (i == p) ? _sizeWorld - Solver<K>::_ldistribution[i - 1] : Solver<K>::_ldistribution[i] - Solver<K>::_ldistribution[i - 1];
            ps = new int[tmp];
            for(unsigned int j = 0; j < tmp; ++j)
                ps[j] = Solver<K>::_ldistribution[i - 1] + j;
        }
#ifndef HPDDM_CONTIGUOUS
        else if(T == 1) {
            if(_rankWorld == p - 1 || _rankWorld > p - 1 + (p - 1) * ((_sizeWorld - p) / p))
                tmp = _sizeWorld - (p - 1) * (_sizeWorld / p);
            else
                tmp = _sizeWorld / p;
            ps = new int[tmp];
            if(_rankWorld < p)
                ps[0] = _rankWorld;
            else {
                if(tmp == _sizeWorld / p)
                    ps[0] = (_rankWorld - p) / ((_sizeWorld - p) / p);
                else
                    ps[0] = p - 1;
            }
            unsigned int offset = ps[0] * (_sizeWorld / p - 1) + p - 1;
            std::iota(ps + 1, ps + tmp, offset + 1);
            std::iota(Solver<K>::_ldistribution, Solver<K>::_ldistribution + p, 0);
        }
#endif
        else {
            if(T != 0)
                opt["master_topology"] = 0;
            if(_rankWorld < (p - 1) * (_sizeWorld / p))
                tmp = _sizeWorld / p;
            else
                tmp = _sizeWorld - (p - 1) * (_sizeWorld / p);
            ps = new int[tmp];
            unsigned int offset;
            if(tmp != _sizeWorld / p)
                offset = _sizeWorld - tmp;
            else
                offset = (_sizeWorld / p) * (_rankWorld / (_sizeWorld / p));
            std::iota(ps, ps + tmp, offset);
            for(unsigned short i = 0; i < p; ++i)
                Solver<K>::_ldistribution[i] = i * (_sizeWorld / p);
        }
        MPI_Group_incl(world, p, Solver<K>::_ldistribution, &master);
        MPI_Group_incl(world, tmp, ps, &split);
        delete [] ps;

        MPI_Comm_create(comm, master, &(Solver<K>::_communicator));
        if(Solver<K>::_communicator != MPI_COMM_NULL)
            MPI_Comm_rank(Solver<K>::_communicator, &(Solver<K>::_rank));
        MPI_Comm_create(comm, split, &_scatterComm);

        MPI_Group_free(&master);
        MPI_Group_free(&split);

        if(!exclude)
            MPI_Comm_dup(comm, &_gatherComm);
        else {
            MPI_Group global;
            MPI_Group_excl(world, p - 1, Solver<K>::_ldistribution + 1, &global);
            MPI_Comm_create(comm, global, &_gatherComm);
            MPI_Group_free(&global);
        }
        MPI_Group_free(&world);
    }
}

template<template<class> class Solver, char S, class K>
template<bool U, typename Solver<K>::Distribution D, bool excluded>
inline void CoarseOperator<Solver, S, K>::constructionCollective(const unsigned short* info, unsigned short p, const unsigned short* infoSplit) {
    if(!U) {
        if(excluded)
            _sizeWorld -= p;
        Solver<K>::_gatherCounts = new int[2 * _sizeWorld];
        Solver<K>::_displs = Solver<K>::_gatherCounts + _sizeWorld;

        Solver<K>::_gatherCounts[0] = info[0];
        Solver<K>::_displs[0] = 0;
        for(unsigned int i = 1, j = 1; j < _sizeWorld; ++i)
            if(!excluded || info[i] != 0)
                Solver<K>::_gatherCounts[j++] = info[i];
        std::partial_sum(Solver<K>::_gatherCounts, Solver<K>::_gatherCounts + _sizeWorld - 1, Solver<K>::_displs + 1);
        if(excluded)
            _sizeWorld += p;
        if(D == DMatrix::DISTRIBUTED_SOL) {
            Solver<K>::_gatherSplitCounts = new int[2 * _sizeSplit];
            Solver<K>::_displsSplit = Solver<K>::_gatherSplitCounts + _sizeSplit;
            std::copy_n(infoSplit, _sizeSplit, Solver<K>::_gatherSplitCounts);
            Solver<K>::_displsSplit[0] = 0;
            std::partial_sum(Solver<K>::_gatherSplitCounts, Solver<K>::_gatherSplitCounts + _sizeSplit - 1, Solver<K>::_displsSplit + 1);
        }
    }
    else {
        Solver<K>::_gatherCounts = new int[1];
        *Solver<K>::_gatherCounts = _local;
    }
}

template<template<class> class Solver, char S, class K>
template<char T, bool U, bool excluded>
inline void CoarseOperator<Solver, S, K>::constructionMap(unsigned short p, const unsigned short* info) {
    if(T == 0) {
        if(!U) {
            unsigned int accumulate = 0;
            for(unsigned short i = 0; i < p - 1; accumulate += Solver<K>::_ldistribution[i++])
                Solver<K>::_ldistribution[i] = std::accumulate(info + i * (_sizeWorld / p), info + (i + 1) * (_sizeWorld / p), 0);
            Solver<K>::_ldistribution[p - 1] = Solver<K>::_n - accumulate;
        }
        else {
            if(p == 1)
                *Solver<K>::_ldistribution = Solver<K>::_n;
            else {
                std::fill_n(Solver<K>::_ldistribution, p - 1, _local * (_sizeWorld / p - excluded));
                Solver<K>::_ldistribution[p - 1] = Solver<K>::_n - _local * (_sizeWorld / p - excluded) * (p - 1);
            }
        }
    }
    else if(T == 1) {
        Solver<K>::_idistribution = new int[Solver<K>::_n];
        unsigned int j = 0;
        if(!excluded)
            for(unsigned int i = 0; i < p * (_sizeWorld / p); ++i) {
                unsigned int offset;
                if(i % (_sizeWorld / p) == 0) {
                    j = i / (_sizeWorld / p);
                    offset = U ? (_sizeWorld / p) * _local * j : (std::accumulate(info, info + j, 0) + std::accumulate(info + p, info + p + j * (_sizeWorld / p - 1), 0));
                }
                else {
                    j = p - 1 + i - i / (_sizeWorld / p);
                    offset  = U ? _local * (1 + i  / (_sizeWorld / p)) : std::accumulate(info, info + 1 + i / (_sizeWorld / p), 0);
                    offset += U ? (j - p) * _local : std::accumulate(info + p, info + j, 0);
                }
                std::iota(Solver<K>::_idistribution + offset, Solver<K>::_idistribution + offset + (U ? _local : info[j]), U ? _local * j : std::accumulate(info, info + j, 0));
                if(i % (_sizeWorld / p) != 0)
                    j = offset + (U ? _local : info[j]);
            }
        std::iota(Solver<K>::_idistribution + j, Solver<K>::_idistribution + Solver<K>::_n, j);
        if(!U) {
            unsigned int accumulate = 0;
            for(unsigned short i = 0; i < p - 1; accumulate += Solver<K>::_ldistribution[i++])
                Solver<K>::_ldistribution[i] = std::accumulate(info + p + i * (_sizeWorld / p - 1), info + p + (i + 1) * (_sizeWorld / p - 1), info[i]);
            Solver<K>::_ldistribution[p - 1] = Solver<K>::_n - accumulate;
        }
        else {
            std::fill_n(Solver<K>::_ldistribution, p - 1, _local * (_sizeWorld / p - excluded));
            Solver<K>::_ldistribution[p - 1] = Solver<K>::_n - _local * (_sizeWorld / p - excluded) * (p - 1);
        }
    }
    else if(T == 2) {
        if(!U) {
            unsigned int accumulate = 0;
            for(unsigned short i = 0; i < p - 1; accumulate += Solver<K>::_ldistribution[i++])
                Solver<K>::_ldistribution[i] = std::accumulate(info + Solver<K>::_ldistribution[i], info + Solver<K>::_ldistribution[i + 1], 0);
            Solver<K>::_ldistribution[p - 1] = Solver<K>::_n - accumulate;
        }
        else {
            for(unsigned short i = 0; i < p - 1; ++i)
                Solver<K>::_ldistribution[i] = (Solver<K>::_ldistribution[i + 1] - Solver<K>::_ldistribution[i] - excluded) * _local;
            Solver<K>::_ldistribution[p - 1] = Solver<K>::_n - (Solver<K>::_ldistribution[p - 1] - (excluded ? p - 1 : 0)) * _local;
        }
    }
}

template<template<class> class Solver, char S, class K>
template<unsigned short U, unsigned short excluded, class Operator>
inline std::pair<MPI_Request, const K*>* CoarseOperator<Solver, S, K>::construction(Operator&& v, const MPI_Comm& comm) {
    static_assert(Solver<K>::_numbering == 'C' || Solver<K>::_numbering == 'F', "Unknown numbering");
    static_assert(Operator::_pattern == 's' || Operator::_pattern == 'c', "Unknown pattern");
    constructionCommunicator<static_cast<bool>(excluded)>(comm);
    if(excluded > 0 && Solver<K>::_communicator != MPI_COMM_NULL) {
        int result;
        MPI_Comm_compare(v._p.getCommunicator(), Solver<K>::_communicator, &result);
        if(result != MPI_CONGRUENT)
            std::cerr << "The communicators for the coarse operator don't match those of the domain decomposition" << std::endl;
    }
    if(Operator::_pattern == 'c')
        v.adjustConnectivity(_scatterComm);
    Solver<K>::initialize();
    Option& opt = *Option::get();
    if(U == 2 && _local == 0)
        _offset = true;
    switch(static_cast<int>(opt["master_topology"])) {
#ifndef HPDDM_CONTIGUOUS
        case  1: return constructionMatrix<1, U, excluded>(v);
#endif
        case  2: return constructionMatrix<2, U, excluded>(v);
        default: return constructionMatrix<0, U, excluded>(v);
    }
}

template<template<class> class Solver, char S, class K>
template<char T, unsigned short U, unsigned short excluded, class Operator>
inline std::pair<MPI_Request, const K*>* CoarseOperator<Solver, S, K>::constructionMatrix(Operator& v) {
    unsigned short* const info = new unsigned short[(U != 1 ? 3 : 1) + v.getConnectivity()];
    const std::vector<unsigned short>& sparsity = v.getPattern();
    info[0] = sparsity.size(); // number of intersections
    int rank;
    MPI_Comm_rank(v._p.getCommunicator(), &rank);
    const unsigned short first = (S == 'S' ? std::distance(sparsity.cbegin(), std::upper_bound(sparsity.cbegin(), sparsity.cend(), rank)) : 0);
    int rankSplit;
    MPI_Comm_size(_scatterComm, &_sizeSplit);
    MPI_Comm_rank(_scatterComm, &rankSplit);
    unsigned short* infoNeighbor;

    K*     sendMaster;
    unsigned int size;
    int* I;
    int* J;
    K*   C;

    const Option& opt = *Option::get();
    unsigned short p = static_cast<int>(opt["master_p"]);
    if(U != 1) {
        infoNeighbor = new unsigned short[info[0]];
        info[1] = (excluded == 2 ? 0 : _local); // number of eigenvalues
        MPI_Request* rq = v._p.getRq();
        if(excluded == 0) {
            if(T != 2)
                for(unsigned short i = 0; i < info[0]; ++i) {
                    if(!(T == 1 && sparsity[i] < p) &&
                       !(T == 0 && (sparsity[i] % (_sizeWorld / p) == 0) && sparsity[i] < p * (_sizeWorld / p)))
                        MPI_Isend(info + 1, 1, MPI_UNSIGNED_SHORT, sparsity[i], 1, v._p.getCommunicator(), rq++);
                }
            else
                for(unsigned short i = 0; i < info[0]; ++i) {
                    if(!std::binary_search(Solver<K>::_ldistribution, Solver<K>::_ldistribution + p, sparsity[i]))
                        MPI_Isend(info + 1, 1, MPI_UNSIGNED_SHORT, sparsity[i], 1, v._p.getCommunicator(), rq++);
                }
        }
        else if(excluded < 2) {
            for(unsigned short i = 0; i < info[0]; ++i)
                MPI_Isend(info + 1, 1, MPI_UNSIGNED_SHORT, sparsity[i], 1, v._p.getCommunicator(), rq++);
        }
        if(rankSplit != 0) {
            for(unsigned short i = 0; i < info[0]; ++i)
                MPI_Irecv(infoNeighbor + i, 1, MPI_UNSIGNED_SHORT, sparsity[i], 1, v._p.getCommunicator(), rq + i);
            unsigned int tmp = (S != 'S' ? _local : 0);
            for(unsigned short i = 0; i < info[0]; ++i) {
                int index;
                MPI_Waitany(info[0], rq, &index, MPI_STATUS_IGNORE);
                if(!(S == 'S' && sparsity[index] < rank))
                    tmp += infoNeighbor[index];
            }
            if(S != 'S')
                size = _local * tmp;
            else {
                info[0] -= first;
                size = _local * tmp + _local * (_local + 1) / 2;
            }
            info[2] = size;
            if(_local) {
                sendMaster = new K[size];
                if(excluded == 0)
                    std::copy_n(sparsity.cbegin() + first, info[0], info + (U != 1 ? 3 : 1));
                else {
                    if(T == 0 || T == 2) {
                        for(unsigned short i = 0; i < info[0]; ++i) {
                            info[(U != 1 ? 3 : 1) + i] = sparsity[i + first] + 1;
                            for(unsigned short j = 0; j < p - 1 && info[(U != 1 ? 3 : 1) + i] >= (T == 0 ? (_sizeWorld / p) * (j + 1) : Solver<K>::_ldistribution[j + 1]); ++j)
                                ++info[(U != 1 ? 3 : 1) + i];
                        }
                    }
                    else if(T == 1) {
                        for(unsigned short i = 0; i < info[0]; ++i)
                            info[(U != 1 ? 3 : 1) + i] = p + sparsity[i + first];
                    }
                }
            }
        }
        int nbRq = std::distance(v._p.getRq(), rq);
        MPI_Waitall(nbRq, rq - nbRq, MPI_STATUSES_IGNORE);
    }
    else {
        infoNeighbor = nullptr;
        if(rankSplit != 0) {
            if(S == 'S') {
                info[0] -= first;
                size = _local * _local * info[0] + _local * (_local + 1) / 2;
            }
            else
                size = _local * _local * (1 + info[0]);
            sendMaster = new K[size];
            std::copy_n(sparsity.cbegin() + first, info[0], info + (U != 1 ? 3 : 1));
        }
    }
    unsigned short** infoSplit;
    unsigned int*    offsetIdx;
    unsigned short*  infoWorld;

    unsigned int offset;
#ifdef HPDDM_CSR_CO
    unsigned int nrow;
#ifdef HPDDM_LOC2GLOB
    int* loc2glob;
#endif
#endif

    if(rankSplit != 0)
        MPI_Gather(info, (U != 1 ? 3 : 1) + v.getConnectivity(), MPI_UNSIGNED_SHORT, NULL, 0, MPI_DATATYPE_NULL, 0, _scatterComm);
    else {
        size = 0;
        infoSplit = new unsigned short*[_sizeSplit];
        *infoSplit = new unsigned short[_sizeSplit * ((U != 1 ? 3 : 1) + v.getConnectivity()) + (U != 1) * _sizeWorld];
        MPI_Gather(info, (U != 1 ? 3 : 1) + v.getConnectivity(), MPI_UNSIGNED_SHORT, *infoSplit, (U != 1 ? 3 : 1) + v.getConnectivity(), MPI_UNSIGNED_SHORT, 0, _scatterComm);
        for(unsigned int i = 1; i < _sizeSplit; ++i)
            infoSplit[i] = *infoSplit + i * ((U != 1 ? 3 : 1) + v.getConnectivity());
        if(S == 'S' && Operator::_pattern == 's')
            **infoSplit -= first;
        offsetIdx = new unsigned int[std::max(_sizeSplit - 1, 2 * p)];
        if(U != 1) {
            infoWorld = *infoSplit + _sizeSplit * (3 + v.getConnectivity());
            int* recvcounts = reinterpret_cast<int*>(offsetIdx);
            int* displs = recvcounts + p;
            displs[0] = 0;
            if(T == 2) {
                std::adjacent_difference(Solver<K>::_ldistribution + 1, Solver<K>::_ldistribution + p, recvcounts);
                recvcounts[p - 1] = _sizeWorld - Solver<K>::_ldistribution[p - 1];
            }
            else {
                std::fill_n(recvcounts, p - 1, _sizeWorld / p);
                recvcounts[p - 1] = _sizeWorld - (p - 1) * (_sizeWorld / p);
            }
            std::partial_sum(recvcounts, recvcounts + p - 1, displs + 1);
            for(unsigned int i = 0; i < _sizeSplit; ++i)
                infoWorld[displs[Solver<K>::_rank] + i] = infoSplit[i][1];
#ifdef HPDDM_CSR_CO
            nrow = std::accumulate(infoWorld + displs[Solver<K>::_rank], infoWorld + displs[Solver<K>::_rank] + _sizeSplit, 0);
#endif
            MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, infoWorld, recvcounts, displs, MPI_UNSIGNED_SHORT, Solver<K>::_communicator);
            if(T == 1) {
                unsigned int i = (p - 1) * (_sizeWorld / p);
                for(unsigned short k = p - 1, j = 1; k-- > 0; i -= _sizeWorld / p, ++j) {
                    recvcounts[k] = infoWorld[i];
                    std::copy_backward(infoWorld + k * (_sizeWorld / p), infoWorld + (k + 1) * (_sizeWorld / p), infoWorld + (k + 1) * (_sizeWorld / p) + j);
                }
                std::copy_n(recvcounts, p - 1, infoWorld + 1);
            }
            offset = std::accumulate(infoWorld, infoWorld + _rankWorld, 0);
            Solver<K>::_n = std::accumulate(infoWorld + _rankWorld, infoWorld + _sizeWorld, offset);
            if(Solver<K>::_numbering == 'F')
                ++offset;
            unsigned short tmp = 0;
            for(unsigned short i = 0; i < info[0]; ++i) {
                infoNeighbor[i] = infoWorld[sparsity[i]];
                if(!(S == 'S' && i < first))
                    tmp += infoNeighbor[i];
            }
            for(unsigned short k = 1; k < _sizeSplit; size += infoSplit[k++][2])
                offsetIdx[k - 1] = size;
            if(excluded < 2)
                size += _local * tmp + (S == 'S' ? _local * (_local + 1) / 2 : _local * _local);
            if(S == 'S')
                info[0] -= first;
        }
        else {
            Solver<K>::_n = (_sizeWorld - (excluded == 2 ? p : 0)) * _local;
            offset = (_rankWorld - (excluded == 2 ? rank : 0)) * _local + (Solver<K>::_numbering == 'F');
#ifdef HPDDM_CSR_CO
            nrow = (_sizeSplit - (excluded == 2)) * _local;
#endif
            if(S == 'S') {
                for(unsigned short i = 1; i < _sizeSplit; size += infoSplit[i++][0])
                    offsetIdx[i - 1] = size * _local * _local + (i - 1) * _local * (_local + 1) / 2;
                info[0] -= first;
                size = (size + info[0]) * _local * _local + _local * (_local + 1) / 2 * (_sizeSplit - (excluded == 2));
            }
            else {
                for(unsigned short i = 1; i < _sizeSplit; size += infoSplit[i++][0])
                    offsetIdx[i - 1] = (i - 1 + size) * _local * _local;
                size = (size + info[0] + _sizeSplit - (excluded == 2)) * _local * _local;
            }
        }
#ifdef HPDDM_CSR_CO
        I = new int[nrow + 1 + size];
        I[0] = (Solver<K>::_numbering == 'F');
        J = I + nrow + 1;
#ifdef HPDDM_LOC2GLOB
#ifndef HPDDM_CONTIGUOUS
        loc2glob = new int[nrow];
#else
        loc2glob = new int[2];
#endif
#endif
#else
        I = new int[2 * size];
        J = I + size;
#endif
        C = new K[size];
    }
    const vectorNeighbor& M = v._p.getMap();

    MPI_Request* rqSend = v._p.getRq();
    MPI_Request* rqRecv;

    K** sendNeighbor = v._p.getBuffer();
    K** recvNeighbor;
    int coefficients = (U == 1 ? _local * (info[0] + (S != 'S')) : std::accumulate(infoNeighbor + first, infoNeighbor + sparsity.size(), (S == 'S' ? 0 : _local)));
    if(Operator::_pattern == 's') {
#if HPDDM_ICOLLECTIVE
        rqRecv = rqSend + (S != 'S' ? info[0] : first);
#else
        rqRecv = (rankSplit == 0 ? new MPI_Request[_sizeSplit - 1 + info[0]] : rqSend + (S != 'S' ? info[0] : first));
#endif
        unsigned int accumulate = 0;
        for(unsigned short i = 0; i < (S != 'S' ? info[0] : first); ++i)
            if(U == 1 || infoNeighbor[i])
                accumulate += _local * M[i].second.size();
        if(U == 1 || _local)
            for(unsigned short i = 0; i < info[0]; ++i)
                accumulate += (U == 1 ? _local : infoNeighbor[i + first]) * M[i + first].second.size();
        if(excluded < 2)
            *sendNeighbor = new K[accumulate];
        accumulate = 0;
        for(unsigned short i = 0; i < (S != 'S' ? info[0] : first); ++i) {
            sendNeighbor[i] = *sendNeighbor + accumulate;
            if(U == 1 || infoNeighbor[i])
                accumulate += _local * M[i].second.size();
        }
        recvNeighbor = (U == 1 || _local ? sendNeighbor + (S != 'S' ? info[0] : first) : nullptr);
        if(U == 1 || _local) {
            for(unsigned short i = 0; i < info[0]; ++i) {
                recvNeighbor[i] = *sendNeighbor + accumulate;
                MPI_Irecv(recvNeighbor[i], (U == 1 ? _local : infoNeighbor[i + first]) * M[i + first].second.size(), Wrapper<K>::mpi_type(), M[i + first].first, 2, v._p.getCommunicator(), rqRecv + i);
                accumulate += (U == 1 ? _local : infoNeighbor[i + first]) * M[i + first].second.size();
            }
        }
        else
            std::fill_n(rqRecv, info[0], MPI_REQUEST_NULL);
    }
    else {
#if HPDDM_ICOLLECTIVE
        rqRecv = rqSend + M.size();
#else
        rqRecv = (rankSplit == 0 ? new MPI_Request[_sizeSplit - 1 + M.size()] : (rqSend + M.size()));
#endif
        recvNeighbor = (U == 1 || _local) ? sendNeighbor + M.size() : nullptr;
    }
    K* work = nullptr;
    if(Operator::_pattern == 's' && excluded < 2) {
        const K* const* const& EV = v._p.getVectors();
        const int n = v._p.getDof();
        v.initialize(n * (U == 1 || info[0] == 0 ? _local : std::max(static_cast<unsigned short>(_local), *std::max_element(infoNeighbor + first, infoNeighbor + sparsity.size()))), work, S != 'S' ? info[0] : first);
        v.template applyToNeighbor<S, U == 1>(sendNeighbor, work, rqSend, infoNeighbor);
        if(S != 'S') {
            unsigned short before = 0;
            for(unsigned short j = 0; j < info[0] && sparsity[j] < rank; ++j)
                before += (U == 1 ? _local : infoNeighbor[j]);
            K* const pt = (rankSplit != 0 ? sendMaster + before : C + before);
            Blas<K>::gemm(&(Wrapper<K>::transc), "N", &_local, &_local, &n, &(Wrapper<K>::d__1), work, &n, *EV, &n, &(Wrapper<K>::d__0), pt, &coefficients);
            Wrapper<K>::template imatcopy<'R'>(_local, _local, pt, coefficients, coefficients);
            if(rankSplit == 0)
                for(unsigned short j = 0; j < _local; ++j) {
#ifndef HPDDM_CSR_CO
                    std::fill_n(I + before + j * coefficients, _local, offset + j);
#endif
                    std::iota(J + before + j * coefficients, J + before + j * coefficients + _local, offset);
                }
        }
        else {
            if(rankSplit != 0)
                if(coefficients >= _local) {
                    Blas<K>::gemm(&(Wrapper<K>::transc), "N", &_local, &_local, &n, &(Wrapper<K>::d__1), *EV, &n, work, &n, &(Wrapper<K>::d__0), sendMaster, &_local);
                    for(unsigned short j = _local; j-- > 0; )
                        std::copy_backward(sendMaster + j * (_local + 1), sendMaster + (j + 1) * _local, sendMaster - (j * (j + 1)) / 2 + j * coefficients + (j + 1) * _local);
                }
                else
                    for(unsigned short j = 0; j < _local; ++j) {
                        int local = _local - j;
                        Blas<K>::gemv(&(Wrapper<K>::transc), &n, &local, &(Wrapper<K>::d__1), EV[j], &n, work + n * j, &i__1, &(Wrapper<K>::d__0), sendMaster - (j * (j - 1)) / 2 + j * (coefficients + _local), &i__1);
                    }
            else {
                if(coefficients >= _local)
                    Blas<K>::gemm(&(Wrapper<K>::transc), "N", &_local, &_local, &n, &(Wrapper<K>::d__1), *EV, &n, work, &n, &(Wrapper<K>::d__0), C, &_local);
                for(unsigned short j = _local; j-- > 0; ) {
#ifndef HPDDM_CSR_CO
                    std::fill_n(I + j * (coefficients + _local) - (j * (j - 1)) / 2, _local - j, offset + j);
#endif
                    std::iota(J + j * (coefficients + _local - 1) - (j * (j - 1)) / 2 + j, J + j * (coefficients + _local - 1) - (j * (j - 1)) / 2 + _local, offset + j);
                    if(coefficients >= _local)
                        std::copy_backward(C + j * (_local + 1), C + (j + 1) * _local, C - (j * (j + 1)) / 2 + j * coefficients + (j + 1) * _local);
                    else {
                        int local = _local - j;
                        Blas<K>::gemv(&(Wrapper<K>::transc), &n, &local, &(Wrapper<K>::d__1), EV[j], &n, work + n * j, &i__1, &(Wrapper<K>::d__0), C - (j * (j - 1)) / 2 + j * (coefficients + _local), &i__1);
                    }
                }
            }
        }
    }
    else if(Operator::_pattern != 's' && excluded < 2)
        v.template applyToNeighbor<S, U == 1>(sendNeighbor, work, rqSend, U == 1 ? nullptr : infoNeighbor, recvNeighbor, rqRecv);
    std::pair<MPI_Request, const K*>* ret = nullptr;
    if(rankSplit != 0) {
        if(U == 1 || _local) {
            if(Operator::_pattern == 's') {
                unsigned int* offsetArray = new unsigned int[info[0]];
                if(S != 'S')
                    offsetArray[0] = M[0].first > rank ? _local : 0;
                else if(info[0] > 0)
                    offsetArray[0] = _local;
                for(unsigned short k = 1; k < info[0]; ++k) {
                    offsetArray[k] = offsetArray[k - 1] + (U == 1 ? _local : infoNeighbor[k - 1 + first]);
                    if(S != 'S' && sparsity[k - 1] < rank && sparsity[k] > rank)
                        offsetArray[k] += _local;
                }
                for(unsigned short k = 0; k < info[0]; ++k) {
                    int index;
                    MPI_Waitany(info[0], rqRecv, &index, MPI_STATUS_IGNORE);
                    v.template assembleForMaster<S, U == 1>(sendMaster + offsetArray[index], recvNeighbor[index], coefficients + (S == 'S' ? _local - 1 : 0), index + first, work, infoNeighbor + first + index);
                }
                delete [] offsetArray;
            }
            else {
                for(unsigned short k = 0; k < M.size(); ++k) {
                    int index;
                    MPI_Waitany(M.size(), rqRecv, &index, MPI_STATUS_IGNORE);
                    v.template assembleForMaster<S, U == 1>(sendMaster, recvNeighbor[index], coefficients, index, work, infoNeighbor);
                }
            }
            if(excluded > 0) {
                ret = new std::pair<MPI_Request, const K*>(MPI_REQUEST_NULL, sendMaster);
#if HPDDM_ICOLLECTIVE
                MPI_Igatherv(sendMaster, size, Wrapper<K>::mpi_type(), NULL, 0, 0, MPI_DATATYPE_NULL, 0, _scatterComm, &(ret->first));
#else
                MPI_Isend(sendMaster, size, Wrapper<K>::mpi_type(), 0, 3, _scatterComm, &(ret->first));
#endif
            }
            else {
#if HPDDM_ICOLLECTIVE
                MPI_Request rq;
                MPI_Igatherv(sendMaster, size, Wrapper<K>::mpi_type(), NULL, 0, 0, MPI_DATATYPE_NULL, 0, _scatterComm, &rq);
                MPI_Wait(&rq, MPI_STATUS_IGNORE);
#else
                MPI_Send(sendMaster, size, Wrapper<K>::mpi_type(), 0, 3, _scatterComm);
#endif
                delete [] sendMaster;
            }
        }
        delete [] info;
        _sizeRHS = _local;
        if(U != 1)
            delete [] infoNeighbor;
        if(U == 0)
            Solver<K>::_displs = &_rankWorld;
        int nbRq = std::distance(v._p.getRq(), rqSend);
        MPI_Waitall(nbRq, rqSend - nbRq, MPI_STATUSES_IGNORE);
        delete [] work;
    }
    else {
        unsigned short rankRelative = (T == 0 || T == 2) ? _rankWorld : p + _rankWorld * ((_sizeWorld / p) - 1) - 1;
        unsigned int* offsetPosition;
        unsigned int idx;
        if(excluded < 2) {
            std::for_each(offsetIdx, offsetIdx + _sizeSplit - 1, [&](unsigned int& i) { i += coefficients * _local + (S == 'S' ? (_local * (_local + 1)) / 2 : 0); });
            idx = Operator::_pattern == 's' ? info[0] : M.size();
        }
        else
            idx = 0;
#if HPDDM_ICOLLECTIVE
        int* counts = new int[2 * _sizeSplit];
        counts[0] = 0;
        counts[_sizeSplit] = (U == 1 ? (S == 'S' ? _local * infoSplit[0][0] * _local + _local * (_local + 1) / 2 : (_local * infoSplit[0][0] + _local) * _local) : infoSplit[0][2]);
        for(unsigned short k = 1; k < _sizeSplit; ++k) {
            counts[k] = offsetIdx[k - 1];
            counts[_sizeSplit + k] = (U == 1 ? (counts[_sizeSplit + k - 1] + (S == 'S' ? (_local * infoSplit[k][0] * _local + _local * (_local + 1) / 2) : ((_local * infoSplit[k][0] + _local) * _local))) : infoSplit[k][2]);
        }
        MPI_Request rq;
        MPI_Igatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, C, counts + _sizeSplit, counts, Wrapper<K>::mpi_type(), 0, _scatterComm, &rq);
#endif
        if(U != 1) {
#if !HPDDM_ICOLLECTIVE
            for(unsigned short k = 1; k < _sizeSplit; ++k) {
                if(infoSplit[k][2])
                    MPI_Irecv(C + offsetIdx[k - 1], infoSplit[k][2], Wrapper<K>::mpi_type(), k, 3, _scatterComm, rqRecv + idx + k - 1);
                else
                    rqRecv[idx + k - 1] = MPI_REQUEST_NULL;
            }
#endif
            offsetPosition = new unsigned int[_sizeSplit];
            offsetPosition[0] = std::accumulate(infoWorld, infoWorld + rankRelative, static_cast<unsigned int>(Solver<K>::_numbering == 'F'));
            if(T == 0 || T == 2)
                for(unsigned int k = 1; k < _sizeSplit; ++k)
                    offsetPosition[k] = offsetPosition[k - 1] + infoSplit[k - 1][1];
            else if(T == 1)
                for(unsigned int k = 1; k < _sizeSplit; ++k)
                    offsetPosition[k] = offsetPosition[k - 1] + infoWorld[rankRelative + k - 1];
        }
#if !HPDDM_ICOLLECTIVE
        else {
            for(unsigned short k = 1; k < _sizeSplit; ++k)
                MPI_Irecv(C + offsetIdx[k - 1], S == 'S' ? _local * infoSplit[k][0] * _local + _local * (_local + 1) / 2 : (_local * infoSplit[k][0] + _local) * _local, Wrapper<K>::mpi_type(), k, 3, _scatterComm, rqRecv + idx + k - 1);
        }
#endif
#pragma omp parallel for shared(I, J, infoWorld, infoSplit, rankRelative, offsetIdx, offsetPosition) schedule(dynamic, 64)
        for(unsigned int k = 1; k < _sizeSplit; ++k) {
            if(U == 1 || infoSplit[k][2]) {
                unsigned int tmp = U == 1 ? (rankRelative + k - (excluded == 2 ? (T == 1 ? p : 1 + rank) : 0)) * _local + (Solver<K>::_numbering == 'F') : offsetPosition[k];
                unsigned int offsetSlave = static_cast<unsigned int>(Solver<K>::_numbering == 'F');
                if(U != 1 && infoSplit[k][0])
                    offsetSlave = std::accumulate(infoWorld, infoWorld + infoSplit[k][3], offsetSlave);
                unsigned short i = 0;
                int* colIdx = J + offsetIdx[k - 1];
                if(S != 'S')
                    while(i < infoSplit[k][0] && infoSplit[k][(U != 1 ? 3 : 1) + i] < rankRelative + k - (U == 1 && excluded == 2 ? (T == 1 ? p : 1 + rank) : 0)) {
                        if(U != 1) {
                            if(i > 0)
                                offsetSlave = std::accumulate(infoWorld + infoSplit[k][2 + i], infoWorld + infoSplit[k][3 + i], offsetSlave);
                        }
                        else
                            offsetSlave = infoSplit[k][1 + i] * _local + (Solver<K>::_numbering == 'F');
                        std::iota(colIdx, colIdx + (U == 1 ? _local : infoWorld[infoSplit[k][3 + i]]), offsetSlave);
                        colIdx += (U == 1 ? _local : infoWorld[infoSplit[k][3 + i]]);
                        ++i;
                    }
                std::iota(colIdx, colIdx + (U == 1 ? _local : infoSplit[k][1]), tmp);
                colIdx += (U == 1 ? _local : infoSplit[k][1]);
                while(i < infoSplit[k][0]) {
                    if(U != 1) {
                        if(i > 0)
                            offsetSlave = std::accumulate(infoWorld + infoSplit[k][2 + i], infoWorld + infoSplit[k][3 + i], offsetSlave);
                    }
                    else
                        offsetSlave = infoSplit[k][1 + i] * _local + (Solver<K>::_numbering == 'F');
                    std::iota(colIdx, colIdx + (U == 1 ? _local : infoWorld[infoSplit[k][3 + i]]), offsetSlave);
                    colIdx += (U == 1 ? _local : infoWorld[infoSplit[k][3 + i]]);
                    ++i;
                }
#ifndef HPDDM_CSR_CO
                int* rowIdx = I + std::distance(J, colIdx);
                std::fill(I + offsetIdx[k - 1], rowIdx, tmp);
#else
                offsetSlave = (U == 1 ? (k - (excluded == 2)) * _local : offsetPosition[k] - offsetPosition[1] + (excluded == 2 ? 0 : _local));
                I[offsetSlave + 1] = colIdx - J - offsetIdx[k - 1];
#if defined(HPDDM_LOC2GLOB) && !defined(HPDDM_CONTIGUOUS)
                loc2glob[offsetSlave] = tmp;
#endif
#endif
                unsigned int coefficientsSlave = colIdx - J - offsetIdx[k - 1];
                for(i = 1; i < (U == 1 ? _local : infoSplit[k][1]); ++i) {
                    if(S == 'S')
                        --coefficientsSlave;
#ifndef HPDDM_CSR_CO
                    std::fill_n(rowIdx, coefficientsSlave, tmp + i);
                    rowIdx += coefficientsSlave;
#else
                    I[offsetSlave + 1 + i] = coefficientsSlave;
#if defined(HPDDM_LOC2GLOB) && !defined(HPDDM_CONTIGUOUS)
                    loc2glob[offsetSlave + i] = tmp + i;
#endif
#endif
                    std::copy(colIdx - coefficientsSlave, colIdx, colIdx);
                    colIdx += coefficientsSlave;
                }
#if defined(HPDDM_LOC2GLOB) && defined(HPDDM_CONTIGUOUS)
                if(excluded == 2 && k == 1)
                    loc2glob[0] = tmp;
                if(k == _sizeSplit - 1)
                    loc2glob[1] = tmp + (U == 1 ? _local : infoSplit[k][1]) - 1;
#endif
            }
        }
        delete [] offsetIdx;
        if(excluded < 2) {
#ifdef HPDDM_CSR_CO
            for(unsigned short k = 0; k < _local; ++k) {
                I[k + 1] = coefficients + (S == 'S' ? _local - k : 0);
#if defined(HPDDM_LOC2GLOB) && !defined(HPDDM_CONTIGUOUS)
                loc2glob[k] = offset + k;
#endif
            }
#if defined(HPDDM_LOC2GLOB) && defined(HPDDM_CONTIGUOUS)
            loc2glob[0] = offset;
            if(_sizeSplit == 1)
                loc2glob[1] = offset + _local - 1;
#endif
#endif
            unsigned int **offsetArray = new unsigned int*[info[0]];
            *offsetArray = new unsigned int[info[0] * ((Operator::_pattern == 's') + (U != 1))];
            if(Operator::_pattern == 's') {
                if(S != 'S') {
                    offsetArray[0][0] = sparsity[0] > _rankWorld ? _local : 0;
                    if(U != 1)
                        offsetArray[0][1] = std::accumulate(infoWorld, infoWorld + sparsity[0], static_cast<unsigned int>(Solver<K>::_numbering == 'F'));
                }
                else {
                    if(info[0] > 0) {
                        offsetArray[0][0] = _local;
                        if(U != 1)
                            offsetArray[0][1] = std::accumulate(infoWorld, infoWorld + sparsity[first], static_cast<unsigned int>(Solver<K>::_numbering == 'F'));
                    }
                }
                for(unsigned short k = 1; k < info[0]; ++k) {
                    if(U != 1) {
                        offsetArray[k] = *offsetArray + 2 * k;
                        offsetArray[k][1] = std::accumulate(infoWorld + sparsity[first + k - 1], infoWorld + sparsity[first + k], offsetArray[k - 1][1]);
                        offsetArray[k][0] = offsetArray[k - 1][0] + infoNeighbor[k - 1 + first];
                    }
                    else {
                        offsetArray[k] = *offsetArray + k;
                        offsetArray[k][0] = offsetArray[k - 1][0] + _local;
                    }
                    if((S != 'S') && sparsity[k - 1] < _rankWorld && sparsity[k] > _rankWorld)
                        offsetArray[k][0] += _local;
                }
            }
            else {
                if(U != 1) {
                    if(S != 'S')
                        offsetArray[0][0] = std::accumulate(infoWorld, infoWorld + sparsity[0], static_cast<unsigned int>(Solver<K>::_numbering == 'F'));
                    else if(info[0] > 0)
                        offsetArray[0][0] = std::accumulate(infoWorld, infoWorld + sparsity[first], static_cast<unsigned int>(Solver<K>::_numbering == 'F'));
                    for(unsigned short k = 1; k < info[0]; ++k) {
                        offsetArray[k] = *offsetArray + k;
                        offsetArray[k][0] = std::accumulate(infoWorld + sparsity[first + k - 1], infoWorld + sparsity[first + k], offsetArray[k - 1][0]);
                    }
                }
                info[0] = M.size();
            }
            if(U == 1 || _local)
                for(unsigned int k = 0; k < info[0]; ++k) {
                    int index;
                    MPI_Waitany(info[0], rqRecv, &index, MPI_STATUS_IGNORE);
                    if(Operator::_pattern == 's')
                        v.template applyFromNeighborMaster<S, Solver<K>::_numbering, U == 1>(recvNeighbor[index], index + first, I + offsetArray[index][0], J + offsetArray[index][0], C + offsetArray[index][0], coefficients + (S == 'S') * (_local - 1), offset, U == 1 ? nullptr : (offsetArray[index] + 1), work, U == 1 ? nullptr : infoNeighbor + first + index);
                    else
                        v.template applyFromNeighborMaster<S, Solver<K>::_numbering, U == 1>(recvNeighbor[index], index, I, J, C, coefficients, offset, U == 1 ? nullptr : *offsetArray, work, U == 1 ? nullptr : infoNeighbor);
                }
            delete [] *offsetArray;
            delete [] offsetArray;
        }
        delete [] info;
#if !HPDDM_ICOLLECTIVE
        MPI_Waitall(_sizeSplit - 1, rqRecv + idx, MPI_STATUSES_IGNORE);
#else
        MPI_Wait(&rq, MPI_STATUS_IGNORE);
#endif
#if HPDDM_ICOLLECTIVE
        delete [] counts;
#endif
        if(U != 1) {
            delete [] infoNeighbor;
            delete [] offsetPosition;
        }
        delete [] work;
        std::string filename = opt.prefix("master_filename", true);
        if(filename.size() > 0) {
            if(excluded == 2)
                filename += "_excluded";
            std::ofstream output { filename + "_" + S + "_" + Solver<K>::_numbering + "_" + to_string(T) + "_" + to_string(Solver<K>::_rank) + ".txt" };
#ifndef HPDDM_CSR_CO
            for(unsigned int i = 0; i < size; ++i)
                output << std::setw(9) << I[i] + (Solver<K>::_numbering == 'C') << std::setw(9) << J[i] + (Solver<K>::_numbering == 'C') << " " << std::scientific << C[i] << std::endl;
#else
            unsigned int accumulate = 0;
            for(unsigned int i = 0; i < nrow; ++i) {
                accumulate += I[i];
                for(unsigned int j = 0; j < I[i + 1]; ++j)
                    output << std::setw(9) <<
#ifndef HPDDM_LOC2GLOB
                    i + 1 <<
#elif !defined(HPDDM_CONTIGUOUS)
                    loc2glob[i] + (Solver<K>::_numbering == 'C') <<
#else
                    loc2glob[0] + i + (Solver<K>::_numbering == 'C') <<
#endif
                    std::setw(9) << J[accumulate + j - (Solver<K>::_numbering == 'F')] + (Solver<K>::_numbering == 'C') << " " << std::scientific << C[accumulate + j - (Solver<K>::_numbering == 'F')] << std::endl;
            }
#endif
        }
#ifdef HPDDM_CSR_CO
#ifndef DHYPRE
        std::partial_sum(I, I + nrow + 1, I);
#endif
#ifndef HPDDM_LOC2GLOB
        Solver<K>::template numfact<S>(nrow, I, J, C);
#else
        Solver<K>::template numfact<S>(nrow, I, loc2glob, J, C);
#endif
#else
        Solver<K>::template numfact<S>(size, I, J, C);
#endif

#ifdef DMKL_PARDISO
        if(S == 'S' || p != 1)
            delete [] C;
#else
        delete [] C;
#endif
        delete [] rqRecv;
    }
    if(excluded < 2)
        delete [] *sendNeighbor;
    if(U != 2) {
        switch(Solver<K>::_distribution) {
            case DMatrix::CENTRALIZED:
                _scatterComm = _gatherComm;
                break;
            case DMatrix::DISTRIBUTED_SOL:
                break;
            case DMatrix::DISTRIBUTED_SOL_AND_RHS:
                _gatherComm = _scatterComm;
                break;
        }
    }
    else {
        unsigned short* pt;
        unsigned short size;
        switch(Solver<K>::_distribution) {
            case DMatrix::CENTRALIZED:
                if(rankSplit != 0)
                    infoWorld = new unsigned short[_sizeWorld];
                pt = infoWorld;
                size = _sizeWorld;
                break;
            case DMatrix::DISTRIBUTED_SOL:
                size = _sizeWorld + _sizeSplit;
                pt = new unsigned short[size];
                if(rankSplit == 0) {
                    std::copy_n(infoWorld, _sizeWorld, pt);
                    for(unsigned int i = 0; i < _sizeSplit; ++i)
                        pt[_sizeWorld + i] = infoSplit[i][1];
                }
                break;
            case DMatrix::DISTRIBUTED_SOL_AND_RHS:
                unsigned short* infoMaster;
                if(rankSplit == 0) {
                    infoMaster = infoSplit[0];
                    for(unsigned int i = 0; i < _sizeSplit; ++i)
                        infoMaster[i] = infoSplit[i][1];
                }
                else
                    infoMaster = new unsigned short[_sizeSplit];
                pt = infoMaster;
                size = _sizeSplit;
                break;
        }
        MPI_Bcast(pt, size, MPI_UNSIGNED_SHORT, 0, _scatterComm);
        if(Solver<K>::_distribution == DMatrix::CENTRALIZED || Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL_AND_RHS) {
            if(Solver<K>::_distribution == DMatrix::CENTRALIZED)
                constructionCommunicatorCollective<(excluded > 0)>(pt, size, _gatherComm, &_scatterComm);
            else
                constructionCommunicatorCollective<false>(pt, size, _scatterComm);
            _gatherComm = _scatterComm;
        }
        else if(Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL) {
            constructionCommunicatorCollective<(excluded > 0)>(pt, _sizeWorld, _gatherComm);
            constructionCommunicatorCollective<false>(pt + _sizeWorld, _sizeSplit, _scatterComm);
        }
        if(rankSplit != 0 || Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL)
            delete [] pt;
    }
    if(rankSplit == 0) {
        if(Solver<K>::_distribution == DMatrix::CENTRALIZED) {
            if(_rankWorld == 0) {
                _sizeRHS = Solver<K>::_n;
                if(U == 1)
                    constructionCollective<true, DMatrix::CENTRALIZED, excluded == 2>();
                else if(U == 2) {
                    Solver<K>::_gatherCounts = new int[1];
                    if(_local == 0) {
                        _local = *Solver<K>::_gatherCounts = *std::find_if(infoWorld, infoWorld + _sizeWorld, [](const unsigned short& nu) { return nu != 0; });
                        _sizeRHS += _local;
                    }
                    else
                        *Solver<K>::_gatherCounts = _local;
                }
                else
                    constructionCollective<false, DMatrix::CENTRALIZED, excluded == 2>(infoWorld, p - 1);
            }
            else {
                if(U == 0)
                    Solver<K>::_displs = &_rankWorld;
                _sizeRHS = _local;
            }
        }
        else {
            constructionMap<T, U == 1, excluded == 2>(p, U == 1 ? nullptr : infoWorld);
            if(Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL) {
                if(_rankWorld == 0)
                    _sizeRHS = Solver<K>::_n;
                else
                    _sizeRHS = Solver<K>::_ldistribution[Solver<K>::_rank];
            }
            else
                _sizeRHS = Solver<K>::_ldistribution[Solver<K>::_rank];
            if(U == 1)
                constructionCollective<true, DMatrix::DISTRIBUTED_SOL, excluded == 2>();
            else if(U == 2) {
                Solver<K>::_gatherCounts = new int[1];
                if(_local == 0) {
                    _local = *Solver<K>::_gatherCounts = *std::find_if(infoWorld, infoWorld + _sizeWorld, [](const unsigned short& nu) { return nu != 0; });
                    _sizeRHS += _local;
                }
                else
                    *Solver<K>::_gatherCounts = _local;
            }
            else {
                unsigned short* infoMaster = infoSplit[0];
                for(unsigned int i = 0; i < _sizeSplit; ++i)
                    infoMaster[i] = infoSplit[i][1];
                constructionCollective<false, DMatrix::DISTRIBUTED_SOL, excluded == 2>(infoWorld, p - 1, infoMaster);
            }
        }
        delete [] *infoSplit;
        delete [] infoSplit;
        if(excluded == 2) {
            if(Solver<K>::_distribution == DMatrix::CENTRALIZED && _rankWorld == 0)
                _sizeRHS += _local;
            else if(Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL || Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL_AND_RHS)
                _sizeRHS += _local;
        }
    }
    return ret;
}

template<template<class> class Solver, char S, class K>
template<bool excluded>
inline void CoarseOperator<Solver, S, K>::callSolver(K* const rhs, const unsigned short& mu) {
    if(_scatterComm != MPI_COMM_NULL) {
        if(Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL) {
            if(Solver<K>::_displs) {
                if(_rankWorld == 0) {
                    transfer<false>(Solver<K>::_gatherCounts, _sizeWorld, mu, rhs);
                    std::for_each(Solver<K>::_gatherCounts, Solver<K>::_displs + _sizeWorld, [&](int& i) { i /= mu; });
                }
                else if(_gatherComm != MPI_COMM_NULL)
                    MPI_Gatherv(rhs, mu * _local, Wrapper<K>::mpi_type(), NULL, 0, 0, Wrapper<K>::mpi_type(), 0, _gatherComm);
                if(Solver<K>::_communicator != MPI_COMM_NULL) {
                    Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL>(rhs, mu);
                    std::for_each(Solver<K>::_gatherSplitCounts, Solver<K>::_displsSplit + _sizeSplit, [&](int& i) { i *= mu; });
                    transfer<true>(Solver<K>::_gatherSplitCounts, mu, _sizeSplit, rhs);
                }
                else
                    MPI_Scatterv(NULL, 0, 0, Wrapper<K>::mpi_type(), rhs, mu * _local, Wrapper<K>::mpi_type(), 0, _scatterComm);
            }
            else {
                if(_rankWorld == 0) {
                    MPI_Gather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, mu * *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), 0, _gatherComm);
                    Wrapper<K>::template cycle<'T'>(_sizeWorld, mu, rhs, *Solver<K>::_gatherCounts);
                }
                else if(_gatherComm != MPI_COMM_NULL)
                    MPI_Gather(rhs, mu * _local, Wrapper<K>::mpi_type(), NULL, 0, MPI_DATATYPE_NULL, 0, _gatherComm);
                if(Solver<K>::_communicator != MPI_COMM_NULL) {
                    Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL>(rhs + (_offset || excluded ? *Solver<K>::_gatherCounts : 0), mu);
                    Wrapper<K>::template cycle<'T'>(mu, _sizeSplit, rhs, *Solver<K>::_gatherCounts);
                    MPI_Scatter(rhs, mu * *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 0, _scatterComm);
                }
                else
                    MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL, rhs, mu * _local, Wrapper<K>::mpi_type(), 0, _scatterComm);
            }
        }
        else if(Solver<K>::_distribution == DMatrix::CENTRALIZED) {
            if(Solver<K>::_displs) {
                if(_rankWorld == 0)
                    transfer<false>(Solver<K>::_gatherCounts, _sizeWorld, mu, rhs);
                else if(_gatherComm != MPI_COMM_NULL)
                    MPI_Gatherv(rhs, mu * _local, Wrapper<K>::mpi_type(), NULL, 0, 0, Wrapper<K>::mpi_type(), 0, _gatherComm);
                if(Solver<K>::_communicator != MPI_COMM_NULL)
                    Solver<K>::template solve<DMatrix::CENTRALIZED>(rhs, mu);
                if(_rankWorld == 0)
                    transfer<true>(Solver<K>::_gatherCounts, mu, _sizeWorld, rhs);
                else if(_gatherComm != MPI_COMM_NULL)
                    MPI_Scatterv(NULL, 0, 0, Wrapper<K>::mpi_type(), rhs, mu * _local, Wrapper<K>::mpi_type(), 0, _gatherComm);
            }
            else {
                if(_rankWorld == 0) {
                    MPI_Gather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, mu * *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), 0, _gatherComm);
                    Wrapper<K>::template cycle<'T'>(_sizeWorld, mu, rhs, *Solver<K>::_gatherCounts);
                }
                else
                    MPI_Gather(rhs, mu * _local, Wrapper<K>::mpi_type(), NULL, 0, MPI_DATATYPE_NULL, 0, _gatherComm);
                if(Solver<K>::_communicator != MPI_COMM_NULL)
                    Solver<K>::template solve<DMatrix::CENTRALIZED>(rhs + (_offset || excluded ? _local : 0), mu);
                if(_rankWorld == 0) {
                    Wrapper<K>::template cycle<'T'>(mu, _sizeWorld, rhs, *Solver<K>::_gatherCounts);
                    MPI_Scatter(rhs, mu * *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 0, _scatterComm);
                }
                else
                    MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL, rhs, mu * _local, Wrapper<K>::mpi_type(), 0, _scatterComm);
            }
        }
        else if(Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL_AND_RHS) {
            if(Solver<K>::_displs) {
                if(Solver<K>::_communicator != MPI_COMM_NULL) {
                    transfer<false>(Solver<K>::_gatherSplitCounts, _sizeSplit, mu, rhs);
                    Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL_AND_RHS>(rhs, mu);
                    transfer<true>(Solver<K>::_gatherSplitCounts, mu, _sizeSplit, rhs);
                }
                else {
                    MPI_Gatherv(rhs, mu * _local, Wrapper<K>::mpi_type(), NULL, 0, 0, Wrapper<K>::mpi_type(), 0, _gatherComm);
                    MPI_Scatterv(NULL, 0, 0, Wrapper<K>::mpi_type(), rhs, mu * _local, Wrapper<K>::mpi_type(), 0, _scatterComm);
                }
            }
            else {
                    if(Solver<K>::_communicator != MPI_COMM_NULL) {
                        MPI_Gather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, mu * *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), 0, _gatherComm);
                        Wrapper<K>::template cycle<'T'>(_sizeSplit, mu, rhs, *Solver<K>::_gatherCounts);
                        Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL_AND_RHS>(rhs + (_offset || excluded ? *Solver<K>::_gatherCounts : 0), mu);
                        Wrapper<K>::template cycle<'T'>(mu, _sizeSplit, rhs, *Solver<K>::_gatherCounts);
                        MPI_Scatter(rhs, mu * *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 0, _scatterComm);
                    }
                    else {
                        MPI_Gather(rhs, mu * _local, Wrapper<K>::mpi_type(), NULL, 0, MPI_DATATYPE_NULL, 0, _gatherComm);
                        MPI_Scatter(NULL, 0, MPI_DATATYPE_NULL, rhs, mu * _local, Wrapper<K>::mpi_type(), 0, _scatterComm);
                    }
            }
        }
    }
    else if(Solver<K>::_communicator != MPI_COMM_NULL) {
        switch(Solver<K>::_distribution) {
            case DMatrix::CENTRALIZED:             Solver<K>::template solve<DMatrix::CENTRALIZED>(rhs, mu); break;
            case DMatrix::DISTRIBUTED_SOL:         Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL>(rhs, mu); break;
            case DMatrix::DISTRIBUTED_SOL_AND_RHS: Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL_AND_RHS>(rhs, mu); break;
        }
    }
}

#if HPDDM_ICOLLECTIVE
template<template<class> class Solver, char S, class K>
template<bool excluded>
inline void CoarseOperator<Solver, S, K>::IcallSolver(K* const rhs, MPI_Request* rq) {
    if(_scatterComm != MPI_COMM_NULL) {
        if(Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL) {
            if(Solver<K>::_displs) {
                if(_rankWorld == 0)                   MPI_Igatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, Solver<K>::_gatherCounts, Solver<K>::_displs, Wrapper<K>::mpi_type(), 0, _gatherComm, rq);
                else if(_gatherComm != MPI_COMM_NULL) MPI_Igatherv(rhs, _local, Wrapper<K>::mpi_type(), NULL, 0, 0, MPI_DATATYPE_NULL, 0, _gatherComm, rq);
                if(Solver<K>::_communicator != MPI_COMM_NULL) {
                    MPI_Wait(rq, MPI_STATUS_IGNORE);
                    Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL>(rhs);
                    MPI_Iscatterv(rhs, Solver<K>::_gatherSplitCounts, Solver<K>::_displsSplit, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 0, _scatterComm, rq + 1);
                }
                else
                    MPI_Iscatterv(NULL, 0, 0, MPI_DATATYPE_NULL, rhs, _local, Wrapper<K>::mpi_type(), 0, _scatterComm, rq + 1);
            }
            else {
                if(_rankWorld == 0)                   MPI_Igather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), 0, _gatherComm, rq);
                else if(_gatherComm != MPI_COMM_NULL) MPI_Igather(rhs, _local, Wrapper<K>::mpi_type(), NULL, 0, MPI_DATATYPE_NULL, 0, _gatherComm, rq);
                if(Solver<K>::_communicator != MPI_COMM_NULL) {
                    MPI_Wait(rq, MPI_STATUS_IGNORE);
                    Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL>(rhs + (_offset || excluded ? *Solver<K>::_gatherCounts : 0));
                    MPI_Iscatter(rhs, *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 0, _scatterComm, rq + 1);
                }
                else
                    MPI_Iscatter(NULL, 0, MPI_DATATYPE_NULL, rhs, _local, Wrapper<K>::mpi_type(), 0, _scatterComm, rq + 1);
            }
        }
        else if(Solver<K>::_distribution == DMatrix::CENTRALIZED) {
            if(Solver<K>::_displs) {
                if(_rankWorld == 0)                   MPI_Igatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, Solver<K>::_gatherCounts, Solver<K>::_displs, Wrapper<K>::mpi_type(), 0, _gatherComm, rq);
                else if(_gatherComm != MPI_COMM_NULL) MPI_Igatherv(rhs, _local, Wrapper<K>::mpi_type(), NULL, 0, 0, MPI_DATATYPE_NULL, 0, _gatherComm, rq);
                if(Solver<K>::_communicator != MPI_COMM_NULL) {
                    MPI_Wait(rq, MPI_STATUS_IGNORE);
                    Solver<K>::template solve<DMatrix::CENTRALIZED>(rhs);
                }
                if(_rankWorld == 0)                   MPI_Iscatterv(rhs, Solver<K>::_gatherCounts, Solver<K>::_displs, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 0, _scatterComm, rq + 1);
                else if(_gatherComm != MPI_COMM_NULL) MPI_Iscatterv(NULL, 0, 0, MPI_DATATYPE_NULL, rhs, _local, Wrapper<K>::mpi_type(), 0, _scatterComm, rq + 1);
            }
            else {
                if(_rankWorld == 0)                   MPI_Igather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), 0, _gatherComm, rq);
                else                                  MPI_Igather(rhs, _local, Wrapper<K>::mpi_type(), NULL, 0, MPI_DATATYPE_NULL, 0, _gatherComm, rq);
                if(Solver<K>::_communicator != MPI_COMM_NULL) {
                    MPI_Wait(rq, MPI_STATUS_IGNORE);
                    Solver<K>::template solve<DMatrix::CENTRALIZED>(rhs + (_offset || excluded ? _local : 0));
                }
                if(_rankWorld == 0)                   MPI_Iscatter(rhs, *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 0, _gatherComm, rq + 1);
                else                                  MPI_Iscatter(NULL, 0, MPI_DATATYPE_NULL, rhs, _local, Wrapper<K>::mpi_type(), 0, _gatherComm, rq + 1);
            }
        }
        else if(Solver<K>::_distribution == DMatrix::DISTRIBUTED_SOL_AND_RHS) {
            if(Solver<K>::_displs) {
                if(Solver<K>::_communicator != MPI_COMM_NULL) {
                    MPI_Igatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, Solver<K>::_gatherSplitCounts, Solver<K>::_displsSplit, Wrapper<K>::mpi_type(), 0, _gatherComm, rq);
                    MPI_Wait(rq, MPI_STATUS_IGNORE);
                    Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL_AND_RHS>(rhs);
                    MPI_Iscatterv(rhs, Solver<K>::_gatherSplitCounts, Solver<K>::_displsSplit, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 0, _scatterComm, rq + 1);
                }
                else {
                    MPI_Igatherv(rhs, _local, Wrapper<K>::mpi_type(), NULL, 0, 0, MPI_DATATYPE_NULL, 0, _gatherComm, rq);
                    MPI_Iscatterv(NULL, 0, 0, MPI_DATATYPE_NULL, rhs, _local, Wrapper<K>::mpi_type(), 0, _scatterComm, rq + 1);
                }
            }
            else {
                    if(Solver<K>::_communicator != MPI_COMM_NULL) {
                        MPI_Igather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, rhs, *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), 0, _gatherComm, rq);
                        MPI_Wait(rq, MPI_STATUS_IGNORE);
                        Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL_AND_RHS>(rhs + (_offset || excluded ? *Solver<K>::_gatherCounts : 0));
                        MPI_Iscatter(rhs, *Solver<K>::_gatherCounts, Wrapper<K>::mpi_type(), MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 0, _scatterComm, rq + 1);
                    }
                    else {
                        MPI_Igather(rhs, _local, Wrapper<K>::mpi_type(), NULL, 0, MPI_DATATYPE_NULL, 0, _gatherComm, rq);
                        MPI_Iscatter(NULL, 0, MPI_DATATYPE_NULL, rhs, _local, Wrapper<K>::mpi_type(), 0, _scatterComm, rq + 1);
                    }
            }
        }
    }
    else if(Solver<K>::_communicator != MPI_COMM_NULL) {
        switch(Solver<K>::_distribution) {
            case DMatrix::CENTRALIZED:             Solver<K>::template solve<DMatrix::CENTRALIZED>(rhs, mu); break;
            case DMatrix::DISTRIBUTED_SOL:         Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL>(rhs, mu); break;
            case DMatrix::DISTRIBUTED_SOL_AND_RHS: Solver<K>::template solve<DMatrix::DISTRIBUTED_SOL_AND_RHS>(rhs, mu); break;
        }
    }
}
#endif // HPDDM_ICOLLECTIVE
} // HPDDM
#endif // _HPDDM_COARSE_OPERATOR_IMPL_
