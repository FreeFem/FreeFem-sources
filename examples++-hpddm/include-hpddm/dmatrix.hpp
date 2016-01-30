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

#ifndef _HPDDM_DMATRIX_
#define _HPDDM_DMATRIX_

#include <map>

namespace HPDDM {
/* Class: DMatrix
 *  A class for handling all communications and computations involving a distributed matrix. */
class DMatrix {
    public:
        /* Enum: Distribution
         *
         *  Defines the distribution of both right-hand sides and solution vectors.
         *
         * CENTRALIZED             - Neither are distributed, both are centralized on the root of <DMatrix::communicator>.
         * DISTRIBUTED_SOL         - Right-hand sides are centralized, while solution vectors are distributed on <DMatrix::communicator>.
         * DISTRIBUTED_SOL_AND_RHS - Both are distributed on <DMatrix::communicator>. */
        enum Distribution : char {
            CENTRALIZED, DISTRIBUTED_SOL, DISTRIBUTED_SOL_AND_RHS
        };
        /* Function: splitCommunicator
         *
         *  If requested, splits a communicator into one made of master processes and another one made of slave processes.
         *
         * Parameters:
         *    in             - Original communicator.
         *    out            - Output communicator which may be left untouched.
         *    exclude        - True if the master processes have to be excluded from the original communicator.
         *    p              - Number of master processes.
         *    T              - Master processes distribution topology. */
        static bool splitCommunicator(const MPI_Comm& in, MPI_Comm& out, const bool& exclude, unsigned short& p, const unsigned short& T) {
            int size, rank;
            MPI_Comm_size(in, &size);
            MPI_Comm_rank(in, &rank);
            if(p > size / 2 && size > 1) {
                p = size / 2;
                if(rank == 0)
                    std::cout << "WARNING -- the number of master processes was set to a value greater than MPI_Comm_size, the value has been reset to " << p << std::endl;
            }
            else if(p < 1)
                p = 1;
            if(exclude) {
                MPI_Group oldGroup, newGroup;
                MPI_Comm_group(in, &oldGroup);
                int* pm = new int[p];
                if(T == 1)
                    std::iota(pm, pm + p, 0);
                else if(T == 2) {
                    float area = size * size / (2.0 * p);
                    *pm = 0;
                    for(unsigned short i = 1; i < p; ++i)
                        pm[i] = static_cast<int>(size - std::sqrt(std::max(size * size - 2 * size * pm[i - 1] - 2 * area + pm[i - 1] * pm[i - 1], 1.0f)) + 0.5);
                }
                else
                    for(unsigned short i = 0; i < p; ++i)
                        pm[i] = i * (size / p);
                bool excluded = std::binary_search(pm, pm + p, rank);
                if(excluded)
                    MPI_Group_incl(oldGroup, p, pm, &newGroup);
                else
                    MPI_Group_excl(oldGroup, p, pm, &newGroup);
                MPI_Comm_create(in, newGroup, &out);
                MPI_Group_free(&oldGroup);
                MPI_Group_free(&newGroup);
                delete [] pm;
                return excluded;
            }
            else {
                MPI_Comm_dup(in, &out);
                return false;
            }
        }
    protected:
        /* Typedef: pair_type
         *  std::pair of unsigned integers. */
        typedef std::pair<unsigned int, unsigned int>   pair_type;
        /* Typedef: map_type
         *
         *  std::map of std::vector<T> indexed by unsigned short integers.
         *
         * Template Parameter:
         *    T              - Class. */
        template<class T>
        using map_type = std::map<unsigned short, std::vector<T>>;
        /* Variable: mapRecv
         *  Values that have to be received to match the distribution of the direct solver and of the user. */
        map_type<pair_type>*        _mapRecv;
        /* Variable: mapSend
         *  Values that have to be sent to match the distribution of the direct solver and of the user. */
        map_type<pair_type>*        _mapSend;
        /* Variable: mapOwn
         *  Values that have to remain on this process to match the distribution of the direct solver and of the user. */
        std::vector<pair_type>*      _mapOwn;
        /* Variable: ldistribution
         *  User distribution. */
        int*                  _ldistribution;
        /* Variable: idistribution */
        int*                  _idistribution;
        /* Variable: gatherCounts */
        int*                   _gatherCounts;
        /* Variable: gatherSplitCounts */
        int*              _gatherSplitCounts;
        /* Variable: displs */
        int*                         _displs;
        /* Variable: displsSplit */
        int*                    _displsSplit;
        /* Variable: communicator
         *  MPI communicator on which the matrix is distributed. */
        MPI_Comm               _communicator;
        /* Variable: n
         *  Size of the coarse operator. */
        int                               _n;
        /* Variable: rank
         *  Rank of the current master process in <Coarse operator::communicator>. */
        int                            _rank;
        /* Variable: distribution
         *  <Distribution> used for right-hand sides and solution vectors. */
        Distribution           _distribution;
        /* Function: initializeMap
         *
         *  Initializes <Coarse operator::mapRecv>, <Coarse operator::mapSend>, and <Coarse operator::mapOwn>.
         *
         * Template Parameter:
         *    isRHS          - True if this function is first called to redistribute a right-hand side, false otherwise.
         *
         * Parameters:
         *    info           - Local dimension returned by the direct solver.
         *    isol_loc       - Numbering of the direct solver.
         *    sol_loc        - Vector following the numbering of the direct solver.
         *    sol            - Vector following the numbering of the user. */
        template<bool isRHS, class K>
        void initializeMap(const int& info, const int* const isol_loc, K* const sol_loc, K* const sol) {
            int size;
            MPI_Comm_size(_communicator, &size);
            int* lsol_loc_glob = new int[size];
            lsol_loc_glob[_rank] = info;
            MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, lsol_loc_glob, 1, MPI_INT, _communicator);
            int* isol_loc_glob = new int[_n];
            int* disp_lsol_loc_glob = new int[size];
            disp_lsol_loc_glob[0] = 0;
            std::partial_sum(lsol_loc_glob, lsol_loc_glob + size - 1, disp_lsol_loc_glob + 1);
            MPI_Allgatherv(const_cast<int*>(isol_loc), info, MPI_INT, isol_loc_glob, lsol_loc_glob, disp_lsol_loc_glob, MPI_INT, _communicator);
            delete [] disp_lsol_loc_glob;
            std::vector<std::pair<unsigned short, unsigned int>> mapping(_n);
            std::vector<std::pair<unsigned short, unsigned int>> mapping_user(_n);
            _mapRecv = new map_type<pair_type>;
            _mapSend = new map_type<pair_type>;
            _mapOwn = new std::vector<pair_type>;
            unsigned int offset = 0;
            for(unsigned short i = 0; i < size; ++i)
                for(unsigned int j = 0; j < lsol_loc_glob[i]; ++j)
                    mapping[isol_loc_glob[offset++] - 1] = { i, j };
            offset = 0;
            if(_idistribution) {
                for(unsigned short i = 0; i < size; ++i)
                    for(unsigned int j = 0; j < _ldistribution[i]; ++j)
                        mapping_user[_idistribution[offset++]] = { i, j };
            }
            else {
                for(unsigned short i = 0; i < size; ++i)
                    for(unsigned int j = 0; j < _ldistribution[i]; ++j)
                        mapping_user[offset++] = { i, j };
            }
            map_type<K> map_recv;
            map_type<K> map_send;
            offset = std::accumulate(_ldistribution, _ldistribution + _rank, 0);
            if(!isRHS) {
                for(unsigned int i = 0; i < info; ++i) {
                    std::pair<unsigned short, unsigned int> tmp = mapping_user[isol_loc[i] - 1];
                    if(tmp.first != _rank) {
                        map_send[tmp.first].emplace_back(sol_loc[i]);
                        (*_mapSend)[tmp.first].emplace_back(i, tmp.second);
                    }
                }
                if(_idistribution)
                    for(unsigned int i = offset; i < offset + _ldistribution[_rank]; ++i) {
                        std::pair<unsigned short, unsigned int> tmp = mapping[_idistribution[i]];
                        if(tmp.first != _rank)
                            map_recv[tmp.first].resize(map_recv[tmp.first].size() + 1);
                        else {
                            _mapOwn->emplace_back(mapping_user[_idistribution[i]].second, tmp.second);
                            sol[mapping_user[_idistribution[i]].second] = sol_loc[tmp.second];
                        }
                    }
                else
                    for(unsigned int i = offset; i < offset + _ldistribution[_rank]; ++i) {
                        std::pair<unsigned short, unsigned int> tmp = mapping[i];
                        if(tmp.first != _rank)
                            map_recv[tmp.first].resize(map_recv[tmp.first].size() + 1);
                        else {
                            _mapOwn->emplace_back(mapping_user[i].second, tmp.second);
                            sol[mapping_user[i].second] = sol_loc[tmp.second];
                        }
                    }
            }
            else {
                for(unsigned int i = 0; i < info; ++i) {
                    unsigned short tmp = mapping_user[isol_loc[i] - 1].first;
                    if(tmp != _rank)
                        map_recv[tmp].resize(map_recv[tmp].size() + 1);
                }
                if(_idistribution)
                    for(unsigned int i = offset; i < offset + _ldistribution[_rank]; ++i) {
                        std::pair<unsigned short, unsigned int> tmp = mapping[_idistribution[i]];
                        if(tmp.first != _rank) {
                            map_send[tmp.first].emplace_back(sol[i - offset]);
                            (*_mapSend)[tmp.first].emplace_back(i - offset, tmp.second);
                        }
                        else {
                            _mapOwn->emplace_back(mapping[_idistribution[i]].second, mapping_user[_idistribution[i]].second);
                            sol_loc[mapping[_idistribution[i]].second] = sol[mapping_user[_idistribution[i]].second];
                        }
                    }
                else
                    for(unsigned int i = offset; i < offset + _ldistribution[_rank]; ++i) {
                        std::pair<unsigned short, unsigned int> tmp = mapping[i];
                        if(tmp.first != _rank) {
                            map_send[tmp.first].emplace_back(sol[i - offset]);
                            (*_mapSend)[tmp.first].emplace_back(i - offset, tmp.second);
                        }
                        else {
                            _mapOwn->emplace_back(mapping[i].second, mapping_user[i].second);
                            sol_loc[mapping[i].second] = sol[mapping_user[i].second];
                        }
                    }
            }
            MPI_Request* rqSend = new MPI_Request[map_send.size() + map_recv.size()];
            MPI_Request* rqRecv = rqSend + map_send.size();
            unsigned int i = 0;
            for(typename map_type<K>::iterator it = map_send.begin(); it != map_send.end(); ++it)
                MPI_Isend(it->second.data(), it->second.size(), Wrapper<K>::mpi_type(), it->first, 4, _communicator, rqSend + i++);
            for(typename map_type<K>::iterator it = map_recv.begin(); it != map_recv.end(); ++it)
                MPI_Irecv(it->second.data(), it->second.size(), Wrapper<K>::mpi_type(), it->first, 4, _communicator, rqSend + i++);
            for(unsigned int i = 0; i < map_recv.size(); ++i) {
                int index;
                MPI_Waitany(map_recv.size(), rqRecv, &index, MPI_STATUS_IGNORE);
                typename map_type<K>::const_iterator it_index = std::next(map_recv.cbegin(), index);
                if(!isRHS) {
                    unsigned int offset = std::accumulate(lsol_loc_glob, lsol_loc_glob + it_index->first, 0);
                    for(unsigned int j = offset, accumulate = 0; j < offset + lsol_loc_glob[it_index->first]; ++j)
                        if(mapping_user[isol_loc_glob[j] - 1].first == _rank) {
                            (*_mapRecv)[it_index->first].emplace_back(mapping_user[isol_loc_glob[j] - 1].second, accumulate);
                            sol[mapping_user[isol_loc_glob[j] - 1].second] = (it_index->second)[accumulate++];
                        }
                }
                else {
                    unsigned int offset = std::accumulate(_ldistribution, _ldistribution + it_index->first, 0);
                    if(_idistribution) {
                        for(unsigned int j = offset, accumulate = 0; j < offset + _ldistribution[it_index->first]; ++j)
                            if(mapping[_idistribution[j]].first == _rank) {
                                (*_mapRecv)[it_index->first].emplace_back(mapping[_idistribution[j]].second, accumulate);
                                sol_loc[mapping[_idistribution[j]].second] = (it_index->second)[accumulate++];
                            }
                    }
                    else {
                        for(unsigned int j = offset, accumulate = 0; j < offset + _ldistribution[it_index->first]; ++j)
                            if(mapping[j].first == _rank) {
                                (*_mapRecv)[it_index->first].emplace_back(mapping[j].second, accumulate);
                                sol_loc[mapping[j].second] = (it_index->second)[accumulate++];
                            }
                    }
                }
            }
            MPI_Waitall(map_send.size(), rqSend, MPI_STATUSES_IGNORE);
            delete [] rqSend;
            delete [] isol_loc_glob;
            delete [] lsol_loc_glob;
            delete [] _ldistribution;
            delete [] _idistribution;
        }
        /* Function: redistribute
         *
         *  Transfers a vector numbered by the user to match the numbering of the direct solver, and vice versa.
         *
         * Template Parameter:
         *    P              - Renumbering identifier.
         *
         * Parameters:
         *    vec            - Vector following the numbering of the direct solver.
         *    res            - Vector following the numbering of the user.
         *    fuse           - Number of fused reductions (optional). */
        template<char P, class K>
        void redistribute(K* const vec, K* const res, const unsigned short& fuse = 0) {
            map_type<pair_type>* map_recv_index;
            map_type<pair_type>* map_send_index;
            if(P == 0 || P == 1) {
                map_recv_index = _mapRecv;
                map_send_index = _mapSend;
            }
            else {
                map_recv_index = _mapSend;
                map_send_index = _mapRecv;
            }
            map_type<K> map_recv;
            map_type<K> map_send;
            int gatherCount = fuse > 0 ? *_gatherCounts - fuse : 1;
            MPI_Request* rqSend = new MPI_Request[map_send_index->size() + map_recv_index->size()];
            MPI_Request* rqRecv = rqSend + map_send_index->size();
            unsigned short i = 0;
            for(map_type<pair_type>::const_reference q : *map_send_index) {
                map_send[q.first].reserve(q.second.size());
                for(std::vector<pair_type>::const_reference p : q.second) {
                    if(P == 0)
                        map_send[q.first].emplace_back(vec[p.first]);
                    else if(P == 1)
                        map_send[q.first].emplace_back(res[p.first + fuse * (p.first / gatherCount)]);
                    else if(P == 2)
                        map_send[q.first].emplace_back(res[p.first]);
                }
                MPI_Isend(map_send[q.first].data(), q.second.size(), Wrapper<K>::mpi_type(), q.first, 5, _communicator, rqSend + i++);
            }
            for(map_type<pair_type>::const_reference q : *map_recv_index) {
                map_recv[q.first].reserve(q.second.size());
                MPI_Irecv(map_recv[q.first].data(), q.second.size(), Wrapper<K>::mpi_type(), q.first, 5, _communicator, rqSend + i++);
            }
            for(std::vector<pair_type>::const_reference p : *_mapOwn) {
                if(P == 0)
                    res[p.first] = vec[p.second];
                else if(P == 1)
                    vec[p.first] = res[p.second + fuse * (p.second / gatherCount)];
                else if(P == 2)
                    vec[p.second + fuse * (p.second / gatherCount)] = res[p.first];
            }
            for(i = 0; i < map_recv.size(); ++i) {
                int index;
                MPI_Waitany(map_recv.size(), rqRecv, &index, MPI_STATUS_IGNORE);
                typename map_type<K>::const_iterator it_index = std::next(map_recv.cbegin(), index);
                K* pt = map_recv[it_index->first].data();
                for(std::vector<pair_type>::const_reference p : (*map_recv_index)[it_index->first]) {
                    if(P == 0)
                        res[p.first] = pt[p.second];
                    else if(P == 1)
                        vec[p.first] = pt[p.second];
                    else if(P == 2)
                        vec[p.first + fuse * (p.first / gatherCount)] = *pt++;
                }
            }
            MPI_Waitall(map_send.size(), rqSend, MPI_STATUSES_IGNORE);
            delete [] rqSend;
        }
        /* Function: initialize
         *
         *  Initializes <DMatrix::rank> and <DMatrix::distribution>.
         *
         * Parameters:
         *    solver         - Name of the solver.
         *    list           - Supported <DMatrix::Distribution>s. */
        void initialize(char const* solver, std::initializer_list<Distribution> list) {
            Option& opt = *Option::get();
            _distribution = static_cast<Distribution>(opt["master_distribution"]);
            std::initializer_list<Distribution>::iterator it = std::find(list.begin(), list.end(), _distribution);
            if(it == list.end()) {
                opt["master_distribution"] = _distribution = *list.begin();
                if(_communicator != MPI_COMM_NULL && _rank == 0) {
                    std::cout << "WARNING -- " << solver << " only supports the following distribution" << (list.size() > 1 ? "s" : "") << ":";
                    for(std::initializer_list<Distribution>::iterator it = list.begin(); it != list.end(); ++it) {
                        switch(*it) {
                            case CENTRALIZED: std::cout << " 'centralized'"; break;
                            case DISTRIBUTED_SOL: std::cout << " 'distributed sol'"; break;
                            case DISTRIBUTED_SOL_AND_RHS: std::cout << " 'distributed sol and rhs'"; break;
                        }
                        if(list.size() > 1 && it != list.end() - 1)
                            std::cout << ",";
                    }
                    switch(_distribution) {
                        case CENTRALIZED: std::cout << " (forcing the distribution to 'centralized')" << std::endl; break;
                        case DISTRIBUTED_SOL: std::cout << " (forcing the distribution to 'distributed sol')" << std::endl; break;
                        case DISTRIBUTED_SOL_AND_RHS: std::cout << " (forcing the distribution to 'distributed sol and rhs')" << std::endl; break;
                    }
                }
            }
        }
    public:
        DMatrix() : _mapRecv(), _mapSend(), _mapOwn(), _ldistribution(), _idistribution(), _gatherCounts(), _gatherSplitCounts(), _displs(), _displsSplit(), _communicator(MPI_COMM_NULL), _n(), _rank(), _distribution() { }
        DMatrix(const DMatrix&) = delete;
        ~DMatrix() {
            if(!_mapRecv) {
                delete [] _ldistribution;
                delete [] _idistribution;
            }
            delete _mapRecv;
            delete _mapSend;
            delete _mapOwn;
            delete [] _gatherCounts;
            delete [] _gatherSplitCounts;
        }
};
} // HPDDM
#endif // _HPDDM_DMATRIX_
