/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY  :   add interface with partionning library METSI and ParMETIS
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht, P. Jolivet
// E-MAIL  : frederic.hecht@sorbonne-universite.fr
// E-MAIL  :  P. Jolivet <pierre.jolivet@enseeiht.fr>

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep: parmetis metis mpi
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include <ff++.hpp>
#include <cmath>
#include <parmetis.h>
#ifndef IDX_T
#define IDX_T MPI_INT
#endif

template<class Type, class Mesh>
class ParMETIS_Op : public E_F0mps {
public:
    Expression part;
    Expression pTh;
    Expression lparts;
    static const int n_name_param = 2;
    static basicAC_F0::name_and_type name_param[];
    Expression nargs[n_name_param];
    ParMETIS_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : part(param1), pTh(param2), lparts(param3) {
        args.SetNameParam(n_name_param, name_param, nargs);
    }
    AnyType operator()(Stack stack) const;
};
template<class Type, class Mesh>
basicAC_F0::name_and_type ParMETIS_Op<Type, Mesh>::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"worker", &typeid(long)}
};
template<class Type, class Mesh>
class ParMETIS : public OneOperator {
public:
    ParMETIS() : OneOperator(atype<long>(), atype<KN<Type>*>(), atype<const Mesh*>(), atype<long>()) { }
    E_F0* code(const basicAC_F0& args) const {
        return new ParMETIS_Op<Type, Mesh>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
    }
};
template<class Type, class Mesh>
AnyType ParMETIS_Op<Type, Mesh>::operator()(Stack stack) const {
    KN<Type>* ptKN = GetAny<KN<Type>*>((*part)(stack));
    idx_t nparts = GetAny<long>((*lparts)(stack));
    Type* pt = *ptKN;
    long n = ptKN->n;
    idx_t* ptInt = sizeof(idx_t) <= sizeof(Type) ? reinterpret_cast<idx_t*>(pt) : new idx_t[n];
    std::fill_n(ptInt, n, 0);
    MPI_Comm comm = nargs[0] ? *((MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack))) : MPI_COMM_WORLD;
    int worker = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : 0;
    MPI_Comm workComm = comm;
    if(worker == 0)
        MPI_Comm_size(comm, &worker);
    else {
        int size;
        MPI_Comm_size(comm, &size);
        worker = std::min(size, worker);
        MPI_Group worldGroup, workGroup;
        MPI_Comm_group(workComm, &worldGroup);
        int ranges[1][3];
        ranges[0][0] = 0;
        ranges[0][1] = worker - 1;
        ranges[0][2] = 1;
        MPI_Group_range_incl(worldGroup, 1, ranges, &workGroup);
        MPI_Comm_create(comm, workGroup, &workComm);
        MPI_Group_free(&worldGroup);
    }
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank < worker) {
        idx_t* vtxdist = new idx_t[worker + 1];
        vtxdist[0] = 0;
        for(int i = 1; i < worker; ++i)
            vtxdist[i] = vtxdist[i - 1] + n / worker;
        vtxdist[worker] = n;
        int loc = vtxdist[rank + 1] - vtxdist[rank];
        idx_t* xadg = new idx_t[loc + 1];
        const Mesh& Th(*GetAny<const Mesh*>((*pTh)(stack)));
        idx_t nve = Mesh::Rd::d + 1;
        std::vector<idx_t> adjncy;
        adjncy.reserve(loc * nve);
        xadg[0] = 0;
        for(int k = vtxdist[rank]; k < vtxdist[rank + 1]; ++k) {
            for(idx_t j = 0; j < nve; ++j) {
                int l = j;
                idx_t m = Th.ElementAdj(k, l);
                if(k != m && m > 0)
                    adjncy.push_back(m);
            }
            xadg[k + 1 - vtxdist[rank]] = adjncy.size();
        }
#if 0
        for(int i = 0; i < worker; ++i) {
            MPI_Barrier(workComm);
            if(i == rank) {
                for(int j = 0; j < worker + 1; ++j) {
                    std::cout << vtxdist[j] << " ";
                }
                std::cout << std::endl;
                for(int j = 0; j < loc + 1; ++j) {
                    std::cout << xadg[j] << " ";
                }
                std::cout << std::endl;
                for(int j = 0; j < adjncy.size(); ++j) {
                    std::cout << adjncy[j] << " ";
                }
                std::cout << std::endl;
            }
            MPI_Barrier(workComm);
        }
#endif
        idx_t wgtflag = 0;
        idx_t ncon = 1;
        idx_t edgecut;
        real_t* tpwgts = new real_t[nparts];
        for(int i = 0; i < nparts; ++i)
            tpwgts[i] = 1.0 / static_cast<real_t>(nparts);
        real_t ubvec = 1.05;
        idx_t* part = ptInt + vtxdist[rank];
        ParMETIS_V3_PartKway(vtxdist, xadg, adjncy.data(), NULL, NULL, &wgtflag, &wgtflag, &ncon, &nparts, tpwgts, &ubvec, &wgtflag, &edgecut, part, &workComm);
        delete [] tpwgts;
        delete [] xadg;
        delete [] vtxdist;
    }
    MPI_Allreduce(MPI_IN_PLACE, ptInt, n, IDX_T, MPI_SUM, comm);
    for(int i = n; i-- > 0; )
        pt[i] = ptInt[i];
    if(sizeof(idx_t) > sizeof(Type))
        delete [] ptInt;
    if(nargs[1] && workComm != MPI_COMM_NULL)
        MPI_Comm_free(&workComm);
    return 0L;
}

static void Load_Init() {
    Global.Add("parmetis", "(", new ParMETIS<long, Mesh>);
    Global.Add("parmetis", "(", new ParMETIS<double, Mesh>);
    Global.Add("parmetis", "(", new ParMETIS<long, Mesh3>);
    Global.Add("parmetis", "(", new ParMETIS<double, Mesh3>);
}

LOADFUNC(Load_Init)
