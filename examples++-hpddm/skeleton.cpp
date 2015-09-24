#ifndef _ALL_IN_ONE_
#include "ff++.hpp"
#include <vector>
#include <cmath>
#endif

class Skeleton_Op : public E_F0mps {
    public:
        Expression interface;
        Expression index;
        Expression restriction;
        Expression outInterface;
        static const int n_name_param = 3;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        Skeleton_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3, Expression param4) : interface(param1), index(param2), restriction(param3), outInterface(param4) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }
        AnyType operator()(Stack stack) const;
};

basicAC_F0::name_and_type Skeleton_Op::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"interface", &typeid(KN<long>*)},
    {"redundancy", &typeid(bool)}
};

class Skeleton : public OneOperator {
    public:
        Skeleton() : OneOperator(atype<long>(), atype<KN<double>*>(), atype<KN<long>*>(), atype<KN<Matrice_Creuse<double> >*>(), atype<KN<KN<long> >*>()) {}

        E_F0* code(const basicAC_F0& args) const
        {
            return new Skeleton_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
        }
};

AnyType Skeleton_Op::operator()(Stack stack) const {
    KN<double>* in = GetAny<KN<double>*>((*interface)(stack));
    KN<KN<long> >* out = GetAny<KN<KN<long> >*>((*outInterface)(stack));
    KN<long>* arrayNeighbor = GetAny<KN<long>*>((*index)(stack));
    KN<Matrice_Creuse<double> >* interpolation = GetAny<KN<Matrice_Creuse<double> >*>((*restriction)(stack));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    unsigned short n = arrayNeighbor->n;
    MPI_Request* rq = new MPI_Request[2 * n];
    std::vector<unsigned char*> send(n);
    std::vector<unsigned char*> recv(n);
    unsigned short neighborAfter = 0;
    if(out->n != n)
        out->resize(n);
    for(unsigned short i = 0; i < n; ++i) {
        MatriceMorse<double>* pt = static_cast<MatriceMorse<double>*>(&(*interpolation->operator[](i).A));
        send[i] = new unsigned char[pt->n];
        recv[i] = new unsigned char[pt->n];
        unsigned int dest = arrayNeighbor->operator[](i);
        if(dest < mpirank) {
            unsigned int col = 0;
            for(unsigned int j = 0; j < pt->n; ++j) {
                if(pt->lg[j + 1] != pt->lg[j]) {
                    if(std::abs(in->operator[](pt->cl[col++]) - 1.0) < 0.1)
                        send[i][j] = '1';
                    else
                        send[i][j] = '0';
                }
                else
                    send[i][j] = '0';
            }
            MPI_Isend(send[i], pt->n, MPI_UNSIGNED_CHAR, dest, 0, *comm, rq + i);
            ++neighborAfter;
        }
        else
            MPI_Irecv(recv[i], pt->n, MPI_UNSIGNED_CHAR, dest, 0, *comm, rq + i);
    }
    for(unsigned short i = 0; i < neighborAfter; ++i) {
        MatriceMorse<double>* pt = static_cast<MatriceMorse<double>*>(&(*interpolation->operator[](i).A));
        // cout << mpirank << " receives from " << arrayNeighbor->operator[](i) << ", " << pt->n << "." << endl;
        MPI_Irecv(recv[i], pt->n, MPI_UNSIGNED_CHAR, arrayNeighbor->operator[](i), 0, *comm, rq + n + i);
    }
    for(unsigned short i = neighborAfter; i < n; ++i) {
        int index;
        MPI_Waitany(n - neighborAfter, rq + neighborAfter, &index, MPI_STATUS_IGNORE);
        unsigned short dest = neighborAfter + index;
        MatriceMorse<double>* pt = static_cast<MatriceMorse<double>*>(&(*interpolation->operator[](dest).A));
        KN<long>& resOut = out->operator[](dest);
        resOut.resize(pt->nbcoef);
        unsigned int nnz = 0;
        unsigned int col = 0;
        for(unsigned int j = 0; j < pt->n; ++j) {
            if(pt->lg[j + 1] != pt->lg[j]) {
                if(std::abs(in->operator[](pt->cl[col]) - 1.0) < 0.1 && recv[dest][j] == '1') {
                    send[dest][j] = '1';
                    resOut[(int)nnz++] = pt->cl[col++];
                }
                else {
                    ++col;
                    send[dest][j] = '0';
                }
            }
            else
                send[dest][j] = '0';
        }
        // cout << mpirank << " sends to " << arrayNeighbor->operator[](dest) << ", " << pt->n << "." << endl;
        MPI_Isend(send[dest], pt->n, MPI_UNSIGNED_CHAR, arrayNeighbor->operator[](dest), 0, *comm, rq + n + dest);
        resOut.resize(nnz);
    }
    for(unsigned short i = 0; i < neighborAfter; ++i) {
        int index;
        MPI_Waitany(neighborAfter, rq + n, &index, MPI_STATUS_IGNORE);
        KN<long>& resOut = out->operator[](index);
        MatriceMorse<double>* pt = static_cast<MatriceMorse<double>*>(&(*interpolation->operator[](index).A));
        resOut.resize(pt->nbcoef);
        unsigned int nnz = 0;
        unsigned int col = 0;
        for(unsigned int j = 0; j < pt->n; ++j) {
            if(recv[index][j] == '1') {
                if(pt->lg[j + 1] != pt->lg[j])
                    resOut[(int)nnz++] = pt->cl[col++];
            }
            else if(pt->lg[j + 1] != pt->lg[j])
                ++col;
        }
        resOut.resize(nnz);
    }
    MPI_Waitall(neighborAfter, rq, MPI_STATUSES_IGNORE);
    MPI_Waitall(n - neighborAfter, rq + n + neighborAfter, MPI_STATUSES_IGNORE);
    for(unsigned short i = 0; i < n; ++i) {
        delete [] recv[i];
        delete [] send[i];
        // cout << mpirank << " <=> " << arrayNeighbor->operator[](i) << " : " << out->operator[](i).n << endl;
    }
    delete [] rq;
    KN<long>* interfaceNb = nargs[1] ? GetAny<KN<long>* >((*nargs[1])(stack)) : (KN<long>*) 0;
    if(interfaceNb) {
        std::vector<unsigned int> vec;
        vec.reserve(in->n);
        for(int i = 0; i < in->n; ++i) {
            if(in->operator[](i) != 0.0)
                vec.emplace_back(i);
        }
        std::sort(vec.begin(), vec.end());
        if(interfaceNb->n != vec.size())
            interfaceNb->resize(vec.size());
        for(  signed int i = 0; i < vec.size(); ++i)
            interfaceNb->operator[](i) = vec[i];
        for(unsigned short i = 0; i < n; ++i) {
            KN<long>& res = out->operator[](i);
            for(  signed int j = 0; j < res.n; ++j) {
                std::vector<unsigned int>::const_iterator idx = std::lower_bound(vec.cbegin(), vec.cend(), (unsigned int)res[j]);
                if(idx == vec.cend() || res[j] < *idx) {
                    std::cout << "Problem !" << std::endl;
                    res[j] = -1;
                }
                else
                    res[j] = std::distance(vec.cbegin(), idx);
            }
        }
        bool redundancy = nargs[2] ? GetAny<bool>((*nargs[2])(stack)) : 1;
        if(!redundancy) {
            std::vector<std::pair<unsigned short, unsigned int> >* array = new std::vector<std::pair<unsigned short, unsigned int> >[interfaceNb->n];
            for(unsigned short i = 0; i < n; ++i) {
                KN<long>& res = out->operator[](i);
                for(  signed int j = 0; j < res.n; ++j)
                    array[res[j]].push_back(std::make_pair(i, res[j]));
            }
            for(unsigned int i = 0; i < interfaceNb->n; ++i) {
                if(array[i].size() > 1) {
                    if(mpirank > arrayNeighbor->operator[](array[i].back().first))
                        array[i].erase(array[i].begin());
                    else if(mpirank < arrayNeighbor->operator[](array[i].front().first))
                        array[i].pop_back();
                }
                else if (array[i].size() < 1)
                    std::cout << "Problem !" << std::endl;
            }
            std::vector<long>* copy = new std::vector<long>[n];
            for(unsigned short i = 0; i < n; ++i)
                copy[i].reserve(out->operator[](i).n);
            for(unsigned int i = 0; i < interfaceNb->n; ++i) {
                for(std::vector<std::pair<unsigned short, unsigned int> >::const_iterator it = array[i].cbegin(); it != array[i].cend(); ++it) {
                    copy[it->first].push_back(it->second);
                }
            }
            for(unsigned short i = 0; i < n; ++i) {
                unsigned int sizeVec = copy[i].size();
                if(sizeVec != out->operator[](i).n) {
                    out->operator[](i).resize(sizeVec);
                    long* pt = (static_cast<KN_<long> >(out->operator[](i)));
                    std::reverse_copy(copy[i].begin(), copy[i].end(), pt);
                }
            }
            delete [] copy;
            delete [] array;
        }
    }
    return 0L;
}

#ifndef _ALL_IN_ONE_
static void Init_Skeleton() {
    Global.Add("buildSkeleton", "(", new Skeleton);
}

LOADFUNC(Init_Skeleton)
#endif
