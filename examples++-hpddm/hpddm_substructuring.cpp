//ff-c++-LIBRARY-dep: cxx11 hpddm [umfpack|mumps parmetis ptscotch scotch scalapack] [mkl|blas] mpi pthread mpifc fc
//ff-c++-cpp-dep:

#define HPDDM_SCHWARZ 0
#define HPDDM_FETI    1
#define HPDDM_BDD     1

#include "common.hpp"

namespace Substructuring {
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


template<class Type, class K>
class initDDM_Op : public E_F0mps {
    public:
        Expression A;
        Expression Mat;
        Expression o;
        Expression R;
        static const int n_name_param = 2;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        initDDM_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3, Expression param4) : A(param1), Mat(param2), o(param3), R(param4) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type initDDM_Op<Type, K>::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"deflation", &typeid(FEbaseArrayKn<K>*)},
};
template<class Type, class K>
class initDDM : public OneOperator {
    public:
        initDDM() : OneOperator(atype<Type*>(), atype<Type*>(), atype<Matrice_Creuse<K>*>(), atype<KN<long>*>(), atype<KN<KN<long>>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new initDDM_Op<Type, K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
        }
};
template<class Type, class K>
AnyType initDDM_Op<Type, K>::operator()(Stack stack) const {
    Type* ptA = GetAny<Type*>((*A)(stack));
    Matrice_Creuse<K>* pA = GetAny<Matrice_Creuse<K>*>((*Mat)(stack));
    MatriceMorse<K>* mA = pA->A ? static_cast<MatriceMorse<K>*>(&(*pA->A)) : nullptr;
    KN<long>* ptO = GetAny<KN<long>*>((*o)(stack));
    KN<KN<long>>* ptR = GetAny<KN<KN<long>>*>((*R)(stack));
    if(ptO)
        ptA->HPDDM::template Subdomain<K>::initialize(mA ? new HPDDM::MatrixCSR<K>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique) : 0, STL<long>(*ptO), *ptR, nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0);
    FEbaseArrayKn<K>* deflation = nargs[1] ? GetAny<FEbaseArrayKn<K>*>((*nargs[1])(stack)) : 0;
    K** const& v = ptA->getVectors();
    if(deflation && deflation->N > 0 && !v) {
        K** ev = new K*[deflation->N];
        *ev = new K[deflation->N * deflation->get(0)->n];
        for(int i = 0; i < deflation->N; ++i) {
            ev[i] = *ev + i * deflation->get(0)->n;
            std::copy(&(*deflation->get(i))[0], &(*deflation->get(i))[deflation->get(i)->n], ev[i]);
        }
        ptA->setVectors(ev);
        ptA->Type::super::super::initialize(deflation->N);
    }
    return ptA;
}

template<class Type, class K>
class attachCoarseOperator_Op : public E_F0mps {
    public:
        Expression comm;
        Expression A;
        static const int n_name_param = 4;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        attachCoarseOperator_Op(const basicAC_F0& args, Expression param1, Expression param2) : comm(param1), A(param2) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type attachCoarseOperator_Op<Type, K>::name_param[] = {
    {"R", &typeid(FEbaseArrayKn<K>*)},
    {"threshold", &typeid(HPDDM::underlying_type<K>)},
    {"timing", &typeid(KN<double>*)},
    {"ret", &typeid(Pair<K>*)}
};
template<class Type, class K>
class attachCoarseOperator : public OneOperator {
    public:
        attachCoarseOperator() : OneOperator(atype<long>(), atype<pcommworld>(), atype<Type*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new attachCoarseOperator_Op<Type, K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
        }
};
template<class Type, class K>
AnyType attachCoarseOperator_Op<Type, K>::operator()(Stack stack) const {
    pcommworld ptComm = GetAny<pcommworld>((*comm)(stack));
    MPI_Comm comm = *(MPI_Comm*)ptComm;
    Type* ptA = GetAny<Type*>((*A)(stack));
    FEbaseArrayKn<K>* R = nargs[0] ? GetAny<FEbaseArrayKn<K>*>((*nargs[0])(stack)) : 0;
    Pair<K>* pair = nargs[3] ? GetAny<Pair<K>*>((*nargs[3])(stack)) : 0;
    unsigned short nu = R ? static_cast<unsigned short>(R->N) : 0;
    HPDDM::Option& opt = *HPDDM::Option::get();
    HPDDM::underlying_type<K> threshold = opt.val("geneo_threshold", 0.0);
    KN<double>* timing = nargs[2] ? GetAny<KN<double>*>((*nargs[2])(stack)) : 0;
    std::pair<MPI_Request, const K*>* ret = nullptr;
    bool adaptive = opt.set("geneo_nu") || threshold > 0.0;
    if(!adaptive)
        ptA->setDeficiency(nu);
    double t = MPI_Wtime();
    if(R) {
        if(adaptive)
            ptA->computeSchurComplement();
        ptA->callNumfactPreconditioner();
        if(timing)
            (*timing)[3] = MPI_Wtime() - t;
        if(adaptive) {
            if(opt.set("geneo_nu"))
                nu = opt["geneo_nu"];
#if defined(MUMPSSUB) || defined(PASTIXSUB) || defined(MKL_PARDISOSUB)
            t = MPI_Wtime();
            ptA->solveGEVP();
            if(timing)
                (*timing)[5] = MPI_Wtime() - t;
#else
            cout << "Please change your solver" << endl;
#endif
        }
        K** const ev = ptA->getVectors();
        if(!R && !ev)
            cout << "Problem !" << endl;
        R->resize(0);
        MPI_Barrier(MPI_COMM_WORLD);
        if(timing)
            t = MPI_Wtime();
        if(ptA->exclusion(comm)) {
            if(pair)
                pair->p = ptA->template buildTwo<1>(comm);
            else
                ret = ptA->template buildTwo<1>(comm);
        }
        else {
            if(pair)
                pair->p = ptA->template buildTwo<0>(comm);
            else
                ret = ptA->template buildTwo<0>(comm);
        }
        if(timing)
            (*timing)[4] = MPI_Wtime() - t;
        if(pair)
            if(pair->p) {
                int flag;
                MPI_Test(&(pair->p->first), &flag, MPI_STATUS_IGNORE);
            }
        t = MPI_Wtime();
        ptA->callNumfact();
    }
    else {
        MPI_Barrier(MPI_COMM_WORLD);
        ret = ptA->template buildTwo<2>(comm);
    }
    if(timing)
        (*timing)[2] = MPI_Wtime() - t;
    if(ret)
        delete ret;
    if(pair)
        if(pair->p) {
            if(timing)
                t = MPI_Wtime();
            MPI_Wait(&(pair->p->first), MPI_STATUS_IGNORE);
            if(timing)
                (*timing)[timing->n - 1] = MPI_Wtime() - t;
            delete [] pair->p->second;
            pair->destroy();
            pair = nullptr;
        }
    return 0L;
}

template<class Type, class K>
class solveDDM_Op : public E_F0mps {
    public:
        Expression A;
        Expression x;
        Expression rhs;
        static const int n_name_param = 5;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        solveDDM_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : A(param1), x(param2), rhs(param3) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type solveDDM_Op<Type, K>::name_param[] = {
    {"eps", &typeid(HPDDM::underlying_type<K>)},
    {"iter", &typeid(long)},
    {"timing", &typeid(KN<double>*)},
    {"excluded", &typeid(bool)},
    {"ret", &typeid(Pair<K>*)}
};
template<class Type, class K>
class solveDDM : public OneOperator {
    public:
        solveDDM() : OneOperator(atype<long>(), atype<Type*>(), atype<KN<K>*>(), atype<KN<K>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new solveDDM_Op<Type, K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
        }
};
template<class Type, class K>
AnyType solveDDM_Op<Type, K>::operator()(Stack stack) const {
    KN<K>* ptX = GetAny<KN<K>*>((*x)(stack));
    KN<K>* ptRHS = GetAny<KN<K>*>((*rhs)(stack));
    Type* ptA = GetAny<Type*>((*A)(stack));
    if(ptX->n != ptRHS->n)
        return 0L;
    HPDDM::Option& opt = *HPDDM::Option::get();
    HPDDM::underlying_type<K> eps = nargs[0] ? GetAny<HPDDM::underlying_type<K>>((*nargs[0])(stack)) : -1.0;
    if(nargs[0])
        std::cerr << "Please do not use the legacy option \"-eps\", set instead \"-hpddm_tol\", cf. \"-hpddm_help\"" << std::endl;
    if(std::abs(eps + 1.0) > 1.0e-6)
        opt["tol"] = eps;
    int iter = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : -1;
    if(iter != -1)
        opt["max_it"] = iter;
    KN<double>* timing = nargs[2] ? GetAny<KN<double>*>((*nargs[2])(stack)) : 0;
    bool excluded = nargs[3] && GetAny<bool>((*nargs[3])(stack));
    if(excluded)
        opt["master_exclude"];
    double timer = MPI_Wtime();
    if(mpisize == 1) {
        ptA->computeSchurComplement();
        ptA->callNumfactPreconditioner();
        if(timing)
            (*timing)[3] = MPI_Wtime() - timer;
        timer = MPI_Wtime();
        ptA->callNumfact();
        if(timing)
            (*timing)[2] = MPI_Wtime() - timer;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(!excluded && timing && mpisize > 1)
        (*timing)[timing->n - 1] += MPI_Wtime() - timer;
    timer = MPI_Wtime();
    int rank;
    MPI_Comm_rank(ptA->getCommunicator(), &rank);
    if(rank != mpirank || rank != 0)
        opt.remove("verbosity");
    timer = MPI_Wtime();
    if(!excluded)
        HPDDM::IterativeMethod::solve(*ptA, (K*)*ptRHS, (K*)*ptX, 1, MPI_COMM_WORLD);
    else
        HPDDM::IterativeMethod::solve<true>(*ptA, (K*)nullptr, (K*)nullptr, 1, MPI_COMM_WORLD);
    timer = MPI_Wtime() - timer;
    if(!excluded && verbosity > 0 && rank == 0)
        std::cout << scientific << " --- system solved (in " << timer << ")" << std::endl;
    return 0L;
}

template<class Type, class K>
class renumber_Op : public E_F0mps {
    public:
        Expression A;
        Expression Mat;
        Expression interface;
        static const int n_name_param = 4;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        renumber_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : A(param1), Mat(param2), interface(param3) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type renumber_Op<Type, K>::name_param[] = {
    {"R", &typeid(FEbaseArrayKn<K>*)},
    {"effort", &typeid(KN<K>*)},
    {"rho", &typeid(KN<K>*)},
    {"timing", &typeid(KN<double>*)}
};
template<class Type, class K>
class renumber : public OneOperator {
    public:
        renumber() : OneOperator(atype<long>(), atype<Type*>(), atype<Matrice_Creuse<K>*>(), atype<KN<long>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new renumber_Op<Type, K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
        }
};
template<class Type, class K>
AnyType renumber_Op<Type, K>::operator()(Stack stack) const {
    Type* ptA = GetAny<Type*>((*A)(stack));
    KN<long>* ptInterface = GetAny<KN<long>*>((*interface)(stack));
    FEbaseArrayKn<K>* deflation = nargs[0] ? GetAny<FEbaseArrayKn<K>*>((*nargs[0])(stack)) : 0;
    KN<K>* ptEffort = nargs[1] ? GetAny<KN<K>*>((*nargs[1])(stack)) : 0;
    KN<K>* rho = nargs[2] ? GetAny<KN<K>*>((*nargs[2])(stack)) : 0;
    KN<double>* timing = nargs[3] ? GetAny<KN<double>*>((*nargs[3])(stack)) : 0;
    double t = MPI_Wtime();
    K** ev;
    if(deflation && deflation->N > 0) {
        ev = new K*[deflation->N];
        *ev = new K[deflation->N * deflation->get(0)->n];
        for(int i = 0; i < deflation->N; ++i) {
            ev[i] = *ev + i * deflation->get(0)->n;
            std::copy(static_cast<K*>(*(deflation->get(i))), static_cast<K*>(*(deflation->get(i))) + deflation->get(i)->n, ev[i]);
        }
        ptA->setVectors(ev);
        ptA->Type::super::super::initialize(deflation->N);
    }

    ptA->renumber(STL<long>(*ptInterface), ptEffort ? static_cast<K*>(*ptEffort) : nullptr);
    MatriceMorse<K>* mA = static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*Mat)(stack))->A));
    if(mA) {
        const HPDDM::MatrixCSR<K>* dA = ptA->getMatrix();
        mA->lg = dA->_ia;
    }

    HPDDM::Option& opt = *HPDDM::Option::get();
    char scaling = opt.val<char>("substructuring_scaling", 0);
    if(scaling == 2 && rho) {
        ptA->renumber(STL<long>(*ptInterface), *rho);
        ptA->buildScaling(scaling, *rho);
    }
    else
        ptA->buildScaling(scaling);

    if(timing) {
        (*timing)[1] = MPI_Wtime() - t;
        t = MPI_Wtime();
    }

    if(deflation && deflation->N > 0)
        for(int i = 0; i < deflation->N; ++i)
            std::copy(ev[i], ev[i] + deflation->get(i)->n, static_cast<K*>(*(deflation->get(i))));
    return 0L;
}

template<class Type>
long nbMult(Type* const& A) {
    return A->getMult();
}
template<class Type>
double nbDof(Type* const& A) {
    return static_cast<double>(A->getAllDof());
}
template<class Type, class K>
long originalNumbering(Type* const& A, KN<K>* const& in, KN<long>* const& interface) {
    A->originalNumbering(STL<long>(*interface), *in);
    return 0;
}

template<class T, class U, class K>
class InvSubstructuring {
    public:
        const T t;
        const U u;
        InvSubstructuring(T v, U w) : t(v), u(w) {}
        void solve(U out) const {
            if(out->n != u->n)
                return;
            if(mpisize == 1) {
                (*t).computeSchurComplement();
                (*t).callNumfactPreconditioner();
                (*t).callNumfact();
            }
            HPDDM::Option& opt = *HPDDM::Option::get();
            if(mpirank != 0)
                opt.remove("verbosity");
            HPDDM::IterativeMethod::solve(*t, (K*)*u, (K*)*out, 1, MPI_COMM_WORLD);
        }
        static U init(U Ax, InvSubstructuring<T, U, K> A) {
            A.solve(Ax);
            return Ax;
        }
};

template<template<class, char> class Type, class K, char S>
void add() {
    Dcl_Type<Type<K, S>*>(Initialize<Type<K, S>>, Delete<Type<K, S>>);

    TheOperators->Add("<-", new initDDM<Type<K, S>, K>);
    Global.Add("attachCoarseOperator", "(", new attachCoarseOperator<Type<K, S>, K>);
    Global.Add("DDM", "(", new solveDDM<Type<K, S>, K>);
    Global.Add("renumber", "(", new renumber<Type<K, S>, K>);
    Global.Add("nbDof", "(", new OneOperator1_<double, Type<K, S>*>(nbDof));
    Global.Add("nbMult", "(", new OneOperator1_<long, Type<K, S>*>(nbMult));
    Global.Add("originalNumbering", "(", new OneOperator3_<long, Type<K, S>*, KN<K>*, KN<long>*>(originalNumbering));
    addInv<Type<K, S>, InvSubstructuring, KN<K>, K>();
    Global.Add("statistics", "(", new OneOperator1_<bool, Type<K, S>*>(statistics<Type<K, S>>));
    Global.Add("exchange", "(", new exchangeInOut<Type<K, S>, K>);
}
}

template<class K, char S>
using HpFetiPrec = HpFeti<HPDDM::FetiPrcndtnr::DIRICHLET, K, S>;
static void Init_Substructuring() {
    Init_Common();
    Global.Add("buildSkeleton", "(", new Substructuring::Skeleton);
#if defined(DSUITESPARSE) || defined(DHYPRE)
    const char ds = 'G';
#else
    const char ds = 'S';
#endif
    const char zs = 'G';
    Substructuring::add<HpBdd, double, ds>();
    zzzfff->Add("dbdd", atype<HpBdd<double, ds>*>());
#ifndef DHYPRE
    // Substructuring::add<HpBdd, float, ds>();
    // zzzfff->Add("sbdd", atype<HpBdd<float, ds>*>());
    Substructuring::add<HpBdd, std::complex<double>, zs>();
    zzzfff->Add("zbdd", atype<HpBdd<std::complex<double>, zs>*>());
    // Substructuring::add<HpBdd, std::complex<float>, zs>();
    // zzzfff->Add("cbdd", atype<HpBdd<std::complex<float>, zs>*>());
#endif
    Substructuring::add<HpFetiPrec, double, ds>();
    zzzfff->Add("dfeti", atype<HpFetiPrec<double, ds>*>());
#ifndef DHYPRE
    // Substructuring::add<HpFetiPrec, float, ds>();
    // zzzfff->Add("sfeti", atype<HpFetiPrec<float, ds>*>());
    Substructuring::add<HpFetiPrec, std::complex<double>, zs>();
    zzzfff->Add("zfeti", atype<HpFetiPrec<std::complex<double>, zs>*>());
    // Substructuring::add<HpFetiPrec, std::complex<float>, zs>();
    // zzzfff->Add("cfeti", atype<HpFetiPrec<std::complex<float>, zs>*>());
#endif
    aType t;
    int r;
    if(!zzzfff->InMotClef("dpair", t, r)) {
        // Dcl_Type<Pair<float>*>(InitP<Pair<float>>, Destroy<Pair<float>>);
        // zzzfff->Add("spair", atype<Pair<double>*>());
        Dcl_Type<Pair<double>*>(InitP<Pair<double>>, Destroy<Pair<double>>);
        zzzfff->Add("dpair", atype<Pair<double>*>());
        // Dcl_Type<Pair<std::complex<float>>*>(InitP<Pair<std::complex<float>>>, Destroy<Pair<std::complex<float>>>);
        // zzzfff->Add("cpair", atype<Pair<std::complex<float>>*>());
        Dcl_Type<Pair<std::complex<double>>*>(InitP<Pair<std::complex<double>>>, Destroy<Pair<std::complex<double>>>);
        zzzfff->Add("zpair", atype<Pair<std::complex<double>>*>());
    }
}

LOADFUNC(Init_Substructuring)
