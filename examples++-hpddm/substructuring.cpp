#ifndef _ALL_IN_ONE_
#include <math.h>
#include <mpi.h>
#include <ff++.hpp>
#include "AFunction_ext.hpp"
#ifdef SCHWARZ
#undef SCHWARZ
#endif
#include <HPDDM.hpp>

template<class T>
class STL {
    T* const _it;
    const int _size;
    public:
        STL(const KN<T>& v) : _it(v), _size(v.size()) { };
        int size() const {
            return _size;
        }
        T* begin() const {
            return _it;
        }
        T* end() const {
            return _it + _size;
        }
};
template<class K>
class Pair {
    public:
        Pair() : p() { };
        std::pair<MPI_Request, const K*>* p;
        void init() {
        }
        void destroy() {
        }
};
#endif

#ifdef FETI
template<class K, char S>
using HpFetiPrec = HpFeti<HPDDM::FetiPrcndtnr::DIRICHLET, K, S>;
#endif

extern KN<String>* pkarg;

namespace Substructuring {
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
        initDDM() : OneOperator(atype<Type*>(), atype<Type*>(), atype<Matrice_Creuse<K>*>(), atype<KN<long>*>(), atype<KN<KN<long> >*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new initDDM_Op<Type, K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
        }
};
template<class Type, class K>
AnyType initDDM_Op<Type, K>::operator()(Stack stack) const {
    const char** argv = new const char*[pkarg->n];
    for(int i = 0; i < pkarg->n; ++i)
        argv[i] = (*((*pkarg)[i].getap()))->data();
    HPDDM::Option::get()->parse(pkarg->n, argv, mpirank == 0);
    delete [] argv;
    Type* ptA = GetAny<Type*>((*A)(stack));
    MatriceMorse<K>* mA = static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*Mat)(stack))->A));
    KN<long>* ptO = GetAny<KN<long>*>((*o)(stack));
    KN<KN<long> >* ptR = GetAny<KN<KN<long> >*>((*R)(stack));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    if(ptO && mA) {
        HPDDM::MatrixCSR<K>* dA = new HPDDM::MatrixCSR<K>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->HPDDM::template Subdomain<K>::initialize(dA, STL<long>(*ptO), *ptR, comm);
    }
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
    {"threshold", &typeid(typename HPDDM::Wrapper<K>::ul_type)},
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
    if(opt["geneo_nu"] != 0)
        nu = opt["geneo_nu"];
    typename HPDDM::Wrapper<K>::ul_type threshold = opt.val("geneo_threshold", 0.0);
    KN<double>* timing = nargs[2] ? GetAny<KN<double>*>((*nargs[2])(stack)) : 0;
    std::pair<MPI_Request, const K*>* ret = nullptr;
    double t = MPI_Wtime();
    if(R) {
        if(opt["geneo_nu"] > 0 || threshold > 0.0)
            ptA->computeSchurComplement();
        else if(R->N)
            nu = opt["geneo_nu"] = R->N;
        ptA->callNumfactPreconditioner();
        if(timing)
            (*timing)[3] = MPI_Wtime() - t;
        if(nu != R->N || threshold > 0.0) {
#if defined(MUMPSSUB) || defined(PASTIXSUB) || defined(MKL_PARDISOSUB)
            t = MPI_Wtime();
            ptA->solveGEVP(nu, threshold);
            if(timing)
                (*timing)[5] = MPI_Wtime() - t;
            opt["geneo_nu"] = nu;
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
            delete pair->p;
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
    {"eps", &typeid(typename HPDDM::Wrapper<K>::ul_type)},
    {"iter", &typeid(long)},
    {"timing", &typeid(KN<double>*)},
    {"excluded", &typeid(long)},
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
    HPDDM::Option& opt = *HPDDM::Option::get();
    typename HPDDM::Wrapper<K>::ul_type eps = nargs[0] ? GetAny<typename HPDDM::Wrapper<K>::ul_type>((*nargs[0])(stack)) : -1.0;
    if(std::abs(eps + 1.0) > 1.0e-6)
        opt["tol"] = eps;
    int iter = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : -1;
    if(iter != -1)
        opt["max_it"] = iter;
    KN<double>* timing = nargs[2] ? GetAny<KN<double>*>((*nargs[2])(stack)) : 0;
    long excluded = nargs[3] ? GetAny<long>((*nargs[3])(stack)) : 0;
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
    if(excluded == 2 && timing && mpisize > 1)
        (*timing)[timing->n - 1] += MPI_Wtime() - timer;
    timer = MPI_Wtime();
    int rank;
    MPI_Comm_rank(ptA->getCommunicator(), &rank);
    if(rank != 0 || excluded)
        opt.remove("verbosity");
    timer = MPI_Wtime();
    if(excluded == 1)
        HPDDM::IterativeMethod::PCG<1>(*ptA, (K*)*ptX, (K*)*ptRHS, MPI_COMM_WORLD);
    else
        HPDDM::IterativeMethod::PCG<0>(*ptA, (K*)*ptX, (K*)*ptRHS, MPI_COMM_WORLD);
    timer = MPI_Wtime() - timer;
    if(excluded != 1) {
        if(rank == 0)
            std::cout << scientific << " --- system solved (in " << timer << ")" << std::endl;
        typename HPDDM::Wrapper<K>::ul_type storage[2];
        ptA->computeError(*ptX, *ptRHS, storage);
        if(rank == 0)
            std::cout << scientific << " --- error = " << storage[1] << " / " << storage[0] << std::endl;
    }
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
    std::vector<unsigned int> it(static_cast<long*>(*ptInterface), static_cast<long*>(*ptInterface + ptInterface->n));
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

    ptA->renumber(it, ptEffort ? static_cast<K*>(*ptEffort) : nullptr);
    MatriceMorse<K>* mA = static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*Mat)(stack))->A));
    if(mA) {
        const HPDDM::MatrixCSR<K>* dA = ptA->getMatrix();
        mA->lg = dA->_ia;
    }

    HPDDM::Option& opt = *HPDDM::Option::get();
    unsigned short scaling = opt["substructuring_scaling"];
    if(scaling == 2 && rho) {
        ptA->renumber(it, *rho);
        ptA->buildScaling(scaling, *rho);
    }
    else {
        ptA->buildScaling(scaling);
        opt["substructuring_scaling"] = scaling;
    }

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
    std::vector<unsigned int> it(static_cast<long*>(*interface), static_cast<long*>(*interface + interface->n));
    A->originalNumbering(it, *in);
    return 0;
}

template<template<class, char> class Type, class K, char S>
void add() {
    Dcl_Type<Type<K, S>*>(Initialize<Type<K, S> >, Delete<Type<K, S> >);

    TheOperators->Add("<-", new initDDM<Type<K, S>, K>);
    Global.Add("attachCoarseOperator", "(", new attachCoarseOperator<Type<K, S>, K>);
    Global.Add("renumber", "(", new renumber<Type<K, S>, K>);
    Global.Add("DDM", "(", new solveDDM<Type<K, S>, K>);
    Global.Add("nbDof", "(", new OneOperator1_<double, Type<K, S>*>(nbDof));
    Global.Add("nbMult", "(", new OneOperator1_<long, Type<K, S>*>(nbMult));
    Global.Add("originalNumbering", "(", new OneOperator3_<long, Type<K, S>*, KN<K>*, KN<long>*>(originalNumbering));
}
}

#ifndef _ALL_IN_ONE_
static void Init_Substructuring() {
#include "init.hpp"
}

LOADFUNC(Init_Substructuring)
#endif
