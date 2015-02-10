//ff-c++-LIBRARY-dep:   hpddm  umfpack amd mumps scalapack blas [mkl]   mpifc  fc mpi  pthread   [petsc]
//ff-c++-cpp-dep: 

#ifdef  WITH_petsc
#ifndef __APPLE__
#define WITH_PETSC
#endif
#endif

#ifdef WITH_mkl
#ifdef __APPLE__
//#include "mkl.h"
// INTEL_MKL_VERSION 
//#if INTEL_MKL_VERSION <= 110103
//#undef INTEL_MKL_VERSION
// MKL to old ..... 
#define HPDDM_MKL            0
#define HAVE_cblas_axpby
#pragma message("schwarz plugin without MKL because MKL is too old " )
#else
#define HPDDM_MKL            1
#endif
#endif
#ifndef WITH_PETSC 
#pragma message("schwarz plugin compile without PETSc")
#endif

#include <math.h>
#include <mpi.h>
#include <ff++.hpp>
#include "AFunction_ext.hpp"
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

#ifdef WITH_PETSC
#include "PETSc.hpp"
#endif

namespace Schwarz {
template<class Type, class K>
class initDDM_Op : public E_F0mps {
    public:
        Expression A;
        Expression Mat;
        Expression o;
        Expression R;
        static const int n_name_param = 3;
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
    {"scaling", &typeid(KN<typename HPDDM::Wrapper<K>::ul_type>*)},
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
    MatriceMorse<K>* mA = static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*Mat)(stack))->A));
    KN<long>* ptO = GetAny<KN<long>*>((*o)(stack));
    KN<KN<long>>* ptR = GetAny<KN<KN<long>>*>((*R)(stack));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    if(ptO && mA) {
        HPDDM::MatrixCSR<K>* dA = new HPDDM::MatrixCSR<K>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->HPDDM::template Subdomain<K>::initialize(dA, STL<long>(*ptO), *ptR, comm);
    }
    FEbaseArrayKn<K>* deflation = nargs[2] ? GetAny<FEbaseArrayKn<K>*>((*nargs[2])(stack)) : 0;
    K** const& v = ptA->getVectors();
    if(deflation && deflation->N > 0 && !v) {
        K** ev = new K*[deflation->N];
        *ev = new K[deflation->N * deflation->get(0)->n];
        for(int i = 0; i < deflation->N; ++i) {
            ev[i] = *ev + i * deflation->get(0)->n;
            std::copy(&(*deflation->get(i))[0], &(*deflation->get(i))[deflation->get(i)->n], ev[i]);
        }
        ptA->setVectors(ev);
        ptA->Type::super::initialize(deflation->N);
    }
    KN<typename HPDDM::Wrapper<K>::ul_type>* ptD = nargs[1] ? GetAny<KN<typename HPDDM::Wrapper<K>::ul_type>*>((*nargs[1])(stack)) : 0;
    if(ptD)
        ptA->initialize(*ptD);
    else
        std::cerr << "Something is really wrong here !" << std::endl;
    return ptA;
}

template<class Type, class K>
class attachCoarseOperator_Op : public E_F0mps {
    public:
        Expression comm;
        Expression A;
        static const int n_name_param = 7;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        attachCoarseOperator_Op(const basicAC_F0& args, Expression param1, Expression param2) : comm(param1), A(param2) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type attachCoarseOperator_Op<Type, K>::name_param[] = {
    {"A", &typeid(Matrice_Creuse<K>*)},
    {"B", &typeid(Matrice_Creuse<K>*)},
    {"pattern", &typeid(Matrice_Creuse<K>*)},
    {"threshold", &typeid(typename HPDDM::Wrapper<K>::ul_type)},
    {"parameters", &typeid(KN<long>*)},
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
    MatriceMorse<K>* mA = nargs[0] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[0])(stack))->A)) : 0;
    KN<long>* ptParm = nargs[4] ? GetAny<KN<long>*>((*nargs[4])(stack)) : 0;
    if(ptParm && ptParm->n != 5) {
        if(ptParm->n == 1) {
            ptParm->resize(5);
            (*ptParm)[HPDDM::P] = 1;
        }
        else if(ptParm->n == 2)
            ptParm->resize(5);
        else
            cout << "Input array must be of size 1, 2, or 5 !" << endl;
        if(ptParm->n == 5) {
            (*ptParm)[HPDDM::TOPOLOGY] = 0;
            (*ptParm)[HPDDM::DISTRIBUTION] = HPDDM::DMatrix::NON_DISTRIBUTED;
            (*ptParm)[HPDDM::STRATEGY] = 3;
        }
    }
    KN<double>* timing = nargs[5] ? GetAny<KN<double>*>((*nargs[5])(stack)) : 0;
    Pair<K>* pair = nargs[6] ? GetAny<Pair<K>*>((*nargs[6])(stack)) : 0;
    typename HPDDM::Wrapper<K>::ul_type threshold = nargs[3] ? GetAny<typename HPDDM::Wrapper<K>::ul_type>((*nargs[3])(stack)) : 0;
    std::vector<unsigned short> parm(5);
    if(ptParm) {
        parm[HPDDM::P]            = (*ptParm)[HPDDM::P];
        parm[HPDDM::TOPOLOGY]     = (*ptParm)[HPDDM::TOPOLOGY];
        parm[HPDDM::DISTRIBUTION] = (*ptParm)[HPDDM::DISTRIBUTION];
        parm[HPDDM::STRATEGY]     = (*ptParm)[HPDDM::STRATEGY];
    }
    else {
        parm[HPDDM::P]            = 1;
        parm[HPDDM::TOPOLOGY]     = 0;
        parm[HPDDM::DISTRIBUTION] = HPDDM::DMatrix::NON_DISTRIBUTED;
        parm[HPDDM::STRATEGY]     = 3;
    }
    unsigned short nu = ptParm ? (*ptParm)[HPDDM::NU] : 20;
    nu = std::max(nu, static_cast<unsigned short>(1));
    std::pair<MPI_Request, const K*>* ret = nullptr;
    if(mA) {
        HPDDM::MatrixCSR<K> dA(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        MatriceMorse<K>* mB = nargs[1] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[1])(stack))->A)) : nullptr;
        MatriceMorse<K>* mP = nargs[2] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[2])(stack))->A)) : nullptr;
        const HPDDM::MatrixCSR<K>* const dP = mP ? new HPDDM::MatrixCSR<K>(mP->n, mP->m, mP->nbcoef, mP->a, mP->lg, mP->cl, mP->symetrique) : nullptr;
        double t = MPI_Wtime();
        if(mB) {
            HPDDM::MatrixCSR<K> dB(mB->n, mB->m, mB->nbcoef, mB->a, mB->lg, mB->cl, mB->symetrique);
            ptA->template solveGEVP<EIGENSOLVER>(&dA, nu, threshold, &dB, dP);
        }
        else
            ptA->template solveGEVP<EIGENSOLVER>(&dA, nu, threshold, nullptr, dP);
        mA->nbcoef = dA._nnz;
        mA->a = dA._a;
        mA->lg = dA._ia;
        mA->cl = dA._ja;
        if(timing)
            (*timing)[3] = MPI_Wtime() - t;
        delete dP;
        parm[HPDDM::NU]           = nu;
        if(ptA->exclusion(comm)) {
            if(pair)
                pair->p = ptA->template buildTwo<1>(comm, parm);
            else
                ret = ptA->template buildTwo<1>(comm, parm);
        }
        else {
            if(pair)
                pair->p = ptA->template buildTwo<0>(comm, parm);
            else
                ret = ptA->template buildTwo<0>(comm, parm);
        }
    }
    else {
        if(!threshold) {
            parm[HPDDM::NU]       = nu;
            ret = ptA->template buildTwo<2>(comm, parm);
        }
    }
    delete ret;
    return 0L;
}

template<class Type, class K>
class solveDDM_Op : public E_F0mps {
    public:
        Expression A;
        Expression x;
        Expression rhs;
        static const int n_name_param = 9;
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
    {"dim", &typeid(long)},
    {"iter", &typeid(long)},
    {"timing", &typeid(KN<double>*)},
    {"solver", &typeid(long)},
    {"pipelined", &typeid(long)},
    {"excluded", &typeid(bool)},
    {"ret", &typeid(Pair<K>*)},
    {"O", &typeid(Matrice_Creuse<K>*)}
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
    typename HPDDM::Wrapper<K>::ul_type eps = nargs[0] ? GetAny<typename HPDDM::Wrapper<K>::ul_type>((*nargs[0])(stack)) : 1e-8;
    unsigned short dim = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : 50;
    unsigned short iter = nargs[2] ? GetAny<long>((*nargs[2])(stack)) : 50;
    KN<double>* timing = nargs[3] ? GetAny<KN<double>*>((*nargs[3])(stack)) : 0;
    Pair<K>* pair = nargs[7] ? GetAny<Pair<K>*>((*nargs[7])(stack)) : 0;
    if(pair)
        if(pair->p) {
            int flag;
            MPI_Test(&(pair->p->first), &flag, MPI_STATUS_IGNORE);
        }
    MatriceMorse<K>* mA = nargs[8] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[8])(stack))->A)) : 0;
    long solver = nargs[4] ? GetAny<long>((*nargs[4])(stack)) : 0;
    double timer = MPI_Wtime();
    if(mA && (solver == 6 || solver == 7)) {
        HPDDM::MatrixCSR<K> dA(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->callNumfact(&dA);
    }
    else
        ptA->callNumfact(nullptr);
    ptA->setType(solver == 1 || solver == 6 || solver == 11);
    if(timing) (*timing)[1] = MPI_Wtime() - timer;
    long pipelined = nargs[5] ? GetAny<long>((*nargs[5])(stack)) : 0;
    bool excluded = nargs[6] ? GetAny<bool>((*nargs[6])(stack)) : false;
    if(pair)
        if(pair->p) {
            if(timing)
                timer = MPI_Wtime();
            MPI_Wait(&(pair->p->first), MPI_STATUS_IGNORE);
            if(timing)
                (*timing)[4] = MPI_Wtime() - timer;
            delete [] pair->p->second;
            delete pair->p;
        }
    timer = MPI_Wtime();
    int rank;
    MPI_Comm_rank(ptA->getCommunicator(), &rank);
    if(!excluded) {
        if(solver == 1)
            HPDDM::IterativeMethod::CG(*ptA, (K*)*ptX, (K*)*ptRHS, iter, eps, MPI_COMM_WORLD, rank == 0 && !excluded ? 1 : 0);
        else
            switch(pipelined) {
#if (OMPI_MAJOR_VERSION > 1 || (OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION >= 7)) || MPICH_NUMVERSION >= 30000000
                case 1:  HPDDM::IterativeMethod::GMRES<HPDDM::PIPELINED>(*ptA, (K*)*ptX, (K*)*ptRHS, dim, iter, eps, MPI_COMM_WORLD, rank == 0 && !excluded ? 1 : 0); break;
#endif
#if defined(DPASTIX) || defined(DMKL_PARDISO)
                case 2:  HPDDM::IterativeMethod::GMRES<HPDDM::FUSED>(*ptA, (K*)*ptX, (K*)*ptRHS, dim, iter, eps, MPI_COMM_WORLD, rank == 0 && !excluded ? 1 : 0); break;
#endif
                default: HPDDM::IterativeMethod::GMRES<HPDDM::CLASSICAL>(*ptA, (K*)*ptX, (K*)*ptRHS, dim, iter, eps, MPI_COMM_WORLD, rank == 0 && !excluded ? 1 : 0); break;
            }
    }
    else {
        if(solver == 1)
            HPDDM::IterativeMethod::CG(*ptA, (K*)nullptr, (K*)nullptr, iter, eps, MPI_COMM_WORLD, rank == 0 && !excluded ? 1 : 0);
        else
            switch(pipelined) {
#if (OMPI_MAJOR_VERSION > 1 || (OMPI_MAJOR_VERSION == 1 && OMPI_MINOR_VERSION >= 7)) || MPICH_NUMVERSION >= 30000000
                case 1:  HPDDM::IterativeMethod::GMRES<HPDDM::PIPELINED, true>(*ptA, (K*)nullptr, (K*)nullptr, dim, iter, eps, MPI_COMM_WORLD, rank == 0 && !excluded ? 1 : 0); break;
#endif
#if defined(DPASTIX) || defined(DMKL_PARDISO)
                case 2:  HPDDM::IterativeMethod::GMRES<HPDDM::FUSED, true>(*ptA, (K*)nullptr, (K*)nullptr, dim, iter, eps, MPI_COMM_WORLD, rank == 0 && !excluded ? 1 : 0); break;
#endif
                default: HPDDM::IterativeMethod::GMRES<HPDDM::CLASSICAL, true>(*ptA, (K*)nullptr, (K*)nullptr, dim, iter, eps, MPI_COMM_WORLD, rank == 0 && !excluded ? 1 : 0); break;
            }
    }
    timer = MPI_Wtime() - timer;
    if(!excluded) {
        if(rank == 0)
            std::cout << scientific << " --- system solved (in " << timer << ")" << std::endl;
        typename HPDDM::Wrapper<K>::ul_type storage[2];
        ptA->computeError(*ptX, *ptRHS, storage);
        if(rank == 0)
            std::cout << scientific << " --- error = " << storage[1] << " / " << storage[0] << std::endl;
    }
    return 0L;
}

double distributedDot(KN<double>* const& A, KN<double>* const& in, KN<double>* const& out) {
    double* tmp = new double[in->n];
    HPDDM::Wrapper<double>::diagv(in->n, *A, *in, tmp);
    KN_<double> KN(tmp, in->n);
    double dot = (KN, *out);
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    delete [] tmp;
    return dot;
}

template<class T, class U>
class GMV {
    public:
        const T t;
        const U u;
        GMV(T v, U w) : t(v), u(w) {}
        void prod(U x) const;
};
template<class Type, class U>
class GMV<Type*, U> {
    public:
        const Type* t;
        const U u;
        GMV(Type* v, U w) : t(v), u(w) {}
        void prod(U x) const { t->GMV(*(this->u), *x); };
};
template<class R, class A, class B> R Build(A a, B b) {
    return R(a, b);
}
template<class Type, class K>
KN<K>* GlobalMV(KN<K>* Ax, GMV<Type*, KN<K>*> A) {
    A.prod(Ax);
    return Ax;
}

template<template<class, char> class Type, class K, char S>
void add() {
    Dcl_Type<Type<K, S>*>(Initialize<Type<K, S>>, Delete<Type<K, S>>);

    TheOperators->Add("<-", new initDDM<Type<K, S>, K>);
    Global.Add("attachCoarseOperator", "(", new attachCoarseOperator<Type<K, S>, K>);
    Global.Add("DDM", "(", new solveDDM<Type<K, S>, K>);
    Dcl_Type<GMV<Type<K, S>*, KN<K>*>>();
    TheOperators->Add("*", new OneOperator2<GMV<Type<K, S>*, KN<K>*>, Type<K, S>*, KN<K>*>(Build));
    TheOperators->Add("=", new OneOperator2<KN<K>*, KN<K>*, GMV<Type<K, S>*, KN<K>*>>(GlobalMV));
    if(std::is_same<K, double>::value)
        Global.Add("dscalprod", "(", new OneOperator3_<double, KN<double>*, KN<double>*, KN<double>*>(distributedDot));
}
}

class Init {
    public:
        Init();
};

LOADINIT(Init);

Init::Init() {
    const char ds = 'S';
    const char zs = 'G';
    Schwarz::add<HpSchwarz, double, ds>();
    zzzfff->Add("dschwarz", atype<HpSchwarz<double, ds>*>());
    Schwarz::add<HpSchwarz, std::complex<double>, zs>();
    zzzfff->Add("zschwarz", atype<HpSchwarz<std::complex<double>, zs>*>());
#ifdef WITH_PETSC
    int argc = pkarg->n;
    char** argv = new char*[argc];
    for(int i = 0; i < argc; ++i)
        argv[i] = const_cast<char*>((*(*pkarg)[i].getap())->c_str());
    PetscInitialize(&argc, &argv, 0, "");
    delete [] argv;
    ff_atend(finalizePETSc);
    Dcl_Type<DistributedCSR*>(Initialize<DistributedCSR>, Delete<DistributedCSR>);
    Dcl_Type<GMV<DistributedCSR*, KN<double>*>>();
    zzzfff->Add("dmatrix", atype<DistributedCSR*>());
    TheOperators->Add("<-", new OneOperator1_<long, DistributedCSR*>(initEmptyCSR));
    TheOperators->Add("<-", new initCSR<double>);
    TheOperators->Add("*", new OneOperator2<GMV<DistributedCSR*, KN<double>*>, DistributedCSR*, KN<double>*>(Build));
    TheOperators->Add("=", new OneOperator2<KN<double>*, KN<double>*, GMV<DistributedCSR*, KN<double>*>>(GlobalMV));
    Global.Add("set", "(", new setOptions<double>());
    Global.Add("distributedNumbering", "(", new OneOperator2_<long, DistributedCSR*, KN<double>*>(distributedNumbering));
    Global.Add("renumber", "(", new OneOperator2_<long, DistributedCSR*, FEbaseArrayKn<double>*>(renumberCSR));
    Global.Add("solvePETSc", "(", new OneOperator2_<long, DistributedCSR*, KN<double>*>(solvePETSc));
    Dcl_Type<DistributedCSR_inv>();
    TheOperators->Add("^", new OneBinaryOperatorPETSc());
    Dcl_Type<Inv<DistributedCSR_inv, KN<double>*>>();
    TheOperators->Add("*", new OneOperator2<Inv<DistributedCSR_inv, KN<double>*>, DistributedCSR_inv, KN<double>*>(Build));
    TheOperators->Add("=", new OneOperator2<KN<double>*, KN<double>*, Inv<DistributedCSR_inv, KN<double>*>>(InvPETSc));
#endif
    Dcl_Type<Pair<double>*>(InitP<Pair<double>>, Destroy<Pair<double>>);
    zzzfff->Add("dpair", atype<Pair<double>*>());
    Dcl_Type<Pair<std::complex<double>>*>(InitP<Pair<std::complex<double>>>, Destroy<Pair<std::complex<double>>>);
    zzzfff->Add("zpair", atype<Pair<std::complex<double>>*>());
}
