//ff-c++-LIBRARY-dep: cxx11 hpddm petsc umfpack mpi
//ff-c++-cpp-dep: 

#include "common.hpp"
#include "petsc.h"

#if PETSC_VERSION_LT(3,7,0)
#define FFPetscOptionsGetInt(a,b,c,d) PetscOptionsGetInt(a,b,c,d)
#define FFPetscOptionsGetReal(a,b,c,d) PetscOptionsGetReal(a,b,c,d)
#define FFPetscOptionsInsert(a,b,c) PetscOptionsInsert(a,b,c)
#else
#define FFPetscOptionsGetInt(a,b,c,d) PetscOptionsGetInt(NULL,a,b,c,d)
#define FFPetscOptionsGetReal(a,b,c,d) PetscOptionsGetReal(NULL,a,b,c,d)
#define FFPetscOptionsInsert(a,b,c) PetscOptionsInsert(NULL,a,b,c)
#endif

#if PETSC_VERSION_LT(3,6,0)
#define MatCreateVecs MatGetVecs
#endif

namespace PETSc {
class DistributedCSR {
    public:
        HpSchwarz<PetscScalar>*     _A;
        Mat                     _petsc;
        Vec                         _x;
        KSP                       _ksp;
        unsigned int*             _num;
        int*                       _ia;
        int*                       _ja;
        PetscScalar*                _c;
        unsigned int            _first;
        unsigned int             _last;
        unsigned int           _global;
        bool                     _free;
        DistributedCSR() : _A(), _num(), _ia(), _ja(), _c(), _first(), _last() { };
        ~DistributedCSR() {
            if(_A) {
                VecDestroy(&_x);
                MatDestroy(&_petsc);
                KSPDestroy(&_ksp);
                delete [] _num;
                _num = nullptr;
                if(_free) {
                    delete []  _ia;
                    delete []  _ja;
                    delete []   _c;
                    _ia = _ja = nullptr;
                    _c = nullptr;
                }
                _A->clearBuffer();
                delete _A;
                _A = nullptr;
            }
        }
};

void finalizePETSc() {
    PETSC_COMM_WORLD = MPI_COMM_WORLD;
    PetscFinalize();
}

long initEmptyCSR(DistributedCSR* const&) {
    return 0;
}

template<class Type>
class initCSR_Op : public E_F0mps {
    public:
        Expression A;
        Expression K;
        Expression O;
        Expression R;
        Expression D;
        static const int n_name_param = 3;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        initCSR_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3, Expression param4, Expression param5) : A(param1), K(param2), O(param3), R(param4), D(param5) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};

template<class Type>
basicAC_F0::name_and_type initCSR_Op<Type>::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"bs", &typeid(long)},
    {"rhs", &typeid(KN<PetscScalar>*)}
};

template<class Type>
class initCSR : public OneOperator {
    public:
        initCSR() : OneOperator(atype<DistributedCSR*>(), atype<DistributedCSR*>(), atype<Matrice_Creuse<PetscScalar>*>(), atype<KN<long>*>(), atype<KN<KN<long>>*>(), atype<KN<double>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new initCSR_Op<Type>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]), t[4]->CastTo(args[4]));
        }
};

template<class Type>
AnyType initCSR_Op<Type>::operator()(Stack stack) const {
    DistributedCSR* ptA = GetAny<DistributedCSR*>((*A)(stack));
    KN<KN<long>>* ptR = GetAny<KN<KN<long>>*>((*R)(stack));
    KN<long>* ptO = GetAny<KN<long>*>((*O)(stack));
    KN<double>* ptD = GetAny<KN<double>*>((*D)(stack));
    long bs = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : 1;
    ptA->_A = new HpSchwarz<PetscScalar>;
    MatriceMorse<PetscScalar> *mA = static_cast<MatriceMorse<PetscScalar>*>(&(*GetAny<Matrice_Creuse<PetscScalar>*>((*K)(stack))->A));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    if(ptO && ptA) {
        HPDDM::MatrixCSR<PetscScalar>* dA = new HPDDM::MatrixCSR<PetscScalar>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->_A->Subdomain<PetscScalar>::initialize(dA, STL<long>(*ptO), *ptR, comm);
    }
    if(comm)
        PETSC_COMM_WORLD = *comm;
    ptA->_A->HpSchwarz<PetscScalar>::initialize(*ptD);
    if(!ptA->_num)
        ptA->_num = new unsigned int[ptA->_A->getDof()];
    double timing = MPI_Wtime();
    ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, ptA->_global);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- global numbering created (in " << MPI_Wtime() - timing << ")" << endl;
    timing = MPI_Wtime();
    ptA->_free = ptA->_A->distributedCSR(ptA->_num, ptA->_first, ptA->_last, ptA->_ia, ptA->_ja, ptA->_c);
    MatCreate(PETSC_COMM_WORLD, &(ptA->_petsc));
    if(bs > 1)
        MatSetBlockSize(ptA->_petsc, bs);
    MatSetSizes(ptA->_petsc, ptA->_last - ptA->_first, ptA->_last - ptA->_first, ptA->_global, ptA->_global);
    if(mpisize > 1) {
        MatSetType(ptA->_petsc, MATMPIAIJ);
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, ptA->_ia, ptA->_ja, ptA->_c);
    }
    else {
        MatSetType(ptA->_petsc, MATSEQAIJ);
        MatSeqAIJSetPreallocationCSR(ptA->_petsc, ptA->_ia, ptA->_ja, ptA->_c);
    }
    if(mA->symetrique)
        MatSetOption(ptA->_petsc, MAT_SYMMETRIC, PETSC_TRUE);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- global CSR created (in " << MPI_Wtime() - timing << ")" << endl;
    KN<PetscScalar>* rhs = nargs[2] ? GetAny<KN<PetscScalar>*>((*nargs[2])(stack)) : 0;
    if(rhs)
        ptA->_A->HPDDM::template Subdomain<PetscScalar>::exchange(*rhs);
    timing = MPI_Wtime();
    KSPCreate(PETSC_COMM_WORLD, &(ptA->_ksp));
    KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
    double eps = 1e-8;
    FFPetscOptionsGetReal(NULL, "-eps", &eps, NULL);
    int it = 100;
    FFPetscOptionsGetInt(NULL, "-iter", &it, NULL);
    KSPSetTolerances(ptA->_ksp, eps, PETSC_DEFAULT, PETSC_DEFAULT, it);
    KSPSetFromOptions(ptA->_ksp);
    KSPSetUp(ptA->_ksp);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- PETSc preconditioner built (in " << MPI_Wtime() - timing << ")" << endl;
    MatCreateVecs(ptA->_petsc, &(ptA->_x), nullptr);
    ptA->_A->setBuffer(1);
    return ptA;
}

template<class Type>
class setOptions_Op : public E_F0mps {
    public:
        Expression A;
        static const int n_name_param = 2;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        setOptions_Op(const basicAC_F0& args, Expression param1) : A(param1) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};

template<class Type>
basicAC_F0::name_and_type setOptions_Op<Type>::name_param[] = {
    {"sparams", &typeid(std::string*)},
    {"nearnullspace", &typeid(FEbaseArrayKn<PetscScalar>*)}
};

template<class Type>
class setOptions : public OneOperator {
    public:
        setOptions() : OneOperator(atype<long>(), atype<DistributedCSR*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new setOptions_Op<Type>(args, t[0]->CastTo(args[0]));
        }
};

template<class Type>
AnyType setOptions_Op<Type>::operator()(Stack stack) const {
    DistributedCSR* ptA = GetAny<DistributedCSR*>((*A)(stack));
    std::string* options = nargs[0] ? GetAny<std::string*>((*nargs[0])(stack)) : NULL;
    if(options) {
        std::vector<std::string> elems;
        std::stringstream ss(*options);
        std::string item;
        while (std::getline(ss, item, ' '))
            elems.push_back(item);
        int argc = elems.size() + 1;
        char** data = new char*[argc];
        data[0] = new char[options->size() + argc]();
        data[1] = data[0] + 1;
        for(int i = 0; i < argc - 1; ++i) {
            if(i > 0)
                data[i + 1] = data[i] + elems[i - 1].size() + 1;
            strcpy(data[i + 1], elems[i].c_str());
        }
        FFPetscOptionsInsert(&argc, &data, NULL);
        delete [] *data;
        delete [] data;
    }
    FEbaseArrayKn<PetscScalar>* ptNS = nargs[1] ? GetAny<FEbaseArrayKn<PetscScalar>*>((*nargs[1])(stack)) : 0;
    int dim = 0;
    if(ptNS)
        dim = ptNS->N;
    if(dim) {
        Vec x;
        MatCreateVecs(ptA->_petsc, &x, NULL);
        Vec* ns;
        VecDuplicateVecs(x, dim, &ns);
        for(unsigned short i = 0; i < dim; ++i) {
            PetscScalar* x;
            VecGetArray(ns[i], &x);
            ptA->_A->distributedVec<0>(ptA->_num, ptA->_first, ptA->_last, *(ptNS->get(i)), x, ptA->_A->getDof());
            VecRestoreArray(ns[i], &x);
        }
        PetscScalar* dots = new PetscScalar[dim];
        for(unsigned short i = 0; i < dim; ++i) {
            if(i > 0) {
                VecMDot(ns[i], i, ns, dots);
                for(int j = 0; j < i; ++j)
                    dots[j] *= -1.0;
                VecMAXPY(ns[i], i, dots, ns);
            }
            VecNormalize(ns[i], NULL);
        }
        delete [] dots;
        MatNullSpace sp;
        MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, dim, ns, &sp);
        MatSetNearNullSpace(ptA->_petsc, sp);
        MatNullSpaceDestroy(&sp);
        if(ns)
            VecDestroyVecs(dim, &ns);
        delete [] ns;
    }
    double timing = MPI_Wtime();
    KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
    double eps = 1e-8;
    FFPetscOptionsGetReal(NULL, "-eps", &eps, NULL);
    int it = 100;
    FFPetscOptionsGetInt(NULL, "-iter", &it, NULL);
    KSPSetTolerances(ptA->_ksp, eps, PETSC_DEFAULT, PETSC_DEFAULT, it);
    KSPSetFromOptions(ptA->_ksp);
    KSPSetUp(ptA->_ksp);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- PETSc preconditioner built (in " << MPI_Wtime() - timing << ")" << endl;
    return 0L;
}

template<class T, class U, class K>
class InvPETSc {
    public:
        const T t;
        const U u;
        InvPETSc(T v, U w) : t(v), u(w) {}
        void solve(U out) const {
            Vec y;
            double timing = MPI_Wtime();
            MatCreateVecs((*t)._petsc, nullptr, &y);
            K* x;
            VecGetArray((*t)._x, &x);
            (*t)._A->template distributedVec<0>((*t)._num, (*t)._first, (*t)._last, static_cast<PetscScalar*>(*u), x, (*t)._A->getDof());
            VecRestoreArray((*t)._x, &x);
            std::fill(static_cast<PetscScalar*>(*out), static_cast<PetscScalar*>(*out) + out->n, 0);
            timing = MPI_Wtime();
            KSPSolve((*t)._ksp, (*t)._x, y);
            if(verbosity > 0 && mpirank == 0)
                cout << " --- system solved with PETSc (in " << MPI_Wtime() - timing << ")" << endl;
            VecGetArray(y, &x);
            (*t)._A->template distributedVec<1>((*t)._num, (*t)._first, (*t)._last, *out, x, (*t)._A->getDof());
            VecRestoreArray(y, &x);
            VecDestroy(&y);
            (*t)._A->HPDDM::template Subdomain<K>::exchange(*out);
        };
        static U init(U Ax, InvPETSc<T, U, K> A) {
            A.solve(Ax);
            return Ax;
        }
};

template<class T, class U, class K>
class ProdPETSc {
    public:
        const T t;
        const U u;
        ProdPETSc(T v, U w) : t(v), u(w) {}
        void prod(U x) const {
            Vec z;
            Vec y;
            MatCreateVecs((*t)._petsc, &z, &y);
            K* w;
            VecGetArray(z, &w);
            (*t)._A->template distributedVec<0>(t->_num, t->_first, t->_last, *u, w, t->_A->getDof());
            VecRestoreArray(z, &w);
            MatMult(t->_petsc, z, y);
            VecGetArray(y, &w);
            (*t)._A->template distributedVec<1>(t->_num, t->_first, t->_last, *x, w, t->_A->getDof());
            VecRestoreArray(y, &w);
            VecDestroy(&y);
            VecDestroy(&z);
            (*t)._A->HPDDM::template Subdomain<K>::exchange(*x);
        }
        static U mv(U Ax, ProdPETSc<T, U, K> A) {
            A.prod(Ax);
            return Ax;
        }
};
}

static void Init_PETSc() {
    int argc = pkarg->n;
    char** argv = new char*[argc];
    for(int i = 0; i < argc; ++i)
        argv[i] = const_cast<char*>((*(*pkarg)[i].getap())->c_str());
    PetscInitialize(&argc, &argv, 0, "");
    delete [] argv;
    ff_atend(PETSc::finalizePETSc);
    Dcl_Type<PETSc::DistributedCSR*>(Initialize<PETSc::DistributedCSR>, Delete<PETSc::DistributedCSR>);
    zzzfff->Add("dmatrix", atype<PETSc::DistributedCSR*>());
    TheOperators->Add("<-", new OneOperator1_<long, PETSc::DistributedCSR*>(PETSc::initEmptyCSR));
    TheOperators->Add("<-", new PETSc::initCSR<PetscScalar>);
    Global.Add("set", "(", new PETSc::setOptions<PetscScalar>());
    addProd<PETSc::DistributedCSR, PETSc::ProdPETSc, KN<PetscScalar>, PetscScalar>();
    addInv<PETSc::DistributedCSR, PETSc::InvPETSc, KN<PetscScalar>, PetscScalar>();
}

LOADFUNC(Init_PETSc)
