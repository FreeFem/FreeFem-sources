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
template<class HpddmType>
class DistributedCSR {
    public:
        HpddmType*                  _A;
        Mat                     _petsc;
        Vec                         _x;
        ISLocalToGlobalMapping   _rmap;
        VecScatter            _scatter;
        Vec                     _isVec;
        KSP                       _ksp;
        unsigned int*             _num;
        int*                       _ia;
        int*                       _ja;
        PetscScalar*                _c;
        unsigned int            _first;
        unsigned int             _last;
        unsigned int           _global;
        bool                     _free;
        DistributedCSR() : _A(), _num(), _ia(), _ja(), _c(), _first(), _last(), _free(false) { };
        ~DistributedCSR() {
            if(_A) {
                if(std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value) {
                    VecDestroy(&_x);
                    MatDestroy(&_petsc);
                    KSPDestroy(&_ksp);
                }
                else {
                    ISLocalToGlobalMappingDestroy(&_rmap);
                    VecDestroy(&_isVec);
                    VecScatterDestroy(&_scatter);
                }
                delete [] _num;
                _num = nullptr;
                if(_free) {
                    delete []  _ia;
                    delete []  _ja;
                    delete []   _c;
                    _ia = _ja = nullptr;
                    _c = nullptr;
                }
                if(std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value)
                    _A->clearBuffer();
                delete _A;
                _A = nullptr;
            }
        }
};
template<class HpddmType, typename std::enable_if<std::is_same<HpddmType, DistributedCSR<HpSchwarz<PetscScalar>>>::value>::type* = nullptr>
void initPETScStructure(HpddmType* ptA, MatriceMorse<PetscScalar>* mA, long& bs, KN<typename std::conditional<std::is_same<HpddmType, DistributedCSR<HpSchwarz<PetscScalar>>>::value, double, long>::type>* ptD, KN<PetscScalar>* rhs) {
    ptA->_A->initialize(*ptD);
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
    ptA->_A->setBuffer();
    if(verbosity > 0 && mpirank == 0)
        cout << " --- global CSR created (in " << MPI_Wtime() - timing << ")" << endl;
    if(rhs)
        ptA->_A->exchange(*rhs);
}
template<class HpddmType, typename std::enable_if<!std::is_same<HpddmType, DistributedCSR<HpSchwarz<PetscScalar>>>::value>::type* = nullptr>
void initPETScStructure(HpddmType* ptA, MatriceMorse<PetscScalar>* mA, long& bs, KN<typename std::conditional<std::is_same<HpddmType, DistributedCSR<HpSchwarz<PetscScalar>>>::value, double, long>::type>* ptD, KN<PetscScalar>* rhs) {
    const HPDDM::MatrixCSR<PetscScalar>* M = ptA->_A->getMatrix();
    if(!M->_sym)
        std::cout << "Please assemble a symmetric CSR" << std::endl;
    double timing = MPI_Wtime();
    ptA->_A->template renumber<false>(STL<long>(*ptD), nullptr);
    ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, ptA->_global);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- global numbering created (in " << MPI_Wtime() - timing << ")" << endl;
    timing = MPI_Wtime();
    PetscInt* indices;
    PetscMalloc(sizeof(PetscInt) * M->_n / bs, &indices);
    for(unsigned int i = 0; i < M->_n; i += bs)
        indices[i / bs] = ptA->_num[i] / bs;
    ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, bs, M->_n / bs, indices, PETSC_OWN_POINTER, &(ptA->_rmap));
    MatCreateIS(PETSC_COMM_WORLD, bs, PETSC_DECIDE, PETSC_DECIDE, ptA->_global, ptA->_global, ptA->_rmap, NULL, &(ptA->_petsc));
    Mat local;
    MatISGetLocalMat(ptA->_petsc, &local);
    MatSetType(local, MATSEQSBAIJ);
    std::vector<std::vector<std::pair<int, PetscScalar>>> transpose(M->_n);
    for(int i = 0; i < transpose.size(); ++i)
        for(int j = M->_ia[i]; j < M->_ia[i + 1]; ++j) {
            transpose[M->_ja[j]].emplace_back(i, M->_a[j]);
            if(bs > 1 && (i - M->_ja[j] <= (i % bs)) && M->_ja[j] != i)
                transpose[i].emplace_back(M->_ja[j], M->_a[j]);
        }
    int nnz = 0;
    for(int i = 0; i < transpose.size(); ++i) {
        std::sort(transpose[i].begin(), transpose[i].end(), [](const std::pair<int, PetscScalar>& lhs, const std::pair<int, PetscScalar>& rhs) { return lhs.first < rhs.first; });
        nnz += transpose[i].size();
    }
    int* ia = new int[M->_n / bs + 1];
    int* ja = new int[nnz / (bs * bs)];
    PetscScalar* a = new PetscScalar[nnz];
    ia[0] = 0;
    for(int i = 0; i < transpose.size(); ++i) {
        for(int j = 0; j < transpose[i].size(); ++j) {
            if(i % bs == 0 && j % bs == 0)
                ja[ia[i / bs] + j / bs] = transpose[i][j].first / bs;
            a[ia[i / bs] * (bs * bs) + j % bs + (j / bs) * (bs * bs) + (i % bs) * bs] = transpose[i][j].second;
        }
        if(i % bs == 0)
            ia[i / bs + 1] = ia[i / bs] + transpose[i].size() / bs;
    }
    MatSeqSBAIJSetPreallocationCSR(local, bs, ia, ja, a);
    MatAssemblyBegin(ptA->_petsc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(ptA->_petsc, MAT_FINAL_ASSEMBLY);
    delete [] a;
    delete [] ja;
    delete [] ia;
    IS to, from;
    PetscInt nr;
    Vec rglobal;
    ISLocalToGlobalMappingGetSize(ptA->_rmap, &nr);
    ISCreateStride(PETSC_COMM_SELF, nr, 0, 1, &to);
    ISLocalToGlobalMappingApplyIS(ptA->_rmap, to, &from);
    MatCreateVecs(ptA->_petsc, &rglobal, NULL);
    VecCreate(PETSC_COMM_SELF, &(ptA->_isVec));
    VecSetType(ptA->_isVec, VECMPI);
    VecSetSizes(ptA->_isVec, PETSC_DECIDE, nr);
    VecScatterCreate(rglobal, from, ptA->_isVec, to, &(ptA->_scatter));
    VecDestroy(&rglobal);
    ISDestroy(&from);
    ISDestroy(&to);
    // ISLocalToGlobalMappingView(ptA->_rmap, PETSC_VIEWER_STDOUT_WORLD);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- global CSR created (in " << MPI_Wtime() - timing << ")" << endl;
}
template<class Type>
long globalNumbering(Type* const& A, KN<long>* const& numbering) {
    numbering->resize(A->_A->getMatrix()->_n);
    if(A->_num)
        for(int i = 0; i < numbering->n; ++i)
            numbering->operator[](i) = A->_num[i];
    return 0L;
}
template<class Type, class K>
long originalNumbering(Type* const& A, KN<K>* const& in, KN<long>* const& interface) {
    A->_A->originalNumbering(STL<long>(*interface), *in);
    return 0;
}

void finalizePETSc() {
    PETSC_COMM_WORLD = MPI_COMM_WORLD;
    PetscFinalize();
}

template<class Type>
long initEmptyCSR(Type* const&) {
    return 0;
}

template<class HpddmType>
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

template<class HpddmType>
basicAC_F0::name_and_type initCSR_Op<HpddmType>::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"bs", &typeid(long)},
    {"rhs", &typeid(KN<PetscScalar>*)}
};

template<class HpddmType>
class initCSR : public OneOperator {
    public:
        initCSR() : OneOperator(atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<Matrice_Creuse<PetscScalar>*>(), atype<KN<long>*>(), atype<KN<KN<long>>*>(), atype<KN<typename std::conditional<std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value, double, long>::type>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new initCSR_Op<HpddmType>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]), t[4]->CastTo(args[4]));
        }
};

template<class HpddmType>
AnyType initCSR_Op<HpddmType>::operator()(Stack stack) const {
    DistributedCSR<HpddmType>* ptA = GetAny<DistributedCSR<HpddmType>*>((*A)(stack));
    KN<KN<long>>* ptR = GetAny<KN<KN<long>>*>((*R)(stack));
    KN<long>* ptO = GetAny<KN<long>*>((*O)(stack));
    KN<typename std::conditional<std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value, double, long>::type>* ptD = GetAny<KN<typename std::conditional<std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value, double, long>::type>*>((*D)(stack));
    long bs = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : 1;
    MatriceMorse<PetscScalar> *mA = static_cast<MatriceMorse<PetscScalar>*>(&(*GetAny<Matrice_Creuse<PetscScalar>*>((*K)(stack))->A));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    KN<PetscScalar>* rhs = nargs[2] ? GetAny<KN<PetscScalar>*>((*nargs[2])(stack)) : 0;
    ptA->_A = new HpddmType;
    if(comm)
        PETSC_COMM_WORLD = *comm;
    if(ptO && ptA) {
        HPDDM::MatrixCSR<PetscScalar>* dA = new HPDDM::MatrixCSR<PetscScalar>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->_A->HPDDM::template Subdomain<PetscScalar>::initialize(dA, STL<long>(*ptO), *ptR, comm);
    }
    if(!ptA->_num)
        ptA->_num = new unsigned int[ptA->_A->getMatrix()->_n];
    initPETScStructure(ptA, mA, bs, ptD, rhs);
    if(ptO && ptA) {
        mA->a = ptA->_A->getMatrix()->_a;
        mA->lg = ptA->_A->getMatrix()->_ia;
        mA->cl = ptA->_A->getMatrix()->_ja;
    }
    double timing = MPI_Wtime();
    KSPCreate(PETSC_COMM_WORLD, &(ptA->_ksp));
    KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
    double eps = 1.0e-8;
    FFPetscOptionsGetReal(NULL, "-eps", &eps, NULL);
    int it = 100;
    FFPetscOptionsGetInt(NULL, "-iter", &it, NULL);
    KSPSetTolerances(ptA->_ksp, eps, PETSC_DEFAULT, PETSC_DEFAULT, it);
    KSPSetFromOptions(ptA->_ksp);
    KSPSetUp(ptA->_ksp);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- PETSc preconditioner built (in " << MPI_Wtime() - timing << ")" << endl;
    MatCreateVecs(ptA->_petsc, &(ptA->_x), nullptr);
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
        setOptions() : OneOperator(atype<long>(), atype<Type*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new setOptions_Op<Type>(args, t[0]->CastTo(args[0]));
        }
};

template<class Type>
AnyType setOptions_Op<Type>::operator()(Stack stack) const {
    Type* ptA = GetAny<Type*>((*A)(stack));
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
    if(std::is_same<Type, PETSc::DistributedCSR<HpSchwarz<PetscScalar>>>::value) {
        FEbaseArrayKn<PetscScalar>* ptNS = nargs[1] ? GetAny<FEbaseArrayKn<PetscScalar>*>((*nargs[1])(stack)) : 0;
        int dim = ptNS ? ptNS->N : 0;
        if(dim) {
            Vec x;
            MatCreateVecs(ptA->_petsc, &x, NULL);
            Vec* ns;
            VecDuplicateVecs(x, dim, &ns);
            for(unsigned short i = 0; i < dim; ++i) {
                PetscScalar* x;
                VecGetArray(ns[i], &x);
                ptA->_A->template distributedVec<0>(ptA->_num, ptA->_first, ptA->_last, *(ptNS->get(i)), x, ptA->_A->getDof());
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
    static_assert(std::is_same<K, PetscScalar>::value, "Wrong types");
    public:
        const T t;
        const U u;
        InvPETSc(T v, U w) : t(v), u(w) {}
        void solve(U out) const {
            Vec y;
            double timing = MPI_Wtime();
            MatCreateVecs((*t)._petsc, nullptr, &y);
            PetscScalar* x;
            if(std::is_same<typename std::remove_reference<decltype(*t.A->_A)>::type, HpSchwarz<PetscScalar>>::value) {
                VecGetArray((*t)._x, &x);
                (*t)._A->template distributedVec<0>((*t)._num, (*t)._first, (*t)._last, static_cast<PetscScalar*>(*u), x, (*t)._A->getDof());
                VecRestoreArray((*t)._x, &x);
            }
            else {
                PetscScalar zero = 0.0;
                VecSet((*t)._x, zero);
#if 0
                Mat_IS* is = (Mat_IS*)(*t)._petsc->data;
                VecGetArray(is->y, &x);
                std::copy_n(static_cast<PetscScalar*>(*u), (*t)._A->getMatrix()->_n, x);
                VecRestoreArray(is->y, &x);
                VecScatterBegin(is->rctx,is->y,(*t)._x,ADD_VALUES,SCATTER_REVERSE);
                VecScatterEnd(is->rctx,is->y,(*t)._x,ADD_VALUES,SCATTER_REVERSE);
#else
                VecGetArray((*t)._isVec, &x);
                std::copy_n(static_cast<PetscScalar*>(*u), (*t)._A->getMatrix()->_n, x);
                VecRestoreArray((*t)._isVec, &x);
                VecScatterBegin((*t)._scatter, (*t)._isVec, (*t)._x, ADD_VALUES, SCATTER_REVERSE);
                VecScatterEnd((*t)._scatter, (*t)._isVec, (*t)._x, ADD_VALUES, SCATTER_REVERSE);
#endif
            }
            timing = MPI_Wtime();
            KSPSolve((*t)._ksp, (*t)._x, y);
            if(verbosity > 0 && mpirank == 0)
                cout << " --- system solved with PETSc (in " << MPI_Wtime() - timing << ")" << endl;
            if(std::is_same<typename std::remove_reference<decltype(*t.A->_A)>::type, HpSchwarz<PetscScalar>>::value) {
                VecGetArray(y, &x);
                (*t)._A->template distributedVec<1>((*t)._num, (*t)._first, (*t)._last, *out, x, (*t)._A->getDof());
                VecRestoreArray(y, &x);
            }
            else {
#if 0
                Mat_IS* is = (Mat_IS*)(*t)._petsc->data;
                VecScatterBegin(is->rctx,y,is->y,INSERT_VALUES,SCATTER_FORWARD);
                VecScatterEnd(is->rctx,y,is->y,INSERT_VALUES,SCATTER_FORWARD);
                VecGetArray(is->y, &x);
                std::copy_n(x, (*t)._A->getMatrix()->_n, (PetscScalar*)*out);
                VecRestoreArray(is->y, &x);
#else
                VecScatterBegin((*t)._scatter, y, (*t)._isVec, INSERT_VALUES, SCATTER_FORWARD);
                VecScatterEnd((*t)._scatter, y, (*t)._isVec, INSERT_VALUES, SCATTER_FORWARD);
                VecGetArray((*t)._isVec, &x);
                std::copy_n(x, (*t)._A->getMatrix()->_n, (PetscScalar*)*out);
                VecRestoreArray((*t)._isVec, &x);
#endif
            }
            VecDestroy(&y);
            if(std::is_same<typename std::remove_reference<decltype(*t.A->_A)>::type, HpSchwarz<PetscScalar>>::value)
                (*t)._A->HPDDM::template Subdomain<PetscScalar>::exchange(*out);
        };
        static U init(U Ax, InvPETSc<T, U, K> A) {
            A.solve(Ax);
            return Ax;
        }
};

template<class T, class U, class K>
class ProdPETSc {
    static_assert(std::is_same<K, PetscScalar>::value, "Wrong types");
    public:
        const T t;
        const U u;
        ProdPETSc(T v, U w) : t(v), u(w) {}
        void prod(U x) const {
            Vec z;
            Vec y;
            MatCreateVecs((*t)._petsc, &z, &y);
            PetscScalar* w;
            VecGetArray(z, &w);
            (*t)._A->template distributedVec<0>(t->_num, t->_first, t->_last, *u, w, t->_A->getDof());
            VecRestoreArray(z, &w);
            MatMult(t->_petsc, z, y);
            VecGetArray(y, &w);
            (*t)._A->template distributedVec<1>(t->_num, t->_first, t->_last, *x, w, t->_A->getDof());
            VecRestoreArray(y, &w);
            VecDestroy(&y);
            VecDestroy(&z);
            (*t)._A->HPDDM::template Subdomain<PetscScalar>::exchange(*x);
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
    Dcl_Type<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>*>(Initialize<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>>, Delete<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>>);
    zzzfff->Add("dmatrix", atype<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>*>());
    TheOperators->Add("<-", new OneOperator1_<long, PETSc::DistributedCSR<HpSchwarz<PetscScalar>>*>(PETSc::initEmptyCSR<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>>));
    TheOperators->Add("<-", new PETSc::initCSR<HpSchwarz<PetscScalar>>);
    Global.Add("set", "(", new PETSc::setOptions<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>>());
    addProd<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PETSc::ProdPETSc, KN<PetscScalar>, PetscScalar>();
    addInv<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PETSc::InvPETSc, KN<PetscScalar>, PetscScalar>();

    Dcl_Type<PETSc::DistributedCSR<HpBdd<PetscScalar>>*>(Initialize<PETSc::DistributedCSR<HpBdd<PetscScalar>>>, Delete<PETSc::DistributedCSR<HpBdd<PetscScalar>>>);
    zzzfff->Add("dbddc", atype<PETSc::DistributedCSR<HpBdd<PetscScalar>>*>());
    TheOperators->Add("<-", new OneOperator1_<long, PETSc::DistributedCSR<HpBdd<PetscScalar>>*>(PETSc::initEmptyCSR<PETSc::DistributedCSR<HpBdd<PetscScalar>>>));
    TheOperators->Add("<-", new PETSc::initCSR<HpBdd<PetscScalar>>);

    Global.Add("globalNumbering", "(", new OneOperator2_<long, PETSc::DistributedCSR<HpSchwarz<PetscScalar>>*, KN<long>*>(PETSc::globalNumbering<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>>));
    Global.Add("globalNumbering", "(", new OneOperator2_<long, PETSc::DistributedCSR<HpBdd<PetscScalar>>*, KN<long>*>(PETSc::globalNumbering<PETSc::DistributedCSR<HpBdd<PetscScalar>>>));
    Global.Add("originalNumbering", "(", new OneOperator3_<long, PETSc::DistributedCSR<HpBdd<PetscScalar>>*, KN<PetscScalar>*, KN<long>*>(PETSc::originalNumbering));
    Global.Add("set", "(", new PETSc::setOptions<PETSc::DistributedCSR<HpBdd<PetscScalar>>>());
    addInv<PETSc::DistributedCSR<HpBdd<PetscScalar>>, PETSc::InvPETSc, KN<PetscScalar>>();
}

LOADFUNC(Init_PETSc)
