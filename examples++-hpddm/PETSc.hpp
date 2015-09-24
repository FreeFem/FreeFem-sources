#include <mpi.h>
#include <iostream>

extern KN<String>* pkarg;

#include "petsc.h"
#if PETSC_VERSION_LT(3,6,0)
#define MatCreateVecs MatGetVecs
#endif

class DistributedCSR {
    public:
        HpSchwarz<>*                _A;
        Mat                     _petsc;
        Vec                         _x;
        KSP                       _ksp;
        unsigned int*             _num;
        int*                       _ia;
        int*                       _ja;
        double*                     _c;
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
    {"rhs", &typeid(KN<double>*)}
};

template<class Type>
class initCSR : public OneOperator {
    public:
        initCSR() : OneOperator(atype<DistributedCSR*>(), atype<DistributedCSR*>(), atype<Matrice_Creuse<double>*>(), atype<KN<long>*>(), atype<KN<KN<long>>*>(), atype<KN<double>*>()) { }

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
    ptA->_A = new HpSchwarz<>;
    MatriceMorse<double> *mA = static_cast<MatriceMorse<double>*>(&(*GetAny<Matrice_Creuse<double>*>((*K)(stack))->A));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    if(ptO && ptA) {
        HPDDM::MatrixCSR<double>* dA = new HPDDM::MatrixCSR<double>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->_A->Subdomain<double>::initialize(dA, STL<long>(*ptO), *ptR, comm);
    }
    if(comm)
        PETSC_COMM_WORLD = *comm;
    ptA->_A->HpSchwarz<>::initialize(*ptD);
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
    KN<double>* rhs = nargs[2] ? GetAny<KN<double>*>((*nargs[2])(stack)) : 0;
    if(rhs)
        ptA->_A->HPDDM::template Subdomain<double>::exchange(*rhs);
    timing = MPI_Wtime();
    KSPCreate(PETSC_COMM_WORLD, &(ptA->_ksp));
    KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
    double eps = 1e-8;
    PetscOptionsGetReal(NULL, "-eps", &eps, NULL);
    int it = 100;
    PetscOptionsGetInt(NULL, "-iter", &it, NULL);
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
    {"nearnullspace", &typeid(FEbaseArrayKn<double>*)}
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
        PetscOptionsInsert(&argc, &data, NULL);
        delete [] *data;
        delete [] data;
    }
    FEbaseArrayKn<double>* ptNS = nargs[1] ? GetAny<FEbaseArrayKn<double>*>((*nargs[1])(stack)) : 0;
    int dim = 0;
    if(ptNS)
        dim = ptNS->N;
    if(dim) {
        Vec x;
        MatCreateVecs(ptA->_petsc, &x, NULL);
        Vec* ns;
        VecDuplicateVecs(x, dim, &ns);
        for(unsigned short i = 0; i < dim; ++i) {
            double* x;
            VecGetArray(ns[i], &x);
            ptA->_A->distributedVec<0>(ptA->_num, ptA->_first, ptA->_last, *(ptNS->get(i)), x, ptA->_A->getDof());
            VecRestoreArray(ns[i], &x);
        }
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
    PetscOptionsGetReal(NULL, "-eps", &eps, NULL);
    int it = 100;
    PetscOptionsGetInt(NULL, "-iter", &it, NULL);
    KSPSetTolerances(ptA->_ksp, eps, PETSC_DEFAULT, PETSC_DEFAULT, it);
    KSPSetFromOptions(ptA->_ksp);
    KSPSetUp(ptA->_ksp);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- PETSc preconditioner built (in " << MPI_Wtime() - timing << ")" << endl;
    return 0L;
}

class DistributedCSR_inv  { public:
  DistributedCSR* A;
  DistributedCSR_inv(DistributedCSR* B) : A(B) { assert(A); }
  operator DistributedCSR& () const { return *A; }
  operator DistributedCSR* () const { return A; }
};

class OneBinaryOperatorPETSc : public OneOperator {
    public:
        OneBinaryOperatorPETSc() : OneOperator(atype<DistributedCSR_inv>(), atype<DistributedCSR*>(), atype<long>()) { }
        E_F0* code(const basicAC_F0 & args) const {
            Expression p = args[1];
            if(!p->EvaluableWithOutStack())
                CompileError("A^p, The p must be a constant == -1, sorry");
            long pv = GetAny<long>((*p)(NullStack));
            if(pv != -1) {
                char buf[100];
                sprintf(buf, "A^%ld, The pow must be == -1, sorry", pv);
                CompileError(buf);
            }
            return new E_F_F0<DistributedCSR_inv, DistributedCSR*>(Build<DistributedCSR_inv, DistributedCSR*>, t[0]->CastTo(args[0]));
        }
};


template<class T, class U>
class Inv {
    public:
        const T t;
        const U u;
        Inv(T v, U w) : t(v), u(w) {}
        void solve(U out) const {
            Vec y;
            double timing = MPI_Wtime();
            MatCreateVecs((*t)._petsc, nullptr, &y);
            double* x;
            VecGetArray((*t)._x, &x);
            (*t)._A->template distributedVec<0>((*t)._num, (*t)._first, (*t)._last, static_cast<double*>(*u), x, (*t)._A->getDof());
            VecRestoreArray((*t)._x, &x);
            std::fill(static_cast<double*>(*out), static_cast<double*>(*out) + out->n, 0);
            timing = MPI_Wtime();
            KSPSolve((*t)._ksp, (*t)._x, y);
            if(verbosity > 0 && mpirank == 0)
                cout << " --- system solved with PETSc (in " << MPI_Wtime() - timing << ")" << endl;
            VecGetArray(y, &x);
            (*t)._A->template distributedVec<1>((*t)._num, (*t)._first, (*t)._last, *out, x, (*t)._A->getDof());
            VecRestoreArray(y, &x);
            VecDestroy(&y);
            (*t)._A->HPDDM::template Subdomain<double>::exchange(*out);
        };
};

KN<double>* InvPETSc(KN<double>* Ax, Inv<DistributedCSR_inv, KN<double>*> A) {
    A.solve(Ax);
    return Ax;
}

template<class T, class U>
class GMV {
    public:
        const T t;
        const U u;
        GMV(T v, U w) : t(v), u(w) {}
        void prod(U x) const;
};

template<class U>
class GMV<HpSchwarz<>*, U> {
    public:
        const HpSchwarz<>* t;
        const U u;
        GMV(HpSchwarz<>* v, U w) : t(v), u(w) {}
        void prod(U x) const { t->GMV(*(this->u), *x); };
};

template<class R, class A, class B> R Build(A a, B b) {
    return R(a, b);
}

template<class K>
KN<K>* GlobalMV(KN<K>* Ax, GMV<HpSchwarz<>*, KN<K>*> A) {
    A.prod(Ax);
    return Ax;
}

template<class U>
class GMV<DistributedCSR*, U> {
    public:
        const DistributedCSR* t;
        const U u;
        GMV(DistributedCSR* v, U w) : t(v), u(w) {}
        void prod(U x) const {
            Vec z;
            Vec y;
            MatCreateVecs(t->_petsc, &z, &y);
            double* w;
            VecGetArray(z, &w);
            t->_A->distributedVec<0>(t->_num, t->_first, t->_last, *u, w, t->_A->getDof());
            VecRestoreArray(z, &w);
            MatMult(t->_petsc, z, y);
            VecGetArray(y, &w);
            t->_A->template distributedVec<1>(t->_num, t->_first, t->_last, *x, w, t->_A->getDof());
            VecRestoreArray(y, &w);
            VecDestroy(&y);
            VecDestroy(&z);
            t->_A->Subdomain<double>::exchange(*x);
        };
};
KN<double>* GlobalMV(KN<double>* Ax, GMV<DistributedCSR*, KN<double>*> A) {
    A.prod(Ax);
    return Ax;
}
