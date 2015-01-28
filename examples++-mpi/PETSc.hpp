#include <mpi.h>
#include <iostream>

#include "petsc.h"

class DistributedCSR {
    public:
        HpSchwarz<>*                _A;
        Mat                     _petsc;
        Vec                         _x;
        unsigned int*             _num;
        int*                       _ia;
        int*                       _ja;
        double*                     _c;
        unsigned int            _first;
        unsigned int             _last;
        unsigned int           _global;
        bool                     _free;
        DistributedCSR() : _A(), _num(), _ia(), _ja(), _c(), _first(), _last() {} ;
        ~DistributedCSR() {
            VecDestroy(&_x);
            MatDestroy(&_petsc);
            delete [] _num;
            if(_free) {
                delete []  _ia;
                delete []  _ja;
                delete []   _c;
            }
            delete _A;
        }
};

long initEmptyCSR(DistributedCSR* const& A) {
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
        static const int n_name_param = 1;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        initCSR_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3, Expression param4, Expression param5) : A(param1), K(param2), O(param3), R(param4), D(param5) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};

template<class Type>
basicAC_F0::name_and_type initCSR_Op<Type>::name_param[] = {
    {"bs", &typeid(long)}
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
    Matrice_Creuse<double>* ptK = GetAny<Matrice_Creuse<double>*>((*K)(stack));
    MatriceMorse<double> *mK = static_cast<MatriceMorse<double>*>(&(*ptK->A));
    KN<KN<long>>* ptR = GetAny<KN<KN<long>>*>((*R)(stack));
    KN<long>* ptO = GetAny<KN<long>*>((*O)(stack));
    KN<double>* ptD = GetAny<KN<double>*>((*D)(stack));
    long bs = nargs[0] ? GetAny<long>((*nargs[0])(stack)) : 1;
    ptA->_A = new HpSchwarz<>;
    std::vector<KN<long>*> vec(ptR->n);
    for(unsigned short i = 0; i < ptR->n; ++i)
        vec[i] = &(ptR->operator[](i));
    if(ptO && ptA) {
        MatriceMorse<double>* mA = static_cast<MatriceMorse<double>*>(&(*ptK->A));
        HPDDM::MatrixCSR<double>* dA = new HPDDM::MatrixCSR<double>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->_A->Subdomain<double>::initialize(dA, static_cast<long*>(*ptO), static_cast<long*>(*ptO) + ptO->N(), vec);
    }
    ptA->_A->HpSchwarz<>::initialize(*ptD);
    if(!ptA->_num)
        ptA->_num = new unsigned int[ptA->_A->getDof()];
    double timing = MPI_Wtime();
    ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, ptA->_global);
    if(mpirank == 0)
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
    if(mK->symetrique)
        MatSetOption(ptA->_petsc, MAT_SYMMETRIC, PETSC_TRUE);
    if(mpirank == 0)
        cout << " --- global CSR created (in " << MPI_Wtime() - timing << ")" << endl;
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
                data[i + 1] = data[i] + elems[i - 1].size()+1;
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
            ptA->_A->distributedRHS<0>(ptA->_num, ptA->_first, ptA->_last, *(ptNS->get(i)), x, ptA->_A->getDof());
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
    return 0L;
}

long distributedNumbering(DistributedCSR* const& A, KN<double>* const& vec) {
    if(!A->_num)
        A->_num = new unsigned int[vec->n];
    double timing = MPI_Wtime();
    A->_A->distributedNumbering(A->_num, A->_first, A->_last, A->_global);
    if(mpirank == 0)
        cout << " --- global numbering created (in " << MPI_Wtime() - timing << ")" << endl;
    for(int i = 0; i < vec->n; ++i)
        vec->operator[](i) = A->_num[i];
    return 0;
}
long renumberCSR(DistributedCSR* const& A, FEbaseArrayKn<double>* const& nullspace) {
    double timing = MPI_Wtime();
    A->_free = A->_A->distributedCSR(A->_num, A->_first, A->_last, A->_ia, A->_ja, A->_c);
    MatCreate(PETSC_COMM_WORLD, &(A->_petsc));
    MatSetSizes(A->_petsc, A->_last - A->_first, A->_last - A->_first, A->_global, A->_global);
    MatSetBlockSize(A->_petsc, (nullspace->N == 3) ? 2 : ((nullspace->N == 6) ? 3 : 1));
    MatSetType(A->_petsc, MATMPIAIJ);
    MatMPIAIJSetPreallocationCSR(A->_petsc, A->_ia, A->_ja, A->_c);
    MatSetOption(A->_petsc, MAT_SPD, PETSC_TRUE);
    if(mpirank == 0)
        cout << " --- global CSR created (in " << MPI_Wtime() - timing << ")" << endl;
    int dim = nullspace->N;
    if(dim) {
        Vec x;
        MatCreateVecs(A->_petsc, &x, NULL);
        timing = MPI_Wtime();
        Vec* ns;
        VecDuplicateVecs(x, dim, &ns);
        for(unsigned short i = 0; i < dim; ++i) {
            double* x;
            VecGetArray(ns[i], &x);
            A->_A->distributedRHS<0>(A->_num, A->_first, A->_last, *(nullspace->get(i)), x, A->_A->getDof());
            VecRestoreArray(ns[i], &x);
        }
        MatNullSpace sp;
        MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, dim, ns, &sp);
        MatSetNearNullSpace(A->_petsc, sp);
        if(mpirank == 0)
            cout << " --- near null space set (in " << MPI_Wtime() - timing << ")" << endl;
        MatNullSpaceDestroy(&sp);
        if(ns)
            VecDestroyVecs(dim, &ns);
        delete [] ns;
#if 0
        PetscViewer viewer;
        const Vec* vecs;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "A", FILE_MODE_WRITE, &viewer);
        MatView(A->_petsc, viewer);
        MatNullSpace mnull;
        MatGetNearNullSpace(A->_petsc, &mnull);
        MatNullSpaceGetVecs(mnull, NULL, &dim, &vecs);
        for(unsigned short i = 0; i < dim; ++i)
            VecView(vecs[i], viewer);
        PetscViewerDestroy(&viewer);
#endif
    }
    return 0;
}
long solvePETSc(DistributedCSR* const& A, KN<double>* const& in) {
    KSP ksp;
    Vec y;
    double timing = MPI_Wtime();
    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A->_petsc, A->_petsc);
    double eps = 1e-8;
    PetscOptionsGetReal(NULL, "-eps", &eps, NULL);
    int it = 100;
    PetscOptionsGetInt(NULL, "-iter", &it, NULL);
    KSPSetTolerances(ksp, eps, PETSC_DEFAULT, PETSC_DEFAULT, it);
    KSPSetFromOptions(ksp);
    KSPSetUp(ksp);
    if(mpirank == 0)
        cout << " --- PETSc preconditioner built (in " << MPI_Wtime() - timing << ")" << endl;
    MatCreateVecs(A->_petsc, &(A->_x), &y);
    double* x;
    VecGetArray(A->_x, &x);
    A->_A->distributedRHS<0>(A->_num, A->_first, A->_last, *in, x, A->_A->getDof());
    VecRestoreArray(A->_x, &x);
    std::fill(static_cast<double*>(*in), static_cast<double*>(*in) + in->n, 0);
    timing = MPI_Wtime();
    KSPSolve(ksp, A->_x, y);
    if(mpirank == 0)
        cout << " --- system solved with PETSc (in " << MPI_Wtime() - timing << ")" << endl;
    VecGetArray(y, &x);
    A->_A->distributedRHS<1>(A->_num, A->_first, A->_last, *in, x, A->_A->getDof());
#if 0
    Vec z;
    MatCreateVecs(A->_petsc, &z, NULL);
    MatMult(A->_petsc, y, z);
    VecAXPY(z, -1, A->_x);
    double err[2];
    VecNorm(z, NORM_2, err);
    VecNorm(A->_x, NORM_2, err + 1);
    if(mpirank == 0)
        std::cout << " --- PETSc error = " << err[0] << " / " << err[1] << std::endl;
#endif
    VecRestoreArray(y, &x);
    VecDestroy(&y);
    KSPDestroy(&ksp);
    A->_A->Subdomain<double>::exchange(*in);
    return 0;
}

class DistributedCSR_inv  { public:
  DistributedCSR* A;
  DistributedCSR_inv(DistributedCSR* AA) : A(AA) {assert(A);}
  operator DistributedCSR& () const {return *A;}
  operator DistributedCSR* () const {return A;}
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
            KSP ksp;
            Vec y;
            double timing = MPI_Wtime();
            KSPCreate(PETSC_COMM_WORLD, &ksp);
            KSPSetOperators(ksp, (*t)._petsc, (*t)._petsc);
            double eps = 1e-8;
            PetscOptionsGetReal(NULL, "-eps", &eps, NULL);
            int it = 100;
            PetscOptionsGetInt(NULL, "-iter", &it, NULL);
            KSPSetTolerances(ksp, eps, PETSC_DEFAULT, PETSC_DEFAULT, it);
            KSPSetFromOptions(ksp);
            KSPSetUp(ksp);
            if(mpirank == 0)
                cout << " --- PETSc preconditioner built (in " << MPI_Wtime() - timing << ")" << endl;
            MatCreateVecs((*t)._petsc, &((*t)._x), &y);
            double* x;
            VecGetArray((*t)._x, &x);
            (*t)._A->template distributedRHS<0>((*t)._num, (*t)._first, (*t)._last, static_cast<double*>(*u), x, (*t)._A->getDof());
            VecRestoreArray((*t)._x, &x);
            std::fill(static_cast<double*>(*out), static_cast<double*>(*out) + out->n, 0);
            timing = MPI_Wtime();
            KSPSolve(ksp, (*t)._x, y);
            if(mpirank == 0)
                cout << " --- system solved with PETSc (in " << MPI_Wtime() - timing << ")" << endl;
            VecGetArray(y, &x);
            (*t)._A->template distributedRHS<1>((*t)._num, (*t)._first, (*t)._last, *out, x, (*t)._A->getDof());
            VecRestoreArray(y, &x);
            VecDestroy(&y);
            KSPDestroy(&ksp);
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
            t->_A->distributedRHS<0>(t->_num, t->_first, t->_last, *u, w, t->_A->getDof());
            VecRestoreArray(z, &w);
            MatMult(t->_petsc, z, y);
            VecGetArray(y, &w);
            t->_A->template distributedRHS<1>(t->_num, t->_first, t->_last, *x, w, t->_A->getDof());
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
