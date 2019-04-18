#include "petsc.h"

#include "PETSc.hpp"

typedef PETSc::DistributedCSR<HpSchwarz<PetscScalar>> Dmat;
typedef PETSc::DistributedCSR<HpSchwarz<PetscReal>> DmatR;
typedef PETSc::DistributedCSR<HpSchwarz<PetscComplex>> DmatC;
typedef PETSc::DistributedCSR<HpSchur<PetscScalar>> Dbddc;
typedef PETSc::DistributedCSR<HpSchur<PetscReal>> DbddcR;
typedef PETSc::DistributedCSR<HpSchur<PetscComplex>> DbddcC;

namespace PETSc {
template<class Type>
struct _n_User;
template<class Type>
using User = _n_User<Type>*;
template<class HpddmType, typename std::enable_if<std::is_same<HpddmType, Dmat>::value>::type* = nullptr>
void initPETScStructure(HpddmType* ptA, PetscInt bs, PetscBool symmetric, KN<typename std::conditional<std::is_same<HpddmType, Dmat>::value, double, long>::type>* ptD, KN<PetscScalar>* rhs) {
    double timing = MPI_Wtime();
    PetscInt global;
    if(ptD) {
        ptA->_A->initialize(*ptD);
        unsigned int g;
        ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, g);
        global = g;
        if(verbosity > 0 && mpirank == 0)
            cout << " --- global numbering created (in " << MPI_Wtime() - timing << ")" << endl;
    }
    else
        global = PETSC_DECIDE;
    timing = MPI_Wtime();
    int* ia = nullptr;
    int* ja = nullptr;
    PetscScalar* c = nullptr;
    bool free = ptA->_A->getMatrix()->_ia ? ptA->_A->distributedCSR(ptA->_num, ptA->_first, ptA->_last, ia, ja, c) : false;
    MatCreate(PETSC_COMM_WORLD, &ptA->_petsc);
    if(bs > 1)
        MatSetBlockSize(ptA->_petsc, bs);
    MatSetSizes(ptA->_petsc, ptA->_last - ptA->_first, ptA->_last - ptA->_first, global, global);
    MatSetType(ptA->_petsc, MATMPIAIJ);
    if(ia)
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, reinterpret_cast<PetscInt*>(ia), reinterpret_cast<PetscInt*>(ja), c);
    else
        MatSetUp(ptA->_petsc);
    MatSetOption(ptA->_petsc, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
    MatSetOption(ptA->_petsc, MAT_SYMMETRIC, symmetric);
    if(free) {
        delete [] ia;
        delete [] ja;
        delete [] c;
    }
    ptA->_A->setBuffer();
    if(verbosity > 0 && mpirank == 0)
        cout << " --- global CSR created (in " << MPI_Wtime() - timing << ")" << endl;
    if(rhs)
        ptA->_A->exchange(*rhs);
}
template<class HpddmType, typename std::enable_if<!std::is_same<HpddmType, Dmat>::value>::type* = nullptr>
void initPETScStructure(HpddmType* ptA, PetscInt& bs, PetscBool symmetric, KN<typename std::conditional<std::is_same<HpddmType, Dmat>::value, double, long>::type>* ptD, KN<PetscScalar>* rhs) {
    const HPDDM::MatrixCSR<PetscScalar>* M = ptA->_A->getMatrix();
    if(!M->_sym)
        cout << "Please assemble a symmetric CSR" << endl;
    double timing = MPI_Wtime();
    ptA->_A->template renumber<false>(STL<long>(*ptD), nullptr);
    unsigned int global;
    ptA->_A->distributedNumbering(ptA->_num, ptA->_first, ptA->_last, global);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- global numbering created (in " << MPI_Wtime() - timing << ")" << endl;
    timing = MPI_Wtime();
    PetscInt* indices;
    PetscMalloc(sizeof(PetscInt) * M->_n / bs, &indices);
    for(unsigned int i = 0; i < M->_n; i += bs)
        indices[i / bs] = ptA->_num[i] / bs;
    ISLocalToGlobalMapping rmap;
    ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, bs, M->_n / bs, indices, PETSC_OWN_POINTER, &rmap);
    MatCreateIS(PETSC_COMM_WORLD, bs, PETSC_DECIDE, PETSC_DECIDE, global, global, rmap, NULL, &ptA->_petsc);
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
    PetscInt* ia = new PetscInt[M->_n / bs + 1];
    PetscInt* ja = new PetscInt[nnz / (bs * bs)];
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
    MatSetOption(ptA->_petsc, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
    MatSetOption(ptA->_petsc, MAT_SYMMETRIC, symmetric);
    MatAssemblyBegin(ptA->_petsc, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(ptA->_petsc, MAT_FINAL_ASSEMBLY);
    delete [] a;
    delete [] ja;
    delete [] ia;
    IS to, from;
    PetscInt nr;
    Vec rglobal;
    ISLocalToGlobalMappingGetSize(rmap, &nr);
    ISCreateStride(PETSC_COMM_SELF, nr, 0, 1, &to);
    ISLocalToGlobalMappingApplyIS(rmap, to, &from);
    MatCreateVecs(ptA->_petsc, &rglobal, NULL);
    Vec isVec;
    VecCreate(PETSC_COMM_SELF, &isVec);
    VecSetType(isVec, VECMPI);
    VecSetSizes(isVec, PETSC_DECIDE, nr);
    VecScatterCreate(rglobal, from, isVec, to, &ptA->_scatter);
    VecDestroy(&isVec);
    VecDestroy(&rglobal);
    ISDestroy(&from);
    ISDestroy(&to);
    // ISLocalToGlobalMappingView(rmap, PETSC_VIEWER_STDOUT_WORLD);
    ISLocalToGlobalMappingDestroy(&rmap);
    if(verbosity > 0 && mpirank == 0)
        cout << " --- global CSR created (in " << MPI_Wtime() - timing << ")" << endl;
}
template<class Type>
long globalNumbering(Type* const& A, KN<long>* const& numbering) {
    if(A) {
        numbering->resize(A->_A->getMatrix()->_n);
        if(A->_num)
            for(int i = 0; i < numbering->n; ++i)
                numbering->operator[](i) = A->_num[i];
    }
    return 0L;
}
template<class Type>
class changeOperator : public OneOperator {
    public:
        const int c;
        class changeOperator_Op : public E_F0mps {
            public:
                Expression A;
                Expression B;
                const int c;
                static const int n_name_param = 2;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                changeOperator_Op(const basicAC_F0& args, int d) : A(0), B(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    A = to<Type*>(args[0]);
                    if(c == 0)
                        B = to<Matrice_Creuse<PetscScalar>*>(args[1]);
                    else
                        B = to<Type*>(args[1]);
                }

                AnyType operator()(Stack stack) const;
                operator aType() const { return atype<long>(); }
        };
        E_F0* code(const basicAC_F0& args) const { return new changeOperator_Op(args, c); }
        changeOperator() : OneOperator(atype<long>(), atype<Type*>(), atype<Matrice_Creuse<PetscScalar>*>()), c(0) { }
        changeOperator(int) : OneOperator(atype<long>(), atype<Type*>(), atype<Type*>()), c(1) { }
};
template<class Type>
basicAC_F0::name_and_type changeOperator<Type>::changeOperator_Op::name_param[] = {
    {"restriction", &typeid(Matrice_Creuse<double>*)},
    {"nest", &typeid(Type*)}
};
template<class Type>
AnyType changeOperator<Type>::changeOperator_Op::operator()(Stack stack) const {
    Type* ptA = GetAny<Type*>((*A)(stack));
    if(c == 0) {
        Matrice_Creuse<PetscScalar>* mat = GetAny<Matrice_Creuse<PetscScalar>*>((*B)(stack));
        if(ptA && mat) {
            MatriceMorse<PetscScalar>* mN = nullptr;
            if(mat->A)
                mN = static_cast<MatriceMorse<PetscScalar>*>(&*(mat->A));
            PetscBool assembled;
            MatAssembled(ptA->_petsc, &assembled);
            if(mN) {
                Matrice_Creuse<double>* pList = nargs[0] ? GetAny<Matrice_Creuse<double>*>((*nargs[0])(stack)) : 0;
                HPDDM::MatrixCSR<void>* dL = nullptr;
                if(pList && pList->A) {
                    MatriceMorse<double>* mList = static_cast<MatriceMorse<double>*>(&*(pList->A));
#ifndef VERSION_MATRICE_CREUSE
                    // code faux .. F.H
                    // remove last term if it zero.
                    if(mList->n == mList->nbcoef - 1 && std::abs(mList->a[mList->nbcoef - 1]) < 1.0e-12) {
                        mList->lg[mList->n] -= 1;
                        mList->nbcoef -= 1;
                    }
                    ffassert(mList->n == mList->nbcoef);
                    ffassert(mList->m == mN->n);
#else
                    ffassert(mList->n == mList->nnz);
                    ffassert(mList->m == mN->n);
#endif
                    dL = new_HPDDM_MatrixCSRvoid(mList,false); //->n, mN->n, mList->n, mList->lg, mList->cl, false);
                }
                HPDDM::MatrixCSR<PetscScalar>* dM = new_HPDDM_MatrixCSR<PetscScalar>(mN);//->n, mN->m, mN->nbcoef, mN->a, mN->lg, mN->cl, mN->symetrique);
                HPDDM::MatrixCSR<PetscScalar>* dN;
                if(!dL)
                    dN = dM;
                else {
                    unsigned int* perm = new unsigned int[dM->_n]();
                    for(unsigned int i = 0; i < dL->_n; ++i)
                        perm[dL->_ja[i]] = i + 1;
                    dN = new HPDDM::MatrixCSR<PetscScalar>(dM, dL, perm);
                    delete [] perm;
                    delete dM;
                    dM = nullptr;
                }
                if(ptA->_A)
                    ptA->_A->setMatrix(dN);
                int* ia = nullptr;
                int* ja = nullptr;
                PetscScalar* c = nullptr;
                bool free = true;
                if(ptA->_cnum)
                    free = HPDDM::template Subdomain<PetscScalar>::distributedCSR(ptA->_num, ptA->_first, ptA->_last, ia, ja, c, dN, ptA->_num + dN->_n);
                else
                    free = ptA->_A->distributedCSR(ptA->_num, ptA->_first, ptA->_last, ia, ja, c);
                if(assembled) {
                    MatZeroEntries(ptA->_petsc);
                    for(PetscInt i = 0; i < ptA->_last - ptA->_first; ++i) {
                        PetscInt row = ptA->_first + i;
                        MatSetValues(ptA->_petsc, 1, &row, ia[i + 1] - ia[i], reinterpret_cast<PetscInt*>(ja + ia[i]), c + ia[i], INSERT_VALUES);
                    }
                }
                else
                    MatMPIAIJSetPreallocationCSR(ptA->_petsc, reinterpret_cast<PetscInt*>(ia), reinterpret_cast<PetscInt*>(ja), c);
                if(free) {
                    delete [] ia;
                    delete [] ja;
                    delete [] c;
                }
            }
            else if(!assembled)
                MatMPIAIJSetPreallocationCSR(ptA->_petsc, NULL, NULL, NULL);
            MatAssemblyBegin(ptA->_petsc, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(ptA->_petsc, MAT_FINAL_ASSEMBLY);
            Type* ptNest = nargs[1] ? GetAny<Type*>((*nargs[1])(stack)) : 0;
            if(ptNest) {
                PetscBool assembled;
                MatAssembled(ptNest->_petsc, &assembled);
                if(!assembled) {
                    PetscInt M, N;
                    Mat** mat;
                    MatNestGetSubMats(ptNest->_petsc, &M, &N, &mat);
                    PetscBool assemble = PETSC_TRUE;
                    for(PetscInt i = 0; i < M && assemble; ++i) {
                        for(PetscInt j = 0; j < N && assemble; ++j) {
                            if(mat[i][j]) {
                                PetscBool assembled;
                                MatAssembled(mat[i][j], &assembled);
                                if(!assembled)
                                    assemble = PETSC_FALSE;
                            }
                        }
                    }
                    if(assemble) {
                        MatAssemblyBegin(ptNest->_petsc, MAT_FINAL_ASSEMBLY);
                        MatAssemblyEnd(ptNest->_petsc, MAT_FINAL_ASSEMBLY);
                    }
                }
            }
            if(ptA->_ksp) {
                KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
                if(std::is_same<Type, Dmat>::value && !ptA->_S.empty()) {
                    KSPSetFromOptions(ptA->_ksp);
                    PetscBool assembled;
                    MatAssembled(ptA->_petsc, &assembled);
                    if(assembled) {
                        if(ptNest)
                            MatAssembled(ptNest->_petsc, &assembled);
                        if(assembled)
                            KSPSetUp(ptA->_ksp);
                    }
                    PC pc;
                    KSPGetPC(ptA->_ksp, &pc);
                    PCType type;
                    PCGetType(pc, &type);
                    PetscBool isFieldSplit;
                    PetscStrcmp(type, PCFIELDSPLIT, &isFieldSplit);
                    if(isFieldSplit)
                        PETSc::setCompositePC(pc, ptA->_S);
                }
            }
        }
    }
    else {
        Type* ptB = GetAny<Type*>((*B)(stack));
        if(ptB->_petsc) {
            MatType type;
            PetscBool isType;
            MatGetType(ptB->_petsc, &type);
            PetscStrcmp(type, MATNEST, &isType);
            if(!isType && ptA->_S.empty() && ptB->_S.empty()) {
                Mat backup = ptA->_petsc;
                ptA->_petsc = nullptr;
                ptA->dtor();
                ptA->_petsc = backup;
                PetscMPIInt flag;
                MPI_Comm_compare(PetscObjectComm((PetscObject)ptA->_petsc), PetscObjectComm((PetscObject)ptB->_petsc), &flag);
                if(flag != MPI_CONGRUENT && flag != MPI_IDENT) {
                    MatDestroy(&ptA->_petsc);
                    MatConvert(ptB->_petsc, MATSAME, MAT_INITIAL_MATRIX, &ptA->_petsc);
                    MatDestroy(&ptB->_petsc);
                }
                else {
                    MatHeaderReplace(ptA->_petsc, &ptB->_petsc);
                }
                if(ptB->_ksp)
                    KSPDestroy(&ptB->_ksp);
                ptA->_A = ptB->_A;
                ptB->_A = nullptr;
                ptA->_num = ptB->_num;
                ptA->_cnum = ptB->_cnum;
                ptA->_first = ptB->_first;
                ptA->_cfirst = ptB->_cfirst;
                ptA->_last = ptB->_last;
                ptA->_clast = ptB->_clast;
                ptB->_num = nullptr;
                ptA->_exchange = ptB->_exchange;
                if(ptB->_exchange) {
                    ptA->_exchange[1] = ptB->_exchange[1];
                    ptB->_exchange = nullptr;
                }
            }
        }
    }
    return 0L;
}
template<class Type, class K>
long originalNumbering(Type* const& A, KN<K>* const& in, KN<long>* const& interface) {
    if(A)
        A->_A->originalNumbering(STL<long>(*interface), *in);
    return 0L;
}
template<class K>
long renumber(KN<K>* const& in, KN<long>* const& numbering, KN<K>* const& out) {
    PetscInt low;
    out->resize(in->n);
    Vec x, y;
    VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, in->n, PETSC_DECIDE, static_cast<PetscScalar*>(*in), &x);
    VecGetOwnershipRange(x, &low, NULL);
    VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, in->n, PETSC_DECIDE, static_cast<PetscScalar*>(*out), &y);

    IS from, to;
    PetscInt* idx = new PetscInt[in->n];
    for(int i = 0; i < in->n; ++i)
        idx[i] = numbering->operator[](i);
    ISCreateGeneral(PETSC_COMM_SELF, in->n, idx, PETSC_USE_POINTER, &to);
    ISCreateStride(PETSC_COMM_SELF, in->n, low, 1, &from);
    VecScatter scatter;
    VecScatterCreate(x, NULL, y, to, &scatter);
    VecScatterBegin(scatter, x, y, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(scatter, x, y, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&scatter);
    ISDestroy(&to);
    ISDestroy(&from);
    delete [] idx;

    VecDestroy(&x);
    VecDestroy(&y);
    return 0L;
}
void finalizePETSc() {
    PETSC_COMM_WORLD = MPI_COMM_WORLD;
    PetscBool isFinalized;
    PetscFinalized(&isFinalized);
    if(!isFinalized)
        PetscFinalize();
}
template<class Type>
long initEmptyCSR(Type* const&) {
    return 0L;
}

template<class HpddmType>
class initCSRfromDMatrix_Op : public E_F0mps {
    public:
        Expression A;
        Expression B;
        Expression K;
        static const int n_name_param = 3;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        initCSRfromDMatrix_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : A(param1), B(param2), K(param3) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class HpddmType>
basicAC_F0::name_and_type initCSRfromDMatrix_Op<HpddmType>::name_param[] = {
    {"rhs", &typeid(KN<PetscScalar>*)},
    {"clean", &typeid(bool)},
    {"symmetric", &typeid(bool)},
};
template<class HpddmType>
class initCSRfromDMatrix : public OneOperator {
    public:
        initCSRfromDMatrix() : OneOperator(atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<Matrice_Creuse<PetscScalar>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new initCSRfromDMatrix_Op<HpddmType>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
        }
};
template<class HpddmType>
AnyType initCSRfromDMatrix_Op<HpddmType>::operator()(Stack stack) const {
    DistributedCSR<HpddmType>* ptA = GetAny<DistributedCSR<HpddmType>*>((*A)(stack));
    DistributedCSR<HpddmType>* ptB = GetAny<DistributedCSR<HpddmType>*>((*B)(stack));
    Matrice_Creuse<PetscScalar>* ptK = GetAny<Matrice_Creuse<PetscScalar>*>((*K)(stack));
    if(ptB->_A && ptK->A) {
        ptA->_A = new HpddmType(static_cast<const HPDDM::Subdomain<PetscScalar>&>(*ptB->_A));
        MatriceMorse<PetscScalar>* mA = static_cast<MatriceMorse<PetscScalar>*>(&(*ptK->A));
        HPDDM::MatrixCSR<PetscScalar>* dA = new_HPDDM_MatrixCSR<PetscScalar>(mA);//->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->_A->setMatrix(dA);
        ptA->_num = new unsigned int[mA->n];
        std::copy_n(ptB->_num, dA->_n, ptA->_num);
        ptA->_first = ptB->_first;
        ptA->_last = ptB->_last;
        PetscInt bs;
        MatGetBlockSize(ptB->_petsc, &bs);
        KN<PetscScalar>* rhs = nargs[0] ? GetAny<KN<PetscScalar>*>((*nargs[0])(stack)) : nullptr;
        initPETScStructure(ptA, bs, nargs[2] ? (GetAny<bool>((*nargs[2])(stack)) ? PETSC_TRUE : PETSC_FALSE) : PETSC_FALSE, static_cast<KN<double>*>(nullptr), rhs);
    }
    if(nargs[1] && GetAny<bool>((*nargs[1])(stack)))
        ptK->destroy();
    return ptA;
}

template<class HpddmType>
class initRectangularCSRfromDMatrix : public OneOperator {
    public:
        const int c;
        class initRectangularCSRfromDMatrix_Op : public E_F0mps {
            public:
                Expression A;
                Expression B;
                Expression C;
                Expression K;
                const int c;
                static const int n_name_param = 1;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                initRectangularCSRfromDMatrix_Op(const basicAC_F0& args, int d) : A(0), B(0), C(0), K(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    A = to<DistributedCSR<HpddmType>*>(args[0]);
                    B = to<DistributedCSR<HpddmType>*>(args[1]);
                    C = to<DistributedCSR<HpddmType>*>(args[2]);
                    if(c != 1)
                        K = to<Matrice_Creuse<PetscScalar>*>(args[3]);
                }

                AnyType operator()(Stack stack) const;
                operator aType() const { return atype<DistributedCSR<HpddmType>*>(); }
        };
        E_F0* code(const basicAC_F0& args) const { return new initRectangularCSRfromDMatrix_Op(args, c); }
        initRectangularCSRfromDMatrix() : OneOperator(atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<Matrice_Creuse<PetscScalar>*>()), c(0) { }
        initRectangularCSRfromDMatrix(int) : OneOperator(atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>()), c(1) { }
};
template<class HpddmType>
basicAC_F0::name_and_type initRectangularCSRfromDMatrix<HpddmType>::initRectangularCSRfromDMatrix_Op::name_param[] = {
    {"clean", &typeid(bool)}
};
template<class HpddmType>
AnyType initRectangularCSRfromDMatrix<HpddmType>::initRectangularCSRfromDMatrix_Op::operator()(Stack stack) const {
    DistributedCSR<HpddmType>* ptA = GetAny<DistributedCSR<HpddmType>*>((*A)(stack));
    DistributedCSR<HpddmType>* ptB = GetAny<DistributedCSR<HpddmType>*>((*B)(stack));
    DistributedCSR<HpddmType>* ptC = GetAny<DistributedCSR<HpddmType>*>((*C)(stack));
    Matrice_Creuse<PetscScalar>* ptK = (c != 1 ? GetAny<Matrice_Creuse<PetscScalar>*>((*K)(stack)) : nullptr);
    if(ptB->_A && ptC->_A) {
        ptA->_first = ptB->_first;
        ptA->_last = ptB->_last;
        ptA->_cfirst = ptC->_first;
        ptA->_clast = ptC->_last;
        PetscInt bsB, bsC;
        MatGetBlockSize(ptB->_petsc, &bsB);
        MatGetBlockSize(ptC->_petsc, &bsC);
        MatCreate(PETSC_COMM_WORLD, &ptA->_petsc);
        if(bsB == bsC && bsB > 1)
            MatSetBlockSize(ptA->_petsc, bsB);
        MatSetSizes(ptA->_petsc, ptB->_last - ptB->_first, ptC->_last - ptC->_first, PETSC_DECIDE, PETSC_DECIDE);
        MatSetType(ptA->_petsc, MATMPIAIJ);
        if(c != 1) {
            int* ia = nullptr;
            int* ja = nullptr;
            PetscScalar* a = nullptr;
            bool free = true;
            if(ptK->A) {
                MatriceMorse<PetscScalar>* mA = static_cast<MatriceMorse<PetscScalar>*>(&(*ptK->A));
                ff_HPDDM_MatrixCSR<PetscScalar> dA(mA);//->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
                ptA->_num = new unsigned int[mA->n + mA->m];
                ptA->_cnum = ptA->_num + mA->n;
                std::copy_n(ptB->_num, mA->n, ptA->_num);
                std::copy_n(ptC->_num, mA->m, ptA->_cnum);
                free = HPDDM::template Subdomain<PetscScalar>::distributedCSR(ptA->_num, ptA->_first, ptA->_last, ia, ja, a, &dA, ptA->_num + mA->n);
            }
            else
                ia = new int[ptB->_last - ptB->_first + 1]();
            MatMPIAIJSetPreallocationCSR(ptA->_petsc, reinterpret_cast<PetscInt*>(ia), reinterpret_cast<PetscInt*>(ja), a);
            if(free) {
                delete [] ia;
                delete [] ja;
                delete [] a;
            }
        }
        else {
            MatSetUp(ptA->_petsc);
            ptA->_num = new unsigned int[ptB->_A->getMatrix()->_n + ptC->_A->getMatrix()->_m];
            ptA->_cnum = ptA->_num + ptB->_A->getMatrix()->_n;
            std::copy_n(ptB->_num, ptB->_A->getMatrix()->_n, ptA->_num);
            std::copy_n(ptC->_num, ptC->_A->getMatrix()->_m, ptA->_cnum);
        }
        MatSetOption(ptA->_petsc, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
        ptA->_exchange = new HPDDM::template Subdomain<PetscScalar>*[2];
        ptA->_exchange[0] = new HPDDM::template Subdomain<PetscScalar>(*ptB->_A);
        ptA->_exchange[0]->setBuffer();
        ptA->_exchange[1] = new HPDDM::template Subdomain<PetscScalar>(*ptC->_A);
        ptA->_exchange[1]->setBuffer();
    }
    if(c != 1 && nargs[0] && GetAny<bool>((*nargs[0])(stack)))
        ptK->destroy();
    return ptA;
}

template<class HpddmType>
class initCSRfromMatrix_Op : public E_F0mps {
    public:
        Expression A;
        Expression K;
        Expression size;
        static const int n_name_param = 6;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        initCSRfromMatrix_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : A(param1), K(param2), size(param3) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class HpddmType>
basicAC_F0::name_and_type initCSRfromMatrix_Op<HpddmType>::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"bs", &typeid(long)},
    {"symmetric", &typeid(bool)},
    {"clean", &typeid(bool)},
    {"bsr", &typeid(bool)},
    {"prune", &typeid(bool)}
};
template<class HpddmType>
class initCSRfromMatrix : public OneOperator {
    public:
        initCSRfromMatrix() : OneOperator(atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<Matrice_Creuse<PetscScalar>*>(), atype<KN<long>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new initCSRfromMatrix_Op<HpddmType>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
        }
};
template<class HpddmType>
AnyType initCSRfromMatrix_Op<HpddmType>::operator()(Stack stack) const {
    if(nargs[0])
        PETSC_COMM_WORLD = *static_cast<MPI_Comm*>(GetAny<pcommworld>((*nargs[0])(stack)));
    DistributedCSR<HpddmType>* ptA = GetAny<DistributedCSR<HpddmType>*>((*A)(stack));
    Matrice_Creuse<PetscScalar>* ptK = GetAny<Matrice_Creuse<PetscScalar>*>((*K)(stack));
    KN<long>* ptSize = GetAny<KN<long>*>((*size)(stack));
    MatriceMorse<PetscScalar>* mK = static_cast<MatriceMorse<PetscScalar>*>(&(*(ptK->A)));
    PetscInt bs = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : 1;
    MatCreate(PETSC_COMM_WORLD, &ptA->_petsc);
    if(bs > 1)
        MatSetBlockSize(ptA->_petsc, bs);
    bool bsr = nargs[4] && GetAny<bool>((*nargs[4])(stack));
    bool prune = bsr ? (nargs[5] && GetAny<bool>((*nargs[5])(stack))) : false;
    if(prune) {
        ptA->_S.resize(1);
        MatCreate(PETSC_COMM_WORLD, &ptA->_S[0]);
    }
    if(mpisize > 1) {
        ffassert(ptSize->n >= 3 + mK->n);
        ptA->_first = ptSize->operator()(0);
        ptA->_last = ptSize->operator()(1);
        MatSetSizes(ptA->_petsc, mK->n * (bsr ? bs : 1), mK->n * (bsr ? bs : 1), ptSize->operator()(2) * (bsr ? bs : 1), ptSize->operator()(2) * (bsr ? bs : 1));
        if(prune)
            MatSetSizes(ptA->_S[0], mK->n * bs, mK->n * bs, ptSize->operator()(2) * bs, ptSize->operator()(2) * bs);
        ptA->_num = new unsigned int[ptSize->n - 3];
        for(int i = 3; i < ptSize->n; ++i)
            ptA->_num[i - 3] = ptSize->operator()(i);
    }
    else {
        ptA->_first = 0;
        ptA->_last = mK->n;
        ptA->_num = nullptr;
        MatSetSizes(ptA->_petsc, mK->n * (bsr ? bs : 1), mK->n * (bsr ? bs : 1), mK->n * (bsr ? bs : 1), mK->n * (bsr ? bs : 1));
        if(prune)
            MatSetSizes(ptA->_S[0], mK->n * bs, mK->n * bs, mK->n * bs, mK->n * bs);
    }
    bool clean = nargs[3] && GetAny<bool>((*nargs[3])(stack));
    if(clean)
        ptSize->resize(0);
    MatSetType(ptA->_petsc, bsr ? MATMPIBAIJ : MATMPIAIJ);
#ifndef VERSION_MATRICE_CREUSE
    if(bsr)
        MatMPIBAIJSetPreallocationCSR(ptA->_petsc, bs, reinterpret_cast<PetscInt*>(mK->lg), reinterpret_cast<PetscInt*>(mK->cl), mK->a);
    else
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, reinterpret_cast<PetscInt*>(mK->lg), reinterpret_cast<PetscInt*>(mK->cl), mK->a);
    if(prune) {
        MatSetType(ptA->_S[0], MATMPIAIJ);
        PetscInt* pI = new PetscInt[mK->n * bs + 1];
        PetscInt* pJ = new PetscInt[mK->nbcoef * bs];
        PetscScalar* pC = new PetscScalar[mK->nbcoef * bs];
        for(int i = 0; i < mK->n; ++i) {
            for(int k = 0; k < bs; ++k) {
                for(int j = mK->lg[i]; j < mK->lg[i + 1]; ++j) {
                    pJ[j - mK->lg[i] + bs * mK->lg[i] + k * (mK->lg[i + 1] - mK->lg[i])] = mK->cl[j] * bs + k;
                    pC[j - mK->lg[i] + bs * mK->lg[i] + k * (mK->lg[i + 1] - mK->lg[i])] = mK->a[j * bs * bs + k * bs + k];
                }
                pI[i * bs + k] = mK->lg[i] * bs + k * (mK->lg[i + 1] - mK->lg[i]);
            }
        }
        pI[mK->n * bs] = mK->nbcoef * bs;
#else
    mK->CSR();
    if(bsr)
        MatMPIBAIJSetPreallocationCSR(ptA->_petsc, bs, reinterpret_cast<PetscInt*>(mK->p), reinterpret_cast<PetscInt*>(mK->j), mK->aij);
    else
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, reinterpret_cast<PetscInt*>(mK->p), reinterpret_cast<PetscInt*>(mK->j), mK->aij);
    if(prune) {
        MatSetType(ptA->_S[0], MATMPIAIJ);
        PetscInt* pI = new PetscInt[mK->n * bs + 1];
        PetscInt* pJ = new PetscInt[mK->nnz * bs];
        PetscScalar* pC = new PetscScalar[mK->nnz * bs];
        for(int i = 0; i < mK->n; ++i) {
            for(int k = 0; k < bs; ++k) {
                for(int j = mK->p[i]; j < mK->p[i + 1]; ++j) {
                    pJ[j - mK->p[i] + bs * mK->p[i] + k * (mK->p[i + 1] - mK->p[i])] = mK->j[j] * bs + k;
                    pC[j - mK->p[i] + bs * mK->p[i] + k * (mK->p[i + 1] - mK->p[i])] = mK->aij[j * bs * bs + k * bs + k];
                }
                pI[i * bs + k] = mK->p[i] * bs + k * (mK->p[i + 1] - mK->p[i]);
            }
        }
        pI[mK->n * bs] = mK->nnz * bs;
#endif
        MatMPIAIJSetPreallocationCSR(ptA->_S[0], pI, pJ, pC);
        delete [] pC;
        delete [] pJ;
        delete [] pI;
    }
    MatSetOption(ptA->_petsc, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
    MatSetOption(ptA->_petsc, MAT_SYMMETRIC, nargs[2] ? (GetAny<bool>((*nargs[2])(stack)) ? PETSC_TRUE : PETSC_FALSE) : PETSC_FALSE);
    if(clean)
        ptK->destroy();
    if(prune) {
        KSPCreate(PETSC_COMM_WORLD, &ptA->_ksp);
        KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_S[0]);
    }
    return ptA;
}
template<class HpddmType>
class initCSRfromArray_Op : public E_F0mps {
    public:
        Expression A;
        Expression K;
        static const int n_name_param = 5;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        initCSRfromArray_Op(const basicAC_F0& args, Expression param1, Expression param2) : A(param1), K(param2) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class HpddmType>
basicAC_F0::name_and_type initCSRfromArray_Op<HpddmType>::name_param[] = {
    {"columns", &typeid(KN<long>*)},
    {"communicator", &typeid(pcommworld)},
    {"bs", &typeid(long)},
    {"symmetric", &typeid(bool)},
    {"clean", &typeid(bool)}
};
template<class HpddmType>
class initCSRfromArray : public OneOperator {
    public:
        initCSRfromArray() : OneOperator(atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<KN<Matrice_Creuse<PetscScalar>>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new initCSRfromArray_Op<HpddmType>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
        }
};
template<class HpddmType>
AnyType initCSRfromArray_Op<HpddmType>::operator()(Stack stack) const {
    if(nargs[1])
        PETSC_COMM_WORLD = *static_cast<MPI_Comm*>(GetAny<pcommworld>((*nargs[1])(stack)));
    int size;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    DistributedCSR<HpddmType>* ptA = GetAny<DistributedCSR<HpddmType>*>((*A)(stack));
    KN<Matrice_Creuse<PetscScalar>>* ptK = GetAny<KN<Matrice_Creuse<PetscScalar>>*>((*K)(stack));
    KN<long>* ptJ = nargs[0] ? GetAny<KN<long>*>((*nargs[0])(stack)) : nullptr;
    if(!ptJ || ptJ->n == ptK->n) {
        double timing = MPI_Wtime();
        std::vector<std::pair<int, int>> v;
        if(ptJ) {
            v.reserve(ptJ->n);
            for(int i = 0; i < ptJ->n; ++i)
                v.emplace_back(std::make_pair(ptJ->operator[](i), i));
            std::sort(v.begin(), v.end());
        }
        PetscInt bs = nargs[2] ? GetAny<long>((*nargs[2])(stack)) : 1;
        int rank;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
        int* dims = new int[size]();
        std::vector<std::pair<int, int>>::const_iterator it = std::lower_bound(v.cbegin(), v.cend(), std::make_pair(rank, 0), [](const std::pair<int, int>& lhs, const std::pair<int, int>& rhs) { return lhs.first < rhs.first; });
        int n = 0;
        if(!ptJ) {
            dims[rank] = ptK->operator[](mpisize / ptK->n > 1 ? 0 : rank).M();
            n = ptK->operator[](mpisize / ptK->n > 1 ? 0 : rank).N();
        }
        else if(it->first == rank) {
            dims[rank] = ptK->operator[](it->second).M();
            n = ptK->operator[](it->second).N();
        }
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, dims, 1, MPI_INT, PETSC_COMM_WORLD);
        if(ptJ)
            size = v.size();
        else
            size = ptK->n;
        PetscInt* ia = new PetscInt[n + 1];
        for(int i = 0; i < n + 1; ++i) {
            ia[i] = 0;
            for(int k = 0; k < size; ++k) {
                MatriceMorse<PetscScalar> *mA = static_cast<MatriceMorse<PetscScalar>*>(&(*(ptK->operator[](k)).A));
#ifndef VERSION_MATRICE_CREUSE
                if(i < mA->n + 1)
                    ia[i] += mA->lg[i];
#else
                mA->CSR();
                if(i < mA->n + 1)
                 ia[i] += mA->p[i];
#endif
            }
        }
        PetscInt* ja = new PetscInt[ia[n]];
        PetscScalar* c = new PetscScalar[ia[n]];
        int nnz = 0;
        int offset = (ptJ ? std::accumulate(dims, dims + v.begin()->first, 0) : 0);
 #ifndef VERSION_MATRICE_CREUSE
        for(int i = 0; i < n; ++i) {
            int backup = offset;
            for(int k = 0; k < size; ++k) {
                MatriceMorse<PetscScalar> *mA = static_cast<MatriceMorse<PetscScalar>*>(&(*(ptK->operator[](ptJ ? v[k].second : k)).A));
                if(i < mA->n) {
                    int j = mA->lg[i];
                    for( ; j < mA->lg[i + 1] && mA->cl[j] < dims[ptJ ? v[k].first : k]; ++j)
                        ja[nnz + j - mA->lg[i]] = mA->cl[j] + offset;
                    std::copy_n(mA->a + mA->lg[i], j - mA->lg[i], c + nnz);
                    nnz += j - mA->lg[i];
                    if(k < (ptJ ? v.size() : size) - 1)
                        offset += (ptJ ? std::accumulate(dims + v[k].first, dims + v[k + 1].first, 0) : dims[k]);
                }
            }
            offset = backup;
        }
#else
        for(int i = 0; i < n; ++i) {
            int backup = offset;
            for(int k = 0; k < size; ++k) {
                MatriceMorse<PetscScalar> *mA = static_cast<MatriceMorse<PetscScalar>*>(&(*(ptK->operator[](ptJ ? v[k].second : k)).A));
                mA->CSR();
                if(i < mA->n) {
                    int j = mA->p[i];
                    for( ; j < mA->p[i + 1] && mA->p[j] < dims[ptJ ? v[k].first : k]; ++j)
                        ja[nnz + j - mA->p[i]] = mA->j[j] + offset;
                    std::copy_n(mA->aij + mA->p[i], j - mA->p[i], c + nnz);
                    nnz += j - mA->p[i];
                    if(k < (ptJ ? v.size() : size) - 1)
                        offset += (ptJ ? std::accumulate(dims + v[k].first, dims + v[k + 1].first, 0) : dims[k]);
                }
            }
            offset = backup;
        }

#endif
        delete [] dims;
        MatCreate(PETSC_COMM_WORLD, &ptA->_petsc);
        if(bs > 1)
            MatSetBlockSize(ptA->_petsc, bs);
        if(!ptJ && (mpisize / ptK->n > 1))
            MatSetSizes(ptA->_petsc, n, n, PETSC_DECIDE, PETSC_DECIDE);
        else
            MatSetSizes(ptA->_petsc, n, ptK->operator[](ptJ ? it->second : rank).M(), PETSC_DECIDE, PETSC_DECIDE);
        ptA->_first = 0;
        ptA->_last = n;
        if(nargs[4] && GetAny<bool>((*nargs[4])(stack))) {
            int* cl = nullptr;
            int* lg = nullptr;
            PetscScalar* a = nullptr;
            for(int k = 0; k < size; ++k) {
                MatriceMorse<PetscScalar> *mA = static_cast<MatriceMorse<PetscScalar>*>(&(*(ptK->operator[](ptJ ? v[k].second : k)).A));
#ifndef VERSION_MATRICE_CREUSE
                if(k == 0) {
                    cl = mA->cl;
                    lg = mA->lg;
                    a = mA->a;
                }
                else {
                    if(mA->cl == cl)
                        mA->cl = nullptr;
                    if(mA->lg == lg)
                        mA->lg = nullptr;
                    if(mA->a == a)
                        mA->a = nullptr;
                }

#else
                if(k == 0) { 
                    mA->CSR();
                    cl = mA->j;
                    lg = mA->p;
                    a = mA->aij;
                }
                else {
                    if(mA->j == cl)
                        mA->j = nullptr;
                    if(mA->p == lg)
                        mA->p = nullptr;
                    if(mA->aij == a)
                        mA->aij = nullptr;
                }

#endif
            }
            ptK->resize(0);
        }
        MatSetType(ptA->_petsc, MATMPIAIJ);
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
        MatSetOption(ptA->_petsc, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
        MatSetOption(ptA->_petsc, MAT_SYMMETRIC, nargs[3] ? (GetAny<bool>((*nargs[3])(stack)) ? PETSC_TRUE : PETSC_FALSE) : PETSC_FALSE);
        delete [] ia;
        delete [] ja;
        delete [] c;
        if(verbosity > 0 && mpirank == 0)
            cout << " --- global CSR created (in " << MPI_Wtime() - timing << ")" << endl;
        timing = MPI_Wtime();
    }
    return ptA;
}

template<class HpddmType>
class initCSR : public OneOperator {
    public:
        const int c;
        class E_initCSR : public E_F0mps {
            public:
                Expression A;
                Expression K;
                Expression O;
                Expression R;
                Expression D;
                const int c;
                static const int n_name_param = 6;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                E_initCSR(const basicAC_F0& args, int d) : A(0), K(0), O(0), R(0), D(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    A = to<DistributedCSR<HpddmType>*>(args[0]);
                    if(c == 1)
                        K = to<long>(args[1]);
                    else
                        K = to<Matrice_Creuse<PetscScalar>*>(args[1]);
                    O = to<KN<long>*>(args[2]);
                    R = to<KN<KN<long>>*>(args[3]);
                    D = to<KN<typename std::conditional<std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value, double, long>::type>*>(args[4]);
                }

                AnyType operator()(Stack stack) const;
                operator aType() const { return atype<DistributedCSR<HpddmType>*>(); }
        };
        E_F0* code(const basicAC_F0 & args) const { return new E_initCSR(args, c); }
        initCSR() : OneOperator(atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<Matrice_Creuse<PetscScalar>*>(), atype<KN<long>*>(), atype<KN<KN<long>>*>(), atype<KN<typename std::conditional<std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value, double, long>::type>*>()), c(0) { }
        initCSR(int) : OneOperator(atype<DistributedCSR<HpddmType>*>(), atype<DistributedCSR<HpddmType>*>(), atype<long>(), atype<KN<long>*>(), atype<KN<KN<long>>*>(), atype<KN<typename std::conditional<std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value, double, long>::type>*>()), c(1) { }
};
template<class HpddmType>
basicAC_F0::name_and_type initCSR<HpddmType>::E_initCSR::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"bs", &typeid(long)},
    {"rhs", &typeid(KN<PetscScalar>*)},
    {"clean", &typeid(bool)},
    {"symmetric", &typeid(bool)},
    {"restriction", &typeid(Matrice_Creuse<double>*)}
};
template<class HpddmType>
AnyType initCSR<HpddmType>::E_initCSR::operator()(Stack stack) const {
    DistributedCSR<HpddmType>* ptA = GetAny<DistributedCSR<HpddmType>*>((*A)(stack));
    KN<KN<long>>* ptR = GetAny<KN<KN<long>>*>((*R)(stack));
    KN<long>* ptO = GetAny<KN<long>*>((*O)(stack));
    KN<typename std::conditional<std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value, double, long>::type>* ptD = GetAny<KN<typename std::conditional<std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value, double, long>::type>*>((*D)(stack));
    PetscInt bs = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : 1;
    int dof = 0;
    MatriceMorse<PetscScalar>* mA = nullptr;
    if(c == 0)
        mA = static_cast<MatriceMorse<PetscScalar>*>(&(*GetAny<Matrice_Creuse<PetscScalar>*>((*K)(stack))->A));
    else
        dof = GetAny<long>((*K)(stack));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    KN<PetscScalar>* rhs = nargs[2] ? GetAny<KN<PetscScalar>*>((*nargs[2])(stack)) : 0;
    ptA->_A = new HpddmType;
    if(comm)
        PETSC_COMM_WORLD = *comm;
    if(ptO && (mA || dof)) {
        HPDDM::MatrixCSR<PetscScalar>* dA;
        if(mA)
            dA = new_HPDDM_MatrixCSR<PetscScalar>(mA);//->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        else
            dA = new HPDDM::MatrixCSR<PetscScalar>(dof, dof, 0, nullptr, nullptr, nullptr, false);
        Matrice_Creuse<double>* pList = nargs[5] ? GetAny<Matrice_Creuse<double>*>((*nargs[5])(stack)) : 0;
        HPDDM::MatrixCSR<void>* dL = nullptr;
        if(std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value && pList && (mA || dof)) {
            int n = 0;
            ptA->_exchange = new HPDDM::template Subdomain<PetscScalar>*[2]();
            ptA->_exchange[0] = new HPDDM::template Subdomain<PetscScalar>();
            ptA->_exchange[0]->initialize(dA, STL<long>(*ptO), *ptR, comm);
            ptA->_exchange[0]->setBuffer();
            if(pList->A) {
                MatriceMorse<double>* mList = static_cast<MatriceMorse<double>*>(&*(pList->A));
#ifndef VERSION_MATRICE_CREUSE
                if(mList->n == mList->nbcoef - 1 && std::abs(mList->a[mList->nbcoef - 1]) < 1.0e-12) {
                    mList->lg[mList->n] -= 1;
                    mList->nbcoef -= 1;
                }
                ffassert(mList->n == mList->nbcoef);
#else
                mList->CSR();
                ffassert(mList->n == mList->nnz);
#endif
                ffassert(mList->m == (mA ? mA->n : dof));
                n = mList->n;
#ifndef VERSION_MATRICE_CREUSE
                dL = new HPDDM::MatrixCSR<void>(n, mA ? mA->n : dof, n, mList->lg, mList->cl, false);
                double* D = new double[n];
                for(int i = 0; i < n; ++i)
                    D[i] = ptD->operator[](mList->cl[i]);
#else
                dL = new HPDDM::MatrixCSR<void>(n, mA ? mA->n : dof, n, mList->p, mList->j, false);
                double* D = new double[n];
                for(int i = 0; i < n; ++i)
                    D[i] = ptD->operator[](mList->j[i]);
#endif
                ptD->resize(n);
                for(int i = 0; i < n; ++i)
                    ptD->operator[](i) = D[i];
                delete [] D;
            }
            else {
                dL = new HPDDM::MatrixCSR<void>(0, mA ? mA->n : dof, 0, nullptr, nullptr, false);
                ptD->destroy();
            }
        }
        ptA->_A->HPDDM::template Subdomain<PetscScalar>::initialize(dA, STL<long>(*ptO), *ptR, comm, dL);
        delete dL;
        ptA->_num = new unsigned int[ptA->_A->getMatrix()->_n];
        initPETScStructure(ptA, bs, nargs[4] ? (GetAny<bool>((*nargs[4])(stack)) ? PETSC_TRUE : PETSC_FALSE) : PETSC_FALSE, ptD, rhs);
        if(!std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value) {
#ifndef VERSION_MATRICE_CREUSE
            mA->lg = ptA->_A->getMatrix()->_ia;
#else
            mA->p = ptA->_A->getMatrix()->_ia;
#endif
        }
        if(nargs[3] && GetAny<bool>((*nargs[3])(stack))) {
            ptO->resize(0);
            ptR->resize(0);
            if(c == 0)
                GetAny<Matrice_Creuse<PetscScalar>*>((*K)(stack))->destroy();
        }
    }
    return ptA;
}
template<class HpddmType>
class initCSRfromBlockMatrix : public E_F0 {
    public:
        typedef DistributedCSR<HpddmType>* Result;
        PetscInt N, M;
        int init;
        Expression emat;
        Expression** e_M;
        int** t_M;
        initCSRfromBlockMatrix(const basicAC_F0& args, int i = 0) {
            init = i;
            M = 0;
            args.SetNameParam();
            emat = args[0];
            const E_Array& eM = *dynamic_cast<const E_Array*>((Expression)args[1]);
            N = eM.size();
            int err = 0;
            for(int i = 0; i < N; ++i) {
                const E_Array* e = dynamic_cast<const E_Array*>((Expression)eM[i]);
                if(!e)
                    ++err;
                else {
                    if(i == 0)
                        M = e->size();
                    else if(M != e->size())
                        ++err;
                }
            }
            if(err)
                CompileError("Wrong format of block matrix");
            assert(N && M);
            e_M = new Expression*[N];
            t_M = new int*[N];
            for(int i = 0; i < N; ++i) {
                const E_Array l = *dynamic_cast<const E_Array*>((Expression)eM[i]);
                e_M[i] = new Expression[M];
                t_M[i] = new int[M];
                for(int j = 0; j < M; ++j) {
                    C_F0 c_M(l[j]);
                    Expression e = c_M.LeftValue();
                    aType r = c_M.left();
                    if(r == atype<long>() && e->EvaluableWithOutStack()) {
                        long c = GetAny<long>((*e)(NullStack));
                        if(c == 0) {
                            e_M[i][j] = 0;
                            t_M[i][j] = 0;
                        }
                        else if(atype<PetscScalar>()->CastingFrom(r)) {
                            e_M[i][j] = to<PetscScalar>(c_M);
                            t_M[i][j] = 7;
                        }
                        else
                            CompileError("Empty block matrix");
                    }
                    else if(r == atype<DistributedCSR<HpddmType>*>()) {
                        e_M[i][j] = e;
                        t_M[i][j] = 1;
                    }
                    else if(r == atype<OpTrans<DistributedCSR<HpddmType>>>()) {
                        e_M[i][j] = e;
                        t_M[i][j] = 2;
                    }
                    else if(atype<KN_<PetscScalar>>()->CastingFrom(r)) {
                        e_M[i][j] = to<KN_<PetscScalar>>(c_M);
                        t_M[i][j] = 3;
                    }
                    else if(atype<Transpose<KN_<PetscScalar>>>()->CastingFrom(r)) {
                        e_M[i][j] = to<Transpose<KN_<PetscScalar>>>(c_M);
                        t_M[i][j] = 4;
                    }
                    else if(atype<PetscScalar>()->CastingFrom(r)) {
                        e_M[i][j] = to<PetscScalar>(c_M);
                        t_M[i][j] = 7;
                    }
                    else
                        CompileError("Unsupported type in submatrix");
                }
            }
        }
        ~initCSRfromBlockMatrix() {
            if(e_M) {
                for(int i = 0; i < N; ++i) {
                    delete [] e_M[i];
                    delete [] t_M[i];
                }
                delete [] e_M;
                delete [] t_M;
                N = 0;
                M = 0;
                e_M = nullptr;
                t_M = nullptr;
            }
        }
        static ArrayOfaType typeargs() { return ArrayOfaType(atype<Result>(), atype<E_Array>()); }
        static E_F0* f(const basicAC_F0& args) { return new initCSRfromBlockMatrix(args, 0); }
        AnyType operator()(Stack s) const {
            Mat* a = new Mat[N * M]();
            std::vector<int> zeros;
            zeros.reserve(std::min(N, M));
            for(int i = 0; i < N; ++i) {
                for(int j = 0; j < M; ++j) {
                    Expression e = e_M[i][j];
                    int t = t_M[i][j];
                    if(e) {
                        AnyType e_ij = (*e)(s);
                        if(t == 1) {
                            DistributedCSR<HpddmType>* pt = GetAny<DistributedCSR<HpddmType>*>(e_ij);
                            a[i * M + j] = pt->_petsc;
                        }
                        else if(t == 2) {
                            DistributedCSR<HpddmType>* pt = GetAny<DistributedCSR<HpddmType>*>(e_ij);
                            Mat B;
                            if(std::is_same<PetscScalar, PetscReal>::value)
                                MatCreateTranspose(pt->_petsc, &B);
                            else
                                MatCreateHermitianTranspose(pt->_petsc, &B);
                            a[i * M + j] = B;
                        }
                        else if(t == 7) {
                            PetscScalar r = GetAny<PetscScalar>(e_ij);
                            if(std::abs(r) > 1.0e-16)
                                ExecError("Nonzero scalar in submatrix");
                            if(i == j)
                                zeros.emplace_back(i);
                        }
                        else if(t == 3 || t == 4) {
                            KN<PetscScalar> x;
                            Mat B;
                            MatCreate(PETSC_COMM_WORLD, &B);
                            if(t == 4) {
                                x = GetAny<KN_<PetscScalar>>(e_ij);
                                MatSetSizes(B, PETSC_DECIDE, x.n, 1, PETSC_DECIDE);
                            }
                            else {
                                x = GetAny<Transpose<KN_<PetscScalar>>>(e_ij);
                                MatSetSizes(B, x.n, PETSC_DECIDE, PETSC_DECIDE, 1);
                            }
                            MatSetType(B, MATMPIDENSE);
                            MatMPIDenseSetPreallocation(B, PETSC_NULL);
                            PetscScalar* array;
                            MatDenseGetArray(B, &array);
                            if(array)
                                std::copy_n(static_cast<PetscScalar*>(x), x.n, array);
                            MatDenseRestoreArray(B, &array);
                            MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
                            MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
                            a[i * M + j] = B;
                        }
                        else {
                            ExecError("Unknown type in submatrix");
                        }
                    }
                    else if(i == j)
                        zeros.emplace_back(i);
                }
            }
            for(int i = 0; i < zeros.size(); ++i) {
                int posX = -1, posY = -1;
                for(int j = 0; j < M && posX == -1; ++j) {
                    if(j != zeros[i] && a[zeros[i] * M + j])
                        posX = j;
                }
                for(int j = 0; j < N && posY == -1; ++j) {
                    if(j != zeros[i] && a[zeros[i] + j * M])
                        posY = j;
                }
                if(posX == -1 && posY == -1)
                    ExecError("Zero row and zero column");
                else if(posX == -1 || posY == -1) {
                    PetscInt x, X, y, Y;
                    if(posX != -1) {
                        MatGetSize(a[zeros[i] * M + posX], &X, &Y);
                        MatGetLocalSize(a[zeros[i] * M + posX], &x, &y);
                    }
                    if(posY != -1) {
                        MatGetSize(a[zeros[i] + posY * M], posX == -1 ? &X : NULL, &Y);
                        MatGetLocalSize(a[zeros[i] + posY * M], posX == -1 ? &x : NULL, &y);
                        if(posX == -1) {
                            x = y;
                            X = Y;
                        }
                    }
                    else {
                        y = x;
                        Y = X;
                    }
                    MatCreate(PETSC_COMM_WORLD, a + zeros[i] * M + zeros[i]);
                    MatSetSizes(a[zeros[i] * M + zeros[i]], x, y, X, Y);
                    MatSetType(a[zeros[i] * M + zeros[i]], MATMPIAIJ);
                    MatMPIAIJSetPreallocation(a[zeros[i] * M + zeros[i]], 0, NULL, 0, NULL);
                    MatAssemblyBegin(a[zeros[i] * M + zeros[i]], MAT_FINAL_ASSEMBLY);
                    MatAssemblyEnd(a[zeros[i] * M + zeros[i]], MAT_FINAL_ASSEMBLY);
                }
            }
            Result sparse_mat = GetAny<Result>((*emat)(s));
            if(sparse_mat->_petsc)
                sparse_mat->~DistributedCSR();
            MatCreateNest(PETSC_COMM_WORLD, N, NULL, M, NULL, a, &sparse_mat->_petsc);
            delete [] a;
            return sparse_mat;
        }
};
template<class HpddmType>
class assignBlockMatrix : public initCSRfromBlockMatrix<HpddmType> {
    public:
        assignBlockMatrix(const basicAC_F0& args) : initCSRfromBlockMatrix<HpddmType>(args, 1) { }
        static E_F0* f(const basicAC_F0 & args) { return new initCSRfromBlockMatrix<HpddmType>(args, 1); }
};

template<class Type>
class setOptions_Op : public E_F0mps {
    public:
        Expression A;
        static const int n_name_param = 9;
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
    {"nearnullspace", &typeid(FEbaseArrayKn<PetscScalar>*)},
    {"fields", &typeid(KN<double>*)},
    {"names", &typeid(KN<String>*)},
    {"prefix", &typeid(std::string*)},
    {"schurPreconditioner", &typeid(KN<Matrice_Creuse<PetscScalar>>*)},
    {"schurList", &typeid(KN<double>*)},
    {"subksp", &typeid(Type*)},
    {"MatNullSpace", &typeid(KNM<PetscScalar>*)}
};
template<class Type, char> class LinearSolve;
template<class Type> class NonlinearSolve;
template<class Type>
struct _n_User {
    typename Type::MatF_O* op;
};
template<class Type>
struct _n_User<NonlinearSolve<Type>> {
    typename NonlinearSolve<Type>::VecF_O* op;
    typename LinearSolve<Type, 'N'>::MatF_O* r;
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
    Type* ptSub = nargs[7] ? GetAny<Type*>((*nargs[7])(stack)) : 0;
    if(ptSub && ptA->_ksp) {
        PC pc;
        KSPGetPC(ptA->_ksp, &pc);
        PCType type;
        PCGetType(pc, &type);
        PetscBool isFieldSplit;
        PetscStrcmp(type, PCFIELDSPLIT, &isFieldSplit);
        if(!isFieldSplit)
            return 0L;
        else {
            PetscInt nsplits;
            KSP* subksp;
            PCFieldSplitGetSubKSP(pc, &nsplits, &subksp);
            Mat** mat;
            PetscInt M, N;
            MatNestGetSubMats(ptA->_petsc, &M, &N, &mat);
            ffassert(M == N);
            for(int i = 0; i < M; ++i) {
                if(mat[i][i] == ptSub->_petsc) {
                    Mat A;
                    KSPGetOperators(subksp[i], &A, NULL);
                    KSPDestroy(&ptSub->_ksp);
                    ptSub->_ksp = subksp[i];
                    PetscObjectReference((PetscObject)subksp[i]);
                    KSPSetOperators(subksp[i], ptSub->_petsc, ptSub->_petsc);
                    ptA = ptSub;
                    break;
                }
            }
            PetscFree(subksp);
        }
    }
    std::string* options = nargs[0] ? GetAny<std::string*>((*nargs[0])(stack)) : NULL;
    bool isFieldSplit = insertOptions(options);
    if(ptA && ptA->_petsc) {
        PetscBool assembled;
        MatAssembled(ptA->_petsc, &assembled);
        if(!ptA->_ksp) {
            KSPCreate(PETSC_COMM_WORLD, &ptA->_ksp);
            KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
        }
        if(std::is_same<Type, Dmat>::value && isFieldSplit) {
            KN<double>* fields = nargs[2] ? GetAny<KN<double>*>((*nargs[2])(stack)) : 0;
            KN<String>* names = nargs[3] ? GetAny<KN<String>*>((*nargs[3])(stack)) : 0;
            KN<Matrice_Creuse<PetscScalar>>* mS = nargs[5] ? GetAny<KN<Matrice_Creuse<PetscScalar>>*>((*nargs[5])(stack)) : 0;
            KN<double>* pL = nargs[6] ? GetAny<KN<double>*>((*nargs[6])(stack)) : 0;
            setFieldSplitPC(ptA, ptA->_ksp, fields, names, mS, pL);
        }
        if(std::is_same<Type, Dmat>::value) {
            FEbaseArrayKn<PetscScalar>* ptNS = nargs[1] ? GetAny<FEbaseArrayKn<PetscScalar>*>((*nargs[1])(stack)) : 0;
            KNM<PetscScalar>* ptPETScNS = nargs[8] ? GetAny<KNM<PetscScalar>*>((*nargs[8])(stack)) : 0;
            PetscInt dim = ptNS ? ptNS->N : 0;
            int dimPETSc = ptPETScNS ? ptPETScNS->M() : 0;
            if(dim || dimPETSc) {
                Vec v;
                MatCreateVecs(ptA->_petsc, &v, NULL);
                Vec* ns;
                VecDuplicateVecs(v, std::max(dim, dimPETSc), &ns);
                VecDestroy(&v);
                if(std::max(dim, dimPETSc) == dimPETSc) {
                    PetscInt m;
                    VecGetLocalSize(ns[0], &m);
                    for(unsigned short i = 0; i < dimPETSc; ++i) {
                        PetscScalar* x;
                        VecGetArray(ns[i], &x);
                        for(int j = 0; j < m; ++j)
                            x[j] = (*ptPETScNS)(j, i);
                        VecRestoreArray(ns[i], &x);
                    }
                }
                else
                    for(unsigned short i = 0; i < dim; ++i) {
                        PetscScalar* x;
                        VecGetArray(ns[i], &x);
                        HPDDM::Subdomain<PetscScalar>::template distributedVec<0>(ptA->_num, ptA->_first, ptA->_last, static_cast<PetscScalar*>(*(ptNS->get(i))), x, ptNS->get(i)->n);
                        VecRestoreArray(ns[i], &x);
                    }
                PetscScalar* dots = new PetscScalar[std::max(dim, dimPETSc)];
                for(unsigned short i = 0; i < std::max(dim, dimPETSc); ++i) {
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
                MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, std::max(dim, dimPETSc), ns, &sp);
                if(!assembled) {
                    MatSetType(ptA->_petsc, MATMPIAIJ);
                    MatSetUp(ptA->_petsc);
                }
                MatSetNearNullSpace(ptA->_petsc, sp);
                MatNullSpaceDestroy(&sp);
                if(ns)
                    VecDestroyVecs(dim, &ns);
                delete [] ns;
            }
        }
        if(nargs[4])
            KSPSetOptionsPrefix(ptA->_ksp, GetAny<std::string*>((*nargs[4])(stack))->c_str());
        KSPSetFromOptions(ptA->_ksp);
        if(std::is_same<Type, Dmat>::value && nargs[2] && nargs[5] && nargs[6]) {
            if(assembled) {
                if(ptSub)
                    MatAssembled(GetAny<Type*>((*A)(stack))->_petsc, &assembled);
                if(assembled)
                    KSPSetUp(ptA->_ksp);
            }
            PC pc;
            KSPGetPC(ptA->_ksp, &pc);
            PETSc::setCompositePC(pc, ptA->_S);
        }
    }
    return 0L;
}

template<class Type>
class IterativeMethod_Op : public E_F0mps {
    public:
        Expression A;
        Expression rhs;
        Expression x;
        static const int n_name_param = 2;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        IterativeMethod_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : A(param1), rhs(param2), x(param3) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type>
basicAC_F0::name_and_type IterativeMethod_Op<Type>::name_param[] = {
    {"sparams", &typeid(std::string*)},
    {"prefix", &typeid(std::string*)}
};
template<class Type>
class IterativeMethod : public OneOperator {
    public:
        IterativeMethod() : OneOperator(atype<long>(), atype<Type*>(), atype<KN<PetscScalar>*>(), atype<KN<PetscScalar>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new IterativeMethod_Op<Type>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
        }
};
template<class Type>
AnyType IterativeMethod_Op<Type>::operator()(Stack stack) const {
    Type* ptA = GetAny<Type*>((*A)(stack));
    KN<PetscScalar>* ptRhs = GetAny<KN<PetscScalar>*>((*rhs)(stack));
    KN<PetscScalar>* ptX = GetAny<KN<PetscScalar>*>((*x)(stack));
    std::string* options = nargs[0] ? GetAny<std::string*>((*nargs[0])(stack)) : NULL;
    std::string* prefix = nargs[1] ? GetAny<std::string*>((*nargs[1])(stack)) : NULL;
    PetscInt bs;
    MatType type;
    MatGetType(ptA->_petsc, &type);
    PetscBool isNotBlock;
    PetscStrcmp(type, MATMPIAIJ, &isNotBlock);
    if(isNotBlock)
        bs = 1;
    else
        MatGetBlockSize(ptA->_petsc, &bs);
    if(!ptA->_ksp) {
        KSPCreate(PETSC_COMM_WORLD, &ptA->_ksp);
        KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
    }
    HPDDM::PETScOperator op(ptA->_ksp, ptA->_last - ptA->_first, bs);
    if(prefix)
        op.setPrefix(*prefix);
    HPDDM::Option& opt = *HPDDM::Option::get();
    if(options) {
        opt.parse(*options, mpirank == 0);
        if(mpirank != 0)
            opt.remove("verbosity");
    }
    Vec x, y;
    MatCreateVecs(ptA->_petsc, &x, &y);
    PetscScalar* ptr_x;
    VecGetArray(x, &ptr_x);
    HPDDM::Subdomain<PetscScalar>::template distributedVec<0>(ptA->_num, ptA->_first, ptA->_last, static_cast<PetscScalar*>(*ptRhs), ptr_x, ptRhs->n / bs, bs);
    PetscScalar* ptr_y;
    VecGetArray(y, &ptr_y);
    std::fill_n(ptr_y, op._n, 0.0);
    HPDDM::IterativeMethod::solve(op, ptr_x, ptr_y, 1, PETSC_COMM_WORLD);
    VecRestoreArray(x, &ptr_x);
    VecDestroy(&x);
    HPDDM::Subdomain<PetscScalar>::template distributedVec<1>(ptA->_num, ptA->_first, ptA->_last, static_cast<PetscScalar*>(*ptX), ptr_y, ptX->n / bs, bs);
    VecRestoreArray(y, &ptr_y);
    VecDestroy(&y);
    if(ptA->_A)
        ptA->_A->HPDDM::template Subdomain<PetscScalar>::exchange(static_cast<PetscScalar*>(*ptX));
    return 0L;
}

template<class Type>
class view_Op : public E_F0mps {
    public:
        Expression A;
        static const int n_name_param = 1;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        view_Op(const basicAC_F0& args, Expression param1) : A(param1) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type>
basicAC_F0::name_and_type view_Op<Type>::name_param[] = {
    {"object", &typeid(std::string*)}
};
template<class Type>
class view : public OneOperator {
    public:
        view() : OneOperator(atype<long>(), atype<Type*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new view_Op<Type>(args, t[0]->CastTo(args[0]));
        }
};
template<class Type>
AnyType view_Op<Type>::operator()(Stack stack) const {
    Type* ptA = GetAny<Type*>((*A)(stack));
    std::string* object = nargs[0] ? GetAny<std::string*>((*nargs[0])(stack)) : NULL;
    if(!object || object->compare("mat") == 0)
        MatView(ptA->_petsc, PETSC_VIEWER_STDOUT_WORLD);
    else {
        if(ptA->_ksp) {
            if(object->compare("ksp") == 0)
                KSPView(ptA->_ksp, PETSC_VIEWER_STDOUT_WORLD);
            else if(object->compare("pc") == 0) {
                PC pc;
                KSPGetPC(ptA->_ksp, &pc);
                PCView(pc, PETSC_VIEWER_STDOUT_WORLD);
            }
        }
    }
    return 0L;
}

template<class Type>
class changeNumbering : public OneOperator {
    public:
        const int c;
        class changeNumbering_Op : public E_F0mps {
            public:
                std::vector<std::pair<Expression, Expression>> E;
                Expression A;
                Expression in;
                Expression out;
                const int c;
                static const int n_name_param = 2;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                changeNumbering_Op(const basicAC_F0& args, int d) : A(0), in(0), out(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    if(c != 2) {
                        if(c == 1) {
                            A = to<KN<long>*>(args[0]);
                            in = to<KN<PetscScalar>*>(args[1]);
                        }
                        else {
                            E.reserve(1);
                            E.emplace_back(to<Type*>(args[0]), to<KN<PetscScalar>*>(args[1]));
                        }
                    }
                    else {
                        const E_Array* EA = dynamic_cast<const E_Array*>(args[0].LeftValue());
                        const E_Array* Ex = dynamic_cast<const E_Array*>(args[1].LeftValue());
                        ffassert(EA->size() == Ex->size() && EA->size());
                        E.reserve(EA->size());
                        for(int i = 0; i < EA->size(); ++i)
                            E.emplace_back(to<Type*>((*EA)[i]), to<KN<PetscScalar>*>((*Ex)[i]));
                    }
                    out = to<KN<PetscScalar>*>(args[2]);
                }

                AnyType operator()(Stack stack) const;
                operator aType() const { return atype<long>(); }
        };
        E_F0* code(const basicAC_F0 & args) const { return new changeNumbering_Op(args, c); }
        changeNumbering() : OneOperator(atype<long>(), atype<Type*>(), atype<KN<PetscScalar>*>(), atype<KN<PetscScalar>*>()), c(0) { }
        changeNumbering(int) : OneOperator(atype<long>(), atype<KN<long>*>(), atype<KN<PetscScalar>*>(), atype<KN<PetscScalar>*>()), c(1) { }
        changeNumbering(int, int) : OneOperator(atype<long>(), atype<E_Array>(), atype<E_Array>(), atype<KN<PetscScalar>*>()), c(2) { }
};
template<class Type>
basicAC_F0::name_and_type changeNumbering<Type>::changeNumbering_Op::name_param[] = {
    {"inverse", &typeid(bool)},
    {"exchange", &typeid(bool)}
};
template<class T, typename std::enable_if<std::is_same<T, KN<PetscScalar>>::value>::type* = nullptr>
void resize(T* v, int size) {
    v->resize(size);
}
template<class T, typename std::enable_if<!std::is_same<T, KN<PetscScalar>>::value>::type* = nullptr>
void resize(T* v, int size) { }
template<class T, class U>
void changeNumbering_func(unsigned int* const num, unsigned int first, unsigned int last, PetscInt m, PetscInt n, PetscInt bs, T* ptIn, U* ptOut, bool inverse) {
    PetscScalar* out;
    if(!inverse) {
        resize(ptOut, m ? m * bs : n);
        out = static_cast<PetscScalar*>(*ptOut);
        if(last - first)
            HPDDM::Subdomain<PetscScalar>::template distributedVec<0>(num, first, last, static_cast<PetscScalar*>(*ptIn), out, ptIn->n / bs, bs);
        else
            std::copy_n(static_cast<PetscScalar*>(*ptIn), ptIn->n, out);
    }
    else {
        resize(ptIn, n);
        *ptIn = PetscScalar();
        out = static_cast<PetscScalar*>(*ptOut);
        if(num)
            HPDDM::Subdomain<PetscScalar>::template distributedVec<1>(num, first, last, static_cast<PetscScalar*>(*ptIn), out, ptIn->n / bs, bs);
        else
            std::copy_n(out, ptIn->n, static_cast<PetscScalar*>(*ptIn));
    }
}
template<class Type>
AnyType changeNumbering<Type>::changeNumbering_Op::operator()(Stack stack) const {
    KN<PetscScalar>* ptOut = GetAny<KN<PetscScalar>*>((*out)(stack));
    bool inverse = nargs[0] && GetAny<bool>((*nargs[0])(stack));
    if(c == 0 || c == 2) {
        PetscScalar* pt = *ptOut;
        for(int j = 0; j < E.size(); ++j) {
            Type* ptA = GetAny<Type*>((*(E[j].first))(stack));
            if(ptA && (ptA->_last - ptA->_first || (inverse && ptA->_num))) {
                PetscInt bs;
                MatType type;
                MatGetType(ptA->_petsc, &type);
                PetscBool isNotBlock;
                PetscStrcmp(type, MATMPIAIJ, &isNotBlock);
                if(isNotBlock)
                    bs = 1;
                else
                    MatGetBlockSize(ptA->_petsc, &bs);
                PetscInt m;
                MatGetLocalSize(ptA->_petsc, &m, NULL);
                KN<PetscScalar>* ptIn = GetAny<KN<PetscScalar>*>((*(E[j].second))(stack));
                if(!inverse)
                    ffassert(ptIn->n == ptA->_A->getMatrix()->_n);
                if(c != 2) {
                    if(inverse)
                        ffassert(ptOut->n == bs * m);
                    changeNumbering_func(ptA->_num, ptA->_first, ptA->_last, m, ptA->_A->getMatrix()->_n, bs, ptIn, ptOut, inverse);
                }
                else {
                    KN_<PetscScalar> ptOutShift(pt, bs * m);
                    changeNumbering_func(ptA->_num, ptA->_first, ptA->_last, m, ptA->_A->getMatrix()->_n, bs, ptIn, &ptOutShift, inverse);
                }
                if(inverse && nargs[1]) {
                    bool exchange = GetAny<bool>((*nargs[1])(stack));
                    if(exchange)
                        ptA->_A->HPDDM::template Subdomain<PetscScalar>::exchange(*ptIn);
                }
                if(c == 2) {
                    pt += bs * m;
                    if(j == E.size() - 1)
                        ffassert(ptOut->n == std::distance(static_cast<PetscScalar*>(*ptOut), pt));
                }
            }
        }
    }
    else {
        KN<long>* ptA = GetAny<KN<long>*>((*A)(stack));
        unsigned int* num = reinterpret_cast<unsigned int*>(&((*ptA)[0]));
        KN<PetscScalar>* ptIn = GetAny<KN<PetscScalar>*>((*in)(stack));
        changeNumbering_func(num + 2, num[0], num[1], num[1] - num[0], ptA->n - 2, 1, ptIn, ptOut, inverse);
    }
    return 0L;
}

template<char N, class Type>
long MatMult(Type* const& A, KN<PetscScalar>* const& in, KN<PetscScalar>* const& out) {
    if(A) {
        Vec x, y;
        PetscInt size;
        if(N == 'T' || N == 'H') {
            MatCreateVecs(A->_petsc, &y, &x);
            VecGetLocalSize(x, &size);
            ffassert(in->n == size);
            VecGetLocalSize(y, &size);
            VecPlaceArray(x, *in);
            out->resize(size);
            VecPlaceArray(y, *out);
            if(N == 'T')
                MatMultTranspose(A->_petsc, x, y);
            else
                MatMultHermitianTranspose(A->_petsc, x, y);
        }
        else {
            MatCreateVecs(A->_petsc, &x, &y);
            VecGetLocalSize(x, &size);
            ffassert(in->n == size);
            VecGetLocalSize(y, &size);
            VecPlaceArray(x, *in);
            out->resize(size);
            VecPlaceArray(y, *out);
            MatMult(A->_petsc, x, y);
        }
        VecResetArray(y);
        VecResetArray(x);
        VecDestroy(&y);
        VecDestroy(&x);
    }
    return 0L;
}
template<char P, class Type>
long MatMatMult(Type* const& A, KNM<PetscScalar>* const& in, KNM<PetscScalar>* const& out) {
    if(A) {
        Mat C;
        MatType type;
        MatGetType(A->_petsc, &type);
        PetscBool isType;
        PetscStrcmp(type, MATMPIAIJ, &isType);
        if(mpisize == 1 && isType)
            MatMPIAIJGetLocalMat(A->_petsc, MAT_INITIAL_MATRIX, &C);
        else
            C = A->_petsc;
        Mat x, y;
        PetscInt n, m, bs, N, M;
        MatGetLocalSize(A->_petsc, &n, &m);
        MatGetSize(A->_petsc, &N, &M);
        if(isType)
            bs = 1;
        else {
            MatGetBlockSize(A->_petsc, &bs);
            N *= bs;
            M *= bs;
        }
        MatCreate(PETSC_COMM_WORLD, &x);
        MatSetSizes(x, P == 'T' || P == 'H' ? bs * n : bs * m, PETSC_DECIDE, P == 'T' || P == 'H' ? N : M, in->M());
        MatSetType(x, MATDENSE);
        if(mpisize == 1)
            MatSeqDenseSetPreallocation(x, &(in->operator()(0, 0)));
        else
            MatMPIDenseSetPreallocation(x, &(in->operator()(0, 0)));
        MatAssemblyBegin(x, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(x, MAT_FINAL_ASSEMBLY);
        if(P == 'T' || P == 'H') {
            ffassert(in->N() == bs * n);
            out->resize(bs * m, in->M());
            if(P == 'T')
                MatTransposeMatMult(C, x, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y);
            else {
                Mat w;
                MatDuplicate(x, MAT_COPY_VALUES, &w);
                MatConjugate(w);
                MatTransposeMatMult(C, w, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &y);
                MatDestroy(&w);
                MatConjugate(y);
            }
            PetscScalar* array;
            MatDenseGetArray(y, &array);
            if(array)
                std::copy_n(array, bs * m * in->M(), &(out->operator()(0, 0)));
            MatDenseRestoreArray(y, &array);
        }
        else {
            ffassert(in->N() == bs * m);
            out->resize(bs * n, in->M());
            MatCreate(PETSC_COMM_WORLD, &y);
            MatSetSizes(y, P == 'T' || P == 'H' ? bs * m : bs * n, PETSC_DECIDE, P == 'T' || P == 'H' ? M : N, in->M());
            MatSetType(y, MATDENSE);
            if(mpisize == 1)
                MatSeqDenseSetPreallocation(y, &(out->operator()(0, 0)));
            else
                MatMPIDenseSetPreallocation(y, &(out->operator()(0, 0)));
            MatAssemblyBegin(y, MAT_FINAL_ASSEMBLY);
            MatAssemblyEnd(y, MAT_FINAL_ASSEMBLY);
            MatMatMult(C, x, MAT_REUSE_MATRIX, PETSC_DEFAULT, &y);
        }
        MatDestroy(&y);
        MatDestroy(&x);
        if(C != A->_petsc)
            MatDestroy(&C);
    }
    return 0L;
}
template<class Type>
long MatConvert(Type* const& A, Type* const& B) {
    if(A->_petsc) {
        MatType type;
        PetscBool isType;
        MatGetType(A->_petsc, &type);
        PetscStrcmp(type, MATNEST, &isType);
        if(isType) {
            if(B->_petsc)
                B->~DistributedCSR();
            Mat** mat;
            PetscInt M, N;
            MatNestGetSubMats(A->_petsc, &M, &N, &mat);
            std::vector<std::pair<std::pair<PetscInt, PetscInt>, Mat>> b;
            b.reserve(M * N);
            for(PetscInt i = 0; i < M; ++i) {
                for(PetscInt j = 0; j < N; ++j) {
                    if(mat[i][j]) {
                        MatGetType(mat[i][j], &type);
                        PetscStrcmp(type, MATTRANSPOSEMAT, &isType);
                        if(isType) {
                            b.emplace_back(std::make_pair(std::make_pair(i, j), Mat()));
                            Mat B = mat[i][j];
                            Mat C;
                            if(std::is_same<PetscScalar, PetscReal>::value) {
                                MatTransposeGetMat(B, &b.back().second);
                                MatTranspose(b.back().second, MAT_INITIAL_MATRIX, &C);
                            }
                            else {
#if PETSC_VERSION_RELEASE
                                ffassert(0);
#else
                                MatHermitianTransposeGetMat(B, &b.back().second);
                                MatHermitianTranspose(b.back().second, MAT_INITIAL_MATRIX, &C);
#endif
                            }
                            MatDestroy(&B);
                            MatNestSetSubMat(A->_petsc, i, j, C);
                        }
                    }
                }
            }
            MatConvert(A->_petsc, MATMPIAIJ, MAT_INITIAL_MATRIX, &B->_petsc);
            for(std::pair<std::pair<PetscInt, PetscInt>, Mat> p : b) {
                Mat B = mat[p.first.first][p.first.second];
                MatDestroy(&B);
                if(std::is_same<PetscScalar, PetscReal>::value)
                    MatCreateTranspose(p.second, &B);
                else
                    MatCreateHermitianTranspose(p.second, &B);
                MatNestSetSubMat(A->_petsc, p.first.first, p.first.second, B);
            }
        }
    }
    return 0L;
}
template<class Type>
long MatZeroRows(Type* const& A, KN<double>* const& ptRows) {
    if(A->_petsc) {
        PetscInt bs;
        MatType type;
        MatGetType(A->_petsc, &type);
        PetscBool isNotBlock;
        PetscStrcmp(type, MATMPIAIJ, &isNotBlock);
        if(isNotBlock)
            bs = 1;
        else
            MatGetBlockSize(A->_petsc, &bs);
        PetscInt m;
        MatGetLocalSize(A->_petsc, &m, NULL);
        ffassert(ptRows->n == bs * m);
        std::vector<PetscInt> rows;
        rows.reserve(bs * m);
        PetscInt start;
        MatGetOwnershipRange(A->_petsc, &start, NULL);
        for(int i = 0; i < bs * m; ++i) {
            if(std::abs(ptRows->operator[](i)) > 1.0e-12)
                rows.emplace_back(start + i);
        }
        MatSetOption(A->_petsc, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
        MatSetOption(A->_petsc, MAT_NO_OFF_PROC_ZERO_ROWS, PETSC_TRUE);
        MatZeroRows(A->_petsc, rows.size(), rows.data(), 0.0, NULL, NULL);
    }
    return 0L;
}
template<class Type, char N = 'N'>
class LinearSolve : public OneOperator {
    public:
        typedef KN<PetscScalar> Kn;
        typedef KN_<PetscScalar> Kn_;
        class MatF_O : public RNM_VirtualMatrix<PetscScalar> {
            public:
                Stack stack;
                mutable Kn x;
                C_F0 c_x;
                Expression mat;
                typedef typename RNM_VirtualMatrix<PetscScalar>::plusAx plusAx;
                MatF_O(int n, Stack stk, const OneOperator* op) :
                    RNM_VirtualMatrix<PetscScalar>(n), stack(stk), x(n), c_x(CPValue(x)),
                    mat(op ? CastTo<Kn_>(C_F0(op->code(basicAC_F0_wa(c_x)), (aType)*op)) : 0) { }
                ~MatF_O() {
                    delete mat;
                    Expression zzz = c_x;
                    delete zzz;
                }
                void addMatMul(const Kn_& xx, Kn_& Ax) const {
                    ffassert(xx.N() == Ax.N());
                    x = xx;
                    Ax += GetAny<Kn_>((*mat)(stack));
                    WhereStackOfPtr2Free(stack)->clean();
                }
                plusAx operator*(const Kn& x) const { return plusAx(this, x); }
                bool ChecknbLine(int) const { return true; }
                bool ChecknbColumn(int) const { return true; }
        };
        const int c;
        class E_LinearSolve : public E_F0mps {
            public:
                Expression A;
                Expression x;
                Expression y;
                const OneOperator *codeA, *codeC;
                const int c;
                static const int n_name_param = 2;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                E_LinearSolve(const basicAC_F0& args, int d) : A(0), x(0), y(0), codeA(0), codeC(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    if(c != 2) {
                        if(c == 1) {
                            const Polymorphic* op = dynamic_cast<const Polymorphic*>(args[0].LeftValue());
                            ffassert(op);
                            codeA = op->Find("(", ArrayOfaType(atype<KN<PetscScalar>*>(), false));
                            if(nargs[0]) {
                                op = dynamic_cast<const Polymorphic*>(nargs[0]);
                                ffassert(op);
                                codeC = op->Find("(", ArrayOfaType(atype<Kn*>(), false));
                            }
                        }
                        else {
                            A = to<Type*>(args[0]);
                        }
                        x = to<KN<PetscScalar>*>(args[1]);
                        y = to<KN<PetscScalar>*>(args[2]);
                    }
                    else {
                        A = to<Type*>(args[0]);
                        x = to<KNM<PetscScalar>*>(args[1]);
                        y = to<KNM<PetscScalar>*>(args[2]);
                    }
                }

                AnyType operator()(Stack stack) const;
                operator aType() const { return atype<long>(); }
        };
        E_F0* code(const basicAC_F0 & args) const { return new E_LinearSolve(args, c); }
        LinearSolve() : OneOperator(atype<long>(), atype<Type*>(), atype<KN<PetscScalar>*>(), atype<KN<PetscScalar>*>()), c(0) { }
        LinearSolve(int) : OneOperator(atype<long>(), atype<Polymorphic*>(), atype<KN<PetscScalar>*>(), atype<KN<PetscScalar>*>()), c(1) { }
        LinearSolve(int, int) : OneOperator(atype<long>(), atype<Type*>(), atype<KNM<PetscScalar>*>(), atype<KNM<PetscScalar>*>()), c(2) { }
};
template<class Type, char N>
basicAC_F0::name_and_type LinearSolve<Type, N>::E_LinearSolve::name_param[] = {
    {"precon", &typeid(Polymorphic*)},
    {"sparams", &typeid(std::string*)}
};
template<class Type, class Container>
static PetscErrorCode Op_User(Container A, Vec x, Vec y);
template<class Type, char N>
AnyType LinearSolve<Type, N>::E_LinearSolve::operator()(Stack stack) const {
    if(c != 2) {
        KN<PetscScalar>* in = GetAny<KN<PetscScalar>*>((*x)(stack));
        KN<PetscScalar>* out = GetAny<KN<PetscScalar>*>((*y)(stack));
        if(A) {
            Type* ptA = GetAny<Type*>((*A)(stack));
            Vec x, y;
            MatCreateVecs(ptA->_petsc, &x, &y);
            PetscInt size;
            VecGetLocalSize(y, &size);
            ffassert(in->n == size);
            if(out->n != size) {
                out->resize(size);
                *out = PetscScalar();
            }
            VecPlaceArray(x, *in);
            VecPlaceArray(y, *out);
            if(!ptA->_ksp) {
                KSPCreate(PETSC_COMM_WORLD, &ptA->_ksp);
                KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
            }
            if(N == 'N')
                KSPSolve(ptA->_ksp, x, y);
            else if(N == 'H') {
                VecConjugate(x);
                KSPSolveTranspose(ptA->_ksp, x, y);
                VecConjugate(y);
            }
            VecResetArray(y);
            VecResetArray(x);
            VecDestroy(&y);
            VecDestroy(&x);
        }
        else {
            User<LinearSolve<Type, N>> user = nullptr;
            PetscNew(&user);
            user->op = new LinearSolve<Type, N>::MatF_O(in->n, stack, codeA);
            Mat S;
            MatCreateShell(PETSC_COMM_WORLD, in->n, in->n, PETSC_DECIDE, PETSC_DECIDE, user, &S);
            MatShellSetOperation(S, MATOP_MULT, (void (*)(void))Op_User<LinearSolve<Type, N>, Mat>);
            Vec x, y;
            MatCreateVecs(S, &x, &y);
            if(out->n != in->n) {
                out->resize(in->n);
                *out = PetscScalar();
            }
            VecPlaceArray(x, *in);
            VecPlaceArray(y, *out);
            KSP ksp;
            KSPCreate(PETSC_COMM_WORLD, &ksp);
            KSPSetOperators(ksp, S, S);
            std::string* options = nargs[1] ? GetAny<std::string*>((*nargs[1])(stack)) : NULL;
            insertOptions(options);
            PC pc;
            KSPGetPC(ksp, &pc);
            User<LinearSolve<Type, N>> userPC = nullptr;
            if(codeC) {
                PCSetType(pc, PCSHELL);
                PetscNew(&userPC);
                userPC->op = new LinearSolve<Type, N>::MatF_O(in->n, stack, codeC);
                PCShellSetContext(pc, userPC);
                PCShellSetApply(pc, Op_User<LinearSolve<Type, N>, PC>);
            }
            KSPSetFromOptions(ksp);
            KSPSolve(ksp, x, y);
            VecResetArray(y);
            VecResetArray(x);
            VecDestroy(&y);
            VecDestroy(&x);
            MatDestroy(&S);
            KSPDestroy(&ksp);
            if(codeC) {
                delete userPC->op;
                PetscFree(userPC);
            }
            delete user->op;
            PetscFree(user);
        }
    }
    else {
        KNM<PetscScalar>* in = GetAny<KNM<PetscScalar>*>((*x)(stack));
        KNM<PetscScalar>* out = GetAny<KNM<PetscScalar>*>((*y)(stack));
        if(A) {
            Type* ptA = GetAny<Type*>((*A)(stack));
            PetscInt m, bs;
            MatGetLocalSize(ptA->_petsc, &m, NULL);
            MatType type;
            MatGetType(ptA->_petsc, &type);
            PetscBool isNotBlock;
            PetscStrcmp(type, MATMPIAIJ, &isNotBlock);
            if(isNotBlock)
                bs = 1;
            else
                MatGetBlockSize(ptA->_petsc, &bs);
            ffassert(in->N() == bs * m);
            out->resize(bs * m, in->M());
            if(!ptA->_ksp) {
                KSPCreate(PETSC_COMM_WORLD, &ptA->_ksp);
                KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
            }
            HPDDM::PETScOperator op(ptA->_ksp, m, bs);
            op.apply(&in->operator()(0, 0), &out->operator()(0, 0), in->M());
        }
    }
    return 0L;
}
template<class Type>
class NonlinearSolve : public OneOperator {
    public:
        typedef KN<PetscScalar> Kn;
        typedef KN_<PetscScalar> Kn_;
        class VecF_O {
            public:
                Stack stack;
                mutable Kn x;
                C_F0 c_x;
                Expression mat;
                VecF_O(int n, Stack stk, const OneOperator* op) :
                    stack(stk), x(n), c_x(CPValue(x)),
                    mat(op ? CastTo<long>(C_F0(op->code(basicAC_F0_wa(c_x)), (aType)*op)) : 0) { }
                ~VecF_O() {
                    delete mat;
                    Expression zzz = c_x;
                    delete zzz;
                }
                long apply(const Kn_& xx) const {
                    x = xx;
                    long ret = GetAny<long>((*mat)(stack));
                    WhereStackOfPtr2Free(stack)->clean();
                    return ret;
                }
        };
        const int c;
        class E_NonlinearSolve : public E_F0mps {
            public:
                Expression A;
                Expression J;
                Expression r;
                Expression b;
                Expression x;
                const OneOperator *codeJ, *codeR;
                const int c;
                static const int n_name_param = 3;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                E_NonlinearSolve(const basicAC_F0& args, int d) : A(0), J(0), r(0), x(0), codeJ(0), codeR(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    A = to<Type*>(args[0]);
                    const Polymorphic* op = dynamic_cast<const Polymorphic*>(args[1].LeftValue());
                    ffassert(op);
                    codeJ = op->Find("(", ArrayOfaType(atype<KN<PetscScalar>*>(), false));
                    op = dynamic_cast<const Polymorphic*>(args[2].LeftValue());
                    ffassert(op);
                    codeR = op->Find("(", ArrayOfaType(atype<KN<PetscScalar>*>(), false));
                    if(c == 0)
                        x = to<KN<PetscScalar>*>(args[3]);
                    else {
                        b = to<KN<PetscScalar>*>(args[3]);
                        x = to<KN<PetscScalar>*>(args[4]);
                    }
                }

                AnyType operator()(Stack stack) const;
                operator aType() const { return atype<long>(); }
        };
        E_F0* code(const basicAC_F0 & args) const { return new E_NonlinearSolve(args, c); }
        NonlinearSolve() : OneOperator(atype<long>(), atype<Type*>(), atype<Polymorphic*>(), atype<Polymorphic*>(), atype<KN<PetscScalar>*>()), c(0) { }
        NonlinearSolve(int) : OneOperator(atype<long>(), atype<Type*>(), atype<Polymorphic*>(), atype<Polymorphic*>(), atype<KN<PetscScalar>*>(), atype<KN<PetscScalar>*>()), c(1) { }
};
template<class Type>
basicAC_F0::name_and_type NonlinearSolve<Type>::E_NonlinearSolve::name_param[] = {
    {"sparams", &typeid(std::string*)},
    {"xl", &typeid(KN<PetscScalar>*)},
    {"xu", &typeid(KN<PetscScalar>*)}
};
template<class Type>
PetscErrorCode FormJacobian(SNES snes, Vec x, Mat J, Mat B, void* ctx) {
    User<Type>*            user;
    const PetscScalar*       in;
    PetscScalar*            out;
    PetscErrorCode         ierr;

    PetscFunctionBegin;
    user = reinterpret_cast<User<Type>*>(ctx);
    typename NonlinearSolve<Type>::VecF_O* mat = reinterpret_cast<typename NonlinearSolve<Type>::VecF_O*>((*user)->op);
    VecGetArrayRead(x, &in);
    KN_<PetscScalar> xx(const_cast<PetscScalar*>(in), mat->x.n);
    long ret = mat->apply(xx);
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(ret);
}
template<class Type>
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void* ctx) {
    User<Type>*            user;
    const PetscScalar*       in;
    PetscScalar*            out;
    PetscErrorCode         ierr;

    PetscFunctionBegin;
    user = reinterpret_cast<User<Type>*>(ctx);
    typename LinearSolve<Type, 'N'>::MatF_O* mat = reinterpret_cast<typename LinearSolve<Type, 'N'>::MatF_O*>((*user)->r);
    VecGetArrayRead(x, &in);
    VecGetArray(f, &out);
    KN_<PetscScalar> xx(const_cast<PetscScalar*>(in), mat->N);
    KN_<PetscScalar> yy(out, mat->N);
    yy = *mat * xx;
    VecRestoreArray(f, &out);
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(0);
}
template<class Type>
AnyType NonlinearSolve<Type>::E_NonlinearSolve::operator()(Stack stack) const {
    Type* ptA = GetAny<Type*>((*A)(stack));
    if(ptA->_petsc) {
        KN<PetscScalar>* in = GetAny<KN<PetscScalar>*>((*x)(stack));
        PetscInt first = ptA->_first;
        PetscInt last = ptA->_last;
        if(!ptA->_num)
            MatGetOwnershipRange(ptA->_petsc, &first, &last);
        ffassert(in->n == 0 || in->n == last - first);
        if(in->n == 0)
            in->resize(last - first);
        User<NonlinearSolve<Type>> user = nullptr;
        PetscNew(&user);
        user->op = new NonlinearSolve<Type>::VecF_O(in->n, stack, codeJ);
        user->r = new typename LinearSolve<Type, 'N'>::MatF_O(in->n, stack, codeR);

        SNES snes;
        SNESCreate(PETSC_COMM_WORLD, &snes);
        Vec r, x, f;
        MatCreateVecs(ptA->_petsc, &r, NULL);
        SNESSetFunction(snes, r, FormFunction<NonlinearSolve<Type>>, &user);
        SNESSetJacobian(snes, ptA->_petsc, ptA->_petsc, FormJacobian<NonlinearSolve<Type>>, &user);
        KN<PetscScalar>* fl = nargs[1] ? GetAny<KN<PetscScalar>*>((*nargs[1])(stack)) : nullptr;
        bool setXl = (fl && fl->n == in->n);
        KN<PetscScalar>* fu = nargs[2] ? GetAny<KN<PetscScalar>*>((*nargs[2])(stack)) : nullptr;
        bool setXu = (fu && fu->n == in->n);
        Vec xu, xl = nullptr;
        if(setXl || setXu) {
            MatCreateVecs(ptA->_petsc, &xu, &xl);
            PetscScalar* ptr;
            if(setXl) {
                VecGetArray(xl, &ptr);
                std::copy_n(static_cast<PetscScalar*>(*fl), in->n, ptr);
                VecRestoreArray(xl, &ptr);
            }
            else
                VecSet(xl, PETSC_NINFINITY);
            if(setXu) {
                VecGetArray(xu, &ptr);
                std::copy_n(static_cast<PetscScalar*>(*fu), in->n, ptr);
                VecRestoreArray(xu, &ptr);
            }
            else
                VecSet(xu, PETSC_INFINITY);
            SNESVISetVariableBounds(snes, xl, xu);
        }
        std::string* options = nargs[0] ? GetAny<std::string*>((*nargs[0])(stack)) : NULL;
        insertOptions(options);
        SNESSetFromOptions(snes);
        if(ptA->_ksp)
            SNESSetKSP(snes, ptA->_ksp);
        {
            PetscInt n;
            VecGetSize(r, &n);
            VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, in->n, n, static_cast<PetscScalar*>(*in), &x);
            if(c == 1) {
                KN<PetscScalar>* fb = GetAny<KN<PetscScalar>*>((*b)(stack));
                ffassert(fb->n == last - first);
                VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, fb->n, n, static_cast<PetscScalar*>(*fb), &f);
            }
        }
        SNESSolve(snes, c == 1 ? f : NULL, x);
        VecDestroy(&x);
        if(c == 1)
            VecDestroy(&f);
        VecDestroy(&r);
        if(xl) {
            VecDestroy(&xl);
            VecDestroy(&xu);
        }
        SNESDestroy(&snes);

        delete user->r;
        delete user->op;
        PetscFree(user);
    }
    return 0L;
}
template<class Type, class Container, typename std::enable_if<std::is_same<Container, Mat>::value>::type* = nullptr>
static PetscErrorCode ContainerGetContext(Container A, Type& user) {
    return MatShellGetContext(A, &user);
}
template<class Type, class Container, typename std::enable_if<std::is_same<Container, PC>::value>::type* = nullptr>
static PetscErrorCode ContainerGetContext(Container A, Type& user) {
    return PCShellGetContext(A, (void**)&user);
}
template<class Type, class Container>
static PetscErrorCode Op_User(Container A, Vec x, Vec y) {
    User<Type>             user;
    const PetscScalar*       in;
    PetscScalar*            out;
    PetscErrorCode         ierr;

    PetscFunctionBegin;
    ierr = ContainerGetContext(A, user); CHKERRQ(ierr);
    typename LinearSolve<Type, 'N'>::MatF_O* mat = reinterpret_cast<typename LinearSolve<Type, 'N'>::MatF_O*>(user->op);
    VecGetArrayRead(x, &in);
    VecGetArray(y, &out);
    KN_<PetscScalar> xx(const_cast<PetscScalar*>(in), mat->N);
    KN_<PetscScalar> yy(out, mat->N);
    yy = *mat * xx;
    VecRestoreArray(y, &out);
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(0);
}

template<class Type>
class augmentation_Op : public E_F0mps {
    public:
        Expression A;
        Expression B;
        static const int n_name_param = 3;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        augmentation_Op(const basicAC_F0& args, Expression param1, Expression param2) : A(param1), B(param2) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type>
basicAC_F0::name_and_type augmentation_Op<Type>::name_param[] = {
    {"P", &typeid(Type*)},
    {"alpha", &typeid(PetscScalar)},
    {"pattern", &typeid(long)}
};
template<class Type>
class augmentation : public OneOperator {
    public:
        augmentation() : OneOperator(atype<long>(), atype<Type*>(), atype<Type*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new augmentation_Op<Type>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
        }
};
template<class Type>
AnyType augmentation_Op<Type>::operator()(Stack stack) const {
    Type* ptA = GetAny<Type*>((*A)(stack));
    Type* ptB = GetAny<Type*>((*B)(stack));
    Type* ptR = nargs[0] ? GetAny<Type*>((*nargs[0])(stack)) : nullptr;
    PetscScalar alpha = nargs[1] ? GetAny<PetscScalar>((*nargs[1])(stack)) : PetscScalar(1.0);
    long pattern = nargs[2] ? GetAny<long>((*nargs[2])(stack)) : 0;
    if(ptA->_petsc && ptB->_petsc) {
        if(!ptR || !ptR->_petsc) {
            MatAXPY(ptA->_petsc, alpha, ptB->_petsc, MatStructure(pattern));
        }
        else {
            Mat C;
            MatPtAP(ptB->_petsc, ptR->_petsc, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
            MatAXPY(ptA->_petsc, alpha, C, MatStructure(pattern));
            MatDestroy(&C);
        }
    }
    else if(ptA->_petsc && ptR && ptR->_petsc) {
        Mat C;
        MatTransposeMatMult(ptR->_petsc, ptR->_petsc, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
        MatAXPY(ptA->_petsc, alpha, C, MatStructure(pattern));
        MatDestroy(&C);
    }
    return 0L;
}

template<class T, class U, class K, char trans>
class InvPETSc {
    static_assert(std::is_same<K, PetscScalar>::value, "Wrong types");
    public:
        const T t;
        const U u;
        InvPETSc(T v, U w) : t(v), u(w) {}
        void solve(U out) const {
            if((*t)._petsc) {
                Vec x, y;
                double timing = MPI_Wtime();
                MatCreateVecs((*t)._petsc, &x, &y);
                PetscScalar* ptr;
                PetscInt bs;
                MatType type;
                MatGetType((*t)._petsc, &type);
                PetscBool isNotBlock;
                PetscStrcmp(type, MATMPIAIJ, &isNotBlock);
                if(isNotBlock)
                    bs = 1;
                else
                    MatGetBlockSize((*t)._petsc, &bs);
                if(std::is_same<typename std::remove_reference<decltype(*t.A->_A)>::type, HpSchwarz<PetscScalar>>::value) {
                    VecGetArray(x, &ptr);
                    HPDDM::Subdomain<K>::template distributedVec<0>((*t)._num, (*t)._first, (*t)._last, static_cast<PetscScalar*>(*u), ptr, u->n / bs, bs);
                    VecRestoreArray(x, &ptr);
                    if(t.A->_A)
                        std::fill_n(static_cast<PetscScalar*>(*out), out->n, 0.0);
                }
                else {
                    VecSet(x, PetscScalar());
#if 0
                    Mat_IS* is = (Mat_IS*)(*t)._petsc->data;
                    VecGetArray(is->y, &ptr);
                    std::copy_n(static_cast<PetscScalar*>(*u), (*t)._A->getMatrix()->_n, ptr);
                    VecRestoreArray(is->y, &ptr);
                    VecScatterBegin(is->rctx,is->y,x,ADD_VALUES,SCATTER_REVERSE);
                    VecScatterEnd(is->rctx,is->y,x,ADD_VALUES,SCATTER_REVERSE);
#else
                    Vec isVec;
                    VecCreateMPIWithArray(PETSC_COMM_SELF, bs, (*t)._A->getMatrix()->_n, (*t)._A->getMatrix()->_n, static_cast<PetscScalar*>(*u), &isVec);
                    VecScatterBegin((*t)._scatter, isVec, x, ADD_VALUES, SCATTER_REVERSE);
                    VecScatterEnd((*t)._scatter, isVec, x, ADD_VALUES, SCATTER_REVERSE);
                    VecDestroy(&isVec);
#endif
                }
                timing = MPI_Wtime();
                if(!(*t)._ksp) {
                    KSPCreate(PETSC_COMM_WORLD, &(*t)._ksp);
                    KSPSetOperators((*t)._ksp, (*t)._petsc, (*t)._petsc);
                }
                if(trans == 'N')
                    KSPSolve((*t)._ksp, x, y);
                else {
                    if(!std::is_same<PetscScalar, PetscReal>::value && t.conjugate)
                        VecConjugate(x);
                    KSPSolveTranspose((*t)._ksp, x, y);
                    if(!std::is_same<PetscScalar, PetscReal>::value && t.conjugate)
                        VecConjugate(y);
                }
                if(verbosity > 0 && mpirank == 0)
                    cout << " --- system solved with PETSc (in " << MPI_Wtime() - timing << ")" << endl;
                if(std::is_same<typename std::remove_reference<decltype(*t.A->_A)>::type, HpSchwarz<PetscScalar>>::value) {
                    VecGetArray(y, &ptr);
                    HPDDM::Subdomain<K>::template distributedVec<1>((*t)._num, (*t)._first, (*t)._last, static_cast<PetscScalar*>(*out), ptr, out->n / bs, bs);
                    VecRestoreArray(y, &ptr);
                }
                else {
#if 0
                    Mat_IS* is = (Mat_IS*)(*t)._petsc->data;
                    VecScatterBegin(is->rctx,y,is->y,INSERT_VALUES,SCATTER_FORWARD);
                    VecScatterEnd(is->rctx,y,is->y,INSERT_VALUES,SCATTER_FORWARD);
                    VecGetArray(is->y, &ptr);
                    std::copy_n(ptr, (*t)._A->getMatrix()->_n, (PetscScalar*)*out);
                    VecRestoreArray(is->y, &ptr);
#else
                    Vec isVec;
                    VecCreateMPIWithArray(PETSC_COMM_SELF, bs, (*t)._A->getMatrix()->_n, (*t)._A->getMatrix()->_n, static_cast<PetscScalar*>(*out), &isVec);
                    VecScatterBegin((*t)._scatter, y, isVec, INSERT_VALUES, SCATTER_FORWARD);
                    VecScatterEnd((*t)._scatter, y, isVec, INSERT_VALUES, SCATTER_FORWARD);
                    VecDestroy(&isVec);
#endif
                }
                VecDestroy(&x);
                VecDestroy(&y);
                if(std::is_same<typename std::remove_reference<decltype(*t.A->_A)>::type, HpSchwarz<PetscScalar>>::value && t.A->_A)
                    (*t)._A->HPDDM::template Subdomain<PetscScalar>::exchange(*out);
            }
        };
        static U inv(U Ax, InvPETSc<T, U, K, trans> A) {
            A.solve(Ax);
            return Ax;
        }
        static U init(U Ax, InvPETSc<T, U, K, trans> A) {
            PetscInt n, m;
            MatGetSize(A.t.A->_petsc, &n, &m);
            ffassert(n == m);
            Ax->init(A.u->n);
            return inv(Ax, A);
        }
};

template<class T, class U, class K, char N>
class ProdPETSc {
    static_assert(std::is_same<K, PetscScalar>::value, "Wrong types");
    public:
        const T t;
        const U u;
        ProdPETSc(T v, U w) : t(v), u(w) {}
        void prod(U out) const {
            if((*t)._petsc) {
                PetscBool assembled;
                MatAssembled((*t)._petsc, &assembled);
                if(assembled) {
                    Vec x, y;
                    PetscScalar* ptr;
                    if(N == 'N')
                        MatCreateVecs((*t)._petsc, &x, &y);
                    else
                        MatCreateVecs((*t)._petsc, &y, &x);
                    VecGetArray(x, &ptr);
                    PetscInt bs;
                    MatType type;
                    MatGetType((*t)._petsc, &type);
                    PetscBool isNotBlock;
                    PetscStrcmp(type, MATMPIAIJ, &isNotBlock);
                    if(isNotBlock)
                        bs = 1;
                    else
                        MatGetBlockSize((*t)._petsc, &bs);
                    if(!t->_cnum || N == 'T') {
                        if(t->_num)
                            HPDDM::Subdomain<K>::template distributedVec<0>(t->_num, t->_first, t->_last, static_cast<PetscScalar*>(*u), ptr, u->n / bs, bs);
                    }
                    else
                        HPDDM::Subdomain<K>::template distributedVec<0>(t->_cnum, t->_cfirst, t->_clast, static_cast<PetscScalar*>(*u), ptr, u->n / bs, bs);
                    VecRestoreArray(x, &ptr);
                    if(N == 'T')
                        MatMultTranspose(t->_petsc, x, y);
                    else
                        MatMult(t->_petsc, x, y);
                    VecDestroy(&x);
                    VecGetArray(y, &ptr);
                    if(!t->_A)
                        std::fill_n(static_cast<PetscScalar*>(*out), out->n, 0.0);
                    if(!t->_cnum || N == 'N') {
                        if(t->_num)
                            HPDDM::Subdomain<K>::template distributedVec<1>(t->_num, t->_first, t->_last, static_cast<PetscScalar*>(*out), ptr, out->n / bs, bs);
                    }
                    else
                        HPDDM::Subdomain<K>::template distributedVec<1>(t->_cnum, t->_cfirst, t->_clast, static_cast<PetscScalar*>(*out), ptr, out->n / bs, bs);
                    if(t->_A)
                        (*t)._A->HPDDM::template Subdomain<PetscScalar>::exchange(*out);
                    else if(t->_exchange) {
                        if(N == 'N')
                            (*t)._exchange[0]->exchange(*out);
                        else
                            (*t)._exchange[1]->exchange(*out);
                    }
                    VecRestoreArray(y, &ptr);
                    VecDestroy(&y);
                }
            }
        }
        static U mv(U Ax, ProdPETSc<T, U, K, N> A) {
            *Ax = K();
            A.prod(Ax);
            return Ax;
        }
        static U init(U Ax, ProdPETSc<T, U, K, N> A) {
            PetscInt n, m;
            MatGetSize(A.t->_petsc, &n, &m);
            ffassert(n == m);
            Ax->init(A.u->n);
            return mv(Ax, A);
        }
};
}

template<typename Type>
bool CheckPetscMatrix(Type* ptA) {
    ffassert(ptA);
    PetscValidHeaderSpecific(ptA->_petsc, MAT_CLASSID, 2);
    return true;
}

static void Init_PETSc() {
    if(verbosity > 1 && mpirank == 0)
        cout << " PETSc (" << typeid(PetscScalar).name() << ")" << endl;
    int argc = pkarg->n;
    char** argv = new char*[argc];
    for(int i = 0; i < argc; ++i)
        argv[i] = const_cast<char*>((*(*pkarg)[i].getap())->c_str());
    PetscInitialize(&argc, &argv, 0, "");
    PetscSysInitializePackage();
    MatInitializePackage();
    if(argc > 1) {
        HPDDM::Option& opt = *HPDDM::Option::get();
        opt.parse(argc - 1, argv + 1, mpirank == 0);
        if(mpirank != 0)
            opt.remove("verbosity");
    }
    delete [] argv;
    KSPRegister("hpddm", HPDDM::KSPCreate_HPDDM);
    ff_atend(PETSc::finalizePETSc);
    if(!exist_type<DmatR*>()) {
        Dcl_Type<DmatR*>(Initialize<DmatR>, Delete<DmatR>);
        zzzfff->Add("Mat", atype<DmatR*>());
    }
    if(!exist_type<DmatC*>())
        Dcl_Type<DmatC*>(Initialize<DmatC>, Delete<DmatC>);
    if(!exist_type<DbddcR*>()) {
        Dcl_Type<DbddcR*>(Initialize<DbddcR>, Delete<DbddcR>);
        zzzfff->Add("MatIS", atype<DbddcR*>());
    }
    if(!exist_type<DbddcC*>())
        Dcl_Type<DbddcC*>(Initialize<DbddcC>, Delete<DbddcC>);
    map_type_of_map[make_pair(atype<DmatR*>(), atype<Complex*>())] = atype<DmatC*>();
    map_type_of_map[make_pair(atype<DmatR*>(), atype<double*>())] = atype<DmatR*>();
    map_type_of_map[make_pair(atype<DbddcR*>(), atype<Complex*>())] = atype<DbddcC*>();
    map_type_of_map[make_pair(atype<DbddcR*>(), atype<double*>())] = atype<DbddcR*>();

    TheOperators->Add("<-", new OneOperator1_<long, Dmat*>(PETSc::initEmptyCSR<Dmat>));
    if(std::is_same<PetscInt, int>::value) {
        TheOperators->Add("<-", new PETSc::initCSR<HpSchwarz<PetscScalar>>);
        TheOperators->Add("<-", new PETSc::initCSR<HpSchwarz<PetscScalar>>(1));
        TheOperators->Add("<-", new PETSc::initCSRfromArray<HpSchwarz<PetscScalar>>);
        TheOperators->Add("<-", new PETSc::initCSRfromMatrix<HpSchwarz<PetscScalar>>);
        TheOperators->Add("<-", new PETSc::initCSRfromDMatrix<HpSchwarz<PetscScalar>>);
        TheOperators->Add("<-", new PETSc::initRectangularCSRfromDMatrix<HpSchwarz<PetscScalar>>);
        TheOperators->Add("<-", new PETSc::initRectangularCSRfromDMatrix<HpSchwarz<PetscScalar>>(1));
        TheOperators->Add("<-", new OneOperatorCode<PETSc::initCSRfromBlockMatrix<HpSchwarz<PetscScalar>>>());
        TheOperators->Add("=", new OneOperatorCode<PETSc::assignBlockMatrix<HpSchwarz<PetscScalar>>>());
    }
    Global.Add("set", "(", new PETSc::setOptions<Dmat>());
    addProd<Dmat, PETSc::ProdPETSc, KN<PetscScalar>, PetscScalar, 'N'>();
    addProd<Dmat, PETSc::ProdPETSc, KN<PetscScalar>, PetscScalar, 'T'>();
    addInv<Dmat, PETSc::InvPETSc, KN<PetscScalar>, PetscScalar, 'N'>();
    addInv<Dmat, PETSc::InvPETSc, KN<PetscScalar>, PetscScalar, 'T'>();

    TheOperators->Add("<-", new OneOperator1_<long, Dbddc*>(PETSc::initEmptyCSR<Dbddc>));
    TheOperators->Add("<-", new PETSc::initCSR<HpSchur<PetscScalar>>);

    Global.Add("exchange", "(", new exchangeIn<Dmat, PetscScalar>);
    Global.Add("exchange", "(", new exchangeInOut<Dmat, PetscScalar>);
    Global.Add("changeNumbering", "(", new PETSc::changeNumbering<Dmat>());
    Global.Add("changeNumbering", "(", new PETSc::changeNumbering<Dmat>(1));
    Global.Add("changeNumbering", "(", new PETSc::changeNumbering<Dmat>(1, 2));
    Global.Add("MatMult", "(", new OneOperator3_<long, Dmat*, KN<PetscScalar>*, KN<PetscScalar>*>(PETSc::MatMult<'N'>));
    Global.Add("MatMatMult", "(", new OneOperator3_<long, Dmat*, KNM<PetscScalar>*, KNM<PetscScalar>*>(PETSc::MatMatMult<'N'>));
    Global.Add("MatMultTranspose", "(", new OneOperator3_<long, Dmat*, KN<PetscScalar>*, KN<PetscScalar>*>(PETSc::MatMult<'T'>));
    Global.Add("MatTransposeMatMult", "(", new OneOperator3_<long, Dmat*, KNM<PetscScalar>*, KNM<PetscScalar>*>(PETSc::MatMatMult<'T'>));
    if(!std::is_same<PetscScalar, PetscReal>::value) {
        Global.Add("MatMultHermitianTranspose", "(", new OneOperator3_<long, Dmat*, KN<PetscScalar>*, KN<PetscScalar>*>(PETSc::MatMult<'H'>));
        Global.Add("MatHermitianTransposeMatMult", "(", new OneOperator3_<long, Dmat*, KNM<PetscScalar>*, KNM<PetscScalar>*>(PETSc::MatMatMult<'H'>));
    }
    Global.Add("MatConvert", "(", new OneOperator2_<long, Dmat*, Dmat*>(PETSc::MatConvert));
    Global.Add("MatZeroRows", "(", new OneOperator2_<long, Dmat*, KN<double>*>(PETSc::MatZeroRows));
    Global.Add("KSPSolve", "(", new PETSc::LinearSolve<Dmat>());
    Global.Add("KSPSolve", "(", new PETSc::LinearSolve<Dmat>(1));
    Global.Add("KSPSolve", "(", new PETSc::LinearSolve<Dmat>(1, 1));
    if(!std::is_same<PetscScalar, PetscReal>::value)
        Global.Add("KSPSolveHermitianTranspose", "(", new PETSc::LinearSolve<Dmat, 'H'>());
    Global.Add("SNESSolve", "(", new PETSc::NonlinearSolve<Dmat>());
    Global.Add("SNESSolve", "(", new PETSc::NonlinearSolve<Dmat>(1));
    Global.Add("augmentation", "(", new PETSc::augmentation<Dmat>);
    Global.Add("globalNumbering", "(", new OneOperator2_<long, Dmat*, KN<long>*>(PETSc::globalNumbering<Dmat>));
    Global.Add("globalNumbering", "(", new OneOperator2_<long, Dbddc*, KN<long>*>(PETSc::globalNumbering<Dbddc>));
    Global.Add("changeOperator", "(", new PETSc::changeOperator<Dmat>());
    Global.Add("changeOperator", "(", new PETSc::changeOperator<Dmat>(1));
    Global.Add("IterativeMethod", "(", new PETSc::IterativeMethod<Dmat>);
    Global.Add("view", "(", new PETSc::view<Dmat>);
    Global.Add("originalNumbering", "(", new OneOperator3_<long, Dbddc*, KN<PetscScalar>*, KN<long>*>(PETSc::originalNumbering));
    Global.Add("renumber", "(", new OneOperator3_<long, KN<PetscScalar>*, KN<long>*, KN<PetscScalar>*>(PETSc::renumber));
    Global.Add("set", "(", new PETSc::setOptions<Dbddc>());
    addInv<Dbddc, PETSc::InvPETSc, KN<PetscScalar>, PetscScalar>();
    Init_Common();

    Global.Add("check", "(", new OneOperator1<bool, Dmat*>(CheckPetscMatrix<Dmat>));
}
#ifndef PETScandSLEPc
LOADFUNC(Init_PETSc)
#endif
