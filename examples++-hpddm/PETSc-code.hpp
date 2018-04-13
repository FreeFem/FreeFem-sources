#include "petsc.h"

#include "PETSc.hpp"

typedef PETSc::DistributedCSR<HpSchwarz<PetscScalar>> Dmat;
typedef PETSc::DistributedCSR<HpSchwarz<PetscReal>> DmatR;
typedef PETSc::DistributedCSR<HpSchwarz<PetscComplex>> DmatC;
typedef PETSc::DistributedCSR<HpSchur<PetscScalar>> Dbddc;
typedef PETSc::DistributedCSR<HpSchur<PetscReal>> DbddcR;
typedef PETSc::DistributedCSR<HpSchur<PetscComplex>> DbddcC;

namespace PETSc {
template<class HpddmType, typename std::enable_if<std::is_same<HpddmType, Dmat>::value>::type* = nullptr>
void initPETScStructure(HpddmType* ptA, MatriceMorse<PetscScalar>* mA, PetscInt bs, PetscBool symmetric, KN<typename std::conditional<std::is_same<HpddmType, Dmat>::value, double, long>::type>* ptD, KN<PetscScalar>* rhs) {
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
    bool free = ptA->_A->distributedCSR(ptA->_num, ptA->_first, ptA->_last, ia, ja, c);
    MatCreate(PETSC_COMM_WORLD, &(ptA->_petsc));
    if(bs > 1)
        MatSetBlockSize(ptA->_petsc, bs);
    MatSetSizes(ptA->_petsc, ptA->_last - ptA->_first, ptA->_last - ptA->_first, global, global);
    MatSetType(ptA->_petsc, MATMPIAIJ);
    MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
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
void initPETScStructure(HpddmType* ptA, MatriceMorse<PetscScalar>* mA, PetscInt& bs, PetscBool symmetric, KN<typename std::conditional<std::is_same<HpddmType, Dmat>::value, double, long>::type>* ptD, KN<PetscScalar>* rhs) {
    const HPDDM::MatrixCSR<PetscScalar>* M = ptA->_A->getMatrix();
    if(!M->_sym)
        std::cout << "Please assemble a symmetric CSR" << std::endl;
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
    ISLocalToGlobalMappingCreate(PETSC_COMM_WORLD, bs, M->_n / bs, indices, PETSC_OWN_POINTER, &(ptA->_rmap));
    MatCreateIS(PETSC_COMM_WORLD, bs, PETSC_DECIDE, PETSC_DECIDE, global, global, ptA->_rmap, NULL, &(ptA->_petsc));
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
    if(A) {
        numbering->resize(A->_A->getMatrix()->_n);
        if(A->_num)
            for(int i = 0; i < numbering->n; ++i)
                numbering->operator[](i) = A->_num[i];
    }
    return 0L;
}
template<class Type>
long changeOperator(Type* const& A, Matrice_Creuse<PetscScalar>* const& mat) {
    if(A && mat) {
        MatriceMorse<PetscScalar>* mN = static_cast<MatriceMorse<PetscScalar>*>(&*(mat->A));
        HPDDM::MatrixCSR<PetscScalar>* dN = new HPDDM::MatrixCSR<PetscScalar>(mN->n, mN->m, mN->nbcoef, mN->a, mN->lg, mN->cl, mN->symetrique);
        A->_A->setMatrix(dN);
        int* ia = nullptr;
        int* ja = nullptr;
        PetscScalar* c = nullptr;
        KSPSetOperators(A->_ksp, NULL, NULL);
        bool free = A->_A->distributedCSR(A->_num, A->_first, A->_last, ia, ja, c);
        for(PetscInt i = 0; i < A->_last - A->_first; ++i) {
            PetscInt row = A->_first + i;
            MatSetValues(A->_petsc, 1, &row, ia[i + 1] - ia[i], ja + ia[i], c + ia[i], INSERT_VALUES);
        }
        if(free) {
            delete [] ia;
            delete [] ja;
            delete [] c;
        }
        MatAssemblyBegin(A->_petsc, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(A->_petsc, MAT_FINAL_ASSEMBLY);
        KSPSetOperators(A->_ksp, A->_petsc, A->_petsc);
        if(std::is_same<Type, Dmat>::value && !A->_S.empty()) {
            KSPSetFromOptions(A->_ksp);
            KSPSetUp(A->_ksp);
            PC pc;
            KSPGetPC(A->_ksp, &pc);
            PCType type;
            PCGetType(pc, &type);
            PetscBool isFieldSplit;
            PetscStrcmp(type, PCFIELDSPLIT, &isFieldSplit);
            if(isFieldSplit) {
                setCompositePC(A, pc);
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
        HPDDM::MatrixCSR<PetscScalar>* dA = new HPDDM::MatrixCSR<PetscScalar>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->_A->setMatrix(dA);
        ptA->_num = new unsigned int[mA->n];
        std::copy_n(ptB->_num, dA->_n, ptA->_num);
        ptA->_first = ptB->_first;
        ptA->_last = ptB->_last;
        PetscInt bs;
        MatGetBlockSize(ptB->_petsc, &bs);
        KN<PetscScalar>* rhs = nargs[0] ? GetAny<KN<PetscScalar>*>((*nargs[0])(stack)) : nullptr;
        initPETScStructure(ptA, mA, bs, nargs[2] ? (GetAny<bool>((*nargs[2])(stack)) ? PETSC_TRUE : PETSC_FALSE) : PETSC_FALSE, static_cast<KN<double>*>(nullptr), rhs);
        KSPCreate(PETSC_COMM_WORLD, &(ptA->_ksp));
        KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
    }
    bool clean = nargs[1] && GetAny<bool>((*nargs[1])(stack));
    if(clean)
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
    MatCreate(PETSC_COMM_WORLD, &(ptA->_petsc));
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
    if(bsr)
        MatMPIBAIJSetPreallocationCSR(ptA->_petsc, bs, mK->lg, mK->cl, mK->a);
    else
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, mK->lg, mK->cl, mK->a);
    if(prune) {
        MatSetType(ptA->_S[0], MATMPIAIJ);
        int* pI = new int[mK->n * bs + 1];
        int* pJ = new int[mK->nbcoef * bs];
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
        MatMPIAIJSetPreallocationCSR(ptA->_S[0], pI, pJ, pC);
        delete [] pC;
        delete [] pJ;
        delete [] pI;
    }
    if(clean)
        ptK->destroy();
    MatSetOption(ptA->_petsc, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
    if(nargs[2])
        MatSetOption(ptA->_petsc, MAT_SYMMETRIC, GetAny<bool>((*nargs[2])(stack)) ? PETSC_TRUE : PETSC_FALSE);
    KSPCreate(PETSC_COMM_WORLD, &(ptA->_ksp));
#ifdef PETSC_HAVE_MKL_SPARSE_OPTIMIZE
    if(0 && prune && bsr)
        MatConvert(ptA->_petsc, MATMPIBAIJMKL, MAT_INPLACE_MATRIX, &(ptA->_petsc));
#endif
    KSPSetOperators(ptA->_ksp, ptA->_petsc, prune ? ptA->_S[0] : ptA->_petsc);
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
        int* ia = new int[n + 1];
        for(int i = 0; i < n + 1; ++i) {
            ia[i] = 0;
            for(int k = 0; k < size; ++k) {
                MatriceMorse<PetscScalar> *mA = static_cast<MatriceMorse<PetscScalar>*>(&(*(ptK->operator[](k)).A));
                if(i < mA->n + 1)
                    ia[i] += mA->lg[i];
            }
        }
        int* ja = new int[ia[n]];
        PetscScalar* c = new PetscScalar[ia[n]];
        int nnz = 0;
        int offset = (ptJ ? std::accumulate(dims, dims + v.begin()->first, 0) : 0);
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
        delete [] dims;
        MatCreate(PETSC_COMM_WORLD, &(ptA->_petsc));
        if(bs > 1)
            MatSetBlockSize(ptA->_petsc, bs);
        if(!ptJ && (mpisize / ptK->n > 1))
            MatSetSizes(ptA->_petsc, n, n, PETSC_DECIDE, PETSC_DECIDE);
        else
            MatSetSizes(ptA->_petsc, n, ptK->operator[](ptJ ? it->second : rank).M(), PETSC_DECIDE, PETSC_DECIDE);
        ptA->_first = 0;
        ptA->_last = n;
        bool clean = nargs[4] && GetAny<bool>((*nargs[4])(stack));
        if(clean) {
            int* cl = nullptr;
            int* lg = nullptr;
            PetscScalar* a = nullptr;
            for(int k = 0; k < size; ++k) {
                MatriceMorse<PetscScalar> *mA = static_cast<MatriceMorse<PetscScalar>*>(&(*(ptK->operator[](ptJ ? v[k].second : k)).A));
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
            }
            ptK->resize(0);
        }
        MatSetType(ptA->_petsc, MATMPIAIJ);
        MatMPIAIJSetPreallocationCSR(ptA->_petsc, ia, ja, c);
        delete [] ia;
        delete [] ja;
        delete [] c;
        MatSetOption(ptA->_petsc, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
        if(nargs[3])
            MatSetOption(ptA->_petsc, MAT_SYMMETRIC, GetAny<bool>((*nargs[3])(stack)) ? PETSC_TRUE : PETSC_FALSE);
        if(verbosity > 0 && mpirank == 0)
            cout << " --- global CSR created (in " << MPI_Wtime() - timing << ")" << endl;
        timing = MPI_Wtime();
        KSPCreate(PETSC_COMM_WORLD, &(ptA->_ksp));
        KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
    }
    return ptA;
}

template<class HpddmType>
class initCSR_Op : public E_F0mps {
    public:
        Expression A;
        Expression K;
        Expression O;
        Expression R;
        Expression D;
        static const int n_name_param = 5;
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
    {"rhs", &typeid(KN<PetscScalar>*)},
    {"clean", &typeid(bool)},
    {"symmetric", &typeid(bool)}
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
    PetscInt bs = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : 1;
    MatriceMorse<PetscScalar> *mA = static_cast<MatriceMorse<PetscScalar>*>(&(*GetAny<Matrice_Creuse<PetscScalar>*>((*K)(stack))->A));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    KN<PetscScalar>* rhs = nargs[2] ? GetAny<KN<PetscScalar>*>((*nargs[2])(stack)) : 0;
    ptA->_A = new HpddmType;
    if(comm)
        PETSC_COMM_WORLD = *comm;
    if(ptO && mA) {
        HPDDM::MatrixCSR<PetscScalar>* dA = new HPDDM::MatrixCSR<PetscScalar>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->_A->HPDDM::template Subdomain<PetscScalar>::initialize(dA, STL<long>(*ptO), *ptR, comm);
        ptA->_num = new unsigned int[ptA->_A->getMatrix()->_n];
        initPETScStructure(ptA, mA, bs, nargs[4] ? (GetAny<bool>((*nargs[4])(stack)) ? PETSC_TRUE : PETSC_FALSE) : PETSC_FALSE, ptD, rhs);
        if(!std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value)
            mA->lg = ptA->_A->getMatrix()->_ia;
        KSPCreate(PETSC_COMM_WORLD, &(ptA->_ksp));
        KSPSetOperators(ptA->_ksp, ptA->_petsc, ptA->_petsc);
        bool clean = nargs[3] && GetAny<bool>((*nargs[3])(stack));
        if(clean) {
            ptO->resize(0);
            ptR->resize(0);
            GetAny<Matrice_Creuse<PetscScalar>*>((*K)(stack))->destroy();
        }
    }
    return ptA;
}

template<class Type>
class setOptions_Op : public E_F0mps {
    public:
        Expression A;
        static const int n_name_param = 7;
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
    {"schurList", &typeid(KN<double>*)}
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
    bool fieldsplit = PETSc::insertOptions(options);
    if(std::is_same<Type, Dmat>::value && fieldsplit) {
        KN<double>* fields = nargs[2] ? GetAny<KN<double>*>((*nargs[2])(stack)) : 0;
        KN<String>* names = nargs[3] ? GetAny<KN<String>*>((*nargs[3])(stack)) : 0;
        KN<Matrice_Creuse<PetscScalar>>* mS = nargs[5] ? GetAny<KN<Matrice_Creuse<PetscScalar>>*>((*nargs[5])(stack)) : 0;
        KN<double>* pL = nargs[6] ? GetAny<KN<double>*>((*nargs[6])(stack)) : 0;
        setFieldSplitPC(ptA, ptA->_ksp, fields, names, mS, pL);
    }
    if(std::is_same<Type, Dmat>::value) {
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
                HPDDM::Subdomain<PetscScalar>::template distributedVec<0>(ptA->_num, ptA->_first, ptA->_last, static_cast<PetscScalar*>(*(ptNS->get(i))), x, ptNS->get(i)->n);
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
    if(nargs[4])
        KSPSetOptionsPrefix(ptA->_ksp, GetAny<std::string*>((*nargs[4])(stack))->c_str());
    KSPSetFromOptions(ptA->_ksp);
    KSPSetUp(ptA->_ksp);
    if(std::is_same<Type, Dmat>::value && nargs[2] && nargs[5] && nargs[6]) {
        PC pc;
        KSPGetPC(ptA->_ksp, &pc);
        setCompositePC(ptA, pc);
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
    MatType type;
    MatGetType(ptA->_petsc, &type);
    PetscInt bs = 1;
    PetscBool isBlock;
    PetscStrcmp(type, MATMPIBAIJ, &isBlock);
#ifdef PETSC_HAVE_MKL_SPARSE_OPTIMIZE
    if(!isBlock)
        PetscStrcmp(type, MATMPIBAIJMKL, &isBlock);
#endif
    if(isBlock)
        MatGetBlockSize(ptA->_petsc, &bs);
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
    HPDDM::Subdomain<PetscScalar>::template distributedVec<1>(ptA->_num, ptA->_first, ptA->_last, static_cast<PetscScalar*>(*ptX), ptr_y, ptX->n / bs, bs);
    VecRestoreArray(y, &ptr_y);
    VecDestroy(&y);
    if(ptA->_A)
        ptA->_A->HPDDM::template Subdomain<PetscScalar>::exchange(static_cast<PetscScalar*>(*ptX));
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
            Vec x, y;
            double timing = MPI_Wtime();
            MatCreateVecs((*t)._petsc, &x, &y);
            PetscScalar* ptr;
            PetscBool isBlock;
            if(std::is_same<typename std::remove_reference<decltype(*t.A->_A)>::type, HpSchwarz<PetscScalar>>::value) {
                VecGetArray(x, &ptr);
                MatType type;
                MatGetType((*t)._petsc, &type);
                PetscInt bs = 1;
                PetscStrcmp(type, MATMPIBAIJ, &isBlock);
#ifdef PETSC_HAVE_MKL_SPARSE_OPTIMIZE
                if(!isBlock)
                    PetscStrcmp(type, MATMPIBAIJMKL, &isBlock);
#endif
                if(isBlock)
                    MatGetBlockSize((*t)._petsc, &bs);
                HPDDM::Subdomain<K>::template distributedVec<0>((*t)._num, (*t)._first, (*t)._last, static_cast<PetscScalar*>(*u), ptr, u->n / bs, bs);
                VecRestoreArray(x, &ptr);
                if(t.A->_A)
                    std::fill_n(static_cast<PetscScalar*>(*out), out->n, 0.0);
            }
            else {
                PetscScalar zero = 0.0;
                VecSet(x, zero);
#if 0
                Mat_IS* is = (Mat_IS*)(*t)._petsc->data;
                VecGetArray(is->y, &ptr);
                std::copy_n(static_cast<PetscScalar*>(*u), (*t)._A->getMatrix()->_n, ptr);
                VecRestoreArray(is->y, &ptr);
                VecScatterBegin(is->rctx,is->y,x,ADD_VALUES,SCATTER_REVERSE);
                VecScatterEnd(is->rctx,is->y,x,ADD_VALUES,SCATTER_REVERSE);
#else
                VecGetArray((*t)._isVec, &ptr);
                std::copy_n(static_cast<PetscScalar*>(*u), (*t)._A->getMatrix()->_n, ptr);
                VecRestoreArray((*t)._isVec, &ptr);
                VecScatterBegin((*t)._scatter, (*t)._isVec, x, ADD_VALUES, SCATTER_REVERSE);
                VecScatterEnd((*t)._scatter, (*t)._isVec, x, ADD_VALUES, SCATTER_REVERSE);
#endif
            }
            timing = MPI_Wtime();
            KSPSolve((*t)._ksp, x, y);
            if(verbosity > 0 && mpirank == 0)
                cout << " --- system solved with PETSc (in " << MPI_Wtime() - timing << ")" << endl;
            if(std::is_same<typename std::remove_reference<decltype(*t.A->_A)>::type, HpSchwarz<PetscScalar>>::value) {
                VecGetArray(y, &ptr);
                MatType type;
                MatGetType((*t)._petsc, &type);
                PetscInt bs = 1;
                if(isBlock)
                    MatGetBlockSize((*t)._petsc, &bs);
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
                VecScatterBegin((*t)._scatter, y, (*t)._isVec, INSERT_VALUES, SCATTER_FORWARD);
                VecScatterEnd((*t)._scatter, y, (*t)._isVec, INSERT_VALUES, SCATTER_FORWARD);
                VecGetArray((*t)._isVec, &ptr);
                std::copy_n(ptr, (*t)._A->getMatrix()->_n, static_cast<PetscScalar*>(*out));
                VecRestoreArray((*t)._isVec, &ptr);
#endif
            }
            VecDestroy(&x);
            VecDestroy(&y);
            if(std::is_same<typename std::remove_reference<decltype(*t.A->_A)>::type, HpSchwarz<PetscScalar>>::value && t.A->_A)
                (*t)._A->HPDDM::template Subdomain<PetscScalar>::exchange(*out);
        };
        static U init(U Ax, InvPETSc<T, U, K> A) {
            A.solve(Ax);
            return Ax;
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
            Vec x, y;
            PetscScalar* ptr;
            if(N == 'N')
                MatCreateVecs((*t)._petsc, &x, &y);
            else
                MatCreateVecs((*t)._petsc, &y, &x);
            VecGetArray(x, &ptr);
            MatType type;
            MatGetType((*t)._petsc, &type);
            PetscInt bs = 1;
            PetscBool isBlock;
            PetscStrcmp(type, MATMPIBAIJ, &isBlock);
#ifdef PETSC_HAVE_MKL_SPARSE_OPTIMIZE
            if(!isBlock)
                PetscStrcmp(type, MATMPIBAIJMKL, &isBlock);
#endif
            if(isBlock)
                MatGetBlockSize((*t)._petsc, &bs);
            HPDDM::Subdomain<K>::template distributedVec<0>(t->_num, t->_first, t->_last, static_cast<PetscScalar*>(*u), ptr, u->n / bs, bs);
            VecRestoreArray(x, &ptr);
            if(N == 'T')
                MatMultTranspose(t->_petsc, x, y);
            else
                MatMult(t->_petsc, x, y);
            VecGetArray(y, &ptr);
            if(!t->_A)
                std::fill_n(static_cast<PetscScalar*>(*out), out->n, 0.0);
            HPDDM::Subdomain<K>::template distributedVec<1>(t->_num, t->_first, t->_last, static_cast<PetscScalar*>(*out), ptr, out->n / bs, bs);
            if(t->_A)
                (*t)._A->HPDDM::template Subdomain<PetscScalar>::exchange(*out);
            VecRestoreArray(y, &ptr);
            VecDestroy(&y);
        }
        static U mv(U Ax, ProdPETSc<T, U, K, N> A) {
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

template<typename T>
inline bool exist_type() {
    map<const string,basicForEachType*>::iterator ir = map_type.find(typeid(T).name());
    return ir != map_type.end();
}

static void Init_PETSc() {
    if(verbosity > 1 && mpirank == 0) {
        cout << " PETSc ( " << typeid(PetscScalar).name() << " )" << endl;
    }
    int argc = pkarg->n;
    char** argv = new char*[argc];
    for(int i = 0; i < argc; ++i)
        argv[i] = const_cast<char*>((*(*pkarg)[i].getap())->c_str());
    PetscInitialize(&argc, &argv, 0, "");
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
        zzzfff->Add("dmatrix", atype<DmatR*>());
    }
    if(!exist_type<DmatC*>()) {
        Dcl_Type<DmatC*>(Initialize<DmatC>, Delete<DmatC>);
        zzzfff->Add("zmatrix", atype<DmatC*>());
    }
    if(!exist_type<DbddcR*>()) {
        Dcl_Type<DbddcR*>(Initialize<DbddcR>, Delete<DbddcR>);
        zzzfff->Add("dbddc", atype<DbddcR*>());
    }
    if(!exist_type<DbddcC*>()) {
        Dcl_Type<DbddcC*>(Initialize<DbddcC>, Delete<DbddcC>);
        zzzfff->Add("zbddc", atype<DbddcC*>());
    }
    map_type_of_map[make_pair(atype<DmatR*>(),atype<Complex*>())] = atype<DmatC*>();
    map_type_of_map[make_pair(atype<DmatR*>(),atype<double*>())] = atype<DmatR*>();
    map_type_of_map[make_pair(atype<DbddcR*>(),atype<Complex*>())] = atype<DbddcC*>();
    map_type_of_map[make_pair(atype<DbddcR*>(),atype<double*>())] = atype<DbddcR*>();

    TheOperators->Add("<-", new OneOperator1_<long, Dmat*>(PETSc::initEmptyCSR<Dmat>));
    TheOperators->Add("<-", new PETSc::initCSR<HpSchwarz<PetscScalar>>);
    TheOperators->Add("<-", new PETSc::initCSRfromArray<HpSchwarz<PetscScalar>>);
    TheOperators->Add("<-", new PETSc::initCSRfromMatrix<HpSchwarz<PetscScalar>>);
    TheOperators->Add("<-", new PETSc::initCSRfromDMatrix<HpSchwarz<PetscScalar>>);
    Global.Add("set", "(", new PETSc::setOptions<Dmat>());
    addProd<Dmat, PETSc::ProdPETSc, KN<PetscScalar>, PetscScalar, 'N'>();
    addProd<Dmat, PETSc::ProdPETSc, KN<PetscScalar>, PetscScalar, 'T'>();
    addInv<Dmat, PETSc::InvPETSc, KN<PetscScalar>, PetscScalar>();

    TheOperators->Add("<-", new OneOperator1_<long, Dbddc*>(PETSc::initEmptyCSR<Dbddc>));
    TheOperators->Add("<-", new PETSc::initCSR<HpSchur<PetscScalar>>);

    Global.Add("exchange", "(", new exchangeIn<Dmat, PetscScalar>);
    Global.Add("exchange", "(", new exchangeInOut<Dmat, PetscScalar>);
    Global.Add("globalNumbering", "(", new OneOperator2_<long, Dmat*, KN<long>*>(PETSc::globalNumbering<Dmat>));
    Global.Add("globalNumbering", "(", new OneOperator2_<long, Dbddc*, KN<long>*>(PETSc::globalNumbering<Dbddc>));
    Global.Add("changeOperator", "(", new OneOperator2_<long, Dmat*, Matrice_Creuse<PetscScalar>*>(PETSc::changeOperator<Dmat>));
    Global.Add("IterativeMethod", "(", new PETSc::IterativeMethod<Dmat>);
    Global.Add("originalNumbering", "(", new OneOperator3_<long, Dbddc*, KN<PetscScalar>*, KN<long>*>(PETSc::originalNumbering));
    Global.Add("set", "(", new PETSc::setOptions<Dbddc>());
    addInv<Dbddc, PETSc::InvPETSc, KN<PetscScalar>, PetscScalar>();

    Global.Add("check", "(", new OneOperator1<bool, Dmat*>(CheckPetscMatrix<Dmat>));
}
#ifndef PETScandSLEPc
LOADFUNC(Init_PETSc)
#endif
