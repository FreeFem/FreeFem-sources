#ifndef PETSC_HPP_
#define PETSC_HPP_

#if PETSC_VERSION_LT(3,7,0)
#define FFPetscOptionsInsert(a,b,c) PetscOptionsInsert(a,b,c)
#else
#define FFPetscOptionsInsert(a,b,c) PetscOptionsInsert(NULL,a,b,c)
#endif

#if PETSC_VERSION_LT(3,6,0)
#define MatCreateVecs MatGetVecs
#endif

#define HPDDM_SCHWARZ 0
#define HPDDM_FETI    0
#define HPDDM_BDD     0
#define HPDDM_PETSC   1

#include "common.hpp"

namespace PETSc {
template<class HpddmType>
class DistributedCSR {
    public:
        HpddmType*                  _A;
        Mat                     _petsc;
        std::vector<Mat>            _S;
        ISLocalToGlobalMapping   _rmap;
        VecScatter            _scatter;
        Vec                     _isVec;
        KSP                       _ksp;
        unsigned int*             _num;
        unsigned int            _first;
        unsigned int             _last;
        DistributedCSR() : _A(), _petsc(), _ksp(), _num(), _first(), _last() { _S.clear(); };
        ~DistributedCSR() {
            MatDestroy(&_petsc);
            for(int i = 0; i < _S.size(); ++i)
                MatDestroy(&_S[i]);
            KSPDestroy(&_ksp);
            if(_A) {
                if(!std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value) {
                    ISLocalToGlobalMappingDestroy(&_rmap);
                    VecDestroy(&_isVec);
                    VecScatterDestroy(&_scatter);
                }
                else
                    _A->clearBuffer();
                delete _A;
                _A = nullptr;
            }
            delete [] _num;
            _num = nullptr;
        }
};

template<class HpddmType, typename std::enable_if<std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value>::type* = nullptr>
void globalMapping(HpddmType* const& A, unsigned int*& num, unsigned int& start, unsigned int& end, unsigned int& global, unsigned int* const list) {
    num = new unsigned int[A->getMatrix()->_n];
    A->template globalMapping<'C'>(num, num + A->getMatrix()->_n, start, end, global, A->getScaling(), list);
}
template<class HpddmType, typename std::enable_if<!std::is_same<HpddmType, HpSchwarz<PetscScalar>>::value>::type* = nullptr>
void globalMapping(HpddmType* const& A, unsigned int* const& num, unsigned int& start, unsigned int& end, unsigned int& global, unsigned int* const list) { }
template<class Type>
void setFieldSplitPC(Type* ptA, KSP ksp, KN<double>* const& fields, KN<String>* const& names, KN<Matrice_Creuse<PetscScalar>>* const& mT, KN<double>* const& pL) {
    if(fields) {
        PC pc;
        KSPGetPC(ksp, &pc);
        PCSetType(pc, PCFIELDSPLIT);
        unsigned short* local = new unsigned short[fields->n + ptA->_last - ptA->_first];
        for(int i = 0; i < fields->n; ++i)
            local[i] = fields->operator[](i);
        unsigned short nb = *std::max_element(local, local + fields->n);
        MPI_Allreduce(MPI_IN_PLACE, &nb, 1, MPI_UNSIGNED_SHORT, MPI_MAX, PETSC_COMM_WORLD);
        local += fields->n;
        HPDDM::Subdomain<PetscScalar>::template distributedVec<0>(ptA->_num, ptA->_first, ptA->_last, local - fields->n, local, fields->n);
        unsigned long* counts = new unsigned long[nb]();
        for(unsigned int i = 0; i < ptA->_last - ptA->_first; ++i)
            ++counts[local[i] - 1];
        PetscInt* idx = new PetscInt[*std::max_element(counts, counts + nb)];
        for(unsigned short j = 0; j < nb; ++j) {
            IS is;
            unsigned short* pt = local;
            for(unsigned int i = 0; i < counts[j]; ++pt) {
                if(*pt == j + 1)
                    idx[i++] = ptA->_first + std::distance(local, pt);
            }
            ISCreateGeneral(PETSC_COMM_WORLD, counts[j], idx, PETSC_COPY_VALUES, &is);
            PCFieldSplitSetIS(pc, names && j < names->size() ? (*(names->operator[](j))).c_str() : NULL, is);
            ISDestroy(&is);
        }
        if(mT && mT->n > 0 && pL) {
            int *is, *js;
            PetscScalar *s;
            unsigned int* re = new unsigned int[pL->n];
            unsigned int nbSchur = 1;
            for(int i = 0; i < pL->n; ++i)
                re[i] = ((*pL)[i]) ? nbSchur++ : 0;
            nbSchur--;
            unsigned int* num;
            unsigned int start, end, global;
            ptA->_A->clearBuffer();
            globalMapping(ptA->_A, num, start, end, global, re);
            delete [] re;
            re = new unsigned int[2 * nbSchur];
            unsigned int* numSchur = re + nbSchur;
            for(int i = 0, j = 0; i < pL->n; ++i) {
                if((*pL)[i]) {
                    *numSchur++ = num[i];
                    re[static_cast<int>((*pL)[i]) - 1] = j++;
                }
            }
            numSchur -= nbSchur;
            delete [] num;
            ptA->_S.resize(mT->n);
            for(int k = 0; k < mT->n; ++k) {
                MatriceMorse<PetscScalar>* mS = static_cast<MatriceMorse<PetscScalar>*>(&(*(mT->operator[](k)).A));
                std::vector<std::vector<std::pair<int, PetscScalar>>> tmp(mS->n);
                for(int i = 0; i < mS->n; ++i) {
                    unsigned int row = re[i];
                    tmp[row].reserve(mS->lg[i + 1] - mS->lg[i]);
                    for(int j = mS->lg[i]; j < mS->lg[i + 1]; ++j)
                        tmp[row].emplace_back(re[mS->cl[j]], mS->a[j]);
                    std::sort(tmp[row].begin(), tmp[row].end(), [](const std::pair<int, PetscScalar>& lhs, const std::pair<int, PetscScalar>& rhs) { return lhs.first < rhs.first; });
                }
                is = new int[mS->n + 1];
                js = new int[mS->nbcoef];
                s = new PetscScalar[mS->nbcoef];
                is[0] = 0;
                for(int i = 0; i < mS->n; ++i) {
                    for(int j = 0; j < tmp[i].size(); ++j) {
                        *js++ = tmp[i][j].first;
                        *s++ = tmp[i][j].second;
                    }
                    is[i + 1] = is[i] + tmp[i].size();
                }
                js -= mS->nbcoef;
                s -= mS->nbcoef;
                int* ia, *ja;
                PetscScalar* c;
                ia = ja = nullptr;
                c = nullptr;
                HPDDM::MatrixCSR<PetscScalar>* dN = new HPDDM::MatrixCSR<PetscScalar>(mS->n, mS->m, mS->nbcoef, s, is, js, mS->symetrique, true);
                bool free = ptA->_A->HPDDM::template Subdomain<PetscScalar>::distributedCSR(numSchur, start, end, ia, ja, c, dN);
                MatCreate(PETSC_COMM_WORLD, &ptA->_S[k]);
                MatSetSizes(ptA->_S[k], end - start, end - start, global, global);
                MatSetType(ptA->_S[k], MATMPIAIJ);
                MatMPIAIJSetPreallocationCSR(ptA->_S[k], ia, ja, c);
                MatSetOption(ptA->_S[k], MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
                if(free) {
                    delete [] ia;
                    delete [] ja;
                    delete [] c;
                }
                delete dN;
            }
            delete [] re;
            ptA->_A->setBuffer();
        }
        delete [] idx;
        delete [] counts;
        local -= fields->n;
        delete [] local;
    }
}
template<class Type>
void setCompositePC(Type* A, PC pc) {
    PetscInt nsplits;
    KSP* subksp;
    PCFieldSplitGetSubKSP(pc, &nsplits, &subksp);
    if(A->_S.size() == 1)
        KSPSetOperators(subksp[nsplits - 1], A->_S[0], A->_S[0]);
    else {
        PC pcS;
        KSPGetPC(subksp[nsplits - 1], &pcS);
        for(int i = 0; i < A->_S.size(); ++i)
            PCCompositeAddPC(pcS, PCNONE);
        PCSetUp(pcS);
        for(int i = 0; i < A->_S.size(); ++i) {
            PC subpc;
            PCCompositeGetPC(pcS, i, &subpc);
            PCSetOperators(subpc, A->_S[i], A->_S[i]);
            PCSetFromOptions(subpc);
        }
    }
    PetscFree(subksp);
}
bool insertOptions(std::string* const& options) {
    bool fieldsplit = false;
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
            if(!fieldsplit)
                fieldsplit = (elems[i].compare("fieldsplit") == 0);
        }
        FFPetscOptionsInsert(&argc, &data, NULL);
        delete [] *data;
        delete [] data;
    }
    return fieldsplit;
}
}
#endif
