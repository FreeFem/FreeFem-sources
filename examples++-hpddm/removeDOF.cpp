#ifndef _ALL_IN_ONE_
#include "ff++.hpp"
#include <algorithm>
#include <vector>
#include <cmath>
#endif
#include <utility>

#ifndef EPS
#define EPS 1e-12
#endif

template<class T>
class removeInteraction_Op : public E_F0mps {
    public:
        Expression A;
        Expression R;
        Expression x;
        Expression out;
        static const int n_name_param = 1;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        removeInteraction_Op(const basicAC_F0&  args, Expression param1, Expression param2, Expression param3, Expression param4) : A(param1), R(param2), x(param3), out(param4) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};

template<class T>
basicAC_F0::name_and_type removeInteraction_Op<T>::name_param[] = {
    {"interaction", &typeid(Matrice_Creuse<T>*)}
};

template<class T>
class removeInteraction : public OneOperator {
    public:
        removeInteraction() : OneOperator(atype<long>(), atype<Matrice_Creuse<T>*>(), atype<Matrice_Creuse<double>*>(), atype<KN<T>*>(), atype<KN<T>*>()) {}

        E_F0* code(const basicAC_F0& args) const {
            return new removeInteraction_Op<T>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
        }
};

template<class T>
AnyType removeInteraction_Op<T>::operator()(Stack stack)  const {
    Matrice_Creuse<T>* pA = GetAny<Matrice_Creuse<T>* >((*A)(stack));
    Matrice_Creuse<T>* pR = GetAny<Matrice_Creuse<T>* >((*R)(stack));
    KN<T>* pX = GetAny<KN<T>* >((*x)(stack));
    KN<T>* pOut = GetAny<KN<T>* >((*out)(stack));
    ffassert(pA && pR && pX && pOut);
    pA->Uh = pR->Uh;
    pA->Vh = pR->Vh;
    MatriceMorse<T> *mA = static_cast<MatriceMorse<T>*>(&(*pA->A));
    MatriceMorse<T> *mR = static_cast<MatriceMorse<T>*>(&(*pR->A));
    Matrice_Creuse<T> *pInteraction = nargs[0] ? GetAny<Matrice_Creuse<T>*>((*nargs[0])(stack)) : 0;
    MatriceMorse<T> *mInteraction = pInteraction ? static_cast<MatriceMorse<T>*>(&(*pInteraction->A)) : 0;

    unsigned int n = mR->nbcoef;
    if(pOut->n == 0)
        pOut->set(new T[n], n);
    if(n == mA->n && n == mA->m) {
        for(unsigned int i = 0; i < n; ++i)
            *(*pOut + i) = *(*pX + mR->cl[i]);
        return 0L;
    }
    int* lg = new int[n + 1];
    int* cl;
    T* val;

    std::vector<std::vector<std::pair<unsigned int, T> > > tmpInteraction(pInteraction ? mR->n : 0);
    for(unsigned int i = 0; i < tmpInteraction.size(); ++i)
        tmpInteraction[i].reserve(std::max(mA->lg[mR->cl[i] + 1] - mA->lg[mR->cl[i]] - 1, 0));
    std::vector<int> tmpVec;
    tmpVec.resize(mA->n);
    for(unsigned int i = 0; i < n; ++i)
        tmpVec[mR->cl[i]] = i + 1;
    unsigned int nnzInteraction = 0;
    if(pInteraction)
        for(unsigned int i = 0, j = 0; i < mA->n; ++i)
            if(tmpVec[i] == 0)
                tmpVec[i] = -(++j);
    if(!mA->symetrique) {
        std::vector<std::pair<int, T> > tmp;
        tmp.reserve(mA->nbcoef);
        lg[0] = 0;
        for(unsigned int i = 0; i < n; ++i) {
            for(unsigned int j = mA->lg[mR->cl[i]]; j < mA->lg[mR->cl[i] + 1]; ++j) {
                if(abs(mA->a[j]) > EPS) {
                    int col = tmpVec[mA->cl[j]];
                    if(col > 0)
                        tmp.emplace_back(col - 1, mA->a[j]);
                    else if(pInteraction) {
                        tmpInteraction[i].emplace_back(-col - 1, mA->a[j]);
                        ++nnzInteraction;
                    }
                }
            }
            std::sort(tmp.begin() + lg[i], tmp.end(), [](const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; });
            *(*pOut + i) = *(*pX + mR->cl[i]);
            lg[i + 1] = tmp.size();
        }
        mA->nbcoef = tmp.size();
        cl = new int[tmp.size()];
        val = new T[tmp.size()];
        for(unsigned int i = 0; i < tmp.size(); ++i) {
            cl[i]  = tmp[i].first;
            val[i] = tmp[i].second;
        }
    }
    else {
        std::vector<std::vector<std::pair<unsigned int, T> > > tmp(n);
        for(unsigned int i = 0; i < n; ++i)
            tmp[i].reserve(mA->lg[mR->cl[i] + 1] - mA->lg[mR->cl[i]]);

        unsigned int nnz = 0;
        for(unsigned int i = 0; i < mA->n; ++i) {
            unsigned int row = std::abs(tmpVec[i]) - 1;
            if(tmpVec[i] > 0) {
                for(unsigned int j = mA->lg[mR->cl[row]]; j < mA->lg[mR->cl[row] + 1]; ++j) {
                    if(abs(mA->a[j]) > EPS) {
                        int col = tmpVec[mA->cl[j]];
                        if(col > 0) {
                            if(row < col - 1)
                                tmp[col - 1].emplace_back(row, mA->a[j]);
                            else
                                tmp[row].emplace_back(col - 1, mA->a[j]);
                            ++nnz;
                        }
                        else if(pInteraction) {
                            tmpInteraction[row].emplace_back(-col - 1, mA->a[j]);
                            ++nnzInteraction;
                        }
                    }
                }
                *(*pOut + row) = *(*pX + mR->cl[row]);
            }
            else if(pInteraction) {
                for(unsigned int j = mA->lg[i]; j < mA->lg[i + 1]; ++j) {
                    if(abs(mA->a[j]) > EPS) {
                        int col = tmpVec[mA->cl[j]];
                        if(col > 0) {
                            tmpInteraction[col - 1].emplace_back(row, mA->a[j]);
                            ++nnzInteraction;
                        }
                    }
                }
            }
        }
        mA->nbcoef = nnz;
        cl = new int[nnz];
        val = new T[nnz];
        nnz = 0;
        lg[0] = 0;
        for(unsigned int i = 0; i < n; ++i) {
            std::sort(tmp[i].begin(), tmp[i].end(), [](const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; });
            for(typename std::vector<std::pair<unsigned int, T> >::const_iterator it = tmp[i].begin(); it != tmp[i].end(); ++it) {
                cl[nnz] = it->first;
                val[nnz++] = it->second;
            }
            lg[i + 1] = nnz;
        }
    }
    delete [] mA->cl;
    delete [] mA->lg;
    delete [] mA->a;
    int m = mA->n - n;
    mA->n = n;
    mA->m = n;
    mA->N = n;
    mA->M = n;
    mA->lg = lg;
    mA->cl = cl;
    mA->a = val;
    if(pInteraction) {
        int* lgInteraction = new int[tmpInteraction.size() + 1];
        lgInteraction[0] = 0;
        int* clInteraction = new int[nnzInteraction];
        T* valInteraction = new T[nnzInteraction];
        nnzInteraction = 0;
        for(unsigned int i = 0; i < tmpInteraction.size(); ++i) {
            if(mA->symetrique)
                std::sort(tmpInteraction[i].begin(), tmpInteraction[i].end(), [](const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; });
            for(const std::pair<unsigned int, T>& p : tmpInteraction[i]) {
                clInteraction[nnzInteraction] = p.first;
                valInteraction[nnzInteraction++] = p.second;
            }
            lgInteraction[i + 1] = nnzInteraction;
        }
        delete mInteraction;
        MatriceMorse<T>* mNew = new MatriceMorse<T>(n, m, nnzInteraction, false, valInteraction, lgInteraction, clInteraction, true);
        pInteraction->typemat = TypeSolveMat(TypeSolveMat::GMRES);
        pInteraction->A.master(mNew);
        mNew->dummy = false;
    }
    return 0L;
}

#ifndef _ALL_IN_ONE_
static void Init_Remove() {
    Global.Add("removeInteraction", "(", new removeInteraction<double>);
    Global.Add("removeInteraction", "(", new removeInteraction<std::complex<double> >);
}

LOADFUNC(Init_Remove)
#endif
