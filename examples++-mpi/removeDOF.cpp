//ff-c++-LIBRARY-dep: cxx11 
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
class removeDOF_Op : public E_F0mps {
    public:
        Expression A;
        Expression R;
        Expression x;
        Expression out;
        static const int n_name_param = 3;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        removeDOF_Op(const basicAC_F0&  args, Expression param1, Expression param2, Expression param3, Expression param4) : A(param1), R(param2), x(param3), out(param4) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};

template<class T>
basicAC_F0::name_and_type removeDOF_Op<T>::name_param[] = {
    {"symmetrize", &typeid(bool)},
    {"condensation", &typeid(KN<long>*)},
    {"interaction", &typeid(Matrice_Creuse<T>*)}
};

template<class T>
class removeDOF : public OneOperator {
    public:
        removeDOF() : OneOperator(atype<long>(), atype<Matrice_Creuse<T>*>(), atype<Matrice_Creuse<double>*>(), atype<KN<T>*>(), atype<KN<T>*>()) {}

        E_F0* code(const basicAC_F0& args) const {
            return new removeDOF_Op<T>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
        }
};

template<class T>
AnyType removeDOF_Op<T>::operator()(Stack stack)  const {
    Matrice_Creuse<T>* pA = GetAny<Matrice_Creuse<T>* >((*A)(stack));
    Matrice_Creuse<T>* pR = GetAny<Matrice_Creuse<T>* >((*R)(stack));
    KN<T>* pX = GetAny<KN<T>* >((*x)(stack));
    KN<T>* pOut = GetAny<KN<T>* >((*out)(stack));
    ffassert(pA && pR && pX && pOut);
    pA->Uh = pR->Uh;
    pA->Vh = pR->Vh;
    MatriceMorse<T> *mA = static_cast<MatriceMorse<T>*>(&(*pA->A));
    MatriceMorse<T> *mR = static_cast<MatriceMorse<T>*>(&(*pR->A));
    bool symmetrize = nargs[0] ? GetAny<bool>((*nargs[0])(stack)) : false;
    KN<long>* condensation = nargs[1] ? GetAny<KN<long>* >((*nargs[1])(stack)) : (KN<long>*) 0;

    unsigned int n = condensation ? condensation->n : mR->nbcoef;
    int* lg = new int[n + 1];
    int* cl;
    T* val;
    T* b;
    if(pOut->n == 0) {
        b = new T[n];
        pOut->set(b, n);
    }

    std::vector<signed int> tmpVec;
    if(!condensation) {
        tmpVec.resize(mA->n);
        for(unsigned int i = 0; i < n; ++i)
            tmpVec[mR->cl[i]] = i + 1;
        if(!mA->symetrique) {
            std::vector<std::pair<int, T> > tmp;
            tmp.reserve(mA->nbcoef);

            lg[0] = 0;
            for(unsigned int i = 0; i < n; ++i) {
                for(unsigned int j = mA->lg[mR->cl[i]]; j < mA->lg[mR->cl[i] + 1]; ++j) {
                    unsigned int col = tmpVec[mA->cl[j]];
                    if(col != 0 && abs(mA->a[j]) > EPS) {
                        if(symmetrize) {
                            if(col - 1 <= i)
                                tmp.emplace_back(col - 1, mA->a[j]);
                        }
                        else
                            tmp.emplace_back(col - 1, mA->a[j]);
                    }
                }
                std::sort(tmp.begin() + lg[i], tmp.end(), [](const std::pair<unsigned int, T>& lhs, const std::pair<unsigned int, T>& rhs) { return lhs.first < rhs.first; });
                *(*pOut + i) = *(*pX + mR->cl[i]);
                lg[i + 1] = tmp.size();
            }
            mA->nbcoef = tmp.size();
            if(symmetrize)
                mA->symetrique = true;
            else
                mA->symetrique = false;

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
            for(unsigned int i = 0; i < n; ++i) {
                for(unsigned int j = mA->lg[mR->cl[i]]; j < mA->lg[mR->cl[i] + 1]; ++j) {
                    unsigned int col = tmpVec[mA->cl[j]];
                    if(col != 0 && abs(mA->a[j]) > EPS) {
                        if(i < col - 1)
                            tmp[col - 1].emplace_back(i, mA->a[j]);
                        else
                            tmp[i].emplace_back(col - 1, mA->a[j]);
                        ++nnz;
                    }
                }
                *(*pOut + i) = *(*pX + mR->cl[i]);
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
        mA->n = n;
        mA->m = n;
        mA->N = n;
        mA->M = n;
        mA->lg = lg;
        mA->cl = cl;
        mA->a = val;
    }
    else {
        tmpVec.reserve(mA->n);
        unsigned int i = 0, j = 1;
        for(unsigned int k = 0; k < mA->n; ++k) {
            if(k == *(*condensation + i)) {
                ++i;
                tmpVec.emplace_back(i);
            }
            else {
                tmpVec.emplace_back(-j);
                ++j;
            }
        }


//        if(!mA->symetrique) {
            std::vector<std::pair<int, T> > tmpInterior;
            std::vector<std::pair<int, T> > tmpBoundary;
            std::vector<std::pair<int, T> > tmpInteraction;
            tmpInterior.reserve(mA->nbcoef);
            tmpBoundary.reserve(mA->nbcoef);
            tmpInteraction.reserve(mA->nbcoef);

            lg[0] = 0;
            for(unsigned int i = 0; i < mA->n; ++i) {
                int row = tmpVec[i];
                if(row < 0) {
                    for(unsigned int j = mA->lg[i]; j < mA->lg[i + 1]; ++j) {
                        int col = tmpVec[mA->cl[j]];
                        if(col < 0)
                            tmpInterior.emplace_back(-col - 1, mA->a[j]);
                        else
                            tmpInteraction.emplace_back(col - 1, mA->a[j]);
                    }

                }
                else {
                    for(unsigned int j = mA->lg[i]; j < mA->lg[i + 1]; ++j) {
                        int col = tmpVec[mA->cl[j]];
                        if(col > 0)
                            tmpBoundary.emplace_back(col - 1, mA->a[j]);
                    }
                    // std::sort(tmp.begin() + lg[i], tmp.end());
                    *(*pOut + i) = *(*pX + *(*condensation + i));
                    lg[i + 1] = tmpBoundary.size();
                }
            }
            cl = new int[tmpBoundary.size()];
            val = new T[tmpBoundary.size()];
            for(unsigned int i = 0; i < tmpBoundary.size(); ++i) {
                cl[i]  = tmpBoundary[i].first;
                val[i] = tmpBoundary[i].second;
            }
//        }
        MatriceMorse<T>* m = new MatriceMorse<T>(n, n, tmpBoundary.size(), mA->symetrique, val, lg, cl, true);
        pR->typemat = TypeSolveMat(TypeSolveMat::GMRES);
        pR->A.master(m);
        m->dummy = false;
    }
    return 0L;
}

template<class T>
long symmetrizeCSR(Matrice_Creuse<T>* const& A) {
    MatriceMorse<T>* mA = static_cast<MatriceMorse<T>*>(&(*A->A));
    if(!mA->symetrique) {
        mA->symetrique = true;
        std::vector<int> cl;
        std::vector<T> a;
        a.reserve(mA->nbcoef);
        cl.reserve(mA->nbcoef);
        unsigned int save = mA->lg[0];
        for(unsigned int i = 0; i < mA->n; ++i) {
            for(unsigned int j = save; j < mA->lg[i + 1]; ++j) {
                int col = mA->cl[j];
                if(col <= i) {
                    T val = mA->a[j];
                    if(abs(val) > EPS) {
                        a.emplace_back(val);
                        cl.emplace_back(col);
                    }
                }
                else
                    break;
            }
            save = mA->lg[i + 1];
            mA->lg[i + 1] = cl.size();
        }
        delete [] mA->cl;
        delete [] mA->a;
        int* col = new int[cl.size()];
        T* val = new T[cl.size()];
        for(unsigned int i = 0; i < cl.size(); ++i) {
            col[i] = cl[i];
            val[i] = a[i];
        }
        mA->cl = col;
        mA->a = val;
        mA->nbcoef = cl.size();
    }
    return 0L;
}

#ifndef _ALL_IN_ONE_
/* --FH:   class Init {
    public:
        Init();
	};
Init ...; */
static void Load_Initrm() {
    Global.Add("removeDOF", "(", new removeDOF<double>);
    Global.Add("removeDOF", "(", new removeDOF<std::complex<double>>);
    Global.Add("symmetrizeCSR", "(", new OneOperator1_<long, Matrice_Creuse<double>* >(symmetrizeCSR<double>));
}
 LOADFUNC(Load_Initrm)
#endif

// std::sort(tmp.begin() + lg[i], tmp.end(), [](const std::pair<int, T>& lhs, const std::pair<int, T>& rhs) { return lhs.first < rhs.first; } );

