#ifndef _ALL_IN_ONE_
#include <math.h>
#include <mpi.h>
#include <ff++.hpp>
#include "AFunction_ext.hpp"
#ifdef BDD
#undef BDD
#endif
#ifdef FETI
#undef FETI
#endif
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
#endif

extern KN<String>* pkarg;

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
    {"deflation", &typeid(FEbaseArrayKn<K>*)}
};
template<class Type, class K>
class initDDM : public OneOperator {
    public:
        initDDM() : OneOperator(atype<Type*>(), atype<Type*>(), atype<Matrice_Creuse<K>*>(), atype<KN<long>*>(), atype<KN<KN<long> >*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new initDDM_Op<Type, K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
        }
};
template<class Type, class K>
AnyType initDDM_Op<Type, K>::operator()(Stack stack) const {
    const char** argv = new const char*[pkarg->n];
    for(int i = 0; i < pkarg->n; ++i)
        argv[i] = (*((*pkarg)[i].getap()))->data();
    HPDDM::Option::get()->parse(pkarg->n, argv, mpirank == 0);
    delete [] argv;
    Type* ptA = GetAny<Type*>((*A)(stack));
    MatriceMorse<K>* mA = static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*Mat)(stack))->A));
    KN<long>* ptO = GetAny<KN<long>*>((*o)(stack));
    KN<KN<long> >* ptR = GetAny<KN<KN<long> >*>((*R)(stack));
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
    {"timing", &typeid(KN<double>*)},
    {"ret", &typeid(Pair<K>*)},
    {"solver", &typeid(long)}
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
    KN<double>* timing = nargs[4] ? GetAny<KN<double>*>((*nargs[4])(stack)) : 0;
    Pair<K>* pair = nargs[5] ? GetAny<Pair<K>*>((*nargs[5])(stack)) : 0;
    HPDDM::Option& opt = *HPDDM::Option::get();
    unsigned short nu = opt["geneo_nu"];
    typename HPDDM::Wrapper<K>::ul_type threshold = opt.val("geneo_threshold", 0.0);
    std::pair<MPI_Request, const K*>* ret = nullptr;
    double t;
    if(mA || nargs[6]) {
        nu = std::max(nu, static_cast<unsigned short>(1));
        long nbSolver = nargs[6] ? GetAny<long>((*nargs[6])(stack)) : 0;
        std::vector<const HPDDM::MatrixCSR<K>*> vecAIJ;
        if(mA) {
            HPDDM::MatrixCSR<K> dA(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
            MatriceMorse<K>* mB = nargs[1] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[1])(stack))->A)) : nullptr;
            MatriceMorse<K>* mP = nargs[2] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[2])(stack))->A)) : nullptr;
            t = MPI_Wtime();
            if(dA._n == dA._m) {
                const HPDDM::MatrixCSR<K>* const dP = mP ? new HPDDM::MatrixCSR<K>(mP->n, mP->m, mP->nbcoef, mP->a, mP->lg, mP->cl, mP->symetrique) : nullptr;
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
                delete dP;
            }
            else {
                vecAIJ.emplace_back(&dA);
                nbSolver = 101;
            }
        }
        else {
            ptA->template interaction<false, true>(vecAIJ);
            std::sort(vecAIJ.begin(), vecAIJ.end(), [](const HPDDM::MatrixCSR<K>* lhs, const HPDDM::MatrixCSR<K>* rhs) { return lhs->_m > rhs->_m; });
            nbSolver = 100;
        }
        if(nbSolver != 0) {
            ptA->callNumfact();
            if(!vecAIJ.empty()) {
                int dof = ptA->getDof();
                HPDDM::Eigensolver<K> solver(dof);
                const HPDDM::MatrixCSR<K>& first = *vecAIJ.front();
                nu = std::min(nu, static_cast<unsigned short>(first._m));
                K** ev = new K*[nu];
                *ev = new K[nu * dof];
                for(int i = 0; i < nu; ++i)
                    ev[i] = *ev + i * dof;
                ptA->setVectors(ev);
                ptA->Type::super::initialize(nu);
                int lwork = solver.workspace("S", &first._m);
                K* a;
                typename HPDDM::Wrapper<K>::ul_type* values;
                if(!std::is_same<K, typename HPDDM::Wrapper<K>::ul_type>::value) {
                    a = new K[first._m * (2 * dof + first._m) + lwork];
                    values = new typename HPDDM::Wrapper<K>::ul_type[nu + first._m + std::max(1, first._m * std::max(5 * first._m + 7, 2 * dof + 2 * first._m + 1))];
                }
                else {
                    a = new K[first._m * (2 * dof + first._m + 1) + lwork + nu];
                    values = reinterpret_cast<typename HPDDM::Wrapper<K>::ul_type*>(a + first._m * (2 * dof + first._m) + lwork);
                }
                int* pos = new int[nu + 8 * first._m];
                std::fill(pos, pos + nu, 0);
                std::fill(values, values + nu, 0.0);
                for(const HPDDM::MatrixCSR<K>* A : vecAIJ) {
                    K* u = a + dof * A->_m;
                    K* vt = u + dof * A->_m;
                    K* work = vt + A->_m * A->_m;
                    typename HPDDM::Wrapper<K>::ul_type* s = values + nu;
                    typename HPDDM::Wrapper<K>::ul_type* rwork = s + A->_m;
                    std::fill(a, a + A->_m * dof, K(0.0));
                    for(int i = 0; i < dof; ++i)
                        for(int j = A->_ia[i]; j < A->_ia[i + 1]; ++j)
                            a[i + A->_ja[j] * dof] = A->_a[j];
                    ptA->Type::super::callSolve(a, A->_m);
                    solver.svd("S", &(A->_m), a, s, u, vt, work, &lwork, pos + nu, rwork);
                    for(unsigned int i = 0, j = 0, k = 0; k < nu; ++k) {
                        if(s[i] > values[j])
                            pos[k] = ++i;
                        else
                            pos[k] = -(++j);
                    }
                    for(unsigned int i = nu - 1; i > 0; ) {
                        if(pos[i] < 0) {
                            unsigned int j = i;
                            while(j > 0 && pos[j - 1] < 0)
                                --j;
                            std::copy_backward(values - pos[j] - 1, values - pos[i], values + i + 1);
                            std::copy_backward(ev[-pos[j] - 1], ev[-pos[i] - 1] + dof, ev[i] + dof);
                            i = std::max(j, 1u) - 1;
                        }
                        else
                            --i;
                    }
                    for(unsigned int i = 0; i < nu; ) {
                        if(pos[i] > 0) {
                            unsigned int j = i;
                            while(j < nu - 1 && pos[j + 1] > 0)
                                ++j;
                            std::copy(s + pos[i] - 1, s + pos[j], values + i);
                            std::copy(u + (pos[i] - 1) * dof, u + pos[j] * dof, ev[i]);
                            i = j + 1;
                        }
                        else
                            ++i;
                    }
                }
                delete [] pos;
                if(!std::is_same<K, typename HPDDM::Wrapper<K>::ul_type>::value)
                    delete [] values;
                delete [] a;
                ptA->Type::super::initialize(nu);
                if(nbSolver == 100)
                    std::for_each(vecAIJ.begin(), vecAIJ.end(), std::default_delete<const HPDDM::MatrixCSR<K> >());
            }
            else
                ptA->Type::super::initialize(0);
        }
        if(timing)
            (*timing)[3] = MPI_Wtime() - t;
        opt["geneo_nu"] = nu;
        MPI_Barrier(MPI_COMM_WORLD);
        if(timing)
            t = MPI_Wtime();
        if(ptA->exclusion(comm)) {
            if(pair)
                pair->p = ptA->template buildTwo<1>(comm);
            else
                ret = ptA->template buildTwo<1>(comm);
        }
        else {
            if(pair)
                pair->p = ptA->template buildTwo<0>(comm);
            else
                ret = ptA->template buildTwo<0>(comm);
        }
        if(timing)
            (*timing)[4] = MPI_Wtime() - t;
    }
    else {
        MPI_Barrier(MPI_COMM_WORLD);
        if(!threshold)
            ret = ptA->template buildTwo<2>(comm);
    }
    if(ret)
        delete ret;
    return 0L;
}

template<class Type, class K>
class solveDDM_Op : public E_F0mps {
    public:
        Expression A;
        Expression x;
        Expression rhs;
        static const int n_name_param = 8;
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
    if(ptX->n != ptRHS->n || ptRHS->n < ptA->getDof())
        return 0L;
    HPDDM::Option& opt = *HPDDM::Option::get();
    typename HPDDM::Wrapper<K>::ul_type eps = nargs[0] ? GetAny<typename HPDDM::Wrapper<K>::ul_type>((*nargs[0])(stack)) : -1.0;
    if(std::abs(eps + 1.0) > 1.0e-6)
        opt["tol"] = eps;
    int dim = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : -1;
    if(dim != -1)
        opt["gmres_restart"] = dim;
    int iter = nargs[2] ? GetAny<long>((*nargs[2])(stack)) : -1;
    if(iter != -1)
        opt["max_it"] = iter;
    KN<double>* timing = nargs[3] ? GetAny<KN<double>*>((*nargs[3])(stack)) : 0;
    long solver = nargs[4] ? GetAny<long>((*nargs[4])(stack)) : 0;
    if(std::abs(solver) == 1)
        opt["schwarz_method"] = 3;
    else if(std::abs(solver) == 8) {
        opt["schwarz_method"] = 2;
        if(!opt.set("schwarz_coarse_correction"))
            opt["schwarz_coarse_correction"] = 0;
    }
    else if(std::abs(solver) == 6) {
        opt["schwarz_method"] = 1;
        if(!opt.set("schwarz_coarse_correction"))
            opt["schwarz_coarse_correction"] = 0;
    }
    else if(std::abs(solver) == 3 && !opt.set("schwarz_coarse_correction"))
        opt["schwarz_coarse_correction"] = 0;
    Pair<K>* pair = nargs[6] ? GetAny<Pair<K>*>((*nargs[6])(stack)) : 0;
    if(solver >= 0 && pair)
        if(pair->p) {
            int flag;
            MPI_Test(&(pair->p->first), &flag, MPI_STATUS_IGNORE);
        }
    MatriceMorse<K>* mA = nargs[7] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[7])(stack))->A)) : 0;
    double timer = MPI_Wtime();
    if(solver >= 0 && solver != 36) {
        if(mpisize > 1 && (mA && (solver == 6 || solver == 8))) {
            HPDDM::MatrixCSR<K> dA(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
            ptA->callNumfact(&dA);
        }
        else if((solver != 100 && solver != 101) || ((solver == 100 || solver == 101) && mpisize == 1))
            ptA->callNumfact(nullptr);
        if(timing) (*timing)[1] = MPI_Wtime() - timer;
    }
    else if(solver == 36) {
        if(mpisize > 1) {
            HPDDM::MatrixCSR<K> dA(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
            ptA->callNumfact(&dA);
        }
    }
    bool excluded = nargs[5] ? GetAny<bool>((*nargs[5])(stack)) : false;
    if(excluded)
        opt["master_exclude"];
    if(pair)
        if(pair->p) {
            if(timing)
                timer = MPI_Wtime();
            MPI_Wait(&(pair->p->first), MPI_STATUS_IGNORE);
            if(timing)
                (*timing)[timing->n - 1] = MPI_Wtime() - timer;
            delete [] pair->p->second;
            delete pair->p;
            timer = MPI_Wtime();
        }
    MPI_Barrier(MPI_COMM_WORLD);
    if(solver >= 0 && !excluded && pair && pair->p && timing && mpisize > 1)
        (*timing)[timing->n - 1] += MPI_Wtime() - timer;
    unsigned short mu = ptX->n / ptA->getDof();
    MPI_Allreduce(MPI_IN_PLACE, &mu, 1, MPI_UNSIGNED_SHORT, MPI_MAX, ptA->getCommunicator());
    int rank;
    MPI_Comm_rank(ptA->getCommunicator(), &rank);
    if(rank != 0 || excluded)
        opt.remove("verbosity");
    timer = MPI_Wtime();
    if(!excluded) {
        if(opt["krylov_method"] == 1)
            HPDDM::IterativeMethod::CG(*ptA, (K*)*ptX, (K*)*ptRHS, MPI_COMM_WORLD);
        else
            HPDDM::IterativeMethod::GMRES(*ptA, (K*)*ptX, (K*)*ptRHS, mu, MPI_COMM_WORLD);
    }
    else {
        if(opt["krylov_method"] == 1)
            HPDDM::IterativeMethod::CG(*ptA, (K*)nullptr, (K*)nullptr, MPI_COMM_WORLD);
        else
            HPDDM::IterativeMethod::GMRES<true>(*ptA, (K*)nullptr, (K*)nullptr, mu, MPI_COMM_WORLD);
    }
    timer = MPI_Wtime() - timer;
    if(!excluded) {
        if(rank == 0)
            std::cout << scientific << " --- system solved (in " << timer << ")" << std::endl;
        typename HPDDM::Wrapper<K>::ul_type* storage = new typename HPDDM::Wrapper<K>::ul_type[2 * mu];
        ptA->computeError(*ptX, *ptRHS, storage, mu);
        if(rank == 0) {
            std::cout << scientific << " --- error = " << storage[1] << " / " << storage[0];
            for(unsigned short nu = 1; nu < mu; ++nu)
                std::cout << " (rhs #" << nu << ")" << std::endl << scientific << "             " << storage[2 * nu + 1] << " / " << storage[2 * nu];
            if(mu > 1)
                std::cout << " (rhs #" << mu << ")";
            std::cout << std::endl;
        }
        delete [] storage;
    }
    return 0L;
}

template<class Type, class K>
class changeOperator_Op : public E_F0mps {
    public:
        Expression A;
        Expression mat;
        static const int n_name_param = 0;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        changeOperator_Op(const basicAC_F0& args, Expression param1, Expression param2) : A(param1), mat(param2) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type changeOperator_Op<Type, K>::name_param[] = { };
template<class Type, class K>
class changeOperator : public OneOperator {
    public:
        changeOperator() : OneOperator(atype<long>(), atype<Type*>(), atype<Matrice_Creuse<K>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new changeOperator_Op<Type, K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
        }
};
template<class Type, class K>
AnyType changeOperator_Op<Type, K>::operator()(Stack stack) const {
    MatriceMorse<K>* mN = static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*mat)(stack))->A));
    HPDDM::MatrixCSR<K>* dN = new HPDDM::MatrixCSR<K>(mN->n, mN->m, mN->nbcoef, mN->a, mN->lg, mN->cl, mN->symetrique);
    Type* ptA = GetAny<Type*>((*A)(stack));
    ptA->setMatrix(dN);
    return 0L;
}

class distributedDot_Op : public E_F0mps {
    public:
        Expression A;
        Expression in;
        Expression out;
        static const int n_name_param = 1;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        distributedDot_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : A(param1), in(param2), out(param3) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
basicAC_F0::name_and_type distributedDot_Op::name_param[] = {
    {"communicator", &typeid(pcommworld)}
};
class distributedDot : public OneOperator {
    public:
        distributedDot() : OneOperator(atype<double>(), atype<KN<double>*>(), atype<KN<double>*>(), atype<KN<double>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new distributedDot_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
        }
};
AnyType distributedDot_Op::operator()(Stack stack) const {
    KN<double>* pA = GetAny<KN<double>*>((*A)(stack));
    KN<double>* pin = GetAny<KN<double>*>((*in)(stack));
    KN<double>* pout = GetAny<KN<double>*>((*out)(stack));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    double* tmp = new double[pin->n];
    HPDDM::Wrapper<double>::diag(pin->n, *pA, *pin, tmp);
    KN_<double> KN(tmp, pin->n);
    double dot = (KN, *pout);
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, MPI_DOUBLE, MPI_SUM, comm ? *((MPI_Comm*)comm) : MPI_COMM_WORLD);
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

#ifdef DEBUG_SCHWARZ
template<class Type>
long InitDDM(Type* const& A) {
    A->callNumfact();
    return 0;
}
template<class Type, class K>
long ApplyAS(Type* const& A, KN<K>* const& in, KN<K>* const& out) {
    K* tmp = new K[in->n];
    std::copy_n((K*)*in, in->n, tmp);
    A->apply(tmp, *out);
    delete [] tmp;
    return 0;
}
#endif

template<template<class, char> class Type, class K, char S>
void add() {
    Dcl_Type<Type<K, S>*>(Initialize<Type<K, S> >, Delete<Type<K, S> >);

    TheOperators->Add("<-", new initDDM<Type<K, S>, K>);
    Global.Add("attachCoarseOperator", "(", new attachCoarseOperator<Type<K, S>, K>);
    Global.Add("DDM", "(", new solveDDM<Type<K, S>, K>);
    Global.Add("changeOperator", "(", new changeOperator<Type<K, S>, K>);
    Dcl_Type<GMV<Type<K, S>*, KN<K>*> >();
    TheOperators->Add("*", new OneOperator2<GMV<Type<K, S>*, KN<K>*>, Type<K, S>*, KN<K>*>(Build));
    TheOperators->Add("=", new OneOperator2<KN<K>*, KN<K>*, GMV<Type<K, S>*, KN<K>*> >(GlobalMV));
    if(std::is_same<K, double>::value)
        Global.Add("dscalprod", "(", new distributedDot);
#ifdef DEBUG_SCHWARZ
    Global.Add("initPreconditioner", "(", new OneOperator1_<long, Type<K, S>*>(InitDDM));
    Global.Add("precond", "(", new OneOperator3_<long, Type<K, S>*, KN<K>*, KN<K>*>(ApplyAS));
#endif
}
}

#ifndef _ALL_IN_ONE_
static void Init_Schwarz() {
#include "init.hpp"
}

LOADFUNC(Init_Schwarz)
#endif
