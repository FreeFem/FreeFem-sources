//ff-c++-LIBRARY-dep: cxx11 hpddm [umfpack|mumps parmetis ptscotch scotch scalapack] [mkl|blas] mpi pthread mpifc fc
//ff-c++-cpp-dep:

#define HPDDM_SCHWARZ 1
#define HPDDM_FETI    0
#define HPDDM_BDD     0
#include "common.hpp"

namespace Schwarz {
double getOpt(string* const& ss) {
    return HPDDM::Option::get()->val(*ss);
}
bool isSetOpt(string* const& ss) {
    return HPDDM::Option::get()->set(*ss);
}
template<class Type, class K>
bool destroyRecycling(Type* const& Op) {
    HPDDM::Recycling<K>::get()->destroy(Op->prefix());
    return false;
}


template<class Type, class K>
class initDDM_Op : public E_F0mps {
    public:
        Expression A;
        Expression Mat;
        Expression o;
        Expression R;
        static const int n_name_param = 4;
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
    {"scaling", &typeid(KN<HPDDM::underlying_type<K>>*)},
    {"deflation", &typeid(FEbaseArrayKn<K>*)},
    {"prefix", &typeid(string*)}
};
template<class Type, class K>
class initDDM : public OneOperator {
    public:
        initDDM() : OneOperator(atype<Type*>(), atype<Type*>(), atype<Matrice_Creuse<K>*>(), atype<KN<long>*>(), atype<KN<KN<long>>*>()) { }

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
    KN<KN<long>>* ptR = GetAny<KN<KN<long>>*>((*R)(stack));
    if(ptO && mA) {
        HPDDM::MatrixCSR<K>* dA = new HPDDM::MatrixCSR<K>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->HPDDM::template Subdomain<K>::initialize(dA, STL<long>(*ptO), *ptR, nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0);
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
    if(nargs[1])
        ptA->initialize(*GetAny<KN<HPDDM::underlying_type<K>>*>((*nargs[1])(stack)));
    else
        std::cerr << "Something is really wrong here !" << std::endl;
    if(nargs[3])
        ptA->setPrefix(*(GetAny<string*>((*nargs[3])(stack))));
    return ptA;
}

template<class Type, class K>
class attachCoarseOperator_Op : public E_F0mps {
    public:
        Expression comm;
        Expression A;
        static const int n_name_param = 6;
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
    {"threshold", &typeid(HPDDM::underlying_type<K>)},
    {"timing", &typeid(KN<double>*)},
    {"ret", &typeid(Pair<K>*)}
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
    Pair<K>* pair = nargs[5] ? GetAny<Pair<K>*>((*nargs[5])(stack)) : 0;
    HPDDM::Option& opt = *HPDDM::Option::get();
    unsigned short nu = opt.val<unsigned short>("geneo_nu", 20);
    HPDDM::underlying_type<K> threshold = opt.val("geneo_threshold", 0.0);
    KN<double>* timing = nargs[4] ? GetAny<KN<double>*>((*nargs[4])(stack)) : 0;
    std::pair<MPI_Request, const K*>* ret = nullptr;
    double t;
    if(mA) {
        nu = std::max(nu, static_cast<unsigned short>(1));
        long nbSolver = 0;
        std::vector<const HPDDM::MatrixCSR<K>*> vecAIJ;
        if(mA) {
            HPDDM::MatrixCSR<K> dA(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
            MatriceMorse<K>* mB = nargs[1] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[1])(stack))->A)) : nullptr;
            MatriceMorse<K>* mP = nargs[2] && opt.any_of("schwarz_method", { 1, 2, 4 }) ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[2])(stack))->A)) : nullptr;
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
                int info;
                int lwork = -1;
                {
                    K wkopt;
                    HPDDM::Lapack<K>::gesdd("S", &dof, &first._m, nullptr, &dof, nullptr, nullptr, &dof, nullptr, &first._m, &wkopt, &lwork, nullptr, nullptr, &info);
                }
                K* a;
                HPDDM::underlying_type<K>* values;
                if(!std::is_same<K, HPDDM::underlying_type<K>>::value) {
                    a = new K[first._m * (2 * dof + first._m) + lwork];
                    values = new HPDDM::underlying_type<K>[nu + first._m + std::max(1, first._m * std::max(5 * first._m + 7, 2 * dof + 2 * first._m + 1))];
                }
                else {
                    a = new K[first._m * (2 * dof + first._m + 1) + lwork + nu];
                    values = reinterpret_cast<HPDDM::underlying_type<K>*>(a + first._m * (2 * dof + first._m) + lwork);
                }
                int* pos = new int[nu + 8 * first._m];
                std::fill(pos, pos + nu, 0);
                std::fill(values, values + nu, 0.0);
                for(const HPDDM::MatrixCSR<K>* A : vecAIJ) {
                    K* u = a + dof * A->_m;
                    K* vt = u + dof * A->_m;
                    K* work = vt + A->_m * A->_m;
                    HPDDM::underlying_type<K>* s = values + nu;
                    HPDDM::underlying_type<K>* rwork = s + A->_m;
                    std::fill(a, a + A->_m * dof, K(0.0));
                    for(int i = 0; i < dof; ++i)
                        for(int j = A->_ia[i]; j < A->_ia[i + 1]; ++j)
                            a[i + A->_ja[j] * dof] = A->_a[j];
                    ptA->Type::super::callSolve(a, A->_m);
                    HPDDM::Lapack<K>::gesdd("S", &dof, &(A->_m), a, &dof, s, u, &dof, vt, &(A->_m), work, &lwork, rwork, pos + nu, &info);
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
                if(!std::is_same<K, HPDDM::underlying_type<K>>::value)
                    delete [] values;
                delete [] a;
                ptA->Type::super::initialize(nu);
                if(nbSolver == 100)
                    std::for_each(vecAIJ.begin(), vecAIJ.end(), std::default_delete<const HPDDM::MatrixCSR<K>>());
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
    {"eps", &typeid(HPDDM::underlying_type<K>)},
    {"dim", &typeid(long)},
    {"iter", &typeid(long)},
    {"timing", &typeid(KN<double>*)},
    {"excluded", &typeid(bool)},
    {"ret", &typeid(Pair<K>*)},
    {"O", &typeid(Matrice_Creuse<K>*)},
    {"solver", &typeid(long)}
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
    HPDDM::underlying_type<K> eps = nargs[0] ? GetAny<HPDDM::underlying_type<K>>((*nargs[0])(stack)) : -1.0;
    if(nargs[0])
        std::cerr << "Please do not use the legacy option \"-eps\", set instead \"-hpddm_tol\", cf. \"-hpddm_help\"" << std::endl;
    if(std::abs(eps + 1.0) > 1.0e-6)
        opt["tol"] = eps;
    int dim = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : -1;
    if(nargs[1])
        std::cerr << "Please do not use the legacy option \"-dim\", set instead \"-hpddm_gmres_restart\", cf. \"-hpddm_help\"" << std::endl;
    if(dim != -1)
        opt["gmres_restart"] = dim;
    int iter = nargs[2] ? GetAny<long>((*nargs[2])(stack)) : -1;
    if(nargs[2])
        std::cerr << "Please do not use the legacy option \"-iter\", set instead \"-hpddm_max_it\", cf. \"-hpddm_help\"" << std::endl;
    if(iter != -1)
        opt["max_it"] = iter;
    if(nargs[7])
        std::cerr << "Please do not use the legacy option \"-solver\", set instead \"-hpddm_schwarz_method\" and \"-hpddm_schwarz_coarse_correction\", cf. \"-hpddm_help\"" << std::endl;
    KN<double>* timing = nargs[3] ? GetAny<KN<double>*>((*nargs[3])(stack)) : 0;
    Pair<K>* pair = nargs[5] ? GetAny<Pair<K>*>((*nargs[5])(stack)) : 0;
    if(opt.set("schwarz_coarse_correction") && pair)
        if(pair->p) {
            int flag;
            MPI_Test(&(pair->p->first), &flag, MPI_STATUS_IGNORE);
        }
    MatriceMorse<K>* mA = nargs[6] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[6])(stack))->A)) : 0;
    unsigned short mu = ptX->n / ptA->getDof();
    MPI_Allreduce(MPI_IN_PLACE, &mu, 1, MPI_UNSIGNED_SHORT, MPI_MAX, ptA->getCommunicator());
    const HPDDM::MatrixCSR<K>* A = ptA->getMatrix();
    bool alreadyRenumbered = false;
    if(HPDDM::Wrapper<K>::I == 'F' && SUBDOMAIN<K>::_numbering == 'F' && mu > 1) {
        alreadyRenumbered = true;
        std::for_each(A->_ia, A->_ia + A->_n + 1, [](int& i) { ++i; });
        std::for_each(A->_ja, A->_ja + A->_nnz, [](int& i) { ++i; });
    }
    double timer = MPI_Wtime();
    if(mpisize > 1 && (mA && opt.any_of("schwarz_method", { 1, 2, 4 }))) {
        HPDDM::MatrixCSR<K> dA(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
        ptA->callNumfact(&dA);
    }
    else {
        if(!alreadyRenumbered)
            ptA->callNumfact();
#ifndef SUITESPARSESUB
        else
            ptA->template callNumfact<'F'>();
#endif
    }
    if(timing) (*timing)[1] = MPI_Wtime() - timer;
    if(HPDDM::Wrapper<K>::I == 'F' && SUBDOMAIN<K>::_numbering == 'C' && mu > 1) {
        std::for_each(A->_ia, A->_ia + A->_n + 1, [](int& i) { ++i; });
        std::for_each(A->_ja, A->_ja + A->_nnz, [](int& i) { ++i; });
    }
    bool excluded = nargs[4] ? GetAny<bool>((*nargs[4])(stack)) : false;
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
            pair->destroy();
            pair = nullptr;
            timer = MPI_Wtime();
        }
    MPI_Barrier(MPI_COMM_WORLD);
    if(opt.val<unsigned short>("reuse_preconditioner") <= 1 && !excluded && pair && pair->p && timing && mpisize > 1)
        (*timing)[timing->n - 1] += MPI_Wtime() - timer;
    int rank;
    MPI_Comm_rank(ptA->getCommunicator(), &rank);
    if(rank != 0 || excluded)
        opt.remove("verbosity");
    timer = MPI_Wtime();
    if(!excluded)
        HPDDM::IterativeMethod::solve(*ptA, (K*)*ptRHS, (K*)*ptX, mu, MPI_COMM_WORLD);
    else
        HPDDM::IterativeMethod::solve<true>(*ptA, (K*)nullptr, (K*)nullptr, mu, MPI_COMM_WORLD);
    timer = MPI_Wtime() - timer;
    if(!excluded) {
        if(verbosity > 0 && rank == 0)
            std::cout << std::scientific << " --- system solved (in " << timer << ")" << std::endl;
        HPDDM::underlying_type<K>* storage = new HPDDM::underlying_type<K>[2 * mu];
        ptA->computeError(*ptX, *ptRHS, storage, mu);
        char v = opt.val<char>("verbosity", 0);
        if(v > 0 && rank == 0) {
            std::cout << std::scientific << " --- error = " << storage[1] << " / " << storage[0];
            if(mu > 1)
                std::cout << " (rhs #1)\n";
            else
                std::cout << "\n";
            for(unsigned short nu = (v > 2 ? 1 : v > 1 ? std::max(2, mu - 1) : mu); nu < mu; ++nu)
                std::cout << std::scientific << "             " << storage[2 * nu + 1] << " / " << storage[2 * nu] << " (rhs #" << (nu + 1) << ")\n";
        }
        delete [] storage;
    }
    if(HPDDM::Wrapper<K>::I == 'F' && mu > 1) {
        std::for_each(A->_ja, A->_ja + A->_nnz, [](int& i) { --i; });
        std::for_each(A->_ia, A->_ia + A->_n + 1, [](int& i) { --i; });
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
        changeOperator_Op(const basicAC_F0& args, Expression param1, Expression param2) : A(param1), mat(param2) {
            args.SetNameParam(n_name_param, name_param, nullptr);
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

template<class Type, class K>
class set_Op : public E_F0mps {
    public:
        Expression A;
        static const int n_name_param = 1;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        set_Op(const basicAC_F0& args, Expression param) : A(param) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type set_Op<Type, K>::name_param[] = {
    {"sparams", &typeid(string*)}
};
template<class Type, class K>
class set : public OneOperator {
    public:
        set() : OneOperator(atype<long>(), atype<Type*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new set_Op<Type, K>(args, t[0]->CastTo(args[0]));
        }
};
template<class Type, class K>
AnyType set_Op<Type, K>::operator()(Stack stack) const {
    std::string params = nargs[0] ? *(GetAny<string*>((*nargs[0])(stack))) : "";
    HPDDM::Option::get()->parse(params);
    return 0L;
}

template<class K>
class distributedDot_Op : public E_F0mps {
    public:
        Expression A;
        Expression in;
        Expression out;
        static const int n_name_param = 1;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        distributedDot_Op<K>(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : A(param1), in(param2), out(param3) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class K>
basicAC_F0::name_and_type distributedDot_Op<K>::name_param[] = {
    {"communicator", &typeid(pcommworld)}
};
template<class K>
class distributedDot : public OneOperator {
    public:
        distributedDot() : OneOperator(atype<K>(), atype<KN<double>*>(), atype<KN<K>*>(), atype<KN<K>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new distributedDot_Op<K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]));
        }
};
template<class K, typename std::enable_if<!std::is_same<K, double>::value>::type* = nullptr>
inline K prod(K u, double d, K v) {
    return std::conj(u) * d * v;
}
template<class K, typename std::enable_if<std::is_same<K, double>::value>::type* = nullptr>
inline K prod(K u, double d, K v) {
    return u * d * v;
}
template<class K>
AnyType distributedDot_Op<K>::operator()(Stack stack) const {
    KN<double>* pA = GetAny<KN<double>*>((*A)(stack));
    KN<K>* pin = GetAny<KN<K>*>((*in)(stack));
    KN<K>* pout = GetAny<KN<K>*>((*out)(stack));
    MPI_Comm* comm = nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0;
    K dot = K();
    for(int i = 0; i < pin->n; ++i)
        dot += prod(pin->operator[](i), pA->operator[](i), pout->operator[](i));
    MPI_Allreduce(MPI_IN_PLACE, &dot, 1, HPDDM::Wrapper<K>::mpi_type(), MPI_SUM, comm ? *((MPI_Comm*)comm) : MPI_COMM_WORLD);
    return SetAny<K>(dot);
}

template<class Type, class K>
class distributedMV_Op : public E_F0mps {
    public:
        Expression A;
        Expression Mat;
        Expression in;
        Expression out;
        static const int n_name_param = 0;
        static basicAC_F0::name_and_type name_param[];
        distributedMV_Op<Type, K>(const basicAC_F0& args, Expression param1, Expression param2, Expression param3, Expression param4) : A(param1), Mat(param2), in(param3), out(param4) {
            args.SetNameParam(n_name_param, name_param, nullptr);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type distributedMV_Op<Type, K>::name_param[] = { };
template<class Type, class K>
class distributedMV : public OneOperator {
    public:
        distributedMV() : OneOperator(atype<long>(), atype<Type*>(), atype<Matrice_Creuse<K>*>(), atype<KN<K>*>(), atype<KN<K>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new distributedMV_Op<Type, K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
        }
};
template<class Type, class K>
AnyType distributedMV_Op<Type, K>::operator()(Stack stack) const {
    Type* pA = GetAny<Type*>((*A)(stack));
    KN<K>* pin = GetAny<KN<K>*>((*in)(stack));
    KN<K>* pout = GetAny<KN<K>*>((*out)(stack));
    pout->resize(pin->n);
    unsigned short mu = pin->n / pA->getDof();
    MatriceMorse<K>* mA = static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*Mat)(stack))->A));
    HPDDM::MatrixCSR<K> dA(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
    bool allocate = pA->setBuffer();
    pA->GMV((K*)*pin, (K*)*pout, mu, &dA);
    pA->clearBuffer(allocate);
    return 0L;
}

template<class Type, class K>
class scaledExchange_Op : public E_F0mps {
    public:
        Expression A;
        Expression in;
        static const int n_name_param = 0;
        static basicAC_F0::name_and_type name_param[];
        scaledExchange_Op<Type, K>(const basicAC_F0& args, Expression param1, Expression param2) : A(param1), in(param2) {
            args.SetNameParam(n_name_param, name_param, nullptr);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type scaledExchange_Op<Type, K>::name_param[] = { };
template<class Type, class K>
class scaledExchange : public OneOperator {
    public:
        scaledExchange() : OneOperator(atype<long>(), atype<Type*>(), atype<KN<K>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new scaledExchange_Op<Type, K>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
        }
};
template<class Type, class K>
AnyType scaledExchange_Op<Type, K>::operator()(Stack stack) const {
    Type* pA = GetAny<Type*>((*A)(stack));
    KN<K>* pin = GetAny<KN<K>*>((*in)(stack));
    unsigned short mu = pin->n / pA->getDof();
    pA->Type::template scaledExchange<true>((K*)*pin, mu);
    return 0L;
}

template<class T, class U, class K>
class ProdSchwarz {
    public:
        const T t;
        const U u;
        ProdSchwarz(T v, U w) : t(v), u(w) {}
        void prod(U x) const { bool allocate = t->setBuffer(); t->GMV(*(this->u), *x); t->clearBuffer(allocate); };
        static U mv(U Ax, ProdSchwarz<T, U, K> A) {
            A.prod(Ax);
            return Ax;
        }
};

template<class T, class U, class K>
class InvSchwarz {
    public:
        const T t;
        const U u;
        InvSchwarz(T v, U w) : t(v), u(w) {}
        void solve(U out) const {
            if(out->n != u->n || u->n < (*t).getDof())
                return;
            HPDDM::Option& opt = *HPDDM::Option::get();
            unsigned short mu = u->n / (*t).getDof();
            MPI_Allreduce(MPI_IN_PLACE, &mu, 1, MPI_UNSIGNED_SHORT, MPI_MAX, (*t).getCommunicator());
            const HPDDM::MatrixCSR<K>* A = (*t).getMatrix();
            bool alreadyRenumbered = false;
            if(HPDDM::Wrapper<K>::I == 'F' && SUBDOMAIN<K>::_numbering == 'F' && mu > 1) {
                alreadyRenumbered = true;
                std::for_each(A->_ia, A->_ia + A->_n + 1, [](int& i) { ++i; });
                std::for_each(A->_ja, A->_ja + A->_nnz, [](int& i) { ++i; });
            }
            if(!alreadyRenumbered)
                (*t).callNumfact();
#ifndef SUITESPARSESUB
            else
                (*t).template callNumfact<'F'>();
#endif
            if(HPDDM::Wrapper<K>::I == 'F' && SUBDOMAIN<K>::_numbering == 'C' && mu > 1) {
                std::for_each(A->_ia, A->_ia + A->_n + 1, [](int& i) { ++i; });
                std::for_each(A->_ja, A->_ja + A->_nnz, [](int& i) { ++i; });
            }
            if(mpirank != 0)
                opt.remove("verbosity");
            HPDDM::IterativeMethod::solve(*t, (K*)*u, (K*)*out, mu, MPI_COMM_WORLD);
        }
        static U init(U Ax, InvSchwarz<T, U, K> A) {
            A.solve(Ax);
            return Ax;
        }
};

template<template<class, char> class Type, class K, char S>
void add() {
    Dcl_Type<Type<K, S>*>(Initialize<Type<K, S>>, Delete<Type<K, S>>);

    TheOperators->Add("<-", new initDDM<Type<K, S>, K>);
    Global.Add("attachCoarseOperator", "(", new attachCoarseOperator<Type<K, S>, K>);
    Global.Add("DDM", "(", new solveDDM<Type<K, S>, K>);
    Global.Add("changeOperator", "(", new changeOperator<Type<K, S>, K>);
    Global.Add("set", "(", new set<Type<K, S>, K>);
    addProd<Type<K, S>, ProdSchwarz, KN<K>, K>();
    addInv<Type<K, S>, InvSchwarz, KN<K>, K>();
    Global.Add("dscalprod", "(", new distributedDot<K>);
    Global.Add("dmv", "(", new distributedMV<Type<K, S>, K>);
    Global.Add("scaledExchange", "(", new scaledExchange<Type<K, S>, K>);
    Global.Add("destroyRecycling", "(", new OneOperator1_<bool, Type<K, S>*>(Schwarz::destroyRecycling<Type<K, S>, K>));
}
}

static void Init_Schwarz() {
    Global.Add("getOption", "(", new OneOperator1_<double, string*>(Schwarz::getOpt));
    Global.Add("isSetOption", "(", new OneOperator1_<bool, string*>(Schwarz::isSetOpt));
#if defined(DSUITESPARSE) || defined(DHYPRE)
    const char ds = 'G';
#else
    const char ds = 'S';
#endif
    const char zs = 'G';
    Schwarz::add<HpSchwarz, double, ds>();
    zzzfff->Add("dschwarz", atype<HpSchwarz<double, ds>*>());
#ifndef DHYPRE
    // Schwarz::add<HpSchwarz, float, ds>();
    // zzzfff->Add("sschwarz", atype<HpSchwarz<float, ds>*>());
    Schwarz::add<HpSchwarz, std::complex<double>, zs>();
    zzzfff->Add("zschwarz", atype<HpSchwarz<std::complex<double>, zs>*>());
    // Schwarz::add<HpSchwarz, std::complex<float>, zs>();
    // zzzfff->Add("cschwarz", atype<HpSchwarz<std::complex<float>, zs>*>());
#endif
    // Dcl_Type<Pair<float>*>(InitP<Pair<float>>, Destroy<Pair<float>>);
    // zzzfff->Add("spair", atype<Pair<double>*>());
    Dcl_Type<Pair<double>*>(InitP<Pair<double>>, Destroy<Pair<double>>);
    zzzfff->Add("dpair", atype<Pair<double>*>());
    // Dcl_Type<Pair<std::complex<float>>*>(InitP<Pair<std::complex<float>>>, Destroy<Pair<std::complex<float>>>);
    // zzzfff->Add("cpair", atype<Pair<std::complex<float>>*>());
    Dcl_Type<Pair<std::complex<double>>*>(InitP<Pair<std::complex<double>>>, Destroy<Pair<std::complex<double>>>);
    zzzfff->Add("zpair", atype<Pair<std::complex<double>>*>());
}

LOADFUNC(Init_Schwarz)
