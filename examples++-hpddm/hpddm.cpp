//ff-c++-LIBRARY-dep: cxx11 hpddm [mumps parmetis ptscotch scotch scalapack|umfpack] [mkl|blas] mpi pthread mpifc fc
//ff-c++-cpp-dep:

#define HPDDM_SCHWARZ 1
#define HPDDM_FETI    0
#define HPDDM_BDD     0

#include "common.hpp"

namespace Schwarz {
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
    Type* ptA = GetAny<Type*>((*A)(stack));
    Matrice_Creuse<K>* pA = GetAny<Matrice_Creuse<K>*>((*Mat)(stack));
    MatriceMorse<K>* mA = pA->A ? static_cast<MatriceMorse<K>*>(&(*pA->A)) : nullptr;
    KN<long>* ptO = GetAny<KN<long>*>((*o)(stack));
    KN<KN<long>>* ptR = GetAny<KN<KN<long>>*>((*R)(stack));
    if(ptO)
        ptA->HPDDM::template Subdomain<K>::initialize(mA ? new HPDDM::MatrixCSR<K>(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique) : 0, STL<long>(*ptO), *ptR, nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0);
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
    KN<double>* timing = nargs[4] ? GetAny<KN<double>*>((*nargs[4])(stack)) : 0;
    std::pair<MPI_Request, const K*>* ret = nullptr;
    if(mA) {
        long nbSolver = 0;
        std::vector<const HPDDM::MatrixCSR<K>*> vecAIJ;
        if(mA) {
            HPDDM::MatrixCSR<K> dA(mA->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
            MatriceMorse<K>* mB = nargs[1] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[1])(stack))->A)) : nullptr;
            MatriceMorse<K>* mP = nargs[2] && opt.any_of("schwarz_method", { 1, 2, 4 }) ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[2])(stack))->A)) : nullptr;
            if(dA._n == dA._m) {
                if(timing) { // tic
                    timing->resize(timing->n + 1);
                    (*timing)[timing->n - 1] = MPI_Wtime();
                }
                const HPDDM::MatrixCSR<K>* const dP = mP ? new HPDDM::MatrixCSR<K>(mP->n, mP->m, mP->nbcoef, mP->a, mP->lg, mP->cl, mP->symetrique) : nullptr;
                if(mB) {
                    HPDDM::MatrixCSR<K> dB(mB->n, mB->m, mB->nbcoef, mB->a, mB->lg, mB->cl, mB->symetrique);
                    ptA->template solveGEVP<EIGENSOLVER>(&dA, &dB, dP);
                }
                else
                    ptA->template solveGEVP<EIGENSOLVER>(&dA, nullptr, dP);
                mA->nbcoef = dA._nnz;
                mA->a = dA._a;
                mA->lg = dA._ia;
                mA->cl = dA._ja;
                delete dP;
                if(timing) { // toc
                    (*timing)[timing->n - 1] = MPI_Wtime() - (*timing)[timing->n - 1];
                }
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
                unsigned short nu = std::min(opt.template val<unsigned short>("geneo_nu", 20), static_cast<unsigned short>(first._m));
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
                opt["geneo_nu"] = nu;
            }
            else
                ptA->Type::super::initialize(0);
        }
        MPI_Barrier(comm);
        if(timing) { // tic
            timing->resize(timing->n + 1);
            (*timing)[timing->n - 1] = MPI_Wtime();
        }
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
        if(timing) { // toc
            (*timing)[timing->n - 1] = MPI_Wtime() - (*timing)[timing->n - 1];
        }
    }
    else {
        MPI_Barrier(comm);
        if(!ptA->getVectors())
            ret = ptA->template buildTwo<2>(comm);
        else if(ptA->exclusion(comm))
            ret = ptA->template buildTwo<1>(comm);
        else
            ret = ptA->template buildTwo<0>(comm);
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
        static const int n_name_param = 9;
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
    {"solver", &typeid(long)},
    {"communicator", &typeid(pcommworld)}
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
    const std::string& prefix = ptA->prefix();
    HPDDM::underlying_type<K> eps = nargs[0] ? GetAny<HPDDM::underlying_type<K>>((*nargs[0])(stack)) : -1.0;
    if(nargs[0])
        std::cerr << "Please do not use the legacy option \"-eps\", set instead \"-hpddm_tol\", cf. \"-hpddm_help\"" << std::endl;
    if(std::abs(eps + 1.0) > 1.0e-6)
        opt[prefix + "tol"] = eps;
    int dim = nargs[1] ? GetAny<long>((*nargs[1])(stack)) : -1;
    if(nargs[1])
        std::cerr << "Please do not use the legacy option \"-dim\", set instead \"-hpddm_gmres_restart\", cf. \"-hpddm_help\"" << std::endl;
    if(dim != -1)
        opt[prefix + "gmres_restart"] = dim;
    int iter = nargs[2] ? GetAny<long>((*nargs[2])(stack)) : -1;
    if(nargs[2])
        std::cerr << "Please do not use the legacy option \"-iter\", set instead \"-hpddm_max_it\", cf. \"-hpddm_help\"" << std::endl;
    if(iter != -1)
        opt[prefix + "max_it"] = iter;
    if(nargs[7])
        std::cerr << "Please do not use the legacy option \"-solver\", set instead \"-hpddm_schwarz_method\" and \"-hpddm_schwarz_coarse_correction\", cf. \"-hpddm_help\"" << std::endl;
    KN<double>* timing = nargs[3] ? GetAny<KN<double>*>((*nargs[3])(stack)) : 0;
    Pair<K>* pair = nargs[5] ? GetAny<Pair<K>*>((*nargs[5])(stack)) : 0;
    if(opt.set(prefix + "schwarz_coarse_correction") && pair)
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
    if(timing) { // tic
        timing->resize(timing->n + 1);
        (*timing)[timing->n - 1] = MPI_Wtime();
    }
    if(mpisize > 1 && (mA && opt.any_of(prefix + "schwarz_method", { 1, 2, 4 }))) {
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
    if(timing) { // toc
        (*timing)[timing->n - 1] = MPI_Wtime() - (*timing)[timing->n - 1];
    }
    if(HPDDM::Wrapper<K>::I == 'F' && SUBDOMAIN<K>::_numbering == 'C' && mu > 1) {
        std::for_each(A->_ia, A->_ia + A->_n + 1, [](int& i) { ++i; });
        std::for_each(A->_ja, A->_ja + A->_nnz, [](int& i) { ++i; });
    }
    bool excluded = nargs[4] && GetAny<bool>((*nargs[4])(stack));
    if(excluded)
        opt[prefix + "master_exclude"];
    if(pair)
        if(pair->p) {
            MPI_Wait(&(pair->p->first), MPI_STATUS_IGNORE);
            delete [] pair->p->second;
            pair->destroy();
            pair = nullptr;
        }
    MPI_Comm comm = nargs[8] ? *(MPI_Comm*)GetAny<pcommworld>((*nargs[8])(stack)) : MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(ptA->getCommunicator(), &rank);
    if(rank != mpirank || rank != 0) {
        opt.remove("verbosity");
        if(prefix.size() > 0)
            opt.remove(prefix + "verbosity");
    }
    MPI_Barrier(comm);
    double timer = MPI_Wtime();
    if(timing) { // tic
        timing->resize(timing->n + 1);
        (*timing)[timing->n - 1] = timer;
    }
    if(!excluded)
        HPDDM::IterativeMethod::solve(*ptA, (K*)*ptRHS, (K*)*ptX, mu, comm);
    else
        HPDDM::IterativeMethod::solve<true>(*ptA, (K*)nullptr, (K*)nullptr, mu, comm);
    timer = MPI_Wtime() - timer;
    if(timing) { // toc
        (*timing)[timing->n - 1] = timer;
    }
    if(!excluded && verbosity > 0 && rank == 0)
        std::cout << std::scientific << " --- system solved (in " << timer << ")" << std::endl;
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
    if(nargs[0])
        HPDDM::Option::get()->parse(*(GetAny<string*>((*nargs[0])(stack))));
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

template<class T, class U, class K, char N>
class ProdSchwarz {
    public:
        const T t;
        const U u;
        ProdSchwarz(T v, U w) : t(v), u(w) {}
        void prod(U x) const { bool allocate = t->setBuffer(); t->GMV(*(this->u), *x); t->clearBuffer(allocate); };
        static U mv(U Ax, ProdSchwarz<T, U, K, N> A) {
            A.prod(Ax);
            return Ax;
        }
        static U init(U Ax, ProdSchwarz<T, U, K, N> A) {
            Ax->init(A.u->n);
            return mv(Ax, A);
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
                HPDDM::Option::get()->remove((*t).prefix("verbosity"));
            HPDDM::IterativeMethod::solve(*t, (K*)*u, (K*)*out, mu, MPI_COMM_WORLD);
        }
        static U init(U Ax, InvSchwarz<T, U, K> A) {
            A.solve(Ax);
            return Ax;
        }
};

template<class R>
class IterativeMethod : public OneOperator {
    public:
        typedef KN<R> Kn;
        typedef KN_<R> Kn_;
        class MatF_O : VirtualMatrice<R> {
            public:
                Stack stack;
                mutable Kn x;
                C_F0 c_x;
                Expression mat1, mat;
                typedef typename VirtualMatrice<R>::plusAx plusAx;
                MatF_O(int n, Stack stk, const OneOperator* op) :
                    VirtualMatrice<R>(n), stack(stk), x(n), c_x(CPValue(x)),
                    mat1(op ? op->code(basicAC_F0_wa(c_x)) : 0),
                    mat(op ? CastTo<Kn_>(C_F0(mat1, (aType)*op)) : 0) { }
                ~MatF_O() {
                    if(mat1 != mat)
                        delete mat;
                    delete mat1;
                    Expression zzz = c_x;
                    delete zzz;
                }
                void addMatMul(const Kn_& xx, Kn_& Ax) const {
                    ffassert(xx.N() == Ax.N());
                    x = xx;
                    Ax += GetAny<Kn_>((*mat)(stack));
                    WhereStackOfPtr2Free(stack)->clean();
                }
                void mv(const R* const in, const int& n, R* const out) const {
                    KN_<R> xx((R*)in, n);
                    KN_<R> yy(out, n);
                    yy = R();
                    yy = plusAx(this, xx);
                }
                bool ChecknbLine(int) const { return true; }
                bool ChecknbColumn(int) const { return true; }
        };
        class Operator : public HPDDM::EmptyOperator<R> {
            public:
                MatF_O& mat;
                MatF_O& prec;
                Operator(MatF_O& m, MatF_O& p) : mat(m), prec(p), HPDDM::EmptyOperator<R>(m.x.N()) { }
                void GMV(const R* const in, R* const out, const int& mu = 1) const {
                    mat.mv(in, HPDDM::EmptyOperator<R>::_n, out);
                }
                template<bool = true>
                void apply(const R* const in, R* const out, const unsigned short& mu = 1, R* = nullptr, const unsigned short& = 0) const {
                    if(prec.mat)
                        prec.mv(in, HPDDM::EmptyOperator<R>::_n, out);
                    else
                        std::copy_n(in, HPDDM::EmptyOperator<R>::_n, out);
                }
        };
        class E_LCG : public E_F0mps {
            public:
                static const int n_name_param = 4;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                const OneOperator *A, *C;
                Expression X, B;
                E_LCG(const basicAC_F0& args) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    { const Polymorphic* op = dynamic_cast<const Polymorphic*>(args[0].LeftValue());
                        ffassert(op);
                        A = op->Find("(", ArrayOfaType(atype<Kn*>(), false)); }
                    if(nargs[0]) {
                        const Polymorphic* op = dynamic_cast<const Polymorphic*>(nargs[0]);
                        ffassert(op);
                        C = op->Find("(", ArrayOfaType(atype<Kn*>(), false));
                    }
                    else
                        C = 0;
                    X = to<Kn*>(args[1]);
                    B = to<Kn*>(args[2]);
                }
                virtual AnyType operator()(Stack stack)  const {
                    int ret = -1;
                    try {
                        Kn& x = *GetAny<Kn*>((*X)(stack));
                        int n = x.N();
                        MPI_Comm comm = nargs[3] ? *(MPI_Comm*)GetAny<pcommworld>((*nargs[3])(stack)) : MPI_COMM_WORLD;
                        Kn& b = *GetAny<Kn*>((*B)(stack));
                        MatF_O AA(n, stack, A);
                        MatF_O PP(n, stack, C);
                        Operator Op(AA, PP);
                        if(nargs[1])
                            Op.setPrefix(*(GetAny<string*>((*nargs[1])(stack))));
                        if(nargs[2])
                            HPDDM::Option::get()->parse(*(GetAny<string*>((*nargs[2])(stack))));
                        ret = HPDDM::IterativeMethod::solve(Op, (R*)b, (R*)x, 1, comm);
                    }
                    catch(...) {
                        throw;
                    }
                    return SetAny<long>(ret);
                }
                operator aType() const { return atype<long>(); }
        };
        E_F0* code(const basicAC_F0& args) const { return new E_LCG(args); }
        IterativeMethod() : OneOperator(atype<long>(), atype<Polymorphic*>(), atype<KN<R>*>(), atype<KN<R>*>()) { }
};

template<class R>
basicAC_F0::name_and_type IterativeMethod<R>::E_LCG::name_param[] = {
    {"precon", &typeid(Polymorphic*)},
    {"prefix", &typeid(string*)},
    {"sparams", &typeid(string*)},
    {"comm", &typeid(pcommworld)}
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
    Global.Add("destroyRecycling", "(", new OneOperator1_<bool, Type<K, S>*>(destroyRecycling<Type<K, S>, K>));
    Global.Add("statistics", "(", new OneOperator1_<bool, Type<K, S>*>(statistics<Type<K, S>>));
    Global.Add("exchange", "(", new exchangeIn<Type<K, S>, K>);
    Global.Add("exchange", "(", new exchangeInOut<Type<K, S>, K>);
    Global.Add("IterativeMethod","(",new IterativeMethod<K>());
}
}

static void Init_Schwarz() {
    Init_Common();
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
    aType t;
    int r;
    if(!zzzfff->InMotClef("dpair", t, r)) {
        // Dcl_Type<Pair<float>*>(InitP<Pair<float>>, Destroy<Pair<float>>);
        // zzzfff->Add("spair", atype<Pair<double>*>());
        Dcl_Type<Pair<double>*>(InitP<Pair<double>>, Destroy<Pair<double>>);
        zzzfff->Add("dpair", atype<Pair<double>*>());
        // Dcl_Type<Pair<std::complex<float>>*>(InitP<Pair<std::complex<float>>>, Destroy<Pair<std::complex<float>>>);
        // zzzfff->Add("cpair", atype<Pair<std::complex<float>>*>());
        Dcl_Type<Pair<std::complex<double>>*>(InitP<Pair<std::complex<double>>>, Destroy<Pair<std::complex<double>>>);
        zzzfff->Add("zpair", atype<Pair<std::complex<double>>*>());
    }
}

LOADFUNC(Init_Schwarz)
