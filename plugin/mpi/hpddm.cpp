//ff-c++-LIBRARY-dep: cxx11 hpddm [mumps parmetis metis ptscotch scotch scalapack|umfpack] [mkl|blas] mpi pthread mpifc fc
//ff-c++-cpp-dep:

#define HPDDM_SCHWARZ                   1
#define HPDDM_FETI                      0
#define HPDDM_BDD                       0
#define HPDDM_INEXACT_COARSE_OPERATOR   1

#include "common_hpddm.hpp"

namespace Schwarz {
template<class Type, class K>
class initDDM : public OneOperator {
    public:
        const int c;
        class E_initDDM : public E_F0mps {
            public:
                Expression A;
                Expression Mat;
                Expression R;
                Expression D;
                const int c;
                static const int n_name_param = 3;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                E_initDDM(const basicAC_F0& args, int d) : A(0), Mat(0), R(0), D(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    A = to<Type*>(args[0]);
                    if(c == 0 || c == 2)
                        Mat = to<Matrice_Creuse<K>*>(args[1]);
                    else
                        Mat = to<long>(args[1]);
                    if(c == 0 || c == 1) {
                        R = to<KN<KN<long>>*>(args[2]);
                        D = to<KN<HPDDM::underlying_type<K>>*>(args[3]);
                    }
                }

                AnyType operator()(Stack stack) const;
                operator aType() const { return atype<Type*>(); }
        };
        E_F0* code(const basicAC_F0 & args) const { return new E_initDDM(args, c); }
        initDDM() : OneOperator(atype<Type*>(), atype<Type*>(), atype<Matrice_Creuse<K>*>(), atype<KN<KN<long>>*>(), atype<KN<HPDDM::underlying_type<K>>*>()), c(0) { }
        initDDM(int) : OneOperator(atype<Type*>(), atype<Type*>(), atype<long>(), atype<KN<KN<long>>*>(), atype<KN<HPDDM::underlying_type<K>>*>()), c(1) { }
        initDDM(int, int) : OneOperator(atype<Type*>(), atype<Type*>(), atype<Matrice_Creuse<K>*>()), c(2) { }
        initDDM(int, int, int) : OneOperator(atype<Type*>(), atype<Type*>(), atype<long>()), c(3) { }
};
template<class Type, class K>
basicAC_F0::name_and_type initDDM<Type, K>::E_initDDM::name_param[] = {
    {"communicator", &typeid(pcommworld)},
    {"scaled", &typeid(bool)},
    {"level", &typeid(long)}
};
template<class Type, class K>
AnyType initDDM<Type, K>::E_initDDM::operator()(Stack stack) const {
    Type* ptA = GetAny<Type*>((*A)(stack));
    HPDDM::MatrixCSR<K>* dA;
    if(c == 0 || c == 2) {
        MatriceMorse<K>* mA = static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*Mat)(stack))->A));
        dA = new_HPDDM_MatrixCSR<K>(mA);//->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
    }
    else {
        int dof = GetAny<long>((*Mat)(stack));
        dA = new HPDDM::MatrixCSR<K>(dof, dof, 0, nullptr, nullptr, nullptr, false);
    }
    if(c == 0 || c == 1) {
        KN<KN<long>>* ptR = GetAny<KN<KN<long>>*>((*R)(stack));
        int level = nargs[2] ? std::abs(GetAny< long >((*nargs[2])(stack))) : 0;
        KN<HPDDM::underlying_type<K>>* ptD = GetAny<KN<HPDDM::underlying_type<K>>*>((*D)(stack));
        if(ptR) {
            KN_<KN<long>> sub(ptR->n > 0 && ptR->operator[](0).n > 0 ? (*ptR)(FromTo(1 + level * ptR->operator[](0).n, 1 + (level + 1) * ptR->operator[](0).n - 1)) : KN<KN<long>>());
            ptA->HPDDM::template Subdomain<K>::initialize(dA, STL<long>(ptR->n > 0 ? ptR->operator[](0) : KN<long>()), sub, nargs[0] ? (MPI_Comm*)GetAny<pcommworld>((*nargs[0])(stack)) : 0);
        }
        if(ptD)
            ptA->initialize(*ptD);
        else
            std::cerr << "Something is really wrong here!" << std::endl;
        if(c == 0 && (!nargs[1] || GetAny<bool>((*nargs[1])(stack))))
            ptA->exchange();
    }
    else {
        const MPI_Comm& comm = MPI_COMM_SELF;
        ptA->HPDDM::template Subdomain<K>::initialize(dA, STL<long>( KN<long>()), KN<KN<long>>(), const_cast<MPI_Comm*>(&comm));
    }
    return ptA;
}

template<class Type, class K>
class attachCoarseOperator : public OneOperator {
    public:
        typedef KN<K> Kn;
        typedef KN_<K> Kn_;
        class MatF_O : public RNM_VirtualMatrix<K>, public Type::super::CoarseCorrection {
            public:
                typedef typename Type::super::CoarseCorrection super;
                Stack stack;
                mutable Kn x;
                C_F0 c_x;
                Expression mat;
                typedef typename RNM_VirtualMatrix<K>::plusAx plusAx;
                MatF_O(int n, Stack stk, const OneOperator* op) :
                    RNM_VirtualMatrix<K>(n), stack(stk), x(n), c_x(CPValue(x)),
                    mat(op ? CastTo<Kn_>(C_F0(op->code(basicAC_F0_wa(c_x)), (aType)*op)) : 0) { }
                ~MatF_O() {
                    delete mat;
                    Expression zzz = c_x;
                    delete zzz;
                }
                virtual void operator()(const K* const in, K* const out) {
                    KN_<K> xx(const_cast<K*>(in), this->N);
                    KN_<K> yy(out, this->N);
                    addMatMul(xx, yy);
                }
                void addMatMul(const Kn_& xx, Kn_& Ax) const {
                    ffassert(xx.N() == this->N && Ax.N() == this->M);
                    x = xx;
                    Ax = GetAny<Kn_>((*mat)(stack));
                    WhereStackOfPtr2Free(stack)->clean();
                }
                plusAx operator*(const Kn& x) const { return plusAx(this, x); }
                bool ChecknbLine(int) const { return true; }
                bool ChecknbColumn(int) const { return true; }
        };
        class MatMatF_O : public RNM_VirtualMatrix<K>, public Type::super::CoarseCorrection {
            public:
                typedef typename Type::super::CoarseCorrection super;
                Stack stack;
                mutable Kn x;
                C_F0 c_x;
                mutable long mu;
                C_F0 c_mu;
                Expression mat;
                typedef typename RNM_VirtualMatrix<K>::plusAx plusAx;
                MatMatF_O(int n, Stack stk, const OneOperator* op) :
                    RNM_VirtualMatrix<K>(n), stack(stk), x(0), c_x(CPValue(x)), mu(1), c_mu(CPValue(mu)),
                    mat(op ? CastTo<Kn_>(C_F0(op->code(basicAC_F0_wa({ c_x, c_mu })), (aType)*op)) : 0) { }
                ~MatMatF_O() {
                    delete mat;
                    Expression zzz = c_x;
                    delete zzz;
                    zzz = c_mu;
                    delete zzz;
                }
                virtual void operator()(const K* const in, K* const out) {
                    KN_<K> xx(const_cast<K*>(in), this->N);
                    KN_<K> yy(out, this->N);
                    addMatMul(xx, yy);
                }
                virtual void operator()(const K* const in, K* const out, int n, unsigned short nu) {
                    mu = nu;
                    KN_<K> xx(const_cast<K*>(in), this->N * mu);
                    KN_<K> yy(out, this->N * mu);
                    addMatMul(xx, yy);
                }
                void addMatMul(const Kn_& xx, Kn_& Ax) const {
                    ffassert(xx.N() == this->N * mu && Ax.N() == this->M * mu);
                    K* backup = x;
                    x.set(xx, this->N * mu);
                    Ax = GetAny<Kn_>((*mat)(stack));
                    x.set(backup, 0);
                    WhereStackOfPtr2Free(stack)->clean();
                }
                plusAx operator*(const Kn& x) const { return plusAx(this, x); }
                bool ChecknbLine(int) const { return true; }
                bool ChecknbColumn(int) const { return true; }
        };
        const int c;
        class E_attachCoarseOperator : public E_F0mps {
            public:
                Expression A;
                Expression comm;
                const OneOperator *codeC;
                const OneOperator *codeMatC;
                const int c;
                static const int n_name_param = 7;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                E_attachCoarseOperator(const basicAC_F0& args, int d) : A(0), comm(0), codeC(0), codeMatC(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    comm = to<pcommworld>(args[0]);
                    A = to<Type*>(args[1]);
                    if(c == 1) {
                        const Polymorphic* op = dynamic_cast<const Polymorphic*>(args[2].LeftValue());
                        ffassert(op);
                        codeMatC = op->Find("(", ArrayOfaType(atype<KN<K>*>(), atype<long>(), false));
                        if(!codeMatC)
                            codeC = op->Find("(", ArrayOfaType(atype<KN<K>*>(), false));
                    }
                }

                AnyType operator()(Stack stack) const;
                operator aType() const { return atype<long>(); }
        };
        E_F0* code(const basicAC_F0 & args) const { return new E_attachCoarseOperator(args, c); }
        attachCoarseOperator() : OneOperator(atype<long>(), atype<pcommworld>(), atype<Type*>()), c(0) { }
        attachCoarseOperator(int) : OneOperator(atype<long>(), atype<pcommworld>(), atype<Type*>(), atype<Polymorphic*>()), c(1) { }
};
template<class Type, class K>
basicAC_F0::name_and_type attachCoarseOperator<Type, K>::E_attachCoarseOperator::name_param[] = {
    {"A", &typeid(Matrice_Creuse<K>*)},
    {"B", &typeid(Matrice_Creuse<K>*)},
    {"pattern", &typeid(Matrice_Creuse<K>*)},
    {"threshold", &typeid(HPDDM::underlying_type<K>)},
    {"timing", &typeid(KN<double>*)},
    {"ret", &typeid(Pair<K>*)},
    {"deflation", &typeid(FEbaseArrayKn<K>*)}
};
template<class Type, class K>
AnyType attachCoarseOperator<Type, K>::E_attachCoarseOperator::operator()(Stack stack) const {
    pcommworld ptComm = GetAny<pcommworld>((*comm)(stack));
    MPI_Comm comm = *(MPI_Comm*)ptComm;
    Type* ptA = GetAny<Type*>((*A)(stack));
    if(ptA->_cc) {
        delete ptA->_cc;
        ptA->_cc = nullptr;
    }
    if(c == 0) {
        MatriceMorse<K>* mA = nargs[0] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[0])(stack))->A)) : 0;
        Pair<K>* pair = nargs[5] ? GetAny<Pair<K>*>((*nargs[5])(stack)) : 0;
        FEbaseArrayKn<K>* deflation = nargs[6] ? GetAny<FEbaseArrayKn<K>*>((*nargs[6])(stack)) : 0;
        HPDDM::Option& opt = *HPDDM::Option::get();
        KN<double>* timing = nargs[4] ? GetAny<KN<double>*>((*nargs[4])(stack)) : 0;
        std::pair<MPI_Request, const K*>* ret = nullptr;
        if(mA) {
            ff_HPDDM_MatrixCSR<K> dA(mA);//->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
            MatriceMorse<K>* mB = nargs[1] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[1])(stack))->A)) : nullptr;
            MatriceMorse<K>* mP = nargs[2] && opt.any_of("schwarz_method", { 1, 2, 4 }) ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[2])(stack))->A)) : nullptr;
            if(dA._n == dA._m && !deflation) {
                if(timing) { // tic
                    timing->resize(timing->n + 1);
                    (*timing)[timing->n - 1] = MPI_Wtime();
                }
                const HPDDM::MatrixCSR<K>* const dP = new_HPDDM_MatrixCSR<K>(mP);
                //mP ? new HPDDM::MatrixCSR<K>(mP->n, mP->m, mP->nbcoef, mP->a, mP->lg, mP->cl, mP->symetrique) : nullptr;
                if(mB) {
       //                 HPDDM::MatrixCSR<K> dB(mB->n, mB->m, mB->nbcoef, mB->a, mB->lg, mB->cl, mB->symetrique);
                    ff_HPDDM_MatrixCSR<K> dB(mB);
                    ptA->template solveGEVP<EIGENSOLVER>(&dA, &dB, dP);
                }
                else
                    ptA->template solveGEVP<EIGENSOLVER>(&dA, nullptr, dP);
                set_ff_matrix(mA,dA);
                delete dP;
                if(timing) { // toc
                    (*timing)[timing->n - 1] = MPI_Wtime() - (*timing)[timing->n - 1];
                }
            }
            else if(deflation && deflation->N > 0 && !ptA->getVectors()) {
                K** ev = new K*[deflation->N];
                *ev = new K[deflation->N * deflation->get(0)->n];
                for(int i = 0; i < deflation->N; ++i) {
                    ev[i] = *ev + i * deflation->get(0)->n;
                    std::copy_n(&(*deflation->get(i))[0], deflation->get(i)->n, ev[i]);
                }
                ptA->setVectors(ev);
                ptA->Type::super::initialize(deflation->N);
            }
            if(timing) { // tic
                MPI_Barrier(comm);
                timing->resize(timing->n + 1);
                (*timing)[timing->n - 1] = MPI_Wtime();
            }
            if(ptA->exclusion(comm)) {
                if(pair)
                    pair->p = ptA->template buildTwo<1>(comm, &dA);
                else
                    ret = ptA->template buildTwo<1>(comm, &dA);
            }
            else {
                if(pair)
                    pair->p = ptA->template buildTwo<0>(comm, &dA);
                else
                    ret = ptA->template buildTwo<0>(comm, &dA);
            }
            if(timing) { // toc
                (*timing)[timing->n - 1] = MPI_Wtime() - (*timing)[timing->n - 1];
            }
        }
        else {
            if(timing)
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
    else {
        if(codeMatC)
            ptA->_cc = new attachCoarseOperator<Type, K>::MatMatF_O(ptA->getDof(), stack, codeMatC);
        else if(codeC)
            ptA->_cc = new attachCoarseOperator<Type, K>::MatF_O(ptA->getDof(), stack, codeC);
        else
            ffassert(0);
        return 0L;
    }
}

template<class Type, class K>
class solveDDM_Op : public E_F0mps {
    public:
        Expression A;
        Expression rhs;
        Expression x;
        static const int n_name_param = 5;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        solveDDM_Op(const basicAC_F0& args, Expression param1, Expression param2, Expression param3) : A(param1), rhs(param2), x(param3) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type solveDDM_Op<Type, K>::name_param[] = {
    {"timing", &typeid(KN<double>*)},
    {"excluded", &typeid(bool)},
    {"ret", &typeid(Pair<K>*)},
    {"O", &typeid(Matrice_Creuse<K>*)},
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
    KN<K>* ptRHS = GetAny<KN<K>*>((*rhs)(stack));
    KN<K>* ptX = GetAny<KN<K>*>((*x)(stack));
    Type* ptA = GetAny<Type*>((*A)(stack));
    if(ptX->n != ptRHS->n || ptRHS->n < ptA->getDof())
        return 0L;
    HPDDM::Option& opt = *HPDDM::Option::get();
    const std::string& prefix = ptA->prefix();
    KN<double>* timing = nargs[0] ? GetAny<KN<double>*>((*nargs[0])(stack)) : 0;
    Pair<K>* pair = nargs[2] ? GetAny<Pair<K>*>((*nargs[2])(stack)) : 0;
    if(opt.set(prefix + "schwarz_coarse_correction") && pair)
        if(pair->p) {
            int flag;
            MPI_Test(&(pair->p->first), &flag, MPI_STATUS_IGNORE);
        }
    MatriceMorse<K>* mA = nargs[3] ? static_cast<MatriceMorse<K>*>(&(*GetAny<Matrice_Creuse<K>*>((*nargs[3])(stack))->A)) : 0;
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
        ff_HPDDM_MatrixCSR<K> dA(mA);//->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
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
    bool excluded = nargs[1] && GetAny<bool>((*nargs[1])(stack));
    if(excluded)
        opt[prefix + "level_2_exclude"];
    if(pair)
        if(pair->p) {
            MPI_Wait(&(pair->p->first), MPI_STATUS_IGNORE);
            delete [] pair->p->second;
            pair->destroy();
            pair = nullptr;
        }
    MPI_Comm comm = nargs[4] ? *(MPI_Comm*)GetAny<pcommworld>((*nargs[4])(stack)) : MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(ptA->getCommunicator(), &rank);
    if(rank != mpirank || rank != 0) {
        opt.remove("verbosity");
        if(prefix.size() > 0)
            opt.remove(prefix + "verbosity");
    }
    double timer;
    if(timing) { // tic
        MPI_Barrier(comm);
        timer = MPI_Wtime();
        timing->resize(timing->n + 1);
        (*timing)[timing->n - 1] = timer;
    }
    if(!excluded) {
        const auto& map = ptA->getMap();
        bool allocate = map.size() > 0 && ptA->getBuffer()[0] == nullptr ? ptA->setBuffer() : false;
        ptA->exchange(static_cast<K*>(*ptRHS), mu);
        ptA->clearBuffer(allocate);
        HPDDM::IterativeMethod::solve(*ptA, (K*)*ptRHS, (K*)*ptX, mu, comm);
    }
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
        static const int n_name_param = 1;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        changeOperator_Op(const basicAC_F0& args, Expression param1, Expression param2) : A(param1), mat(param2) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type changeOperator_Op<Type, K>::name_param[] = {
    {"scaled", &typeid(bool)}
};
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
    HPDDM::MatrixCSR<K>* dN = new_HPDDM_MatrixCSR<K>(mN);//mN->n, mN->m, mN->nbcoef, mN->a, mN->lg, mN->cl, mN->symetrique);
    Type* ptA = GetAny<Type*>((*A)(stack));
    ptA->setMatrix(dN);
    if(!nargs[0] || GetAny<bool>((*nargs[0])(stack)))
        ptA->exchange();
    return 0L;
}

template<class Type, class K>
class set_Op : public E_F0mps {
    public:
        Expression A;
        static const int n_name_param = 2;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        set_Op(const basicAC_F0& args, Expression param) : A(param) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }

        AnyType operator()(Stack stack) const;
};
template<class Type, class K>
basicAC_F0::name_and_type set_Op<Type, K>::name_param[] = {
    {"sparams", &typeid(string*)},
    {"prefix", &typeid(string*)}
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
    if(nargs[0]) {
        HPDDM::Option::get()->parse(*(GetAny<string*>((*nargs[0])(stack))));
#ifdef PETSCSUB
        PetscOptionsInsertString(NULL, (GetAny<string*>((*nargs[0])(stack)))->c_str());
#endif
    }
    if(nargs[1]) {
        Type* ptA = GetAny<Type*>((*A)(stack));
        ptA->setPrefix(*(GetAny<string*>((*nargs[1])(stack))));
    }
    return 0L;
}

template<class Type, class K>
class distributedMV_Op : public E_F0mps {
    public:
        Expression A;
        Expression Mat;
        Expression in;
        Expression out;
        static const int n_name_param = 0;
        distributedMV_Op<Type, K>(const basicAC_F0& args, Expression param1, Expression param2, Expression param3, Expression param4) : A(param1), Mat(param2), in(param3), out(param4) {
            args.SetNameParam(n_name_param, nullptr, nullptr);
        }

        AnyType operator()(Stack stack) const;
};
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
    ff_HPDDM_MatrixCSR<K> dA(mA);//->n, mA->m, mA->nbcoef, mA->a, mA->lg, mA->cl, mA->symetrique);
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
            *Ax = K();
            A.prod(Ax);
            return Ax;
        }
        static U init(U Ax, ProdSchwarz<T, U, K, N> A) {
            Ax->init(A.u->n);
            return mv(Ax, A);
        }
};

template<class T, class U, class K, char trans>
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
            const auto& map = (*t).getMap();
            bool allocate = map.size() > 0 && (*t).getBuffer()[0] == nullptr ? (*t).setBuffer() : false;
            (*t).exchange(static_cast<K*>(*u), mu);
            (*t).clearBuffer(allocate);
            HPDDM::IterativeMethod::solve(*t, (K*)*u, (K*)*out, mu, MPI_COMM_WORLD);
            if(HPDDM::Wrapper<K>::I == 'F' && mu > 1) {
                std::for_each(A->_ja, A->_ja + A->_nnz, [](int& i) { --i; });
                std::for_each(A->_ia, A->_ia + A->_n + 1, [](int& i) { --i; });
            }
        }
        static U inv(U Ax, InvSchwarz<T, U, K, trans> A) {
            A.solve(Ax);
            return Ax;
        }
        static U init(U Ax, InvSchwarz<T, U, K, trans> A) {
            Ax->init(A.u->n);
            return inv(Ax, A);
        }
};

template<class R, char S>
class IterativeMethod : public OneOperator {
    public:
        const int c;
        typedef KN<R> Kn;
        typedef KN_<R> Kn_;
        class MatF_O : RNM_VirtualMatrix<R> {
            public:
                Stack stack;
                mutable Kn x;
                C_F0 c_x;
                Expression mat;
                typedef typename RNM_VirtualMatrix<R>::plusAx plusAx;
                MatF_O(int n, Stack stk, const OneOperator* op) :
                    RNM_VirtualMatrix<R>(n), stack(stk), x(n), c_x(CPValue(x)),
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
                void mv(const R* const in, const int& n, const int& m, R* const out) const {
                    ffassert(m == 1);
                    KN_<R> xx((R*)in, n);
                    KN_<R> yy(out, n);
                    yy = R();
                    addMatMul(xx,yy);
                }
                bool ChecknbLine(int) const { return true; }
                bool ChecknbColumn(int) const { return true; }
        };
        class MatMatF_O : RNM_VirtualMatrix<R> {
            public:
                Stack stack;
                mutable KNM<R> x;
                C_F0 c_x;
                Expression mat;
                typedef typename RNM_VirtualMatrix<R>::plusAx plusAx;
                MatMatF_O(int n, int m, Stack stk, const OneOperator* op) :
                    RNM_VirtualMatrix<R>(n), stack(stk), x(n, m), c_x(CPValue(x)),
                    mat(op ? CastTo<KNM_<R>>(C_F0(op->code(basicAC_F0_wa(c_x)), (aType)*op)) : 0) /* */ { }
                ~MatMatF_O() {
                    delete mat;
                    Expression zzz = c_x;
                    delete zzz;
                }
                void addMatMul(const KN_<R>& xx, KN_<R>& Ax) const { }
                void addMatMul(const KNM_<R>& xx, KNM_<R>& Ax) const {
                    ffassert(xx.N() == Ax.N());
                    ffassert(xx.M() == Ax.M());
                    x.resize(xx.N(), xx.M());
                    x = xx;
                    Ax = GetAny<KNM_<R>>((*mat)(stack));
                    WhereStackOfPtr2Free(stack)->clean();
                }
                void mv(const R* const in, const int& n, const int& m, R* const out) const {
                    KNM_<R> xx((R*)in, n, m);
                    KNM_<R> yy(out, n, m);
                    yy = R();
                    addMatMul(xx,yy);
                }
                bool ChecknbLine(int) const { return true; }
                bool ChecknbColumn(int) const { return true; }
        };
        template<class Op>
        class Operator : public HPDDM::EmptyOperator<R> {
            public:
                Op& mat;
                Op& prec;
                Operator(Op& m, Op& p) : HPDDM::EmptyOperator<R>(m.x.N()), mat(m), prec(p) { }
                void GMV(const R* const in, R* const out, const int& mu = 1) const {
                    mat.mv(in, HPDDM::EmptyOperator<R>::_n, mu, out);
                }
                template<bool>
                void apply(const R* const in, R* const out, const unsigned short& mu = 1, R* = nullptr, const unsigned short& = 0) const {
                    if(prec.mat)
                        prec.mv(in, HPDDM::EmptyOperator<R>::_n, mu, out);
                    else
                        std::copy_n(in, HPDDM::EmptyOperator<R>::_n * mu, out);
                }
        };
        template<class Op>
        class SchwarzOperator : public HpSchwarz<R, S> {
            public:
                Op& prec;
                SchwarzOperator(HpSchwarz<R, S>* m, Op& p) : prec(p) { *reinterpret_cast<HpSchwarz<R, S>*>(this) = *m; }
                template<bool>
                void apply(const R* const in, R* const out, const unsigned short& mu = 1, R* = nullptr, const unsigned short& = 0) const {
                    if(prec.mat)
                        prec.mv(in, HPDDM::Subdomain<R>::_dof, mu, out);
                    else
                        std::copy_n(in, HPDDM::Subdomain<R>::_dof * mu, out);
                }
        };
        class E_LCG : public E_F0mps {
            public:
                static const int n_name_param = 4;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                const OneOperator *A, *C;
                Expression Op;
                Expression B, X;
                const int c;
                E_LCG(const basicAC_F0& args, int d) : Op(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    if(c == 0 || c == 2) {
                        const Polymorphic* op = dynamic_cast<const Polymorphic*>(args[0].LeftValue());
                        ffassert(op);
                        A = (c == 0 || c == 1 ? op->Find("(", ArrayOfaType(atype<Kn*>(), false)) : op->Find("(", ArrayOfaType(atype<KNM<R>*>(), false)));
                    }
                    else {
                        Op = to<HpSchwarz<R, S>*>(args[0]);
                    }
                    if(nargs[0]) {
                        const Polymorphic* op = dynamic_cast<const Polymorphic*>(nargs[0]);
                        ffassert(op);
                        C = (c == 0 || c == 1 ? op->Find("(", ArrayOfaType(atype<Kn*>(), false)) : op->Find("(", ArrayOfaType(atype<KNM<R>*>(), false)));
                    }
                    else
                        C = 0;
                    if(c == 0 || c == 1) {
                        B = to<Kn*>(args[1]);
                        X = to<Kn*>(args[2]);
                    }
                    else {
                        B = to<KNM<R>*>(args[1]);
                        X = to<KNM<R>*>(args[2]);
                    }
                }
                virtual AnyType operator()(Stack stack)  const {
                    int ret = -1;
                    try {
                        MPI_Comm comm = nargs[3] ? *(MPI_Comm*)GetAny<pcommworld>((*nargs[3])(stack)) : MPI_COMM_WORLD;
                        if(nargs[2])
                            HPDDM::Option::get()->parse(*(GetAny<string*>((*nargs[2])(stack))));
                        if(c == 0 || c == 1) {
                            Kn& x = *GetAny<Kn*>((*X)(stack));
                            int n = x.N();
                            Kn& b = *GetAny<Kn*>((*B)(stack));
                            MatF_O PP(n, stack, C);
                            if(c == 0) {
                                MatF_O AA(n, stack, A);
                                Operator<MatF_O> Op(AA, PP);
                                if(nargs[1])
                                    Op.setPrefix(*(GetAny<string*>((*nargs[1])(stack))));
                                ret = HPDDM::IterativeMethod::solve(Op, (R*)b, (R*)x, 1, comm);
                            }
                            else {
                                HpSchwarz<R, S>* op = GetAny<HpSchwarz<R, S>*>((*Op)(stack));
                                SchwarzOperator<MatF_O> SchwarzOp(op, PP);
                                ret = HPDDM::IterativeMethod::solve(SchwarzOp, (R*)b, (R*)x, 1, comm);
                            }
                        }
                        else {
                            KNM<R>& x = *GetAny<KNM<R>*>((*X)(stack));
                            int n = x.N();
                            int m = x.M();
                            KNM<R>& b = *GetAny<KNM<R>*>((*B)(stack));
                            MatMatF_O PP(n, m, stack, C);
                            if(c == 2) {
                                MatMatF_O AA(n, m, stack, A);
                                Operator<MatMatF_O> Op(AA, PP);
                                if(nargs[1])
                                    Op.setPrefix(*(GetAny<string*>((*nargs[1])(stack))));
                                ret = HPDDM::IterativeMethod::solve(Op, (R*)b, (R*)x, m, comm);
                            }
                            else {
                                HpSchwarz<R, S>* op = GetAny<HpSchwarz<R, S>*>((*Op)(stack));
                                SchwarzOperator<MatMatF_O> SchwarzOp(op, PP);
                                ret = HPDDM::IterativeMethod::solve(SchwarzOp, (R*)b, (R*)x, m, comm);
                            }
                        }
                    }
                    catch(...) {
                        throw;
                    }
                    return SetAny<long>(ret);
                }
                operator aType() const { return atype<long>(); }
        };
        E_F0* code(const basicAC_F0& args) const { return new E_LCG(args, c); }
        IterativeMethod() : OneOperator(atype<long>(), atype<Polymorphic*>(), atype<KN<R>*>(), atype<KN<R>*>()), c(0) { }
        IterativeMethod(int) : OneOperator(atype<long>(), atype<HpSchwarz<R, S>*>(), atype<KN<R>*>(), atype<KN<R>*>()), c(1) { }
        IterativeMethod(int, int) : OneOperator(atype<long>(), atype<Polymorphic*>(), atype<KNM<R>*>(), atype<KNM<R>*>()), c(2) { }
        IterativeMethod(int, int, int) : OneOperator(atype<long>(), atype<HpSchwarz<R, S>*>(), atype<KNM<R>*>(), atype<KNM<R>*>()), c(3) { }
};

template<class R, char S>
basicAC_F0::name_and_type IterativeMethod<R, S>::E_LCG::name_param[] = {
    {"precon", &typeid(Polymorphic*)},
    {"prefix", &typeid(string*)},
    {"sparams", &typeid(string*)},
    {"communicator", &typeid(pcommworld)}
};

template<class Type>
long globalNumbering(Type* const& A, KN<long>* const& numbering) {
    if(A) {
        numbering->resize(2 + A->getMatrix()->_n);
        unsigned int g;
        unsigned int* num = reinterpret_cast<unsigned int*>(&((*numbering)[0]));
        A->distributedNumbering(num + 2, num[0], num[1], g);
    }
    return 0L;
}

template<class Type, class K>
Type* changeOperatorSimple(Type* const& A, Type* const& B) {
    *A = *B;
    return A;
}

template<template<class, char> class Type, class K, char S, char U = S>
void add() {
    Dcl_Type<Type<K, S>*>(Initialize<Type<K, S>>, Delete<Type<K, S>>);
#ifndef PETSCSUB
    if(std::is_same<K, HPDDM::underlying_type<K>>::value)
#endif
        zzzfff->Add("schwarz", atype<HpSchwarz<K, S>*>());
#ifndef PETSCSUB
    map_type_of_map[make_pair(atype<Type<HPDDM::underlying_type<K>, U>*>(), atype<K*>())] = atype<Type<K, S>*>();
#else
    map_type_of_map[make_pair(atype<Type<K, U>*>(), atype<K*>())] = atype<Type<K, S>*>();
#endif

    TheOperators->Add("<-", new initDDM<Type<K, S>, K>);
    TheOperators->Add("<-", new initDDM<Type<K, S>, K>(1));
    TheOperators->Add("<-", new initDDM<Type<K, S>, K>(1, 1));
    TheOperators->Add("<-", new initDDM<Type<K, S>, K>(1, 1, 1));
    TheOperators->Add("=", new OneOperator2_<Type<K, S>*, Type<K, S>*, Type<K, S>*>(Schwarz::changeOperatorSimple<Type<K, S>, K>));
    Global.Add("attachCoarseOperator", "(", new attachCoarseOperator<Type<K, S>, K>);
    Global.Add("attachCoarseOperator", "(", new attachCoarseOperator<Type<K, S>, K>(1));
    Global.Add("DDM", "(", new solveDDM<Type<K, S>, K>);
    Global.Add("changeOperator", "(", new changeOperator<Type<K, S>, K>);
    Global.Add("set", "(", new set<Type<K, S>, K>);
    addProd<Type<K, S>, ProdSchwarz, KN<K>, K>();
    addInv<Type<K, S>, InvSchwarz, KN<K>, K>();
    addScalarProduct< Type<K, S>, K >( );
#if 0 // if you need this, please make sure you are using the master branch of HPDDM
    addArray<Type<K, S>>();
#endif
    Global.Add("dmv", "(", new distributedMV<Type<K, S>, K>);
    Global.Add("destroyRecycling", "(", new OneOperator1_<bool, Type<K, S>*>(destroyRecycling<Type<K, S>, K>));
    Global.Add("statistics", "(", new OneOperator1_<bool, Type<K, S>*>(statistics<Type<K, S>>));
    Global.Add("exchange", "(", new exchangeIn<Type<K, S>, K>);
    Global.Add("exchange", "(", new exchangeInOut<Type<K, S>, K>);
    Global.Add("IterativeMethod","(",new IterativeMethod<K, S>());
    Global.Add("IterativeMethod","(",new IterativeMethod<K, S>(1));
    Global.Add("IterativeMethod","(",new IterativeMethod<K, S>(1, 1));
    Global.Add("IterativeMethod","(",new IterativeMethod<K, S>(1, 1, 1));
    Global.Add("globalNumbering", "(", new OneOperator2_<long, Type<K, S>*, KN<long>*>(globalNumbering<Type<K, S>>));

    if(!exist_type<Pair<K>*>()) {
        Dcl_Type<Pair<K>*>(InitP<Pair<K>>, Destroy<Pair<K>>);
#ifndef PETSCSUB
        map_type_of_map[make_pair(atype<Pair<HPDDM::underlying_type<K>>*>(), atype<K*>())] = atype<Pair<K>*>();
#else
        map_type_of_map[make_pair(atype<Pair<K>*>(), atype<K*>())] = atype<Pair<K>*>();
#endif
    }
    aType t;
    int r;
    if(!zzzfff->InMotClef("pair", t, r) && std::is_same<K, HPDDM::underlying_type<K>>::value)
        zzzfff->Add("pair", atype<Pair<K>*>());
}
}

static void Init_Schwarz() {
    Init_Common();
#if defined(DSUITESPARSE) || defined(DHYPRE) || defined(PETSCSUB)
    constexpr char ds = 'G';
#else
    constexpr char ds = 'S';
#endif
    constexpr char zs = 'G';
#ifndef PETSCSUB
    Schwarz::add<HpSchwarz, double, ds>();
#ifndef DHYPRE
    Schwarz::add<HpSchwarz, std::complex<double>, zs, ds>();
    // Schwarz::add<HpSchwarz, float, ds>();
    // Schwarz::add<HpSchwarz, std::complex<float>, zs>();
#endif
#else
    Schwarz::add<HpSchwarz, PetscScalar, ds>();
#endif
}

LOADFUNC(Init_Schwarz)
