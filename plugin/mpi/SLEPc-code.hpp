#include "PETSc.hpp"

#ifdef WITH_slepc
#define WITH_SLEPC
#endif
#ifdef WITH_slepccomplex
#define WITH_SLEPC
#endif

#ifdef WITH_SLEPC

#include "slepc.h"

namespace SLEPc {
template<class Type, class K>
struct _m_User;
template<class Type, class K, class SType>
struct _n_User;
template<class Type, class K, class SType>
using User = typename std::conditional<!std::is_same<SType, NEP>::value, _n_User<Type, K, SType>*, _m_User<Type, K>*>::type;
template<class Type, class K, class SType>
static PetscErrorCode MatMult_User(Mat A, Vec x, Vec y);
template<class K, typename std::enable_if<std::is_same<K, double>::value || !std::is_same<PetscScalar, double>::value>::type* = nullptr>
void copy(K* pt, PetscInt n, PetscScalar* xr, PetscScalar* xi) { }
template<class K, typename std::enable_if<!std::is_same<K, double>::value && std::is_same<PetscScalar, double>::value>::type* = nullptr>
void copy(K* pt, PetscInt n, PetscScalar* xr, PetscScalar* xi) {
    for(int i = 0; i < n; ++i)
        pt[i] = K(xr[i], xi[i]);
}
template<class SType, class K, typename std::enable_if<!std::is_same<SType, SVD>::value && (std::is_same<K, double>::value || !std::is_same<PetscScalar, double>::value)>::type* = nullptr>
void assign(K* pt, PetscScalar& kr, PetscScalar& ki) {
    *pt = kr;
}
template<class SType, class K, typename std::enable_if<!std::is_same<SType, SVD>::value && (!std::is_same<K, double>::value && std::is_same<PetscScalar, double>::value)>::type* = nullptr>
void assign(K* pt, PetscScalar& kr, PetscScalar& ki) {
    *pt = K(kr, ki);
}
template<class SType, class K, typename std::enable_if<std::is_same<SType, SVD>::value>::type* = nullptr>
void assign(K* pt, PetscScalar& kr, PetscScalar& ki) { }
template<class K, typename std::enable_if<(std::is_same<PetscScalar, double>::value && std::is_same<K, std::complex<double>>::value)>::type* = nullptr>
void distributedVec(PetscInt* num, PetscInt first, PetscInt last, K* const in, PetscScalar* pt, PetscInt n) { }
template<class K, typename std::enable_if<!(std::is_same<PetscScalar, double>::value && std::is_same<K, std::complex<double>>::value)>::type* = nullptr>
void distributedVec(PetscInt* num, PetscInt first, PetscInt last, K* const in, PetscScalar* pt, PetscInt n) {
    HPDDM::Subdomain<K>::template distributedVec<0>(num, first, last, in, pt, n, 1);
}
template<class Type, class K, class SType>
class eigensolver : public OneOperator {
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
        };
        class ScalarF_O {
            public:
                Stack stack;
                mutable PetscScalar x;
                C_F0 c_x;
                Expression mat;
                ScalarF_O(Stack stk, const OneOperator* op)
                    : stack(stk), x(0), c_x(CPValue(x)),
                    mat(op ? CastTo< long >(C_F0(op->code(basicAC_F0_wa(c_x)), (aType)*op)) : 0) {}
                ~ScalarF_O( ) {
                    delete mat;
                    Expression zzz = c_x;
                    delete zzz;
                }
                long apply(const PetscScalar xx, int* f) const {
                    x = xx;
                    ffassert(mat);
                    int ret = GetAny< int >((*mat)(stack));
                    *f = ret;
                    WhereStackOfPtr2Free(stack)->clean( );
                    return 0;
                }
        };
        const int c;
        class E_eigensolver : public E_F0mps {
            public:
                Expression A;
                Expression B;
                const OneOperator* codeA, *codeB;
                const int c;
                static const int n_name_param = 10;
                static basicAC_F0::name_and_type name_param[];
                Expression nargs[n_name_param];
                E_eigensolver(const basicAC_F0& args, int d) : A(0), B(0), codeA(0), codeB(0), c(d) {
                    args.SetNameParam(n_name_param, name_param, nargs);
                    if(c != 4) {
                        A = to<Type*>(args[0]);
                        if(c == 1 || c == 3) {
                            const Polymorphic* op = dynamic_cast<const Polymorphic*>(args[1].LeftValue());
                            ffassert(op);
                            if(c == 3) {
                                codeA = op->Find("(", ArrayOfaType(atype<PetscScalar>(), false));
                                ffassert(codeA);
                                B = to<Type*>(args[2]);
                                op = dynamic_cast<const Polymorphic*>(args[3].LeftValue());
                                ffassert(op);
                                codeB = op->Find("(", ArrayOfaType(atype<PetscScalar>(), false));
                                ffassert(codeB);
                            }
                            else
                                codeA = op->Find("(", ArrayOfaType(atype<KN<PetscScalar>*>(), false));
                        }
                        else if(c == 0) {
                            B = to<Type*>(args[1]);
                        }
                    }
                    else {
                        A = to<KN<Type>*>(args[0]);
                    }
                }

                AnyType operator()(Stack stack) const;
                operator aType() const { return atype<long>(); }
        };
        E_F0* code(const basicAC_F0 & args) const { return new E_eigensolver(args, c); }
        eigensolver() : OneOperator(atype<long>(), atype<Type*>(), atype<Type*>()), c(0) { }
        eigensolver(int) : OneOperator(atype<long>(), atype<Type*>(), atype<Polymorphic*>()), c(1) { }
        eigensolver(int, int) : OneOperator(atype<long>(), atype<Type*>()), c(2) { }
        eigensolver(int, int, int) : OneOperator(atype<long>(), atype<Type*>(), atype<Polymorphic*>(), atype<Type*>(), atype<Polymorphic*>()), c(3) { }
        eigensolver(int, int, int, int) : OneOperator(atype<long>(), atype<KN<Type>*>()), c(4) { }
};
template<class Type, class K, class SType>
basicAC_F0::name_and_type eigensolver<Type, K, SType>::E_eigensolver::name_param[] = {
    {"sparams", &typeid(std::string*)},
    {"prefix", &typeid(std::string*)},
    {"values", &typeid(KN<typename std::conditional<!std::is_same<SType, SVD>::value, K, PetscReal>::type>*)},
    {!std::is_same<SType, SVD>::value ? "vectors" : "lvectors", &typeid(FEbaseArrayKn<K>*)},
    {!std::is_same<SType, SVD>::value ? "array" : "larray", &typeid(KNM<K>*)},
    {"fields", &typeid(KN<double>*)},
    {"names", &typeid(KN<String>*)},
    {!std::is_same<SType, SVD>::value ? "schurPreconditioner" : "rvectors", !std::is_same<SType, SVD>::value ? &typeid(KN<Matrice_Creuse<HPDDM::upscaled_type<PetscScalar>>>*) : &typeid(FEbaseArrayKn<K>*)},
    {!std::is_same<SType, SVD>::value ? "schurList" : "rarray", !std::is_same<SType, SVD>::value ? &typeid(KN<double>*) : &typeid(KNM<K>*)},
    {"deflation", &typeid(KNM<PetscScalar>*)},
};
template<class Type, class K, class SType>
struct _n_User {
    typename eigensolver<Type, K, SType>::MatF_O* mat;
};
template<class Type, class K>
struct _m_User {
    typename eigensolver<Type, K, NEP>::ScalarF_O* F;
    typename eigensolver<Type, K, NEP>::ScalarF_O* J;
};
PetscErrorCode FormFun(NEP nep, PetscScalar lambda, Mat F, Mat P, void* ctx) {
    User<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, NEP>* user;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, NEP >* >(ctx);
    typename eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, NEP>::ScalarF_O* f =
        reinterpret_cast< typename eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, NEP>::ScalarF_O* >((*user)->F);
    int ret;
    f->apply(lambda, &ret);
    PetscFunctionReturn(ret);
}
PetscErrorCode FormJac(NEP nep, PetscScalar lambda, Mat J, void* ctx) {
    User<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, NEP>* user;

    PetscFunctionBeginUser;
    user = reinterpret_cast< User< PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, NEP >* >(ctx);
    typename eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, NEP>::ScalarF_O* f =
        reinterpret_cast< typename eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, NEP>::ScalarF_O* >((*user)->J);
    int ret;
    f->apply(lambda, &ret);
    PetscFunctionReturn(ret);
}
template<class Type, class K, class SType>
AnyType eigensolver<Type, K, SType>::E_eigensolver::operator()(Stack stack) const {
    if(A && (c == 2 || c == 4 || B || codeA)) {
        KN<Type>* ptTab;
        Type* ptA;
        if(c != 4)
            ptA = GetAny<Type*>((*A)(stack));
        else {
            ptTab = GetAny<KN<Type>*>((*A)(stack));
            ffassert(ptTab && ptTab->N());
            ptA = &(ptTab->operator[](0));
        }
        if(ptA->_petsc) {
            EPS eps;
            SVD svd;
            NEP nep;
            PEP pep;
            if(std::is_same<SType, EPS>::value)
                EPSCreate(PetscObjectComm((PetscObject)ptA->_petsc), &eps);
            else if(std::is_same<SType, SVD>::value)
                SVDCreate(PetscObjectComm((PetscObject)ptA->_petsc), &svd);
            else if(std::is_same<SType, NEP>::value)
                NEPCreate(PetscObjectComm((PetscObject)ptA->_petsc), &nep);
            else
                PEPCreate(PetscObjectComm((PetscObject)ptA->_petsc), &pep);
            Mat S;
            User<Type, K, typename std::conditional<std::is_same<SType, NEP>::value, EPS, SType>::type> user = nullptr;
            User<Type, K, NEP> func = nullptr;
            MatType type;
            PetscBool isType;
            MatGetType(ptA->_petsc, &type);
            PetscStrcmp(type, MATNEST, &isType);
            PetscInt m;
            if(!codeA) {
                Type* ptB = (c == 0 ? GetAny<Type*>((*B)(stack)) : NULL);
                if(std::is_same<SType, EPS>::value)
                    EPSSetOperators(eps, ptA->_petsc, c == 0 && ptB ? ptB->_petsc : NULL);
                else if(std::is_same<SType, SVD>::value)
                    SVDSetOperator(svd, ptA->_petsc);
                else if(std::is_same<SType, PEP>::value) {
                    Mat* tab = new Mat[ptTab->N()];
                    for(int i = 0; i < ptTab->N(); ++i)
                        tab[i] = ptTab->operator[](i)._petsc;
                    PEPSetOperators(pep, ptTab->N(), tab);
                    delete [] tab;
                }
                if(!ptA->_A)
                    MatGetLocalSize(ptA->_petsc, &m, NULL);
            }
            else {
                MatGetLocalSize(ptA->_petsc, &m, NULL);
                PetscInt M;
                MatGetSize(ptA->_petsc, &M, NULL);
                if(!std::is_same<SType, NEP>::value) {
                    PetscNew(&user);
                    user->mat = new typename eigensolver<Type, K, typename std::conditional<std::is_same<SType, NEP>::value, EPS, SType>::type>::MatF_O(m, stack, codeA);
                    MatCreateShell(PetscObjectComm((PetscObject)ptA->_petsc), m, m, M, M, user, &S);
                    MatShellSetOperation(S, MATOP_MULT, (void (*)(void))MatMult_User<Type, K, typename std::conditional<std::is_same<SType, NEP>::value, EPS, SType>::type>);
                    if(std::is_same<SType, EPS>::value)
                        EPSSetOperators(eps, S, NULL);
                    else if(std::is_same<SType, SVD>::value)
                        SVDSetOperator(svd, S);
                }
                else {
                    Type* ptB = GetAny<Type*>((*B)(stack));
                    PetscNew(&func);
                    func->F = new typename eigensolver<Type, K, NEP>::ScalarF_O(stack, codeA);
                    func->J = new typename eigensolver<Type, K, NEP>::ScalarF_O(stack, codeB);
                    NEPSetFunction(nep, ptA->_petsc, ptA->_petsc, FormFun, &func);
                    NEPSetJacobian(nep, ptB->_petsc, FormJac, &func);
                }
            }
            if (nargs[0]) {
                std::string* options = GetAny< std::string* >((*nargs[0])(stack));
                PetscOptionsInsertString(NULL, options->c_str());
            }
            if(nargs[1]) {
                if(std::is_same<SType, EPS>::value)
                    EPSSetOptionsPrefix(eps, GetAny<std::string*>((*nargs[1])(stack))->c_str());
                else if(std::is_same<SType, SVD>::value)
                    SVDSetOptionsPrefix(svd, GetAny<std::string*>((*nargs[1])(stack))->c_str());
                else if(std::is_same<SType, NEP>::value)
                    NEPSetOptionsPrefix(nep, GetAny<std::string*>((*nargs[1])(stack))->c_str());
                else
                    PEPSetOptionsPrefix(pep, GetAny<std::string*>((*nargs[1])(stack))->c_str());
            }
            if(std::is_same<SType, EPS>::value) {
                EPSSetFromOptions(eps);
                ST st;
                EPSGetST(eps, &st);
                if(ptA->_ksp)
                    STSetKSP(st, ptA->_ksp);
                else {
                    KSP ksp;
                    PC pc;
                    STGetKSP(st, &ksp);
                    KSPGetPC(ksp, &pc);
                    PCSetFromOptions(pc);
                    PCType type;
                    PCGetType(pc, &type);
                    PetscBool isFieldSplit;
                    PetscStrcmp(type, PCFIELDSPLIT, &isFieldSplit);
                    if(isFieldSplit) {
                        KN<double>* fields = nargs[5] ? GetAny<KN<double>*>((*nargs[5])(stack)) : 0;
                        KN<String>* names = nargs[6] ? GetAny<KN<String>*>((*nargs[6])(stack)) : 0;
                        KN<Matrice_Creuse<HPDDM::upscaled_type<PetscScalar>>>* mS = nargs[7] ? GetAny<KN<Matrice_Creuse<HPDDM::upscaled_type<PetscScalar>>>*>((*nargs[7])(stack)) : 0;
                        KN<double>* pL = nargs[8] ? GetAny<KN<double>*>((*nargs[8])(stack)) : 0;
                        if(fields && names) {
                            KSP ksp;
                            STGetKSP(st, &ksp);
                            KSPSetOperators(ksp, ptA->_petsc, ptA->_petsc);
                            setFieldSplitPC(ptA, ksp, fields, names, mS, pL);
                            EPSSetUp(eps);
                            if(ptA->_vS && !ptA->_vS->empty()) {
                                PC pc;
                                KSPGetPC(ksp, &pc);
                                PCSetUp(pc);
                                PETSc::setCompositePC(pc, ptA->_vS);
                            }
                        }
                    }
                }
            }
            else if(std::is_same<SType, SVD>::value) {
                SVDSetFromOptions(svd);
                SVDSetUp(svd);
            }
            else if(std::is_same<SType, NEP>::value) {
                NEPSetFromOptions(nep);
                NEPSetUp(nep);
            }
            else {
                PEPSetFromOptions(pep);
                PEPSetUp(pep);
            }
            FEbaseArrayKn<K>* eigenvectors = nargs[3] ? GetAny<FEbaseArrayKn<K>*>((*nargs[3])(stack)) : nullptr;
            Vec* basis = nullptr;
            PetscInt n = 0;
            if(eigenvectors) {
                ffassert(!isType);
                if(eigenvectors->N > 0 && eigenvectors->get(0) && eigenvectors->get(0)->n > 0) {
                    n = eigenvectors->N;
                    basis = new Vec[n];
                    for(int i = 0; i < n; ++i) {
                        MatCreateVecs(ptA->_petsc, &basis[i], NULL);
                        PetscScalar* pt;
                        VecGetArray(basis[i], &pt);
                        if(!(std::is_same<PetscScalar, double>::value && std::is_same<K, std::complex<double>>::value))
                            distributedVec(ptA->_num, ptA->_first, ptA->_last, static_cast<K*>(*(eigenvectors->get(i))), pt, eigenvectors->get(i)->n);
                        VecRestoreArray(basis[i], &pt);
                    }
                }
                eigenvectors->resize(0);
            }
            if(std::is_same<SType, EPS>::value) {
                if(n)
                    EPSSetInitialSpace(eps, n, basis);
                KNM<PetscScalar>* ptDeflation = nargs[9] ? GetAny<KNM<PetscScalar>*>((*nargs[9])(stack)) : NULL;
                if(ptDeflation && ptDeflation->M()) {
                    PetscInt m;
                    MatGetLocalSize(ptA->_petsc, &m, NULL);
                    ffassert(m == ptDeflation->N());
                    Vec* deflation = new Vec[ptDeflation->M()];
                    for(int i = 0; i < ptDeflation->M(); ++i)
                        VecCreateMPIWithArray(PetscObjectComm((PetscObject)ptA->_petsc), 1, ptDeflation->N(), PETSC_DECIDE, &ptDeflation->operator( )(0, i), deflation + i);
                    EPSSetDeflationSpace(eps, ptDeflation->M(), deflation);
                    for(int i = 0; i < ptDeflation->M(); ++i)
                        VecDestroy(deflation + i);
                    delete [] deflation;
                }
                EPSSolve(eps);
            }
            else if(std::is_same<SType, SVD>::value) {
                if(n)
                    SVDSetInitialSpaces(svd, n, basis, 0, NULL);
                SVDSolve(svd);
            }
            else if(std::is_same<SType, NEP>::value) {
                if(n)
                    NEPSetInitialSpace(nep, n, basis);
                NEPSolve(nep);
            }
            else {
                if(n)
                    PEPSetInitialSpace(pep, n, basis);
                PEPSolve(pep);
            }
            for(int i = 0; i < n; ++i)
                VecDestroy(&basis[i]);
            delete [] basis;
            PetscInt nconv;
            if(std::is_same<SType, EPS>::value)
                EPSGetConverged(eps, &nconv);
            else if(std::is_same<SType, SVD>::value)
                SVDGetConverged(svd, &nconv);
            else if(std::is_same<SType, NEP>::value)
                NEPGetConverged(nep, &nconv);
            else
                PEPGetConverged(pep, &nconv);
            if(nconv > 0 && ((nargs[2] || nargs[3] || nargs[4]) || (std::is_same<SType, SVD>::value && (nargs[7] || nargs[8])))) {
                KN<typename std::conditional<!std::is_same<SType, SVD>::value, K, PetscReal>::type>* eigenvalues = nargs[2] ? GetAny<KN<typename std::conditional<!std::is_same<SType, SVD>::value, K, PetscReal>::type>*>((*nargs[2])(stack)) : nullptr;
                KNM<K>* array = nargs[4] ? GetAny<KNM<K>*>((*nargs[4])(stack)) : nullptr;
                FEbaseArrayKn<K>* rvectors = std::is_same<SType, SVD>::value && nargs[7] ? GetAny<FEbaseArrayKn<K>*>((*nargs[7])(stack)) : nullptr;
                KNM<K>* rarray = std::is_same<SType, SVD>::value && nargs[8] ? GetAny<KNM<K>*>((*nargs[8])(stack)) : nullptr;
                if(eigenvalues)
                    eigenvalues->resize(nconv);
                if(eigenvectors && !isType)
                    eigenvectors->resize(nconv);
                if(rvectors && !isType)
                    rvectors->resize(nconv);
                if(array)
                    array->resize(!codeA && !isType && ptA->_A ? ptA->_A->getDof() : m, nconv);
                if(rarray)
                    rarray->resize(!codeA && !isType && ptA->_A ? ptA->_A->getDof() : m, nconv);
                Vec xr, xi;
                PetscInt n;
                if(eigenvectors || array || rvectors || rarray) {
                    MatCreateVecs(ptA->_petsc, &xi, &xr);
                    VecGetLocalSize(xr, &n);
                } else xr = xi = NULL;
                for(PetscInt i = 0; i < nconv; ++i) {
                    PetscScalar kr, ki = 0;
                    PetscReal sigma;
                    if(std::is_same<SType, EPS>::value)
                        EPSGetEigenpair(eps, i, &kr, &ki, (eigenvectors || array) ? xr : NULL, (eigenvectors || array) && std::is_same<PetscScalar, double>::value && std::is_same<K, std::complex<double>>::value ? xi : NULL);
                    else if(std::is_same<SType, SVD>::value)
                        SVDGetSingularTriplet(svd, i, &sigma, xr, xi);
                    else if(std::is_same<SType, NEP>::value)
                        NEPGetEigenpair(nep, i, &kr, &ki, (eigenvectors || array) ? xr : NULL, (eigenvectors || array) && std::is_same<PetscScalar, double>::value && std::is_same<K, std::complex<double>>::value ? xi : NULL);
                    else
                        PEPGetEigenpair(pep, i, &kr, &ki, (eigenvectors || array) ? xr : NULL, (eigenvectors || array) && std::is_same<PetscScalar, double>::value && std::is_same<K, std::complex<double>>::value ? xi : NULL);
                    if(eigenvectors || array || rvectors || rarray) {
                        PetscScalar* tmpr;
                        PetscScalar* tmpi;
                        VecGetArray(xr, &tmpr);
                        K* pt, *pti;
                        if(!std::is_same<SType, SVD>::value) {
                            if(std::is_same<PetscScalar, double>::value && std::is_same<K, std::complex<double>>::value) {
                                VecGetArray(xi, &tmpi);
                                pt = new K[n];
                                copy(pt, n, tmpr, tmpi);
                            }
                            else
                                pt = reinterpret_cast<K*>(tmpr);
                        }
                        else {
                            pt = reinterpret_cast<K*>(tmpr);
                            VecGetArray(xi, reinterpret_cast<PetscScalar**>(&pti));
                        }
                        if(!isType && ptA->_A) {
                            KN<K> cpy(ptA->_A->getDof());
                            cpy = K(0.0);
                            HPDDM::Subdomain<K>::template distributedVec<1>(ptA->_num, ptA->_first, ptA->_last, static_cast<K*>(cpy), pt, static_cast<PetscInt>(cpy.n), 1);
                            ptA->_A->HPDDM::template Subdomain<PetscScalar>::exchange(static_cast<K*>(cpy));
                            if(eigenvectors)
                                eigenvectors->set(i, cpy);
                            if(array && !codeA)
                                (*array)(':', i) = cpy;
                            if(std::is_same<SType, SVD>::value) {
                                /* TODO FIXME: handle rectangular case */
                                HPDDM::Subdomain<K>::template distributedVec<1>(ptA->_num, ptA->_first, ptA->_last, static_cast<K*>(cpy), pti, static_cast<PetscInt>(cpy.n), 1);
                                ptA->_A->HPDDM::template Subdomain<PetscScalar>::exchange(static_cast<K*>(cpy));
                                if(rvectors)
                                    rvectors->set(i, cpy);
                                if(rarray && !codeA)
                                    (*rarray)(':', i) = cpy;
                            }
                        }
                        if(codeA || isType || !ptA->_A) {
                            if(array) {
                                KN<K> cpy(m, pt);
                                (*array)(':', i) = cpy;
                            }
                            if(rarray) {
                                KN<K> cpy(m, pti);
                                (*rarray)(':', i) = cpy;
                            }
                        }
                        if(!std::is_same<SType, SVD>::value && std::is_same<PetscScalar, double>::value && std::is_same<K, std::complex<double>>::value)
                            delete [] pt;
                        else
                            VecRestoreArray(xi, &tmpi);
                        VecRestoreArray(xr, &tmpr);
                    }
                    if(eigenvalues) {
                        if(!std::is_same<SType, SVD>::value)
                            assign<SType>(static_cast<typename std::conditional<!std::is_same<SType, SVD>::value, K*, PetscReal*>::type>(*eigenvalues + i), kr, ki);
                        else
                            eigenvalues->operator[](i) = sigma;
                    }
                }
                if(eigenvectors || array || rvectors || rarray) {
                    VecDestroy(&xr);
                    VecDestroy(&xi);
                }
            }
            if(user) {
                MatDestroy(&S);
                delete user->mat;
                PetscFree(user);
            }
            if(std::is_same<SType, EPS>::value)
                EPSDestroy(&eps);
            else if(std::is_same<SType, SVD>::value)
                SVDDestroy(&svd);
            else if(std::is_same<SType, NEP>::value)
                NEPDestroy(&nep);
            else
                PEPDestroy(&pep);
            return static_cast<long>(nconv);
        }
        else
            return 0L;
    }
    else
        return 0L;
}
template<class Type, class K, class SType>
static PetscErrorCode MatMult_User(Mat A, Vec x, Vec y) {
    User<Type, K, SType>   user;
    const PetscScalar*       in;
    PetscScalar*            out;
    PetscErrorCode         ierr;

    PetscFunctionBegin;
    ierr = MatShellGetContext(A, &user); CHKERRQ(ierr);
    typename SLEPc::eigensolver<Type, K, SType>::MatF_O* mat = reinterpret_cast<typename SLEPc::eigensolver<Type, K, SType>::MatF_O*>(user->mat);
    VecGetArrayRead(x, &in);
    VecGetArray(y, &out);
    KN_<PetscScalar> xx(const_cast<PetscScalar*>(in), mat->N);
    KN_<PetscScalar> yy(out, mat->N);
    yy = *mat * xx;
    VecRestoreArray(y, &out);
    VecRestoreArrayRead(x, &in);
    PetscFunctionReturn(0);
}
void finalizeSLEPc() {
    SlepcFinalize();
}
template<class K, typename std::enable_if<std::is_same<K, double>::value>::type* = nullptr>
void addSLEPc() {
    Global.Add("EPSSolveComplex", "(", new SLEPc::eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, std::complex<double>, EPS>());
    Global.Add("EPSSolveComplex", "(", new SLEPc::eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, std::complex<double>, EPS>(1));
    Global.Add("EPSSolveComplex", "(", new SLEPc::eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, std::complex<double>, EPS>(1, 1));
}
template<class K, typename std::enable_if<!std::is_same<K, double>::value>::type* = nullptr>
void addSLEPc() { }
}

static void Init() {
    //  to load only once
    aType t;
    int r;
#ifdef WITH_slepccomplex
    const char * mmmm= "Petsc Slepc complex";
#else
    const char * mmmm= "Petsc Slepc real";
#endif
    if(!zzzfff->InMotClef(mmmm,t,r))
    {
#ifdef PETScandSLEPc
        Init_PETSc();
#endif
        int argc = pkarg->n;
        char** argv = new char*[argc];
        for(int i = 0; i < argc; ++i)
            argv[i] = const_cast<char*>((*(*pkarg)[i].getap())->c_str());
        PetscBool isInitialized;
        PetscInitialized(&isInitialized);
        if(!isInitialized && mpirank == 0)
            std::cout << "PetscInitialize has not been called, do not forget to load PETSc before loading SLEPc" << std::endl;
        SlepcInitialize(&argc, &argv, 0, "");
        delete [] argv;
        ff_atend(SLEPc::finalizeSLEPc);
        SLEPc::addSLEPc<PetscScalar>();
        Global.Add("EPSSolve", "(", new SLEPc::eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, EPS>());
        Global.Add("EPSSolve", "(", new SLEPc::eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, EPS>(1));
        Global.Add("EPSSolve", "(", new SLEPc::eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, EPS>(1, 1));
        Global.Add("SVDSolve", "(", new SLEPc::eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, SVD>(1, 1));
        Global.Add("NEPSolve", "(", new SLEPc::eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, NEP>(1, 1, 1));
        Global.Add("PEPSolve", "(", new SLEPc::eigensolver<PETSc::DistributedCSR<HpSchwarz<PetscScalar>>, PetscScalar, PEP>(1, 1, 1, 1));
        if(verbosity>1)cout << "*** End:: load PETSc & SELPc "<< typeid(PetscScalar).name() <<"\n\n"<<endl;
        zzzfff->Add(mmmm, atype<Dmat*>());
    }
    else {
        if(verbosity>1)cout << "*** reload and skip load PETSc & SELPc "<< typeid(PetscScalar).name() <<"\n\n"<<endl;
    }
}
#else
static void Init() {
     Init_PETSc();
}
#endif
