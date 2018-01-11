#ifndef _COMMON_
#define _COMMON_

#include <math.h>
#include <mpi.h>
#include <ff++.hpp>
#include <AFunction_ext.hpp>

#ifdef WITH_mkl
#define HPDDM_MKL 1
#endif

#if HPDDM_SCHWARZ || HPDDM_FETI || HPDDM_BDD
#ifdef WITH_mkl
#define MKL_PARDISOSUB
#elif defined(WITH_mumps)
#define MUMPSSUB
#else
#define SUITESPARSESUB
#endif

#ifdef WITH_mumps
#define DMUMPS
#else
#define DSUITESPARSE
#endif
#define MU_ARPACK
#endif

#define HPDDM_NUMBERING 'C'
#undef CBLAS_H

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
        bool empty() const { return _size <= 0; }
        T& operator[](std::size_t idx) { return _it[idx]; }
        const T& operator[](std::size_t idx) const { return _it[idx]; }
        T& back() { return _it[_size - 1]; }
        const T& back() const { return _it[_size - 1]; }
};
template<class K>
class Pair {
    public:
        Pair() : p() { };
        std::pair<MPI_Request, const K*>* p;
        void init() { }
        void destroy() {
            delete p;
            p = nullptr;
        }
};
template<class R, class A, class B> R Build(A a, B b) {
    return R(a, b);
}
template<class Op, class Inv>
class OneBinaryOperatorInv : public OneOperator {
    public:
        OneBinaryOperatorInv() : OneOperator(atype<Inv>(), atype<Op*>(), atype<long>()) { }
        E_F0* code(const basicAC_F0 & args) const {
            Expression p = args[1];
            if(!p->EvaluableWithOutStack())
                CompileError("A^p, The p must be a constant == -1, sorry");
            long pv = GetAny<long>((*p)(NullStack));
            if(pv != -1) {
                char buf[100];
                sprintf(buf, "A^%ld, The pow must be == -1, sorry", pv);
                CompileError(buf);
            }
            return new E_F_F0<Inv, Op*>(Build<Inv, Op*>, t[0]->CastTo(args[0]));
        }
};
template<class Op, template<class, class, class> class Inv, class V, class K = double>
void addInv() {
    class OpInv {
        public:
            Op* A;
            OpInv(Op* B) : A(B) { assert(A); }
            operator Op& () const { return *A; }
            operator Op* () const { return A; }
    };
    Dcl_Type<OpInv>();
    Dcl_Type<Inv<OpInv, V*, K>>();
    TheOperators->Add("^", new OneBinaryOperatorInv<Op, OpInv>());
    TheOperators->Add("*", new OneOperator2<Inv<OpInv, V*, K>, OpInv, V*>(Build));
    TheOperators->Add("=", new OneOperator2<V*, V*, Inv<OpInv, V*, K>>(Inv<OpInv, V*, K>::init));
}
template<class Op, template<class, class, class, char> class Prod, class V, class K = double, char N = 'N'>
void addProd() {
    Dcl_Type<Prod<Op*, V*, K, N>>();
    if(N == 'T') {
        class OpTrans {
            public:
                Op* A;
                OpTrans(Op* B) : A(B) { assert(A); }
                operator Op& () const { return *A; }
                operator Op* () const { return A; }
        };
        Dcl_Type<OpTrans>();
        TheOperators->Add("\'", new OneOperator1<OpTrans, Op*>(Build));
        TheOperators->Add("*", new OneOperator2<Prod<Op*, V*, K, N>, OpTrans, V*>(Build));
    }
    else {
        TheOperators->Add("*", new OneOperator2<Prod<Op*, V*, K, N>, Op*, V*>(Build));
    }
    TheOperators->Add("=", new OneOperator2<V*, V*, Prod<Op*, V*, K, N>>(Prod<Op*, V*, K, N>::mv));
    TheOperators->Add("<-", new OneOperator2<V*, V*, Prod<Op*, V*, K, N>>(Prod<Op*, V*, K, N>::init));
}

extern KN<String>* pkarg;

#if HPDDM_SCHWARZ || HPDDM_FETI || HPDDM_BDD
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
template<class Type>
bool statistics(Type* const& Op) {
    Op->statistics();
    return false;
}
template<template<class, char> class Type, class K, char S>
long exchange(Type<K, S>* const& Op, KN<K>* const& in, KN<K>* const& out) {
    *out = *in;
    const auto& map = Op->getMap();
    bool allocate = map.size() > 0 && Op->getBuffer()[0] == nullptr ? Op->setBuffer() : false;
    Op->HPDDM::template Subdomain<K>::exchange(static_cast<K*>(*out));
    Op->clearBuffer(allocate);
    *out -= *in;
    return 0L;
}

static void Init_Common() {
    aType t;
    int r;
    if(!zzzfff->InMotClef("dpair", t, r)) {
        Global.Add("getOption", "(", new OneOperator1_<double, string*>(getOpt));
        Global.Add("isSetOption", "(", new OneOperator1_<bool, string*>(isSetOpt));
        int argc = pkarg->n;
        const char** argv = new const char*[argc];
        for(int i = 0; i < argc; ++i)
            argv[i] = (*((*pkarg)[i].getap()))->data();
        HPDDM::Option::get()->parse(argc, argv, mpirank == 0);
        delete [] argv;
    }
}
#endif
#endif // _COMMON_
