#ifndef _COMMON_
#define _COMMON_

template<class R, class A, class B> R Build(A a, B b) {
    return R(a, b);
}
template<class RR, class AA = RR, class BB = AA>
struct BinaryOp {
    using first_argument_type  = AA;
    using second_argument_type = BB;
    using result_type          = RR;
    static RR f(Stack s, const AA& a, const BB& b) { return RR(s, a, b); }
};
template<class Op, char trans = 'N'>
class pwr {
    public:
        Op* A;
        pwr(Op* B) : A(B) { assert(A); }
        const typename std::conditional<trans == 'T', std::string*, long>::type c;
        static constexpr char tr = trans;
        mutable bool conjugate;
        pwr(Stack s, Op* const& d, const typename std::conditional<trans == 'T', std::string*, long>::type e) : A(d), c(e), conjugate(false) { }
        operator Op* () const { return A; }
};
template<class RR, class AA = RR, class BB = AA>
struct assign {
using first_argument_type  = AA;
using second_argument_type = BB;
using result_type          = RR;
template<class V, class T>
static T check(T* t, typename std::enable_if<T::tr == 'H'>::type* = 0) {
    t->conjugate = true;
    if(t->c != -1)
        CompileError("A'^p, the p must be a constant == -1, sorry");
    return *t;
}
template<class V, class T>
static T check(T* t, typename std::enable_if<T::tr == 'N'>::type* = 0) {
    if(t->c != -1)
        CompileError("A^p, the p must be a constant == -1 or == \"-T\" or == \"-H\", sorry");
    return *t;
}
template<class V, class T>
static T check(T* t, typename std::enable_if<T::tr == 'T'>::type* = 0) {
    if(t->c->compare("-H") == 0)
        t->conjugate = true;
    if(t->c->compare("-T") != 0 && !t->conjugate)
        CompileError("A^p, the p must be a constant == -1 or == \"-T\" or == \"-H\", sorry");
    return *t;
}
static RR f(Stack stack, const AA& a, const BB& b) {
    ffassert(a);
    check<BB>(&a);
    RR p(a, b);
    return p;
}
};
template<class Op>
class OpTrans {
    public:
        Op* A;
        OpTrans(Op* B) : A(B) { assert(A); }
        operator Op& () const { return *A; }
        operator Op* () const { return A; }
};
template<class Op, template<class, class, class, char> class Inv, class V, class K = double, char trans = 'N'>
void addInv() {
    Dcl_Type<pwr<Op, trans>>();
    Dcl_Type<Inv<pwr<Op, trans>, V*, K, trans>>();
    if(trans == 'T') {
        Dcl_Type<pwr<Op, 'H'>>();
        Dcl_Type<Inv<pwr<Op, 'H'>, V*, K, 'T'>>();
        TheOperators->Add("^", new OneBinaryOperator_st<BinaryOp<pwr<Op, 'H'>, OpTrans<Op>, typename std::conditional<'N' == 'N', long, std::string*>::type>>);
        TheOperators->Add("*", new OneBinaryOperator_st<assign<Inv<pwr<Op, 'H'>, V*, K, 'T'>, pwr<Op, 'H'>, V*>>);
        TheOperators->Add("=", new OneOperator2<V*, V*, Inv<pwr<Op, 'H'>, V*, K, 'T'>>(Inv<pwr<Op, 'H'>, V*, K, 'T'>::inv));
        TheOperators->Add("<-", new OneOperator2<V*, V*, Inv<pwr<Op, 'H'>, V*, K, 'T'>>(Inv<pwr<Op, 'H'>, V*, K, 'T'>::init));
    }
    TheOperators->Add("^", new OneBinaryOperator_st<BinaryOp<pwr<Op, trans>, Op*, typename std::conditional<trans == 'N', long, std::string*>::type>>);
    TheOperators->Add("*", new OneBinaryOperator_st<assign<Inv<pwr<Op, trans>, V*, K, trans>, pwr<Op, trans>, V*>>);
    TheOperators->Add("=", new OneOperator2<V*, V*, Inv<pwr<Op, trans>, V*, K, trans>>(Inv<pwr<Op, trans>, V*, K, trans>::inv));
    TheOperators->Add("<-", new OneOperator2<V*, V*, Inv<pwr<Op, trans>, V*, K, trans>>(Inv<pwr<Op, trans>, V*, K, trans>::init));
}
template<class Op, template<class, class, class, char> class Prod, class V, class K = double, char N = 'N'>
void addProd() {
    Dcl_Type<Prod<Op*, V*, K, N>>();
    if(N == 'T') {
        Dcl_Type<OpTrans<Op>>();
        TheOperators->Add("\'", new OneOperator1<OpTrans<Op>, Op*>(Build));
        TheOperators->Add("*", new OneOperator2<Prod<Op*, V*, K, N>, OpTrans<Op>, V*>(Build));
    }
    else
        TheOperators->Add("*", new OneOperator2<Prod<Op*, V*, K, N>, Op*, V*>(Build));
    TheOperators->Add("=", new OneOperator2<V*, V*, Prod<Op*, V*, K, N>>(Prod<Op*, V*, K, N>::mv));
    TheOperators->Add("<-", new OneOperator2<V*, V*, Prod<Op*, V*, K, N>>(Prod<Op*, V*, K, N>::init));
}

#endif // _COMMON_
