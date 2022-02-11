/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

// Example C++ function "myfunction", dynamically loaded into "ff-c++ dfft.cpp "

/* clang-format off */
//ff-c++-LIBRARY-dep: fftw3
//ff-c++-cpp-dep:
/* clang-format on */

#include "ff++.hpp"
#include "AFunction_ext.hpp"
#include <fftw3.h>

template< class Complex >
class DFFT_1d2dor3d {
public:
    Complex *x;
    int n, m, k;
    int sign;
    DFFT_1d2dor3d(KN< Complex > *xx, long signn, long nn = 1, long kk = 1)
    : x(*xx), n(nn), m(xx->N( ) / (nn * kk)), k(kk), sign(signn) {
        cout << xx << " " << signn << " " << nn << " " << xx->N( ) << " n: " << n << " m:" << m
        << " k:  " << k << endl;
        ffassert(n > 0 && (n * m * k == xx->N( )));
    }
    
    DFFT_1d2dor3d(KNM< Complex > *xx, long signn) : x(*xx), n(xx->M( )), m(xx->N( )), sign(signn) {}
};

DFFT_1d2dor3d< Complex > dfft(KN< Complex > *const &x, const long &sign) {
    return DFFT_1d2dor3d< Complex >(x, sign);
}

DFFT_1d2dor3d< Complex > dfft(KN< Complex > *const &x, const long &nn, const long &sign) {
    return DFFT_1d2dor3d< Complex >(x, sign, nn);
}

DFFT_1d2dor3d< Complex > dfft(KN< Complex > *const &x, const long &nn, const long &kk,
                              const long &sign) {
    return DFFT_1d2dor3d< Complex >(x, sign, nn, kk);
}

DFFT_1d2dor3d< Complex > dfft(KNM< Complex > *const &x, const long &sign) {
    return DFFT_1d2dor3d< Complex >(x, sign);
}

bool ff_execute(fftw_plan *p) {
    if (*p) {
        fftw_execute(*p);
    }
    
    return 0;
}

bool ff_delete(fftw_plan *p) {
    if (*p) {
        fftw_destroy_plan(*p);
    }
    
    *p = 0;
    return 0;
}

KN< Complex > *dfft_eq(KN< Complex > *const &x, const DFFT_1d2dor3d< Complex > &d) {
    ffassert(x->N( ) == d.n * d.m * d.k);
    Complex *px = *x;
    fftw_plan p;
    if (d.k == 1) {
        if (d.n > 1) {
            p = fftw_plan_dft_2d(d.n, d.m, reinterpret_cast< fftw_complex * >(d.x),
                                 reinterpret_cast< fftw_complex * >(px), d.sign, FFTW_ESTIMATE);
        } else {
            p = fftw_plan_dft_1d(d.m, reinterpret_cast< fftw_complex * >(d.x),
                                 reinterpret_cast< fftw_complex * >(px), d.sign, FFTW_ESTIMATE);
        }
    } else {
        if (d.n > 1) {
            p = fftw_plan_dft_3d(d.n, d.m, d.k, reinterpret_cast< fftw_complex * >(d.x),
                                 reinterpret_cast< fftw_complex * >(px), d.sign, FFTW_ESTIMATE);
        } else {
            p = fftw_plan_dft_2d(d.m, d.k, reinterpret_cast< fftw_complex * >(d.x),
                                 reinterpret_cast< fftw_complex * >(px), d.sign, FFTW_ESTIMATE);
        }
    }
    
    fftw_execute(p);
    fftw_destroy_plan(p);
    return x;
}

KN< double > *dfft_eq(KN< double > *const &x, const DFFT_1d2dor3d< double > &d) {
    ffassert(0);
    return x;
}

/*  class Init { public:
 * Init();
 * };*/
// bofbof ..
struct fftw_plan_s {};

// ...

template<>
inline AnyType DeletePtr< fftw_plan * >(Stack, const AnyType &x) {
    fftw_plan *a = PGetAny< fftw_plan >(x);
    
    if (*a) {
        fftw_destroy_plan(*a);
    }
    
    *a = 0;
    return Nothing;
};

fftw_plan *plan__eq(fftw_plan *a, fftw_plan b) {
    if (*a) {
        fftw_destroy_plan(*a);
    }
    
    *a = b;
    return a;
}

fftw_plan *plan_set(fftw_plan *a, fftw_plan b) {
    *a = b;
    return a;
}

fftw_plan plan_dfft(KN< Complex > *const &x, KN< Complex > *const &y, const long &sign) {
    return fftw_plan_dft_1d(x->N( ), reinterpret_cast< fftw_complex * >(&x[0]),
                            reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
}

fftw_plan plan_dfft(KNM< Complex > *const &x, KNM< Complex > *const &y, const long &sign) {
    long m = x->N( ), n = x->M( );
    
    fftw_plan_dft_2d(n, m, reinterpret_cast< fftw_complex * >(&x[0]),
                     reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
    return 0;
}

fftw_plan plan_dfft(KN< Complex > *const &x, KN< Complex > *const &y, const long &n,
                    const long &sign) {
    long nn = n, mm = y->N( ) / nn;
    
    ffassert(mm * nn == y->N( ) && x->N( ) == y->N( ));
    
    return fftw_plan_dft_2d(nn, mm, reinterpret_cast< fftw_complex * >(&x[0]),
                            reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
}

fftw_plan plan_dfft(KN< Complex > *const &x, KN< Complex > *const &y, const long &n, const long &k,
                    const long &sign) {
    int nn = n, mm = y->N( ) / (k * n), kk = k;
    
    ffassert(y->N( ) == nn * mm * kk);
    if (nn > 1) {
        return fftw_plan_dft_3d(nn, mm, kk, reinterpret_cast< fftw_complex * >(&x[0]),
                                reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
    } else {
        return fftw_plan_dft_2d(nn, mm, reinterpret_cast< fftw_complex * >(&x[0]),
                                reinterpret_cast< fftw_complex * >(&y[0]), sign, FFTW_ESTIMATE);
    }
}
template <int NP>
class Mapkk : public E_F0mps {
public:
    typedef Complex R;
    typedef KN_< R > Result;
    ;
    static basicAC_F0::name_and_type *name_param;
    static const int n_name_param = 0;
    Expression expv, expK, expm, expk, exp;
    Expression nargs[n_name_param];
    
    Mapkk(const basicAC_F0 &args) : expv(0),expK(0), expm(0), expk(0), exp(0) {
        args.SetNameParam(n_name_param, name_param, nargs);
        expv = to< KN< R > * >(args[0]);    // a the expression to get the mesh
        expK = to<  R3  * >(args[1]);    // a the expression to get the K fourier variable
        if(NP>2)expm = to< long >(args[2]);// 2d ...
        if(NP==4)
            expk = to< long >(args[3]);
        exp = to< R >(args[NP]);    // a the expression to get the mesh
    }
    
    ~Mapkk( ) {}
    
    static ArrayOfaType typeargs( ) {
        if(NP==2)
            return ArrayOfaType(atype< KN< R > * >( ), atype< R3* >( ), atype< R >( ));
        else if(NP==3)
            return ArrayOfaType(atype< KN< R > * >( ), atype< R3* >( ), atype< long >( ), atype< R >( ));
        else if(NP==4)
            return ArrayOfaType(atype< KN< R > * >( ), atype< R3* >( ), atype< long >( ),atype< long >( ), atype< R >( ));
        else ffassert(0); //
    }
    
    static E_F0 *f(const basicAC_F0 &args) { return new Mapkk<NP>(args); }
    
    AnyType operator( )(Stack s) const;
};
template <int NP>
basicAC_F0::name_and_type *Mapkk<NP>::name_param = 0;
template <int NP>
AnyType Mapkk<NP>::operator( )(Stack s) const {
    // correct July 2015 ... not tested before..
    MeshPoint *mp(MeshPointStack(s)), mps = *mp;
    
    KN< R > *pv = GetAny< KN< R > * >((*expv)(s));
    KN< R > & v(*pv);
    R3 * pK = GetAny< R3 * >((*expK)(s));
    long nn = v.N( );
    long n2 = expm ?GetAny< long >((*expm)(s)):1;
    long n3 = expk  ? GetAny< long >((*expk)(s)): 1;
    if (verbosity > 9) {
        cout << "  map: expm " << expm << " n2 = " << n2 << " n3 =" << n3 << " size array:" << nn << endl;
    }
    long n23 =n2*n3;
    long n1 = nn / n23;
    double k1 = 1. / n1;
    double k2 = 1. / n2;
    double k3 = 1.  / n3;
    
    double k10 = 0., k20 = 0, k30 = 0;;
    if (verbosity > 9) {
        cout << " map: " << n1 << " " << n2 << " " << n3 << " " << nn << " == " << n1 * n2 * n3 << endl;
    }
    
    ffassert(n1 * n2 * n3  == nn);
    long n12 = (n1 + 1) / 2, n22 = (n2 + 1) / 2, n32 = (n3 + 1) / 2;
    int kkk =0; 
    for (long i3 = 0; i3 < n3; ++i3)
    for (long i2 = 0 ; i2 < n2; ++i2)
    for (long i1 = 0; i1 < n1; ++i1,++kkk)
    {
        int ii1 = i1%n12 - (i1/n12)*n12;
        int ii2 = i2%n22 - (i2/n22)*n22;
        int ii3 = i3%n32 - (i3/n32)*n32;
        R3 P(ii1,ii2, ii3);
        *pK = P; // set value of K to P.
        v[kkk] = GetAny< R >((*exp)(s));
        if(verbosity>19) cout <<  "map" << kkk << " " <<ii1 << " " << ii2 << " " << ii3 << " " << v[kkk] << " P=" << P << endl;
    }
    ffassert(kkk==nn);
    
    *mp = mps;
    return 0L;
}

static void Load_Init( ) {
    typedef DFFT_1d2dor3d< Complex > DFFT_C;
    typedef DFFT_1d2dor3d< double > DFFT_R;
    
    cout << " load: init dfft " << endl;
    Dcl_Type< DFFT_C >( );
    Dcl_Type< DFFT_R >( );
    
    Dcl_Type< fftw_plan * >(::InitializePtr< fftw_plan * >, ::DeletePtr< fftw_plan * >);
    Dcl_Type< fftw_plan >( );
    zzzfff->Add("fftwplan", atype< fftw_plan * >( ));
    
    TheOperators->Add("=", new OneOperator2< fftw_plan *, fftw_plan *, fftw_plan >(plan__eq));
    TheOperators->Add("<-", new OneOperator2< fftw_plan *, fftw_plan *, fftw_plan >(plan_set));
    
    Global.Add("plandfft", "(",
               new OneOperator3_< fftw_plan, KN< Complex > *, KN< Complex > *, long >(plan_dfft));
    Global.Add(
               "plandfft", "(",
               new OneOperator4_< fftw_plan, KN< Complex > *, KN< Complex > *, long, long >(plan_dfft));
    Global.Add(
               "plandfft", "(",
               new OneOperator5_< fftw_plan, KN< Complex > *, KN< Complex > *, long, long, long >(plan_dfft));
    Global.Add("plandfft", "(",
               new OneOperator3_< fftw_plan, KNM< Complex > *, KNM< Complex > *, long >(plan_dfft));
    
    Global.Add("execute", "(", new OneOperator1< bool, fftw_plan * >(ff_execute));
    Global.Add("delete", "(", new OneOperator1< bool, fftw_plan * >(ff_delete));
    
    Global.Add("dfft", "(", new OneOperator2_< DFFT_C, KN< Complex > *, long >(dfft));
    Global.Add("dfft", "(", new OneOperator3_< DFFT_C, KN< Complex > *, long, long >(dfft));
    Global.Add("dfft", "(", new OneOperator4_< DFFT_C, KN< Complex > *, long, long, long >(dfft));
    Global.Add("dfft", "(", new OneOperator2_< DFFT_C, KNM< Complex > *, long >(dfft));
    Global.Add("mapk", "(", new OneOperatorCode< Mapkk<2> >( ));
    Global.Add("mapkk", "(", new OneOperatorCode< Mapkk<3> >( ));
    Global.Add("mapkkk", "(", new OneOperatorCode< Mapkk<4> >( ));
    TheOperators->Add("=", new OneOperator2_< KN< Complex > *, KN< Complex > *, DFFT_C >(dfft_eq));
}

LOADFUNC(Load_Init)
