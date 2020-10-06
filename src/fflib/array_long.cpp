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
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

#include "array_tlp.hpp"
#include "array_init.hpp"
const basicForEachType *aatypeknlongp;

/*
void initArrayOperators() {
  ArrayOperator<double>();
  ArrayOperator<Complex>();
  ArrayOperator<long>();
}

void initArrayDCL() {
  ArrayDCL<double>();
  ArrayDCL<Complex>();
  ArrayDCL<long>();
}
*/

aType aaaa_knlp;
void initArrayDCLlong() {
  // ArrayOperator<long>();
  Dcl_Type<Inv_KN_long>(); // Add FH mars 2005
  ArrayDCL<long>();
  aaaa_knlp = atype<KN<long>*>();
}

class OneBinaryOperatorInv_KN_long : public OneOperator {
  public:
    OneBinaryOperatorInv_KN_long(basicForEachType * ti) : OneOperator(atype<Inv_KN_long >(), ti ,atype<long>()) {}
    E_F0 * code(const basicAC_F0 & args) const {
      Expression p=args[1];
      if (!p->EvaluableWithOutStack()) {
        bool bb = p->EvaluableWithOutStack();
        cout << bb << " " << *p << endl;
        CompileError("Inverse: int[int] I, array, with I^p, The p must be a constant == -1, sorry");
      }
      long pv = GetAny<long>((*p)(NullStack));
      if (pv !=-1) {
        char buf[100];
        sprintf(buf, "Inverse: int[int] I, array, I^%ld, The pow must be == -1, sorry", pv);
        CompileError(buf);
      }
      return new E_F_F0<Inv_KN_long, KN_<long> >(Build<Inv_KN_long, KN_<long> >, to< KN_<long> >(args[0]));
    }
};

// Add mars 2010
template<class R> R *set_init_init( R* const & a,const long & n) {
  SHOWVERB(cout << " set_init " << typeid(R).name() << " " << n << endl);
  a->init(n);
  for (int i = 0; i < n; i++)
    (*a)[i].init();
  return a;
}

inline string **get_elements(KN<String> *const &a, long const &b) {
  ffassert(a && b >=0 && b < a->size());
  String &Sret = (*a)[b]; // correction FH feb 2004
  // delete b; la chaine est detruire automatiquement en fin d'instruction FH jan 2010
  return Sret.getap();
}

template<class A> inline AnyType Destroy_KN(Stack,const AnyType &x) {
  KN<A> *a = GetAny<KN<A>*>(x);
  for (int i = 0; i <a->N(); i++)
    (a)[i].destroy();
  a->destroy();
  return Nothing;
}
// end add

template<class A,class B>
struct set_Inv_KN_long : public binary_function<A,B,A> {
  static A f(const A & a, B const & b) {
    int n = a.N();
    KN_<long> I(b.t);
    for (int i = 0; i < I.N(); ++i) {
      int j = I[i];
      if (j >= 0 && j < n)
        a[j] = i;
    }
    return a;
  }
};

template<class A,class B>
struct set_Inv_pKN_longI: public binary_function<A,B,A> {
  static A f(const A & a, B const & b) {
    KN_<long> I(b.t);
    int n = I.max() + 1;
    a->init(n);
    (*a) = -1;
    for (int i = 0; i < I.N(); ++i) {
      int j = I[i];
      if (j >= 0 && j < n)
        (*a)[j] = i;
    }
    return a;
  }
};

long findall(const  KN_<long> & a,  const long &v,  KN<long> * const &  pI)
{
    long nn=0,k=0;;
    KN<long> & I=*pI;
    for(int i=0; i<a.N();++i)
       if( a[i]==v) nn++;
    I.resize(nn);
    for(int i=0; i<a.N();++i)
    if( a[i]==v) I[k++]=i;
    return nn;
}



void initArrayOperatorlong()
{
  typedef long K;
  Dcl_Type< Eye > ();// OK this is the fist array def ..
  Global.Add("eye","(",new OneOperator1<Eye,long>(fEye));
  Global.Add("eye","(",new OneOperator2<Eye,long>(fEye));
  Global.Add("findall", "(", new OneOperator3_< long, KN_<long>, long ,KN<long>*>(findall));// oct 2020 FH.

  ArrayOperator<long, long>();
  // to define inverse permutation // Add FH mars 2005
  TheOperators->Add("^", new OneBinaryOperatorInv_KN_long(atype<KN_<long> >()));
  //- TheOperators->Add("^", new OneBinaryOperatorInv_KN_long(atype<KN<long> *>()));
  aatypeknlongp = atype<KN<long>*>(); // for compilation error with g++ 3.2.2

  Add<KN_<long> >("sort", ".", new OneOperator1_<KN_<K>, KN_<K> >(SortKn<K, KN_<K> >));
  // Add<KN<long> >("sort", ".", new OneOperator1_<KN<K>,KN<K> >(SortKn<K, KN<K> >));
  Add<KN<long> *>("sort", ".", new OneOperator1_<KN<K>*, KN<K>* >(SortpKn<K>));
  Global.Add("sort", "(", new OneOperator2_<KN<K>*, KN<K>*, KN<long>* >(SortpKn2<K,long>));

  // ArrayDCL<long>();
  Dcl_TypeandPtr_<KN_<String>, KN<String> *>(0, 0, 0, ::Destroy<KN<String> >, ::ClearReturnKK_<K,KN<String>, KN_<String> >, ::ClearReturnpKK<String, KN<String> >);
  atype<KN<String>* >()->Add("[", "", new OneOperator2_<string**, KN<String>*, long >(get_elements));
  TheOperators->Add("<-", new OneOperator2_<KN<String> *, KN<String> *, long>(&set_init_init), new InitArrayfromArray<string*, KN<String>*, true>);
  map_type_of_map[make_pair(atype<long>(), atype<string*>())] = atype<KN<String>*>(); // vector
  Add<KN<String> *>("n", ".", new OneOperator1<long, KN<String> *>(get_n));
  extern KN<String> *pkarg;
  Global.New("ARGV", CPValue<KN<String> >(*pkarg));// add FH mars 2010
  Global.Add("toZarray", "(", new OneOperator_2KN_<long>);
  TheOperators->Add("=", new OneBinaryOperator<set_Inv_KN_long<KN_<long>, Inv_KN_long> >);
  TheOperators->Add("<-", new OneBinaryOperator<set_Inv_pKN_longI<KN<long>*, Inv_KN_long> >);

  Add<KN<K> *>("imin", ".", new OneOperator1<long, KN<K> *>(get_imin));
  Add<KN<K> *>("imax", ".", new OneOperator1<long, KN<K> *>(get_imax));
  Add<KNM<K> *>("imin", ".", new OneOperator1<long, KNM<K> *>(get_imin)); // Add april 2018 FH
  Add<KNM<K> *>("imax", ".", new OneOperator1<long, KNM<K> *>(get_imax)); // Add april 2018 FH
  Add<KNM<K> *>("jmin", ".", new OneOperator1<long, KNM<K> *>(get_jmin)); // Add april 2018 FH
  Add<KNM<K> *>("jmax", ".", new OneOperator1<long, KNM<K> *>(get_jmax)); // Add april 2018 FH
  TheOperators->Add("ijmax", new OneOperator3_<NothingType,KNM<K> *, long*, long*>(get_ijmax)); // Add april 2018 FH
  TheOperators->Add("ijmin", new OneOperator3_<NothingType,KNM<K> *, long*, long*>(get_ijmin)); // Add april 2018 FH
  // madd FH. march 2015 ...
  Global.Add("Unique", "(", new Unique<K, K>);
  Global.Add("Unique", "(", new Unique<K, double>);
  // convertion double -> long (via lround)
  Dcl_Type<F_KN_<K, K, double> >();
  Global.Add("lround", "(", new OneOperator1F_KN_<F_KN_<K, K, double>, K, double, KN_<double> >(lround));
  TheOperators->Add("=", new OneBinaryOperator<set_eq_array<KN_<K>, F_KN_<K, K,double> > >); // add FH juin 2005
  TheOperators->Add("+=", new OneBinaryOperator<set_eq_array_add<KN_<K>, F_KN_<K, K, double> > >); // add FH juin 2005
  TheOperators->Add("-=", new OneBinaryOperator<set_eq_array_sub<KN_<K>, F_KN_<K, K, double> > >); // add FH juin 2005
  TheOperators->Add("/=", new OneBinaryOperator<set_eq_array_div<KN_<K>, F_KN_<K, K, double> > >); // add FH juin 2005
  TheOperators->Add("*=", new OneBinaryOperator<set_eq_array_mul<KN_<K>, F_KN_<K, K, double> > >); // add FH juin 2005

  TheOperators->Add("<-", new InitMapfromArray<MyMap<String,String>*, string *, string*, true> );

  TheOperators->Add("<<", new OneBinaryOperator<PrintPnd<KN<String> * > >); // add may 2018 FH
}

// void xxxx() {
// }
