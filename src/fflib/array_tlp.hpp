// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*

 This file is part of Freefem++

 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
//#pragma dont_inline on
//#pragma inline_depth(1)

// TODO: remove this block as soon as autoconf is removed from FreeFem++
#ifndef CMAKE
#include <config.h>
#endif

#ifndef ARRAY_TLP_HPP
#define ARRAY_TLP_HPP

#include <set>
#include <complex>
//#include <type_traits>

#include "AFunction.hpp"
#include <cstdarg>
#include <cstring>
#include "error.hpp"
#include "lex.hpp"
#include "ufunction.hpp"

#include "RNM.hpp"

#include "Operator.hpp"
// for exec routine
#include "rgraph.hpp"
#include "InitFunct.hpp"
#include <queue>
#include "array_resize.hpp"
#include "HeapSort.hpp"
namespace Fem2D {
#include "R3.hpp"
}
template <class T>
struct affectation
{
	using first_argument_type  = T;
	using second_argument_type = T;
	using result_type          = T;
	T& operator()(T& x, const T& y) const {return (x=y);}
};

template <class T>
struct affectation_add
{
	using first_argument_type  = T;
	using second_argument_type = T;
	using result_type          = T;
	T& operator()(T& x, const T& y) const {return (x+=y);}// correct FH 25/10/2013
};

template <class T>
struct affectation_sub
{
	using first_argument_type  = T;
	using second_argument_type = T;
	using result_type          = T;
	T& operator()(T& x, const T& y) const {return (x-=y);}// correct FH 25/10/2013
};



extern Map_type_of_map map_type_of_map ; //  to store te type
extern Map_type_of_map map_pair_of_type ; //  to store te type

extern basicForEachType *  typevarreal,  * typevarcomplex;  //  type of real and complex variable

extern int TheCurrentLine; // unset: by default
extern long mpisize,mpirank;

template<class T> inline T Square (const T &a){return a*a;}



template<class K>
struct Op2_dotproduct {
  using first_argument_type  = Transpose<KN_<K> >;
  using second_argument_type = KN<K> *;
  using result_type          = K;
  static K f( Transpose<KN_<K> > const & a, KN<K> * const& b)
   { return (conj(a.t),*b);} };

template<class K>
struct Op2_dotproduct_ {
  using first_argument_type  = Transpose<KN_<K> >;
  using second_argument_type = KN_<K>;
  using result_type          = K;
  static K f( Transpose<KN_<K> > const & a, KN_<K>  const& b)
   { return (conj(a.t),b);} };


template<class T>
void  HeapSort(T *c,long n,long o)
{ // trie un tableau c de n valeur avec un decalage de o.
   //  le tableau: c[i*o] , pour i = 0 a n-1
    long l,j,r,i;
    T crit;
    c-=o; // on decale de o pour que le tableau commence a o
    if( n <= 1) return;
    l = (n/2 + 1)*o;
    r = n*o;
    while (1) { // label 2
	if(l <= o ) { // label 20
	    crit = c[r];
	    c[r] = c[o];
	    r-=o;
	    if ( r == o ) { c[o]=crit; return;}
	} else  crit = c[l-=o];
	j=l;
	while (1) {// label 4
	    i=j;
	    j=2*j;
	    if  (j>r) {c[i]=crit;break;} // L8 -> G2
	    if ((j<r) && (c[j] < c[j+o])) j+=o; // L5
	    if (crit < c[j]) c[i]=c[j]; // L6+1 G4
	    else {c[i]=crit;break;} //L8 -> G2
	}
    }
}

template<class R,class A> A  SortKn(const A  & ca)
{
    A a(ca);
    if(a.n > 0)
        HeapSort<R>(&a[0],a.n,a.step);
    return a;}

template<class R,class RR,class A,class B> A  SortKn(const A  & ca,const B  & cb)
{
    cout << "SortKn  " << endl;
    const A &a(ca);
    const B &b(cb);
    ffassert(a.n == b.n);
    ffassert(a.step == b.step && b.step ==1);
    if(a.n > 0)
        HeapSort<R,RR>(&a[0],&b[0],a.n);
    cout << b << endl;
return a;}

template<class R,class RR> KN<R> *  SortpKn2( KN<R> * const & pa,KN<RR> * const & pb){
  //  cout << " SortpKn2 " << endl;
    KN<R> &a(*pa);
    KN<RR> &b(*pb);
    ffassert(a.n == b.n);
    ffassert(a.step == b.step && b.step ==1);
    if(a.n > 0)
        HeapSort<R,RR>(&a[0],&b[0],a.n);
   return pa;}

template<class R> KN<R> *  SortpKn( KN<R> * const & pa){
    KN<R> &a(*pa);
    if(a.n > 0)
        HeapSort<R>(&a[0],a.n,a.step);
    return pa;}

template<class R>
class QuantileKN:  public KN_<R> { public:
    QuantileKN(const KN_<R> &a): KN_<R>(a) {}
    QuantileKN(KN<R>  * p): KN_<R>(*p) {}
    operator R *() const {return this->KN_<R>::operator R *() ;}
};


template<class R> R   Quantile(QuantileKN<R>  const & a,const double & q){
    KN<R> b(a);
    HeapSort<R>(b,b.n,b.step);
    long m=lrint(b.n*q);
    if( m >= b.n) m=b.n-1;
    if( m < 0) m=0;
    R qq=b[m];
   // cout <<  "Quantile: m = " << m << " " << b <<endl;
    return qq;}



inline void MyAssert(int i,char * ex,char * file,long line)
{if (i) {
    cout << "CompileError assertion :  " << ex << " in file " << file << "  line = " << line << endl;
     CompileError();}
 }


template<class K>
inline   K * get_element( MyMap<String,K> *  const  &  a,string*  const   & b)
 { K * ret=  &((*a)[*b]); // correction FH feb 2004
  //  cout << "get_element " << *b << " : " << ret << " = "<< * ret << endl;
   // delete b;  modif mars 2006 auto del ptr
    return ret;}

template<class K>
inline   bool exist_element( MyMap<String,K> *  const  &  a,string*  const   & b)
{     return  a->exist(*b);}

template<>
inline   string ** get_element<string*>( MyMap<String,string*> *  const  &  a,string*  const   & b)
 { string** ret=  &((*a)[*b]); // correction FH feb 2004
    if( *ret ==0) *ret = newstring(""); //  string vide ???
     // cout << "get_element " << *b << " : " << ret << " = "<< * ret << endl;
    // delete b;  modif mars 2006 auto del ptr
    return ret;}

inline   string ** get_elements( MyMap<String,String> *  const  &  a,string*  const   & b)
 { String* Sret=  &((*a)[*b]); // correction FH feb 2004
   //  delete b;  modif mars 2006 auto del ptr
    return Sret->getap();}

template<class RR,class A,class B>
RR * get_element_(const A & a,const B & b){
  if( b<0 || a.N() <= b)
   { cerr << " Out of bound  0 <=" << b << " < "  << a.N() << " array type = " << typeid(A).name() << endl;
     ExecError("Out of bound in operator []");}
    return  &((a)[b]);}


template<class RR,class A,class B>
RR * get_elementp_(const A & a,const B & b){
  if( b<0 || a->N() <= b)
   { cerr << " Out of bound  0 <=" << b << " < "  << a->N() << " array type = " << typeid(A).name() << endl;
     ExecError("Out of bound in operator []");}
    return  &((*a)[b]);}

template<class K>
KN_<K> fSubArray(const KN_<K> & a,const SubArray & b)
 { return a(b);}
template<class K>
KN_<K> fSubArrayp( KN<K>  * const & a,const SubArray & b)
 { return (*a)(b);}

template<class K>
KNM_<K> fSubArraybb(const KNM_<K> & a,const SubArray & b,const SubArray & c)
{ return a(b,c);}
template<class K>
KNM_<K> fSubArraypbb( KNM<K> * const & a,const SubArray & b,const SubArray & c)
{ return (*a)(b,c);}

template<class K>
KN_<K> fSubArrayib(const KNM_<K> & a,const long &i,const SubArray & b)
{ return a(i,b);}
template<class K>
KN_<K> fSubArraybi(const KNM_<K> & a,const SubArray & b,const long &i)
{ return a(b,i);}

template<class K>
KN_<K> fSubArraypib( KNM<K> *const & a,const long &i,const SubArray & b)
{ return (*a)(i,b);}
template<class K>
KN_<K> fSubArraypbi( KNM<K> *const & a,const SubArray & b,const long &i)
{ return (*a)(b,i);}


template<class A>
A fSubArrayc(const A & a,const char & b)
 { return a;}

template<class RR,class A,class B,class C>
RR * get_elementp2_(const A & a,const B & b,const C & c){
  if( b<0 || a->N() <= b || c<0 || a->M() <= c  )
   { cerr << " Out of bound  0 <=" << b << " < "  << a->N() << " " << c << " < "  << a->M()
           << " array type = " << typeid(A).name() << endl;
     ExecError("Out of bound in operator (,)");}
    return  &((*a)(b,c));}

template<class RR,class A,class B,class C>
RR get_element_is(const A &  a,const B & b,const C & c){
 //  cout << b << " .... " << ((*a)(SubArray(1,b),c)) << endl;;
    return  ((*a)(b,'.')(c));}

template<class RR,class A,class B,class C>
RR get_element_si(const A &  a,const B & b,const C & c){
 //  cout << c << " .... " << ((*a)(b,SubArray(1,c) )) << endl;;
     return  ((*a)('.',c)(b));}

template<class A,class B,class C>
struct check_get_element_lineorcol
{
  static  void check(const A &  a,const B & b,const C & c){
        if(c == ':' && (b<0 || a->N() <= b))
            ExecError("Out of bound");
        if(b == ':' && (c<0 || a->M() <= c))
            ExecError("Out of bound");

    }
};
template<class A,class B>
struct check_get_element_lineorcol<A,B, char>
{
    static  void check(const A &  a,const B & b,const char & c){
        if( (b<0 || a->N() <= b))
            ExecError("Out of bound");
    }
};
template<class A,class B>
struct check_get_element_lineorcol<A, char,B>
{
    static  void check(const A &  a,const char & b,const B & c){
        if( (c<0 || a->M() <= c))
            ExecError("Out of bound");
    }
};

template<class RR,class A,class B,class C>
RR get_element_lineorcol(const A &  a,const B & b,const C & c){
 //  cout << b << " .... " << ((*a)(SubArray(1,b),c)) << endl;;
     check_get_element_lineorcol<A,B,C>::check(a,b,c);
     return  ((*a)(b,c));
    }

template<class RR,class A,class B,class C>
RR get_element_is_(const A &  a,const B & b,const C & c){
    //  cout << b << " .... " << ((*a)(SubArray(1,b),c)) << endl;;
return  ((a)(b,'.')(c));}

template<class RR,class A,class B,class C>
RR get_element_si_(const A &  a,const B & b,const C & c){
    //  cout << c << " .... " << ((*a)(b,SubArray(1,c) )) << endl;;
return  ((a)('.',c)(b));}

template<class RR,class A,class B,class C>
RR get_element_lineorcol_(const A &  a,const B & b,const C & c){
    //  cout << b << " .... " << ((*a)(SubArray(1,b),c)) << endl;;
return  ((a)(b,c));}

template<class RR,class A,class B,class C>
RR * get_elementp2__(const A & a,const B & b,const C & c){
    if( b<0 || a.N() <= b || c<0 || a.M() <= c  )
      { cerr << " Out of bound  0 <=" << b << " < "  << a.N() << " " << c << " < "  << a.M()
	  << " array type = " << typeid(A).name() << endl;
      ExecError("Out of bound in operator (,)");}
return  &((a)(b,c));}


template<typename T>
struct myremove_pointer
{
    typedef T type;
};

template<typename T>
struct myremove_pointer<T*>
{
    typedef typename myremove_pointer<T>::type type;
};

template<class Map,class Key, class Value,bool isinit>
class  InitMapfromArray : public OneOperator {
public:
   typedef typename  myremove_pointer<Map>::type MMap ;
  //  typedef typename KNR::K RR;
    typedef Map  A;
    typedef Map  R;
    typedef E_Array B;

    class CODE : public  E_F0 { public:
        Expression a0;
        int N;
        Expression * tab;
        int * what;//  0  RR, 1 KN<RR>,
        const  bool mi;
        /*
         static KN_<RR> &set(KN<RR> * a,KN<RR> *& p,int n){
            if(isinit) a->init(n);
            else a->resize(n);
            p =a;
            return *a;}*/

       /* static KN_<RR> &set(KN_<RR> & a,KN<RR> *& p,int n){p=0;return a;}*/

        CODE(Expression a,const E_Array & tt)
        : a0(a),N(tt.size()),
        tab(new Expression [N]),
        what(new int[N])  ,
        mi(tt.MeshIndependent())

        {
            assert(&tt);
            //      int err=0;
            for (int i=0;i<N;i++)
                if(i%2==0 && atype<Key>()->CastingFrom(tt[i].right() ) )
                {
                    tab[i]=atype<Key>()->CastTo(tt[i]);
                    what[i]=1;
                }
               else if(i%2==1 && atype<Value>()->CastingFrom(tt[i].right() ) )
            {
                tab[i]=atype<Value>()->CastTo(tt[i]);
                what[i]=1;
            }

                 else
                    CompileError(" InitMapfromArray: we are waiting for Key or Value  type");
        }
        AnyType operator()(Stack stack)  const
        {

            A  aa=GetAny<A>((*a0)(stack));
            if(isinit) aa->init();
            for(int i=0,j=0; i<N/2; ++i)
            {
                Key k= GetAny<Key>((*(tab[j++]))(stack));
              //  String sk(*k);
                Value v= GetAny<Value>((*(tab[j++]))(stack));

                //cout << "InitMapfromArray  "<< *k << " " << (string) sk << " "<<  v << endl;
                aa->insert(*k,v);
            }
          return SetAny<R>(aa);
        }
        bool MeshIndependent() const     {return  mi;} //
        ~CODE() { delete [] tab; delete[] what;}
        operator aType () const { return atype<R>();}
    }; // end sub class CODE


public:
    E_F0 * code(const basicAC_F0 & args) const
    {   if(verbosity>9999)
        cout << "\n code InitMapfromArray:" << *args[0].left() << " " << *args[1].left()
        << "( "<< *t[0] << " " << *t[1] <<")" <<endl;
        return  new CODE(t[0]->CastTo(args[0]),*dynamic_cast<const E_Array*>( t[1]->CastTo(args[1]).LeftValue()));}
    InitMapfromArray(int preff=0):   OneOperator(atype<R>(),atype<A>(),atype<B>())  {
        pref=preff;
    }

};


template<class CR, class KNRR,bool isinit>
class  InitArrayfromArray : public OneOperator {
public:
    typedef typename  myremove_pointer<KNRR>::type KNR ;
    typedef typename KNR::K RR;
     typedef KNRR  A;
    typedef KNRR  R;
    typedef E_Array B;

    class CODE : public  E_F0 { public:
       Expression a0;
       int N;
       Expression * tab;
    int * what;//  0  RR, 1 KN<RR>,
       const  bool mi;

    static KN_<RR> &set(KN<RR> * a,KN<RR> *& p,int n){
                if(isinit) a->init(n);
                else a->resize(n);
        p =a;
        return *a;}

    static KN_<RR> &set(KN_<RR> & a,KN<RR> *& p,int n){p=0;return a;}

    CODE(Expression a,const E_Array & tt)
      : a0(a),N(tt.size()),
	tab(new Expression [N]),
	what(new int[N])  ,
	mi(tt.MeshIndependent())

      {
        assert(&tt);
	//      int err=0;
        for (int i=0;i<N;i++)
	if(atype<CR>()->CastingFrom(tt[i].right() ) )
	  {
          tab[i]=atype<CR>()->CastTo(tt[i]);
	    what[i]=0;
	  }
	else if(atype<KN_<RR> >()->CastingFrom(tt[i].right() ) )
	  {
	    tab[i]=atype<KN_<RR> >()->CastTo(tt[i].RightExp());
	    what[i]=1;
	  }
	else
	  CompileError(" InitArrayfromArray: we are waiting for scalar or vector of scalar");
    }
    AnyType operator()(Stack stack)  const
    {
        // a verifier ... FH nov 2015.....
	//extern void xxxx();
	//xxxx();
      KN<RR> * pa=0;
      A  aa=GetAny<A>((*a0)(stack));
       KN<AnyType> v(N);
      KN<int>  nn(N+1);
      for (int i=0;i<N;i++)
        v[i]= (*(tab[i]))(stack);

      int n=0;
      for (int i=0;i<N;i++)
	{
	  if (what[i]==0) nn[i]=1;
	  else if (what[i]==1) nn[i]=GetAny<KN_<RR> >(v[i]).size();
          n += nn[i];
	}
        if(verbosity>10000)cout << " InitArrayfromArray  aa = " <<aa << " n="<< n << endl;
        KN_<RR> a =set(aa,pa,n);
        if(verbosity>10000) cout << " InitArrayfromArray a.N() "<< a.N() << " " << n << " "
                             << pa << " " << isinit  <<endl;
        ffassert(a.N()>=n);
      for (int i=0,j=0 ;i<N; j += nn[i++])
      {
       // cout << " ### " << i << " " << j << endl;
        if (what[i]==0)
          a[j]= GetAny<CR>(v[i]);
        else if (what[i]==1)
          a(SubArray(nn[i],j)) = GetAny<KN_<RR> >((*(tab[i]))(stack));// correct bug nov 2014
      }
        //  (due to resize=> pointer  change Fh
      return SetAny<R>(aa);
    }
    bool MeshIndependent() const     {return  mi;} //
    ~CODE() { delete [] tab; delete[] what;}
    operator aType () const { return atype<R>();}
  }; // end sub class CODE


    public:
    E_F0 * code(const basicAC_F0 & args) const
    {   if(verbosity>9999)
          cout << "\n code InitArrayfromArray:" << *args[0].left() << " " << *args[1].left()
               << "( "<< *t[0] << " " << *t[1] <<")" <<endl;
         return  new CODE(t[0]->CastTo(args[0]),*dynamic_cast<const E_Array*>( t[1]->CastTo(args[1]).LeftValue()));}
    InitArrayfromArray(int preff=0):   OneOperator(atype<R>(),atype<A>(),atype<B>())  {
        pref=preff;
       // cout << "\n @@@ R " << *atype<R>()<< " A " << *atype<A>() << " B " << * atype<B>() << " " << preff <<endl;
    }

};
// correct compile PB 12/03/24 FH
template<bool t> struct InitMatfromAArray_Barray{ typedef E_Array Btype; };
template<> struct InitMatfromAArray_Barray<true>{ typedef TransE_Array Btype; };

template<class RR,bool isinit,bool Trans=false>
class  InitMatfromAArray : public OneOperator {
public:
    typedef KNM<RR> * A;
    typedef KNM<RR> * R;
    

    typedef typename InitMatfromAArray_Barray<Trans>::Btype B;

    class CODE : public  E_F0 { public:
       Expression a0;
       int N;
       int M;
       Expression ** tab;
       const  bool mi;

    CODE(Expression a,const E_Array & tt)
      : a0(a),N(tt.size()),M(0),
	tab(new Expression* [N]),
	mi(tt.MeshIndependent())

      {
        assert(&tt);
	//       int err=0;
        for (int i=0;i<N;i++)
         {
          const E_Array  *li =  dynamic_cast<const E_Array *>(tt[i].LeftValue());
          if (li)
	  {
	     const E_Array & lli = *li;
	     // -- check ---
	     if( i == 0) {
	         M = lli.size(); ffassert( M>0 );
	        for (int i=0;i<N;i++) tab[i] = new Expression [M];
	       }
	     else {
	        if ( M != li->size() ) {
	        cout << " line " << i << " the size of the column change " << M << " to " << li->size() << endl;
	        CompileError(" Is not a matrix,  M is not constant" ); } }

	    for (int j=0;j<M;j++)
              tab[i][j]=atype<RR>()->CastTo(  lli[j]);
	   }
	 else  // li == 0
	  CompileError(" we are waiting for  vector of scalar [  , , ,  ] ");
	 }
        
    }

    AnyType operator()(Stack stack)  const
    {
      A  a=GetAny<A>((*a0)(stack));
        int NN=N,MM=M;
        if(Trans) swap(NN,MM);
      if (isinit)
        a->init(NN,MM);
      else
	a->resize(NN,MM);

       for (int i =0;i<N;++i)
       for (int j =0;j<M;++j)
           if( Trans)
              (*a)(j,i)=  RNM::conj(GetAny< RR >( (*(tab[i][j]))(stack))) ;
            else
             (*a)(i,j)=   GetAny< RR >( (*(tab[i][j]))(stack)) ;
      return SetAny<R>(a);
    }
    bool MeshIndependent() const     {return  mi;} //
    ~CODE() { for (int i=0;i<N;i++) delete [] tab[i]; delete [] tab; }
    operator aType () const { return atype<R>();}
  }; // end sub class CODE


    public:
    E_F0 * code(const basicAC_F0 & args) const
     {
         const E_Array * ea;
         if(Trans)
         {
             const TransE_Array *tea =dynamic_cast<const TransE_Array*>( t[1]->CastTo(args[1]).LeftValue());
             ea =tea->v;
         }
         else
          ea =dynamic_cast<const E_Array*>( t[1]->CastTo(args[1]).LeftValue());
        return  new CODE(t[0]->CastTo(args[0]),*ea) ;}
    InitMatfromAArray():   OneOperator(atype<R>(),atype<A>(),atype<B>())  {}

};

template<typename RR>
class  SetArrayofKNfromKN : public OneOperator {
public:
    typedef KN_<RR>  A; // Warning  B type of  1 parameter
    typedef KN_<RR>  R;
    typedef E_Array B; //   A type of 2 parameter

    class CODE : public  E_F0 { public:
       Expression a0;
       int N;
       Expression * tab;
       int * what;//  0  RR, 1 KN<RR>,
       const  bool mi;

    CODE(Expression a,const E_Array & tt)
      : a0(a),N(tt.size()),
	tab(new Expression [N]),
	what(new int[N])  ,
	mi(tt.MeshIndependent())
      {
        assert(&tt);
	//      int err=0;
        for (int i=0;i<N;i++)
	if(atype<RR*>()->CastingFrom(tt[i].left() ) )
	  {
          tab[i]=atype<RR*>()->CastTo(tt[i]);
	    what[i]=0;
	  }
	else if(atype<KN_<RR> >()->CastingFrom(tt[i].right() ) )
	  {
	    tab[i]=atype<KN_<RR> >()->CastTo(tt[i].RightExp());
	    what[i]=1;
	  }
	else
	  CompileError("SetArrayofKNfromKN: we are waiting for scalar or vector of scalar");
    }

    AnyType operator()(Stack stack)  const
    {
      A  a=GetAny<A>((*a0)(stack));
      KN<AnyType> v(N);
      KN<int>  nn(N+1);
      for (int i=0;i<N;i++)
        v[i]= (*(tab[i]))(stack);

      int n=0;
      for (int i=0;i<N;i++)
	{
	  if (what[i]==0) nn[i]=1;
	  else if (what[i]==1) nn[i]=GetAny<KN_<RR> >(v[i]).size();
          n += nn[i];
	}
      ffassert(n == a.size());
      for (int i=0,j=0 ;i<N; j += nn[i++])

        if (what[i]==0)
         * GetAny<RR*>(v[i]) = a[j];
        else if (what[i]==1) { // hack FH
           KN_<RR> tab(GetAny<KN_<RR> >(v[i]));
           tab  =a(SubArray(nn[i],j));
           }
      return SetAny<R>(a);
    }
    bool MeshIndependent() const     {return  mi;} //
    ~CODE() { delete [] tab; delete[] what;}
    operator aType () const { return atype<R>();}
  }; // end sub class CODE


    public: // warning hack  A and B
    E_F0 * code(const basicAC_F0 & args) const
     { return  new CODE(t[1]->CastTo(args[1]),*dynamic_cast<const E_Array*>( t[0]->CastTo(args[0]).RightValue()));}
    SetArrayofKNfromKN():   OneOperator(atype<R>(),atype<B>(),atype<A>())  {} // warning with A and B

};


template<class K> long get_MyMap_n(MyMap<String,K> *p) {return p ? (p->m ? p->m->size():0):0 ;} 
template<class K> long get_n(KN<K> * p){ return p->N();}//
template<class K> long get__n(KN_<K>  p){ return p.N();}//

template<class K> long get_n(KNM<K> * p){ return p->N();}//
template<class K> long get_m(KNM<K> * p){ return p->M();}//
template<class K> long get__n(KNM_<K> p){ return p.N();}//
template<class K> long get__m(KNM_<K> p){ return p.M();}//

template<class K> K get_max(KN<K> * p){ return p->max();}//
template<class K> K get_min(KN<K> * p){ return p->min();}//

template<class K> long get_imax(KN<K> * p){ long  i =0; for(long k=1;k< p->N(); ++k) if(  (*p)[k]>(*p)[i] ) i=k;  return i ;}
template<class K> long get_imin(KN<K> * p){ long  i =0; for(long k=1;k< p->N(); ++k) if(  (*p)[k]<(*p)[i] ) i=k;  return i ;}
template<class K> long get_imin(KNM<K> * p){
    long  i =0,j=0;
         for(long k=0;k< p->N(); ++k)
             for(long l=0;l< p->M(); ++l) if(  (*p)(k,l) <(*p)(i,j) )  i=k,j=l;
        return i ;}

template<class K> long get_jmin(KNM<K> * p){
    long  i =0,j=0;
    for(long k=0;k< p->N(); ++k)
        for(long l=0;l< p->M(); ++l) if(  (*p)(k,l) <(*p)(i,j) )  i=k,j=l;
    return j ;}

template<class K> long get_imax(KNM<K> * p){
    long  i =0,j=0;
    for(long k=0;k< p->N(); ++k)
        for(long l=0;l< p->M(); ++l) if(  (*p)(k,l) >(*p)(i,j) )  i=k,j=l;
    return i ;}
/*template<class K> K get_max(KNM<K> * p){
    K m = (*p)(0,0), mm;
    for(long k=0;k< p->N(); ++k)
        for(long l=0;l< p->M(); ++l) if(  (mm=(*p)(k,l)) > m ) m=mm ;
    return m ;}
template<class K> K get_min(KNM<K> * p){
    K m = (*p)(0,0), mm;
    for(long k=0;k< p->N(); ++k)
        for(long l=0;l< p->M(); ++l) if(  (mm=(*p)(k,l)) < m ) m=mm ;
    return m ;}
*/
template<class K> long get_jmax(KNM<K> * p){
    long  i =0,j=0;
    for(long k=0;k< p->N(); ++k)
        for(long l=0;l< p->M(); ++l) if(  (*p)(k,l) >(*p)(i,j) )  i=k,j=l;
    return j ;}
template<class K> NothingType get_ijmax(KNM<K> * const & p,long * const & pi,long * const & pj ){
    long  i =0,j=0;
    for(long k=0;k< p->N(); ++k)
        for(long l=0;l< p->M(); ++l) if(  (*p)(k,l) >(*p)(i,j) )  i=k,j=l;
    *pi = i;
    *pj = j;
    return NothingType() ;}
template<class K> NothingType get_ijmin(KNM<K> * const & p,long * const & pi,long * const & pj )
{
    long  i =0,j=0;
    for(long k=0;k< p->N(); ++k)
        for(long l=0;l< p->M(); ++l) if(  (*p)(k,l) <(*p)(i,j) )  i=k,j=l;
    *pi = i;
    *pj = j;
    return NothingType() ;}

template<class K> K get_sum(KN<K> * p){ return p->sum();}
template<class K> double get_l2(KN<K> * p){ return p->l2();}
template<class K> double get_l1(KN<K> * p){ return p->l1();}
template<class K> double get_linfty(KN<K> * p){ return p->linfty();}

template<class K> K get_max(KNM<K> * p){ return p->max();}
template<class K> K get_min(KNM<K> * p){ return p->min();}
template<class K> K get_sum(KNM<K> * p){ return p->sum();}
template<class K> double get_l2(KNM<K> * p){ return p->l2();}
template<class K> double get_l1(KNM<K> * p){ return p->l1();}
template<class K> double get_linfty(KNM<K> * p){ return p->linfty();}

template<class K,class T > K get_sum0(const T & p){ return p.sum();}
template<class K,class T > K get_max0(const T &p){ return p.max();}
template<class K,class T > K get_min0(const T &p){ return p.min();}
template<class K,class T> K  get_l2_0(const T &p){ return p.l2();}
template<class K,class T> K  get_l1_0(const T &p){ return p.l1();}
template<class K,class T> K  get_linfty_0(const T &p){ return p.linfty();}



 class ostream_precis { public:
 ostream_precis(ostream * ff) :f(ff) {}
  ostream * f;
   operator long () const {return f->precision();}
 };

template<class A,class B> B castto(const A & a){ return a;}

/*
template<class K>
AnyType ClearReturnpKN(Stack stack, const AnyType & a)
{
    KN<K> * m = GetAny<K>(a);
    Add2StackOfPtr2FreeRC(stack, (K*) (*m) );
    if(verbosity>1)
	cout << "AddIncrement:: increment + Add2StackOfPtr2FreeRC " << endl;
    return new KN<K>(true, *m);
}*/

template<typename K,typename KK>
AnyType ClearReturnpKK(Stack stack, const AnyType & a)
{
    // a ne faire que pour les variables local au return...
    //  pour l'instant on copie pour fqire mqrche
    // a repense  FH  mqi 2009....
    KK * m = GetAny<KK * >(a);
  //   KN<K> *cm=new KN<K>(true, *m); bug quant KN est une variable global
   // KN<K> *cm=new KN<K>( *m); // on duplique le tableau comme en C++  (dur dur ?????? FH)
    m->increment();
    Add2StackOfPtr2FreeRC(stack,m);
    if(verbosity>400)
	cout << "ClearReturnpKK:: increment + Add2StackOfPtr2FreeRC nb ref  " <<  -m->next  << endl;
    return m;
}
template<typename K,typename KK,typename KK_>
AnyType ClearReturnpKK_(Stack stack, const AnyType & a)
{
   // il faut faire un copie du tableau
    KK_ * m = GetAny<KK_ * >(a);
    KK *cm=new KK(*m);

    Add2StackOfPtr2Free(stack,cm);// detruire la copie
    if(verbosity>400)
	cout << "ClearReturnpKK_:: copie  Add2StackOfPtr2Free "  << endl;
    return (KK_ *) cm;
}
template<typename K,typename KK,typename KK_>
AnyType ClearReturnKK_(Stack stack, const AnyType & a)
{
    // il faut faire un copie du tableau
    KK_  m = GetAny<KK_>(a);
    KK *cm=new KK(m);

    Add2StackOfPtr2Free(stack,cm);// detruire la copie
    if(verbosity>400)
	cout << "ClearReturnKK_:: copie  Add2StackOfPtr2Free   "  << endl;
    return SetAny<KK_>(*cm);
}
template<typename K,typename KK_,typename KK>
AnyType CopieKK_pKK(Stack stack,const AnyType &a) {
    KK_  m = GetAny<KK_>(a);
    KK *cm=new KK(m);
    if(verbosity>400)
	cout << "CopieKK_pKK:: copie  Add2StackOfPtr2Free   "<< cm   << endl;
    Add2StackOfPtr2Free(stack,cm);// detruire la copie
return cm;}


template<typename KK,typename KK_>
AnyType UnRefpKN(Stack,const AnyType &a) {
    KK_ a_(*PGetAny<KK >(a));
    return  SetAny<KK_ >(a_);}

template<class A> inline AnyType DestroyKN(Stack,const AnyType &x){
    KN<A> * a=PGetAny<KN<A> >(x);
    SHOWVERB(cout << "DESTROY " <<typeid(A).name() << " " << a <<  endl);
    for(int i=0;i<a->N(); ++i)
        (*a)[i].destroy();
    a->destroy();
    return  Nothing;
}

template<typename K>
AnyType TransposeKNM(Stack,const AnyType &a) {
    Transpose< KNM<K> * > ta(GetAny<Transpose< KNM<K> * > >(a));
    return  SetAny<KNM_<K> >((*ta.t).t());}

template<class K>
void ArrayDCL()
{
  //  Dcl_TypeandPtr<KN<K> >(0,0,0,::Destroy<KN<K> >, 0 ,  ::ClearReturnKN<K> );

    //Dcl_Type<KN<K> *>(0,::Destroy<KN<K> >,   ::ClearReturnpKK<K,KN<K> > );
    //Dcl_TypeandPtr<KN_<K> >(0,0,0,0,::ClearReturnKK_<K,KN<K>,KN_<K> >,::ClearReturnpKK_<K,KN<K>,KN_<K> >);
    Dcl_TypeandPtr_<KN_<K> ,KN<K>*  > (0,0,::InitP<KN<K> >,::Destroy<KN<K> >, ::ClearReturnKK_<K,KN<K>,KN_<K> >,::ClearReturnpKK<K,KN<K> >); // add init 0

    //  Dcl_Type<KN<Complex> *>(0,::Destroy<KN<Complex> >);
   // Dcl_Type<KN<K> *>(0,::Destroy<KN<K> >); // Modif 17102005
   // attention un exp KN<> * right est un KN<> et non un KN<> *

  //  Dcl_Type<KNM<K> *>(0,::Destroy<KNM<K> > ,::ClearReturnpKK<K,KNM<K> >);
    Dcl_TypeandPtr_<KNM_<K> ,KNM<K>*  > (0,0,0,::Destroy<KNM<K> >, ::ClearReturnKK_<K,KNM<K>,KNM_<K> >,::ClearReturnpKK<K,KNM<K> >);
    Dcl_Type<  KN<KNM<K> >* >(InitP<KN<KNM<K>>>,::DestroyKN<KNM<K> >,::ClearReturnpKK<KNM<K>,KN<KNM<K> > >);
    Dcl_Type<  KN<KN<K> >* >(InitP<KN<KN<K>>>,::DestroyKN<KN<K> >,::ClearReturnpKK<KN<K>,KN<KN<K> > >);

    Dcl_Type< outProduct_KN_<K>* >();
    Dcl_Type< Transpose<KN_<K> > > ();
    Dcl_Type< Transpose< KNM<K> *> >();

    Dcl_Type<Add_KN_<K> >();

    Dcl_Type<DotStar_KN_<K> >();
    Dcl_Type<DotSlash_KN_<K> >();
    Dcl_Type<Sub_KN_<K> >();
    Dcl_Type<Mulc_KN_<K> >();
    Dcl_Type<Divc_KN_<K> >();
    Dcl_Type<Mul_KNM_KN_<K> >();
    Dcl_Type<Add_Mulc_KN_<K> *>();
    Dcl_Type<if_arth_KN_<K> *>();
    // for    B(I) and B(I^-1)
    Dcl_Type<pair<KN_<K>,Inv_KN_long> *>();
    Dcl_Type<pair<KN_<K>,KN_<long> > *>();

    map_type[typeid(KN<K> * ).name()]->AddCast(
    new E_F1_funcT<KN<K>*,KN_<K> >(CopieKK_pKK<K,KN_<K>,KN<K> > )
	 );
// add  august 2009 FH  to see full  matrix as a array
    map_type[typeid(KN_<K>  ).name()]->AddCast(
						     new E_F1_funcT<KN_<K>,KNM<K>* >(UnRef<KN_<K>,KNM<K> *> ));


     map_type[typeid(KN_<K> ).name()]->AddCast(
    //   new E_F1_funcT<KN_<K>,KN_<K>*>(UnRefpKN_<K> ),
       new E_F1_funcT<KN_<K>,KN<K>*>(UnRefpKN<KN<K>,KN_<K> >  )
	//  inutil cas KN<K> est right expression de KN<K>*
//       new E_F1_funcT<KN_<K>,KN<K> >(Cast<KN_<K>,KN<K> >)

       );
    map_type[typeid(KNM_<K> ).name()]->AddCast(
					      new E_F1_funcT<KNM_<K>,KNM<K>*>(UnRefpKN<KNM<K>,KNM_<K> > ) ,
                                              new E_F1_funcT<KNM_<K>,Transpose< KNM<K> *> >(TransposeKNM<K>)
					      );

    //   ,new E_F1_funcT<KN_<K>,K>(ValueToKN_<K>),
    //   new E_F1_funcT<KN_<K>,K*>(PtrToKN_<K>)
/*
     // Ajoute FH
     map_type[typeid(KN<K> ).name()]->AddCast(
       new E_F1_funcT<KN<K>,KN<K>*>(UnRef<KN<K> >)
    //   ,new E_F1_funcT<KN_<K>,K>(ValueToKN_<K>),
    //   new E_F1_funcT<KN_<K>,K*>(PtrToKN_<K>)
       ); */
    map_type_of_map[make_pair(atype<long>(),atype<K>())]=atype<KN<K>*>(); // vector
    map_pair_of_type[make_pair(atype<long>(),atype<long>())] =atype<pair<long,long> >();
    map_type_of_map[make_pair(atype<pair<long,long> >(),atype<K>())]=atype<KNM<K>*>(); // matrix
    map_type_of_map[make_pair(atype<long>(),atype<KN_<K> >())]=atype<KN<KN<K> >*>();// tableau de tableau
    map_type_of_map[make_pair(atype<long>(),atype<KNM_<K> >())]=atype<KN<KNM<K> >*>();// tableau de matrix

    // Add FH for loop 2015
    atype<KN<K>*>()->SetTypeLoop(atype<K*>(),atype<long*>());
    atype<KNM<K>*>()->SetTypeLoop(atype<K*>(),atype<long*>(),atype<long*>());
    atype<KN_<K> >()->SetTypeLoop(atype<K*>(),atype<long*>());
    atype<KNM_<K> >()->SetTypeLoop(atype<K*>(),atype<long*>(),atype<long*>());

      // map string  of array ..  Test FH Mai 2022 ...
  
    map_type[typeid(MyMap<String,K>*).name()] = new ForEachType<MyMap<String,K>*>(Initialize<MyMap<String,K> >,Delete<MyMap<String,K> >) ;

    map_type_of_map[make_pair(atype<string*>(),atype<K>())]=atype<MyMap<String,K>*>();


    typedef MyMap<String,KN<K> > MyMapofArray;
    map_type[typeid(MyMapofArray*).name()] = new ForEachType<MyMapofArray*>(Initialize<MyMapofArray >,Delete<MyMapofArray >) ;
    ffassert(TypeArray(atype<K>(),atype<string*>()));
    // make_pair(atype<string*>(),atype<KN_<K> >())
    //  real[string][int]
    map_type_of_map[make_pair(TypeArray(atype<K>(),atype<string*>()),atype<long >())]=atype<MyMapofArray*>();
    ;// map of   tableau
  /*
    atype<MyMap<String,K>*>()->Add("[","",new OneOperator2_<K*,MyMap<String,K>*,string*>(get_element<K>));
    TheOperators->Add("&",new OneOperator2_<bool,MyMap<String,K>*,string*>(exist_element<K>));
   TheOperators->Add("<<",new OneBinaryOperator<PrintP<MyMap<String,K>*> >);
   Add<MyMap<String,K>* >("n",".",new OneOperator1<long,MyMap<String,K> *>(get_MyMap_n));
   */

}



template<class A,class B> pair<A,B> * pBuild(const A & a,const B & b)
  { return new pair<A,B>(a,b);}

// add mars 2006
template<class K,class L,class OP>
struct set_A_BI {
  using first_argument_type  = KN_<K>;
  using second_argument_type = pair<KN_<K>, KN_<L> > *;
  using result_type          = KN_<K>;
  static KN_<K> f(const KN_<K>   & a, pair<KN_<K>, KN_<L> > * const & b)  {
    KN_<K> x(a);
    OP op;
     const KN_<K> & y(b->first);
    const KN_<L> & I(b->second);
    L  N = x.N();
    L n = y.N();

    L maxI=I(SubArray(N)).max() ;
    L minI=I(SubArray(N)).min() ;

    if( maxI >= n || I.N()  < N)
       { cerr << " Out of Bound x=y(I)  :  0 <= " << minI << " < "<< maxI << "< " << n  << endl;
         cerr << " or I.N() " << I.N() << " > " << N << endl;
         ExecError("Out of Bound error");
       }

    for(int i=0;i<N;i++)
      if(I[i]>=0)
      op(x(i),y(I[i]));
    delete b;
    return a;

  }
};
//  add oct 2019  To:  real[int] b = a(I); // where a and I is also a array..
template<class K,class L,class OP>
struct init_A_BI {
  using first_argument_type  = KN<K>*;
  using second_argument_type = pair<KN_<K>, KN_<L> > *;
  using result_type          = KN<K>*;
    static KN<K>* f( KN<K>  * const  & a, pair<KN_<K>, KN_<L> > * const & b)  {
        KN<K> * px(a);
        OP op;
        const KN_<K> & y(b->first);
        const KN_<L> & I(b->second);
        L n = y.N();
        px->init(I.N());
        KN<K> & x=*px;
        L  N = x.N();

        L maxI=I.max() ;
        L minI=I.min() ;
        
        if( maxI >= n )
        { cerr << " Out of Bound x=y(I)  :  0 <= " << minI << " < "<< maxI << "< " << n  << endl;
            cerr << " or I.N() " << I.N() << " > " << N << endl;
            ExecError("Out of Bound error");
        }
        
        for(int i=0;i<N;i++)
            if(I[i]>=0)
                op(x(i),y(I[i]));
        delete b;
        return a;
        
    }
};
template<class K,class L,class OP>
struct set_AI_B {
  using first_argument_type  = pair<KN_<K>, KN_<L> > *;
  using second_argument_type = KN_<K>;
  using result_type          = NothingType;
  static NothingType  f( pair<KN_<K>, KN_<L> > * const & b,const KN_<K>   & a)  {
    KN_<K> x(a);
    OP op;
     const KN_<K> & y(b->first);
    const KN_<L> & I(b->second);
    L  N = I.N();
    L n = y.N();

    L maxI=I(SubArray(N)).max() ;
    L minI=I(SubArray(N)).min() ;

    if(  maxI >= n || x.N()  < N )
       { cerr << " Out of Bound x(I)=y  :  0 <= " << minI << " < "<< maxI << "< " << n  << endl;
         cerr << " or x.N() " << I.N() << " > " << N << endl;
         ExecError("Out of Bound error");
       }

    for(int i=0;i<N;i++)
      if(I[i] >=0)
      op(y(I[i]),x[i]);
    delete b;
    return  NothingType();

  }
};

template<class K>
struct Op3_paac: public ternary_function<KN_<K>,KN_<K>,K,if_arth_KN_<K>*> {
static if_arth_KN_<K>* f(Stack s,const KN_<K> & a,const KN_<K> & b,const  K & c )  {
    //K cc(c);
    KN_<K> kc(new(NewAllocTmp(s,sizeof(c))) K(c),1,0);
  return new if_arth_KN_<K>(a,b,kc);}
};
template<class K>
struct Op3_paca: public ternary_function<KN_<K>,K,KN_<K>,if_arth_KN_<K>*> {
    static if_arth_KN_<K>* f(Stack s,const KN_<K> & a,const  K & b,const KN_<K> & c )  {
	//K bb(b);
	KN_<K> kb(new(NewAllocTmp(s,sizeof(b))) K(b),1,0);
    return new if_arth_KN_<K>(a,kb,c);}
};

template<class K>
struct Op3_pacc: public ternary_function<KN_<K>,K,K,if_arth_KN_<K>*> {
    static if_arth_KN_<K>* f(Stack s,const KN_<K> & a,const K & b,const  K & c )  {
	K cc(c),bb(b);
	KN_<K> kc(new(NewAllocTmp(s,sizeof(c))) K(c),1,0),
	       kb(new(NewAllocTmp(s,sizeof(b))) K(b),1,0);
    return new if_arth_KN_<K>(a,kb,kc);}
};
template<class K> KNM_<K> Transp(KNM_<K>  M){ return M.t();} // Add FH July 2015
template<class K>
struct SetArray2{
  using first_argument_type  = K;
  using second_argument_type = K;
  using result_type          = SetArray<K>;
  static SetArray<K> f(const K & a,const K & b)  {
    // cout << "SubArray: " << a << " " << b << endl;
    //     SetArray(long nn,R oo=R(),R sstep=R(1)): o(oo),n(nn),step(sstep) {}
    long  n= long(abs((b-a)));
      
      
    ffassert(n>=0);
    K s= (n==0) ? 1 : (b-a)/K(n);// correction oct 2019  FH Thanks to P. Ventura
    n++;
    if(verbosity>100)
      cout << "    SetArray " << n << " " << a << " " << s << endl;
    return SetArray<K>(n,a,s);} };

template<class K>
struct SetArray3: public ternary_function<K,K,K,SetArray<K> > {
  static SetArray<K> f(Stack s,const K & a,const K &b,const K & c)  {
    // cout << "SubArray: " << a << " " << b << " " <<  c << endl;
    long n= long(1+abs((c-a)/b));
    if(verbosity>100)
      cout << "    SetArray " << n << " :  "  << " " << a << " " << b << " " << c << endl;
    return SetArray<K>(n,a,b);} };

template<class R,class A>  R * set_init_array( R* const & a,const A & b){
    SHOWVERB( cout << " set_init " << typeid(R).name() << " " << &b << endl);
    a->init(b.size());
    *a=b;
return a;}
template<class R,class A>  R * set_array( R* const & a,const A & b){
    SHOWVERB( cout << " set_init " << typeid(R).name() << " " << &b << endl);
    a->resize(b.size());
    *a=b;
return a;}
// missing FH august 2009
template<class R,class A>  R * set_arrayp( R* const & a,const A & b){
    SHOWVERB( cout << " set_init " << typeid(R).name() << " " << &b << endl);
    a->resize(b->size());
    *a=*b;
return a;}
template<class R,class A>  R  set_array_( R const & a,const A & b){
    SHOWVERB( cout << " set_array_ " << typeid(R).name() << " " << &b << endl);
    ffassert(a.N()==b.size());
    R aa=a;
    aa=b;
return a;}
// xxxxxxxxxxx
template<class K>  KNM<K> * set_initmat_t(KNM<K> * a,Transpose<KNM<K> * > b ){
    ConjKNM_<K>  tb=b.t->t(); ;
     a->init(tb.N(),tb.M());
    *a=tb;
    return a;}
template<class K>  KNM<K> * set_initmat(KNM<K> * a,KNM<K> *  b ){

    a->init(b->N(),b->M());
    *a=*b;
    return a;}
template<class K>  KNM<K> * set_mat_t(KNM<K> * a,Transpose<KNM<K> * > b ){
    ConjKNM_<K>  tb=b.t->t(); ;
    a->resize(tb.N(),tb.M());// correction mars 2013
    *a=tb;
    return a;}
template<class K>  KNM<K> * addto_mat_t(KNM<K> * a,Transpose<KNM<K> * > b ){
    ConjKNM_<K>  tb=b.t->t(); ;
    *a+=tb;
    return a;}
template<class K>  KNM<K> * subto_mat_t(KNM<K> * a,Transpose<KNM<K> * > b ){
    ConjKNM_<K>  tb=b.t->t(); ;
    *a-=tb;
    return a;}
template<class K>  KNM<K> * set_mat(KNM<K> * a,KNM<K> *  b ){

    a->resize(b->N(),b->M());
    *a=*b;
    return a;}

template<class K,class KK=K>
class  OneOperator_2KN_ : public OneOperator {public:
    class Op : public E_F0 {
       public:
	int N;
	Expression *tab;

	Op( const  E_Array &bb) : N(bb.size()), tab(new Expression[N])
	{
	  for(int i=0;i<N;++i)
	    tab[i]=atype<KK>()->CastTo( bb[i]);
	}
	AnyType operator()(Stack s)  const {
	    K * p = Add2StackOfPtr2FreeA<K>(s,new K[N]); //   mark to be delete ..
            KN_<K> A(p,N); // FH:  Correct jan 2019  bug is: stupide A what delete if KN  type
	    for(int i=0;i<N;++i)
		A[i]= GetAny<KK>( (*tab[i])(s));
	    return SetAny<KN_<K> >(A);}
    };
    E_F0 * code(const basicAC_F0 & a) const
    {  const  E_Array * b = dynamic_cast<const E_Array *>(a[0].LeftValue());
	ffassert(b);
        return new Op(*b);}
    
    OneOperator_2KN_(): OneOperator(atype<KN_<K> >(),atype<E_Array>()) { pref=-1;}
};
template<class K, class L>
class Unique_Op : public E_F0mps {
    public:
        Expression ar;
        Expression va;
        static const int n_name_param = 1;
        static basicAC_F0::name_and_type name_param[];
        Expression nargs[n_name_param];
        Unique_Op<K, L>(const basicAC_F0& args, Expression param1, Expression param2) : ar(param1), va(param2) {
            args.SetNameParam(n_name_param, name_param, nargs);
        }
        AnyType operator()(Stack stack) const;
};
template<class K, class L>
basicAC_F0::name_and_type Unique_Op<K, L>::name_param[] = {
    {"remove", &typeid(long)},
};
template<class K, class L>
class Unique : public OneOperator {
    public:
        Unique() : OneOperator(atype<long>(), atype<KN<K>*>(), atype<KN<L>*>()) { }

        E_F0* code(const basicAC_F0& args) const {
            return new Unique_Op<K, L>(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]));
        }
};
template<class K, class L>
AnyType Unique_Op<K, L>::operator()(Stack stack) const {
    KN<K>* array = GetAny<KN<K>*>((*ar)(stack));
    KN<L>* val = GetAny<KN<L>*>((*va)(stack));
    std::set<L> vals;
    for(int i = 0; i < array->n; ++i)
    //    if(!nargs[0] || (*array)[i] != remove)
            vals.insert((*array)[i]);
    if( nargs[0]) // remove
    {
        long remove = GetAny<long>((*nargs[0])(stack)) ;
        typename std::set<L>::iterator lr = vals.find(remove);
        if( lr != vals.end()) vals.erase(lr);
    }
    val->resize(vals.size());
    int i = 0;
    for(typename std::set<L>::iterator it = vals.begin(); it != vals.end(); ++it)
        (*val)[i++] = *it;
    return static_cast<long>(vals.size());
}





template<class K>
class E_ForAllLoopRNM
{  public:
    typedef KNM_<K> Tab;
    typedef  ForAllLoopOpBase DataL;
    const DataL *data;
    E_ForAllLoopRNM(const DataL *t): data(t){}
    AnyType f(Stack s) const {
        Tab t= GetAny<KNM_<K> >(data->tab(s));
        long * i=   GetAny<long * >(data->i(s));
        long * j=  GetAny<long * >(data->j(s));
        K * v   =   GetAny<K * >(data->v(s));
        if(verbosity>1000) {
        cout << i << " " << j << " " << v << " " << data->epl <<     endl;
        cout << " i " << (char*) (void *) i -  (char*)(void*) s ;
        cout << " j " <<  (char*)(void *) j -  (char*)(void*) s ;
        cout << " vij " <<  (char*) (void *) v -  (char*)(void*) s ;
        cout << endl;
        }

        ffassert(i && j && v);
        for ( *i=0;*i<t.N();++*i)
            for ( *j=0;*j<t.M();++*j)
            {
                *v = t(*i,*j);
                data->code(s);
                t(*i,*j)= *v;
            }
      //  data->end(s);
        return Nothing  ;
     }

};

template<class K>
class E_ForAllLoopRN
{  public:
    typedef KN_<K> Tab;
    typedef  ForAllLoopOpBase DataL;
    const DataL *data;
    E_ForAllLoopRN(const DataL *t): data(t){}
    AnyType f(Stack s) const {
        Tab t= GetAny<KN_<K> >(data->tab(s));
        long * i=   GetAny<long * >(data->i(s));
        K * v   =   GetAny<K * >(data->v(s));
  //      cout << i << " " << j << " " << v << " " << data->epl <<     endl;
         if(verbosity>1000) {
        cout << " i " << (char*) (void *) i -  (char*)(void*) s ;
        cout << " vi " <<  (char*) (void *) v -  (char*)(void*) s ;
        cout << endl;
         }

        ffassert(i && v);
        for ( *i=0;*i<t.N();++*i)
            {
                *v = t[*i];
                data->code(s);
                t[*i]= *v;
            }
     //   data->end(s);
        return Nothing  ;
    }

};

template<class V>
class E_ForAllLoopMapSI
{  public:
    typedef String K;
    typedef string *KK;
    typedef V *VV;

    typedef MyMap<K,V>  *Tab;
    typedef typename  MyMap<K,V>::iterator TabI ;

    typedef  ForAllLoopOpBase DataL;
    const DataL *data;
    E_ForAllLoopMapSI(const DataL *t): data(t){}
    AnyType f(Stack s) const {
        Tab t= GetAny<Tab >(data->tab(s));
        KK * i   =   GetAny<KK*>(data->i(s));
        VV  v   =   GetAny<VV >(data->v(s));
        //      cout << i << " " << j << " " << v << " " << data->epl <<     endl;
        if(verbosity>1000) {
        cout << " i " << (char*) (void *) i -  (char*)(void*) s ;
        cout << " vi " <<  (char*) (void *) v -  (char*)(void*) s ;
        cout << endl;
        }

        ffassert(i && v);
        if(t->m)
            for (TabI ii=t->m->begin();ii != t->m->end();++ii)
            {
                String  kk = ii->first;
                V  vv = ii->second;

                *i =  kk;
                *v =  vv;
                // for  Windows otherwise trap ???? FH. march 2016
                if(verbosity>99999) cout << " " << i << " "<< v  << " "  << kk << " " <<  vv << endl;
                 data->code(s);

                ii->second  = *v;
                *i=0;
            }
      //  data->end(s);
        return Nothing  ;
    }

};
template <class K>
KN_<K> getdiag_(KNM_<K> A) {
    int n = A.N(), m=A.M();
    int nn= min(n,m);
    int si= A.shapei.step;
    //int ni =A.shapei.next;
    int sj= A.shapej.step;
    
    KN_<K> d(A,SubArray(nn,0,sj+si));
    return d;
}
template <class K>
KN_<K> getdiag(KNM<K> *pA) {
    ffassert(pA);
    return getdiag_(*pA);
}

template <class K>
KN_<K> asarray(KNM<K> *pA) {
    ffassert( pA->IsVector1());
    return *pA; }
// Add Oct 2019
template<class K, bool init>
KNM<K> * set_Eye(KNM<K> *pA,const  Eye eye)
{
    int n = eye.n, m=eye.m, nn= min(n,m);
    if( init) pA->init(n,m);
    else pA->resize(n,m);
    *pA = K();
    KN_<K> d(*pA,SubArray(nn,0,n+1));
    d= K(1.);
    return  pA;
}
extern aType aaaa_knlp;

template<class K,class Z>
void ArrayOperator()
{
    //  juin 2009  remove type KN_< > *
    // and set  KN<> * 9left expression) qnd KN_<> is the associated expression..
    // =>  lot of change because  KN<>* and KN_< > can generqte ambuguity.
    // so remove all to code with KN<>* type.
    // the remove cqde are in comment :
    //  the comment begin //-
    // and the if(0) in comment /* */


     Dcl_Type< Resize<KN<K> > > ();
     Dcl_Type< Resize<KNM<K> > > ();

   //-  typedef KN<Z> ZN;

    // add  dec 2009.  ne marche pas ( incompatible  avec MatrixBlock) a comprendre ????? FH.
    //  //   xxxxxxxxxx  2010 feb.   retest .. FH
    //   il y a plusieurs problems
    //    1)   [1,2,3.] ->  tableau de quel type  int, real , complex ????
    //
     //   map_type[typeid(KN_<K>).name()]->AddCast(new OneOperator_2KN_<K>);
    // fin add
    // ----

     atype<KN<K>* >()->Add("[","",new OneOperator2_<K*,KN<K>*,Z >(get_elementp_<K,KN<K>*,Z>));
     atype<KN<K>* >()->Add("(","",new OneOperator2_<K*,KN<K>*,Z >(get_elementp_<K,KN<K>*,Z>));
     atype<KN_<K> >()->Add("(","",new OneOperator2_<KN_<K>,KN_<K>,char >(fSubArrayc<KN_<K> >));
     atype<KN_<K> >()->Add("(","",new OneOperator2_<KN_<K>,KN_<K>,SubArray>(fSubArray<K> ));
     atype<KN<K>*>()->Add("(","",new OneOperator2_<KN_<K>,KN<K>*,SubArray>(fSubArrayp<K> ));
     atype<KN<K>* >()->Add("(","",new OneOperator2_<KN<K>*,KN<K>*,char >(fSubArrayc<KN<K>* >));
//

     atype<KNM_<K> >()->Add("(","",new OneOperator3_<KNM_<K>,KNM_<K>,SubArray,SubArray>(fSubArraybb<K> ));
     atype<KNM<K>* >()->Add("(","",new OneOperator3_<KNM_<K>,KNM<K>*,SubArray,SubArray>(fSubArraypbb<K> ));
    /*
     atype<KN_<K> >()->Add("(","",new OneOperator3_<KN_<K>,KNM_<K>,SubArray,long>(fSubArraybi<K> ));
     atype<KN_<K> >()->Add("(","",new OneOperator3_<KN_<K>,KNM_<K>,long,SubArray>(fSubArrayib<K> ));
     atype<KN_<K> >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,SubArray,long>(fSubArraypbi<K> ));
     atype<KN_<K> >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,long,SubArray>(fSubArraypib<K> ));
     */
//

    atype<KN_<K> >()->Add("[","",new OneOperator2_<K*,KN_<K>,Z >(get_element_<K,KN_<K>,Z>));
    atype<KN_<K> >()->Add("(","",new OneOperator2_<K*,KN_<K>,Z >(get_element_<K,KN_<K>,Z>));


     atype<KNM<K>* >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,Z,SubArray >(get_element_is<KN_<K>,KNM<K>*,Z,SubArray>));
     atype<KNM<K>* >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,SubArray,Z >(get_element_si<KN_<K>,KNM<K>*,SubArray,Z>));
     atype<KNM<K>* >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,Z,char >(get_element_lineorcol<KN_<K>,KNM<K>*,Z,char>));
     atype<KNM<K>* >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,char,Z >(get_element_lineorcol<KN_<K>,KNM<K>*,char,Z>));
     atype<KNM<K>* >()->Add("(","",new OneOperator3_<K*,KNM<K>*,Z,Z >(get_elementp2_<K,KNM<K>*,Z,Z>));

    atype<KNM_<K> >()->Add("(","",new OneOperator3_<KN_<K>,KNM_<K>,Z,SubArray >(get_element_is_<KN_<K>,KNM_<K>,Z,SubArray>));
    atype<KNM_<K> >()->Add("(","",new OneOperator3_<KN_<K>,KNM_<K>,SubArray,Z >(get_element_si_<KN_<K>,KNM_<K>,SubArray,Z>));
    atype<KNM_<K> >()->Add("(","",new OneOperator3_<KN_<K>,KNM_<K>,Z,char >(get_element_lineorcol_<KN_<K>,KNM_<K>,Z,char>));
    atype<KNM_<K> >()->Add("(","",new OneOperator3_<KN_<K>,KNM_<K>,char,Z >(get_element_lineorcol_<KN_<K>,KNM_<K>,char,Z>));
    atype<KNM_<K> >()->Add("(","",new OneOperator3_<K*,KNM_<K>,Z,Z >(get_elementp2__<K,KNM_<K>,Z,Z>));


     Add<KN<K> *>("sum",".",new OneOperator1<K,KN<K> *>(get_sum));
     Add<KN<K> *>("min",".",new OneOperator1<K,KN<K> *>(get_min));
     Add<KN<K> *>("max",".",new OneOperator1<K,KN<K> *>(get_max));

     Add<KN<K> *>("l2",".",new OneOperator1<double,KN<K> *>(get_l2));
     Add<KN<K> *>("l1",".",new OneOperator1<double,KN<K> *>(get_l1));
     Add<KN<K> *>("linfty",".",new OneOperator1<double,KN<K> *>(get_linfty));
// add july 2009
    Add<KNM<K> *>("sum",".",new OneOperator1<K,KNM<K> *>(get_sum));
    Add<KNM<K> *>("min",".",new OneOperator1<K,KNM<K> *>(get_min));
    Add<KNM<K> *>("max",".",new OneOperator1<K,KNM<K> *>(get_max));
    Add<KNM<K> *>("l2",".",new OneOperator1<double,KNM<K> *>(get_l2));
    Add<KNM<K> *>("l1",".",new OneOperator1<double,KNM<K> *>(get_l1));
    Add<KNM<K> *>("linfty",".",new OneOperator1<double,KNM<K> *>(get_linfty));
// end add

     Add<KN_<K> >("sum",".",new OneOperator1_<K,KN_<K> >(get_sum0<K,KN_<K> >));
     Add<KN_<K> >("min",".",new OneOperator1_<K,KN_<K> >(get_min0<K,KN_<K> >));
     Add<KN_<K> >("max",".",new OneOperator1_<K,KN_<K> >(get_max0<K,KN_<K> >));
     Add<KN_<K> >("l2",".",new OneOperator1_<double,KN_<K> >(get_l2_0<double,KN_<K> >));
     Add<KN_<K> >("l1",".",new OneOperator1_<double,KN_<K> >(get_l1_0<double,KN_<K> >));
     Add<KN_<K> >("linfty",".",new OneOperator1_<double,KN_<K> >(get_linfty_0<double,KN_<K> >));

// add july 2009
    Add<KNM_<K> >("sum",".",new OneOperator1_<K,KNM_<K> >(get_sum0<K,KNM_<K> >));
    Add<KNM_<K> >("min",".",new OneOperator1_<K,KNM_<K> >(get_min0<K,KNM_<K> >));
    Add<KNM_<K> >("max",".",new OneOperator1_<K,KNM_<K> >(get_max0<K,KNM_<K> >));
    Add<KNM_<K> >("l2",".",new OneOperator1_<double,KNM_<K> >(get_l2_0<double,KNM_<K> >));
    Add<KNM_<K> >("l1",".",new OneOperator1_<double,KNM_<K> >(get_l1_0<double,KNM_<K> >));
    Add<KNM_<K> >("linfty",".",new OneOperator1_<double,KNM_<K> >(get_linfty_0<double,KNM_<K> >));
// end add


/*
     Add<KN<K> >("sum",".",   new OneOperator1_<K,KN<K> >(get_sum0<K,KN<K> >));
     Add<KN<K> >("min",".",   new OneOperator1_<K,KN<K> >(get_min0<K,KN<K> >));
     Add<KN<K> >("max",".",   new OneOperator1_<K,KN<K> >(get_max0<K,KN<K> >));
     Add<KN<K> >("l2",".",    new OneOperator1_<double,KN<K> >(get_l2_0<double,KN<K> >));
     Add<KN<K> >("l1",".",    new OneOperator1_<double,KN<K> >(get_l1_0<double,KN<K> >));
     Add<KN<K> >("linfty",".",new OneOperator1_<double,KN<K> >(get_linfty_0<double,KN<K> >));
*/

     Add<KN<K> *>("resize",".",new OneOperator1< Resize<KN<K> >,KN<K> *>(to_Resize));
     Add<KNM<K> *>("resize",".",new OneOperator1< Resize<KNM<K> >,KNM<K> *>(to_Resize));

     Add<Resize<KN<K> > >("(","",new OneOperator2_<KN<K> *,Resize<KN<K> > , Z   >(resize1));
     Add<Resize<KNM<K> > >("(","",new OneOperator3_<KNM<K> *,Resize<KNM<K> > , Z, Z  >(resize2));


     TheOperators->Add("<-",
 //      new OneOperator2_<KN<K> *,KN<K> *,R3>(&set_initR3),
 //      new OneOperator2_<KN<K> *,KN<K> *,R3*>(&set_initR3),

       new OneOperator2_<KN<K> *,KN<K> *,Z>(&set_init),
       new InitArrayfromArray<K,KN<K>*,true>
    //   new OneOperator2_<KN<K> *,KN<K> *,KN<K> >(&set_init),
    //   new OneOperator2_<KN<K> *,KN<K> *,KN_<K> >(&set_init)		????
     //  new OneOperator2_<KN<K> *,KN<K> *,KN<K> * >(&set_initp)
       );
    TheOperators->Add("<-",
                      new OneOperator2_<KN< KN<K> > *,KN< KN<K> > * ,Z  >(&set_init));
    TheOperators->Add("<-",
                      new OneOperator2_<KN< KNM<K> > *,KN< KNM<K> > * ,Z  >(&set_init));

     TheOperators->Add("<-",
        new OneOperator3_<KNM<K> *,KNM<K> *,Z,Z>(&set_init2),
        new InitMatfromAArray<K,true,false>,
        new InitMatfromAArray<K,true,true>
       );

     Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,Z>(&set_init));
   //  Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,KN<K> >(&set_init));
     //Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,KN_<K> >(&set_init));
    // Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,KN<K> * >(&set_initp));
     Add<KNM<K> *>("<-","(",new OneOperator3_<KNM<K> *,KNM<K> *,Z,Z>(&set_init2));
     TheOperators->Add("<-",new OneOperator2<KNM<K> *,KNM<K> *,Transpose<KNM<K> * > >(&set_initmat_t));// may 2011 FH..
    TheOperators->Add("=",new OneOperator2<KNM<K> *,KNM<K> *, Eye  >(set_Eye<K,false>));// may 2011 FH..
    TheOperators->Add("<-",new OneOperator2<KNM<K> *,KNM<K> *, Eye  >(set_Eye<K,true>));// may 2011 FH..
    TheOperators->Add("=",new OneOperator2<KNM<K> *,KNM<K> *,Transpose<KNM<K> * > >(&set_mat_t));// may 2011 FH..
    TheOperators->Add("+=",new OneOperator2<KNM<K> *,KNM<K> *,Transpose<KNM<K> * > >(&addto_mat_t));// oct 2011 FH..
    TheOperators->Add("-=",new OneOperator2<KNM<K> *,KNM<K> *,Transpose<KNM<K> * > >(&subto_mat_t));// oct 2011 FH..
//     TheOperators->Add("-=",new OneOperator2<KNM<K> *,KNM<K> *,Transpose<KNM<K> * > >(&subto_mat_t<-1>));// oct 2011 FH..
   //  Add<KNM<K> *>("<-","(",new OneOperator2<KNM<K> *,KNM<K> *,KNM<K> *  >(&set_initmat));// may 2011 FH..
   //  Add<KNM<K> *>("=","(",new OneOperator2<KNM<K> *,KNM<K> *,Transpose<KNM<K> * > >(&set_mat_t));// may 2011 FH..
  //   Add<KNM<K> *>("=","(",new OneOperator2<KNM<K> *,KNM<K> *,KNM<K> *  >(&set_mat));// may 2011 FH..

    // Add<KNM<K> *>("=","(",new OneOperator2_<KNM<K> *,KNM<K> *,Transpose<KNM<K> * > >(&set_tt));

     Add<KN<K> *>("<-","(",new InitArrayfromArray<K,KN<K>*,true>);
     Add<KNM<K> *>("<-","(",new InitMatfromAArray<K,true,false>);
     Add<KNM<K> *>("<-","(",new InitMatfromAArray<K,true,true>);
     Add<KN<K> *>("n",".",new OneOperator1<Z,KN<K> *>(get_n));
     Add<KN_<K> >("n",".",new OneOperator1<Z,KN_<K> >(get__n));
     Add<KNM<K> *>("n",".",new OneOperator1<Z,KNM<K> *>(get_n));
     Add<KNM<K> *>("m",".",new OneOperator1<Z,KNM<K> *>(get_m));
    Add<KNM_<K> >("n",".",new OneOperator1<Z,KNM_<K> >(get__n));
    Add<KNM_<K> >("m",".",new OneOperator1<Z,KNM_<K> >(get__m));
 //ajout ars 2012 FH
     Add<KN<KN<K> > *>("n",".",new OneOperator1<long,KN<KN<K> > *>(get_n));
     Add<KN<KNM<K> > *>("n",".",new OneOperator1<long,KN<KNM<K> > *>(get_n));
     atype<KN<KN<K> > * >()->Add("[","",new OneOperator2_<KN<K>*,KN<KN<K> >*,Z >(get_elementp_<KN<K>,KN<KN<K> >*,Z>));
    atype<KN<KNM<K> > * >()->Add("[","",new OneOperator2_<KNM<K>*,KN<KNM<K> >*,Z >(get_elementp_<KNM<K>,KN<KNM<K> >*,Z>));
    Dcl_Type< Resize<KN<KN<K> > > > ();
    Dcl_Type< Resize<KN<KNM<K> > > >();
    Add<KN<KN<K> > * >("resize",".",new OneOperator1< Resize<KN<KN<K> > >,KN<KN<K> > *>(to_Resize));
    Add<KN<KNM<K> > * >("resize",".",new OneOperator1< Resize<KN<KNM<K> > >,KN<KNM<K> > *>(to_Resize));
    Add<Resize<KN<KN<K> > > >("(","",new OneOperator2_<KN<KN<K> >  *,Resize<KN<KN<K> > > , long   >(resize1));
    Add<Resize<KN<KNM<K> > > >("(","",new OneOperator2_<KN<KNM<K> >  *,Resize<KN<KNM<K> > > , long   >(resize1));
    Dcl_Type<Mul_KNMh_KN_<K> >();
//     AddOpeqarray<set_eqarray,KN,K>("=");

     TheOperators->Add("=", new InitArrayfromArray<K,KN<K>*,false>(10));
     TheOperators->Add("=", new InitArrayfromArray<K,KN_<K>,false>(1));// ???????? FH nov 2015 ..
     TheOperators->Add("=", new InitMatfromAArray<K,false,false>,
                            new InitMatfromAArray<K,false,true>
       );
     TheOperators->Add("=", new SetArrayofKNfromKN<K>
       );
 if(0) //  a  change il faut regle un PB ambiguite ...
     TheOperators->Add("=",
        new OneBinaryOperator<set_eqarray<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Divc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,KN_<K> > > , // Add FH juin 2005
        new OneBinaryOperator<set_eqarraypd<KN<K> ,Add_Mulc_KN_<K>* > > , // Add FH aug 2005
        new OneBinaryOperator<set_eqarraypd<KN<K> ,if_arth_KN_<K>* > >
      //  new OneBinaryOperator<set_eqarrayp<KN<K> ,KN<K>* > >   // test aug 2009
      );
  // add august 2007

     TheOperators->Add("<-",
		      // new OneBinaryOperator<set_eqarray<KN<K> ,K > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,Add_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,DotStar_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,DotSlash_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,Sub_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,Mulc_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,Divc_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,Mul_KNM_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,KN_<K> > > , // Add FH juin 2005
		       new OneBinaryOperator<init_eqarraypd<KN<K> ,Add_Mulc_KN_<K>* > > , // Add FH aug 2005
		       new OneBinaryOperator<init_eqarraypd<KN<K> ,if_arth_KN_<K>* > >
		      // new OneBinaryOperator<init_eqarrayp<KN<K> ,KN<K>* > >
		       );


     TheOperators->Add("=",
        new OneBinaryOperator<set_eqarray<KNM<K>  ,K > > ,
         new OneBinaryOperator<set_eqarrayp<KNM<K>  , KNM<K> *  > >

       );

     TheOperators->Add("=",
        new OneBinaryOperator<set_eq_array<KN_<K> ,K > > ,
        new OneBinaryOperator<set_eq_array<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array<KN_<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array<KN_<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array<KN_<K> ,Divc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array<KN_<K> ,Mul_KNM_KN_<K> > > ,
	new OneBinaryOperator<set_eq_arraypd<KN_<K> ,if_arth_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arraypd<KN_<K> ,Add_Mulc_KN_<K>* > >  , // Add FH aug 2005
	new OneBinaryOperator<set_eq_array<KN_<K> ,KN_<K> > >, // add FH juin 2005
        new OneBinaryOperator<set_eq_arraypd<KN_<K> ,KN<K>* > >

      //-  new OneBinaryOperator<set_eq_arrayp<KN_<K> ,KN<K>* > >
      );
//  ajoute mars 2010  FH
    TheOperators->Add("<-",
		      new OneBinaryOperator<init_eqarray<KNM<K> ,KNM_<K> > >
		      );

    TheOperators->Add("=",
		      new OneBinaryOperator<set_eqarray<KNM<K>  ,KNM_<K> > >

		    //  new OneBinaryOperator<set_eq_array<KNM_<K> ,K > > ,
		    //  new OneBinaryOperator<set_eq_array<KNM_<K> ,KNM_<K> > >,
		    //  new OneBinaryOperator<set_eq_arraypd<KNM_<K> ,KNM<K>* > >
		      );

//  end add ...
/*if(0)
     TheOperators->Add("+=",
        new OneBinaryOperator<set_eqarray_add<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_add<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarraypd_add<KN<K> ,if_arth_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,KN_<K> > > , // Add FH juin 2005
        new OneBinaryOperator<set_eqarrayp_add<KN<K> ,KN<K>* > >

      );
*/
     TheOperators->Add("+=",
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,K > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Divc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_add<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arraypd_add<KN_<K> ,if_arth_KN_<K>* > > ,
	new OneBinaryOperator<set_eq_array_add<KN_<K> ,KN_<K> > >  // add FH juin 2005

       // new OneBinaryOperator<set_eq_arrayp_add<KN_<K> ,KN<K>* > >
      );
/*    if(0)
     TheOperators->Add("-=",
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_sub<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarraypd_sub<KN<K> ,if_arth_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,KN_<K> > > , // Add FH juin 2005
        new OneBinaryOperator<set_eqarrayp_sub<KN<K> ,KN<K>* > >
      );*/

     TheOperators->Add("-=",
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,K > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Divc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_sub<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arraypd_sub<KN_<K> ,if_arth_KN_<K>* > > ,
       //- new OneBinaryOperator<set_eq_arrayp_sub<KN_<K> ,KN<K>* > >
	new OneBinaryOperator<set_eq_array_sub<KN_<K> ,KN_<K> > >   // Add FH juin 2005
      );

/*    if(0)
    TheOperators->Add("*=",
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,K > >  ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_mul<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarrayp_mul<KN<K> ,KN<K>* > >
      );*/

      TheOperators->Add("*=",
                        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,K > >  );
     TheOperators->Add(".*=",
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Divc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_mul<KN_<K> ,Add_Mulc_KN_<K>* > > ,
       //- new OneBinaryOperator<set_eq_arrayp_mul<KN_<K> ,KN<K>* > >
	new OneBinaryOperator<set_eq_array_mul<KN_<K> ,KN_<K> > >
      );
// FH correction  01 nov 2005 FH  copy paste mistake eq_ exchange  ok  v2.0-3
/*    if(0)
     TheOperators->Add("/=",
        new OneBinaryOperator<set_eqarray_div<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_div<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,KN_<K> > >
     );*/

     TheOperators->Add("/=",
                       new OneBinaryOperator<set_eq_array_div<KN_<K> ,K > > );
    TheOperators->Add("./=",
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,Divc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_div<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,KN_<K> > >
     );
// end correction
     TheOperators->Add("+",
       new OneBinaryOperator<Op2_add0<Add_KN_<K>,KN_<K>,KN_<K> > >,
     //-  new OneBinaryOperator<Op2_add0<Add_KN_<K>,KN_<K>,KN_<K> > >(knrp,knrp),
       new OneBinaryOperator<Op2_add__n<Add_Mulc_KN_<K>,Mulc_KN_<K>,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_add__n<Add_Mulc_KN_<K>,KN_<K>,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_add__n<Add_Mulc_KN_<K>,Mulc_KN_<K> ,KN_<K> > >
       );

     TheOperators->Add("-",
       new OneBinaryOperator<Op2_sub0<Sub_KN_<K>,KN_<K> ,KN_<K> > >,
     //-  new OneBinaryOperator<Op2_sub0<Sub_KN_<K>,KN_<K> ,KN_<K> > >(knrp,knrp),
       new OneUnaryOperator<Op1_sub<Mulc_KN_<K>,KN_<K> > >,
       new OneBinaryOperator<Op2_sub__n<Add_Mulc_KN_<K>,Mulc_KN_<K>,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_sub__n<Add_Mulc_KN_<K>,KN_<K>,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_sub__n<Add_Mulc_KN_<K>,Mulc_KN_<K> ,KN_<K> > >
       );

     TheOperators->Add("*",
     //-  new OneBinaryOperator<Op2_mulpc<Mulc_KN_<K>,KN<K>*,K> >,
     //-  new OneBinaryOperator<Op2_mulcp<Mulc_KN_<K>,K,KN<K>*> >,
       new OneBinaryOperator<Op2_mulc<Mulc_KN_<K>,KN_<K>,K> >,
       new OneBinaryOperator<Op2_mulc<Mulc_KN_<K>,K,KN_<K> > >,
       new OneBinaryOperator<Op2_mulpcp<Mul_KNM_KN_<K>,KNM<K>*,KN<K>*> >,// A*b zzzzzzz
       
      // new OneBinaryOperator<Op2_mulp<Mul_KNM_KN_<K>,KNM_<K>,KN_<K>> >, // - add #1 mqi 2009
      // new OneBinaryOperator<Op2_dotproduct<K> >,
       new OneBinaryOperator<Op2_dotproduct_<K> >
     //-  ,new OneBinaryOperator<Op2_pbuild<outProduct_KN_<K>,KN<K>*,Transpose<KN_<K> > > >
       // ,new OneBinaryOperatorBug<Transpose<KN_<K> >,KNM<K>*  >
        ,new OneBinaryOperatorBug<Transpose<KN_<K> >,KNM_<K> >
        ,new OneBinaryOperatorBug<Transpose<KN_<K> >,Transpose<KNM<K>* >  >

       ,new OneBinaryOperator<Op2_pbuild<outProduct_KN_<K>,KN_<K>,Transpose<KN_<K> > > >
       ,new OneBinaryOperator<Op2_pbuild<outProduct_KN_<K>,Mulc_KN_<K>,Transpose<KN_<K> > > >

       );
    TheOperators->Add("*", new OneBinaryOperator<Op2_2p_<Mul_KNMh_KN_<K>, Transpose<KNM<K>*>, KN<K>*> >); // A'*b
    TheOperators->Add("=", new OneBinaryOperator<init_eqarray<KN<K>, Mul_KNMh_KN_<K> > >);
    TheOperators->Add("<-", new OneBinaryOperator<init_eqarray<KN<K>, Mul_KNMh_KN_<K> > >);

    TheOperators->Add("/",
                      new OneBinaryOperator<Op2_divc<Divc_KN_<K>,K,KN_<K> > >,
                      new OneBinaryOperator<Op2_divc<Mulc_KN_<K>,KN_<K>,K > >

                      );

//  nouvel operateur
     TheOperators->Add("+=",
        new OneBinaryOperator<set_eqarraypd_add<KNM<K> ,outProduct_KN_<K>* > >
       );

     TheOperators->Add("-=",
        new OneBinaryOperator<set_eqarraypd_sub<KNM<K> ,outProduct_KN_<K>* > >
       );

     TheOperators->Add("=",
        new OneBinaryOperator<set_eqarraypd<KNM<K> ,outProduct_KN_<K>* > >
       );
//   tested ok ...  FH
     TheOperators->Add("?:",
       new OneTernaryOperator3<Op3_p<if_arth_KN_<K>, KN_<K> > > ,
       new OneTernaryOperator3<Op3_paac<K > > ,
       new OneTernaryOperator3<Op3_pacc<K > > ,
       new OneTernaryOperator3<Op3_paca<K > >

       );
// end ...

// add mars 2006
// atype<KN_<K> >()->Add("(","",new OneOperator2_< pair<KN_<K>,KN_<long> > * ,KN_<K>  , KN_<long>  >(pBuild< KN_<K>   , KN_<long>  >));
 atype<KN_<K> >()->Add("(","",new OneOperator2_< pair<KN_<K>,KN_<long> > * ,KN_<K>  , KN_<long>  >(pBuild< KN_<K>   , KN_<long>  >,atype<KN_<K>  >(), atype<KN_<long> >() ));
 atype<KN<K> *>()->Add("(","",new OneOperator2_< pair<KN_<K>,KN_<long> > * ,KN_<K>  , KN_<long>  >(pBuild< KN_<K>   , KN_<long>  >,atype<KN<K> * >(), atype<KN_<long> >() ));
 //atype<KN_<K> >()->Add("(","",new OneOperator2_< pair<KN_<K>,KN_<long> > * ,KN_<K>  , KN_<long>  >(pBuild< KN_<K>   , KN_<long>  >,atype<KN_<K>  >(), knlp ));
 //atype<KN<K> *>()->Add("(","",new OneOperator2_< pair<KN_<K>,KN_<long> > * ,KN_<K>  , KN_<long>  >(pBuild< KN_<K>   , KN_<long>  >,atype<KN<K> * >(), knlp ));
 Add<KNM<K> *>("diag",".",new OneOperator1<KN_<K> ,KNM<K> *>(getdiag<K>) );
  Add<KNM_<K> >("diag",".",new OneOperator1<KN_<K> ,KNM_<K> >(getdiag_<K>) );
 Add<KNM<K> *>("asarray",".",new OneOperator1<KN_<K> ,KNM<K> *>(asarray<K>) );

 TheOperators->Add("<-",
                   new OneBinaryOperator<init_A_BI< K,Z,affectation<K>  > > );
 TheOperators->Add("=",
        new OneBinaryOperator<set_A_BI< K,Z,affectation<K>  > > ,
        new OneBinaryOperator<set_AI_B< K,Z,affectation<K>  > >
 );
 TheOperators->Add("+=",
        new OneBinaryOperator<set_A_BI< K,Z,affectation_add<K>  > > ,
        new OneBinaryOperator<set_AI_B< K,Z,affectation_add<K>  > >
 );
 TheOperators->Add("-=",
        new OneBinaryOperator<set_A_BI< K,Z,affectation_sub<K>  > > ,
        new OneBinaryOperator<set_AI_B< K,Z,affectation_sub<K>  > >
 );
// fin
  TheOperators->Add("\'",
      // new OneOperator1<Transpose<KN_<K> >,KN<K> *>(&Build<Transpose<KN_<K> >,KN<K> *>),
       new OneOperator1<Transpose<KN_<K> >,KN_<K> >(&Build<Transpose<KN_<K> >,KN_<K> >),
       new OneOperator1<Transpose<KNM<K> * >, KNM<K> * >(&Build<Transpose<KNM<K> * >,KNM<K> * >,10)  ,
       new OneOperator1<KNM_<K> , KNM_<K> >(Transp<K>, 1)
  );

     TheOperators->Add(".*",
       new OneBinaryOperator<Op2_build<DotStar_KN_<K>,KN_<K>,KN_<K> > > //-,
     //-  new OneBinaryOperator<Op2_build<DotStar_KN_<K>,KN_<K>,KN_<K> > >(knrp,knrp),
     //-  new OneBinaryOperator<Op2_build<DotStar_KN_<K>,KN_<K>,KN_<K> > >(knr_,knrp),
      //- new OneBinaryOperator<Op2_build<DotStar_KN_<K>,KN_<K>,KN_<K> > >(knrp,knr_)

      );


     TheOperators->Add("./",
       new OneBinaryOperator<Op2_build<DotSlash_KN_<K>,KN_<K>,KN_<K> > > //-,
     //- new OneBinaryOperator<Op2_build<DotSlash_KN_<K>,KN_<K>,KN_<K> > >(knrp,knrp),
     //-  new OneBinaryOperator<Op2_build<DotSlash_KN_<K>,KN_<K>,KN_<K> > >(knr_,knrp),
     //-  new OneBinaryOperator<Op2_build<DotSlash_KN_<K>,KN_<K>,KN_<K> > >(knrp,knr_)
      );

     TheOperators->Add("<<",
    //   new OneBinaryOperator<PrintPnd<KN<K>*> >,
       new OneBinaryOperator<Print<KNM_<K> > >,
       new OneBinaryOperator<Print<KN_<K> > >
       );
    TheOperators->Add("<<",
                      new OneBinaryOperator< PrintPnd< KN< KNM<K> >* > >,
                      new OneBinaryOperator< PrintPnd< KN< KN<K> >* > >
                      );
  // remove f.H ne marche pas
  //  TheOperators->Add("<<",
  //     new OneBinaryOperator<Op_WriteKN<K> >,
  //     new OneBinaryOperator<Op_WriteKNM<K> >
  //   );

     TheOperators->Add(">>",
        new OneBinaryOperator<Op_ReadKN<K> >,
        new OneBinaryOperator<Op_ReadKNM<K> >
      );
     atype<MyMap<String,K>*>()->Add("[","",new OneOperator2_<K*,MyMap<String,K>*,string*>(get_element<K>));
     TheOperators->Add("&",new OneOperator2_<bool,MyMap<String,K>*,string*>(exist_element<K>));
    TheOperators->Add("<<",new OneBinaryOperator<PrintP<MyMap<String,K>*> >);
    Add<MyMap<String,K>* >("n",".",new OneOperator1<Z,MyMap<String,K> *>(get_MyMap_n));
    
    // Add Mai 2009
    Dcl_Type<SetArray<K> >();
    TheOperators->Add("::",

		      new OneBinaryOperator<SetArray2<K> >,
		      new OneTernaryOperator3<SetArray3<K> >);
    TheOperators->Add("<-",
		      new OneOperator2_<KN<K> *,KN<K> *,SetArray<K> >(&set_init_array));

    TheOperators->Add("=",
		      new OneOperator2_<KN<K> *,KN<K> *,SetArray<K> >(&set_array),
		      new OneOperator2_<KN<K> *,KN<K> *,KN<K> * >(&set_arrayp),  //  to reomve ambiguity aug 2009
		      new OneOperator2_<KN_<K> ,KN_<K> ,SetArray<K> >(-1,&set_array_) // missing aug 2009 a(:)=1:3 less prioritaire
    );

    atype<MyMap<String,K>*>()->SetTypeLoop(atype<K*>(),atype<string**>());


    TheOperators->Add("{}",new ForAllLoop<E_ForAllLoopMapSI<K> >);

    TheOperators->Add("{}",new ForAllLoop<E_ForAllLoopRNM<K> >);
    TheOperators->Add("{}",new ForAllLoop<E_ForAllLoopRN<K> >);

    TheOperators->Add("<-",new InitMapfromArray<MyMap<String,K>*,string *,K,true> );

}

template<class R,class A,class B=A,class BB=B>
class  OneOperator1F_KN_ : public OneOperator {
    aType r; //  return type
    typedef  A (*func)( B ) ;
    func  f;
    public:
    E_F0 * code(const basicAC_F0 & args) const
     { return  new Op(f,t[0]->CastTo(args[0]));}
    OneOperator1F_KN_(func  ff):
      OneOperator(map_type[typeid(R).name()],map_type[typeid(BB).name()]),f(ff){}

 class Op :public  E_F0 { public:
  typedef  A (*func)(B ) ;
  func f;
  Expression a;
  Op(func ff,Expression aa) : f(ff),a(aa) {}
  AnyType operator()(Stack s)  const  {return SetAny<R>( R(f, GetAny<BB>( (*a)(s)) ) );}
   bool EvaluableWithOutStack() const
      {return a->EvaluableWithOutStack() ;} //
   bool MeshIndependent() const
      {return a->MeshIndependent();} //

};

};
template<class K,class KK>
void ArrayOperatorF()
{
     Dcl_Type<F_KN_<K,K,K,KK> >();


     Global.Add("exp","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(exp));
     Global.Add("log","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(log));
     Global.Add("log10","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(log10));
     Global.Add("sqrt","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(sqrt));
     Global.Add("sin","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(sin));
     Global.Add("cos","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(cos));
     Global.Add("tan","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(tan));
     Global.Add("cosh","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(cosh));
     Global.Add("sinh","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(sinh));
     Global.Add("tanh","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(tanh));
    // Global.Add("acos","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(acos));
    // Global.Add("asin","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(asin));
    // Global.Add("atan","(",new OneOperator1F_KN_<F_KN_<K,K,K,KK>,K,KK,KN_<K> >(atan));

     TheOperators->Add("=",new OneBinaryOperator<set_eq_array<KN_<K> ,F_KN_<K,K,K,KK> > > ); // add FH juin 2005
     TheOperators->Add("+=",new OneBinaryOperator<set_eq_array_add<KN_<K> ,F_KN_<K,K,K,KK> > > ); // add FH juin 2005
     TheOperators->Add("-=",new OneBinaryOperator<set_eq_array_sub<KN_<K> ,F_KN_<K,K,K,KK> > > ); // add FH juin 2005
     TheOperators->Add("/=",new OneBinaryOperator<set_eq_array_div<KN_<K> ,F_KN_<K,K,K,KK> > > ); // add FH juin 2005
     TheOperators->Add("*=",new OneBinaryOperator<set_eq_array_mul<KN_<K> ,F_KN_<K,K,K,KK> > > ); // add FH juin 2005

    TheOperators->Add("<-", new OneOperator2_<KN<K> *,KN<K> *,F_KN_<K,K,K,KK> >(set_init_N));

}
// Add nov 2019  version 4.4-3 FH
template<class K> struct KN_rmeps {KN_<K> v;
    KN_rmeps(KN_<K> vv):v(vv) {}
} ;
template<class K> KN_rmeps<K> build_rmeps(KN_<K> v){ return KN_rmeps<K>(v);}
#endif
