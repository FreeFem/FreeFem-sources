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

#include "config-wrapper.h"

#include <complex>
#include "AFunction.hpp"
#include <cstdarg>
#include <cstring>
#include "error.hpp"
#include "lex.hpp"

#include "RNM.hpp"

#include "Operator.hpp"
// for exec routine 
#include "rgraph.hpp"
#include "InitFunct.hpp"
#include <queue>
#include "array_resize.hpp"

template <class T>
struct affectation: binary_function<T, T, T>
{
	T& operator()(T& x, const T& y) const {return (x=y);}
};

template <class T>
struct affectation_add: binary_function<T, T, T>
{
	T& operator()(T& x, const T& y) const {return (x=+y);}
};

template <class T>
struct affectation_sub: binary_function<T, T, T>
{
	T& operator()(T& x, const T& y) const {return (x=-y);}
};



extern Map_type_of_map map_type_of_map ; //  to store te type 
extern Map_type_of_map map_pair_of_type ; //  to store te type 

extern basicForEachType *  typevarreal,  * typevarcomplex;  //  type of real and complex variable

extern int TheCurrentLine; // unset: by default
extern long mpisize,mpirank;

template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
template<class T> inline T Min (const T &a,const T & b){return a < b ? a : b;}
template<class T> inline T Abs (const T &a){return a <0 ? -a : a;}
template<class T> inline T Max (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}
template<class T> inline T Min (const T &a,const T & b,const T & c){return Min(Min(a,b),c);}
template<class T> inline T Square (const T &a){return a*a;}


 
template<class K> 
struct Op2_dotproduct: public binary_function<Transpose<KN_<K> >,KN<K> *,K> { 
  static K f( Transpose<KN_<K> > const & a, KN<K> * const& b)  
   { return (conj(a.t),*b);} }; 

template<class K> 
struct Op2_dotproduct_: public binary_function<Transpose<KN_<K> >,KN_<K> ,K> { 
  static K f( Transpose<KN_<K> > const & a, KN_<K>  const& b)  
   { return (conj(a.t),b);} }; 
   
template<class A,class B>  A Build(B b) {  return A(b);}

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

template<class R,class A> A  SortKn(const A  & ca){ 
    A a(ca);
    HeapSort<R>(&a[0],a.n,a.step);
    return a;}

template<class R> KN<R> *  SortpKn( KN<R> * const & pa){ 
    KN<R> &a(*pa);
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
    
template<>
inline   string ** get_element<string*>( MyMap<String,string*> *  const  &  a,string*  const   & b)
 { string** ret=  &((*a)[*b]); // correction FH feb 2004
    if( *ret ==0) *ret = new string(""); //  string vide ???
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
     
template<class RR,class A,class B,class C>  
RR get_element_lineorcol(const A &  a,const B & b,const C & c){ 
 //  cout << b << " .... " << ((*a)(SubArray(1,b),c)) << endl;;
    return  ((*a)(b,c));}

    

template<class RR,bool isinit>
class  InitArrayfromArray : public OneOperator { 
public:
    typedef KN<RR> * A;
    typedef KN<RR> * R;
    typedef E_Array B;
    
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
	if(atype<RR>()->CastingFrom(tt[i].right() ) ) 
	  {
          tab[i]=atype<RR>()->CastTo(tt[i]);
	    what[i]=0;
	  }
	else if(atype<KN_<RR> >()->CastingFrom(tt[i].right() ) ) 
	  {
	    tab[i]=atype<KN_<RR> >()->CastTo(tt[i].RightExp());
	    what[i]=1;
	  }      
	else 
	  CompileError(" we are waiting for scalar or vector of scalar");
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
      if (isinit) 
        a->init(n);
      else
	a->resize(n);
      
      for (int i=0,j=0 ;i<N; j += nn[i++])
	
        if (what[i]==0)
          (*a)[j]= GetAny<RR>(v[i]);
        else if (what[i]==1) 
          (*a)(SubArray(nn[i],j)) = GetAny<KN_<RR> >(v[i]);
      return SetAny<R>(a);
    } 
    bool MeshIndependent() const     {return  mi;} // 
    ~CODE() { delete [] tab; delete[] what;}
    operator aType () const { return atype<R>();}    
  }; // end sub class CODE
  
  
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new CODE(t[0]->CastTo(args[0]),*dynamic_cast<const E_Array*>( t[1]->CastTo(args[1]).LeftValue()));} 
    InitArrayfromArray():   OneOperator(atype<R>(),atype<A>(),atype<B>())  {}
  
};

template<class RR,bool isinit>
class  InitMatfromAArray : public OneOperator { 
public:
    typedef KNM<RR> * A;
    typedef KNM<RR> * R;
    typedef E_Array B;
    
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
      if (isinit) 
        a->init(N,M);
      else
	a->resize(N,M);
      
       for (int i =0;i<N;++i)
       for (int j =0;j<M;++j)
          (*a)(i,j)=   GetAny< RR >( (*(tab[i][j]))(stack)) ; 
      return SetAny<R>(a);
    } 
    bool MeshIndependent() const     {return  mi;} // 
    ~CODE() { for (int i=0;i<N;i++) delete [] tab[i]; delete [] tab; }
    operator aType () const { return atype<R>();}    
  }; // end sub class CODE
  
  
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new CODE(t[0]->CastTo(args[0]),*dynamic_cast<const E_Array*>( t[1]->CastTo(args[1]).LeftValue()));} 
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
	  CompileError(" we are waiting for scalar or vector of scalar");
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
   


template<class K> long get_n(KN<K> * p){ return p->N();}
template<class K> long get_n(KNM<K> * p){ return p->N();}

template<class K> long get_m(KNM<K> * p){ return p->M();}
template<class K> K get_max(KN<K> * p){ return p->max();}
template<class K> K get_min(KN<K> * p){ return p->min();}
template<class K> K get_sum(KN<K> * p){ return p->sum();}
template<class K> double get_l2(KN<K> * p){ return p->l2();}
template<class K> double get_l1(KN<K> * p){ return p->l1();}
template<class K> double get_linfty(KN<K> * p){ return p->linfty();}

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
 
template<class K>
void ArrayDCL()
{
    Dcl_TypeandPtr<KN<K> >(0,0,0,::Destroy<KN<K> >);
  //  Dcl_Type<KN<Complex> *>(0,::Destroy<KN<Complex> >);
   // Dcl_Type<KN<K> *>(0,::Destroy<KN<K> >); // Modif 17102005 
   // attention un exp KN<> * right est un KN<> et non un KN<> *

    Dcl_Type<KNM<K> *>(0,::Destroy<KNM<K> >);
   //  Dcl_Type< Transpose<KN<K> *> > ();  remove mars 2006 FH 
    Dcl_Type< outProduct_KN_<K>* >();
    Dcl_Type< Transpose<KN_<K> > > ();
    Dcl_Type< Transpose< KNM<K> *> >();
    //Dcl_Type< Transpose<KN<Complex> > > ();
    Dcl_TypeandPtr<KN_<K> >(0,0,0,0);
    //Dcl_TypeandPtr<KN_<Complex> >(0,0,0,0);

    Dcl_Type<Add_KN_<K> >();
    
    Dcl_Type<DotStar_KN_<K> >();
    Dcl_Type<DotSlash_KN_<K> >();
    Dcl_Type<Sub_KN_<K> >();
    Dcl_Type<Mulc_KN_<K> >();
    Dcl_Type<Mul_KNM_KN_<K> >();
    Dcl_Type<Add_Mulc_KN_<K> *>();
    Dcl_Type<if_arth_KN_<K> *>();
    // for    B(I) and B(I^-1)
    Dcl_Type<pair<KN_<K>,Inv_KN_long> *>();
    Dcl_Type<pair<KN_<K>,KN_<long> > *>();
    

     map_type[typeid(KN_<K> ).name()]->AddCast(
       new E_F1_funcT<KN_<K>,KN_<K>*>(UnRef<KN_<K> >),
     //  new E_F1_funcT<KN_<K>,KN<K>*>(UnRef<KN_<K>,KN<K>* >), inutil cas KN<K> est right expression de KN<K>* 
       new E_F1_funcT<KN_<K>,KN<K> >(Cast<KN_<K>,KN<K> >)
       
       );
    //   ,new E_F1_funcT<KN_<K>,K>(ValueToKN_<K>),
    //   new E_F1_funcT<KN_<K>,K*>(PtrToKN_<K>)       
       
     // Ajoute FH   
     map_type[typeid(KN<K> ).name()]->AddCast(
       new E_F1_funcT<KN<K>,KN<K>*>(UnRef<KN<K> >)
    //   ,new E_F1_funcT<KN_<K>,K>(ValueToKN_<K>),
    //   new E_F1_funcT<KN_<K>,K*>(PtrToKN_<K>)       
       ); 
    map_type_of_map[make_pair(atype<long>(),atype<K>())]=atype<KN<K>*>(); // vector
    map_pair_of_type[make_pair(atype<long>(),atype<long>())] =atype<pair<long,long> >();   
    map_type_of_map[make_pair(atype<pair<long,long> >(),atype<K>())]=atype<KNM<K>*>(); // matrix                                               
}



template<class A,class B> pair<A,B> * pBuild(const A & a,const B & b)
  { return new pair<A,B>(a,b);}

// add mars 2006
template<class K,class L,class OP>
struct set_A_BI: public binary_function<KN_<K>,pair<KN_<K>, KN_<L> > *,KN_<K> > {
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

template<class K,class L,class OP>
struct set_AI_B: public binary_function<pair<KN_<K>, KN_<L> > * ,KN_<K>, NothingType > {
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

  
extern aType aaaa_knlp;
template<class K,class Z>
void ArrayOperator()
{
     Dcl_Type< Resize<KN<K> > > ();
     Dcl_Type< Resize<KNM<K> > > ();
     aType knrp = atype<KN<K> *>();
     aType knr_ = atype<KN_<K> >();
     typedef KN<Z> ZN;
      
     aType knlp=  aaaa_knlp ;
     //atype<KN<Z> *>();
     // aType knl_ = atype<KN_<Z> >();
    
     atype<KN<K>* >()->Add("[","",new OneOperator2_<K*,KN<K>*,Z >(get_elementp_<K,KN<K>*,Z>));
     atype<KN<K>* >()->Add("(","",new OneOperator2_<K*,KN<K>*,Z >(get_elementp_<K,KN<K>*,Z>));
     atype<KN_<K> >()->Add("(","",new OneOperator2_<KN_<K>,KN_<K>,char >(fSubArrayc<KN_<K> >));
     atype<KN_<K> >()->Add("(","",new OneOperator2_<KN_<K>,KN_<K>,SubArray>(fSubArray<K> ));
     atype<KN<K>*>()->Add("(","",new OneOperator2_<KN_<K>,KN<K>*,SubArray>(fSubArrayp<K> ));
     atype<KN<K>* >()->Add("(","",new OneOperator2_<KN<K>*,KN<K>*,char >(fSubArrayc<KN<K>* >));

     atype<KNM<K>* >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,Z,SubArray >(get_element_is<KN_<K>,KNM<K>*,Z,SubArray>));
     atype<KNM<K>* >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,SubArray,Z >(get_element_si<KN_<K>,KNM<K>*,SubArray,Z>));
     atype<KNM<K>* >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,Z,char >(get_element_lineorcol<KN_<K>,KNM<K>*,Z,char>));
     atype<KNM<K>* >()->Add("(","",new OneOperator3_<KN_<K>,KNM<K>*,char,Z >(get_element_lineorcol<KN_<K>,KNM<K>*,char,Z>));

     atype<KNM<K>* >()->Add("(","",new OneOperator3_<K*,KNM<K>*,Z,Z >(get_elementp2_<K,KNM<K>*,Z,Z>));

     Add<KN<K> *>("sum",".",new OneOperator1<K,KN<K> *>(get_sum));
     Add<KN<K> *>("min",".",new OneOperator1<K,KN<K> *>(get_min));
     Add<KN<K> *>("max",".",new OneOperator1<K,KN<K> *>(get_max));
     Add<KN<K> *>("l2",".",new OneOperator1<double,KN<K> *>(get_l2));
     Add<KN<K> *>("l1",".",new OneOperator1<double,KN<K> *>(get_l1));
     Add<KN<K> *>("linfty",".",new OneOperator1<double,KN<K> *>(get_linfty));
     
     
     Add<KN_<K> >("sum",".",new OneOperator1_<K,KN_<K> >(get_sum0<K,KN_<K> >));
     Add<KN_<K> >("min",".",new OneOperator1_<K,KN_<K> >(get_min0<K,KN_<K> >));
     Add<KN_<K> >("max",".",new OneOperator1_<K,KN_<K> >(get_max0<K,KN_<K> >));
     Add<KN_<K> >("l2",".",new OneOperator1_<double,KN_<K> >(get_l2_0<double,KN_<K> >));
     Add<KN_<K> >("l1",".",new OneOperator1_<double,KN_<K> >(get_l1_0<double,KN_<K> >));
     Add<KN_<K> >("linfty",".",new OneOperator1_<double,KN_<K> >(get_linfty_0<double,KN_<K> >));
    
     Add<KN<K> >("sum",".",   new OneOperator1_<K,KN<K> >(get_sum0<K,KN<K> >));
     Add<KN<K> >("min",".",   new OneOperator1_<K,KN<K> >(get_min0<K,KN<K> >));
     Add<KN<K> >("max",".",   new OneOperator1_<K,KN<K> >(get_max0<K,KN<K> >));
     Add<KN<K> >("l2",".",    new OneOperator1_<double,KN<K> >(get_l2_0<double,KN<K> >));
     Add<KN<K> >("l1",".",    new OneOperator1_<double,KN<K> >(get_l1_0<double,KN<K> >));
     Add<KN<K> >("linfty",".",new OneOperator1_<double,KN<K> >(get_linfty_0<double,KN<K> >));
     

     Add<KN<K> *>("resize",".",new OneOperator1< Resize<KN<K> >,KN<K> *>(to_Resize));
     Add<KNM<K> *>("resize",".",new OneOperator1< Resize<KNM<K> >,KNM<K> *>(to_Resize));
     
     Add<Resize<KN<K> > >("(","",new OneOperator2_<KN<K> *,Resize<KN<K> > , Z   >(resize1));
     Add<Resize<KNM<K> > >("(","",new OneOperator3_<KNM<K> *,Resize<KNM<K> > , Z, Z  >(resize2));

     TheOperators->Add("<-", 
       new OneOperator2_<KN<K> *,KN<K> *,Z>(&set_init),
       new InitArrayfromArray<K,true>,
       new OneOperator2_<KN<K> *,KN<K> *,KN<K> >(&set_init),
       new OneOperator2_<KN<K> *,KN<K> *,KN_<K> >(&set_init)		       
     //  new OneOperator2_<KN<K> *,KN<K> *,KN<K> * >(&set_initp)
       );
     TheOperators->Add("<-", 
        new OneOperator3_<KNM<K> *,KNM<K> *,Z,Z>(&set_init2),
        new InitMatfromAArray<K,true>
       );
       
     Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,Z>(&set_init));
   //  Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,KN<K> >(&set_init));
     //Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,KN_<K> >(&set_init));
    // Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,KN<K> * >(&set_initp));
     Add<KNM<K> *>("<-","(",new OneOperator3_<KNM<K> *,KNM<K> *,Z,Z>(&set_init2));
     Add<KN<K> *>("<-","(",new InitArrayfromArray<K,true>);
     Add<KNM<K> *>("<-","(",new InitMatfromAArray<K,true>);
     Add<KN<K> *>("n",".",new OneOperator1<Z,KN<K> *>(get_n));
     Add<KNM<K> *>("n",".",new OneOperator1<Z,KNM<K> *>(get_n));
     Add<KNM<K> *>("m",".",new OneOperator1<Z,KNM<K> *>(get_m));
     
//     AddOpeqarray<set_eqarray,KN,K>("=");

     TheOperators->Add("=", new InitArrayfromArray<K,false>
       );
     TheOperators->Add("=", new InitMatfromAArray<K,false>
       );
     TheOperators->Add("=", new SetArrayofKNfromKN<K>
       );
     
     TheOperators->Add("=",
        new OneBinaryOperator<set_eqarray<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,KN_<K> > > , // Add FH juin 2005         
        new OneBinaryOperator<set_eqarraypd<KN<K> ,Add_Mulc_KN_<K>* > > , // Add FH aug 2005     
        new OneBinaryOperator<set_eqarraypd<KN<K> ,if_arth_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarrayp<KN<K> ,KN<K>* > >       
      );
  // add august 2007 
     TheOperators->Add("<-",
		      // new OneBinaryOperator<set_eqarray<KN<K> ,K > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,Add_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,DotStar_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,DotSlash_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,Sub_KN_<K> > > ,
		       new OneBinaryOperator<init_eqarray<KN<K> ,Mulc_KN_<K> > > ,
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
        new OneBinaryOperator<set_eq_array<KN_<K> ,KN_<K> > > , // add FH juin 2005
        new OneBinaryOperator<set_eq_array<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array<KN_<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arraypd<KN_<K> ,Add_Mulc_KN_<K>* > > , // Add FH aug 2005     
        new OneBinaryOperator<set_eq_arrayp<KN_<K> ,KN<K>* > >       
      );

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
     TheOperators->Add("+=",
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,K > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_add<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arraypd_add<KN_<K> ,if_arth_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp_add<KN_<K> ,KN<K>* > >        
      );
      
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
      );
      
     TheOperators->Add("-=",
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,K > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_sub<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arraypd_sub<KN_<K> ,if_arth_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp_sub<KN_<K> ,KN<K>* > >        
      );
      
     TheOperators->Add("*=",
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,K > >  ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_mul<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarrayp_mul<KN<K> ,KN<K>* > >       
      );
 
      TheOperators->Add("*=",
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,K > >  ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_mul<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp_mul<KN_<K> ,KN<K>* > >       
      );
// FH correction  01 nov 2005 FH  copy paste mistake eq_ exchange  ok  v2.0-3 
     TheOperators->Add("/=",
        new OneBinaryOperator<set_eqarray_div<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_div<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarrayp_div<KN<K> ,KN<K>* > >        
     );

     TheOperators->Add("/=",
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,K > > ,
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN_<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_div<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp_div<KN_<K> ,KN<K>* > >        
     );
// end correction 
     TheOperators->Add("+",
       new OneBinaryOperator<Op2_add0<Add_KN_<K>,KN_<K>,KN_<K> > >,
       new OneBinaryOperator<Op2_add0<Add_KN_<K>,KN_<K>,KN_<K> > >(knrp,knrp),
       new OneBinaryOperator<Op2_add__n<Add_Mulc_KN_<K>,Mulc_KN_<K>,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_addp_n<Add_Mulc_KN_<K>,KN<K>*,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_add_pn<Add_Mulc_KN_<K>,Mulc_KN_<K> ,KN<K>* > >
       );

     TheOperators->Add("-",
       new OneBinaryOperator<Op2_sub0<Sub_KN_<K>,KN_<K> ,KN_<K> > >,
       new OneBinaryOperator<Op2_sub0<Sub_KN_<K>,KN_<K> ,KN_<K> > >(knrp,knrp),
       new OneUnaryOperator<Op1_subp<Mulc_KN_<K>,KN<K>*> >,
       new OneBinaryOperator<Op2_sub__n<Add_Mulc_KN_<K>,Mulc_KN_<K>,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_subp_n<Add_Mulc_KN_<K>,KN<K>*,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_sub_pn<Add_Mulc_KN_<K>,Mulc_KN_<K> ,KN<K>* > >
       );
     TheOperators->Add("*",
       new OneBinaryOperator<Op2_mulpc<Mulc_KN_<K>,KN<K>*,K> >,
       new OneBinaryOperator<Op2_mulcp<Mulc_KN_<K>,K,KN<K>*> >,
       new OneBinaryOperator<Op2_mulc<Mulc_KN_<K>,KN_<K>,K> >,
       new OneBinaryOperator<Op2_mulc<Mulc_KN_<K>,K,KN_<K> > >,
       new OneBinaryOperator<Op2_mulpcp<Mul_KNM_KN_<K>,KNM<K>*,KN<K>*> >,
       new OneBinaryOperator<Op2_dotproduct<K> >,
       new OneBinaryOperator<Op2_dotproduct_<K> > 
       ,new OneBinaryOperator<Op2_pbuild<outProduct_KN_<K>,KN<K>*,Transpose<KN_<K> > > >
       ,new OneBinaryOperator<Op2_pbuild<outProduct_KN_<K>,KN_<K>,Transpose<KN_<K> > > > 
       ,new OneBinaryOperator<Op2_pbuild<outProduct_KN_<K>,Mulc_KN_<K>,Transpose<KN_<K> > > > 
             
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
// not tested
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
 atype<KN_<K> >()->Add("(","",new OneOperator2_< pair<KN_<K>,KN_<long> > * ,KN_<K>  , KN_<long>  >(pBuild< KN_<K>   , KN_<long>  >,atype<KN_<K>  >(), knlp ));
 atype<KN<K> *>()->Add("(","",new OneOperator2_< pair<KN_<K>,KN_<long> > * ,KN_<K>  , KN_<long>  >(pBuild< KN_<K>   , KN_<long>  >,atype<KN<K> * >(), knlp ));
 
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
       new OneOperator1<Transpose<KN_<K> >,KN<K> *>(&Build<Transpose<KN_<K> >,KN<K> *>),
       new OneOperator1<Transpose<KN_<K> >,KN_<K> >(&Build<Transpose<KN_<K> >,KN_<K> >),
       new OneOperator1<Transpose<KNM<K> * >, KNM<K> * >(&Build<Transpose<KNM<K> * >,KNM<K> * >)            
  );
       
     TheOperators->Add(".*",
       new OneBinaryOperator<Op2_build<DotStar_KN_<K>,KN_<K>,KN_<K> > >,
       new OneBinaryOperator<Op2_build<DotStar_KN_<K>,KN_<K>,KN_<K> > >(knrp,knrp),
       new OneBinaryOperator<Op2_build<DotStar_KN_<K>,KN_<K>,KN_<K> > >(knr_,knrp),
       new OneBinaryOperator<Op2_build<DotStar_KN_<K>,KN_<K>,KN_<K> > >(knrp,knr_)
       
      );

      
     TheOperators->Add("./",
       new OneBinaryOperator<Op2_build<DotSlash_KN_<K>,KN_<K>,KN_<K> > >,
      new OneBinaryOperator<Op2_build<DotSlash_KN_<K>,KN_<K>,KN_<K> > >(knrp,knrp),
       new OneBinaryOperator<Op2_build<DotSlash_KN_<K>,KN_<K>,KN_<K> > >(knr_,knrp),
       new OneBinaryOperator<Op2_build<DotSlash_KN_<K>,KN_<K>,KN_<K> > >(knrp,knr_)
      );
      
     TheOperators->Add("<<",
       new OneBinaryOperator<PrintPnd<KN<K>*> >,
       new OneBinaryOperator<PrintPnd<KNM<K>*> >,
       new OneBinaryOperator<Print<KN_<K> > >
       ); 
     
       
     TheOperators->Add(">>",
        new OneBinaryOperator<Op_ReadKN<K> >
      );            

     map_type[typeid(MyMap<String,K>*).name()] = new ForEachType<MyMap<String,K>*>(Initialize<MyMap<String,K> >,Delete<MyMap<String,K> >) ;
         
     map_type_of_map[make_pair(atype<string*>(),atype<K>())]=atype<MyMap<String,K>*>(); 
     
     atype<MyMap<String,K>*>()->Add("[","",new OneOperator2_<K*,MyMap<String,K>*,string*>(get_element<K>));

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
     Dcl_Type<F_KN_<K,K,KK> >();


     Global.Add("exp","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(exp));
     Global.Add("log","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(log));
     Global.Add("log10","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(log10));
     Global.Add("sqrt","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(sqrt));
     Global.Add("sin","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(sin));
     Global.Add("cos","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(cos));
     Global.Add("tan","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(tan));
     Global.Add("cosh","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(cosh));
     Global.Add("sinh","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(sinh));
     Global.Add("tanh","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(tanh));
    // Global.Add("acos","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(acos));
    // Global.Add("asin","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(asin));
    // Global.Add("atan","(",new OneOperator1F_KN_<F_KN_<K,K,KK>,K,KK,KN_<K> >(atan));

     TheOperators->Add("=",new OneBinaryOperator<set_eq_array<KN_<K> ,F_KN_<K,K,KK> > > ); // add FH juin 2005
     TheOperators->Add("+=",new OneBinaryOperator<set_eq_array_add<KN_<K> ,F_KN_<K,K,KK> > > ); // add FH juin 2005
     TheOperators->Add("-=",new OneBinaryOperator<set_eq_array_sub<KN_<K> ,F_KN_<K,K,KK> > > ); // add FH juin 2005
     TheOperators->Add("/=",new OneBinaryOperator<set_eq_array_div<KN_<K> ,F_KN_<K,K,KK> > > ); // add FH juin 2005
     TheOperators->Add("*=",new OneBinaryOperator<set_eq_array_mul<KN_<K> ,F_KN_<K,K,KK> > > ); // add FH juin 2005
  
}
