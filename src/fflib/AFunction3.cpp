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
/*
struct SubArray2: public binary_function<long,long,SubArray> { 
  static SubArray f(const long & a,const long & b)  { 
   // cout << "SubArray: " << a << " " << b << endl;
    return SubArray(b-a+1,a);} }; 
struct SubArray3: public ternary_function<long,long,long,SubArray> { 
  static SubArray f(const long & a,const long & b,const long & c)  {  
  // cout << "SubArray: " << a << " " << b << " " <<  c << endl;
   return SubArray((b-a+1)/c,a,c);} }; 
*/

 
template<class K> 
struct Op2_dotproduct: public binary_function<Transpose<KN<K> *>,KN<K> *,K> { 
  static K f( Transpose<KN<K> *> const & a, KN<K> * const& b)  
   { return (conj(*a.t),*b);} }; 

template<class K> 
struct Op2_dotproduct_: public binary_function<Transpose<KN_<K> >,KN_<K> ,K> { 
  static K f( Transpose<KN_<K> > const & a, KN_<K>  const& b)  
   { return (conj(a.t),b);} }; 
   
template<class A,class B>  A Build(B b) {  return A(b);}
  
  

/*
long Exit(long i) {throw(ErrorExec("Exit",i));return 0;}
bool Assert(bool b) {if (!b) throw(ErrorExec("exec assert",1));return true;}
*/
  
inline void MyAssert(int i,char * ex,char * file,long line)
{if (i) {
    cout << "CompileError assertion :  " << ex << " in file " << file << "  line = " << line << endl; 
     CompileError();}
 }

/*  
  C_F0 *pOne=0,*pZero=0,*pminusOne=0;
// const C_F0 & One(*pOne), &Zero(*pZero);
 
 Polymorphic * TheOperators=0, //=new Polymorphic(), 
             * TheRightOperators=0; //=new Polymorphic();

TableOfIdentifier Global;

 long E_Border::Count =0;

typedef list<TableOfIdentifier *> ListOfTOfId;

  ListOfTOfId tables_of_identifier;

  const int AC_F0::MaxSize=1024; // maximal number of parameters
*/  
/*  
template<class R>
class  OneOperator0 : public OneOperator {
 class E_F0_F :public  E_F0 { public:
  typedef  R (*func)( ) ; 
  func f;
  E_F0_F(func ff)  : f(ff) {}
  AnyType operator()(Stack )  const  {return SetAny<R>( f()) ;}  
     operator aType () const { return atype<R>();} 

};

  //  aType r; //  return type
    typedef  R (*func)() ; 
    func  f;
    public: 
    E_F0 * code(const basicAC_F0 & ) const 
     { return  new E_F0_F(f);} 
    OneOperator0(func  ff): OneOperator(map_type[typeid(R).name()]),f(ff){}
};
*/
/*
template<class R>
class  OneOperatorConst : public OneOperator {
    E_F0 * e;
    public: 
    E_F0 * code(const basicAC_F0 & ) const  { return  e;} 
    OneOperatorConst(E_F0 * ee):  OneOperator(map_type[typeid(R).name()]),e(ee){}
};

class  OneOperator_array : public OneOperator {public:
    E_F0 * code(const basicAC_F0 & a) const 
     { return  new E_Array(a);} 
    OneOperator_array(): OneOperator(atype<E_Array>(),true) {}
};
class  OneOperator_border : public OneOperator {public:
    E_F0 * code(const basicAC_F0 & a) const 
     { if (a.size()==1 && a[0].left()==atype<E_Array>() ) 
        return new E_Border(dynamic_cast<const E_Array*>(a[0].LeftValue()));
        else     
        return  new E_Border(a);} 
    OneOperator_border(): OneOperator(atype<const E_Border *>(),true) {}
};

class  OneOperator_border_label : public OneOperator {public:
  class Op : public E_F0 {public:
   const  E_Border *b;
      Op( const  E_Border *bb) : b(bb) {}
      AnyType operator()(Stack)  const { return SetAny<long>(b->label);}
   };
    E_F0 * code(const basicAC_F0 & a) const 
     {  const  E_Border * b = dynamic_cast<const E_Border *>(a[0].LeftValue());
        return new Op(b);} 
    OneOperator_border_label(): OneOperator(atype<long>(),atype<const E_Border *>()) {}
};



template<class RR> RR LIncremantation(RR* a){ return ++(*a);}
template<class RR> RR RIncremantation(RR* a){ return (*a)++;}
template<class RR> RR LDecremantation(RR* a){ return --(*a);}
template<class RR> RR RDecremantation(RR* a){ return (*a)--;}

template<class RR,class B>
 RR * New_form_string(string * s) {B * r=  new B(s);delete *s;return r;}
 
 


typedef MyMap<String,double> mapSd ;
*/
template<class K>
inline   K * get_element( MyMap<String,K> *  const  &  a,string*  const   & b)
 { K * ret=  &((*a)[*b]); // correction FH feb 2004
  //  cout << "get_element " << *b << " : " << ret << " = "<< * ret << endl;
    delete b;
    return ret;}
    
template<>
inline   string ** get_element<string*>( MyMap<String,string*> *  const  &  a,string*  const   & b)
 { string** ret=  &((*a)[*b]); // correction FH feb 2004
    if( *ret ==0) *ret = new string(""); //  string vide ???
     // cout << "get_element " << *b << " : " << ret << " = "<< * ret << endl;
    delete b;
    return ret;}

inline   string ** get_elements( MyMap<String,String> *  const  &  a,string*  const   & b)
 { String* Sret=  &((*a)[*b]); // correction FH feb 2004
    delete b;
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
/* 
template<class RR> RR Abs(RR a) { return a<0?-a:a;}

template<class R,class A,class B>
R *MakePtrWithDel( A  const & a)
{ R *r= new B(a->c_str());
  delete a;
  return r;}

template<class R,class RR> 
struct Op1_new_pstring: public unary_function<string*,R> { 
  static R f(string * const & a)  {R r =  new RR(a->c_str()); delete a;return r;} }; 

template<class R,class RR> 
struct Op2_set_pstring: public binary_function<R,string*,R> { 
  static R  f(R const & p,string * const & a)  {*p =  new RR(a->c_str());
   if ( !*p || !**p) { cerr << " Error openning file " << *a << endl; ExecError("Error openning file");}
   delete a;return p;} }; 

template<class R,class RR> 
struct Op2_set_pstringiomode: public ternary_function<R,string*,ios::openmode,R> { 
  static R  f(R const & p,string * const & a,const ios::openmode & mode) 
   {*p =  new RR(a->c_str(),mode); delete a;return p;} }; 

AnyType FWhile(Stack s ,Expression test,Expression ins)
{ 
  AnyType a;
  while ( GetAny<bool>((*test)(s)))
     try  { a=(*ins)(s);}
     catch ( E_exception & e) { 
       if (e.code == E_exception::e_break) break;
       else if  (e.code == E_exception::e_continue) continue;
       }
  return a;
}
    
AnyType FFor(Stack s ,Expression i0,Expression i1,Expression i2,Expression ins)
{ 
  AnyType a;
  for ( (*i0)(s);GetAny<bool>((*i1)(s));(*i2)(s))
   {
     try  {a=(*ins)(s);}
     catch ( E_exception & e) { 
        if (verbosity>50)
        cerr << "FFor " << e.what() << e.code << endl; 
       if (e.code == E_exception::e_break) break;
       else if  (e.code == E_exception::e_continue) continue;
       }
   }
  return a;
}

AnyType FIf(Stack s ,Expression test,Expression i1,Expression i2,Expression )
 {  AnyType a;
   if (GetAny<bool>((*test)(s))) 
      a=(*i1)(s);
   	else if (i2) 
   	  a=(*i2)(s);   	
   	return a;
 }



aType TypeArray(aType b,aType a)
{ // type of  b[a]
   aType r=map_type_of_map[make_pair(a->right(),b->right())];
   if (!r) {
      cerr << "Sorry is not possible to make a map "<< *b->right() << " [" << *a->right() << "]" << endl;
      cerr << " list: " << endl;
      Map_type_of_map::const_iterator i;
      for(i=map_type_of_map.begin();i!=map_type_of_map.end();i++)
        cerr << "\t " << *i->first.second << " [" << *i->first.first << "]" << "=" << *i->second << endl;        
      CompileError();
   }
   return r;
}

aType TypeTemplate(aType b,aType a)
{ // type of  b[a]
   aType r=map_type_of_map[make_pair(b,a)];
   if (!r) {
      cerr << "Sorry is not possible to make a map "<< *b << "<" << *a << ">" << endl;
      cerr << " list: " << endl;
      Map_type_of_map::const_iterator i;
      for(i=map_type_of_map.begin();i!=map_type_of_map.end();i++)
        cerr << "\t " << *i->first.second << " <" << *i->first.first << ">" << "=" << *i->second << endl;        
      CompileError();
   }
   return r;
}
aType TypeArray(aType c,aType b,aType a)
{
   // type of  c[ b, a] 
   aType ba=map_pair_of_type[make_pair(b->right(),a->right())];
   if (!ba) {
      cerr << "Sorry is not possible to make a type of pair  "<< *b->right() << ", " << *c->right() << " " << endl;
      cerr << " list: " << endl;
      Map_type_of_map::const_iterator i;
      for(i=map_pair_of_type.begin();i!=map_pair_of_type.end();i++)
        cerr << "\t (" << *i->first.second << " , " << *i->first.first << ") " << "=" << *i->second << endl;        
      CompileError();
   }
   return TypeArray(c,ba);
}


inline  void ShowOn_cerr(const pair<const char * ,const OneOperator *> & i)
{ 
   cerr << "\t" <<  *i.first << ":" <<  endl;
   i.second->Show(cerr);
}
*/

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
	mi(tt.MeshIndependent()),
	tab(new Expression [N]),
	what(new int[N])  
      {
        assert(&tt);
      int err=0;
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
	mi(tt.MeshIndependent()),
	tab(new Expression [N]),
	what(new int[N])  
      {
        assert(&tt);
      int err=0;
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
   



/*
void ShowKeyWord(ostream & f ) 
 {
   zzzfff->dump(f);
 
 }

ostream* dumptable(ostream* f)
{

  *f << " the keywords " << endl;
  ShowKeyWord(*f);
  *f << " the types " << endl; 
  ShowType(*f);
   ListOfTOfId::const_iterator i=tables_of_identifier.begin();
   for(;i!=tables_of_identifier.end();++i)
    { 
      cout << " --------- table of identifier ---------\n";
      TableOfIdentifier * ti=*i;
      TableOfIdentifier::const_iterator mc=ti->m.begin();
      TableOfIdentifier::const_iterator end=ti->m.end();
      for (;mc != end;mc++)
       {
         *f  << "  - " << mc->first << ",  type :" <<  *mc->second.first << endl;
         const Polymorphic * op =dynamic_cast<const Polymorphic *>(mc->second.second) ;
         if ( op )  *f << *op << endl;
       }
      
     }
  
  return f;
}


long exec(string *s)
    {
      int r=execute(s->c_str());
      delete s;
      return r;}
      

*/ 
template<class K> long get_n(KN<K> * p){ return p->N();}
template<class K> long get_n(KNM<K> * p){ return p->N();}

template<class K> long get_m(KNM<K> * p){ return p->M();}
template<class K> K get_max(KN<K> * p){ return p->max();}
template<class K> K get_min(KN<K> * p){ return p->min();}
template<class K> K get_sum(KN<K> * p){ return p->sum();}
template<class K> K get_sum0(const KN_<K> & p){ return p.sum();}
template<class K> K get_max0(const KN_<K> &p){ return p.max();}
template<class K> K get_min0(const KN_<K> &p){ return p.min();}


/*
template<class F>  class dot1_F0
 { public:
   typedef typename F::argument_type pA;
   typedef typename F::Result R;
   
   F f;
   pA k;
   R operator()() const { return f(*k); }
   dot1_F0(F ff,pA kk):f(ff),k(kk) {}
 };
 
 
typedef  mem_fun_t<long,ostream> ostream_dot_f;

typedef dot1_F0<mem_fun_t<long,ostream> > ostream_prec;

//template <class S,class T> 
 template<class S,class T,S (T::*f)()>
 dot1_F0<mem_fun_t<S,T> >  Build(T*k) { return dot1_F0(f,k);}
 */
 
 class ostream_precis { public:
 ostream_precis(ostream * ff) :f(ff) {}
  ostream * f;
   operator long () const {return f->precision();}
 };
/*
 ostream_precis ostream_precision(ostream **f){ return ostream_precis(*f);}
  ostream_precis ostream_precision(ostream *f){ return ostream_precis(f);}
long get_precis( ostream_precis  pf) { return pf.f->precision();}
 long set_precis( ostream_precis  pf, long  l) { return pf.f->precision(l);}

 class istream_good { public:
  istream_good(istream * ff) :f(ff) {}
  istream * f;
  operator bool () const {return f->good();}
 };
 inline istream_good to_istream_good(istream **f){ return istream_good(*f);}
 inline istream_good to_istream_good(istream *f){ return istream_good(f);}
  
  inline long get_good( istream_good  pf) { return pf.f->good();}
  inline bool get_eof(istream ** p){ return (**p).eof();}
 
*/
 
template<class K>
void ArrayDCL()
{
    Dcl_TypeandPtr<KN<K> >(0,0,0,::Destroy<KN<K> >);
  //  Dcl_Type<KN<Complex> *>(0,::Destroy<KN<Complex> >);
    Dcl_Type<KN<K> *>(0,::Destroy<KN<K> >);
    Dcl_Type<KNM<K> *>(0,::Destroy<KNM<K> >);
    Dcl_Type< Transpose<KN<K> *> > ();
    
    Dcl_Type< Transpose<KN_<K> > > ();
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
    

     map_type[typeid(KN_<K> ).name()]->AddCast(
       new E_F1_funcT<KN_<K>,KN_<K>*>(UnRef<KN_<K> >),
       new E_F1_funcT<KN_<K>,KN<K>*>(UnRef<KN_<K>,KN<K> >)
    //   ,new E_F1_funcT<KN_<K>,K>(ValueToKN_<K>),
    //   new E_F1_funcT<KN_<K>,K*>(PtrToKN_<K>)       
       ); 
    map_type_of_map[make_pair(atype<long>(),atype<K>())]=atype<KN<K>*>(); // vector
    map_pair_of_type[make_pair(atype<long>(),atype<long>())] =atype<pair<long,long> >();   
    map_type_of_map[make_pair(atype<pair<long,long> >(),atype<K>())]=atype<KNM<K>*>(); // matrix                                               
}


template<class T> struct  Resize{ T *v;
  Resize( T * vv) : v(vv) {}
 }; 

template<class T> T *resize1(const Resize<T> & t,const long &n)
 {  
  t.v->resize(n);
  return  t.v;
 }
 
template<class T> T *resize2(const Resize<T> & t,const long &n, const long & m)
 {  
  t.v->resize(n,m);
  return  t.v;
 }

template<class T> Resize<T> to_Resize( T *v){ return Resize<T>(v);}

template<class K>
void ArrayOperator()
{
     Dcl_Type< Resize<KN<K> > > ();
     Dcl_Type< Resize<KNM<K> > > ();

     atype<KN<K>* >()->Add("[","",new OneOperator2_<K*,KN<K>*,long >(get_elementp_<K,KN<K>*,long>));
     atype<KN<K>* >()->Add("(","",new OneOperator2_<K*,KN<K>*,long >(get_elementp_<K,KN<K>*,long>));
     atype<KN_<K> >()->Add("(","",new OneOperator2_<KN_<K>,KN_<K>,char >(fSubArrayc<KN_<K> >));
     atype<KN_<K> >()->Add("(","",new OneOperator2_<KN_<K>,KN_<K>,SubArray>(fSubArray<K> ));
     atype<KN<K>*>()->Add("(","",new OneOperator2_<KN_<K>,KN<K>*,SubArray>(fSubArrayp<K> ));
     atype<KN<K>* >()->Add("(","",new OneOperator2_<KN<K>*,KN<K>*,char >(fSubArrayc<KN<K>* >));

     atype<KNM<K>* >()->Add("(","",new OneOperator3_<K*,KNM<K>*,long,long >(get_elementp2_<K,KNM<K>*,long,long>));

     Add<KN<K> *>("sum",".",new OneOperator1<K,KN<K> *>(get_sum));
     Add<KN<K> *>("min",".",new OneOperator1<K,KN<K> *>(get_min));
     Add<KN<K> *>("max",".",new OneOperator1<K,KN<K> *>(get_max));
     Add<KN_<K> >("sum",".",new OneOperator1_<K,KN_<K> >(get_sum0));
     Add<KN_<K> >("min",".",new OneOperator1_<K,KN_<K> >(get_min0));
     Add<KN_<K> >("max",".",new OneOperator1_<K,KN_<K> >(get_max0));

     Add<KN<K> *>("resize",".",new OneOperator1< Resize<KN<K> >,KN<K> *>(to_Resize));
     Add<KNM<K> *>("resize",".",new OneOperator1< Resize<KNM<K> >,KNM<K> *>(to_Resize));
     
     Add<Resize<KN<K> > >("(","",new OneOperator2_<KN<K> *,Resize<KN<K> > , long   >(resize1));
     Add<Resize<KNM<K> > >("(","",new OneOperator3_<KNM<K> *,Resize<KNM<K> > , long, long  >(resize2));

     TheOperators->Add("<-", 
       new OneOperator2_<KN<K> *,KN<K> *,long>(&set_init),
		       new InitArrayfromArray<K,true>
       );
     TheOperators->Add("<-", 
        new OneOperator3_<KNM<K> *,KNM<K> *,long,long>(&set_init2)
       );
       
     Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,long>(&set_init));
     Add<KNM<K> *>("<-","(",new OneOperator3_<KNM<K> *,KNM<K> *,long,long>(&set_init2));
     Add<KN<K> *>("<-","(",new InitArrayfromArray<K,true>);
     Add<KN<K> *>("n",".",new OneOperator1<long,KN<K> *>(get_n));
     Add<KNM<K> *>("n",".",new OneOperator1<long,KNM<K> *>(get_n));
     Add<KNM<K> *>("m",".",new OneOperator1<long,KNM<K> *>(get_m));
     
//     AddOpeqarray<set_eqarray,KN,K>("=");

     TheOperators->Add("=", new InitArrayfromArray<K,false>
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
     
     TheOperators->Add("=",
        new OneBinaryOperator<set_eqarray<KNM<K>  ,K > > 
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

     TheOperators->Add("/=",
        new OneBinaryOperator<set_eq_array_div<KN<K> ,K > > ,
        new OneBinaryOperator<set_eq_array_div<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_div<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp_div<KN_<K> ,KN<K>* > >        
     );

     TheOperators->Add("/=",
        new OneBinaryOperator<set_eqarray_div<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Mul_KNM_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_div<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarrayp_div<KN<K> ,KN<K>* > >        
     );

     TheOperators->Add("+",
       new OneBinaryOperator<Op2_addp<Add_KN_<K>,KN<K>*,KN<K>*> >,
       new OneBinaryOperator<Op2_add__n<Add_Mulc_KN_<K>,Mulc_KN_<K>,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_addp_n<Add_Mulc_KN_<K>,KN<K>*,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_add_pn<Add_Mulc_KN_<K>,Mulc_KN_<K> ,KN<K>* > >
       );

     TheOperators->Add("-",
       new OneBinaryOperator<Op2_subp<Sub_KN_<K>,KN<K>*,KN<K>*> >,
       new OneUnaryOperator<Op1_subp<Mulc_KN_<K>,KN<K>*> >,
       new OneBinaryOperator<Op2_sub__n<Add_Mulc_KN_<K>,Mulc_KN_<K>,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_subp_n<Add_Mulc_KN_<K>,KN<K>*,Mulc_KN_<K> > >,
       new OneBinaryOperator<Op2_sub_pn<Add_Mulc_KN_<K>,Mulc_KN_<K> ,KN<K>* > >
       );
     TheOperators->Add("*",
       new OneBinaryOperator<Op2_mulpc<Mulc_KN_<K>,KN<K>*,K> >,
       new OneBinaryOperator<Op2_mulcp<Mulc_KN_<K>,K,KN<K>*> >,
       new OneBinaryOperator<Op2_mulpcp<Mul_KNM_KN_<K>,KNM<K>*,KN<K>*> >,
       new OneBinaryOperator<Op2_dotproduct<K> >,
       new OneBinaryOperator<Op2_dotproduct_<K> >
       
       );

     TheOperators->Add("?:",
       new OneTernaryOperator3<Op3_p<if_arth_KN_<K>, KN_<K> > >       
       );

  TheOperators->Add("\'",       
       new OneOperator1<Transpose<KN<K> *>,KN<K> *>(&Build<Transpose<KN<K> *>,KN<K> *>),
       new OneOperator1<Transpose<KN_<K> >,KN_<K> >(&Build<Transpose<KN_<K> >,KN_<K> >)       
  );
       
     TheOperators->Add(".*",
       new OneBinaryOperator<Op2_dotstarp<DotStar_KN_<K>,KN<K>*,KN<K>*> >
      );
      
     TheOperators->Add("./",
       new OneBinaryOperator<Op2_dotstarp<DotSlash_KN_<K>,KN<K>*,KN<K>*> >
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
void initArrayOperators()
{
     ArrayOperator<double>();
     ArrayOperator<Complex>();
     ArrayOperator<long>();

}
void  initArrayDCL()
{
    ArrayDCL<double>();
    ArrayDCL<Complex>();
    ArrayDCL<long>();
}
