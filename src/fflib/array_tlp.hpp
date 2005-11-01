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

 
template<class K> 
struct Op2_dotproduct: public binary_function<Transpose<KN<K> *>,KN<K> *,K> { 
  static K f( Transpose<KN<K> *> const & a, KN<K> * const& b)  
   { return (conj(*a.t),*b);} }; 

template<class K> 
struct Op2_dotproduct_: public binary_function<Transpose<KN_<K> >,KN_<K> ,K> { 
  static K f( Transpose<KN_<K> > const & a, KN_<K>  const& b)  
   { return (conj(a.t),b);} }; 
   
template<class A,class B>  A Build(B b) {  return A(b);}
  
  
  
inline void MyAssert(int i,char * ex,char * file,long line)
{if (i) {
    cout << "CompileError assertion :  " << ex << " in file " << file << "  line = " << line << endl; 
     CompileError();}
 }


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

