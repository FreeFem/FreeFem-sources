//#pragma dont_inline on
//#pragma inline_depth(1)

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

Map_type_of_map map_type_of_map ; //  to store te type 
Map_type_of_map map_pair_of_type ; //  to store te type 

double  VersionNumber(); 
int TheCurrentLine=-1; // unset: by default
long mpisize=0,mpirank=0;
queue<pair<const E_Routine*,int> > debugstack;

bool showCPU= false;
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
template<class T> inline T Min (const T &a,const T & b){return a < b ? a : b;}
template<class T> inline T Abs (const T &a){return a <0 ? -a : a;}
template<class T> inline T Max (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}
template<class T> inline T Min (const T &a,const T & b,const T & c){return Min(Min(a,b),c);}
template<class T> inline T Square (const T &a){return a*a;}

struct SubArray2: public binary_function<long,long,SubArray> { 
  static SubArray f(const long & a,const long & b)  { 
   // cout << "SubArray: " << a << " " << b << endl;
    return SubArray(b-a+1,a);} }; 
struct SubArray3: public ternary_function<long,long,long,SubArray> { 
  static SubArray f(const long & a,const long & b,const long & c)  {  
  // cout << "SubArray: " << a << " " << b << " " <<  c << endl;
   return SubArray((b-a+1)/c,a,c);} }; 

#ifdef OLDCPP
template<class T> inline
complex<T> polar(const T& r, const T& theta)
{
	return complex<T>(r * cos(theta), r * sin(theta));
}
#endif


double preal( Complex * const& p){return real(*p);}
 
template<class K> 
struct Op2_dotproduct: public binary_function<Transpose<KN<K> *>,KN<K> *,K> { 
  static K f( Transpose<KN<K> *> const & a, KN<K> * const& b)  
   { return (*a.t,*b);} }; 

template<class K> 
struct Op2_dotproduct_: public binary_function<Transpose<KN_<K> >,KN_<K> ,K> { 
  static K f( Transpose<KN_<K> > const & a, KN_<K>  const& b)  
   { return (a.t,b);} }; 
   
template<class A,class B>  A Build(B b) {  return A(b);}
  
  


long Exit(long i) {throw(ErrorExec("Exit",i));return 0;}
bool Assert(bool b) {if (!b) throw(ErrorExec("exec assert",1));return true;}
  
inline void MyAssert(int i,char * ex,char * file,long line)
{if (i) {
    cout << "CompileError assertion :  " << ex << " in file " << file << "  line = " << line << endl; 
     CompileError();}
 }

  C_F0 *pOne=0,*pZero=0,*pminusOne=0;
// const C_F0 & One(*pOne), &Zero(*pZero);
 
 Polymorphic * TheOperators=new Polymorphic(), 
             * TheRightOperators=new Polymorphic();

TableOfIdentifier  Global;

 long E_Border::Count =0;

typedef list<TableOfIdentifier *> ListOfTOfId;

  ListOfTOfId tables_of_identifier;

  const int AC_F0::MaxSize=1024; // maximal number of parameters
  
  
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


map<const string,basicForEachType *> map_type;

template<class RR> RR LIncremantation(RR* a){ return ++(*a);}
template<class RR> RR RIncremantation(RR* a){ return (*a)++;}
template<class RR> RR LDecremantation(RR* a){ return --(*a);}
template<class RR> RR RDecremantation(RR* a){ return (*a)--;}

template<class RR,class B>
 RR * New_form_string(string * s) {B * r=  new B(s);delete *s;return r;}
 
 


typedef MyMap<String,double> mapSd ;
template<class K>
inline   K * get_element( MyMap<String,K> *  const  &  a,string*  const   & b)
 { K * ret=  &((*a)[*b]); // correction FH feb 2004
  //  cout << "get_element " << *b << " : " << ret << " = "<< * ret << endl;
    delete b;
    return ret;}

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

OneOperator::pair_find OneOperator::Find(const ArrayOfaType & at)const
 { 
      const OneOperator *w=0,*oo;
      int nn=0,p=0;
 /*     for (oo=this;oo;oo=oo->next)
        if (oo->pref>=p && oo->WithOutCast(at)) 
          {
           if(p<oo->pref) {nn=0;p=oo->pref;}
           nn++;
           w=oo;}
      if (nn) return make_pair(w,nn);*/
      for (int ncast=0;ncast<=n;ncast++) // loop on the number of cast 
       {
         p=0;
         for (oo=this;oo;oo=oo->next)
          if (oo->pref>=p && oo->WithCast(at,ncast)) 
          { 
           if(p<oo->pref) {nn=0;p=oo->pref;}
            nn++;
            w=oo;}
         if (nn) return make_pair(w,nn);
       }
      for (oo=this;oo;oo=oo->next)
        if (oo->WithCast(at)) 
          {nn++;
           w=oo;}
       return make_pair(w,nn);       
}

OneOperator::pair_find OneOperator::FindWithOutCast(const ArrayOfaType & at)const
 { 
      const OneOperator *w=0,*oo;
      int n=0;
      for (oo=this;oo;oo=oo->next)
        if (oo->WithOutCast(at)) 
          {n++;
           w=oo;}
      return make_pair(w,n);
}

OneOperator* OneOperator::FindSameR(const ArrayOfaType & at)
 { 
     if (!this) return 0;
      OneOperator *w=0,*oo,*r;
      int n=0;
      for (oo=this;oo;oo=oo->next)
        //if (oo->WithOutCast(at)) 
        if  (at==*oo)  n++,r=oo;        
        else if (oo->WithOutCast(at)) n++,r=oo;
     // if (n>1) cout << "FindSameR " << n << endl;
      return n==1 ? r : 0;
}
      
void OneOperator::Show(ostream &f) const
{         
   const OneOperator *w=0,*oo;
   for (oo=this;oo;oo=oo->next)
     f << "\t (" <<  *oo << ")\n";
 }   

void OneOperator::Show(const ArrayOfaType & at,ostream &f) const
{         
         const OneOperator *w=0,*oo;
         int n=0,np=0;
         for (oo=this;oo;oo=oo->next)
           if (oo->WithOutCast(at)) {n++;f << "\t (" <<  *oo << ")\n";}
         if(n==0) 
          for (oo=this;oo;oo=oo->next)
           if (oo->WithCast(at)) {
              n++;
              if (oo->pref) np++;
              if (oo->pref) 
                f <<   " c(" << oo->pref << ") \t (" <<  *oo << ")\n" ;
                else f <<  " \t c(" <<  *oo << ")\n";
              }
         if (n==0) 
          {
           f << " List of choises "<< endl;
           Show(f);            
          }
         else if (np != 1) 
           f << " We have ambigu•ty " << n << endl; 
 }   
       
const  OneOperator * Polymorphic::Find(const char *op, const  ArrayOfaType &at) const
  {
    const_iterator i=m.find(op);
    if (i!=m.end())  
      { 
       OneOperator::pair_find r=i->second->Find(at);
       if (r.second==1) return r.first;
       }    
    return 0;
  }
const  OneOperator * Polymorphic::FindWithOutCast(const char *op, const  ArrayOfaType &at) const
  {
    const_iterator i=m.find(op);
    if (i!=m.end())  
      { 
       OneOperator::pair_find r=i->second->FindWithOutCast(at);
       if (r.second==1) return r.first;
       }    
    return 0;
  }
 
  
void Polymorphic::Show(const char *op,const ArrayOfaType & at,ostream &f)  const
    {
    const_iterator i=m.find(op);
    if (i==m.end()) f << " unknow " << op << " operator " << endl;
    else i->second->Show(at,f);
  }

C_F0::C_F0(const Polymorphic * poly,const char *op,const basicAC_F0 & p)
{
  ArrayOfaType at(p); 
  if (poly) { // a Polymorphic => polymorphisme
      const  OneOperator *  ff=poly->Find(op,at);
      if (ff) { 
        /* cout << endl;
         poly->Show(op,at,cout);
         cout << op << ": (in " << at << ") => " << " " << *ff<< "\n\n";*/
         f=ff->code(p);
         r=*ff;}
      else
        { 
           cerr << " error operator " << op << " " << at << endl;
           poly->Show(op,at,cerr);
           CompileError();
        }
   }
  else { 
      //  no polymorphisme
     cerr << " const Polymorphic * poly,const char *op,const basicAC_F0 & p)   " << endl;
     cerr  << op << " " << p << endl;
     CompileError();          
  }
}
       



//  operator without parameter
C_F0::C_F0(const Polymorphic * pop,const char *op) 
{
  basicAC_F0  p;
  p=0;
  *this= C_F0(pop,op,p);
}
//  operator unaire
C_F0::C_F0(const Polymorphic * pop,const char *op,const C_F0 & aa) 
{
  basicAC_F0  p;
  C_F0 a(aa);
  p=a;
  *this= C_F0(pop,op,p);
}
//  operator binaire
C_F0::C_F0(const Polymorphic * pop,const char *op,const  C_F0 & a,const  C_F0 & b) 
{
  C_F0 tab[2]={a,b};
  basicAC_F0  p;
  p=make_pair<int,C_F0*>(2,tab);
  *this= C_F0(pop,op,p);
}

//  operator trinaire
C_F0::C_F0(const Polymorphic * pop,const char *op,const  C_F0 & a,const  C_F0 & b,const  C_F0 & c) 
{
  C_F0 tab[3]={a,b,c};
  basicAC_F0  p;
  p=make_pair<int,C_F0*>(3,tab);
  *this= C_F0(pop,op,p);
}


 OneOperator::~OneOperator(){ 
       OneOperator * d=next;
       next=0; 
       while(d) 
        { 
         OneOperator * dd=d->next;
         d->next=0;
         delete d;
         d=dd;
        }
  }

inline  void ShowOn_cerr(const pair<const char * ,const OneOperator *> & i)
{ 
   cerr << "\t" <<  *i.first << ":" <<  endl;
   i.second->Show(cerr);
}

void Polymorphic::Addp(const char * op,Value pp,...) const 
{
  pair<iterator,bool>  p=m.insert(make_pair<const Key,Value>(op,pp));
  Value f= p.first->second;	
  if (!p.second)  // not insert => old 
    *f += *pp;
  va_list ap;
  va_start(ap,pp);
  for(pp=va_arg(ap,OneOperator * );pp;pp=va_arg(ap,OneOperator * ))
    *f += *pp;
/*  if ( ! strlen(op) )  
   { // no polymorphisme
     if(m.size() !=1 ||  !f->Simple()) {
       cerr << " no polymorphisme and polymorphisme are mixed " << endl;
    //   for_each(m.begin,m.end(),ShowOn_cerr);
       CompileError();       
     }   
   } */
}

void Polymorphic::Add(const char * op,Value *pp) const
{
  if (*pp)
   {
    pair<iterator,bool>  p=m.insert(make_pair<const Key,Value>(op,*pp));
    Value f= p.first->second;	
    if (!p.second)  // not insert => old 
      *f += **pp;
    pp++;
    for(;*pp;pp++)
     *f += **pp;
   /*if ( ! strlen(op) )  
     { // no polymorphisme
      if(m.size() !=1 ||  !f->Simple()) {
       cerr << " no polymorphisme and polymorphisme are mixed " << endl;
      //   for_each(m.begin,m.end(),ShowOn_cerr);
       CompileError();       
     }  } */ }
     
}


 int  FindType(const char * name)  
   {
   C_F0 r;

     ListOfTOfId::const_iterator i=tables_of_identifier.begin();
      for(;i!=tables_of_identifier.end();++i)
      { 
      TableOfIdentifier * ti=*i;     
      r = ti->Find(name);
      if (r.NotNull()) return r.TYPEOFID();
    }
     return 0;
   } 

  C_F0 Find(const char * name)   
{
   C_F0 r;
   ListOfTOfId::const_iterator i=tables_of_identifier.begin();
   for(;i!=tables_of_identifier.end();++i)
    { 
      TableOfIdentifier * ti=*i;     
      r = ti->Find(name);
      if (r.NotNull()) return r;
    }
    cerr << " The Identifier " << name << " does not exist " << endl;
    CompileError();
    return r;
}

C_F0 TableOfIdentifier::destroy()
{
 int k=0;
// cout << "\n\t List of destroy variables " << m.size() << " : " ; 
 for (pKV * i=listofvar;i;i=i->second.next)
   {
     if  (i->second.del && i->second.first->ExistDestroy() )
    // cout  << i->first << ", " ;
     assert(i->second.first);
     if (i->second.del && i->second.first->ExistDestroy() ) k++;
   }
// cout << endl;
 ListOfInst *l=new ListOfInst(k);
 for (pKV * i=listofvar;i;i=i->second.next)
     if (i->second.del && i->second.first->ExistDestroy()) 
       l->Add(i->second.first->Destroy(i->second) );
   
 return C_F0(l);     
}


Expression basicForEachType::Destroy(const C_F0 & e) const 
{
    return destroy ? NewExpression(destroy,e) : (Expression)  e;
} 

 void basicForEachType::SetArgs(const ListOfId *lid) const
{ SHOWVERB(cout << "SetArgs::\n ") ;throwassert(lid==0 || lid->size()==0);}


const  Type_Expr &   TableOfIdentifier::New(Key k,const Type_Expr & v,bool del)
 {
    pair<iterator,bool>  p=m.insert(pKV(k,Value(v,listofvar,del)));
    listofvar = &*m.find(k);
    if (!p.second) 
     {
       cerr << " The identifier " << k << " exist \n";
       cerr << " \t  the existing type is " << *p.first->second.first << endl;
       cerr << " \t  the new  type is " << *v.first << endl;
       CompileError();
     }
    return v;
 }
 void  TableOfIdentifier::Add(Key k,Key op,OneOperator *p0,OneOperator *p1,
      OneOperator *p2,OneOperator *p3,OneOperator *p4,OneOperator *p5,OneOperator *p6)
{
  iterator i= m.find(k);
  if (i==m.end()) // new
    {
     Value poly0=Value(atype<Polymorphic*>(),new Polymorphic(),listofvar);     
    i=m.insert(make_pair<const Key,Value>(k,poly0)).first;
    listofvar= &*i;
    }
    const Polymorphic * p= dynamic_cast<const Polymorphic *>(i->second.second);
  if ( !p) { 
    cerr << k << " is not a Polymorphic id " << endl;
    CompileError();
  }
  p->Add(op,p0,p1,p2,p3,p4,p5,p6);
}

 ArrayOfaType::ArrayOfaType(const ListOfId * l)
  : n(l->size()),t(new aType[n]),ellipse(false)
 {
    for (int i=0;i<n;i++) 
      {
      t[i]=(*l)[i].r;  
       if ( ! t[i])
        {
           cerr << " Argument " << i << " '"<< (*l)[i].id << "' without type\n";
           CompileError("DCL routine: Argument without type ");
         }
      }
 }

bool ArrayOfaType::WithOutCast( const ArrayOfaType & a) const 
 {
   if ( ( !ellipse && (a.n != n))  || (ellipse && n > a.n) ) return false;
   for (int i=0;i<n;i++)
       if (! a.t[i]->SametypeRight(t[i])) 
        return false;       
 // cerr << " TRUE " << endl;    
   return true;
 }

  
bool ArrayOfaType::WithCast( const ArrayOfaType & a,int nbcast) const 
 {  
   if (  ( !ellipse && (a.n != n))  || (ellipse && n > a.n) ) return false;
   for (int i=0;i<n;i++)
     if ( a.t[i]->SametypeRight(t[i])) ;
     else if (! t[i]->CastingFrom(a.t[i])) return false; 
     else if ( --nbcast <0) return false;
   return true;
 }    
 
void basicForEachType::AddCast(CastFunc f1,CastFunc f2,CastFunc f3,CastFunc f4,
  CastFunc f5,CastFunc f6,CastFunc f7,CastFunc f8)
 {
   CastFunc ff[]={f1,f2,f3,f4,f5,f6,f7,f8,0};
   for (int i=0;ff[i];i++)
    {
      throwassert(this == *ff[i] );
      if (casting->FindSameR(*ff[i]))
       {
         cerr << " The casting to " << *ff[i] << " exist " << endl;
         cerr << " List of cast " << endl;
         casting->Show(cerr);
         CompileError();
       }
      if (casting)  *casting += *ff[i];
      else casting = ff[i];
/*      
   if( ! mapofcast.insert(make_pair<const aType,CastFunc>(ff[i]->a,ff[i])).second) 
     {
        cerr << " The casting to "<< *this << " from " << ff[i]->a << " exist " << endl;
        cerr << " List of cast " << endl;
        for_each(mapofcast.begin(),mapofcast.end(),CerrCast);
        CompileError();
     } */
  }
}

 ostream & operator<<(ostream & f,const OneOperator & a)     
{
//   for(const OneOperator * tt=&a;tt;tt=tt->next)
     f << "\t  " << * (a.r) << " :  "  <<(const ArrayOfaType &) a;
   return f;   
}

 ostream & operator<<(ostream & f,const Polymorphic & a)
{
  Polymorphic::const_iterator i;
  if(!&a) return f << "Null " << endl;
  for (i=a.m.begin();i!=a.m.end();i++)
   {
    f << "   operator" << i->first << " : " << endl;
    i->second->Show(f);
   }
  return f;
}  
 ostream & operator<<(ostream & f,const ArrayOfaType & a)
   { 
     for (int i=0;i<a.n;i++) 
       f <<  (i ? ", " : " ") << *a.t[i];
       if (a.ellipse ) f << "... ";
       else            f << " "; 
      return f;}
    ostream & operator<<(ostream & f,const TableOfIdentifier & t )
 {
   TableOfIdentifier::const_iterator i;
   for(i=t.m.begin();i!=t.m.end();i++)
    {
      TableOfIdentifier::Value v=i->second;
      f << i->first << ":  " << *v.first << " <- " ; 
      const Polymorphic * p=dynamic_cast<const Polymorphic *>(v.second);
      if(p) f << "Polymorphic " << *p << endl;
      else  f << " Simple @" <<  v.second << endl;
    }
    return f;
 }
 
Expression NewExpression(Function1 f,Expression a)
{
  throwassert(f);
  return new E_F0_Func1(f,a);
}
Expression NewExpression(Function2 f,Expression a,Expression b)
{
  throwassert(f);
  return new E_F0_Func2(f,a,b);
 
}
 
 void ShowType(ostream & f)
 {
 
   map<const string,basicForEachType *>::const_iterator i;
   for(i=map_type.begin();i!=map_type.end();i++)
     {
       f << " --"<< i->first <<" = " ; 
       i->second->Show(f) ;
       f << endl;
     }

 }

 void basicForEachType::Show(ostream & f) const {
       f << " " <<* this << endl;
       if (casting) f <<*casting ;
       if (ti.m.size())
        { 
          TableOfIdentifier::const_iterator mc=ti.m.begin();
          TableOfIdentifier::const_iterator end=ti.m.end();
          for (;mc != end;mc++)
          {
            f  << "    " << mc->first << ",  type :" <<  *mc->second.first << endl;
            const Polymorphic * op =dynamic_cast<const Polymorphic *>(mc->second.second) ;
            if ( op )  f << *op << endl;
          }
        }
   }
 
//size_t CodeAlloc::nb=0;
//size_t CodeAlloc::lg=0;
/*
FuncForEachType::FuncForEachType(const basicForEachType * t)
: rtype(t),basicForEachType(t->,sizeof((void*)()),
         0,0,0,0){ throwassert(t);}
*/

E_Routine::E_Routine(const Routine * routine,const basicAC_F0 & args)
 : rt(routine->tret),code(routine->ins),clean(routine->clean),
   nbparam(args.size()),param(new Expression[nbparam]),name(routine->name)
{    
   assert(routine->ins); 
   for (int i=0;i<args.size();i++)
        param[i]=routine->param[i].r->CastTo(args[i]);
};

AnyType E_Routine::operator()(Stack s)  const  {
   debugstack.push(make_pair<const E_Routine*,int>(this,TheCurrentLine));
   const int lgsave=BeginOffset*sizeof(void*);
   char  save[lgsave];
   AnyType ret=Nothing;
   memcpy(save,s,lgsave); // save register 
   AnyType *listparam=new AnyType[nbparam];
   for (int i=0;i<nbparam;i++)
     listparam[i]= (*param[i])(s); // set of the parameter 
   Stack_Ptr<AnyType>(s,ParamPtrOffset) = listparam; 
   try {  (*code)(s);  }
   catch( E_exception & e) { 
           (*clean)(s); 
          // cout << " catch " << e.what() << " clean & throw " << endl;
           if (e.type() == E_exception::e_return)  
              ret = e.r;
           else 
              ErrorExec("E_exception: break or contine not in loop ",1);
        }
        
   delete [] listparam; 
    memcpy(s,save,lgsave);  // restore register
    TheCurrentLine=debugstack.front().second;
    debugstack.pop();
 
   return ret;
}
class E_F0para :public E_F0 { public:
  const int i;
  AnyType operator()(Stack s)  const  {
  //  return  (* Stack_Ptr<Expression>(s,ParamPtrOffset)[i])(s);
    return Stack_Ptr<AnyType>(s,ParamPtrOffset)[i];
  }
   E_F0para(int ii) : i(ii){}
};        

template<class RR>
class  InitArrayfromArray : public OneOperator { public:
    typedef KN<RR> * A;
    typedef KN<RR> * R;
    typedef E_Array B;
    
    class CODE : public  E_F0 { public:
       Expression a0;
       int N;
       Expression * tab;
       const  bool mi;
    CODE(Expression a,const E_Array & tt)  :
     a0(a),N(tt.size()),mi(tt.MeshIndependent()),tab(new Expression [N]) 
      {
        assert(&tt);
        
        for (int i=0;i<N;i++)
          tab[i]=atype<RR>()->CastTo(tt[i]);
      }
    AnyType operator()(Stack s)  const 
    {
      A  a=GetAny<A>((*a0)(s));
      a->init(N);
      for (int i=0;i<N;i++)
        (*a)[i]= GetAny<RR>((*(tab[i]))(s));
      return SetAny<R>(a);} 
    bool MeshIndependent() const 
      {return  mi;} // 
    operator aType () const { return atype<R>();}    
    }; 
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new CODE(t[0]->CastTo(args[0]),*dynamic_cast<const E_Array*>( t[1]->CastTo(args[1]).LeftValue()));} 
    InitArrayfromArray():   OneOperator(atype<R>(),atype<A>(),atype<B>())  {}
};
   
Routine::Routine(aType tf,aType tr,const char * iden, const ListOfId *l,Block * & cb)
    : OneOperator(tr,l),offset(cb->OffSet(sizeof(void*))),
     tfunc(tf),tret(tr),name(iden),param(*l),
      currentblock(new Block(cb)),ins(0),clean(0) 
     {
       cb = currentblock; 
       for (int i=0;i<param.size();i++)
           currentblock->NewID(param[i].r,param[i].id,C_F0(new E_F0para(i),param[i].r),!param[i].ref);
     }
   Block * Routine::Set(C_F0 instrs) 
       { 
         ins=instrs;
         clean = (C_F0) currentblock->close(currentblock);
         return    currentblock;} 
         
 
E_F0 * Routine::code(const basicAC_F0 & args) const 
{
   
   return new E_Routine(this,args);
}

extern void ShowKeyWord(ostream & f ); 
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
      

 
template<class K> long get_n(KN<K> * p){ return p->N();}
template<class K> long get_n(KNM<K> * p){ return p->N();}
template<class K> long get_m(KNM<K> * p){ return p->M();}
template<class K> double get_max(KN<K> * p){ return p->max();}
template<class K> double get_min(KN<K> * p){ return p->min();}
template<class K> K get_sum(KN<K> * p){ return p->sum();}


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
 
 ostream_precis ostream_precision(ostream **f){ return ostream_precis(*f);}
  ostream_precis ostream_precision(ostream *f){ return ostream_precis(f);}
long get_precis( ostream_precis  pf) { return pf.f->precision();}
 long set_precis( ostream_precis  pf, long  l) { return pf.f->precision(l);}
 
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
    Dcl_Type<Add_Mulc_KN_<K> *>();

     map_type[typeid(KN_<K> ).name()]->AddCast(
       new E_F1_funcT<KN_<K>,KN_<K>*>(UnRef<KN_<K> >),
       new E_F1_funcT<KN_<K>,KN<K>*>(UnRef<KN_<K>,KN<K> >)
       ); 
    map_type_of_map[make_pair(atype<long>(),atype<K>())]=atype<KN<K>*>(); // vector
    map_pair_of_type[make_pair(atype<long>(),atype<long>())] =atype<pair<long,long> >();   
    map_type_of_map[make_pair(atype<pair<long,long> >(),atype<K>())]=atype<KNM<K>*>(); // matrix                                               
}

template<class K>
void ArrayOperator()
{

     atype<KN<K>* >()->Add("[","",new OneOperator2_<K*,KN<K>*,long >(get_elementp_<K,KN<K>*,long>));
     atype<KN<K>* >()->Add("(","",new OneOperator2_<K*,KN<K>*,long >(get_elementp_<K,KN<K>*,long>));
     atype<KN_<K> >()->Add("(","",new OneOperator2_<KN_<K>,KN_<K>,char >(fSubArrayc<KN_<K> >));
     atype<KN_<K> >()->Add("(","",new OneOperator2_<KN_<K>,KN_<K>,SubArray>(fSubArray<K> ));
     atype<KN<K>*>()->Add("(","",new OneOperator2_<KN_<K>,KN<K>*,SubArray>(fSubArrayp<K> ));
     atype<KN<K>* >()->Add("(","",new OneOperator2_<KN<K>*,KN<K>*,char >(fSubArrayc<KN<K>* >));

     atype<KNM<K>* >()->Add("(","",new OneOperator3_<K*,KNM<K>*,long,long >(get_elementp2_<K,KNM<K>*,long,long>));

     Add<KN<K> *>("sum",".",new OneOperator1<K,KN<K> *>(get_sum));

     TheOperators->Add("<-", 
       new OneOperator2_<KN<K> *,KN<K> *,long>(&set_init),
       new InitArrayfromArray<K>
       );
     TheOperators->Add("<-", 
        new OneOperator3_<KNM<K> *,KNM<K> *,long,long>(&set_init2)
       );
       
     Add<KN<K> *>("<-","(",new OneOperator2_<KN<K> *,KN<K> *,long>(&set_init));
     Add<KNM<K> *>("<-","(",new OneOperator3_<KNM<K> *,KNM<K> *,long,long>(&set_init2));
     Add<KN<K> *>("<-","(",new InitArrayfromArray<K>);
     Add<KN<K> *>("n",".",new OneOperator1<long,KN<K> *>(get_n));
     Add<KNM<K> *>("n",".",new OneOperator1<long,KNM<K> *>(get_n));
     Add<KNM<K> *>("m",".",new OneOperator1<long,KNM<K> *>(get_m));
     
     TheOperators->Add("=",
        new OneBinaryOperator<set_eqarray<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd<KN<K> ,Add_Mulc_KN_<K>* > > ,
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
        new OneBinaryOperator<set_eq_array<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp<KN_<K> ,KN<K>* > >       
      );

     TheOperators->Add("+=",
        new OneBinaryOperator<set_eqarray_add<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_add<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_add<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarrayp_add<KN<K> ,KN<K>* > >        
      );
     TheOperators->Add("+=",
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,K > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_add<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_add<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp_add<KN_<K> ,KN<K>* > >        
      );
      
     TheOperators->Add("-=",
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_sub<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_sub<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarrayp_sub<KN<K> ,KN<K>* > >        
      );
      
     TheOperators->Add("-=",
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,K > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,DotStar_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,DotSlash_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_sub<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_sub<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp_sub<KN_<K> ,KN<K>* > >        
      );
      
     TheOperators->Add("*=",
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,K > >  ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_mul<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eqarraypd_mul<KN<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eqarrayp_mul<KN<K> ,KN<K>* > >       
      );
 
      TheOperators->Add("*=",
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,K > >  ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_mul<KN_<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_mul<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp_mul<KN_<K> ,KN<K>* > >       
      );

     TheOperators->Add("/=",
        new OneBinaryOperator<set_eq_array_div<KN<K> ,K > > ,
        new OneBinaryOperator<set_eq_array_div<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eq_array_div<KN<K> ,Mulc_KN_<K> > > ,
        new OneBinaryOperator<set_eq_arraypd_div<KN_<K> ,Add_Mulc_KN_<K>* > > ,
        new OneBinaryOperator<set_eq_arrayp_div<KN_<K> ,KN<K>* > >        
     );

     TheOperators->Add("/=",
        new OneBinaryOperator<set_eqarray_div<KN<K> ,K > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Add_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Sub_KN_<K> > > ,
        new OneBinaryOperator<set_eqarray_div<KN<K> ,Mulc_KN_<K> > > ,
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
       new OneBinaryOperator<Op2_dotproduct<K> >,
       new OneBinaryOperator<Op2_dotproduct_<K> >
       
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
void Init_map_type()
{
  //  cout << sizeof(string) << endl;
    map_type[typeid(AnyType).name()] = new ForTypeAnyType();
    map_type[typeid(void).name()] = new ForTypeVoid();

    Dcl_Type<Expression>(0);    
    Dcl_TypeandPtr<double>(0);
    Dcl_TypeandPtr<long>(0);
    Dcl_TypeandPtr<bool>(0);
    Dcl_TypeandPtr<Complex>(0);
    Dcl_Type<char*>();
    Dcl_Type<const char *>();
    Dcl_Type<char>();
    Dcl_TypeandPtr<string*>(0,::Delete<string>,::InitializePtr<string*>,::DeletePtr<string*>);
    Dcl_TypeandPtr<ostream*>(0,0,::InitializePtr<ostream*>,::DeletePtr<ostream*>);
    Dcl_TypeandPtr<istream*>(0,0,::InitializePtr<istream*>,::DeletePtr<istream*>);
    Dcl_Type< ostream_precis > ();
    
    Dcl_Type<Polymorphic*>();
    
//    Dcl_Type<C_F0>();
    map_type[typeid(C_F0).name()] =  new TypeLineFunction; 
    Dcl_Type<E_Array>();
    Dcl_Type<const E_Border *>();
    Dcl_Type<const E_BorderN *>();

    
    
    Dcl_Type<SubArray>();
    Dcl_Type<pair<long,long> >();
    
    ArrayDCL<double>();
    ArrayDCL<Complex>();
    
    Dcl_Type<ios::openmode>();
    
//  les types des variables 
    
  zzzfff->Add("real",atype<double*>());
  zzzfff->Add("int",atype<long*>());
  zzzfff->Add("complex",atype<Complex*>());
  zzzfff->Add("bool",atype<bool*>());
  zzzfff->Add("string",atype<string**>());
  zzzfff->Add("ifstream",atype<istream**>());
  zzzfff->Add("ofstream",atype<ostream**>());
  zzzfff->AddF("func",atype<C_F0>());
  
  
    
//  end of know types

     map_type[typeid(bool).name()]->AddCast(
       new E_F1_funcT<bool,bool*>(UnRef<bool>),
       new E_F1_funcT<bool,long>(Cast<bool,long>),
       new E_F1_funcT<bool,double>(Cast<bool,double>)
       );                                         
                        

     map_type[typeid(long).name()]->AddCast(
       new E_F1_funcT<long,long*>(UnRef<long>),
       new E_F1_funcT<long,double>(Cast<long,double>),
       new E_F1_funcT<long,bool>(Cast<long,bool>),
       new E_F1_funcT<long,ostream_precis>(Cast<long,ostream_precis>)

       );
       
       
     map_type[typeid(double).name()]->AddCast(
       new E_F1_funcT<double,double*>(UnRef<double>),
       new E_F1_funcT<double,long>(Cast<double,long>),
       new E_F1_funcT<double,bool>(Cast<double,bool>)
       
       ); 
                                               
     map_type[typeid(Complex).name()]->AddCast(
       new E_F1_funcT<Complex,Complex*>(UnRef<Complex>),
       new E_F1_funcT<Complex,long>(Cast<Complex,long>),
       new E_F1_funcT<Complex,double>(Cast<Complex,double>)
       );                                         

     map_type[typeid(string*).name()]->AddCast(
       new E_F1_funcT<string*,string**>(UnRefCopyPtr<string>),
       new E_F1_funcT<string*,long>(FCast<string*,long,toString>),
       new E_F1_funcT<string*,double>(FCast<string*,double,toString>),
       new E_F1_funcT<string*,bool>(FCast<string*,bool,toString>),
       new E_F1_funcT<string*,Complex>(FCast<string*,Complex,toString>)
        ); 
     //   a changer ---------------  modif                                         
        map_type[typeid(string*).name()]->AddCast(					       
          new E_F1_funcT<string*,char *>(FCast<string*,char *,toStringC>),
          new E_F1_funcT<string*,const char *>(FCast<string* ,const char *,toStringCconst>)
       );
       
     map_type[typeid(long).name()]->AddCast(new OneOperator_border_label);
      
     Global.New("verbosity",CPValue<long>(verbosity));
     
     Global.New("cout",CConstant<ostream*>(&cout));
     Global.New("cin",CConstant<istream*>(&cin));
     Global.New("append",CConstant<ios::openmode>(ios::app));
     Global.New("endl",CConstant<char*>("\n"));
     Global.New("true",CConstant<bool>(true));
     Global.New("false",CConstant<bool>(false));
     Global.New("pi",CConstant<double>(3.14159265358979323846264338328));
     Global.New("version",CConstant<double>(VersionNumber()));
      
     Global.New("CPUTime",CConstant<bool*>(&showCPU));
     // def de Zero et One
     pZero = new  C_F0(CConstant<double>(0.0));
     pOne = new  C_F0(CConstant<double>(1.0));
     pminusOne = new  C_F0(CConstant<double>(-1.0));
     
     TheOperators->Add(":",
       new OneOperatorConst<char>(new EConstant<char>(':')),
       new OneBinaryOperator<SubArray2>,
       new OneTernaryOperator3<SubArray3>);
     
       
     TheOperators->Add("+",
       new OneBinaryOperator<Op2_add<long,long,long> >,
       new OneBinaryOperator<Op2_add<double,double,double> >,
       new OneBinaryOperator<Op2_add<double,double,long> >,
       new OneBinaryOperator<Op2_add<double,long,double> >,
       new OneBinaryOperator<Op2_add<long,bool,bool> >,       
       new OneBinaryOperator<Op2_add<long,long,bool> >,       
       new OneBinaryOperator<Op2_add<long,bool,long> >,       
       new OneBinaryOperator<Op2_add<Complex,Complex,Complex> >,
       new OneBinaryOperator<Op2_add<Complex,Complex,double> >,
       new OneBinaryOperator<Op2_add<Complex,double,Complex> >,
       new OneBinaryOperator<Op2_add<Complex,Complex,long> >,
       new OneBinaryOperator<Op2_add<Complex,long,Complex> > ,      
       new OneBinaryOperator<Op2_padd<string,string*,string*> >       
       );
     TheOperators->Add("-",
       new OneBinaryOperator<Op2_sub<long,long,long> >,
       new OneBinaryOperator<Op2_sub<double,double,double> >,
       new OneBinaryOperator<Op2_sub<double,double,long> >,
       new OneBinaryOperator<Op2_sub<double,long,double> >,
       new OneBinaryOperator<Op2_sub<long,bool,bool> >,              
       new OneBinaryOperator<Op2_sub<Complex,Complex,Complex> >,
       new OneBinaryOperator<Op2_sub<Complex,Complex,double> >,
       new OneBinaryOperator<Op2_sub<Complex,double,Complex> >,
       new OneBinaryOperator<Op2_sub<Complex,Complex,long> >,
       new OneBinaryOperator<Op2_sub<Complex,long,Complex> >       
       );
     TheOperators->Add("*",
       new OneBinaryOperator<Op2_mul<long,long,long> >,
       new OneBinaryOperator<Op2_mul<double,double,double> >,
       new OneBinaryOperator<Op2_mul<double,double,long> >,
       new OneBinaryOperator<Op2_mul<double,long,double> >,
       new OneBinaryOperator<Op2_mul<Complex,Complex,Complex> >,
       new OneBinaryOperator<Op2_mul<Complex,Complex,double> >,
       new OneBinaryOperator<Op2_mul<Complex,double,Complex> >,
       new OneBinaryOperator<Op2_mul<Complex,Complex,long> >,
       new OneBinaryOperator<Op2_mul<Complex,long,Complex> >       
       );
     TheOperators->Add("/",
       new OneBinaryOperator<Op2_div<long,long,long> >,
       new OneBinaryOperator<Op2_div<double,double,double> >,
       new OneBinaryOperator<Op2_div<double,double,long> >,
       new OneBinaryOperator<Op2_div<double,long,double> >,
       new OneBinaryOperator<Op2_div<Complex,Complex,Complex> >,
       new OneBinaryOperator<Op2_div<Complex,Complex,double> >,
       new OneBinaryOperator<Op2_div<Complex,double,Complex> >,
       new OneBinaryOperator<Op2_div<Complex,Complex,long> >,
       new OneBinaryOperator<Op2_div<Complex,long,Complex> >       
       );

     TheOperators->Add("%",
       new OneBinaryOperator<Op2_mod<long,long,long> >
       );


    TheOperators->Add("+",
       new OneUnaryOperator<Op1_plus<double> >,
       new OneUnaryOperator<Op1_plus<long> >,
       new OneUnaryOperator<Op1_plus<Complex> >);
     
    TheOperators->Add("-",
       new OneUnaryOperator<Op1_neg<double> >,
       new OneUnaryOperator<Op1_neg<long> >,
       new OneUnaryOperator<Op1_neg<Complex> >);
              
     TheOperators->Add("^",
       new OneBinaryOperator<Op2_pow<long,long,long> >,
    //   new OneBinaryOperator<Op2_pow<double,long,double> >,
       new OneBinaryOperator<Op2_pow<double,double,double> >,
       new OneBinaryOperator<Op2_pow<double,double,long> >,
    //   new OneBinaryOperator<Op2_pow<Complex,Complex,double> >,
    //  new OneBinaryOperator<Op2_pow<Complex,double,Complex> >,
       new OneBinaryOperator<Op2_pow<Complex,Complex,Complex> >
     );
     
     TheOperators->Add("<",
       new OneBinaryOperator<Op2_lt<long,long> >,
       new OneBinaryOperator<Op2_lt<double,double> >,
       new OneBinaryOperator<Op2_plt<string*,string*> >
     );
     TheOperators->Add("<=",
       new OneBinaryOperator<Op2_le<long,long> >,
       new OneBinaryOperator<Op2_le<double,double> >,
       new OneBinaryOperator<Op2_ple<string*,string*> >
     );
     TheOperators->Add(">",
       new OneBinaryOperator<Op2_gt<long,long> >,
       new OneBinaryOperator<Op2_gt<double,double> >,
       new OneBinaryOperator<Op2_pgt<string*,string*> >
     );
     TheOperators->Add(">=",
       new OneBinaryOperator<Op2_ge<long,long> >,
       new OneBinaryOperator<Op2_ge<double,double> >,
       new OneBinaryOperator<Op2_pge<string*,string*> >
     );
     TheOperators->Add("==",
       new OneBinaryOperator<Op2_eq<long,long> >,
       new OneBinaryOperator<Op2_eq<double,double> >,
       new OneBinaryOperator<Op2_eq<Complex,Complex> >,
       new OneBinaryOperator<Op2_peq<string*,string*> >
     );

     TheOperators->Add("!=",
       new OneBinaryOperator<Op2_ne<long,long> >,
       new OneBinaryOperator<Op2_ne<double,double> >,
       new OneBinaryOperator<Op2_ne<Complex,Complex> >,
       new OneBinaryOperator<Op2_pne<string*,string*> >
     );
     
     TheOperators->Add("!",
       new OneUnaryOperator<Op1_not<bool > >
     );
     
     TheOperators->Add("&&", new OneBinaryOperator<Op2_and > );
     TheOperators->Add("&", new OneBinaryOperator<Op2_and > );
     TheOperators->Add("||", new OneBinaryOperator<Op2_or> );
     TheOperators->Add("|", new OneBinaryOperator<Op2_or> );
       
      // Unary_Op_Comparaision
     
     TheOperators->Add("=",
       new OneBinaryOperator<set_eq<bool> >,
       new OneBinaryOperator<set_eq<long> >,
       new OneBinaryOperator<set_eq<double> >,
       new OneBinaryOperator<set_eq<Complex> >,
       new OneBinaryOperator<set_peq<string*> >
       ); 
       
     ArrayOperator<double>();
     ArrayOperator<Complex>();
     
 

     TheOperators->Add("+=",
       new OneBinaryOperator<set_eq_add<long> >,
       new OneBinaryOperator<set_eq_add<double> >,
       new OneBinaryOperator<set_eq_add<Complex> >
      );


     TheOperators->Add("-=",
       new OneBinaryOperator<set_eq_sub<long> >,
       new OneBinaryOperator<set_eq_sub<double> >,
       new OneBinaryOperator<set_eq_sub<Complex> >
      );



     TheOperators->Add("*=",
       new OneBinaryOperator<set_eq_mul<long> >,
       new OneBinaryOperator<set_eq_mul<double> >,
       new OneBinaryOperator<set_eq_mul<Complex> >     
      );


     TheOperators->Add("/=",
       new OneBinaryOperator<set_eq_div<long> >,
       new OneBinaryOperator<set_eq_div<double> >,
       new OneBinaryOperator<set_eq_div<Complex> >     
     );

     TheOperators->Add("+",
    //   new OneBinaryOperator<Op2_addp<const E_BorderN *,const E_BorderN *,const E_BorderN * > >,  
       new AddBorderOperator
       );

      

      

     TheOperators->Add(">>",
       new OneBinaryOperator<Op_Read<bool> >,
       new OneBinaryOperator<Op_Read<long> >,
       new OneBinaryOperator<Op_Read<double> >,
       new OneBinaryOperator<Op_Read<Complex> >
       );
     
     TheOperators->Add("<<",
       new OneBinaryOperator<Print<bool> >,
       new OneBinaryOperator<Print<long> >,
       new OneBinaryOperator<Print<double> >,
       new OneBinaryOperator<Print<Complex> >,
       new OneBinaryOperator<PrintP<string*> >
       );

     
     TheRightOperators->Add("++",       
       new OneOperator1<long,long*>(&RIncremantation<long>));
     TheRightOperators->Add("--",       
       new OneOperator1<long,long*>(&RDecremantation<long>));
     TheOperators->Add("++",       
       new OneOperator1<long,long*>(&LIncremantation<long>));
     TheOperators->Add("--",       
       new OneOperator1<long,long*>(&LDecremantation<long>));
//   init        
     TheOperators->Add("<-", 
       new OneOperator2_<string**,string**,string*>(&set_copy),
       new OneOperator2_<double*,double*,double>(&set_copy),
       new OneOperator2_<long*,long*,long>(&set_copy),
       new OneOperator2_<bool*,bool*,bool>(&set_copy),
       new OneOperator2_<Complex*,Complex*,Complex>(&set_copy),
       new OneOperator2_<istream**,istream**,istream*>(&set_copy),
       new OneOperator2_<ostream**,ostream**,ostream*>(&set_copy),
//       new OneUnaryOperator<Op1_new_pstring<istream*,ifstream> >,
//       new OneUnaryOperator<Op1_new_pstring<ostream*,ofstream> >,
       new OneBinaryOperator<Op2_set_pstring<istream**,ifstream> >,
       new OneBinaryOperator<Op2_set_pstring<ostream**,ofstream> >,
       new OneTernaryOperator3<Op2_set_pstringiomode<ostream**,ofstream> >       
       );  
       
     atype<istream* >()->AddCast( new E_F1_funcT<istream*,istream**>(UnRef<istream* >)); 
     atype<ostream* >()->AddCast( new E_F1_funcT<ostream*,ostream**>(UnRef<ostream* >)); 
   
//     Add<istream**>("<-","(", new OneUnaryOperator<Op1_new_pstring<istream*,ifstream> >);
     Add<ostream**>("<-","(", new OneUnaryOperator<Op1_new_pstring<ostream*,ofstream> >);
     
     Add<KN<double> *>("min",".",new OneOperator1<double,KN<double> *>(get_min));
     Add<KN<double> *>("max",".",new OneOperator1<double,KN<double> *>(get_max));
    // Polymorphic * precis =new Polymorphic();
    //  Add<ostream*>("precision",".",precis);
     Add<ostream**>("precision",".",new OneOperator1<ostream_precis,ostream**>(ostream_precision));
     Add<ostream*>("precision",".",new OneOperator1<ostream_precis,ostream*>(ostream_precision));
     Add<ostream_precis>("(","",new OneOperator1<long,ostream_precis>(get_precis),
                                new OneOperator2<long,ostream_precis,long>(set_precis));
                                
      
     TheOperators->Add("[]",new OneOperator_array );
     TheOperators->Add("[border]",new OneOperator_border );
     
      
     Global.Add("cos","(",new OneOperator1<double>(cos));
//     Global.Add("square","(",new OneOperator1_<double>(Square));
     Global.Add("square","(",new OneOperator1<double,double,E_F_F0<double,const double &> >(Square));
 //    Global.Add("square","(",new OneUnaryOperator<bUnary_Op<F_1<double,const double&,Square,double> >,
 //      bUnary_Op<F_1<double,const double&,Square,double> >
 //     >);
//   bUnary_Op
     Global.Add("sin","(",new OneOperator1<double>(sin));
     Global.Add("tan","(",new OneOperator1<double>(tan));
     Global.Add("atan","(",new OneOperator1<double>(atan));
     Global.Add("sinh","(",new OneOperator1<double>(sinh));
     Global.Add("cosh","(",new OneOperator1<double>(cosh));
     Global.Add("tanh","(",new OneOperator1<double>(tanh));
     Global.Add("atanh","(",new OneOperator1<double>(atanh));
     Global.Add("asin","(",new OneOperator1<double>(asin));
     Global.Add("acos","(",new OneOperator1<double>(acos));
     Global.Add("asinh","(",new OneOperator1<double>(asinh));
     Global.Add("acosh","(",new OneOperator1<double>(acosh));
     Global.Add("exp","(",new OneOperator1<double>(exp));
     Global.Add("log","(",new OneOperator1<double>(log));
     Global.Add("log10","(",new OneOperator1<double>(log10));
     Global.Add("pow","(",new OneOperator2<double,double>(pow));
//     Global.Add("pow","(",new OneOperator2<double,double,long>(pow));
     Global.Add("max","(",new OneOperator2_<double,double>(Max<double> ));
     Global.Add("min","(",new OneOperator2_<double,double>(Min<double> ));
     Global.Add("max","(",new OneOperator2_<long,long>(Max));
     Global.Add("min","(",new OneOperator2_<long,long>(Min));
     Global.Add("atan2","(",new OneOperator2<double>(atan2));
     Global.Add("atan","(",new OneOperator2<double>(atan2));
     Global.Add("sqrt","(",new OneOperator1<double>(sqrt));
     Global.Add("abs","(",new OneOperator1<double>(Abs));
     Global.Add("abs","(",new OneOperator1<long>(Abs));
     Global.Add("cos","(",new OneOperator1_<Complex>(cos));
     Global.Add("sin","(",new OneOperator1_<Complex>(sin));
     Global.Add("sinh","(",new OneOperator1_<Complex>(sinh));
     Global.Add("cosh","(",new OneOperator1_<Complex>(cosh));
     Global.Add("log","(",new OneOperator1_<Complex>(log));
     //     Global.Add("log10","(",new OneOperator1_<Complex>(log10));
     //     Global.Add("tan","(",new OneOperator1_<Complex>(tan));
     Global.Add("exp","(",new OneOperator1_<Complex>(exp));
     Global.Add("pow","(",new OneOperator2_<Complex,Complex>(pow));
     Global.Add("sqrt","(",new OneOperator1_<Complex>(sqrt));
     Global.Add("conj","(",new OneOperator1_<Complex>(conj));
     Global.Add("imag","(",new OneOperator1_<double,Complex>(imag));
    // Global.Add("real","(",new OneOperator1_<double,Complex>(real));
     Add<Complex>(atype<double>()->name(),".",new OneOperator1_<double,Complex>(real));
     Add<Complex*>(atype<double>()->name(),".",new OneOperator1_<double,Complex*>(preal));
    
     Global.Add("abs","(",new OneOperator1_<double,Complex>(abs));

     Global.Add("arg","(",new OneOperator1_<double,Complex>(arg));
     Global.Add("norm","(",new OneOperator1_<double,Complex>(norm));
     Global.Add("exit","(",new OneOperator1<long>(Exit));     
     Global.Add("assert","(",new OneOperator1<bool>(Assert));     
     
     Global.Add("clock","(",new OneOperator0<double>(CPUtime));
     Global.Add("dumptable","(",new OneOperator1<ostream*,ostream*>(dumptable));
     Global.Add("exec","(",new OneOperator1<long,string* >(exec));
    
     Global.Add("polar","(",new OneOperator2_<Complex,double,double>(polar));
  
     
   //  NEW_TYPE( mapSd  );
  //   NEW_TYPE_I( mapSd  );
     
          
     tables_of_identifier.push_back(&Global);

}

void basicAC_F0::SetNameParam(int n,name_and_type *l , Expression * e) const 
{
 int k=0;
 if ( !n && !named_parameter)  return;
 
  for (int i=0;i<n;i++)
  {
     C_F0  ce=find(l[i].name) ;
     if (ce.LeftValue()==0)
       e[i]=0;
     else  {
       e[i]= map_type[l[i].type->name()]->CastTo(ce);
       k++;
       }
  }
  
 if (!named_parameter) return;
  
  if (k!=  named_parameter->size()) 
   {
      cout << " Sorry some name parameter are not used!  found" <<  k << " == " << named_parameter->size() <<endl;
      for(const_iterator ii=named_parameter->begin(); ii != named_parameter->end();ii++)
       {
        for (int i=0;i<n;i++)
          if (!strcmp(l[i].name,ii->first)) 
            goto L1;
         cout << "\t the parameter is '" << ii->first << "' is unused " << endl;
        L1:;
       }
    if ( n) {
    cerr << " The named parameter can be " << endl;
    for (int i=0;i<n;i++)
       cerr << "\t" << l[i].name << " =  <" << l[i].type->name() << ">\n";
    }
    CompileError("Unused named parameter");  
   }
}

static addingInitFunct TheaddingInitFunct(-10000,Init_map_type); 

void ShowDebugStack()
 {
   if (mpisize)
   cerr << "  current line = " << TheCurrentLine  
        << " mpirank " << mpirank << " / " << mpisize <<endl; 
   else
   cerr << "  current line = " << TheCurrentLine  << endl; 
   while ( debugstack.size() )
     {
        
        cerr << " call " << debugstack.front().first->name<< "  at  line " 
             <<debugstack.front().second << endl; 
        debugstack.pop();     
     }
 }
 

  int  E_F0::Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const 
     {
      int rr = find(m);
      if (rr) return rr;
      if( (verbosity / 100)% 10 == 1) 
      	 cout << "\n new expression : " << n  << " mi=" << MeshIndependent()<< " " << typeid(*this).name()  
      	      << " :" << *this << endl;        
       return insert(this,l,m,n);
     }   
