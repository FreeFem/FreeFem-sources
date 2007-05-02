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

#include "array_init.hpp"

extern Map_type_of_map map_type_of_map ; //  to store te type 
extern Map_type_of_map map_pair_of_type ; //  to store te type 

extern basicForEachType *  typevarreal,  * typevarcomplex;  //  type of real and complex variable

extern int TheCurrentLine; // unset: by default
extern long mpisize,mpirank;
// FH  for g++ 3.4  the prototypage  have change
double  VersionNumber(); 
double Imag(const  complex<double> & z){ return imag(z);}
double Real(const  complex<double> & z){ return real(z);}


// FH

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

 
template<class A,class B>  A Build(B b) {  return A(b);}
  
  


long Exit(long i) {throw(ErrorExec("Exit",i));return 0;}
bool Assert(bool b) {if (!b) throw(ErrorExec("exec assert",1));return true;}
  
inline void MyAssert(int i,char * ex,char * file,long line)
{if (i) {
    cout << "CompileError assertion :  " << ex << " in file " << file << "  line = " << line << endl; 
     CompileError();}
 }


  
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



template<class RR> RR LIncremantation(RR* a){ return ++(*a);}
template<class RR> RR RIncremantation(RR* a){ return (*a)++;}
template<class RR> RR LDecremantation(RR* a){ return --(*a);}
template<class RR> RR RDecremantation(RR* a){ return (*a)--;}

template<class RR,class B>
 RR * New_form_string(string * s) {B * r=  new B(s);delete *s;return r;}
 
 


inline   string ** get_elements( MyMap<String,String> *  const  &  a,string*  const   & b)
 { String* Sret=  &((*a)[*b]); // correction FH feb 2004
    delete b;
    return Sret->getap();}

template<class RR> RR Abs(RR a) { return a<0?-a:a;}

template<class R,class A,class B>
R *MakePtrWithDel( A  const & a)
{ R *r= new B(a->c_str());
  delete a;
  return r;}

template<class R,class RR> 
struct Op1_new_pstring: public unary_function<string*,R> { 
  static R f(string * const & a)  {R r =  new RR(a->c_str());
    // delete a;  (stack ptr) FH mars 2006
    return r;} }; 

template<class R,class RR> 
struct Op2_set_pstring: public binary_function<R,string*,R> { 
  static R  f(R const & p,string * const & a)  {*p =  new RR(a->c_str());
   if ( !*p || !**p) { 
       cerr << " Error openning file " << *a << endl; 
       ExecError("Error openning file");}
  //  delete a; modif mars 2006 FH
   return p;} }; 

template<class R,class RR> 
struct Op2_set_pstringiomode: public ternary_function<R,string*,ios::openmode,R> { 
  static R  f(R const & p,string * const & a,const ios::openmode & mode) 
   {*p =  new RR(a->c_str(),mode);
     // delete a;   modif mars 2006 FH
    return p;} }; 

AnyType FWhile(Stack s ,Expression test,Expression ins)
{ 
  bool sptrclean=true;
  AnyType a;
  StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);
  while ( GetAny<bool>((*test)(s)))
     try  { 
        a=(*ins)(s);
        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
       }
     catch ( E_exception & e) { 
        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
       if (e.code == E_exception::e_break) break;
       else if  (e.code == E_exception::e_continue) continue;
       }
  return a;
}
    
AnyType FFor(Stack s ,Expression i0,Expression i1,Expression i2,Expression ins)
{ 
  bool sptrclean=true;
  AnyType a;
     StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);
  for ( (*i0)(s);GetAny<bool>((*i1)(s));(*i2)(s))
   {
     try  {
        a=(*ins)(s);
        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
       }
     catch ( E_exception & e) { 
        if (verbosity>50)
          cerr << "FFor " << e.what() << e.code << endl; 
        if(sptrclean) sptrclean=sptr->clean(); // modif FH mars 2006  clean Ptr
       if (e.code == E_exception::e_break) break;
       else if  (e.code == E_exception::e_continue) continue;
       }
   }
  return a;
}

AnyType TTry(Stack s ,Expression ins,Expression ccatch,Expression fin,Expression notused)
{ 
  assert(notused == 0);
  AnyType a;
     try  {a=(*ins)(s);}
     catch ( E_exception & e) { 
        throw e; 
       }
     catch(...) {
        if(verbosity> 2) cerr << "Try:: catch (...) exception " << endl;
        a=(*ccatch)(s);         
     }
   
   a=(*fin)(s);
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
    //  delete s;    modif mars 2006 FH
      return r;}
      


 
 class ostream_precis { public:
 ostream_precis(ostream * ff) :f(ff) {}
  ostream * f;
   operator long () const {return f->precision();}
 };
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
 
typedef ios_base& ( * ostream_manipulateur )(ios_base&);

ios_base&  default1(ios_base& f)
{
    f.flags( (ios_base::fmtflags) 0 ) ; // (/*ios_base::scientific | */ios_base::fixed) );
    return f;
}

template< ostream_manipulateur pf> 
inline ostream **set_os(ostream **f)
{ 
    **f << pf  ; return f;
}
template< ostream_manipulateur pf> 
inline ostream *set_os1(ostream *f)
{ 
    *f << pf  ; return f;
}


template<class R>
class  OneOperator_0 : public OneOperator {
  class E_F0_F :public  E_F0mps { public:
    typedef  R (*func)( ) ; 
    func f;
    E_F0_F(func ff)  : f(ff) {}
    AnyType operator()(Stack )  const {return SetAny<R>( f()) ;}  
    operator aType () const { return atype<R>();} 

  };

  typedef  R (*func)() ; 
  func  f;
public: 
  E_F0 * code(const basicAC_F0 & ) const 
  { return  new E_F0_F(f);} 
  OneOperator_0(func  ff): OneOperator(map_type[typeid(R).name()]),f(ff){}
};


void init_by_array(unsigned long init_key[], int key_length);
long genrand_int32(void);
void init_genrand(unsigned long);
long genrandint (long  s) { init_genrand( (unsigned long ) s); return 0;}
long genrandint32 () {return (long)  genrand_int32();}

template<class A,class B,bool RO=true>
struct  MIMul {
  static bool MeshIndependent(Expression a,Expression b)   
   { 
    bool mia= a->MeshIndependent() ;
    bool mib= b->MeshIndependent();
    if ( mia && mib) return true;
    else
      {
        if (mia && a->EvaluableWithOutStack() )
          {
            A va = GetAny<A>((*a)(0));
           // cout << " va = " << va << endl;
            if ( va == A() )
             { 
             //  cout << " va = " << va << endl; 
               return true;
             }
          }
        if (mib && b->EvaluableWithOutStack() )
          {
            B vb = GetAny<B>((*b)(0));
            // cout << " vb = " << vb << endl;
            if ( vb == B() ) 
            { //cout << " vb = " << vb << endl;
             return true; }
           }
         return false;
        }
      
    }
  static bool ReadOnly() { return RO;}
    
};
// add frev 2007
class TransE_Array:  public E_F0 {  public:  
     const E_Array * v;
    int size() const {return v->size();}
    size_t nbitem() const {return v->size();}
    bool MeshIndependent(){return v->MeshIndependent();}
    TransE_Array(const E_Array * e): v(e) {ffassert(e);}
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
};


// add frev 2007 
class opTrans : public OneOperator{
public:
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    opTrans():   OneOperator(atype<TransE_Array>(),atype<E_Array>()  ) {}
    E_F0 * code(const basicAC_F0 & args) const {
	  return new TransE_Array(dynamic_cast<const E_Array*>((Expression) args[0])); } 
};

/*
class opTTrans : public OneOperator{
public:
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    opTTrans():   OneOperator(atype<E_Array>(),atype<TransE_Array>()  ) {}
    E_F0 * code(const basicAC_F0 & args) const {
	return dynamic_cast<const TransE_Array*>((Expression) args[0])->v; } 
};
*/ 

class opDot : public OneOperator{
public:
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    bool MeshIndependent() const { return false;}
    
    opDot(aType A, aType B): OneOperator(atype<C_F0>(),A,B) {}
    opDot(): OneOperator(atype<C_F0>(),atype<TransE_Array >(),atype<E_Array>()  ) {}
    
    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const;       
};

class opFormal : public OneOperator{
public:
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    bool MeshIndependent() const { return false;}
    C_F0  (*thecode2)(const basicAC_F0 &args);     
    opFormal(aType A,C_F0  (c2)(const basicAC_F0 &args) ): OneOperator(atype<C_F0>(),A),thecode2(c2) {}       
    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const { return (*thecode2)(args);}    
};
// fin frev 2007
// avril 2007

C_F0  formalMatTrace(const basicAC_F0 &args)       
{
    bool ta =args[0].left()==atype<TransE_Array>();
    const TransE_Array * tea=0;
    const E_Array * ea=0;
	if( ta)  tea = dynamic_cast<const TransE_Array*>((Expression) args[0]);
    else ea = dynamic_cast<const E_Array*>((Expression) args[0]);
    assert( ea || tea );
    const E_Array & a=  ta ? *tea->v : *ea;
    int ma =1;
    int na=a.size();
    if(na <1 ) CompileError(" trace  [ ...]  ");
    bool maa= a[0].left()==atype<E_Array>();
    if(maa) {
	ma= a[0].LeftValue()->nbitem();
	for (int i=1;i<na;i++)
	    if( ma != a[i].LeftValue()->nbitem()) 
		CompileError(" first matrix with variable number of columm");
        
    }

    int na1=na,ma1=ma;
    if(ta) RNM::Exchange(na1,ma1);
    if(na1 != ma1) CompileError(" trace:  no square matrix ");
    KNM<CC_F0> A(na1,ma1);
    
    if(maa)
	for (int i=0;i<na;++i)
	{
	    const E_Array * li=  dynamic_cast<const E_Array *>(a[i].LeftValue());
	    ffassert(li);
	    for (int j=0; j<ma;++j)
		if(!ta)  A(i,j) = (*li)[j];
		else     A(j,i) = (*li)[j];
	} 
	    else
		for (int i=0;i<na;++i)
		    if(!ta)  A(i,0) = a[i];
		    else     A(0,i) = a[i];
    
    
    CC_F0 s;
    s= A(0,0);
    for (int i=0;i<na1;++i)
	s = C_F0(TheOperators,"+",s,A(i,i));
    return  s;
     
}

C_F0  formalMatDet(const basicAC_F0 &args)       
{
    bool ta =args[0].left()==atype<TransE_Array>();
    const TransE_Array * tea=0;
    const E_Array * ea=0;
    if( ta)  tea = dynamic_cast<const TransE_Array*>((Expression) args[0]);
    else ea = dynamic_cast<const E_Array*>((Expression) args[0]);
    assert( ea || tea );
    const E_Array & a=  ta ? *tea->v : *ea;
    int ma =1;
    int na=a.size();
    if(na <1 ) CompileError(" trace  [ ...]  ");
    bool maa= a[0].left()==atype<E_Array>();
    if(maa) {
	ma= a[0].LeftValue()->nbitem();
	for (int i=1;i<na;i++)
	    if( ma != a[i].LeftValue()->nbitem()) 
		CompileError("  matrix with variable number of columm");
        
    }
    
    int na1=na,ma1=ma;
    if(ta) RNM::Exchange(na1,ma1);
    if(na1 != ma1) CompileError(" trace:  no square matrix ");
    KNM<CC_F0> A(na1,ma1);
    
    if(maa)
	for (int i=0;i<na;++i)
	{
	    const E_Array * li=  dynamic_cast<const E_Array *>(a[i].LeftValue());
	    ffassert(li);
	    for (int j=0; j<ma;++j)
		if(!ta)  A(i,j) = (*li)[j];
		else     A(j,i) = (*li)[j];
	} 
	    else
		for (int i=0;i<na;++i)
		    if(!ta)  A(i,0) = a[i];
		    else     A(0,i) = a[i];
    
    
    if(na1==1)
      return  A(0,0);
    else if( na1==2 )
    {
	C_F0 s1(TheOperators,"*",A(0,0),A(1,1));
	C_F0 s2(TheOperators,"*",A(0,1),A(1,0));
	return C_F0(TheOperators,"-",s1,s2);
    }
    else
    {
	CompileError("FH: sorry only det of 1x1 and 2x2 matrix ");
    }
    return  C_F0(); 
    
}




// fiun avril 2007
void Init_map_type()
{
   TheOperators=new Polymorphic(), 
   TheRightOperators=new Polymorphic();
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
    Dcl_Type< istream_good > ();
    Dcl_Type< NothingType > ();
    
    Dcl_Type<Polymorphic*>();
    
//    Dcl_Type<C_F0>();
    map_type[typeid(C_F0).name()] =  new TypeLineFunction; 
    Dcl_Type<E_Array>();
    Dcl_Type<TransE_Array >();// add
    Dcl_Type<const E_Border *>();
    Dcl_Type<const E_BorderN *>();

    
    
    Dcl_Type<SubArray>();
    Dcl_Type<pair<long,long> >();
    
    initArrayDCLlong();
    initArrayDCLdouble();
    initArrayDCLComplex();
        
    Dcl_Type<ios::openmode>();
    
//  les types des variables 
    
  zzzfff->Add("real",typevarreal=atype<double*>());
  zzzfff->Add("int",atype<long*>());
  zzzfff->Add("complex",typevarcomplex=atype<Complex*>());
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
       new OneBinaryOperator_st<Op2_padd<string,string*,string*> >  // a changer to do FH string * mars 2006    
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
       new OneBinaryOperator<Op2_mul<double,double,double>,MIMul<double,double> >,
       new OneBinaryOperator<Op2_mul<double,double,long>, MIMul<double,long> >,
       new OneBinaryOperator<Op2_mul<double,long,double>,MIMul<long,double> >,
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
       new OneBinaryOperator<Op2_plt<string*,string*> >  //  FH string * mars 2006 
     );
     TheOperators->Add("<=",
       new OneBinaryOperator<Op2_le<long,long> >,
       new OneBinaryOperator<Op2_le<double,double> >,
       new OneBinaryOperator<Op2_ple<string*,string*> >  //  FH string * mars 2006 
     );
     TheOperators->Add(">",
       new OneBinaryOperator<Op2_gt<long,long> >,
       new OneBinaryOperator<Op2_gt<double,double> >,
       new OneBinaryOperator<Op2_pgt<string*,string*> >  //  string * mars 2006 
     );
     TheOperators->Add(">=",
       new OneBinaryOperator<Op2_ge<long,long> >,
       new OneBinaryOperator<Op2_ge<double,double> >,
       new OneBinaryOperator<Op2_pge<string*,string*> >  //  FH string * mars 2006 
     );
     TheOperators->Add("==",
       new OneBinaryOperator<Op2_eq<long,long> >,
       new OneBinaryOperator<Op2_eq<double,double> >,
       new OneBinaryOperator<Op2_eq<Complex,Complex> >,
       new OneBinaryOperator<Op2_peq<string*,string*> >  //   FH string * mars 2006 
     );

     TheOperators->Add("!=",
       new OneBinaryOperator<Op2_ne<long,long> >,
       new OneBinaryOperator<Op2_ne<double,double> >,
       new OneBinaryOperator<Op2_ne<Complex,Complex> >,
       new OneBinaryOperator<Op2_pne<string*,string*> >  //  FH string * mars 2006 
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
       new OneBinaryOperator<set_eq<bool> ,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<set_eq<long> ,OneBinaryOperatorMIWO>,
       new OneBinaryOperator<set_eq<double> ,OneBinaryOperatorMIWO>,
       new OneBinaryOperator<set_eq<Complex> ,OneBinaryOperatorMIWO>,
       new OneBinaryOperator<set_peq<string> ,OneBinaryOperatorMIWO>  // FH string * mars 2006 
       ); 

     TheOperators->Add("?:",
       new Operator_Aritm_If<long >,
       new Operator_Aritm_If<double >,
       new Operator_Aritm_If<Complex >,
       new Operator_Aritm_If<string* >  // (OK???)  to do FH string * mars 2006 
       ); 
       
/*
     ArrayOperator<double>();
     ArrayOperator<Complex>();
     ArrayOperator<long>();
*/
//      initArrayOperators()   ;  
     initArrayOperatorlong();
     initArrayOperatordouble();
     initArrayOperatorComplex();
     initStringOperator();


     TheOperators->Add("+=",
       new OneBinaryOperator<set_eq_add<long>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<set_eq_add<double>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<set_eq_add<Complex>,OneBinaryOperatorMIWO >
      );


     TheOperators->Add("-=",
       new OneBinaryOperator<set_eq_sub<long>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<set_eq_sub<double>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<set_eq_sub<Complex>,OneBinaryOperatorMIWO >
      );



     TheOperators->Add("*=",
       new OneBinaryOperator<set_eq_mul<long> ,OneBinaryOperatorMIWO>,
       new OneBinaryOperator<set_eq_mul<double>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<set_eq_mul<Complex>,OneBinaryOperatorMIWO >     
      );


     TheOperators->Add("/=",
       new OneBinaryOperator<set_eq_div<long>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<set_eq_div<double>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<set_eq_div<Complex>,OneBinaryOperatorMIWO >     
     );

     TheOperators->Add("+",
    //   new OneBinaryOperator<Op2_addp<const E_BorderN *,const E_BorderN *,const E_BorderN * > >,  
       new AddBorderOperator
       );

      // add frev 2007
      TheOperators->Add("\'", new opTrans); 
      
     // TheOperators->Add("\'", new opTTrans); 
      TheOperators->Add("*",new opDot(atype<TransE_Array >(),atype<E_Array>() )   );  // a faire mais dur 
      TheOperators->Add("*",new opDot(atype<E_Array >(),atype<E_Array>() )   );  // a faire mais dur 
      TheOperators->Add("*",new opDot(atype<E_Array >(),atype<TransE_Array>() )   );  // a faire mais dur 
      TheOperators->Add("*",new opDot(atype<TransE_Array >(),atype<TransE_Array>() )   );  // a faire mais dur 
     // car le type de retour depent des objets du tableau 
     // il faut refechir  .....  FH 
     // il faut definir le type d'un tableau bof, bof (atype<C_F0>())
     TheOperators->Add(">>",
       new OneBinaryOperator<Op_Read<bool>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<Op_Read<long>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<Op_Read<double>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<Op_Read<Complex>,OneBinaryOperatorMIWO >,
       new OneBinaryOperator<Op_ReadP<string>,OneBinaryOperatorMIWO >
       );
     
     TheOperators->Add("<<",
       new OneBinaryOperator<Print<bool> >,
       new OneBinaryOperator<Print<long> >,
       new OneBinaryOperator<Print<double> >,
       new OneBinaryOperator<Print<Complex> >,
       new OneBinaryOperator<PrintP<string*> >  //  FH string * mars 2006 
       );

     
     TheRightOperators->Add("++",       
       new OneOperator1<long,long*, E_F_F0<long,long*,false> >(&RIncremantation<long>));
     TheRightOperators->Add("--",       
       new OneOperator1<long,long*, E_F_F0<long,long*,false> >(&RDecremantation<long>));
     TheOperators->Add("++",       
       new OneOperator1<long,long*, E_F_F0<long,long*,false> >(&LIncremantation<long>));
     TheOperators->Add("--",       
       new OneOperator1<long,long*, E_F_F0<long,long*,false> >(&LDecremantation<long>));
//   init        
     TheOperators->Add("<-", 
       new OneOperator2<string**,string**,string*>(&set_copyp_new<string>),  //  FH string * mars 2006 
       new OneOperator2_<double*,double*,double>(&set_copyp),  // 
       new OneOperator2_<long*,long*,long>(&set_copyp),
       new OneOperator2_<bool*,bool*,bool>(&set_copyp), //  mars 2006
       new OneOperator2_<Complex*,Complex*,Complex>(&set_copy),
       new OneOperator2_<istream**,istream**,istream*>(&set_copy),
       new OneOperator2_<ostream**,ostream**,ostream*>(&set_copy),
//       new OneUnaryOperator<Op1_new_pstring<istream*,ifstream> >,
//       new OneUnaryOperator<Op1_new_pstring<ostream*,ofstream> >,
       new OneBinaryOperator<Op2_set_pstring<istream**,ifstream> >,  //  FH string * mars 2006 
       new OneBinaryOperator<Op2_set_pstring<ostream**,ofstream> >,  //  FH string * mars 2006 
       new OneTernaryOperator3<Op2_set_pstringiomode<ostream**,ofstream> >      //  FH string * mars 2006   
       );  
       
     atype<istream* >()->AddCast( new E_F1_funcT<istream*,istream**>(UnRef<istream* >)); 
     atype<ostream* >()->AddCast( new E_F1_funcT<ostream*,ostream**>(UnRef<ostream* >)); 
   
//     Add<istream**>("<-","(", new OneUnaryOperator<Op1_new_pstring<istream*,ifstream> >);
     Add<ostream**>("<-","(", new OneUnaryOperator<Op1_new_pstring<ostream*,ofstream> >);  //  FH string * mars 2006 
     
    // Polymorphic * precis =new Polymorphic();
    //  Add<ostream*>("precision",".",precis);
     Add<ostream**>("precision",".",new OneOperator1<ostream_precis,ostream**>(ostream_precision));
     Add<ostream*>("precision",".",new OneOperator1<ostream_precis,ostream*>(ostream_precision));
     
     
     Add<ostream_precis>("(","",new OneOperator1<long,ostream_precis>(get_precis),
                                new OneOperator2<long,ostream_precis,long>(set_precis));
//  add v 1.41   
     Add<istream**>("good",".",new OneOperator1<istream_good,istream**>(to_istream_good));
     Add<istream*>("good",".",new OneOperator1<istream_good,istream*>(to_istream_good));
     Add<istream*>("good",".",new OneOperator1<istream_good,istream*>(to_istream_good));    
     Add<istream_good>("(","",new OneOperator1<long,istream_good>(get_good));

     Add<istream**>("eof",".",new OneOperator1<bool,istream**>(get_eof));
// add v 2.8 
     Add<ostream**>("scientific",".",new OneOperator1<ostream**,ostream**>(set_os<scientific>));
     Add<ostream**>("fixed",".",new OneOperator1<ostream**,ostream**>(set_os<fixed>));
     Add<ostream**>("showbase",".",new OneOperator1<ostream**,ostream**>(set_os<showbase>));
     Add<ostream**>("noshowbase",".",new OneOperator1<ostream**,ostream**>(set_os<noshowbase>));
     Add<ostream**>("showpos",".",new OneOperator1<ostream**,ostream**>(set_os<showpos>));
     Add<ostream**>("noshowpos",".",new OneOperator1<ostream**,ostream**>(set_os<noshowpos>));
     Add<ostream**>("default",".",new OneOperator1<ostream**,ostream**>(set_os<default1>));
     
     Add<ostream*>("scientific",".",new OneOperator1<ostream*,ostream*>(set_os1<scientific>));
     Add<ostream*>("fixed",".",new OneOperator1<ostream*,ostream*>(set_os1<fixed>));
     Add<ostream*>("showbase",".",new OneOperator1<ostream*,ostream*>(set_os1<showbase>));
     Add<ostream*>("noshowbase",".",new OneOperator1<ostream*,ostream*>(set_os1<noshowbase>));
     Add<ostream*>("showpos",".",new OneOperator1<ostream*,ostream*>(set_os1<showpos>));
     Add<ostream*>("noshowpos",".",new OneOperator1<ostream*,ostream*>(set_os1<noshowpos>));
     Add<ostream*>("default",".",new OneOperator1<ostream*,ostream*>(set_os1<default1>));
     
// add 2.16
     TheOperators->Add("trace",new opFormal(atype<E_Array>(),formalMatTrace ));
     TheOperators->Add("det",new opFormal(atype<E_Array>(),formalMatDet ));
// end add     
                                
      
     TheOperators->Add("[]",new OneOperator_array );
     TheOperators->Add("[border]",new OneOperator_border );
     
      
     Global.Add("cos","(",new OneOperator1<double>(cos));
//     Global.Add("square","(",new OneOperator1_<double>(Square));
     Global.Add("square","(",new OneOperator1<double,double,E_F_F0<double,const double &> >(Square));

     Global.Add("floor","(",new OneOperator1<double>(floor)); // add march 2006
     Global.Add("ceil","(",new OneOperator1<double>(ceil));  // add march 2006
     Global.Add("rint","(",new OneOperator1<double>(rint));  // add june 2006
     
     Global.Add("sin","(",new OneOperator1<double>(sin));
     Global.Add("tan","(",new OneOperator1<double>(tan));
     Global.Add("atan","(",new OneOperator1<double>(atan));
     Global.Add("sinh","(",new OneOperator1<double>(sinh));
     Global.Add("cosh","(",new OneOperator1<double>(cosh));
     Global.Add("tanh","(",new OneOperator1<double>(tanh));
#ifdef HAVE_ATANH
     Global.Add("atanh","(",new OneOperator1<double>(atanh));
#endif
     Global.Add("asin","(",new OneOperator1<double>(asin));
     Global.Add("acos","(",new OneOperator1<double>(acos));
#ifdef HAVE_ASINH
     Global.Add("asinh","(",new OneOperator1<double>(asinh));
#endif
#ifdef HAVE_ACOSH
     Global.Add("acosh","(",new OneOperator1<double>(acosh));
#endif
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
     TheOperators->Add("\'",new OneOperator1_<Complex>(conj));       
     
     
     Global.Add("imag","(",new OneOperator1_<double,Complex>(Imag));
     //  Big probleme  real is a type
     Add<double>("<--","(",new OneOperator1_<double,Complex>(Real));
    // Global.Add("real","(",new OneOperator1_<double,Complex>(Real));
    // Add<double>(typevarreal->right()->name(),".",new OneOperator1_<double,Complex>(Real));
    // Global.Add(typevarreal->right()->name(),".",new OneOperator1_<double,Complex>(Real));
    // Add<double*>(typevarreal->left()->name(),".",new OneOperator1_<double,Complex*>(preal));
    
     Global.Add("abs","(",new OneOperator1_<double,Complex>(abs));

     Global.Add("arg","(",new OneOperator1_<double,Complex>(arg));
     Global.Add("norm","(",new OneOperator1_<double,Complex>(norm));
     Global.Add("exit","(",new OneOperator1<long>(Exit));     
     Global.Add("assert","(",new OneOperator1<bool>(Assert));     
     
     Global.Add("clock","(",new OneOperator0<double>(CPUtime));
     Global.Add("dumptable","(",new OneOperator1<ostream*,ostream*>(dumptable));
     Global.Add("exec","(",new OneOperator1<long,string* >(exec));  //FH string * mars 2006 
    
     Global.Add("polar","(",new OneOperator2_<Complex,double,double>(polar));
 // rand generator ---
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;                                                            
  init_by_array(init, length);
  extern long genrand_int31(void);   
  extern double genrand_real1(void);
  extern double genrand_real2(void);
  extern double genrand_real3(void);
  extern double  genrand_res53(void) ;
  
  Global.Add("randint32","(",new OneOperator_0<long>(genrandint32));
  Global.Add("randint31","(",new OneOperator_0<long>(genrand_int31));
  Global.Add("randreal1","(",new OneOperator_0<double>(genrand_real1));
  Global.Add("randreal2","(",new OneOperator_0<double>(genrand_real2));
  Global.Add("randreal3","(",new OneOperator_0<double>(genrand_real3));
  Global.Add("randres53","(",new OneOperator_0<double>(genrand_res53));
  Global.Add("randinit","(",new OneOperator1<long>(genrandint));
  
       

typedef MyMap<String,String> MyMapSS;
     map_type[typeid(MyMapSS*).name()] = new ForEachType<MyMapSS*>(Initialize<MyMapSS >,Delete<MyMapSS >) ;         
     map_type_of_map[make_pair(atype<string*>(),atype<string*>())]=atype<MyMapSS*>();      
     atype<MyMapSS*>()->Add("[","",new OneOperator2_<string**,MyMapSS*,string*>(get_elements));
     
          
     tables_of_identifier.push_back(&Global);

}
int ShowAlloc(const char *s,size_t & lg); 




 void ClearMem()
 {
     size_t lg;
     ShowAlloc("ClearMem: begin" , lg);
     delete pZero;
     delete pOne;
     delete pminusOne;
      
     tables_of_identifier.clear();
     for (map<const string,basicForEachType *>::iterator i=map_type.begin();i!=map_type.end();++i)
        delete i->second;
        
     map_type.clear();
     map_type_of_map.clear();
     map_pair_of_type.clear();
     Global.clear();
     if(TheOperators) 
       TheOperators->clear();
     if(TheRightOperators)
       TheRightOperators->clear();
        
     CodeAlloc::clear();
     ShowAlloc("ClearMem: end" , lg); 

 } 
static addingInitFunct TheaddingInitFunct(-10000,Init_map_type); 

C_F0  opDot::code2(const basicAC_F0 &args) const      
{
    bool ta =args[0].left()==atype<TransE_Array>();
    bool tb = args[1].left()==atype<TransE_Array>();
    const TransE_Array * tea=0;
    const TransE_Array * teb=0;
    const E_Array * ea=0;
    const E_Array * eb=0;// E_F0
	if( ta)  tea = dynamic_cast<const TransE_Array*>((Expression) args[0]);
    else ea = dynamic_cast<const E_Array*>((Expression) args[0]);
    if( tb)  teb = dynamic_cast<const TransE_Array*>((Expression) args[1]);
    else eb = dynamic_cast<const E_Array*>((Expression) args[1]);
    assert( ea || tea );
    assert( eb || teb );
    const E_Array & a=  ta ? *tea->v : *ea;
    const E_Array & b=  tb ? *teb->v : *eb;
    int ma =1;
    int mb =1;
    int na=a.size();
    int nb=b.size();
    if(na <1 && nb < 1) CompileError(" empty array  [ ...]'*[ ...  ]  ");
    bool mab= b[0].left()==atype<E_Array>();
    bool maa= a[0].left()==atype<E_Array>();
    if(maa) {
	ma= a[0].LeftValue()->nbitem();
	for (int i=1;i<na;i++)
	    if( ma != a[i].LeftValue()->nbitem()) 
		CompileError(" first matrix with variable number of columm");
        
    }
    if(mab) {
	mb= b[1].LeftValue()->nbitem();
	for (int i=1;i<nb;i++)
	    if( mb != b[i].LeftValue()->nbitem()) 
		CompileError(" second matrix with variable number of columm");
    }
    int na1=na,ma1=ma,nb1=nb,mb1=mb;
    if(ta) RNM::Exchange(na1,ma1);
    if(tb) RNM::Exchange(nb1,mb1);
    
    KNM<CC_F0> A(na1,ma1), B(nb1,mb1);
    if ( A.M() != B.N())
    {
	cout << "   formal prod array or matrix : [ .. ] * [ .. ]   " << endl;
	cout << " first  array :  matrix " << maa << " trans " << ta << " " << na << "x" << ma <<endl;
	cout << " second array :  matrix " << mab << " trans " << tb << " " << nb << "x" << mb <<endl;	
	CompileError(" no same size  [ ...]'*[ ...  ] sorry ");
    }
    
    if(maa)
	for (int i=0;i<na;++i)
	{
	    const E_Array * li=  dynamic_cast<const E_Array *>(a[i].LeftValue());
	    ffassert(li);
	    for (int j=0; j<ma;++j)
		if(!ta)  A(i,j) = (*li)[j];
		else     A(j,i) = (*li)[j];
	} 
    else
	for (int i=0;i<na;++i)
	    if(!ta)  A(i,0) = a[i];
	    else     A(0,i) = a[i];
	 
    if(mab)
	for (int i=0;i<nb;++i)
	{
	    const E_Array * li=  dynamic_cast<const E_Array *>(b[i].LeftValue());
	    ffassert(li);
	    for (int j=0; j<mb;++j)
		if(!tb)  B(i,j) = (*li)[j];
		else     B(j,i) = (*li)[j];
	} 
    else
	for (int i=0;i<nb;++i)
	    if(!tb)  B(i,0) = b[i];
	    else     B(0,i) = b[i];
    
    KNM<CC_F0> C(na1,mb1);
    CC_F0 s,abi;
    for (int i=0;i<na1;++i)
	for (int j=0;j<mb1;++j)
	{
	    s= C_F0(TheOperators,"*",A(i,0),B(0,j));
	    for (int k=1;k<ma1;++k) {
		abi = C_F0(TheOperators,"*",A(i,k),B(k,j));
		s = C_F0(TheOperators,"+",s,abi);}
	    C(i,j)=s;
	};
    if( na1==1 && mb1 ==1)
	return C(0,0);
    else if (( mb1 ==1) || (na1==1))
    {
	AC_F0  v;
	v=C(0,0);
	int i0=na1!=1,j0=mb1!=1, nn= mb1*na1;
	for (int i=1;i<nn;++i)
	    v+=C(i0*i,j0*i);
	C_F0  r(TheOperators,"[]",v);
	if(mb1==1) return r;
	else return C_F0(TheOperators,"\'",r);
    }
    else
    {
	AC_F0  v,cc;
	v=C(0,0);
	for (int i=0;i<na1;++i)	
	{  cc = C(i,0);
	    for (int j=1;j<mb1;++j)
		cc+= C(i,j);
	    C_F0  vi(TheOperators,"[]",cc);
	    if(i==0) v=vi;
	    else v+= vi;
	}
	return C_F0(TheOperators,"[]",v);
    }
/*	  
    if ( !mab && ! maa)
    {
	
	if( na != nb)
	    CompileError(" no same size  [ ...]'*[ ...  ] sorry ");
	
	if( ta && ! tb)
	{
	    s= C_F0(TheOperators,"*",a[0],b[0]);
	    for (int i=1;i<na;++i)
	    {
		abi = C_F0(TheOperators,"*",a[i],b[i]);
		s = C_F0(TheOperators,"+",s,abi);
	    }
	    return s;//Type_Expr(s); //new C_F0(s);   ATTENTION le type est variable ici   FH
	}
	
	if(!ma && mb)
	{  
	}
	
    }*/
    
    cout << "   formal prod array or matrix : [ .. ] * [ .. ]   " << na << "x" << nb << endl;
    cout << "   formal prod array or matrix : [ .. ] * [ .. ]   " <<  endl;
    cout << " first  array :  matrix " << maa << " trans " << ta << " " << na << "x" << ma <<endl;
    cout << " second array :  matrix " << mab << " trans " << tb << " " << nb << "x" << mb <<endl;
    CompileError("  not implemented sorry ..... (FH) to do ???? ");	
    return C_F0();

}

