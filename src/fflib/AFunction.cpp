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

//  for bessel function
// c++11   => __STRICT_ANSI__ => error FH..
#ifdef __STRICT_ANSI__
#undef __STRICT_ANSI__
#endif

// TODO: remove this block as soon as autoconf is removed from FreeFem++
#ifndef CMAKE
#include <config.h>
#endif

#include <cmath>
#include <complex>
//  put here some def dur to c++11
// problem with mixed with using namespace std;
// to correct bug in g++ v 4.8.1 add std
#if defined (_WIN32  ) || (__GNUC__ >=5) || __llvm__
#define NM_STD std::
#else
#define NM_STD
#endif
long isNaN(double x){return NM_STD isnan(x);}
long isInf(double x){return NM_STD isinf(x);}
long isNormal(double x){return std::isnormal(x);}
#ifdef HAVE_JN
double myyn(long n, double x){ return yn((int)n,x);}
double myjn(long n, double x){ return jn((int) n,x);}
#endif
//int  ShowAlloc(const char *s, size_t lg);

// F. Hecht fev. 2015 ...
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
#include "AFunction_ext.hpp"
// Add FH to get memroy used in test .. march 2014
#if __APPLE__
#include <malloc/malloc.h>
#elif HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_TIMES
#include <ctime>
#endif
long storageused()
{
#if HAVE_MSTATS
    struct mstats mem1;
    mem1 = mstats();
    return mem1.bytes_used;
#elif HAVE_MALLINFO
    struct mallinfo mem1;
    mem1=mallinfo();
    return mem1.uordblks;
#else
    return 0;
#endif

}
long storagetotal()
{
#if HAVE_MSTATS
    struct mstats mem1;
    mem1 = mstats();
    return mem1.bytes_total;
#elif HAVE_MALLINFO
    struct mallinfo mem1;
    mem1=mallinfo();
    return mem1.keepcost;
#else
    return 0;
#endif
}
// end add mach 2014 ...
extern Map_type_of_map map_type_of_map ; //  to store te type
extern Map_type_of_map map_pair_of_type ; //  to store te type

extern basicForEachType *  typevarreal,  * typevarcomplex;  //  type of real and complex variable

extern int TheCurrentLine; // unset: by default
extern long mpisize,mpirank;

// rand generator ---
extern long genrand_int31(void);
extern double genrand_real1(void);
extern double genrand_real2(void);
extern double genrand_real3(void);
extern double  genrand_res53(void) ;


// FH  for g++ 3.4  the prototypage  have change
double  VersionNumber();
double Imag(const  complex<double> & z){ return imag(z);}
double Real(const  complex<double> & z){ return real(z);}
const  basicForEachType * basicForEachType::type_C_F0 =0; //  for any type un formal operation .... FH add 09/2012

// FH

template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
template<class T> inline T Min (const T &a,const T & b){return a < b ? a : b;}
template<class T> inline T Abs (const T &a){return a <0 ? -a : a;}
template<class T> inline T Max (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}
template<class T> inline T Min (const T &a,const T & b,const T & c){return Min(Min(a,b),c);}
template<class T> inline T Min (const T &a,const T & b,const T & c,const T & d){return Min(Min(a,b),Min(c,d));}
template<class T> inline T Max (const T &a,const T & b,const T & c,const T & d){return Max(Max(a,b),Max(c,d));}

template<class T> inline T Square (const T &a){return a*a;}

struct SubArray2: public binary_function<long,long,SubArray> {
  static SubArray f(const long & a,const long & b)  {
    return SubArray(b-a+1,a);} };
struct SubArray3: public ternary_function<long,long,long,SubArray> {
  static SubArray f(Stack s,const long & a,const long & b,const long & c)  {
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




long Exit(long i) {throw(ErrorExit("Exit",i));return 0;}
bool Assert(bool b) {if (!b) throw(ErrorExec("exec assert",1));return true;}

inline void MyAssert(int i,char * ex,char * file,long line)
{if (i) {
    cout << "CompileError assertion :  " << ex << " in file " << file << "  line = " << line << endl;
     CompileError();}
 }

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
    E_F0 * code(const basicAC_F0 & a) const {
        if (a.size()==1 && a[0].left()==atype<E_Array>() )
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
 RR * New_form_string(string * s) {B * r=  new B(s);freestring(s);return r;}// correct Mars 2011 remove * if delete

inline   string ** get_elements( MyMap<String,String> *  const  &  a,string*  const   & b)
 { String* Sret=  &((*a)[*b]); // correction FH feb 2004
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
    return r;} };

template<class R,class RR>
struct Op2_set_pstring: public binary_function<R,string*,R> {
  static R  f(R const & p,string * const & a)  {*p =  new RR(a->c_str());
   if ( !*p || !**p) {
       cerr << " Error opening file " << *a << endl;
       ExecError("Error opening file");}
   return p;} };

template<class R,class RR>
struct Op2_set_pstringiomode: public ternary_function<R,string*,ios::openmode,R> {
  static R  f(Stack s,R const & p,string * const & a,const ios::openmode & mode)
   {*p =  new RR(a->c_str(),mode);
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
       else throw e;
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
       else throw e;
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
   if(fin)  a=(*fin)(s);
  return a;
}

AnyType FIf(Stack s ,Expression test,Expression i1,Expression i2,Expression )
 {  AnyType a;
   if (GetAny<bool>((*test)(s)))
     {
       if(i1) a=(*i1)(s);//Add if FH oct 2010
     }
      else if (i2)
    {
     if(i2) a=(*i2)(s); //Add if FH oct 2010
    }

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

// <<dumptable>>
ostream* dumptable(ostream* f)
{

  *f << " the keywords " << endl;
  ShowKeyWord(*f);
  *f << " the types " << endl;
  ShowType(*f);
   ListOfTOfId::const_iterator i=tables_of_identifier.begin();
   for(;i!=tables_of_identifier.end();++i)
    {
      cout << "  --------- table of identifier ---------\n";
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


class OP_setw { public:
    long w;
    OP_setw(long ww) :w(ww) {}
    friend   ostream & operator<<(ostream & f,const OP_setw& op) { return f << setw(op.w);}
};
 OP_setw defOP_setw(long i) {return OP_setw(i);}


 ostream_precis ostream_precision(ostream **f){ return ostream_precis(*f);}
  ostream_precis ostream_precision(ostream *f){ return ostream_precis(f);}
 long get_precis( ostream_precis  pf) { return pf.f->precision();}
 long set_precis( ostream_precis  pf, long  l) { return pf.f->precision(l);}

class ostream_seekp { public:
    ostream_seekp(ostream * ff) :f(ff) {}
    ostream * f;
    operator long () const {return f->tellp();}
};


class istream_seekg { public:
    istream_seekg(istream * ff) :f(ff) {}
    istream * f;
    operator long () const {return f->tellg();}
};

ostream_seekp ff_oseekp(ostream **f){ return ostream_seekp(*f);}
ostream_seekp ff_oseekp(ostream *f){ return ostream_seekp(f);}
istream_seekg ff_iseekg(istream **f){ return istream_seekg(*f);}
istream_seekg ff_iseekg(istream *f){ return istream_seekg(f);}

long ffseekp( ostream_seekp  pf, long  l) {  pf.f->clear();long ll= pf.f->tellp(); return pf.f->seekp(l),ll;}
long fftellp( ostream_seekp  pf) { pf.f->clear(); return pf.f->tellp() ;}
long ffseekg( istream_seekg  pf, long  l) { pf.f->clear(); return pf.f->seekg(l),l;}
long fftellg( istream_seekg  pf) { return pf.f->tellg() ;}

 class istream_good { public:
  istream_good(istream * ff) :f(ff) {}
  istream * f;
  operator bool () const {return f->good();}
 };
 inline istream_good to_istream_good(istream **f){ return istream_good(*f);}
 inline istream_good to_istream_good(istream *f){ return istream_good(f);}

  inline long get_good( istream_good  pf) { return pf.f->good();}
  inline bool get_eof(istream ** p){ return (**p).eof();}
  long filelength(istream ** p) {
      istream *f = *p;
      long where=f->tellg();
      f->seekg (0, f->end);
      long length =f->tellg();
      f->seekg (where);
      return length;
  }
  bool eatspace(istream ** p){
    istream *f = *p;
    int c;
     while ( (c=f->peek())!=EOF)
         if ( !isspace(c)) break;
         else f->get();
     return c != EOF;}

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

inline ostream **set_os_flush(ostream **f)
{
    (**f).flush()  ; return f;
}
inline ostream *set_os_flush(ostream *f)
{
    (*f).flush()  ; return f;
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
            A va = GetAny<A>((*a)(NullStack));
            if ( va == A() )
             {
               return true;
             }
          }
        if (mib && b->EvaluableWithOutStack() )
          {
            B vb = GetAny<B>((*b)(NullStack));
            if ( vb == B() )
            {
             return true; }
           }
         return false;
        }

    }
  static bool ReadOnly() { return RO;}

};


// add frev 2007
class opTrans : public OneOperator{
public:
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    opTrans():   OneOperator(atype<TransE_Array>(),atype<E_Array>()  ) {}
    E_F0 * code(const basicAC_F0 & args) const {
	  return new TransE_Array(dynamic_cast<const E_Array*>((Expression) args[0])); }
};

class opDot : public OneOperator{
public:
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    bool MeshIndependent() const { return false;}

    opDot(aType A, aType B): OneOperator(atype<C_F0>(),A,B) {}
    opDot(): OneOperator(atype<C_F0>(),atype<TransE_Array >(),atype<E_Array>()  ) {}

    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const;
};

class opColumn : public OneOperator{
public:
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    bool MeshIndependent() const { return false;}

    opColumn(aType A, aType B): OneOperator(atype<C_F0>(),A,B) {if( A== basicForEachType::type_C_F0)pref=-100;}
    opColumn(aType A): OneOperator(atype<C_F0>(),ArrayOfaType(A,true)) {pref=-100;}

    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const;
};

class opSum : public OneOperator{
public:
    const char * op;
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    bool MeshIndependent() const { return false;}

    opSum(const char *opp,aType A, aType B): OneOperator(atype<C_F0>(),A,B),op(opp) {}

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
// nov 2007   v[i]
class opVI : public OneOperator{
public:
    AnyType operator()(Stack s)  const {ffassert(0);return 0L;}
    bool MeshIndependent() const { return false;}

    opVI(aType A): OneOperator(atype<C_F0>(),A,atype<long>()) {}

    E_F0 *  code(const basicAC_F0 & ) const {ffassert(0);}
    C_F0  code2(const basicAC_F0 &args) const;
};


C_F0  formalMatCofactor(const basicAC_F0 &args)
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
    if(na <1 ) CompileError(" Cofactor  ([ ...])  ");
    bool maa= a[0].left()==atype<E_Array>();
    if(maa) {
	ma= a[0].LeftValue()->nbitem();
	for (int i=1;i<na;i++)
	    if( ma != (int) a[i].LeftValue()->nbitem())
		CompileError(" a matrix with variable number of columm");

    }

    int na1=na,ma1=ma;
    if(ta) RNM::Exchange(na1,ma1);
    if(na1 != ma1) CompileError(" CoFactor:  no square matrix ");
    if(na1 > 3 || ( na1 <1) ) CompileError(" CoFactor:   square matrix size is more then 3   ");
    KNM<CC_F0> A(na1,na1);
    KNM<CC_F0> C(na1,na1);
    if(maa)
	for (int i=0;i<na;++i)
	{
	    const E_Array * li=  dynamic_cast<const E_Array *>(a[i].LeftValue());
	    ffassert(li);
	    for (int j=0; j<ma;++j)
		if(!ta)  A(i,j) = (*li)[j];
		else     A(j,i) = TryConj((*li)[j]);
	}
    else
        for (int i=0;i<na;++i)
            if(!ta)  A(i,0) = a[i];
            else     A(0,i) = TryConj(a[i]);


    AC_F0  v,cc;
    if(na1==2)
    {
      for(int i=0;i<na1;++i)
        for(int j=0;j<na1;++j)
          if( (i+j) %2 == 0)
            C(i,j) = A(1-i,1-j);
          else
            C(i,j) = C_F0(TheOperators,"-",A(1-i,1-j));
    }
    else if( na1 ==3)
    {
        int i1,i2,j1,j2;
        for(int i=0;i<3;++i)
          for(int j=0;j<3;++j)
          {
              i1 = (i+1)%3;
              i2 = (i+2)%3;
              j1 = (j+1)%3;
              j2 = (j+2)%3;

            C(i,j) = A(i1,j1)*A(i2,j2)-A(i1,j2)*A(i2,j1);
          }
    }
    v=C(0,0);

    for (int i=0;i<na1;++i)
    {  cc = C(i,0);
        for (int j=1;j<ma1;++j)
            cc+= C(i,j);
        C_F0  vi(TheOperators,"[]",cc);
        if(i==0) v=vi;
        else v+= vi;
    }
    return C_F0(TheOperators,"[]",v);


}

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
	    if( ma != (int) a[i].LeftValue()->nbitem())
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
		else     A(j,i) = TryConj((*li)[j]);
	}
	    else
		for (int i=0;i<na;++i)
		    if(!ta)  A(i,0) = a[i];
		    else     A(0,i) = TryConj(a[i]);


    CC_F0 s;
    s= A(0,0); // correction feb. 2016 Thank to O. Pironneau
    for (int i=1;i<na1;++i)
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
	    if( ma != (int) a[i].LeftValue()->nbitem())
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
		else     A(j,i) = TryConj((*li)[j]);
	}
	    else
		for (int i=0;i<na;++i)
		    if(!ta)  A(i,0) = a[i];
		    else     A(0,i) = TryConj(a[i]);


    if(na1==1)
      return  A(0,0);
    else if( na1==2 )
    {
	C_F0 s1(TheOperators,"*",A(0,0),A(1,1));
	C_F0 s2(TheOperators,"*",A(0,1),A(1,0));
	return C_F0(TheOperators,"-",s1,s2);
    }
    else if( na1==3 )
    {
        int i=0,ii=(i+1)%3,iii=(i+2)%3;
        A(i,0)*A(i,0);
        C_F0 det = A(i,0)*A(ii,1)*A(iii,2) - A(i,0)*A(ii,2)*A(iii,1);
        i++;ii=(i+1)%3,iii=(i+2)%3;
        det +=  A(i,0)*A(ii,1)*A(iii,2) - A(i,0)*A(ii,2)*A(iii,1);
        i++;ii=(i+1)%3,iii=(i+2)%3;
        det +=  A(i,0)*A(ii,1)*A(iii,2) - A(i,0)*A(ii,2)*A(iii,1);
        return det;
    }
    else
    {
	CompileError("FH: sorry only det of 1x1 and 2x2 matrix ");
    }
    return  C_F0();

}

//  Add juin  2007
template<class A,class B=A,class R=A>
struct evalE_mul {
    static AnyType eval(Stack s,const E_F0 * ab,const E_F0 * a,const E_F0 * b, bool & meshidenp)
    {
	A aa = GetAny<A>((*a)(s)) ;
	B bb = GetAny<B>((*b)(s)) ;
	R rr(aa*bb);
	bool mia=a->MeshIndependent();
	bool mib=b->MeshIndependent();

	if (( aa == A()) && mia ) meshidenp=true;
	else if(( bb == B()) && mib ) meshidenp=true;
        else meshidenp = mib && mia;
	cout << " meshidenp ??? " << meshidenp << " " << rr << endl;
	return SetAny<R>(static_cast<R>(rr));
    }
};
istream *Getline(istream * f, string ** s)
{
    if( *s==0) *s=newstring();
    getline(*f,**s);
    size_t l = (**s).length();
    if( l > 0 && ((**s)[l-1]=='\r')) (**s).resize(l-1); //
       return f;
}
// Fin Add ne marche pas ....
// fiun avril 2007
//  Hack to Bypass a bug in freefem FH  ...
template<>
class ForEachType<void *>:  public basicForEachType{public:// correction july 2009..... FH  Hoooo....  (Il y a un bug DUR DUR FH  ...)
    ForEachType(Function1 iv=0,Function1 id=0,Function1 OOnReturn=0):basicForEachType(typeid(void *),sizeof(void *),0,0,iv,id,OOnReturn) { }
};

inline double walltime(){
    // add for Pichon mars 2010
    time_t currentWallTime;
    time(&currentWallTime);
    return (double)currentWallTime;
}

inline long fftime()
{
    time_t tloc;
    return time(&tloc);
}
long ffstrtol(string* p)
{
    char * pe;
    const char *pp=p->c_str();
    long r = strtol(pp,&pe,10);
    const char *ppe = pe, *pppe= pp+p->size();
    assert(ppe <= pppe);
    for(const char *ppe = pe; ppe < pppe; ++ppe)
        ffassert(isspace(*ppe));
    return r;

}
long ffstrtol(string* p,long d)
{
    char * pe;
    const char *pp=p->c_str();
    long r = strtol(pp,&pe,d);
    const char *ppe = pe, *pppe= pp+p->size();
    ffassert(ppe <= pppe);
    for(const char *ppe = pe; ppe < pppe; ++ppe)
        ffassert(isspace(*ppe));
    return r;
}

double ffstrtod(string* p)
{
    char * pe;
    const char *pp=p->c_str();
    double r = strtod(pp,&pe);
    const char *ppe = pe, *pppe= pp+p->size();
    ffassert(ppe <= pppe);
    for(const char *ppe = pe; ppe < pppe; ++ppe)
        ffassert(isspace(*ppe));
    return r;

}

long atoi(string* p) {return atol(p->c_str());}// add march 2010
double atof(string* p) {return atof(p->c_str());}// add march 2010
double NaN(string* p) {
return nan(p->c_str());}// add march 2012
double NaN() {return nan("");}// add march 2012
int ShowAlloc(const char *s,size_t & lg);
long ShowAlloc1(string *  s,long * np) { size_t lg; long  n= ShowAlloc(s->c_str(),lg); *np=lg; return n;}
long ShowAlloc1(string *  s) { size_t lg; long  n= ShowAlloc(s->c_str(),lg); return n;}

class E_ForAllLoopMapSS
{  public:
    typedef String K;
    typedef String V;
    typedef string *KK;
    typedef string *VV;

    typedef MyMap<K,V>  *Tab;
    typedef   MyMap<K,V>::iterator TabI ;

    typedef  ForAllLoopOpBase DataL;
    const DataL *data;
    E_ForAllLoopMapSS(const DataL *t): data(t){}
    AnyType f(Stack s) const {
        TabI  ii;
        Tab t= GetAny<Tab >(data->tab(s));

        KK * i   =   GetAny<KK* >(data->i(s));
        VV * v   =   GetAny<VV* >(data->v(s));
        if(verbosity>1000) {
            cout << " i " << (char*) (void *) i -  (char*)(void*) s ;
            cout << " vi " <<  (char*) (void *) v -  (char*)(void*) s ;
            cout << endl;}


        ffassert(i && v);
        if(t->m)
            ii=t->m->begin();
        bool ok = true;
        while(ok)
        {
            if(verbosity>99999) cout << " new n";
            TabI iip=ii++;
            ok = ii != t->m->end();
            String  kk = iip->first;
            String  vv = iip->second;
            const string * pvo=iip->second;

            *i =  kk;
            *v =  vv;
            const string * pv =*v;
            // for  Windows otherwise trap ???? FH. march 2016
            if(verbosity>99999) cout << "  b:" << i << " "<< v  << " "  << kk << " " <<  vv << endl;

            data->code(s);
            if(verbosity>99999) cout << "  a:" << i << " "<< v  << " "  << kk << " " <<  **v << endl;
            if( pvo  == (const string *) iip->second) // no  change m[i]=
            {if( *v != pv ) //  v change
                iip->second  = **v;
            }
            else if( *v != pv )
            {//  v change
                cerr << " Erreur forall change m[i] and mi (stupide please choosse) \n";
                ffassert(0);
            }
            if(verbosity>99999) cout << " A;" << i << " "<< v  << " "  << kk << " " <<  vv << " ok =" << ok << endl;

            *i=0;
            *v=0;
            if(verbosity>99999)
            {
                cout << " 0:" << i << " "<< v  << " "  << kk << " " <<  vv << endl;
                for (TabI iii=t->m->begin();iii!= t->m->end();++iii)
                    cout << " map =" << iii->first << " -> " << iii->second << endl;
                cout << " end of map " << endl;
            }
        }
        return Nothing  ;
    }

};
double projection(const double & aa, const double & bb, const double & x )
{ double a=aa,b=bb; if(a>b) std::swap(a,b); return min(max(a,x),b);}
double dist(const double & aa, const double & bb) { return sqrt( aa*aa+bb*bb);}
double dist(const double & aa, const double & bb,const double & cc) { return sqrt(aa*aa+bb*bb+cc*cc);}
// Add Jan 2017 FH
double diffpos(const double & aa, const double & bb) { return aa<bb ? bb-aa : 0.;}
double invdiffpos(const double & aa, const double & bb) { return aa<bb ? 1./(bb-aa) : 0.;}
double diffnp(const double & aa, const double & bb) { return aa<0. && 0.<bb ? bb-aa : 0.;}//  Corr 08/18  G. Sadaka
double invdiffnp(const double & aa, const double & bb) { return aa<0. && 0.<bb  ? 1./max(bb-aa,1e-30) : 0.;}//  Corr 08/18  G. Sadaka
double invdiff(const double & aa, const double & bb) { double d= aa-bb; return abs(d) < 1e-30 ? d : 1/d;}
double invdiff(const double & aa, const double & bb,const double &eps) { double d= aa-bb; return abs(d) < eps ? d : 1/d;}
extern double ff_tgv; // Add FH jan 2018
double sign(double x){return (x>0.)-(x<0.); }// Add FH jan 2018
long sign(long x){return (x>0)-(x<0); }// Add FH jan 2018
bool ffsignbit(long x){return signbit(x);}
bool ffsignbit(double x){return signbit(x);}
template<typename T>
bool  pswap(T *a,T *b) {swap(*a,*b);return 0; }

void Init_map_type()
{
   TheOperators=new Polymorphic(),
   TheRightOperators=new Polymorphic();
    map_type[typeid(AnyType).name()] = new ForTypeAnyType();
    map_type[typeid(void).name()] = new ForTypeVoid();
       InitLoop();
    Dcl_Type<Expression>(0);
    Dcl_TypeandPtr<double>(0,0,::InitializeDef<double>,0);
    Dcl_TypeandPtr<long>(0,0,::InitializeDef<long>,0);
    Dcl_TypeandPtr<bool>(0,0,::InitializeDef<bool>,0);
    Dcl_TypeandPtr<Complex>(0,0,::InitializeDef<Complex>,0);
    Dcl_Type<void*>(); // add FH ...  for mpi comm world
    Dcl_Type<char*>();
    Dcl_Type<const char *>();
    Dcl_Type<char>();
    Dcl_TypeandPtr<string*>(0,0,::InitializePtr<string*>,::DeletePtr<string*>);
    Dcl_TypeandPtr<ostream*>(0,0,::InitializePtr<ostream*>,::DeletePtr<ostream*>);
    Dcl_TypeandPtr<istream*>(0,0,::InitializePtr<istream*>,::DeletePtr<istream*>);
    Dcl_Type< ostream_precis > ();
    Dcl_Type< ostream_seekp > ();
    Dcl_Type< istream_seekg > ();
    Dcl_Type< istream_good > ();
    Dcl_Type< NothingType > ();
    Dcl_Type<OP_setw>();
    Dcl_Type<Polymorphic*>();

    basicForEachType::type_C_F0 = map_type[typeid(C_F0).name()] =  new TypeLineFunction;
    Dcl_Type<E_Array>();
    Dcl_Type<TransE_Array >();// add
    Dcl_Type<const E_Border *>();
    Dcl_Type<const E_BorderN *>();
    typedef MyMap<String,String> MyMapSS;
    map_type[typeid(MyMapSS*).name()] = new ForEachType<MyMapSS*>(Initialize<MyMapSS >,Delete<MyMapSS >) ;
    map_type_of_map[make_pair(atype<string*>(),atype<string*>())]=atype<MyMapSS*>();


    Dcl_Type<SubArray>();
    Dcl_Type<pair<long,long> >();

    initArrayDCLlong();
    initArrayDCLdouble();
    initArrayDCLComplex();

    Dcl_Type<ios::openmode>();

    // <<known_variable_types>> les types des variables

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
       new E_F1_funcT<long,ostream_precis>(Cast<long,ostream_precis>),
       new E_F1_funcT<long,ostream_seekp>(Cast<long,ostream_seekp>),
       new E_F1_funcT<long,istream_seekg>(Cast<long,istream_seekg>)
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

     Global.New("searchMethod",CPValue<long>(searchMethod)); //pichon
     Global.New("tgv",CPValue<double>(ff_tgv));
     Global.New("lockOrientation",CPValue<bool>(lockOrientation));
     extern long newconvect3;// def in global.cpp
     Global.New("newconvect",CPValue<long>(newconvect3)); //pichon

     // <<cout>> uses [[file:AFunction.hpp::CConstant]]
     Global.New("cout",CConstant<ostream*>(&cout));

     Global.New("cerr",CConstant<ostream*>(&cerr));// add jan 2014 FH.
     Global.New("cin",CConstant<istream*>(&cin));
     Global.New("append",CConstant<ios::openmode>(ios::app));
     Global.New("binary",CConstant<ios::openmode>(ios::binary)); // add FH april 2014
     TheOperators->Add("|",new OneBinaryOperator<Op2_pipe<ios::openmode> >); // add FH april 2014
     Global.New("endl",CConstant<const char*>("\n"));
     Global.New("true",CConstant<bool>(true));
     Global.New("false",CConstant<bool>(false));
     Global.New("pi",CConstant<double>(3.14159265358979323846264338328));
     Global.New("version",CConstant<double>(VersionNumber()));

     Global.New("showCPU",CPValue<bool>(showCPU));
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
       new OneBinaryOperator<Op2_sub<long,long,bool> >,
       new OneBinaryOperator<Op2_sub<long,bool,long> >,
       new OneBinaryOperator<Op2_sub<Complex,Complex,Complex> >,
       new OneBinaryOperator<Op2_sub<Complex,Complex,double> >,
       new OneBinaryOperator<Op2_sub<Complex,double,Complex> >,
       new OneBinaryOperator<Op2_sub<Complex,Complex,long> >,
       new OneBinaryOperator<Op2_sub<Complex,long,Complex> >
       );

     TheOperators->Add("*",
       new OneBinaryOperator<Op2_mul<bool,bool,bool>,OneBinaryOperatorMI,evalE_mul<bool> >,
       new OneBinaryOperator<Op2_mul<long,long,long>,OneBinaryOperatorMI,evalE_mul<long> >,
       new OneBinaryOperator<Op2_mul<double,double,double>,MIMul<double,double>,evalE_mul<double> >,
       new OneBinaryOperator<Op2_mul<double,double,long>, MIMul<double,long>,evalE_mul<double,long,double> >,
       new OneBinaryOperator<Op2_mul<double,long,double>,MIMul<long,double>,evalE_mul<long,double,double>  >,
       new OneBinaryOperator<Op2_mul<Complex,Complex,Complex> >,
       new OneBinaryOperator<Op2_mul<Complex,Complex,double> >,
       new OneBinaryOperator<Op2_mul<Complex,double,Complex> >,
       new OneBinaryOperator<Op2_mul<Complex,Complex,long> >,
       new OneBinaryOperator<Op2_mul<Complex,long,Complex> >,
       new OneBinaryOperator<Op2_mul<long,long,bool> >,
       new OneBinaryOperator<Op2_mul<long,bool,long> >
       );
    // add missing operation jan 2017 FH.
    TheOperators->Add("*",
                      new OneBinaryOperator<Op2_mul<Complex,Complex,bool> >,
                      new OneBinaryOperator<Op2_mul<Complex,bool,Complex> >
                      );
    TheOperators->Add("+",
                      new OneBinaryOperator<Op2_add<Complex,Complex,bool> >,
                      new OneBinaryOperator<Op2_add<Complex,bool,Complex> >
                      );
    TheOperators->Add("-",
                      new OneBinaryOperator<Op2_sub<Complex,Complex,bool> >,
                      new OneBinaryOperator<Op2_sub<Complex,bool,Complex> >
                      );
    TheOperators->Add("/",
                //      new OneBinaryOperator<Op2_sub<Complex,Complex,bool> >,
                      new OneBinaryOperator<Op2_div<Complex,bool,Complex> >
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
		       new OneBinaryOperator<Op2_pow<double,double,double> >,
		       new OneBinaryOperator<Op2_pow<double,double,long> >,
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
       new OneBinaryOperator<set_peqstring ,OneBinaryOperatorMIWO>  // FH string * mars 2006
       );

     TheOperators->Add("?:",
       new Operator_Aritm_If<bool >,
       new Operator_Aritm_If<long >,
       new Operator_Aritm_If<double >,
       new Operator_Aritm_If<Complex >,
       new Operator_Aritm_If<string* >  // (OK???)  to do FH string * mars 2006
       );

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
       new AddBorderOperator
       );
    
      // add frev 2007
      TheOperators->Add("\'", new opTrans);

    TheOperators->Add("*",new opDot(atype<TransE_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add("*",new opDot(atype<E_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add("*",new opColumn(atype<E_Array >() )   );  //  [ ]* C_F0 (all)
      TheOperators->Add("*",new opColumn(basicForEachType::type_C_F0,atype<E_Array >() )   );  //  [ ]* C_F0 (all)
      TheOperators->Add("*",new opColumn(basicForEachType::type_C_F0,atype<TransE_Array >() )   );  //  [ ]* C_F0 (all)
//    type_C_F0
      TheOperators->Add("::",new opColumn(atype<E_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add("*",new opDot(atype<E_Array >(),atype<TransE_Array>() )   );  // a faire mais dur
      TheOperators->Add("*",new opDot(atype<TransE_Array >(),atype<TransE_Array>() )   );  // a faire mais dur

     // car le type de retour depent des objets du tableau
      atype<E_Array >()->Add("[","",new opVI(atype<E_Array >())   );
      atype<TransE_Array >()->Add("[","",new opVI(atype<TransE_Array >())   );
      TheOperators->Add("+",new opSum("+",atype<TransE_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add("+",new opSum("+",atype<E_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add("+",new opSum("+",atype<E_Array >(),atype<TransE_Array>() )   );  // a faire mais dur
      TheOperators->Add("+",new opSum("+",atype<TransE_Array >(),atype<TransE_Array>() )   );  // a faire mais dur
      TheOperators->Add("-",new opSum("-",atype<TransE_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add("-",new opSum("-",atype<E_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add("-",new opSum("-",atype<E_Array >(),atype<TransE_Array>() )   );  // a faire mais dur
      TheOperators->Add("-",new opSum("-",atype<TransE_Array >(),atype<TransE_Array>() )   );  // a faire mais dur
      TheOperators->Add(".*",new opSum("*",atype<TransE_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add(".*",new opSum("*",atype<E_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add(".*",new opSum("*",atype<E_Array >(),atype<TransE_Array>() )   );  // a faire mais dur
      TheOperators->Add(".*",new opSum("*",atype<TransE_Array >(),atype<TransE_Array>() )   );  // a faire mais dur
      TheOperators->Add("./",new opSum("/",atype<TransE_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add("./",new opSum("/",atype<E_Array >(),atype<E_Array>() )   );  // a faire mais dur
      TheOperators->Add("./",new opSum("/",atype<E_Array >(),atype<TransE_Array>() )   );  // a faire mais dur
    // correct in sept. 2009
      TheOperators->Add("./",new opSum("/",atype<TransE_Array >(),atype<TransE_Array>() )   );  // a faire mais dur


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
       new OneBinaryOperator<Op2_set_pstring<istream**,ifstream> >,  //  FH string * mars 2006
       new OneBinaryOperator<Op2_set_pstring<ostream**,ofstream> >,  //  FH string * mars 2006
       new OneTernaryOperator3<Op2_set_pstringiomode<ostream**,ofstream> >  ,    //  FH string * mars 2006
       new OneTernaryOperator3<Op2_set_pstringiomode<istream**,ifstream> >   //  FH string * april  2014
       );

     atype<istream* >()->AddCast( new E_F1_funcT<istream*,istream**>(UnRef<istream* >));
     atype<ostream* >()->AddCast( new E_F1_funcT<ostream*,ostream**>(UnRef<ostream* >));

     Add<ostream**>("<-","(", new OneUnaryOperator<Op1_new_pstring<ostream*,ofstream> >);  //  FH string * mars 2006

     Add<ostream**>("precision",".",new OneOperator1<ostream_precis,ostream**>(ostream_precision));
     Add<ostream*>("precision",".",new OneOperator1<ostream_precis,ostream*>(ostream_precision));

    // add FH jan 2010 ...
    Add<ostream**>("seekp",".",new OneOperator1<ostream_seekp,ostream**>(ff_oseekp));
    Add<ostream*>("seekp",".",new OneOperator1<ostream_seekp,ostream*>(ff_oseekp));

    Add<istream**>("seekg",".",new OneOperator1<istream_seekg,istream**>(ff_iseekg));
    Add<istream*>("seekg",".",new OneOperator1<istream_seekg,istream*>(ff_iseekg));
    Add<ostream**>("tellp",".",new OneOperator1<ostream_seekp,ostream**>(ff_oseekp));
    Add<ostream*>("tellp",".",new OneOperator1<ostream_seekp,ostream*>(ff_oseekp));

    Add<istream**>("tellg",".",new OneOperator1<istream_seekg,istream**>(ff_iseekg));
    Add<istream*>("tellg",".",new OneOperator1<istream_seekg,istream*>(ff_iseekg));

    Add<ostream_seekp>("(","",new OneOperator1<long,ostream_seekp>(fftellp),
		       new OneOperator2<long,ostream_seekp,long>(ffseekp));
    Add<istream_seekg>("(","",new OneOperator1<long,istream_seekg>(fftellg),
		       new OneOperator2<long,istream_seekg,long>(ffseekg));
    // end add  jan 2010 ..
    Add<ostream_precis>("(","",new OneOperator1<long,ostream_precis>(get_precis),
                                new OneOperator2<long,ostream_precis,long>(set_precis));
//  add v 1.41
     Add<istream**>("good",".",new OneOperator1<istream_good,istream**>(to_istream_good));
     Add<istream*>("good",".",new OneOperator1<istream_good,istream*>(to_istream_good));
     Add<istream*>("good",".",new OneOperator1<istream_good,istream*>(to_istream_good));
     Add<istream_good>("(","",new OneOperator1<long,istream_good>(get_good));

     Add<istream**>("eof",".",new OneOperator1<bool,istream**>(get_eof));
// add v 4.5 jan 2020 FH.
        Add<istream**>("eatspace",".",new OneOperator1<bool,istream**>(eatspace));
       Add<istream**>("length",".",new OneOperator1<long,istream**>(filelength));
// add v 2.8
     Add<ostream**>("scientific",".",new OneOperator1<ostream**,ostream**>(set_os<scientific>));
     Add<ostream**>("fixed",".",new OneOperator1<ostream**,ostream**>(set_os<fixed>));
     Add<ostream**>("showbase",".",new OneOperator1<ostream**,ostream**>(set_os<showbase>));
     Add<ostream**>("noshowbase",".",new OneOperator1<ostream**,ostream**>(set_os<noshowbase>));
     Add<ostream**>("showpos",".",new OneOperator1<ostream**,ostream**>(set_os<showpos>));
     Add<ostream**>("noshowpos",".",new OneOperator1<ostream**,ostream**>(set_os<noshowpos>));
     Add<ostream**>("default",".",new OneOperator1<ostream**,ostream**>(set_os<default1>));
     Add<ostream**>("flush",".",new OneOperator1<ostream**,ostream**>(set_os_flush));// ADD may 2010

     Add<ostream*>("scientific",".",new OneOperator1<ostream*,ostream*>(set_os1<scientific>));
     Add<ostream*>("fixed",".",new OneOperator1<ostream*,ostream*>(set_os1<fixed>));
     Add<ostream*>("showbase",".",new OneOperator1<ostream*,ostream*>(set_os1<showbase>));
     Add<ostream*>("noshowbase",".",new OneOperator1<ostream*,ostream*>(set_os1<noshowbase>));
     Add<ostream*>("showpos",".",new OneOperator1<ostream*,ostream*>(set_os1<showpos>));
     Add<ostream*>("noshowpos",".",new OneOperator1<ostream*,ostream*>(set_os1<noshowpos>));
     Add<ostream*>("default",".",new OneOperator1<ostream*,ostream*>(set_os1<default1>));
     Add<ostream*>("flush",".",new OneOperator1<ostream*,ostream*>(set_os_flush));// ADD may 2010

    Global.Add("getline","(",new OneOperator2<istream*,istream*,string **>(Getline));
// add 2.16
     Global.Add("trace","(",new opFormal(atype<E_Array>(),formalMatTrace ));
     Global.Add("det","(",new opFormal(atype<E_Array>(),formalMatDet ));
// end add

    // add 3.20
    Global.Add("Cofactor","(",new opFormal(atype<E_Array>(),formalMatCofactor ));

     TheOperators->Add("[]",new OneOperator_array );
     TheOperators->Add("[border]",new OneOperator_border );

     Global.Add("cos","(",new OneOperator1<double>(cos));
    Global.Add("square","(",new OneOperator1<long,long,E_F_F0<long,const long &> >(Square));// add FH Mai 2011
    Global.Add("square","(",new OneOperator1<double,double,E_F_F0<double,const double &> >(Square));
    Global.Add("square","(",new OneOperator1<Complex,Complex,E_F_F0<Complex,const Complex &> >(Square));// add FH Mai 2011
 //add for Olivier FH July 2014
    Global.Add("sqr","(",new OneOperator1<long,long,E_F_F0<long,const long &> >(Square));//
    Global.Add("sqr","(",new OneOperator1<double,double,E_F_F0<double,const double &> >(Square));
    Global.Add("sqr","(",new OneOperator1<Complex,Complex,E_F_F0<Complex,const Complex &> >(Square));//

     Global.Add("round","(",new OneOperator1<double>(round)); // add june 2007
     Global.Add("lround","(",new OneOperator1<long,double>(lround)); // add june 2007
     Global.Add("floor","(",new OneOperator1<double>(floor)); // add march 2006
     Global.Add("ceil","(",new OneOperator1<double>(ceil));  // add march 2006
    Global.Add("rint","(",new OneOperator1<double>(rint));  // add june 2006
     Global.Add("lrint","(",new OneOperator1<long,double>(lrint));  // add mars  2014

     Global.Add("sin","(",new OneOperator1<double>(sin));
     Global.Add("tan","(",new OneOperator1<double>(tan));
     Global.Add("atan","(",new OneOperator1<double>(atan));
     Global.Add("sinh","(",new OneOperator1<double>(sinh));
     Global.Add("cosh","(",new OneOperator1<double>(cosh));
     Global.Add("tanh","(",new OneOperator1<double>(tanh));

    Global.Add("atoi","(",new OneOperator1<long,string*>(atoi));// add march 2010
    Global.Add("atol","(",new OneOperator1<long,string*>(atoi));// add march 2010
    Global.Add("atof","(",new OneOperator1<double,string*>(atof));// add march 2010

    Global.Add("strtol","(",new OneOperator1<long,string*>(ffstrtol));// add march 2017
    Global.Add("strtol","(",new OneOperator2<long,string*,long>(ffstrtol));// add march 2017
    Global.Add("strtod","(",new OneOperator1<double,string*>(ffstrtod));// add march 2017

     Global.Add("atanh","(",new OneOperator1<double>(atanh));
     Global.Add("asin","(",new OneOperator1<double>(asin));
     Global.Add("acos","(",new OneOperator1<double>(acos));
     Global.Add("asinh","(",new OneOperator1<double>(asinh));
     Global.Add("acosh","(",new OneOperator1<double>(acosh));
#ifdef HAVE_ERFC
     Global.Add("erf","(",new OneOperator1<double>(erf));
     Global.Add("erfc","(",new OneOperator1<double>(erfc));
#endif
#ifdef HAVE_TGAMMA
     Global.Add("tgamma","(",new OneOperator1<double>(tgamma));
     Global.Add("lgamma","(",new OneOperator1<double>(lgamma));
#endif
     //  function de bessel j0, j1, jn, y0, y1, yn -- bessel functions of first and second kind
#ifdef HAVE_JN
      Global.Add("j0","(",new OneOperator1<double>(j0));
      Global.Add("j1","(",new OneOperator1<double>(j1));
      Global.Add("jn","(",new OneOperator2<double,long,double>(myjn));
      Global.Add("y0","(",new OneOperator1<double>(y0));
      Global.Add("y1","(",new OneOperator1<double>(y1));
      Global.Add("yn","(",new OneOperator2<double,long,double>(myyn));
#endif
     Global.Add("exp","(",new OneOperator1<double>(exp));
     Global.Add("log","(",new OneOperator1<double>(log));
     Global.Add("log10","(",new OneOperator1<double>(log10));
     Global.Add("pow","(",new OneOperator2<double,double>(pow));
//     Global.Add("pow","(",new OneOperator2<double,double,long>(pow));
     Global.Add("max","(",new OneOperator2_<double,double>(Max<double> ));
     Global.Add("min","(",new OneOperator2_<double,double>(Min<double> ));
    Global.Add("diffpos","(",new OneOperator2_<double,double>(diffpos )); // jan 2018 FH
    Global.Add("invdiffpos","(",new OneOperator2_<double,double>(invdiffpos )); // jan 2018 FH
    Global.Add("diffnp","(",new OneOperator2_<double,double>(diffnp )); // jan 2018 FH
    Global.Add("invdiffnp","(",new OneOperator2_<double,double>(invdiffnp )); // jan 2018 FH
    Global.Add("invdiff","(",new OneOperator2_<double,double>(invdiff )); // jan 2018 FH
    Global.Add("invdiff","(",new OneOperator3_<double,double,double>(invdiff )); // jan 2018 FH

     Global.Add("max","(",new OneOperator2_<long,long>(Max));
     Global.Add("min","(",new OneOperator2_<long,long>(Min));
    Global.Add("max","(",new OneOperator3_<double,double>(Max<double> ));
    Global.Add("min","(",new OneOperator3_<double,double>(Min<double> ));
    Global.Add("max","(",new OneOperator3_<long,long>(Max));
    Global.Add("min","(",new OneOperator3_<long,long>(Min));
    Global.Add("max","(",new OneOperator4_<long,long>(Max));
    Global.Add("min","(",new OneOperator4_<long,long>(Min));
    Global.Add("max","(",new OneOperator4_<double,double>(Max));
    Global.Add("min","(",new OneOperator4_<double,double>(Min));

    Global.Add("atan2","(",new OneOperator2<double>(atan2));
    Global.Add("fmod","(",new OneOperator2<double>(fmod));// add sep 2017
    Global.Add("fdim","(",new OneOperator2<double>(fdim));// add sep 2017
    Global.Add("fmax","(",new OneOperator2<double>(fmax));// add sep 2017
    Global.Add("fmin","(",new OneOperator2<double>(fmin));// add sep 2017

    Global.Add("hypot","(",new OneOperator2<double>(hypot));// add Jan 2014

     Global.Add("atan","(",new OneOperator2<double>(atan2));
     Global.Add("sqrt","(",new OneOperator1<double>(sqrt,2));
     Global.Add("abs","(",new OneOperator1<double>(Abs));
     Global.Add("abs","(",new OneOperator1<long>(Abs));
     Global.Add("cos","(",new OneOperator1_<Complex>(cos));
     Global.Add("sin","(",new OneOperator1_<Complex>(sin));
     Global.Add("sinh","(",new OneOperator1_<Complex>(sinh));
     Global.Add("cosh","(",new OneOperator1_<Complex>(cosh));
     Global.Add("tanh","(",new OneOperator1_<Complex>(tanh));// Add June 2016 FH..
     Global.Add("log","(",new OneOperator1_<Complex>(log));
     Global.Add("tan","(",new OneOperator1_<Complex>(tan));
     Global.Add("exp","(",new OneOperator1_<Complex>(exp));

    Global.Add("pow","(",new OneBinaryOperator<Op2_pow<Complex,Complex,Complex> >);
                //new OneOperator2_<Complex,Complex>(pow ));
     Global.Add("sqrt","(",new OneOperator1_<Complex>(sqrt,0));
     Global.Add("conj","(",new OneOperator1_<Complex>(conj,0));
     Global.Add("conj","(",new OneOperator1_<double>(RNM::conj,1));
     TheOperators->Add("\'",new OneOperator1_<Complex>(conj,0));
     TheOperators->Add("\'",new OneOperator1_<double>(RNM::conj,1));       //  add F.  Feb 2010  of conj of varf..


     Global.Add("imag","(",new OneOperator1_<double,Complex>(Imag));
     //  Big probleme  real is a type
     Add<double>("<--","(",new OneOperator1_<double,Complex>(Real));

     Global.Add("abs","(",new OneOperator1_<double,Complex>(abs));

     Global.Add("arg","(",new OneOperator1_<double,Complex>(arg));
     Global.Add("norm","(",new OneOperator1_<double,Complex>(norm));
     Global.Add("exit","(",new OneOperator1<long>(Exit));
     Global.Add("assert","(",new OneOperator1<bool>(Assert));

     Global.Add("clock","(",new OneOperator0<double>(CPUtime));
    Global.Add("time","(",new OneOperator0<double>(walltime));// add mars 2010 for Pichon.
    Global.Add("ltime","(",new OneOperator0<long>(fftime));// add mars 2014 ( the times unix fonction)
    Global.Add("storageused","(",new OneOperator0<long>(storageused));
    Global.Add("storagetotal","(",new OneOperator0<long>(storagetotal));

     Global.Add("dumptable","(",new OneOperator1<ostream*,ostream*>(dumptable));
     Global.Add("exec","(",new OneOperator1<long,string* >(exec));  //FH string * mars 2006
     Global.Add("system","(",new OneOperator1<long,string* >(exec));  //FH string fevr 2011

     Global.Add("polar","(",new OneOperator2_<Complex,double,double>(polar));
 // rand generator ---
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
  init_by_array(init, length);

  Global.Add("randint32","(",new OneOperator_0<long>(genrandint32));
  Global.Add("randint31","(",new OneOperator_0<long>(genrand_int31));
  Global.Add("randreal1","(",new OneOperator_0<double>(genrand_real1));
  Global.Add("randreal2","(",new OneOperator_0<double>(genrand_real2));
  Global.Add("randreal3","(",new OneOperator_0<double>(genrand_real3));
  Global.Add("randres53","(",new OneOperator_0<double>(genrand_res53));
  Global.Add("randinit","(",new OneOperator1<long>(genrandint));

  //  NaN and Inf
  Global.Add("ShowAlloc","(",new OneOperator1<long,string*>(ShowAlloc1));// debuging
  Global.Add("ShowAlloc","(",new OneOperator2<long,string*,long*>(ShowAlloc1));// debuging
  Global.Add("NaN","(",new OneOperator0<double>(NaN));

  Global.Add("NaN","(",new OneOperator1<double,string*   >(NaN));
  Global.Add("isNaN","(",new OneOperator1<long,double>(isNaN));
  Global.Add("copysign","(",new OneOperator2<double>(copysign));// Add jan 2018 FH
  Global.Add("sign","(",new OneOperator1<double>(sign));// Add jan 2018 FH
  Global.Add("sign","(",new OneOperator1<long>(sign));// Add jan 2018 FH
  Global.Add("signbit","(",new OneOperator1<bool,long>(ffsignbit));// Add jan 2018 FH
  Global.Add("signbit","(",new OneOperator1<bool,double>(ffsignbit));// Add jan 2018 FH

  Global.Add("isInf","(",new OneOperator1<long,double>(isInf));
  Global.Add("isNormal","(",new OneOperator1<long,double>(isNormal));
  Global.Add("chtmpdir","(",new OneOperator0<long>(ffapi::chtmpdir));
  Global.Add("projection","(",new OneOperator3_<double,double   >(projection));
  Global.Add("dist","(",new OneOperator2_<double,double>(dist));
  Global.Add("dist","(",new OneOperator3_<double,double>(dist));
  Global.Add("swap","(",new OneOperator2<bool,double*>(pswap));
  Global.Add("swap","(",new OneOperator2<bool,long*>(pswap));
  Global.Add("swap","(",new OneOperator2<bool,bool*>(pswap));
  Global.Add("swap","(",new OneOperator2<bool,Complex*>(pswap));
  Global.Add("swap","(",new OneOperator2<bool,string**>(pswap));


  atype<MyMapSS*>()->Add("[","",new OneOperator2_<string**,MyMapSS*,string*>(get_elements));

  atype<MyMapSS*>()->SetTypeLoop(atype<string**>(),atype<string**>());


  tables_of_identifier.push_back(&Global);

  TheOperators->Add("<<",new OneBinaryOperator<PrintP<MyMapSS*> >);

  TheOperators->Add("{}",new ForAllLoop<E_ForAllLoopMapSS >);
  // add setw feb 2015 FH
  Global.Add("setw","(",new OneOperator1<OP_setw,long>(defOP_setw));
  TheOperators->Add("<<", new OneBinaryOperator<Print<OP_setw> >);

}




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

// <<addingInitFunct>>
static addingInitFunct TheaddingInitFunct(-10000,Init_map_type);

C_F0  opVI::code2(const basicAC_F0 &args) const
{
    Expression p=args[1];
    if ( ! p->EvaluableWithOutStack() )
    {
        CompileError(" [...][p], The p must be a constant , sorry");
    }
    int pv = GetAny<long>((*p)(NullStack));
    bool ta =args[0].left()==atype<TransE_Array>();
    const TransE_Array * tea=0;
    const E_Array * ea=0;
    if( ta)  tea = dynamic_cast<const TransE_Array*>((Expression) args[0]);
    else ea = dynamic_cast<const E_Array*>((Expression) args[0]);
    assert( ea || tea );
    const E_Array & a=  ta ? *tea->v : *ea;
    if(tea){
        AC_F0  v  ;
        for (int i=0;i<a.size();++i)
        {
            const E_Array * li=  dynamic_cast<const E_Array *>(a[i].LeftValue());
            ffassert(li && (li->size() >pv));

            const C_F0 vi = TryConj( (*li)[pv]);
            if(i==0) v=vi;
            else v+= vi;

        }

        return C_F0(TheOperators,"[]",v);

    }

    if(!(pv >=0 && pv <a.size()))
    {
        cerr << "\n\nerror [ ... ][" << pv <<" ] " << " the  size of [ ...]  is "<< a.size() << endl;
        lgerror(" bound of  [ .., .. , ..][ . ] operation  ");
    }

    return (* a.v)[pv];
}

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
	    if( ma != (int) a[i].LeftValue()->nbitem())
		CompileError(" first matrix with variable number of columm");

    }
    if(mab) {
	mb= b[1].LeftValue()->nbitem();
	for (int i=1;i<nb;i++)
	    if( mb != (int) b[i].LeftValue()->nbitem())
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
		else     A(j,i) = TryConj((*li)[j]);
	}
    else
	for (int i=0;i<na;++i)
	    if(!ta)  A(i,0) = a[i];
	    else     A(0,i) = TryConj(a[i]);

    if(mab)
	for (int i=0;i<nb;++i)
	{
	    const E_Array * li=  dynamic_cast<const E_Array *>(b[i].LeftValue());
	    ffassert(li);
	    for (int j=0; j<mb;++j)
		if(!tb)  B(i,j) = (*li)[j];
		else     B(j,i) = TryConj((*li)[j]);
	}
    else
	for (int i=0;i<nb;++i)
	    if(!tb)  B(i,0) = b[i];
	    else     B(0,i) = TryConj(b[i]);

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
    else if ( mb1 ==1 ) // || (na1==1)) // correct du car ' on conj encore r . mars 2010
    {
	AC_F0  v;
	v=C(0,0);
	int i0=na1!=1,j0=mb1!=1, nn= mb1*na1;
	for (int i=1;i<nn;++i)
	    v+=C(i0*i,j0*i);
	C_F0  r(TheOperators,"[]",v);
	if(mb1==1) return r;
	else return C_F0(TheOperators,"\'",r);// Bug car on conj encore r . mars 2010
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

    cout << "   formal prod array or matrix : [ .. ] * [ .. ]   " << na << "x" << nb << endl;
    cout << "   formal prod array or matrix : [ .. ] * [ .. ]   " <<  endl;
    cout << " first  array :  matrix " << maa << " trans " << ta << " " << na << "x" << ma <<endl;
    cout << " second array :  matrix " << mab << " trans " << tb << " " << nb << "x" << mb <<endl;
    CompileError("  not implemented sorry ..... (FH) to do ???? ");
    return C_F0();

}
C_F0  opColumn::code2(const basicAC_F0 &args) const
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

    if( (eb || teb) && ( ea || tea ) )
    {
        const E_Array & a=  ta ? *tea->v : *ea;
        const E_Array & b=  tb ? *teb->v : *eb;
        int ma =1;
        int mb =1;
        int na=a.size();
        int nb=b.size();
        if(na <1 && nb < 1) CompileError(" empty array  [ ...]':[ ...  ]  ");
        bool mab= b[0].left()==atype<E_Array>();
        bool maa= a[0].left()==atype<E_Array>();
        if(maa) {
            ma= a[0].LeftValue()->nbitem();
            for (int i=1;i<na;i++)
                if( ma != (int) a[i].LeftValue()->nbitem())
                    CompileError(" first matrix with variable number of columm");

        }
        if(mab) {
            mb= b[1].LeftValue()->nbitem();
            for (int i=1;i<nb;i++)
                if( mb != (int) b[i].LeftValue()->nbitem())
                    CompileError(" second matrix with variable number of columm");
        }
        int na1=na,ma1=ma,nb1=nb,mb1=mb;
        if(ta) RNM::Exchange(na1,ma1);
        if(tb) RNM::Exchange(nb1,mb1);

        KNM<CC_F0> A(na1,ma1), B(nb1,mb1);
        if ( (na1!=nb1 ) || (ma1 != mb1) || (na1 * ma1 ==0)  )
        {
            cout << "\n   formal  array or matrix : [ .. ] : [ .. ]   " << endl;
            cout << " first  array :  matrix " << maa << " trans " << ta << " " << na << "x" << ma <<endl;
            cout << " second array :  matrix " << mab << " trans " << tb << " " << nb << "x" << mb <<endl;
            CompileError(" no same size  [ ...] : [ ...  ] sorry ");
        }

        if(maa)
            for (int i=0;i<na;++i)
            {
                const E_Array * li=  dynamic_cast<const E_Array *>(a[i].LeftValue());
                ffassert(li);
                for (int j=0; j<ma;++j)
                    if(!ta)  A(i,j) = (*li)[j];
                    else     A(j,i) = TryConj((*li)[j]);
            }
        else
            for (int i=0;i<na;++i)
                if(!ta)  A(i,0) = a[i];
                else     A(0,i) = TryConj(a[i]);

        if(mab)
            for (int i=0;i<nb;++i)
            {
                const E_Array * li=  dynamic_cast<const E_Array *>(b[i].LeftValue());
                ffassert(li);
                for (int j=0; j<mb;++j)
                    if(!tb)  B(i,j) = (*li)[j];
                    else     B(j,i) = TryConj((*li)[j]);
            }
        else
            for (int i=0;i<nb;++i)
                if(!tb)  B(i,0) = b[i];
                else     B(0,i) = TryConj(b[i]);

        CC_F0 s,aibi;

        for (int i=0;i<na1;++i)
            for (int j=0;j<ma1;++j)
            {
                aibi = C_F0(TheOperators,"*",A(i,j),B(i,j));
                if( (i==0) && (j==0))
                    s = aibi;
                else
                    s = C_F0(TheOperators,"+",s,aibi);
            };
        return s;
    }
    else if ( ea || tea )
    { // modif 2 /08/  2013  FH .. bug in [ a0,a1,... ]'*b
        //  [a0,a1,... ]*b  or [ a0,a1,... ]'*b  => [ a0*b',a1*b',
        const E_Array & a=  ta ? *tea->v : *ea;
        int na=a.size();
        AC_F0  v;
        v = 0; // empty
        C_F0 b =ta ? TryConj(args[1]) :args[1];
        for (int i=0;i<na;++i)
        v += C_F0(TheOperators,"*",a[i],b)  ;
        return ta ? C_F0(TheOperators,"\'",C_F0(TheOperators,"[]",v)) :   C_F0(TheOperators,"[]",v);

    }
    else if(eb || teb)
    {  // modif 2 /08/  2013  FH .. bug in a*[ b0,b1,... ]'
        const E_Array & b=  tb ? *teb->v : *eb;
        int nb=b.size();
        C_F0 a =tb ? TryConj(args[0]) :args[0];
        AC_F0  v;
        v = 0; // empty
        for (int i=0;i<nb;++i)
            v += C_F0(TheOperators,"*",a,b[i]) ;
        return tb ? C_F0(TheOperators,"\'",C_F0(TheOperators,"[]",v)) :   C_F0(TheOperators,"[]",v);

    }
    else ffassert(0);
    return C_F0();
}


C_F0  opSum::code2(const basicAC_F0 &args) const
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
    int na=a.size();
    int nb=b.size();
    if(na != nb) CompileError(" formal   [ [...] [] ] : [ [..], [..] , ... ]  ");


    AC_F0  v;
    v = 0; // empty
	for (int i=0;i<na;++i)
	    v += C_F0(TheOperators,op,ta ? TryConj(a[i]) : a[i],tb ? TryConj(b[i]): b[i]) ;
	return C_F0(TheOperators,"[]",v);

}
//  to be sure new and delele be in see dll for windows
string  *newstring(){string * p=new string();
    if(verbosity>999999) cout << p << "=newstring() " <<endl;;
      return p ;
}
string  *newstring(const string & c){
    string * p=new string(c);
    if(verbosity>999999) cout <<  p << "=newstring((string) "<< c << ") \n";
   return p;}
string  *newstring(const char * c){
     string * p=new string(c);
    if(verbosity>999999) cout << p << "=newstring((char*)  "<< c << " ) \n";
    return p;}
void   freestring(const string * c){
    if(verbosity>999999) cout << c << "freestring(  "<< *c <<") \n";
    delete c;}
