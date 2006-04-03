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

class SubString {  public:
  string * s;
  long i,j;
  SubString(string **ss,long ii,long jj) :s(*ss),i(ii),j(jj)
     {}
  SubString(string **ss,const SubArray & sb) :s(*ss),i(sb.start),j(sb.n)
     {ffassert(sb.step==1);}
};


extern Map_type_of_map map_type_of_map ; //  to store te type 
extern Map_type_of_map map_pair_of_type ; //  to store te type 

extern basicForEachType *  typevarreal,  * typevarcomplex;  //  type of real and complex variable

extern int TheCurrentLine; // unset: by default
extern long mpisize,mpirank;


long get_size(string * p)  { ffassert(p ) ;  return p->size() ;}
long get_sizep(string ** p)  { ffassert(p && *p) ;  return (*p)->size() ;}



string ** get_replace(string ** pp,long i,long j, string *rr) 
    {
     ffassert(pp && *pp ) ; 
      string & s=**pp;
      s.replace(i,j,*rr);
      return pp;
    } 
// a( : ) = "sqsd";

struct set_substring: public binary_function<SubString,string *,SubString> {

  static SubString f(SubString const  & a, string  *const  & b)  
  {  
     a.s->replace(a.i,a.j,*b);
     return a;}
};    

SubString fSubString( string ** const & a,const SubArray & b)
 { return SubString(a,b);}
 
template<bool B>
 struct String_find {
   string *p;
   String_find(string *pp) :p(pp){ffassert(p);}
   String_find(string **pp) :p(*pp){ffassert(p);}
   
   long find(string *f) const 
      { return p->find(*f);}
   long find(string *f,long i) const 
      { return p->find(*f,i);}
   
 }; 

// spacialisation find -> rfind  (bofbof ??) 
template<>
 struct String_find<false> {
   string *p;
   String_find(string *pp) :p(pp){ffassert(p);}
   String_find(string **pp) :p(*pp){ffassert(p);}
   
   long find(string *f) const 
      { return p->rfind(*f);}
   long find(string *f,long i) const 
      { return p->rfind(*f,i);}
   
 }; 

template<bool B>
String_find<B> to_String_find(string * p)
  { return String_find<B>(p);}

template<bool B>
String_find<B> to_String_findp(string ** p)
  { return String_find<B>(*p);}

template<bool B>
 long string_find(String_find<B> sf,string *s)
 {  return sf.find(s); }
template<bool B>
 long string_find( String_find<B> const & sf,string *const &s,long const &i)
 {  return sf.find(s,i); }
 
istream* getlinep( istream* f, string ** s)
 {
   getline(*f,**s);
   return f;
 }
void     initStringOperator()
{

  Dcl_Type<SubString>();
  aType  tstringp =atype<string*>();
  aType  tstringpp =atype<string**>();
  
   Dcl_Type< String_find<true> > ();
   Dcl_Type< String_find<false> > ();

// ajoute pour les chaines
 Add<string **>("size",".",new OneOperator1<long,string **>(get_sizep) );
 Add<string **>("length",".",new OneOperator1<long,string **>(get_sizep) );
 
 Add<string *>("size",".",new OneOperator1<long,string *>(get_size) );
 Add<string *>("length",".",new OneOperator1<long,string *>(get_size) );
 
  Add<string *>("find",".",new OneOperator1<String_find<true>,string * >(to_String_find<true> ));
  Add<string *>("rfind",".",new OneOperator1<String_find<false>,string * >(to_String_find<false> ));
  Add<string **>("find",".",new OneOperator1<String_find<true>,string ** >(to_String_findp<true> ));
  Add<string **>("rfind",".",new OneOperator1<String_find<false>,string ** >(to_String_findp<false> ));
  
  
 
  Add<String_find<true> >("(","",new OneOperator2<long,String_find<true>,string *>(string_find ) );
  Add<String_find<false> >("(","",new OneOperator2<long,String_find<false>,string *>(string_find ) );

  Add<String_find<true> >("(","",new OneOperator3_<long,String_find<true>,string *,long>(string_find ) );
  Add<String_find<false> >("(","",new OneOperator3_<long,String_find<false>,string *,long>(string_find ) );


 TheOperators->Add("=",
         new OneBinaryOperator<set_substring > );

 Add<string**>("(","",new OneOperator2_<SubString,string **,SubArray>(fSubString)); 

 TheOperators->Add("getline",new OneOperator2<istream*,istream*,string**>(getlinep));
  
 
 
// Add<string**>("[","",new OneOperator2_<SubString,string **,SubArray>(fSubString)); 
 
 //Add<string **>("rfind",".",new OneOperator4_<long,string **,long,long,string *>(get_replace) );
 
}
