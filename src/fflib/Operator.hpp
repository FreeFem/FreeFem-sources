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

#ifndef Operator_hpp_
#define Operator_hpp_

#if defined(__GNUC__) && __GNUC__+0 >= 3
inline double pow(double x,long l) { return pow(x,(double)l);}
#endif


template<class R,class A=R> 
struct Op1_neg: public unary_function<A,R> { 
  static R f(const A & a)  { return - (R)a;} }; 
  
template<class R,class A=R> 
struct Op1_plus: public unary_function<A,R> { 
  static R f(const A & a)  { return + (R)a;} }; 
  
template<class A> 
struct Op1_not: public unary_function<A,bool> { 
  static bool f(const A & a)  { return ! (bool)a;} }; 
  
template<class R,class A=R,class B=A> 
struct Op2_add: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return ((R)a + (R)b);} }; 

template<class R,class A=R,class B=A> 
struct Op2_sub: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return ((R)a - (R)b);} }; 

template<class R,class A=R,class B=A> 
struct Op2_mul: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { 
  // cout << a << " * " << b <<" => "  << ((R)a * (R)b) << endl;
  return ((R)a * (R)b);} }; 

template<class R,class A=R,class B=A> 
struct Op2_div: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  {
     if (b == B()) 
       {cerr <<  a << "/" << b << " : " <<  typeid(A).name()  << " " << typeid(B).name() 
             << " " << typeid(R).name() << endl;ExecError(" Div by 0");}
     return ((R)a / (R)b);} }; 

template<class R>
struct Op2_pipe: public binary_function<R,R,R> {
    static R f(const R & a,const R & b)  {   return (a | b);} };

template<class R,class A=R,class B=A> 
struct Op2_mod: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return ((R)a % (R)b);} }; 

template<class A,class B=A> 
struct Op2_lt: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { 
//    cout << a << " < " << b << " = " << ( a<b) << endl;
    return a  < b;} }; 

template<class A,class B=A> 
struct Op2_le: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { return a  <= b;} }; 


template<class A,class B=A> 
struct Op2_gt: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { return a  > b;} }; 


template<class A,class B=A> 
struct Op2_ge: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { return a  >= b;} }; 

template<class A,class B=A> 
struct Op2_eq: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b) 
   { //cout << a << " == " << b << " => " <<( a  == b) << endl;
   return a  == b;} }; 

template<class A,class B=A> 
struct Op2_ne: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { return a  != b;} }; 

struct Op2_and: public binary_function<bool,bool,bool> { 
  static bool f(const bool & a,const bool & b)  { return a  && b;} }; 
  
struct Op2_or: public binary_function<bool,bool,bool> { 
  static bool f(const bool & a,const bool & b)  { return a  || b;} }; 


template<class R,class A,class B> 
struct Op2_padd: public binary_function<A,B,R*> { 
  static R * f(Stack s,const A & a,const B & b)  { 
   R* r= Add2StackOfPtr2Free(s, a || b ? new R ((a ? *a : nullptr) + (b ? *b : nullptr)) : nullptr);
  // delete a,delete b;
  return r;} }; 


template<class A,class B=A> 
struct Op2_plt: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  < *b;
  //delete a,delete b;
  return r;} }; 

template<class A,class B=A> 
struct Op2_ple: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  <= *b;
  // delete a,delete b;
  return r;} }; 


template<class A,class B=A> 
struct Op2_pgt: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  > *b;
 // delete a,delete b;
  return r;} }; 


template<class A,class B=A> 
struct Op2_pge: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  >= *b;
  // delete a,delete b;
  return r;} }; 

template<class A,class B=A> 
struct Op2_peq: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  == *b;
 //  delete a,delete b;
  return r;} }; 

template<class A,class B=A> 
struct Op2_pne: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  {  bool r=*a  != *b;
  // delete a,delete b;
  return r;} }; 


template<class R,class A=R,class B=A>
struct Op2_pow: public binary_function<A,B,R> {
  static R f(const A & a,const B & b)  { return R(pow(a,b));}};
  


template<class A>
struct Op_Read : public binary_function<istream*,A*,istream*> {
  static istream *  f(istream  * const  & f,A  * const  &  a)  
   {
       if( !f || !*f) ExecError("Fatal Error: file not open in read value (Op_Read)");
       *f >> *a;
       if(!f->good()) ExecError("Fatal Error: file  not good in read array (Op_ReadKN)");
     return f;
   }
};

template<class A>
struct Op_ReadP : public binary_function<istream*,A**,istream*> {
  static istream *  f(istream  * const  & f,A  ** const  &  a)  
   {
     assert(a);
     if( ! *a)  *a= new A ;
     if( !f || !*f) ExecError("Fatal Error: file not open in read value (Op_Read)");
     *f >> **a;
     if(!f->good()) ExecError("Fatal Error: file  not good in read array (Op_ReadKN)");
     return f;
   }
};

template<class A>
struct Op_ReadKN : public binary_function<istream*,KN<A>*,istream*> {
  static istream *  f(istream  * const  & f,KN<A>* const  &  a)  
   { 
     if( !f || !*f) ExecError("Fatal Error: file not open in read array (Op_ReadKN)");
     int n;char c;
     *f >> n;
     if(!f->good()) ExecError("Fatal Error: file  not good in read array (Op_ReadKN)");
     
     if(n !=a->N()) {
        cerr << " length on the array " << a->N() << " != " << n << " length in file " << endl;
        ExecError("Fatal Error: incompatible length in read array (Op_ReadKN)");
       assert(n==a->N());
       }
      while (f->get(c) && (c!='\n' && c!='\r' ) ) ((void) 0); // eat until control (new line

       // buffer problem if reading value are out range
       for (int i=0;i<n;i++) {
           *f >> (*a)[i] ;
           if (!f->good()) {
               cerr << " reading problem with value i:" << i << " .... " <<(*a)[i]<<endl;
               ExecError("Fatal Error: file  not good in read array (Op_ReadKN)");
               cout << "test i: " << i << " -> " << (*a)[i] << " " << f->good() << endl;
           }
        }
     return f;
   }
};
template<class A>
struct Op_ReadKNM : public binary_function<istream*,KNM<A>*,istream*> {
    static istream *  f(istream  * const  & f,KNM<A>* const  &  a)
    {
        if( !f || !*f) ExecError("Fatal Error: file not open in read array (Op_ReadKNM)");
        int n,m;char c;
        *f >> n >> m;
        if(!f->good()) ExecError("Fatal Error: file  not good in read array (Op_ReadKNM)");
        
        if(n !=a->N() || m != a->M()    ) {
            cerr << " length on the array  N " << a->N() << " != " << n << " n in file " << endl;
            cerr << "  or                  M " << a->M() << " != " << m << " m in file " << endl;
            ExecError("Fatal Error: incompatible length in read array (Op_ReadKNM)");
            assert(n==a->N());
        }
        while (f->get(c) &&  (c!='\n' && c!='\r' ) ) ((void) 0); // eat until control (new line
        
        for (int i=0;i<n;i++)
            for (int j=0;j<m;j++) {
                *f >> (*a)(i,j) ;
                if (!f->good()) {
                    cerr << " reading problem with value i:" << i << " .... " <<(*a)[i]<<endl;
                    ExecError("Fatal Error: file  not good in read array (Op_ReadKN)");
                    cout << "test (i,j): (" << i <<","<<j<<")"<< " -> " << (*a)(i,j) << " " << f->good() << endl;
                }
            }
        return f;
    }
};


template<class A>
struct Op_WriteKNM : public binary_function<ostream*,KNM<A>*,ostream*> {
    static ostream *  f(ostream  * const  & f,KNM<A>* const  &  a) {
        *f << *a;
        return f;
    }
};


template<class A>
struct Op_WriteKN : public binary_function<ostream*,KN<A>*,ostream*> {
    static ostream *  f(ostream  * const  & f,KN<A>* const  &  a) {
         *f << *a;
        return f;
    }
};

/*  FH je supprime  mauvais code nov 2019
template<>
struct Op_WriteKNM<double> : public binary_function<ostream*,KNM<double>*,ostream*> {
    static ostream *  f(ostream  * const  & f,KNM<double>* const  &  a)
    {
       int n =a->N(), m = a->M();
        
        for (int i=0;i<n;i++)
            for (int j=0;j<m;j++)
                if(abs((*a)(i,j))>1.e300 || abs((*a)(i,j))<1.e-300)
                    (*a)(i,j)=0.;
        *f << *a ;
        return f;
    }
};

template<>
struct Op_WriteKN<double> : public binary_function<ostream*,KN<double>*,ostream*> {
    static ostream *  f(ostream  * const  & f,KN<double>* const  &  a)
    {
        int n =a->N();
        for (int i=0;i<n;i++)
            if( (abs((*a)[i])>1.e300) || (abs((*a)[i])<1.e-300) )
                (*a)[i]=0.;
        *f << *a ;
        return f;
    }
};
*/
template<class A>
struct Print: public binary_function<ostream*,A,ostream*> {
  static ostream* f(ostream* const  & a,const A & b)  { *a << b;  return a;}
};

//  ---------------------------------------------
template<class A>
struct set_eq: public binary_function<A*,A,A*> {
  static A* f(A* const  & a,const A & b)  { *a = b; return a;}
};

template<class A,class B>
struct set_eqq: public binary_function<A,B,A> {
  static A f(const A & a,const B & b)  {A aa(a); aa = b; return aa;}
};

template<class A,class B>
struct set_eqq_add: public binary_function<A,B,A> {
  static A f(const A & a,const B & b)  {A aa(a); aa += b; return aa;}
};

template<class A>
struct set_eq_add: public binary_function<A*,A,A*> {
  static A* f(A* const  & a,const A & b)  { *a += b; return a;}
};

template<class A>
struct set_eq_sub: public binary_function<A*,A,A*> {
  static A* f(A* const  & a,const A & b)  { *a -= b; return a;}
};

template<class A>
struct set_eq_mul: public binary_function<A*,A,A*> {
  static A* f(A* const  & a,const A & b)  { *a *= b; return a;}
};

template<class A>
struct set_eq_div: public binary_function<A*,A,A*> {
  static A* f(A* const  & a,const A & b)  { *a /= b; return a;}
};


struct set_peqstring: public binary_function<string**,string*,string**> {
  static string** f(string** const  & a, string * const & b)  {
    if(*a != b )
    { 
	//cerr << " set_peq " << *a << endl;
	freestring(*a);
	//cerr << " set_peq " << *a << " " << " = " << b << " " <<  endl;
	*a = newstring(*b); //(stack ptr) FH mars 2006
	//cerr << " set_peq " << *a << " " << **a << " = " << *b << " " << b <<  endl;
    }
     return a;}
};

//  ---------------------------------------------
template<class A,class B>
struct set_eqarrayp: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {  *a = *b; return a;}
};

template<class A,class B>
struct set_eqarrayp_add: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  { assert(SameShape(*a,*b)); *a += *b; return a;}
};

template<class A,class B>
struct set_eqarrayp_sub: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  { assert(SameShape(*a,*b)); *a -= *b; return a;}
};

template<class A,class B>
struct set_eqarrayp_mul: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  { assert(SameShape(*a,*b)); *a *= *b; return a;}
};


template<class A,class B>
struct set_eqarrayp_div: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  { assert(SameShape(*a,*b)); *a /= *b; return a;}
};

//  ---------------------------------------------
template<class A,class B>
struct set_eqarraypd: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a = *b;
   delete b;
   return a;}
};


template<class A,class B>
struct set_eqarraypd_add: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a += *b; 
  delete b;
  return a;}
};

template<class A,class B>
struct set_eqarraypd_sub: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a -= *b; 
  delete b;
  return a;}
};

template<class A,class B>
struct set_eqarraypd_mul: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a *= *b;
   delete b;
   return a;}
};

template<class A,class B>
struct set_eqarraypd_div: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a /= *b;
   delete b;
   return a;}
};

//  ---------------------------------------------
template<class A,class B>
struct set_eqarray: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  
  {  *a = b;
     return a;}
};

//  --------------------------------------------- august 2007 FH 
template<class A,class B>
struct init_eqarray: public binary_function<A*,B,A*> {
    static A* f(A* const  & a, B const & b)  
    {  a->init(); *a = b;
	return a;}
};
/*
template<class A,class B>
struct init_eqarray_call: public binary_function<A*,B,A*> {
    static A* f(A* const  & a, B const & b)
    {  a->init();
       b.call(*a);
        return a;}
};*/
//  ---------------------------------------------
template<class A,class B>
struct init_eqarraypd: public binary_function<A*,B,A*> {
    static A* f(A* const  & a, B const & b)  {a->init();  *a = *b;
    delete b;
    return a;}
};

//  ---------------------------------------------
template<class A,class B>
struct init_eqarrayp: public binary_function<A*,B,A*> {
    static A* f(A* const  & a, B const & b)  {  a->init(); *a = *b; return a;}
};

// ----------------------------------------  fin modif august 2007

template<class A,class B>
struct set_eqarray_add: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,b));  *a += b; return a;}
};

template<class A,class B>
struct set_eqarray_sub: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,b));  *a -= b; return a;}
};

template<class A,class B>
struct set_eqarray_mul: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,b));  *a *= b; return a;}
};

template<class A,class B>
struct set_eqarray_div: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,b));  *a /= b; return a;}
};

template<class A,class B>
struct set_eq_array: public binary_function<A,B,A> {
  static A f(const A & a, B const & b)  
  {  A aa=a;aa = b;
     return a;}
};

/*
template<class A,class B>
struct set_eq_array_call: public binary_function<A,B,A> {
    static A f(const A & a, B const & b)
    {  A aa=a;
        b.call(aa);//A aa=a;aa = b;
        return a;}
};
 */
template<class A,class B>
struct set_eq_array_add: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,b));  A aa(a);  aa += b; return a;}
};

template<class A,class B>
struct set_eq_array_sub: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,b));  A aa(a);  aa -= b; return a;}
};

template<class A,class B>
struct set_eq_array_mul: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,b));  A aa(a);  aa *= b; return a;}
};

template<class A,class B>
struct set_eq_array_div: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,b));  A aa(a);  aa /= b; return a;}
};
template<class A,class B>
struct set_eq_arrayp: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {  A aa(a);  aa = *b; return a;}
};
//  ---------------------------------------------
template<class A,class B>
struct set_eq_arraypd: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b));A aa(a);  aa = *b;
  delete b;
  return a;}
};

template<class A,class B>
struct set_eq_arrayp_sub: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  { assert(SameShape(a,*b));  A aa(a);  aa -= *b; return a;}
};

template<class A,class B>
struct set_eq_arrayp_mul: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  { assert(SameShape(a,*b));  A aa(a);  aa *= *b; return a;}
};

template<class A,class B>
struct set_eq_arrayp_div: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  { assert(SameShape(a,*b));  A aa(a);  aa /= *b; return a;}
};

template<class A,class B>
struct set_eq_arrayp_add: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  { assert(SameShape(a,*b));  A aa(a);  aa += *b; return a;}
};
template<class A,class B>
struct set_eq_arraypd_add: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa += *b;
   delete b;
   return a;}
};
template<class A,class B>
struct set_eq_arraypd_sub: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa -= *b;
   delete b;
   return a;}
};

template<class A,class B>
struct set_eq_arraypd_mul: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa *= *b;
   delete b;
   return a;}
};

template<class A,class B>
struct set_eq_arraypd_div: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa /= *b;
   delete b;
   return a;}
};

template<class A>
struct PrintP: public binary_function<ostream*,A,ostream*> {
    static ostream* f(ostream* const  & a,const A & b)  {  if(b) *a << *b; 
	// a->flush();// ADD FH MAi 2010 to empty the buffer  baf idea  add flush of ostream 
  //delete b; mars 2006 FH 
   return a;}
};
template<class A>
struct PrintPnd: public binary_function<ostream*,A,ostream*> {
  static ostream* f(ostream* const  & a,const A & b)  
    { if(verbosity>9999) cout << "PrintPnd:  " << b << endl;  *a << *b; return a;}
};




template<class R,class A>  R * set_eqP(R* a,A b){ 
   if (*a != b) delete (*a) ;
  ( *a =b); return a;}
template<class R,class A>  R * set_eqdestroy(R* a,A b){ 
   if (*a != b)  (**a).destroy() ;//  le cas debile Th=Th doit marcher
   // cout << " set_eqdestroy " << a << " " << b << endl;
  ( *a =b); return a;}
  
template<class R,class A>  R * set_eqdestroy_incr(R* a,A b){
  if(b) (*b).increment() ;
  if(*a) (**a).destroy() ;//  le cas debile Th=Th doit marcher
   // cout << " set_eqdestroy " << a << " " << b << endl;
  ( *a =b); return a;}
  
template<class R>  R * set_copy( R* const & a,const R & b){ 
 SHOWVERB( cout << " set_copy " << typeid(R).name() << " " << &b << endl);
  memcpy(a,&b,sizeof(R)); return a;}

template<class R>  R ** set_copy_new( R** const & a,const R * & b){ 
    SHOWVERB( cout << " set_copy_new " << typeid(R).name() << " " << &b << endl);
    *a= new R(*b);
   return a;}

template<class R>  R * set_copyp( R* const & a,const R & b){ 
 SHOWVERB( cout << " set_copy " << typeid(R).name() << " " << &b << endl);
  // memcpy(a,&b,sizeof(R));
   *a = b;
   return a;}

template<class R>  R ** set_copyp_new( R**  a,R*  b){ 
 SHOWVERB( cout << " set_copy " << typeid(R).name() << " " << &b << endl);
  //memcpy(a,&b,sizeof(R)); return a;  FH 2007
   // cerr << " set_copyp_new " << typeid(R).name() << " " << b <<  " " << *b ;
  *a = new R(*b);

  //  cerr << "  -> " << *a << endl; 
  return a;
  }

template<class R>  R ** set_copy_incr( R** const & a, R * const & b){ 
   *a=b;
    SHOWVERB( cout << "set_copy_incr  " << b << " dans  "<< a << endl);
   if(b) b->increment();
   return a;}

template<class R,class A>  R * set_init2( R* const & a,const A & b,const A & c){ 
 SHOWVERB( cout << " set_init2 " << typeid(R).name() << " " << &b << " " << &c << endl);
  a->init(b,c); return a;}
template<class R,class A>  R * set_init( R* const & a,const A & b){ 
 SHOWVERB( cout << " set_init " << typeid(R).name() << " " << &b << endl);
  a->init(b); return a;}
template<class R,class A>  R * set_init_N( R* const & a,const A & b){
    SHOWVERB( cout << " set_init_N" << typeid(R).name() << " " << &b << endl);
    a->init(b.N()); *a=b; return a;}

template<class R,class A>  R * set_initp( R* const & a,const A & b){ 
    SHOWVERB( cout << " set_init " << typeid(R).name() << " " << &b << endl);
    a->init(*b); return a;}

template<class R,class A=R,class B=A> 
struct Op2_add0: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return (a + b);} };  

template<class R,class A=R,class B=A> 
struct Op2_build: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return R(a,b);} };  
template<class R,class A=R,class B=A> 

struct Op2_pbuild: public binary_function<A,B,R*> { 
  static R *f(const A & a,const B & b)  { return new R(a,b);} };  
  
template<class R,class A=R,class B=A> 
struct Op2_add__n: public binary_function<A,B,R*> { 
  static R * f(const A & a,const B & b)  { return new R(a + b);} };    
template<class R,class A=R,class B=A> 
struct Op2_addp_n: public binary_function<A,B,R*> { 
  static R* f(const A & a,const B & b)  { return new R(*a + b);} };   
template<class R,class A=R,class B=A> 
struct Op2_add_pn: public binary_function<A,B,R*> { 
  static R* f(const A & a,const B & b)  { return new R(a + *b);} };   


template<class R,class A=R,class B=A> 
struct Op2_sub0: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return (a - b);} }; 
  
template<class R,class A=R> 
struct Op1_subp: public unary_function<A,R> { 
  static R f(const A & a)  { return (- *a );} };   
template<class R,class A=R> 
struct Op1_sub: public unary_function<A,R> { 
static R f(const A & a)  { return (- a );} };   

template<class R,class A=R,class B=A> 
struct Op2_mulcp: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return (a * *b);} }; 

template<class R,class A=R,class B=A> 
struct Op2_mulc: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return (a * b);} };

template<class R,class A=R,class B=A>
struct Op2_divc: public binary_function<A,B,R> {
    static R f(const A & a,const B & b)  { return (a / b);} };

template<class R,class A=R,class B=A> 
struct Op2_mulpc: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return (b * *a);} }; 

template<class R,class A=R,class B=A> 
struct Op2_mulpcp: public binary_function<A,B,R> {
  static R f(const A & a,const B & b)  { return (*a * *b);} };

template<class R,class A=R,class B=A>
struct Op2_2p_: public binary_function<A,B,R> {
    static R f(const A & a,const B & b)  { return R(*a,*b);} };


template<class R,class A=R,class B=A> 
struct Op2_sub__n: public binary_function<A,B,R*> { 
  static R * f(const A & a,const B & b)  { return new R(a - b);} };    
template<class R,class A=R,class B=A> 
struct Op2_subp_n: public binary_function<A,B,R*> { 
  static R* f(const A & a,const B & b)  { return new R(*a - b);} };   
template<class R,class A=R,class B=A> 
struct Op2_sub_pn: public binary_function<A,B,R*> { 
  static R* f(const A & a,const B & b)  { return new R(a - *b);} };   


template<class R,class A=R,class B=A,class C=A> 
struct Op3_p: public ternary_function<A,B,C,R*> { 
  static R* f(Stack s,const A & a,const B & b,const  C & c )  { return new R(a,b,c);} };   

template<class R,class A=R,class B=A> 
struct Op2_p: public binary_function<A,B,R*> { 
  static R* f(const A & a,const B & b)  { return new R(a,b);} };   




template<class T>
class Transpose{ public:
  T  t;
  Transpose( T  v)
   : t(v) {}
  template<class TT> Transpose( TT  v) : t(v) {}  
  template<class TT> Transpose( TT * v) : t(*v) {}  
  operator const T & () const {return t;}
};
#endif
