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
struct Op1_neg {
  using argument_type  = A;
  using result_type    = R;
  static R f(const A & a)  { return - (R)a;} }; 
  
template<class R,class A=R> 
struct Op1_plus{
  using argument_type  = A;
  using result_type    = R;
  static R f(const A & a)  { return + (R)a;} }; 
  
template<class A> 
struct Op1_not{
  using argument_type  = A;
  using result_type    = bool;
  static bool f(const A & a)  { return ! (bool)a;} }; 
  
template<class R,class A=R,class B=A> 
struct Op2_add {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return ((R)a + (R)b);} }; 

template<class R,class A=R,class B=A> 
struct Op2_sub {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return ((R)a - (R)b);} }; 

template<class R,class A=R,class B=A>
struct Op2_DotDiv {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return DotDiv((R)a, (R)b);} };

template<class R,class A=R,class B=A>
struct Op2_DotStar {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return DotStar((R)a, (R)b);} };


template<class R,class A=R,class B=A>
struct Op2_mul {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  {
  // cout << a << " * " << b <<" => "  << ((R)a * (R)b) << endl;
  return ((R)a * (R)b);} }; 
template<class R,class A,class B>
struct Op2_mull {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  {
  // cout << a << " * " << b <<" => "  << ((R)a * (R)b) << endl;
  return (a * b);} };

template<class R,class A=R,class B=A>
struct Op2_div {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  {
     if (b == B())
       {cerr <<  a << "/" << b << " : " <<  typeid(A).name()  << " " << typeid(B).name() 
             << " " << typeid(R).name() << endl;ExecError(" Div by 0");}
     return ((R)a / (R)b);} };

template<class R,class A,class B >
struct Op2_divv {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  {
     if (b == B())
       {cerr <<  a << "/" << b << " : " <<  typeid(A).name()  << " " << typeid(B).name()
             << " " << typeid(R).name() << endl;ExecError(" Div by 0");}
     return (a / b);} };

template<class R>
struct Op2_pipe {
  using first_argument_type  = R;
  using second_argument_type = R;
  using result_type          = R;
    static R f(const R & a,const R & b)  {   return (a | b);} };

template<class R,class A=R,class B=A> 
struct Op2_mod {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return ((R)a % (R)b);} }; 

template<class A,class B=A> 
struct Op2_lt {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { 
//    cout << a << " < " << b << " = " << ( a<b) << endl;
    return a  < b;} }; 

template<class A,class B=A> 
struct Op2_le {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { return a  <= b;} }; 


template<class A,class B=A> 
struct Op2_gt {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { return a  > b;} }; 


template<class A,class B=A> 
struct Op2_ge {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { return a  >= b;} }; 

template<class A,class B=A> 
struct Op2_eq {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b) 
   { //cout << a << " == " << b << " => " <<( a  == b) << endl;
   return a  == b;} }; 

template<class A,class B=A> 
struct Op2_ne {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { return a  != b;} }; 

struct Op2_and {
  using first_argument_type  = bool;
  using second_argument_type = bool;
  using result_type          = bool;
  static bool f(const bool & a,const bool & b)  { return a  && b;} }; 
  
struct Op2_or {
  using first_argument_type  = bool;
  using second_argument_type = bool;
  using result_type          = bool;
  static bool f(const bool & a,const bool & b)  { return a  || b;} }; 


template<class R,class A,class B> 
struct Op2_padd {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R*;
  static R * f(Stack s,const A & a,const B & b)  { 
   R* r= Add2StackOfPtr2Free(s, a || b ? new R ((a ? *a : nullptr) + (b ? *b : nullptr)) : nullptr);
  // delete a,delete b;
  return r;} }; 


template<class A,class B=A> 
struct Op2_plt {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { bool r= *a  < *b;
  //delete a,delete b;
  return r;} }; 

template<class A,class B=A> 
struct Op2_ple {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { bool r= *a  <= *b;
  // delete a,delete b;
  return r;} }; 


template<class A,class B=A> 
struct Op2_pgt {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { bool r= *a  > *b;
 // delete a,delete b;
  return r;} }; 


template<class A,class B=A> 
struct Op2_pge {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { bool r= *a  >= *b;
  // delete a,delete b;
  return r;} }; 

template<class A,class B=A> 
struct Op2_peq {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  { bool r= *a  == *b;
 //  delete a,delete b;
  return r;} }; 

template<class A,class B=A> 
struct Op2_pne {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = bool;
  static bool f(const A & a,const B & b)  {  bool r=*a  != *b;
  // delete a,delete b;
  return r;} }; 


template<class R,class A=R,class B=A>
struct Op2_pow {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return R(pow(a,b));}};
  


template<class A>
struct Op_Read {
  using first_argument_type  = istream*;
  using second_argument_type = A*;
  using result_type          = istream*;
  static istream *  f(istream  * const  & f,A  * const  &  a)  
   {
       if( !f || !*f) ExecError("Fatal Error: file not open in read value (Op_Read)");
       *f >> *a;
       if(!f->good()) ExecError("Fatal Error: file  not good in read array (Op_ReadKN)");
     return f;
   }
};

template<class A>
struct Op_ReadP {
  using first_argument_type  = istream*;
  using second_argument_type = A**;
  using result_type          = istream*;
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
struct Op_ReadKN {
  using first_argument_type  = istream*;
  using second_argument_type = KN<A>*;
  using result_type          = istream*;
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
struct Op_ReadKNM {
  using first_argument_type  = istream*;
  using second_argument_type = KNM<A>*;
  using result_type          = istream*;
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
struct Op_WriteKNM {
  using first_argument_type  = ostream*;
  using second_argument_type = KNM<A>*;
  using result_type          = ostream*;
    static ostream *  f(ostream  * const  & f,KNM<A>* const  &  a) {
        *f << *a;
        return f;
    }
};


template<class A>
struct Op_WriteKN {
  using first_argument_type  = ostream*;
  using second_argument_type = KN<A>*;
  using result_type          = ostream*;
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
struct Print {
  using first_argument_type  = ostream*;
  using second_argument_type = A;
  using result_type          = ostream*;
  static ostream* f(ostream* const  & a,const A & b)  { *a << b;  return a;}
};

//  ---------------------------------------------
template<class A>
struct set_eq {
  using first_argument_type  = A*;
  using second_argument_type = A;
  using result_type          = A*;
  static A* f(A* const  & a,const A & b)  { *a = b; return a;}
};

template<class A,class B>
struct set_eqq {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(const A & a,const B & b)  {A aa(a); aa = b; return aa;}
};

template<class A,class B>
struct set_eqq_add {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(const A & a,const B & b)  {A aa(a); aa += b; return aa;}
};

template<class A>
struct set_eq_add {
  using first_argument_type  = A*;
  using second_argument_type = A;
  using result_type          = A*;
  static A* f(A* const  & a,const A & b)  { *a += b; return a;}
};

template<class A>
struct set_eq_sub {
  using first_argument_type  = A*;
  using second_argument_type = A;
  using result_type          = A*;
  static A* f(A* const  & a,const A & b)  { *a -= b; return a;}
};

template<class A,class B=A>
struct set_eq_mul {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a,const B & b)  { *a *= b; return a;}
};

template<class A,class B=A>
struct set_eq_div {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a,const B & b)  { *a /= b; return a;}
};


struct set_peqstring {
  using first_argument_type  = string**;
  using second_argument_type = string*;
  using result_type          = string**;
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
struct set_eqarrayp {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {  *a = *b; return a;}
};

template<class A,class B>
struct set_eqarrayp_add {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  { assert(SameShape(*a,*b)); *a += *b; return a;}
};

template<class A,class B>
struct set_eqarrayp_sub {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  { assert(SameShape(*a,*b)); *a -= *b; return a;}
};

template<class A,class B>
struct set_eqarrayp_mul {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  { assert(SameShape(*a,*b)); *a *= *b; return a;}
};


template<class A,class B>
struct set_eqarrayp_div {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  { assert(SameShape(*a,*b)); *a /= *b; return a;}
};

//  ---------------------------------------------
template<class A,class B>
struct set_eqarraypd {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a = *b;
   delete b;
   return a;}
};


template<class A,class B>
struct set_eqarraypd_add {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a += *b; 
  delete b;
  return a;}
};

template<class A,class B>
struct set_eqarraypd_sub {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a -= *b; 
  delete b;
  return a;}
};

template<class A,class B>
struct set_eqarraypd_mul {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a *= *b;
   delete b;
   return a;}
};

template<class A,class B>
struct set_eqarraypd_div {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a /= *b;
   delete b;
   return a;}
};

//  ---------------------------------------------
template<class A,class B>
struct set_eqarray {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  
  {  *a = b;
     return a;}
};

//  --------------------------------------------- august 2007 FH 
template<class A,class B>
struct init_eqarray {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
    static A* f(A* const  & a, B const & b)  
    {  a->init(); *a = b;
	return a;}
};
/*
template<class A,class B>
struct init_eqarray_call {
    static A* f(A* const  & a, B const & b)
    {  a->init();
       b.call(*a);
        return a;}
};*/
//  ---------------------------------------------
template<class A,class B>
struct init_eqarraypd {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
    static A* f(A* const  & a, B const & b)  {a->init();  *a = *b;
    delete b;
    return a;}
};

//  ---------------------------------------------
template<class A,class B>
struct init_eqarrayp {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
    static A* f(A* const  & a, B const & b)  {  a->init(); *a = *b; return a;}
};

// ----------------------------------------  fin modif august 2007

template<class A,class B>
struct set_eqarray_add {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,b));  *a += b; return a;}
};

template<class A,class B>
struct set_eqarray_sub {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,b));  *a -= b; return a;}
};

template<class A,class B>
struct set_eqarray_mul {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,b));  *a *= b; return a;}
};

template<class A,class B>
struct set_eqarray_div {
  using first_argument_type  = A*;
  using second_argument_type = B;
  using result_type          = A*;
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,b));  *a /= b; return a;}
};

template<class A,class B>
struct set_eq_array{
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
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
struct set_eq_array_add  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {assert(SameShape(a,b));  A aa(a);  aa += b; return a;}
};

template<class A,class B>
struct set_eq_array_sub  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {assert(SameShape(a,b));  A aa(a);  aa -= b; return a;}
};

template<class A,class B>
struct set_eq_array_mul  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {assert(SameShape(a,b));  A aa(a);  aa *= b; return a;}
};

template<class A,class B>
struct set_eq_array_div  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {assert(SameShape(a,b));  A aa(a);  aa /= b; return a;}
};
template<class A,class B>
struct set_eq_arrayp  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {  A aa(a);  aa = *b; return a;}
};
//  ---------------------------------------------
template<class A,class B>
struct set_eq_arraypd  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b));A aa(a);  aa = *b;
  delete b;
  return a;}
};

template<class A,class B>
struct set_eq_arrayp_sub  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  { assert(SameShape(a,*b));  A aa(a);  aa -= *b; return a;}
};

template<class A,class B>
struct set_eq_arrayp_mul  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  { assert(SameShape(a,*b));  A aa(a);  aa *= *b; return a;}
};

template<class A,class B>
struct set_eq_arrayp_div  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  { assert(SameShape(a,*b));  A aa(a);  aa /= *b; return a;}
};

template<class A,class B>
struct set_eq_arrayp_add  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  { assert(SameShape(a,*b));  A aa(a);  aa += *b; return a;}
};
template<class A,class B>
struct set_eq_arraypd_add  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa += *b;
   delete b;
   return a;}
};
template<class A,class B>
struct set_eq_arraypd_sub  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa -= *b;
   delete b;
   return a;}
};

template<class A,class B>
struct set_eq_arraypd_mul  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa *= *b;
   delete b;
   return a;}
};

template<class A,class B>
struct set_eq_arraypd_div  {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = A;
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa /= *b;
   delete b;
   return a;}
};

template<class A>
struct PrintP {
  using first_argument_type  = ostream*;
  using second_argument_type = A;
  using result_type          = ostream*;
    static ostream* f(ostream* const  & a,const A & b)  {  if(b) *a << *b; 
	// a->flush();// ADD FH MAi 2010 to empty the buffer  baf idea  add flush of ostream 
  //delete b; mars 2006 FH 
   return a;}
};
template<class A>
struct PrintPnd {
  using first_argument_type  = ostream*;
  using second_argument_type = A;
  using result_type          = ostream*;
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
struct Op2_add0 {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return (a + b);} };  

template<class R,class A=R,class B=A> 
struct Op2_build {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return R(a,b);} };  
template<class R,class A=R,class B=A> 

struct Op2_pbuild {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R*;
  static R *f(const A & a,const B & b)  { return new R(a,b);} };  
  
template<class R,class A=R,class B=A> 
struct Op2_add__n {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R*;
  static R * f(const A & a,const B & b)  { return new R(a + b);} };    
template<class R,class A=R,class B=A> 
struct Op2_addp_n {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R*;
  static R* f(const A & a,const B & b)  { return new R(*a + b);} };   
template<class R,class A=R,class B=A> 
struct Op2_add_pn {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R*;
  static R* f(const A & a,const B & b)  { return new R(a + *b);} };   


template<class R,class A=R,class B=A> 
struct Op2_sub0 {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return (a - b);} }; 
  
template<class R,class A=R> 
struct Op1_subp{
  using argument_type  = A;
  using result_type    = R;
  static R f(const A & a)  { return (- *a );} };   
template<class R,class A=R> 
struct Op1_sub{
  using argument_type  = A;
  using result_type    = R;
static R f(const A & a)  { return (- a );} };   

template<class R,class A=R,class B=A> 
struct Op2_mulcp {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return (a * *b);} }; 

template<class R,class A=R,class B=A> 
struct Op2_mulc {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return (a * b);} };

template<class R,class A=R,class B=A>
struct Op2_divc {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
    static R f(const A & a,const B & b)  { return (a / b);} };

template<class R,class A=R,class B=A> 
struct Op2_mulpc {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return (b * *a);} }; 

template<class R,class A=R,class B=A> 
struct Op2_mulpcp {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
  static R f(const A & a,const B & b)  { return (*a * *b);} };

template<class R,class A=R,class B=A>
struct Op2_2p_ {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R;
    static R f(const A & a,const B & b)  { return R(*a,*b);} };


template<class R,class A=R,class B=A> 
struct Op2_sub__n {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R*;
  static R * f(const A & a,const B & b)  { return new R(a - b);} };    
template<class R,class A=R,class B=A> 
struct Op2_subp_n {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R*;
  static R* f(const A & a,const B & b)  { return new R(*a - b);} };   
template<class R,class A=R,class B=A> 
struct Op2_sub_pn {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R*;
  static R* f(const A & a,const B & b)  { return new R(a - *b);} };   


template<class R,class A=R,class B=A,class C=A> 
struct Op3_p: public ternary_function<A,B,C,R*> { 
  static R* f(Stack s,const A & a,const B & b,const  C & c )  { return new R(a,b,c);} };   

template<class R,class A=R,class B=A> 
struct Op2_p {
  using first_argument_type  = A;
  using second_argument_type = B;
  using result_type          = R*;
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
