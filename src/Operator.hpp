
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
  static R * f(const A & a,const B & b)  { R* r= new R (*a  + *b);delete a,delete b;return r;} }; 


template<class A,class B=A> 
struct Op2_plt: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  < *b;delete a,delete b;return r;} }; 

template<class A,class B=A> 
struct Op2_ple: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  <= *b;delete a,delete b;return r;} }; 


template<class A,class B=A> 
struct Op2_pgt: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  > *b;delete a,delete b;return r;} }; 


template<class A,class B=A> 
struct Op2_pge: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  >= *b;delete a,delete b;return r;} }; 

template<class A,class B=A> 
struct Op2_peq: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  { bool r= *a  == *b;delete a,delete b;return r;} }; 

template<class A,class B=A> 
struct Op2_pne: public binary_function<A,B,bool> { 
  static bool f(const A & a,const B & b)  {  bool r=*a  != *b;delete a,delete b;return r;} }; 


template<class R,class A=R,class B=A>
struct Op2_pow: public binary_function<A,B,R> {
  static R f(const A & a,const B & b)  { return R(pow(a,b));}};
  


template<class A>
struct Op_Read : public binary_function<istream*,A*,istream*> {
  static istream *  f(istream  * const  & f,A  * const  &  a)  
   { 
     *f >> *a;
     return f;
   }
};

template<class A>
struct Op_ReadKN : public binary_function<istream*,KN<A>*,istream*> {
  static istream *  f(istream  * const  & f,KN<A>* const  &  a)  
   { 
     if( !f || !*f) ExecError("Fatal Error file: not open in read array (Op_ReadKN)");
     int n;char c;
     *f >> n;
     if(!f->good()) ExecError("Fatal Error file  not good   in read array (Op_ReadKN)");
     
     if(n !=a->N()) {
        cerr << " length on the array " << a->N() << " != " << n << " length in file " << endl;
        ExecError("Fatal Error incompatible lengh in read array (Op_ReadKN)");
       assert(n==a->N());
       }
     while (f->get(c) &&  (c!='\n' && c!='\r' ) ) 0; // eat until control (new line

     for (int i=0;i<n;i++)
       *f >> (*a)[i] ;
     return f;
   }
};



template<class A>
struct Print: public binary_function<ostream*,A,ostream*> {
  static ostream* f(ostream* const  & a,const A & b)  { *a << b;  return a;}
};
//  ---------------------------------------------
template<class A>
struct set_eq: public binary_function<A*,A,A*> {
  static A* f(A* const  & a,const A & b)  { *a = b; return a;}
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

template<class A>
struct set_peq: public binary_function<A*,A,A*> {
  static A* f(A* const  & a,const A & b)  { delete *a; *a = b; return a;}
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
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a = *b; delete b;return a;}
};

template<class A,class B>
struct set_eqarraypd_add: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a += *b; delete b;return a;}
};

template<class A,class B>
struct set_eqarraypd_sub: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a -= *b; delete b;return a;}
};

template<class A,class B>
struct set_eqarraypd_mul: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a *= *b; delete b;return a;}
};

template<class A,class B>
struct set_eqarraypd_div: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  {assert(SameShape(*a,*b));  *a /= *b; delete b;return a;}
};

//  ---------------------------------------------
template<class A,class B>
struct set_eqarray: public binary_function<A*,B,A*> {
  static A* f(A* const  & a, B const & b)  
  {  *a = b;
     return a;}
};



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
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b));A aa(a);  aa = *b; delete b;return a;}
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
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa += *b; delete b;return a;}
};
template<class A,class B>
struct set_eq_arraypd_sub: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa -= *b; delete b;return a;}
};

template<class A,class B>
struct set_eq_arraypd_mul: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa *= *b; delete b;return a;}
};

template<class A,class B>
struct set_eq_arraypd_div: public binary_function<A,B,A> {
  static A f(A const  & a, B const & b)  {assert(SameShape(a,*b)); A aa(a);  aa /= *b; delete b;return a;}
};

template<class A>
struct PrintP: public binary_function<ostream*,A,ostream*> {
  static ostream* f(ostream* const  & a,const A & b)  {  *a << *b;delete b; return a;}
};
template<class A>
struct PrintPnd: public binary_function<ostream*,A,ostream*> {
  static ostream* f(ostream* const  & a,const A & b)  
  {  *a << *b; return a;}
};




template<class R,class A>  R * set_eqP(R* a,A b){ 
   if (*a != b) delete (*a) ;
  ( *a =b); return a;}
template<class R,class A>  R * set_eqdestroy(R* a,A b){ 
   if (*a != b)  (**a).destroy() ;//  le cas debile Th=Th doit marcher
   // cout << " set_eqdestroy " << a << " " << b << endl;
  ( *a =b); return a;}
  
template<class R,class A>  R * set_eqdestroy_incr(R* a,A b){
   (*b).increment() ;
   (**a).destroy() ;//  le cas debile Th=Th doit marcher
   // cout << " set_eqdestroy " << a << " " << b << endl;
  ( *a =b); return a;}
  
template<class R>  R * set_copy( R* const & a,const R & b){ 
 SHOWVERB( cout << " set_copy " << typeid(R).name() << " " << &b << endl);
  memcpy(a,&b,sizeof(R)); return a;}

template<class R>  R ** set_copy_incr( R** const & a, R * const & b){ 
   *a=b;
   b->increment();
   return a;}

template<class R,class A>  R * set_init2( R* const & a,const A & b,const A & c){ 
 SHOWVERB( cout << " set_init2 " << typeid(R).name() << " " << &b << " " << &c << endl);
  a->init(b,c); return a;}
template<class R,class A>  R * set_init( R* const & a,const A & b){ 
 SHOWVERB( cout << " set_init " << typeid(R).name() << " " << &b << endl);
  a->init(b); return a;}

template<class R,class A=R,class B=A> 
struct Op2_addp: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return (*a + *b);} };  

template<class R,class A=R,class B=A> 
struct Op2_dotstarp: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return R(*a,*b);} };  
  
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
struct Op2_subp: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return (*a - *b);} }; 
  
template<class R,class A=R> 
struct Op1_subp: public unary_function<A,R> { 
  static R f(const A & a)  { return (- *a );} };   

template<class R,class A=R,class B=A> 
struct Op2_mulcp: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return (a * *b);} }; 
  
template<class R,class A=R,class B=A> 
struct Op2_mulpc: public binary_function<A,B,R> { 
  static R f(const A & a,const B & b)  { return (b * *a);} }; 

template<class R,class A=R,class B=A> 
struct Op2_sub__n: public binary_function<A,B,R*> { 
  static R * f(const A & a,const B & b)  { return new R(a - b);} };    
template<class R,class A=R,class B=A> 
struct Op2_subp_n: public binary_function<A,B,R*> { 
  static R* f(const A & a,const B & b)  { return new R(*a - b);} };   
template<class R,class A=R,class B=A> 
struct Op2_sub_pn: public binary_function<A,B,R*> { 
  static R* f(const A & a,const B & b)  { return new R(a - *b);} };   











































































































































































  
