//file afonction.h
#ifndef __AFONCTION__
#define __AFONCTION__
#include "showverb.hpp" 


#include <typeinfo>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <cstring>
#include "error.hpp"
#include <map>
#include <deque>
#include <list>
#include <vector>
#include <queue>
#include <complex>
#include <string>
#include <cstdlib>
#include <algorithm>
extern bool showCPU;
#include <time.h> 

inline double CPUtime(){
#ifdef SYSTIMES
  struct tms buf;
  if (times(&buf)!=-1)
    return ((double)buf.tms_utime+(double)buf.tms_stime)/(long) sysconf(_SC_CLK_TCK);
  else
#endif
    return ((double) clock())/CLOCKS_PER_SEC;
}

extern long verbosity;  // level off printing


extern bool  withrgraphique;

//   in the stack we save all the variable 
//   a adresse 0 we have the MeshPointStack to defineP,N, ....
//   a adresse sizeof(void *) 
// 
//
//  Offset in (void *)
const int MeshPointStackOffset =0;
const int ParamPtrOffset = 1;
const int ElemMatPtrOffset = 2;

const int BeginOffset = 3;
//  0 : MeshPoint pointeur 
//  1 : ParamPtrOffset
//  2 : Truc les matrice elementaire


using namespace std;

#include "AnyType.hpp"
#include "String.hpp"


typedef void *Stack;

template<class T>
T * Stack_offset (Stack stack,size_t offset)  {return   (T *) (void *) (((char *) stack)+offset);}

template<class T>
T * & Stack_Ptr (Stack stack,size_t offset)  {return   (T * &)  (((void **) stack)[offset]);}
 void ShowType(ostream & f);

class basicForEachType;
class E_F1_funcT_Type;
class E_F0;  //  une instruction exec time 
class C_F0;  //  une instruction  complie time
class ListOfInst;
class Polymorphic;
class OneOperator;
typedef const E_F0  *  Expression;
class AC_F0;
class basicAC_F0;
typedef complex<double> Complex;

typedef pair<aType,const  E_F0  *>  Type_Expr ;// to store the type and the expression 

 int  FindType(const char * name) ; 
  void lgerror (const char* s) ;  
 void CompileError(string msg="",aType r=0);
 void ExecError(string msg="");
 
struct UnId {
  const char * id;
  aType r;
  Expression  e; 
  deque<UnId> * array; //  to store a array 
  aType re; 
  bool ref; // a ref or non 
  UnId() :id(0),r(0),e(0),array(0),re(0),ref(false) {}
  UnId(const char * idd) :id(idd),r(0),e(0),array(0),re(0),ref(false) {}
  UnId(const char * idd,const C_F0 & ee,aType rr,bool reff) ;  
  UnId(deque<UnId>  * d) : id(0),e(0),array(d),re(0),ref(false) {}
};

typedef deque<UnId> ListOfId;
//  xxx is a type so xxx can't be a parameter 
#define ATYPE(xxx) map_type[typeid(xxx).name()]
/* #define NEW_TYPE(type) map_type[typeid(type).name()] = new ForEachType<type >(0,0)
//#define NEW_TYPE(type) map_type[typeid(type).name()] = new ForEachType<type >()
#define NEW_TYPE_I(type,i,d) map_type[typeid(type).name()] = new ForEachType<type>(i,d)   
#define NEW_TYPE_Ptr(type) map_type[typeid(type*).name()] = new ForEachTypePtr<type  >()
#define NEW_TYPE_PtrND(type) map_type[typeid(type*).name()] = new  ForEachTypePtr<type >(0)
#define NEW_TYPE_PtrNIND(type) map_type[typeid(type*).name()] = new ForEachTypePtr<type >(0,0)
//#define NEW_TYPE_PtrI(type) map_type[typeid(type*).name()] = new ForEachTypePtr<type*>(Initialize<type>)
*/
class CodeAlloc { public:
/*   static size_t nb,lg;   
   static void * mem;
   static void * lgmax;
   void *operator new(size_t ll ) { nb++;lg+=ll; return malloc(ll);} 
   void *operator new[](size_t ll ) {nb++;lg+=ll; return malloc(ll);} 
   void operator delete(void * pp)   {free(pp);}
   void operator delete[](void * pp) { free(pp);} */
};


extern Polymorphic * TheOperators, * TheRightOperators;

//  -------------
extern  C_F0 *pOne,*pZero,*pminusOne;


typedef   AnyType (* Function1)(Stack, const AnyType &);
typedef   AnyType (* Function2)(Stack, const AnyType &,const AnyType &);
typedef   AnyType (* CFunction2)(Stack, const E_F0 *,const E_F0 *);
typedef   AnyType (* CFunction4)(Stack, const E_F0 *,const E_F0 *,const E_F0 *,const E_F0 *);


Expression NewExpression(Function1,Expression);
Expression NewExpression(Function2,Expression,Expression);


inline Type_Expr make_Type_Expr(aType t,const E_F0  * e) {return make_pair(t,e);}
inline Type_Expr make_Type_Expr(const E_F0  * e,aType t) {return make_pair(t,e);}

struct Keyless : binary_function<const char *,const char *, bool>
   { 
    typedef const char * Key;
    bool operator()(const Key& x, const Key& y) const { return strcmp(x,y)<0;} };
    

// un table Iden     
class TableOfIdentifier:CodeAlloc {
  public:
  struct Value;
  typedef const char * Key;
  typedef map<Key,Value,Keyless> maptype;
  typedef pair<const Key,Value> pKV;
  typedef maptype::iterator iterator;
  typedef maptype::const_iterator const_iterator;
  
  struct  Value :public Type_Expr {
    pKV *  next; // link all the variable in reverse order to call delete on each variable 
    bool del; 
    Value(const Type_Expr & vv,pKV * n,bool dd=true) : Type_Expr(vv),next(n),del(dd) {}
    Value(aType t,E_F0  *f,pKV *n,bool dd=true): Type_Expr(t,f),next(n),del(dd) {}
  };//  to store the type and the expression 
  pKV *   listofvar;
  
// struct Keyless : binary_function<Key,Key, bool>
//   { bool operator()(const Key& x, const Key& y) const{ return strcmp(x,y)<0;} };


  maptype m;
  C_F0 Find(Key) const ; 
  C_F0 Find(Key,const basicAC_F0 &) const ; 
  
  const Type_Expr & New(Key k,const Type_Expr &  v,bool del=true);
  void Add(Key k,Key op,OneOperator *p0,OneOperator *p1=0,
      OneOperator *p2=0,OneOperator *p3=0,OneOperator *p4=0,
      OneOperator *p5=0,OneOperator *p6=0)  ;
template<class T>         
  C_F0 NewVar(Key k,aType t,size_t & top,const C_F0 &i) ;
template<class T>         
  C_F0 NewVar(Key k,aType t,size_t & top,const basicAC_F0 &args) ;
template<class T,class U>         
  C_F0 NewVar(Key k,aType t,size_t & top,const basicAC_F0 &args,const U & data) ;
//  C_F0 NewVar(Key k,aType t,size_t & top,const basicAC_F0 &args,const C_F0& f) ;
template<class T>         
  C_F0 NewVar(Key k,aType t,size_t & top) ;
  C_F0 NewID(aType t,Key k, C_F0 & c,size_t & top,bool del=true);   
  C_F0 NewID(aType t,Key k, C_F0 & c,const ListOfId & l,size_t & top,bool del=true);   
template<class T>   
  C_F0 NewFESpace(Key k,aType t,size_t & top,const basicAC_F0 &args);
  friend   ostream & operator<<(ostream & f,const TableOfIdentifier & );
  C_F0 destroy();
  TableOfIdentifier() : listofvar(0) {};
};
//  for all the type of the language 
class basicForEachType : public CodeAlloc {
    const type_info  * ktype;  // the real type_info
  //  const type_info *ktypefunc;// the type of code 
    public:
     const size_t size;

    
    typedef OneOperator * CastFunc;
    typedef map<aType,CastFunc>::const_iterator const_cast_iterator;

    typedef const char * Key;

   // virtual  void print(ostream &f,const void *p) const =0;
                            
    friend ostream & operator<<(ostream & f,const basicForEachType & e) 
      { f << '<' << e.name() << '>' ;return f;}
     void Show(ostream & f) const ;
     const char * name() const  { return this  ?  ktype->name() :"NULL" ;}
     virtual bool CastingFrom(const basicForEachType * t) const ;
     //  modif FH -----  A TESTER  // 
     virtual bool SametypeRight(const basicForEachType * t) const {return  this == t || t == un_ptr_type;}
//     virtual Type_Expr init(const Type_Expr & te) const { return Type_Expr(0,0);}
     virtual int TYPEOFID() const  {return 0;}
//     bool SametypeLeft(const basicForEachType * t) const {return  t == this;}
   //  bool To(const basicForEachType * t) const { throwassert(t && this);return un_ptr_type == this ? t->un_ptr_type == this :  t == this;}
     virtual C_F0 CastTo(const C_F0 & e) const ; 
     virtual void SetArgs(const ListOfId *lid) const ;// { cout << "SetArgs::\n " ;throwassert(lid==0 || lid->size()==0);}
     aType right() const {return un_ptr_type;};
     Expression RightValueExpr(Expression f) const; 
 //    Type_Expr NewVar(Key k,aType t,size_t & top,const C_F0 &i);
     virtual  C_F0 Initialization(const Type_Expr & e) const ;
     virtual Expression Destroy(const C_F0 &) const ;
     virtual bool ExistDestroy() const {return destroy;} 
     virtual Type_Expr SetParam(const C_F0 & c,const ListOfId * l,size_t & top) const;
     // { return make_pair<aType,const E_F0  *>(this,c.left());}

   protected: 
    inline basicForEachType(const type_info  & k ,const size_t ,
                            const E_F1_funcT_Type * p=0,basicForEachType *rr=0,
                            Function1 iv=0,Function1 id=0) ;
/*    inline basicForEachType(const type_info  & k ,const type_info  & kf ,const size_t ,
                            const E_F1_funcT_Type * p=0,basicForEachType *rr=0,
                            Function1 iv=0,Function1 id=0) ;*/

public:    
    const basicForEachType * un_ptr_type;  // type of right exp
   private:
 //   map<aType,CastFunc> mapofcast;
    OneOperator * casting; // list of operator for casting to this type 
    
    const E_F1_funcT_Type * un_ptr;        //  is ptr -> get value function
    
    
    Function1 InitExp;       //  to init the ptr value 
    Function1  destroy;//  the destroy function 
    TableOfIdentifier ti; //  all polymorphisme of the Identifier   
   public:
  // basicForEachType * FunctionType() const;// { return funct_type ? funct_type : (funct_type= new FuncForEachType(this));}
   C_F0  Find(const char * k) const; // {return ti->Find(k);}
   C_F0  Find(const char * k,const basicAC_F0 & args) const; // {return ti->Find(k);}
   void  New(Key k,Type_Expr  v,bool del=true){ti.New(k,v,del);}
  
  void Add(Key k,Key op,OneOperator *p0,OneOperator *p1=0,
      OneOperator *p2=0,OneOperator *p3=0,OneOperator *p4=0,
      OneOperator *p5=0,OneOperator *p6=0)  
     {ti.Add(k,op,p0,p1,p2,p3,p4,p5,p6);}     
 
 	void AddCast(CastFunc f1,CastFunc f2=0,CastFunc f3=0,CastFunc f4=0,
 	             CastFunc f5=0,CastFunc f6=0,CastFunc f7=0,CastFunc f8=0);
    ostream & ShowTable(ostream & f) const { f << ti; return f;}
    
  //  basicForEachType * funct_type;
    
};


template<typename T> 
inline basicForEachType * atype() { 
  basicForEachType * r=map_type[typeid(T).name()];
  if (! r) { cerr << "Error: aType  '" << typeid(T).name() << "', D'ont exist\n";
             ShowType(cerr);
            throw(ErrorExec("exit",1));}
  return r;}


//  --------
//typedef basicForEachType TheType;

//  const basicForEachType * ktype; // compilation time 

//  class for all   exp 
// a left exp is a pointer expression 
//  -------
//  --  exec times le code is just E_F0*(fonction without args)
class C_LF2;
class C_LF1;

//  3 types of function/expression  0,1,2 args  
class E_F0 :public CodeAlloc{ public:

  struct kless : binary_function<Expression,Expression, bool>
   { bool operator()(const Expression& x, const Expression& y) const{ 
     //cout << x << " " << y << x->compare(y) << " ::: ";
      int r1 = x->compare(y);// , r2 = y->compare(x);
     //assert(r1+r2==0);
     return r1<0;} };  
   typedef map<const E_F0 *,int,kless> MapOfE_F0;

    virtual AnyType operator()(Stack)  const =0;
    virtual bool Empty() const {return !this; }
   // virtual E_F0 * destroy(Stack ) const {return 0;}
  //  virtual const E_F0 * Parameter(Stack ) const {return this;}
    virtual size_t nbitem() const {return 1;}
    virtual bool EvaluableWithOutStack() const {return false;} // 
    virtual bool MeshIndependent() const {return true;} // 
    virtual E_F0 * right_E_F0() const { return 0;}
    virtual ~E_F0() {}
    virtual int compare (const E_F0 *t) const { int r= (t==this) ? 0 : ( ( this<t) ?-1 : 1);
     //cout << "cmp " <<  typeid(*this).name() << r << endl; 
     return r;} // to give a order in instuction 
    virtual int Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const;  // build optimisation
    virtual AnyType operator()(Stack stack,AnyType *)  const { return operator()(stack);}  // call optim code
    virtual  operator aType ()  const { assert(0);return 0;}   // the type of the expression
    virtual ostream & dump(ostream &f) const  { f << ' ' << typeid(*this).name() << ' ' << this << ' '  ;return f; }
    // for OPTIMIZATION
    
    int find(const MapOfE_F0 & m) const;
    int insert(Expression  opt,deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const;
 };  
 
inline ostream & operator<<(ostream & f,const E_F0 &e) { if(&e) e.dump(f); else f << " --0-- " ;return f;}
// a  
class E_F0mps : public E_F0 { public:
  virtual bool MeshIndependent() const {return false;} // 
};

class E_F0info : public E_F0 { public:
  // not a real expression just to pass information 
    virtual bool EvaluableWithOutStack() const {return true;} // 
    virtual bool MeshIndependent() const {return true;} // 
    AnyType operator()(Stack s)  const {  
    return SetAny<const E_F0 *>(this);}
    operator aType () const { return atype<Expression>();} 

  
};

class E_F1 : public CodeAlloc{ public: virtual AnyType operator()(Stack,AnyType &)  const =0;}; 
class E_F2 : public CodeAlloc{ public: virtual AnyType operator()(Stack,AnyType &,AnyType &)  const =0;};
class E_FN : public CodeAlloc{ public: virtual AnyType operator()(Stack,size_t N,...)  const =0;};

//   class to play with  polymorphisme 
//   ---------------------------------
class basicAC_F0;
class  ArrayOfaType : public CodeAlloc{ 
  //  class for the type of parameter
   aType tt[4]; 
   protected:
   bool ellipse; 
   int n;
   aType * t; // array of type  
   void operator=(const ArrayOfaType &); // no set  operator
   public:
 //  ArrayOfaType() :n(0),t(0),ellipse(false) {}
   explicit ArrayOfaType(bool ell=false) 
       :n(0),t(0),ellipse(ell) {}
   
   explicit ArrayOfaType(const aType & a,bool ell=false) 
       :n(1),t(tt),ellipse(ell)  {t[0]=a;}
       
   explicit ArrayOfaType(const aType & a,const aType & b,bool ell=false) 
       :n(2),t(tt),ellipse(ell)  {t[0]=a,t[1]=b;}
       
   explicit ArrayOfaType(const aType & a,const aType & b,const aType & c,bool ell=false) 
       :n(3),t(tt),ellipse(ell)  {t[0]=a,t[1]=b;t[2]=c;}
       
   explicit ArrayOfaType(const aType & a,const aType & b,const aType & c,const aType & d,bool ell=false) 
       :n(4),t(tt),ellipse(ell)  {t[0]=a,t[1]=b;t[2]=c;t[3]=d;
       /* cout << * a << *b << * c << * d << " ---------" << endl; */}
       
   ArrayOfaType(const basicAC_F0 & ) ;
   ArrayOfaType(const ArrayOfaType & ); // 
   ArrayOfaType(const ListOfId * l);
   ~ArrayOfaType() { if(t && t != tt) delete [] t;t=0;n=0;}
   bool WithOutCast( const ArrayOfaType & a) const ;  
   bool WithCast( const ArrayOfaType & a,int nbcast=100000) const ;  // return the number of cast 
   // exactly comparaison 
   bool operator==( const ArrayOfaType & a) const { 
     if (a.n != n || a.ellipse !=ellipse) return false;
     for (int i=0;i<n;i++)  
       if (t[i] != a.t[i]) 
         return false; 
     return true;}
   
   friend ostream & operator<<(ostream & f,const ArrayOfaType & a);
};


    
class  OneOperator : public ArrayOfaType {
    friend class MakeVectSpaceN;
    const basicForEachType * r; //  return type 
    OneOperator *next; // to make a list of OneOperator
    public: 
    int pref; //  to try to solve ambiguity for binary operator
    //  10 for bool, 20 for int , 30 for long , 40, for float, 50 double, 60 for complex, 70 string
    //  string+ 1 => string 
    // 1+string => string 
    OneOperator(aType rr) : r(rr),ArrayOfaType(),next(0),pref(0) {throwassert(r);}
    OneOperator(aType rr,aType  a) : r(rr),ArrayOfaType(a,false),next(0),pref(0) {throwassert(rr && a );}
    OneOperator(aType rr,aType  a,aType  b) : r(rr),ArrayOfaType(a,b,false),next(0),pref(0) {
     throwassert(rr && a && b);} 
    OneOperator(aType rr,aType  a,aType  b,aType c) : r(rr),ArrayOfaType(a,b,c,false),next(0),pref(0) {throwassert(rr && a && b && c);} 
    OneOperator(aType rr,aType  a,aType  b,aType c,aType d) : r(rr),ArrayOfaType(a,b,c,d,false),next(0),pref(0) {throwassert(rr && a && b && c);} 
    OneOperator(aType rr,const ArrayOfaType &ta) : r(rr),ArrayOfaType(ta),next(0),pref(0) {throwassert(rr);} 
    OneOperator(aType rr,bool ellipse) : r(rr),ArrayOfaType(ellipse),next(0),pref(0) {throwassert(rr );} 
    OneOperator(aType rr,const ListOfId *l) : r(rr),ArrayOfaType(l),next(0),pref(0) {throwassert(rr );} 
    
    typedef pair<const OneOperator *,int> pair_find;
    void operator+=(OneOperator &a){throwassert(a.next==0);a.next=next;next=&a;} 
    //  a way to make none recurve delete  good   
    virtual ~OneOperator();
    pair_find Find(const ArrayOfaType & at) const ;
    pair_find FindWithOutCast(const ArrayOfaType & at) const ; // for 
    OneOperator * FindSameR(const ArrayOfaType & at)  ; 
       
    void Show(const ArrayOfaType & at,ostream &f=cerr) const;
    void Show(ostream &f=cerr) const;
    operator aType () const { return r;}
    virtual E_F0 * code(const basicAC_F0 &) const =0; 
    const OneOperator * Simple() const { return next||n?0:this;}
    friend ostream & operator<<(ostream & f,const OneOperator & a);
    
};


class Polymorphic:  public E_F0mps {
   //  a list of type 
   //  simple, array or function
   private: 
   typedef const char * Key;
   typedef OneOperator * Value;
 //  struct Keyless : binary_function<Key,Key, bool>
 //  { bool operator()(const Key& x, const Key& y) const{ return strcmp(x,y)<0;} };
   
   typedef map<Key,Value,Keyless> maptype;          //  
   typedef maptype::const_iterator const_iterator;  // 
   typedef maptype::iterator iterator;              // 
   //  remark the map is mutable because 
   //  a expression is const E_F0 *
   // So There is a incompatibility between 
   //   we save an  expression in a variable 
   //   we have to add thing to a polymorphisme expression
   mutable maptype m; //  all polymorphisme of the Identifier
   Expression e; // default expression
  public:    
  Polymorphic() : m(),e(0) {}
  
//  by default Empty and do nothing      
 virtual AnyType operator()(Stack ) const  { return Nothing;}
 virtual bool Empty() const {return true;} //  by default Empty 
 
 const  OneOperator * Find(const char *op, const ArrayOfaType &at) const;
 const  OneOperator * FindWithOutCast(const char *op,const  ArrayOfaType &at) const;
 void Show(const char *op,const ArrayOfaType & at,ostream &f=cerr)const ; 
 void  Add(const char * op,OneOperator * p0  ,OneOperator * p1=0,OneOperator * p2=0,
                           OneOperator * p3=0,OneOperator * p4=0,OneOperator * p5=0,
                           OneOperator * p6=0,OneOperator * p7=0,OneOperator * p8=0,
                           OneOperator * p9=0,OneOperator * pa=0,OneOperator * pb=0,
                           OneOperator * pc=0,OneOperator * pd=0,OneOperator * pe=0
                           ) const
      {Addp(op,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,pa,pb,pc,pd,pe,0);}
 void Add(const char * op,OneOperator ** pp) const ;      
 private:
 void Addp(const char * op,OneOperator * pp,...) const ;
  friend ostream & operator<<(ostream & f,const Polymorphic & a);
  
  
};




//   the type for polymorphisme of id 



//  compile time expression 
class basicAC_F0;
 class C_F0 {
   friend class CC_F0;
  protected: 
  Expression  f; //  the expression code 
  aType r;   // the expression type
  
 public: 
	   //  the constructeur 
	  C_F0() :f(0),r(0) {}
	  C_F0(const C_F0 & c):f(c.f),r(c.r) {}
	  C_F0(const C_F0 & a,const C_F0 & b); // concatenation 
	  C_F0(const Type_Expr & a):r(a.first),f(a.second) {}
	  C_F0(const Polymorphic *,const char *,const basicAC_F0 & );
	  C_F0(const Polymorphic *,const char *, AC_F0 & );
	  //  function, array ..  
	  C_F0(const C_F0 & e,const char *op,const basicAC_F0 & p)  ;
      C_F0(const C_F0 & e,const char *op, AC_F0 & p) ;	  
	  C_F0(const C_F0 & e,const char *op,const C_F0 & ee)  ; 
      C_F0(const C_F0 & e,const char *op,const C_F0 & a,const C_F0 & b) ; 	  
	  C_F0(const C_F0 & e,const char *nm) ; 
      //  without parameter ex f()   
      C_F0(const Polymorphic * pop,const char *op); 
      // unary operator  
      C_F0(const Polymorphic * pop,const char *op,const C_F0 & a); 
      // binary operator  
      C_F0(const Polymorphic * pop,const char *op,const C_F0 & a,const  C_F0  & b); 
      // ternary operator  
      C_F0(const Polymorphic * pop,const char *op,const  C_F0 & a,const  C_F0 & b,const  C_F0 & c); 
	  
	  	  
	  C_F0( Expression ff,aType rr ): f(ff),r(rr) { 
	  // if (!rr && ff)  cerr << "Type Null" << endl;
	    }
	 // operator Expression() const {return f;}
	  AnyType eval(Stack s) const {return (*f)(s);}

	  Expression RightValue() const { return r->RightValueExpr(f);}	  
	  Expression LeftValue() const;
	  
	  aType left() const {return r;}
	  aType right() const {return r->right();}
	  operator  const  E_F0 *  () const {return f;}
	  bool Empty() const {return !f || f->Empty();}
	  bool NotNull() const {return  f;}
	  int  TYPEOFID() const { return r ? r->TYPEOFID(): 0;}
	  int  nbitem() const { return f ? f->nbitem() : 0;}
	  bool EvaluableWithOutStack() const { return f && f->EvaluableWithOutStack();}
	  Expression Destroy() {  return r->Destroy(*this);}
	  operator const Polymorphic * () const {return  dynamic_cast<const Polymorphic *>(f);}
	  bool operator==(const C_F0 & a) const {return f==a.f && r == a.r;}
	  bool operator!=(const C_F0 & a) const {return f!=a.f || r != a.r;}
//          Type_Expr SetParam(const ListOfId * l,size_t & top) const ;
      bool MeshIndependent() const { return f ==0 ? f->MeshIndependent() : false;}
private:
friend class Block;	 
friend class TableOfIdentifier; 
	  C_F0( Expression ff ): f(ff),r(0) {}
 };
//  for bison 
class CListOfInst;
 

 //  a => b
 //  f => t||f
 //  t => t
 //  (a =>b)  <=>  (!a || b )
 
//  warning ------------------
class ForTypeVoid:  public basicForEachType{public:
    ForTypeVoid():basicForEachType(typeid(void),0,0,0,0,0) {}
};

template<class T> 
class ForEachType:  public basicForEachType{public:
    ForEachType(Function1 iv=0,Function1 id=0):basicForEachType(typeid(T),sizeof(T),0,0,iv,id) {
     if (sizeof(T) > sizeof(AnyTypeWithOutCheck) )
      {
        cout << " Sorry the " <<typeid(T).name() << " is to large  ( " << sizeof(T) 
             << " > " << sizeof(AnyTypeWithOutCheck) << " ) " << endl;
       throwassert(sizeof(T) <= sizeof(AnyTypeWithOutCheck) );
      }
    }
};
template<class T> 
class ForEachType<T*>:  public basicForEachType{public:
    ForEachType(Function1 iv=0,Function1 id=0):basicForEachType(typeid(T),sizeof(T),0,0,iv,id) { }
};

template<class A,class B>  AnyType UnRef(Stack,const AnyType &a) ; 
template<class A>  AnyType Initialize(Stack,const AnyType &a) ; 
template<class A>  AnyType Destroy(Stack,const AnyType &a) ; 

//  the type of variable is pointer because we need to write in 
template<class T> 
class ForEachTypePtr:  public basicForEachType { public:
    ForEachTypePtr();
    ForEachTypePtr(Function1 init,Function1 dl);         
    ForEachTypePtr(Function1 dl);
};

template<class T> 
class ForEachTypePtr<T*>:  public basicForEachType { public:
    ForEachTypePtr();
    ForEachTypePtr(Function1 init,Function1 dl);         
    ForEachTypePtr(Function1 dl);
};


template<class T,int RTYPE=1> 
class ForEachTypePtrfspace:  public ForEachTypePtr<T> { public:
    ForEachTypePtrfspace():ForEachTypePtr<T>() {} 
    int TYPEOFID() const {return RTYPE;} 
};


class ForTypeAnyType:  public basicForEachType{public:
    ForTypeAnyType(): basicForEachType(typeid(AnyType),sizeof(AnyType)) {}
      bool CastingFrom(const basicForEachType * ) const {return true;}
	  C_F0 CastTo(const C_F0 & e) const {return e;}     
};


//  for  cast and get value associed to a pointer  

    
template<class A,class B> 
  AnyType Cast(Stack,const AnyType &b) { 
    return   SetAny<A>(static_cast<A>(GetAny<B>(b)));}
    
template<class A,class B,A F(const  B &)> 
  AnyType FCast(Stack,const AnyType &b) { 
    return   SetAny<A>(F(GetAny<B>(b)));}
    
template<class A> 
  AnyType UnRef(Stack,const AnyType &a) { 
    return   SetAny<A>(*GetAny<A*>(a));}

template<class A,class B> 
  AnyType UnRef(Stack,const AnyType &a) { 
    return   SetAny<A>(*GetAny<B*>(a));}
    
    
template<class A> 
  AnyType UnRefCopyPtr(Stack,const AnyType &a) { 
    return   SetAny<A*>(new A(**GetAny<A**>(a))) ;} 
       
    
template<class A> AnyType Initialize(Stack,const AnyType &x){
 A * a=GetAny<A*>(x);
 A *b=new A;// 
  memcpy(a,b,sizeof(A));// bitcopy
  ::operator delete(b); // delete with no destruction 
  return  SetAny<A*>(a);
}

template<class A> AnyType InitializePtr(Stack stack,const AnyType &x){
  A * a=GetAny<A*>(x);
SHOWVERB( cout << " init ptr " << typeid(A*).name() <<  (char *) a  - (char*) stack<< endl);
  *a=0;
  return  x;
}

template<class A> inline AnyType Delete(Stack,const AnyType &x){
  A * a=GetAny<A*>(x);
  SHOWVERB(cout << "DESTROY " <<typeid(A).name() << " " << a <<  endl); 
  (*a).~A(); 
  return  Nothing;
}

template<class A> inline AnyType Destroy(Stack,const AnyType &x){
  A * a=GetAny<A*>(x);
  SHOWVERB(cout << "DESTROY " <<typeid(A).name() << " " << a <<  endl); 
  a->destroy(); 
  return  Nothing;
}

template<class A> inline AnyType DestroyS(Stack,const AnyType &x){
  A a=GetAny<A>(x);
  SHOWVERB(cout << "DESTROY " <<typeid(A).name() << " " << a <<  endl); 
  a.destroy(); 
  return  Nothing;
}

template<class A> inline AnyType InitS(Stack,const AnyType &x){
  A  a=GetAny<A>(x);
  SHOWVERB(cout << "InitS " <<typeid(A).name() << " " << a <<  endl); 
  a.init(); 
  return  Nothing;
}
template<class A> inline AnyType InitP(Stack,const AnyType &x){
  A  *a=GetAny<A*>(x);
  SHOWVERB(cout << "InitP " <<typeid(A).name() << " " << a <<  endl); 
  a->init(); 
  return  Nothing;
}


template<class A> inline AnyType DestroyPtr(Stack,const AnyType &x) {
  const A *  a=GetAny<A*>(x);
  SHOWVERB(cout << "DestroyPtr " << typeid(A).name() << *a  << endl);
   (*a)->destroy(); 
   //  delete *a; 

  return  Nothing; 
};
template<class A> inline AnyType DeletePtr(Stack,const AnyType &x) {
  const A *  a=GetAny<A*>(x);
  SHOWVERB(cout << "DeletePtr " << typeid(A).name() << *a  << endl);
  // (*a)->destroy(); 
    delete *a; 

  return  Nothing; 
};

template<> AnyType inline DestroyPtr<string *>(Stack,const AnyType &x) {
  string **  a=GetAny<string**>(x);
 SHOWVERB( cout << "DestroyPtr " << typeid(string*).name() << *a  << endl);
  delete *a; 
  return  Nothing; 
};



template<class A> AnyType Initialize(Stack,const AnyType &x,const AnyType &y){
 A * a=GetAny<A*>(x);
 A *b=new A(GetAny<A>(x));// 
  memcpy(a,b,sizeof(A));// bitcopy
  ::operator delete(b); // delete with no destruction 
  return  SetAny<A*>(a);
}
 

  
class E_F0_CFunc2 :public  E_F0mps { public:
   CFunction2  f2;
   const  E_F0 *a,*b;
   AnyType operator()(Stack s)  const {return f2(s,a,b);}
   E_F0_CFunc2( CFunction2 ff,const E_F0 *aa,const E_F0 *bb) : f2(ff),a(aa),b(bb){}
   bool EvaluableWithOutStack() const 
      {return a->EvaluableWithOutStack() && b->EvaluableWithOutStack();} // 
    operator aType () const { return atype<void>();}         

};

class E_F0_CFunc4 :public  E_F0mps { public:
   CFunction4  f4;
   const  E_F0 *a,*b,*c,*d;
   AnyType operator()(Stack s)  const {return f4(s,a,b,c,d);}
   E_F0_CFunc4( CFunction4 ff,const E_F0 *aa,const E_F0 *bb,const E_F0 *cc,const E_F0 *dd) 
   : f4(ff),a(aa),b(bb),c(cc),d(dd){}
    operator aType () const { return atype<void>();}         

};



template<class R,class A>
 class E_F1_F :public  E_F1 { public:
  typedef  R (*func)(A) ; 
  func f;
  E_F1_F(func ff) : f(ff) {}
  AnyType operator()(Stack s,AnyType & a)  const 
    {return SetAny<R>(f(GetAny<A>(a)));}  
};

template<class R,class A0,class A1>
 class E_F2_F :public  E_F2 { public:
  typedef  R (*func)(const  A0 &,const  A1&) ; 
  func f;
  E_F2_F(func ff) : f(ff) {}
  AnyType operator()(Stack s,AnyType & a0,AnyType & a1)  const 
    {return SetAny<R>(f(GetAny<A0>(a0),GetAny<A1>(a1)));}  
};

template<class R,class TA0>
 class E_F_F0 :public  E_F0 { public:
   template <class T> struct remove_reference     {typedef T type;};
//   template <class T> struct remove_reference<T&> {typedef T type;};
   template <class T> struct remove_reference<const T&> {typedef T type;};
   typedef typename remove_reference<TA0>::type A0;
   
     class Opt: public E_F_F0  { public :
       size_t ia;  
       Opt(const  E_F_F0 &t,size_t iaa) 
         : E_F_F0(t) , ia(iaa) {assert(iaa<2000000 && iaa >0);}
      AnyType operator()(Stack s)  const 
       {
      // A0 x =  *static_cast<A0 *>(static_cast<void*>(static_cast<char *>(s)+ia));
      // cout << " opt f (" << x << " ) = "   << ": " << ia << endl; 
       return SetAny<R>( f( *static_cast<A0 *>(static_cast<void*>(static_cast<char *>(s)+ia))  ) );}  
         
      };      
 
 
  typedef  R (*func)(  TA0 ) ; 
  func f;
  Expression a;
  E_F_F0(func ff,Expression aa) : f(ff),a(aa) {}
  AnyType operator()(Stack s)  const 
    {return SetAny<R>(f(GetAny<A0>( (*a)(s) )));}  
  bool EvaluableWithOutStack() const 
      {return a->EvaluableWithOutStack();} // 
  bool MeshIndependent() const {return a->MeshIndependent();} // 
      
  int compare (const E_F0 *t) const { 
     int rr;
    // cout << "cmp " << typeid(*this).name() << " and " << typeid(t).name() << endl;
     const  E_F_F0* tt=dynamic_cast<const E_F_F0 *>(t);
     if (tt && f == tt->f) rr = a->compare(tt->a);
     else rr = E_F0::compare(t);
     return rr;
     } // to give a order in instuction 

   int Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const
    {
       int rr = find(m);
       if (rr) return rr;
       return insert(new Opt(*this,a->Optimize(l,m,n)),l,m,n);       
    }
    virtual ostream & dump(ostream &f) const  { f << typeid(*this).name() <<" f= " << f << " a= "<< *a << ' '  ;return f; }
      

};


template<class A0>
 class E_VF_F0 :public  E_F0 { public:
  typedef  void (*func)(  A0 ) ; 
  func f;
  Expression a;
  E_VF_F0(func ff,Expression aa) : f(ff),a(aa) {}
  AnyType operator()(Stack s)  const 
    {f(GetAny<A0>( (*a)(s) ));return Nothing;}  
  bool EvaluableWithOutStack() const 
      {return a->EvaluableWithOutStack();} // 
      
 bool MeshIndependent() const { return a->MeshIndependent();  }    

};

inline int clexico(int i,int j) { return i==0 ? j : i;}

template<class R,class TA0,class TA1>
 class E_F_F0F0 :public  E_F0 { public:
   template <class T> struct remove_reference     {typedef T type;};
   template <class T> struct remove_reference<T&> {typedef T type;};
   typedef typename remove_reference<TA0>::type A0;
   typedef typename remove_reference<TA1>::type A1;
   typedef  R (*func)( A0 , A1 ) ;
    
  func f;
  Expression a0,a1;
  E_F_F0F0(func ff,Expression aa0,Expression aa1) 
    : f(ff),a0(aa0),a1(aa1) {}
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)) , GetAny<A1>((*a1)(s)) ) );}  
   bool EvaluableWithOutStack() const 
      {return a0->EvaluableWithOutStack() && a1->EvaluableWithOutStack();} // 
   bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent();} // 
  int compare (const E_F0 *t) const { 
     int rr;
    // cout << "cmp " << typeid(*this).name() << " and " << typeid(t).name() << endl;
     const  E_F_F0F0* tt=dynamic_cast<const E_F_F0F0 *>(t);
     if (tt && f == tt->f) rr= clexico(a0->compare(tt->a0),a1->compare(tt->a1));
     else rr = E_F0::compare(t);
     return rr;
     } // to give a order in instuction 
      
   int Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const
    {

       int rr = find(m);
       if (rr) return rr;

       return insert(new Opt(*this,a0->Optimize(l,m,n),a1->Optimize(l,m,n)),l,m,n);
    }
     // build optimisation
     class Opt: public E_F_F0F0  { public :
       size_t ia,ib;  
       Opt(const  E_F_F0F0 &t,size_t iaa,size_t ibb) 
         : E_F_F0F0(t) ,
         ia(iaa),ib(ibb) {}
      AnyType operator()(Stack s)  const 
       {
         //A0 aa =*static_cast<A0 *>(static_cast<void*>(static_cast<char *>(s)+ia));
         //A1 bb=*static_cast<A1 *>(static_cast<void*>(static_cast<char *>(s)+ib)) ;
         //cout << ia << " " << ib <<  "f( " << aa << "," << bb  << " )   = "<< f(aa,bb) << endl;
         return SetAny<R>( f( *static_cast<A0 *>(static_cast<void*>(static_cast<char *>(s)+ia)) , 
                             *static_cast<A1 *>(static_cast<void*>(static_cast<char *>(s)+ib)) ) );}  
         
      };     
       
    
};




template<class R,class A0>
 class E_F_F0_ :public  E_F0 { public:
  typedef  R (*func)(const   A0& ) ; 
  func f;
  Expression a;
  E_F_F0_(func ff,Expression aa) : f(ff),a(aa) {}
  AnyType operator()(Stack s)  const 
    {return SetAny<R>(f(GetAny<A0>( (*a)(s) )));}  
   bool EvaluableWithOutStack() const 
      {return a->EvaluableWithOutStack() ;} // 
   bool MeshIndependent() const 
      {return a->MeshIndependent();} // 
    
};

template<class R,class A0>
 class E_F_F0s_ :public  E_F0mps { public:
  typedef  R (*func)(Stack stack,const   A0& ) ; 
  func f;
  Expression a;
  E_F_F0s_(func ff,Expression aa) : f(ff),a(aa) {}
  AnyType operator()(Stack s)  const 
    {return SetAny<R>(f(s,GetAny<A0>( (*a)(s) )));}  
  bool MeshIndependent() const 
      {return true;} // 

    operator aType () const { return atype<R>();}         
    
};

template<class R,class A0,class A1,class E=E_F0>
 class E_F_F0F0_ :public  E { public:
  typedef  R (*func)(const  A0 &,const  A1 & ) ; 
  func f;
  Expression a0,a1;
  E_F_F0F0_(func ff,Expression aa0,Expression aa1) 
    : f(ff),a0(aa0),a1(aa1) {}
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)) , GetAny<A1>((*a1)(s)) ) );} 
    bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent();} // 
 
};

template<class R,class A0,class A1,class A2,class E=E_F0>
 class E_F_F0F0F0_ :public  E { public:
  typedef  R (*func)(const  A0 &,const  A1 & , const A2 &) ; 
  func f;
  Expression a0,a1,a2;
  E_F_F0F0F0_(func ff,Expression aa0,Expression aa1,Expression aa2) 
    : f(ff),a0(aa0),a1(aa1),a2(aa2) {}
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)) , GetAny<A1>((*a1)(s)),GetAny<A2>((*a2)(s))  ) );}  
    virtual size_t nbitem() const {return a2->nbitem(); } 
      bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent()&& a2->MeshIndependent();} // 

};

template<class R,class A0,class A1,class A2,class E=E_F0>
 class E_F_stackF0F0F0_ :public  E_F0mps { public:
  typedef  R (*func)(Stack, const  A0 &,const  A1 & , const A2 &) ; 
  func f;
  Expression a0,a1,a2;
  E_F_stackF0F0F0_(func ff,Expression aa0,Expression aa1,Expression aa2) 
    : f(ff),a0(aa0),a1(aa1),a2(aa2) {}
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f(s, GetAny<A0>((*a0)(s)) , GetAny<A1>((*a1)(s)),GetAny<A2>((*a2)(s))  ) );}  
    virtual size_t nbitem() const {return a2->nbitem(); } 
   bool MeshIndependent() const { return true;}
};

template<class R,class A0,class A1>
 class E_F_F0F0_NC :public  E_F0 { public:
  typedef  R (*func)(  A0 &,const  A1 & ) ; 
  func f;
  Expression a0,a1;
  E_F_F0F0_NC(func ff,Expression aa0,Expression aa1) 
    : f(ff),a0(aa0),a1(aa1) {}
  AnyType operator()(Stack s)  const 
    {return SetAny<R>( f( GetAny<A0>((*a0)(s)) , GetAny<A1>((*a1)(s)) ) );}  
   bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent() ; } // 

};




 class E_F_StackF0F0 :public  E_F0mps { public:
  typedef   AnyType (*func)(Stack,Expression ,Expression ) ; 
  func f;
  Expression a0,a1;
  E_F_StackF0F0(func ff,Expression aa0,Expression aa1) 
    : f(ff),a0(aa0),a1(aa1) { }
  AnyType operator()(Stack s)  const 
    {return  (*f)(s, a0 , a1) ;}  

};


/*
 class E_F_F0F0_<AnyType,AnyType,AnyType> :public  E_F0 { public:
  typedef  AnyType (*func)(const  AnyType &,const  AnyType & ) ; 
  func f;
  Expression a0,a1;
  E_F_F0F0_(func ff,Expression aa0,Expression aa1) 
    : f(ff),a0(aa0),a1(aa1) {}
  AnyType operator()(Stack s)  const 
    {return  f( (*a0)(s) , (*a1)(s) );}  
   bool MeshIndependent() const 
      {return a0->MeshIndependent() && a1->MeshIndependent() ; } // 

};
*/

class E_F2_func :public  E_F2 { public:
   Function2 f;
   AnyType operator()(Stack s,AnyType & a,AnyType & b)  const {return f(s,a,b);}
   E_F2_func( Function2 ff) : f(ff) {}
};

class E_F0_Func1 :public  E_F0 { public:
   Function1  f;
   const  E_F0 *a;
   AnyType operator()(Stack s)  const {return f(s,(*a)(s));}
   E_F0_Func1( Function1 f1,const E_F0 *aa) : f(f1),a(aa){}
   bool EvaluableWithOutStack() const {return a->EvaluableWithOutStack();} // 
   bool MeshIndependent() const {return a->MeshIndependent();} // 
   int compare (const E_F0 *t) const { 
     int rr;
     const  E_F0_Func1* tt=dynamic_cast<const E_F0_Func1 *>(t);
     if (tt && f == tt->f) rr = a->compare(tt->a);
     else rr = E_F0::compare(t);
     if(tt && 0)
      {
       cout << "\n\t\t\t -------- " << (void *) f << " " << (void *) tt->f << " rr=" << a->compare(tt->a) << endl;
       cout << "\t\t\tcmp E_F0_Func1 " << rr <<" << " << *this << " cmp " << *t << " " << tt << ">>\n";
      }
     return rr;
     } // to give a order in instuction 
  // int Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const;  // build optimisation

    virtual ostream & dump(ostream &f) const  { f << "E_F0_Func1 f= " << f << " a= "<< *a << ' '  ;return f; }

};
class E_F0_Func2 :public  E_F0 { public:
   Function2  f;
   const  E_F0 *a,*b;
   AnyType operator()(Stack s)  const {return f(s,(*a)(s),(*b)(s));}
   E_F0_Func2( Function2 f1,const E_F0 *aa,const E_F0 *bb) : f(f1),a(aa),b(bb){}
   bool EvaluableWithOutStack() const 
     {return a->EvaluableWithOutStack() &&b->EvaluableWithOutStack();} // 
   bool MeshIndependent() const {return a->MeshIndependent() && b->MeshIndependent();} // 

};



//  the variable offset / stack (local variable)
template<class R> class Value1:public E_F0
 { 
  size_t offset;
  public:
  AnyType operator()(Stack s) const { return SetAny<R*>(static_cast<R *>(static_cast<void *>(  static_cast<char *>(s)+offset)));}
  Value1(size_t o):offset(o) {}
};

//  the variable globale
template<class R> class GValue:public E_F0
 { 
  mutable R v;
  public:
  AnyType operator()(Stack ) const { return SetAny<R*>(static_cast<R *>(static_cast<void *>(&v)));}
  GValue(R o):v(o) {}
  bool EvaluableWithOutStack() const {return true;} // 
  
};

//  a constante value 
template<class R> int ccompare(const R & a,const R& b){ return a==b ?  0 :( a<b ? -1 : +1);}
template<class R> int ccompare(const complex<R> & a,const complex<R>& b){ 
  int c=ccompare(a.real(),b.real()); 
  return c==0 ? ccompare(a.imag(),b.imag()): c ;}
  
template<class R> class EConstant:public E_F0
 { 
  const R v;
  public:
  AnyType operator()(Stack ) const { /*cout << " ()" << v << endl*/;return SetAny<R>(v);}
  EConstant(const R & o):v(o) { /*cout << "New constant " << o << endl;*/}
  bool EvaluableWithOutStack() const {return true;} //   
  operator aType () const { return atype<R>();} 
  int compare (const E_F0 *t) const { 
        int rr;
        const  EConstant * tt=dynamic_cast<const EConstant *>(t);
            if (tt) rr = ccompare(v,tt->v);
             else rr = E_F0::compare(t);
         return rr;
       } 
   ostream & dump(ostream &f) const { f << " ((" <<typeid(R).name()  << ") " << v << ") " ;return f;}  
};





//  the variable offset / stack (local variable)

 class LocalVariable:public E_F0
 { 
  size_t offset;
  aType t; //  type of the variable just for check  
  public:
  AnyType operator()(Stack s) const { 
    SHOWVERB( cout << "\n\tget var " << offset << " " <<  t->name() << endl);  
   return PtrtoAny(static_cast<void *>(static_cast<char *>(s)+offset),t);}
  LocalVariable(size_t o,aType tt):offset(o),t(tt) {throwassert(tt);     
     SHOWVERB(cout << "\n--------new var " << offset << " " <<  t->name() << endl);
    }
};


class LocalVariableFES : public LocalVariable { public:
  size_t data;
  LocalVariableFES(size_t o,aType tt,const  size_t & d) 
   : LocalVariable(o,tt),data(d) {}
  size_t nbitem() const { return data;}
};

template <class U>
class LocalVariablePlus : public LocalVariable { public:
  U data;
  LocalVariablePlus(size_t o,aType tt,const  U & d) 
   : LocalVariable(o,tt),data(d) {}
};

//  global variable bof bof 
template<class R> class PValue:public E_F0
 { 
  R * p;
  public:
  AnyType operator()(Stack ) const { return SetAny<R*>(p);}
  PValue(R * pp):p(pp) {}
};
template<class R> class PPValue:public E_F0
 { 
  R ** p;
  public:
  AnyType operator()(Stack ) const { return SetAny<R*>(*p);}
  PPValue(R ** pp):p(pp) {}
};


template<class R>
Type_Expr CPValue(R & v)
 {
   throwassert(map_type[typeid(R*).name()]);
  return make_pair(map_type[typeid(R*).name()],new PValue<R>(&v));
 }
template<class R>
Type_Expr CPPValue(R *& v)
 {
   throwassert(map_type[typeid(R*).name()]);
  return make_pair(map_type[typeid(R*).name()],new PPValue<R>(&v));
 }
 
template<class R >
Type_Expr CConstant(const R & v)
 {
  throwassert(map_type[typeid(R).name()]);
  return make_pair(map_type[typeid(R).name()],new  EConstant<R>(v));
 }


class CC_F0 {
  Expression f;
  aType r;
public:
  void operator=(const C_F0& c) { f=c.f;r=c.r;;} 
  void operator=(const AC_F0& a) ; //{ f=new E_Array(a); f= atype<E_Array>();};
  void operator=(long ) {f=0;r=0;}
  void operator=(const CListOfInst& c);//{ C_FO cc=c;f=cc.f;r=cc.r}
  operator C_F0 () const {return C_F0(f,r);}
  bool Empty() const {return !f || f->Empty();}
  aType left() const {return r;}
 // operator const C_F0 &() const {return  *this;}
};


class ListOfInst : public  E_F0mps { 
    int n;
    Expression   *   list;
    int   *   linenumber;
    const int nx;
    public:
    ListOfInst():n(0),list(0),linenumber(0),nx(10){}
    ListOfInst(int nn):n(0),list(0),linenumber(0),nx(nn?nn:10){}
    void Add(const C_F0 & ins); 
    AnyType operator()(Stack s) const; 
    operator aType () const { return n ? (aType) * (list[n-1]) : atype<void>();} 
   
   Expression &operator[](int i){return list[i];}
   bool empty() const {return n==0;}
   int size() const {return n;}
   Expression * ptr() const {return list;}
   int * nlines() const {return linenumber;}
 //  void destroy() { if (list) delete [] list; list=0;}
 
};

class CListOfInst  {  private:
  ListOfInst * f;
  const basicForEachType *r;
  public:
   void operator=(const CC_F0 &a){
     f=new ListOfInst();     
       if( !a.Empty() ) {
         f->Add(a);
         r=a.left(); }}
   CListOfInst & operator+=(const CC_F0 & a);//{ if( !a.Empty()){ f->Add(a);r=a.left();};return *this;} 
    operator C_F0 () const  { return C_F0(f,r);}
   void eval(Stack s) {(*f)(s);}
   int size() const {return f->size();}
   Expression * ptr() const {return f->ptr();}
   int * nlines() const { return f->nlines();}
};


AnyType FWhile(Stack ,const E_F0 * test,const E_F0 * ins);
AnyType FFor(Stack s ,const E_F0 * i0,const E_F0 * i1,const E_F0 * i2,const E_F0 * ins);
AnyType FIf(Stack s ,const E_F0 * test,const E_F0 * i1,const E_F0 * i2,const E_F0 * notuse);



extern TableOfIdentifier Global;
void ShowType(ostream & );

template<class T> 
inline C_F0 to(const C_F0 & a) { return map_type[typeid(T).name()]->CastTo(a);}


/*
inline C_F0 toBool(const C_F0 & a)    {return ATYPE(bool)->CastTo(a);}
inline C_F0 toInt(const C_F0 & a)     {return ATYPE(int)->CastTo(a);}
inline C_F0 toLong(const C_F0 & a)    {return ATYPE(long)->CastTo(a);}
inline C_F0 toDouble(const C_F0 & a)  {return ATYPE(double)->CastTo(a);}
inline C_F0 toComplex(const C_F0 & a)  {return ATYPE(Complex)->CastTo(a);}
*/
inline C_F0 While(C_F0 test,C_F0 ins) {return C_F0(new E_F0_CFunc2(FWhile,to<bool>(test),ins),0);}
inline C_F0 For(C_F0 i0,C_F0 i1,C_F0 i2,C_F0 ins) {return C_F0(new E_F0_CFunc4(FFor,i0,to<bool>(i1),i2,ins),0);}
inline C_F0 FIf(C_F0 i0,C_F0 i1,C_F0 i2) {return C_F0(new E_F0_CFunc4(FIf,to<bool>(i0),i1,i2,0),0);}
inline C_F0 FIf(C_F0 i0,C_F0 i1) {return C_F0(new E_F0_CFunc4(FIf,to<bool>(i0),i1,0,0),0);}
//inline  C_F0 C_F0::PtrValue() const{ 
//   if (!(r && r->un_ptr)) { cerr << "PtrValue: Not a Left value " << *r << endl;CompileError();} 
//   return C_F0(new  E_F0_Func1(r->un_ptr->f,f),r->un_ptr->r);}

 

class basicAC_F0 {
//  version de base d'un tableau d'un parametres  
//   pour les operateurs unaire, binaire, , 
//   pas d'allocation 
  friend class E_Array; // for mapping fonction 
  protected:
  typedef const C_F0 const_C_F0;
 int nb;
 C_F0 *a;
 public:
  typedef map<const char *,C_F0,Keyless> maptype ;
  typedef maptype::iterator iterator;
  typedef maptype::const_iterator const_iterator;
  maptype * named_parameter;
 basicAC_F0 & operator=(int i) {throwassert(i==0);named_parameter=0,nb=0;return *this;} // pas de parametres
 basicAC_F0 & operator=(C_F0 & c) {named_parameter=0;nb=1;a=&c;return *this;}
 basicAC_F0 & operator=(pair<int,C_F0*> p)  {named_parameter=0;nb=p.first;a=p.second;return *this;} 
 const C_F0 & operator [] (int i) const {throwassert(a && i<nb);return a[i];}
 int size() const {return nb;}
 C_F0 * ptr() const  {return a;}
 C_F0  find(const char * k) const  { 
   if (named_parameter) { const_iterator i=named_parameter->find(k) ;
    if (i == named_parameter->end() ) return C_F0();
    else return i->second;}
   else return C_F0();} 
   
 struct name_and_type{ 
  const char * name;
  const type_info * type;
 } ;
 
  void SetNameParam(int n=0,name_and_type *l=0 , Expression * e=0) const ;
   
};


class AC_F0: public basicAC_F0 { //  a Array of C_F0
//    tableau d'un parametres  max 1024 parametres 
//    avec allocation  
 const static  int MaxSize;
 //  no constructor in this class (this class is in a union )
  public:
  AC_F0 & operator=(pair<const char *,const C_F0> p) {
    named_parameter=0; a=new C_F0[MaxSize]; nb=0;Add(p.first,p.second);return *this;}
  AC_F0 & operator+=(pair<const char *,const C_F0> p) {Add(p.first,p.second);return *this;} 

  AC_F0 & operator=(long k) {throwassert(k==0);named_parameter=0;a=new C_F0[MaxSize]; nb=0;return *this;}     
  AC_F0 & operator=(const C_F0& c) {named_parameter=0; a=new C_F0[MaxSize]; nb=0;a[nb++]=c;return *this;} 
  AC_F0 & operator+=(const C_F0& c) { throwassert(a&& nb<MaxSize);a[nb++]=c;return *this;} 
  AC_F0 & Add(const char * nm,const C_F0 &c)  {
     if (!named_parameter) named_parameter=new maptype();
    iterator i=named_parameter->find(nm);
    if(i==named_parameter->end()) named_parameter->insert(make_pair(nm,c));
    else {cerr << " the named in the list all ready exist "<< nm <<endl; CompileError();}
    return *this;}
  int size() const {return nb;}
  const C_F0 & operator [] (int i) const {throwassert(a && i<nb);return a[i];}
  void destroy() {
       if (a) delete []  a;
        nb=0;
        if(named_parameter) 
          delete named_parameter;}
  
}; 

class  basicAC_F0_wa : public basicAC_F0 { public:
 basicAC_F0_wa(const C_F0 & e) {
   named_parameter=0;
   nb=1;
   a= new C_F0[nb];
   a[0]=e;
 }
 basicAC_F0_wa(const C_F0 & e,const C_F0 & ee) {
   named_parameter=0;
   nb=2;
   a= new C_F0[nb];
   a[0]=e;
   a[1]=ee;
 }
 

 basicAC_F0_wa(C_F0 e,const basicAC_F0 & b) { 
   named_parameter=0;
   if (b.named_parameter) named_parameter = new maptype(*b.named_parameter);
   nb=1+b.size();
   a= new C_F0[nb];
   a[0]=e;
   for (int i=1;i<nb;i++) a[i]=b[i-1];}
 ~basicAC_F0_wa(){delete [] a;} 

 basicAC_F0_wa(const basicAC_F0 & b) { 
   named_parameter=0;
   if (b.named_parameter) named_parameter = new maptype(*b.named_parameter);
   nb=b.size();
   a= new C_F0[nb];
   for (int i=0;i<nb;i++) a[i]=b[i];}
   
  private: 
   void operator=(const basicAC_F0 & b);
};



class E_Array  :public E_F0 {  public: 
  basicAC_F0_wa *v;// the value
  E_Array(const basicAC_F0 & aa) : v(new basicAC_F0_wa(aa))  {throwassert(v);}
  AnyType operator()(Stack)  const {
     cerr << " Pas evaluation d'un E_array" << endl;
     throwassert(0);
     return  Nothing;}
 const C_F0 & operator [] (int i) const {throwassert(v );return (*v)[i];}
 int size() const {return v->size();}
 size_t nbitem() const {return v->size();}
 void map(C_F0 (*mapping)(const C_F0 & )) const 
   { for (int i=0;i<v->size();i++) 
      v->a[i]=(*mapping)(v->a[i]);}
  virtual bool MeshIndependent() const {
    for (int i=0;i<v->size();i++) 
      if (v->a[i].MeshIndependent()) return false;
     return false;
   
  } // 
  operator aType () const { return atype<void>();} 
 
  };
  
class E_Border ;
class E_BorderN :public E_F0mps { public: 
   const E_Border * b;
   Expression  n;
   const E_BorderN * next;
   E_BorderN(const E_Border *  bb, C_F0  nn,const E_BorderN * nx=0) ; 
   E_BorderN(const E_BorderN & bb,const E_BorderN * nx)
     : b(bb.b),n(bb.n),next(nx) { }   
  AnyType operator()(Stack)  const {
     return  SetAny<const  E_BorderN *>(this);}  
  operator aType () const { return atype<const  E_BorderN *>();}         
     
   E_BorderN * operator+( const E_BorderN & bb)  const 
   { throwassert(bb.next==0);
     return new E_BorderN(bb,this);}
  long Nbseg(Stack stack) const { return GetAny<long>((*n)(stack));}
  double from(Stack stack) const ;//{ return GetAny<double>((*n)(stack));}
  double to(Stack stack) const ;//{ return GetAny<double>((*b)(stack));}
  double * var(Stack stack) const ;//{ return GetAny<double*>((*n)(stack));}
  void code(Stack stack) const ;
  long label()const  ;
  void Plot(Stack stack) const ;
  void BoundingBox(Stack stack,double  &xmin,double & xmax, double & ymin,double & ymax) const ;
};

class AddBorderOperator: public  OneOperator{
  typedef const E_BorderN * A;
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     {
       A a0=dynamic_cast<A>((Expression) args[0]);
       A a1=dynamic_cast<A>((Expression) args[1]);
       assert( a0 && a1);
       assert(a1->next==0);
       return  new E_BorderN(*a1,a0);} 
    AddBorderOperator(): 
      OneOperator(map_type[typeid(A).name()],map_type[typeid(A).name()],map_type[typeid(A).name()])
      {pref = 0;}

};


class  OneOperator_borderN : public OneOperator {public:
  const  E_Border * theborder;
    E_F0 * code(const basicAC_F0 & a) const 
     { return  new E_BorderN(theborder,a[0]);} 
    OneOperator_borderN(const  E_Border * b)
      : OneOperator(atype<const E_BorderN *>(),atype<long>()),
      theborder(b){}
};

class E_Border  :public Polymorphic  {  public: 
   static basicAC_F0::name_and_type name_param[] ;
   static const int n_name_param =0;
  static long Count;
  Expression xvar,xfrom,xto,xcode;
  basicAC_F0_wa * tab;
  long label;
  E_Border(const E_Array * a) : tab(a? a->v:0) , xvar(0),xfrom(0),xto(0),xcode(0),label(++Count) 
   {
       assert(tab); 
       Add("(",new OneOperator_borderN(this));
   }
  E_Border(const basicAC_F0 & aa) :    
      xvar(to<double*>(aa[0])),
      xfrom(to<double>(aa[1])),
      xto(to<double>(aa[2])),
      xcode(aa[3].LeftValue()),
      tab(0),
      label(++Count)
 {
      Add("(",new OneOperator_borderN(this));}
      AnyType operator()(Stack)  const {
        return  SetAny<const  E_Border *>(this);}
     double length(Stack stack) const { assert(0); /* a faire */ }
};
  
inline  E_BorderN::E_BorderN(const E_Border *bb, C_F0  nn,const E_BorderN * nx)  
     :b(bb),n(::to<long>(nn)),next(nx) { throwassert(b);}

inline  double E_BorderN::from(Stack stack) const { return b->xfrom ? GetAny<double>((*b->xfrom)(stack)): double(0.0);}
inline  double  E_BorderN::to(Stack stack) const { return b->xto? GetAny<double>((*b->xto)(stack)): b->length(stack) ;}
inline  double *  E_BorderN::var(Stack stack) const { return b->xvar ? GetAny<double*>((*b->xvar)(stack)): (double*) 0 ;}
inline  void  E_BorderN::code(Stack stack)const { (*b->xcode)(stack);}
inline  long  E_BorderN::label()const { return b->label;}

inline ArrayOfaType::ArrayOfaType(const basicAC_F0 & aa) : n(aa.size()),t(n ? (n<=4 ? tt : new aType[n]):0),ellipse(false) { 
   for (int i=0;i<n;i++) t[i]=aa[i].left();}
   
inline ArrayOfaType::ArrayOfaType(const ArrayOfaType & aa) : n(aa.n),t(n<=4?tt:new aType[n]),ellipse(aa.ellipse) { 
   for (int i=0;i<n;i++) t[i]=aa.t[i];}   


inline C_F0 TableOfIdentifier::Find(const char * name) const  {
    const_iterator i=m.find(name); 
    if ( i == m.end()) { return C_F0();}
    else return C_F0(i->second);}

inline C_F0 TableOfIdentifier::Find(const char * name,const basicAC_F0 & args) const  {
    const_iterator i=m.find(name); 
    if ( i == m.end()) {cerr<<"No operator " << name<<endl;
    cerr <<*this << endl;CompileError("TableOfIdentifier::Find");return C_F0();}
    else {return C_F0(C_F0(i->second),"(",args);}}
//  Attention il y a moralement un bug
//  les initialisation   x = y   ( passe par l'operateur binaire <-  dans TheOperators
//   les initialisation   x(y)   ( passe par l'operateur unaire <-  du type de x
//   -------

inline size_t align8(size_t &off) 
{ 
 register size_t o= off %8 ;
  off += o ? 8-o : 0;
 return off;
}


template<class T>
inline Type_Expr  NewVariable(aType t,size_t &off) 
{ 
   size_t o= align8(off);//  align    
 //  off += t->un_ptr_type->size;
 // bug    off += t->size;
   off += t->un_ptr_type->size; // coorection 16/09/2003 merci ˆ Richard MICHEL
   return  Type_Expr(t,new T(o,t));
} 

template<class T>
inline Type_Expr  NewVariable(aType t,size_t &off,const basicAC_F0 &args) 
{ 
   size_t o= align8(off);//  align    
   off += t->un_ptr_type->size;
   return  Type_Expr(t,new T(o,t,args));
}

template<class T,class U>
inline Type_Expr  NewVariable(aType t,size_t &off,const U & data) 
{ 
   size_t o= align8(off);//  align    
   off += t->un_ptr_type->size;
   return  Type_Expr(t,new T(o,t,data));
}

template<class T>   
inline  C_F0 TableOfIdentifier::NewVar(Key k,aType t,size_t & top,const C_F0 &i) 
   { 
     return C_F0(TheOperators,"<-",New(k,NewVariable<T>(t,top)),i);}

template<class T>   
inline  C_F0 TableOfIdentifier::NewVar(Key k,aType t,size_t & top,const basicAC_F0 &args) 
   {  
 //      return C_F0(TheOperators,"<-",New(k,NewVariable(t,top)),t->Find("<-",args));}
        return C_F0(TheOperators,"<-",basicAC_F0_wa(New(k,NewVariable<T>(t,top)),args));}
        
template<class T>   
inline  C_F0 TableOfIdentifier::NewFESpace(Key k,aType t,size_t & top,const basicAC_F0 &args) 
   {  
        return C_F0(TheOperators,"<-",basicAC_F0_wa(New(k,NewFESpace<T>(t,top,args)),args));}
        

template<class T,class U>   
inline  C_F0 TableOfIdentifier::NewVar(Key k,aType t,size_t & top,const basicAC_F0 &args,const U & data) 
   {  
 //      return C_F0(TheOperators,"<-",New(k,t->NewVar(top)),t->Find("<-",args));}
        return C_F0(TheOperators,"<-",basicAC_F0_wa(New(k,NewVariable<T,U>(t,top,data)),args));}
   
//inline  C_F0 TableOfIdentifier::NewVar(Key k,aType t,size_t & top,const AC_F0 &args,const C_F0& ) 
//   {   throwassert(0); return C_F0(TheOperators,"<-",New(k,NewVariable(t,top)),t->Find("<-",args));}

template<class T>   
inline  C_F0 TableOfIdentifier::NewVar(Key k,aType t,size_t & top) 
   {  return t->Initialization(New(k,NewVariable<T>(t,top))); }

// save a expression 
inline  C_F0 TableOfIdentifier::NewID(aType r,Key k, C_F0 & c,size_t &top, bool del ) 
   {  New(k,(make_pair<aType,const E_F0  *>(c.left(),c.LeftValue())),del);return 0; }
 //  { return r->Initialization(New(k,r->SetParam(c,ListOfId(),top),del));}

inline  C_F0 TableOfIdentifier::NewID(aType r,Key k, C_F0 & c,const ListOfId & l,size_t & top,bool del) 
   { return r->Initialization(New(k,r->SetParam(c,&l,top),del));}
   


    
typedef list<TableOfIdentifier *> ListOfTOfId;    
    
 extern  list<TableOfIdentifier *> tables_of_identifier;

 


  C_F0 Find(const char * name) ;  
  
inline  C_F0 basicForEachType::Find(const char * k) const
  {  C_F0 r( ti.Find(k));
     //if (r.Empty()) {cerr << " no member " <<k << " in type " << name() << endl; CompileError("  ");}
     return r; }
inline C_F0  basicForEachType::Find(const char * k,const basicAC_F0 & args) const {return ti.Find(k,args);}
inline  C_F0 basicForEachType::Initialization(const Type_Expr & e) const 
  {
     if(!InitExp) 
       { 
          cerr << "Internal Error: Not Way to make the Initialization of this var type " << *this << endl;
          CompileError();
       }
   return C_F0(new  E_F0_Func1(InitExp,e.second),this);        
  }
  

    
//inline  AnyType Args2(const AnyType &,const  AnyType & b) {return b;}
class E_comma : public E_F0 {public:
   Expression a,b;
   E_comma(Expression aa,Expression bb) : a(aa),b(bb) {}
   AnyType operator()(Stack s) const  { (*a)(s); return (*b)(s);}
   bool MeshIndependent() const {
    return a->MeshIndependent() && b->MeshIndependent();}
} ;

inline	  C_F0::C_F0(const C_F0 & a,const C_F0 & b)  
  { // the concatenation 
      if (a.Empty())
         {r=b.r,f=b.f;}
      else if (b.Empty())
         {r=a.r,f=b.f;}
      else 
        {r=b.r;
        f= new E_comma(a.f,b.f);}
  }

inline	  C_F0::C_F0(const C_F0 & e,const char *op, AC_F0 & p)  
   {    *this=C_F0(e,op,(const basicAC_F0 &) p);
        p.destroy();
   }          
inline	  C_F0::C_F0(const Polymorphic * poly,const char *op, AC_F0 & p)
   {    *this=C_F0(poly,op,(const basicAC_F0 &) p);
        p.destroy();
   }          

inline	  C_F0::C_F0(const C_F0 & e,const char *op,const basicAC_F0 & p)  
	   { 
	     const Polymorphic * pop=e;
	     if (pop) 
	      {
	      //cerr << "poly: " <<  *pop << endl;
	      *this=C_F0(pop,op,p);
	      }
	     else { 
	      // cerr << *e.r << " : table  " << endl;
	      // e.r->ShowTable(cerr);
	       C_F0 x=e.r->Find(op);
	       
	       pop=x;
	       if(pop) 
	       	 {  
	       	   basicAC_F0_wa ep(e,p);       
	           *this=C_F0(pop,"",ep); 
	         } 
	       else
	        {
	           cerr << " unknow operator " << op << " on type " << *e.r << endl;
	           CompileError();
	        }}	       
	   }
inline	  C_F0::C_F0(const C_F0 & e,const char *op,const C_F0 & a,const C_F0 & b)  
{
    C_F0 tab[2]={a,b};
    basicAC_F0  p;
    p=make_pair<int,C_F0*>(2,tab);
    *this= C_F0(e,op,p);
}
	   
inline	  C_F0::C_F0(const C_F0 & e,const char *op,const C_F0 & ee)  
{
	     const Polymorphic * pop=e;
	     if (pop) 
	      {
	      *this=C_F0(pop,op,e,ee);
	      }
	     else { 
	      // cerr << *e.r << " : table  " << endl;
	      // e.r->ShowTable(cerr);
	       C_F0 x=e.r->Find(op);       
	       pop=x;
	       if(pop) 	         
	         *this=C_F0(pop,"",e,ee);  
	       else
	        {
	           cerr << " unknow operator " << op << " on type " << *e.r << " " << *ee.r<<  endl;
	           CompileError();
	        }}	       
  
}
	   
inline	  C_F0::C_F0(const C_F0 & e,const char *nm)  
	   { 
	    // cerr << "  C_F0(const C_F0 & e,const char *item) : "   <<  " " << nm << endl;
	     C_F0 x=e.r->Find(nm); 
	     const Polymorphic * pop=x;
	     // cerr << "Find " << *pop << endl;
	     if (pop) 
	       *this=C_F0(pop,".",e); //  unary oper . 
	     else 
	       {
	           cerr << " No operator ." << nm << " for type " << *e.r  <<  endl;
	          lgerror("");	          	        
	       }
	       
	   }
inline const E_F0 * C_F0::LeftValue() const {
    return f;
}

/*inline Type_Expr C_F0::SetParam(const ListOfId * l,size_t & top) const {
    return r->SetParam(*this,l,top);
}*/

inline    void ListOfInst::Add(const C_F0 & ins) { 
       if( (!ins.Empty()) ) {
      if (n%nx==0){ 
                Expression   *  l = new Expression [n+nx];
                int * ln =  new int [n+nx];
      			for (int i=0;i<n;i++) {
      			  l[i]=list[i];
      			  ln[i]=linenumber[i];}
      			delete [] list;
      			delete [] linenumber;
      			list =l;
      			linenumber=ln;
      		    }
      throwassert(list);
      linenumber[n]= TheCurrentLine;		    
      list[n++] = ins;
      }}
      
      
inline   AnyType ListOfInst::operator()(Stack s) const {     
      AnyType r; 
      double s0=CPUtime(),s1=s0,ss0=s0;
      for (int i=0;i<n;i++) 
       {
        TheCurrentLine=linenumber[i];
        r=(*list[i])(s);
        s1=CPUtime();
        if (showCPU)  
          cout << " CPU: "<< i << " " << s1-s0 << "s" << " " << s1-ss0 << "s" << endl;
        s0=CPUtime();
       }
      return r;}
      
aType TypeArray(aType,aType);
aType TypeArray(aType c,aType b,aType a);

void Init_map_type();

class Block { //
   static size_t Max(size_t a,size_t b){ return a < b ? b :a;}
   typedef const char *  Key;
   Block * fatherblock;
   size_t  top,topmax;
   TableOfIdentifier table;
   ListOfTOfId::iterator itabl;    
   public:
   //  list of variable
   size_t OffSet(size_t size) {
      top=align8(top);
     size_t r=top;  top+=size ;topmax=Max(topmax,top);
     return r;}
   Block(Block * f=0):fatherblock(f),top(f?f->top:BeginOffset*sizeof(void*)),topmax(top)
    {     
      itabl=tables_of_identifier.insert(tables_of_identifier.begin(),&table);
    }
   size_t size() const { return Max(topmax,top);}
  void Add(Key k,Key op,OneOperator *p0)  
    { table.Add(k,op,p0);}
   
template<class T>   
   C_F0 NewVar(Key k,aType t,const C_F0 &i) 
     {return table.NewVar<T>(k, t,top,i);}
template<class T>   
   C_F0 NewFESpace(Key k,aType t,const basicAC_F0 &args) 
     {return table.NewFESpace<T>(k, t,top,args);}
template<class T>   
   C_F0 NewVar(Key k,aType t, AC_F0 &args) 
     {C_F0 r= table.NewVar<T>(k, t,top,args);
      args.destroy();
      topmax=Max(topmax,top);
      return r;}
template<class T>   
   C_F0 NewVar(Key k,aType t,const basicAC_F0 &args) 
     {C_F0 r= table.NewVar<T>(k, t,top,args);
      topmax=Max(topmax,top);
      return r;}
template<class T,class U>   
   C_F0 NewVar(Key k,aType t,const basicAC_F0 &args,const U & data) 
     {C_F0 r= table.NewVar<T,U>(k, t,top,args,data);
      topmax=Max(topmax,top);
      return r;}
//   C_F0 NewVar(Key k,aType t,const AC_F0 &args,const C_F0 & f) 
//     {return table.NewVar(k, t,top,args,f);}
template<class T>   
   C_F0 NewVar(Key k,aType t) 
     {C_F0 r= table.NewVar<T>(k, t,top);
      topmax=Max(topmax,top);
      return r;
      }

  // C_F0 NewVar(aType t,Key k,C_F0 f) 
  //   {return table.NewVar(t,k, f);}
  C_F0 NewID(aType t,Key k,C_F0 f,bool del=true) 
     {C_F0 r= table.NewID(t,k, f,top,del);
      topmax=Max(topmax,top);
      return r;}
  C_F0 NewID(aType t,Key k,C_F0 f,const ListOfId & l,bool del=true) 
     {C_F0 r= table.NewID(t,k,f,l,top,del);
      topmax=Max(topmax,top);
      return r;}

   CC_F0  close(Block *& c) {
     tables_of_identifier.erase(itabl);      
     c=fatherblock;
     if (fatherblock) {fatherblock->topmax=topmax;
                       fatherblock->top=top;}
        
     CC_F0 r;
     r = table.destroy();
     delete this;
     return r;}
   C_F0 Find(const char * k) const  {return table.Find(k);}
   
   ~Block(){} 
}; 




template<class R,class A=R,class CODE=E_F_F0<R,A> >
class  OneOperator1 : public OneOperator {
    aType r; //  return type
    typedef typename CODE::func func; // R (*func)(A) ; 
    func  f;
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new CODE(f,t[0]->CastTo(args[0]));} 
    OneOperator1(func  ff): 
      OneOperator(map_type[typeid(R).name()],map_type[typeid(A).name()]),f(ff){}
};


template<class R,class A=R,class B=A,class CODE=E_F_F0F0<R,A,B> >
class  OneOperator2 : public OneOperator {
    aType r; //  return type 
    typedef typename CODE::func func;
    //typedef  R (*func)(A,B) ;
    func f;
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new CODE(f,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]));} 
    OneOperator2(func  ff): 
      OneOperator(map_type[typeid(R).name()],map_type[typeid(A).name()],map_type[typeid(B).name()]),
      f(ff){}
};

/*template<typename C>
struct OneBinaryOperator_Traits {
  typedef C::result_type R;
  typedef C::first_argument_type A;
  typedef C::second_argument_type B;
};*/

template<class A,class B>  struct SameType { static const int OK=0;};
template<class A>  struct SameType<A,A> { static const int OK=0;};
template<>  struct SameType<bool,bool> { static const int OK=10;};
template<>  struct SameType<long,long> { static const int OK=20;};
template<>  struct SameType<double,double> { static const int OK=30;};
template<>  struct SameType<Complex,Complex> { static const int OK=40;};
template<>  struct SameType<string*,string*> { static const int OK=50;};

template <typename Arg1, typename Arg2,typename Arg3, class Result>
struct ternary_function
{
	typedef Arg1   first_argument_type;
	typedef Arg2   second_argument_type;
	typedef Arg3   third_argument_type;
	typedef Result result_type;
};

template<typename T,class CODE >
class  OneTernaryOperator : public OneOperator{
  typedef typename T::result_type R;
  typedef typename T::first_argument_type A;
  typedef typename T::second_argument_type B;
  typedef typename T::third_argument_type C;

    class Op : public E_F0 {
      typedef  typename C::result_type Result;
         Expression a,b,c;
       public:
       AnyType operator()(Stack s)  const 
        {return  SetAny<R>(static_cast<R>(C::f( GetAny<A>((*a)(s)) ,
                                                GetAny<B>((*b)(s)) ,
                                                GetAny<C>((*c)(s)))));}
       Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {} 
       bool MeshIndependent() const {
       return a->MeshIndependent() && b->MeshIndependent() && c->MeshIndependent();}
    };

   public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new CODE(t[0]->CastTo(args[0]),t[1]->CastTo(args[1]),t[2]->CastTo(args[2]));} 
    OneTernaryOperator(): 
      OneOperator(map_type[typeid(R).name()],
                  map_type[typeid(A).name()],
                  map_type[typeid(B).name()],
                  map_type[typeid(C).name()]) {}
};

template<typename T >
class  OneTernaryOperator3 : public OneOperator{
  typedef typename T::result_type R;
  typedef typename T::first_argument_type A;
  typedef typename T::second_argument_type B;
  typedef typename T::third_argument_type C;

    class Op : public E_F0 {
     // typedef  typename C::result_type Result;
         Expression a,b,c;
       public:
       AnyType operator()(Stack s)  const 
        {return  SetAny<R>(static_cast<R>(T::f( GetAny<A>((*a)(s)) ,
                                                GetAny<B>((*b)(s)) ,
                                                GetAny<C>((*c)(s)))));}
       Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {} 
       bool MeshIndependent() const { return a->MeshIndependent() && b->MeshIndependent() &&  c->MeshIndependent();}
       
    };

   public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new Op(t[0]->CastTo(args[0]),t[1]->CastTo(args[1]),t[2]->CastTo(args[2]));} 
    OneTernaryOperator3(): 
      OneOperator(map_type[typeid(R).name()],
                  map_type[typeid(A).name()],
                  map_type[typeid(B).name()],
                  map_type[typeid(C).name()]) {}
};




template<typename C>
class  OneBinaryOperator : public OneOperator{
  typedef  typename C::result_type R;
  typedef typename C::first_argument_type A;
  typedef typename C::second_argument_type B;

    class Op : public E_F0 {
      typedef  typename C::result_type Result;
         Expression a,b;
       public:
       AnyType operator()(Stack s)  const 
        {return  SetAny<R>(static_cast<R>(C::f( GetAny<A>((*a)(s)) , GetAny<B>((*b)(s)))));}
       Op(Expression aa,Expression bb) : a(aa),b(bb) {} 
       bool MeshIndependent() const { return a->MeshIndependent() && b->MeshIndependent();}
       int Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const
         {
          int rr = find(m);
          if (rr) return rr;          
          int Opa = a->Optimize(l,m,n);          
          int Opb =b->Optimize(l,m,n);
          return insert(new Opt(*this,Opa,Opb),l,m,n);       
          } 
   int compare (const E_F0 *t) const { 
     int rr;
     const  Op * tt=dynamic_cast<const Op *>(t);
     if (tt ) rr =   clexico(a->compare(tt->a),b->compare(tt->b));
     else rr = E_F0::compare(t);
    // cout << "cmp E_F0_Func1 " << rr << endl;
     return rr;
     } // to give a order in instuction 
  // int Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const;  // build optimisation

    virtual ostream & dump(ostream &f) const  { 
       f << "Op<" << typeid(C).name() 
         << ">   \n\t\t\t( a= "<< *a<< ")  \n\t\t\t(b= "<< *b << ") "  ;
      return f; }
          
     // build optimisation
     class Opt: public Op  { public :
       size_t ia,ib;  
       Opt(const  Op &t,size_t iaa,size_t ibb) 
         : Op(t) ,
         ia(iaa),ib(ibb) {}
      AnyType operator()(Stack s)  const 
       {
       // cout <<  "Opt2 ::: " << ia << " "<< ib << " f = " 
       //      <<  GetAny<double>(SetAny<R>(C::f( *static_cast<A *>(static_cast<void*>(static_cast<char *>(s)+ia)) , 
        //                     *static_cast<B *>(static_cast<void*>(static_cast<char *>(s)+ib))))) << endl;
              
              
         return SetAny<R>( C::f( *static_cast<A *>(static_cast<void*>(static_cast<char *>(s)+ia)) , 
                             *static_cast<B *>(static_cast<void*>(static_cast<char *>(s)+ib)) ) );}  
                        
         
      };     
          
                
    };
 //   aType r; //  return type 
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { //cout << "A op B \n" ;
       return  new Op(t[0]->CastTo(args[0]),t[1]->CastTo(args[1]));} 
    OneBinaryOperator(): 
      OneOperator(map_type[typeid(R).name()],map_type[typeid(A).name()],map_type[typeid(B).name()])
      {pref = SameType<A,B>::OK ;}
};

/* essai d'unification des classes 

template<class R,class A,R ff(A),class AA=A> 
struct F_1 : unary_function<AA,R>,public E_F0 {
       AnyType operator()(Stack s)  const 
        { return  SetAny<R>( ff(GetAny<A>((*a)(s)))) ;}

};


template<class C>
class bUnary_Op : public C { public:

         Expression a;
       public:
        
       bUnary_Op(Expression aa) : a(aa) {} 
      
       int compare (const E_F0 *t) const { 
           int rr;
            const  bUnary_Op * tt=dynamic_cast<const bUnary_Op *>(t);
            if (tt) rr = a->compare(tt->a);
             else rr = E_F0::compare(t);
             // cout << "cmp E_F0_Func1 " << rr << endl;
         return rr;
       } // to give a order in instuction 
      bool EvaluableWithOutStack() const {return a->EvaluableWithOutStack();} // 
     bool MeshIndependent() const {return a->MeshIndependent();} // 
       
  // int Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const;  // build optimisation

    virtual ostream & dump(ostream &f) const  { 
       f << "Op1<" << typeid(C).name() 
         << ">   \n\t\t\t( a= "<< *a<< ") "  ;
      return f; }       
       
    };
*/ 
template<class C>
class Unary_Op : public E_F0 { public:


  typedef typename C::result_type R;
  typedef typename C::argument_type A; 
  
       Expression a;
       public:
       AnyType operator()(Stack s)  const 
        { return  SetAny<R>( C::f(GetAny<A>((*a)(s)))) ;}
        
       Unary_Op(Expression aa) : a(aa) {} 
      
       int compare (const E_F0 *t) const { 
           int rr;
            const  Unary_Op * tt=dynamic_cast<const Unary_Op *>(t);
            if (tt) rr = a->compare(tt->a);
             else rr = E_F0::compare(t);
             // cout << "cmp E_F0_Func1 " << rr << endl;
         return rr;
       } // to give a order in instuction 
      bool EvaluableWithOutStack() const {return a->EvaluableWithOutStack();} // 
     bool MeshIndependent() const {return a->MeshIndependent();} // 
       
  // int Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const;  // build optimisation

    virtual ostream & dump(ostream &f) const  { 
       f << "Op1<" << typeid(C).name() 
         << ">   \n\t\t\t( a= "<< *a<< ") "  ;
      return f; }
       
       
    };



template<class C,class Op=Unary_Op<C> > 
class  OneUnaryOperator : public OneOperator{
  typedef typename C::result_type R;
  typedef typename C::argument_type A; 
   // aType r; //  return type 
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new Op(t[0]->CastTo(args[0]));} 
    OneUnaryOperator(): 
      OneOperator(map_type[typeid(R).name()],map_type[typeid(A).name()])
      {}
};

template<class R,class A=R>
class  OneOperator1s_ : public OneOperator {
    aType r; //  return type
    typedef  R (*func)(Stack stack, const A &) ; 
    func  f;
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new E_F_F0s_<R,A>(f,t[0]->CastTo(args[0]));} 
    OneOperator1s_(func  ff): 
      OneOperator(map_type[typeid(R).name()],map_type[typeid(A).name()]),f(ff){}
};

template<class R,class A=R>
class  OneOperator1_ : public OneOperator {
    aType r; //  return type
    typedef  R (*func)(const A &) ; 
    func  f;
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new E_F_F0_<R,A>(f,t[0]->CastTo(args[0]));} 
    OneOperator1_(func  ff): 
      OneOperator(map_type[typeid(R).name()],map_type[typeid(A).name()]),f(ff){}
};

template<class R,class A,class B,class E> class E_F_F0F0_;




template<class R,class A=R,class B=A,class CODE=E_F_F0F0_<R,A,B,E_F0> >
class  OneOperator2_ : public OneOperator {
    aType r; //  return type 
    typedef typename  CODE::func  func;
    func f;
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new CODE(f,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]));} 
    OneOperator2_(func  ff): 
      OneOperator(map_type[typeid(R).name()],map_type[typeid(A).name()],map_type[typeid(B).name()]),
      f(ff){}
};

template<class R,class A=R,class B=A,class C=B,class CODE=E_F_F0F0F0_<R,A,B,C,E_F0> >
class  OneOperator3_ : public OneOperator {
    aType r; //  return type 
    typedef typename  CODE::func  func;
    func f;
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new CODE(f,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]),t[2]->CastTo(args[2]));} 
    OneOperator3_(func  ff): 
      OneOperator(map_type[typeid(R).name()],map_type[typeid(A).name()],map_type[typeid(B).name()],map_type[typeid(C).name()]),
      f(ff){}
};

//  la class code doit contenir
/*
  class CODE: public E_F0 {
    typedef  ...  func .. ;
    typedef .. R;
     static ArrayOfaType  typeargs(); // the list of type de l'operateur of the args 
    typedef  ... R;  // return type 
}
*/
//
template<class CODE,int ppref=0>
class  OneOperatorCode : public OneOperator {
    public: 
    E_F0 * code(const basicAC_F0 & args) const  { return CODE::f(args);} 
    OneOperatorCode():  OneOperator(atype<typename CODE::Result>(),CODE::typeargs()) {pref=ppref;}
    OneOperatorCode(aType r,const ArrayOfaType & l):  OneOperator(r,l)  {pref=ppref;}
    OneOperatorCode(aType r,aType a):  OneOperator(r,a)  {pref=ppref;}
    OneOperatorCode(aType r,aType a,aType b):  OneOperator(r,a,b)  {pref=ppref;}
    OneOperatorCode(aType r,aType a,aType b,aType c):  OneOperator(r,a,b,c)  {pref=ppref;}
    
};

template<class A,class B> struct binary_trait{ typedef  A R  ;}; 
template<>  struct binary_trait<int,double> { typedef  double R;}; 
template<>  struct binary_trait<long,double> { typedef  double R;}; 
template<>  struct binary_trait<int,complex<double> > { typedef  complex<double> R;}; 
template<>  struct binary_trait<long,complex<double> > { typedef  complex<double> R;}; 
template<>  struct binary_trait<double,complex<double> > { typedef  complex<double> R ;}; 
template<class A>  struct binary_trait<A,string* > { typedef  string*  R ;}; 

template<class T> 
 ForEachTypePtr<T>::ForEachTypePtr(): 
         basicForEachType(typeid(T*),sizeof(T*),
//         new E_F1_funcT<T,T*>(UnRef<T>),atype<T>(),
         new E_F1_funcT_Type(atype<T>(),this,UnRef<T>),atype<T>(),

         ::Initialize<T>,::Delete<T>){}
         
template<class T> 
 ForEachTypePtr<T>::ForEachTypePtr(Function1 init,Function1 dl): 
         basicForEachType(typeid(T*),sizeof(T*),
//         new E_F1_funcT<T,T*>(UnRef<T>),atype<T>(),
         new E_F1_funcT_Type(atype<T>(),this,UnRef<T>),atype<T>(),

         init,dl){}
         
template<class T> 
 ForEachTypePtr<T>::ForEachTypePtr(Function1 dl): 
         basicForEachType(typeid(T*),sizeof(T*),
         new E_F1_funcT_Type(atype<T>(),this,UnRef<T>),atype<T>(),
         ::Initialize<T>,dl){}
         

template<class T> 
 ForEachTypePtr<T*>::ForEachTypePtr(): 
         basicForEachType(typeid(T**),sizeof(T**),
//         new E_F1_funcT<T*,T**>(UnRef<T*>),atype<T*>(),
         new E_F1_funcT_Type(atype<T*>(),this,UnRef<T*>),atype<T*>(),

         ::InitializePtr<T*>,::DestroyPtr<T*>){}
      
template<class T> 
 ForEachTypePtr<T*>::ForEachTypePtr(Function1 init,Function1 dl): 
         basicForEachType(typeid(T**),sizeof(T**),
        // new E_F1_funcT<T*,T**>(UnRef<T*>),atype<T*>(),
         new E_F1_funcT_Type(atype<T*>(),this,UnRef<T*>),atype<T*>(),
         init,dl){}
        
template<class T> 
 ForEachTypePtr<T*>::ForEachTypePtr(Function1 dl): 
         basicForEachType(typeid(T**),sizeof(T**),
//         new E_F1_funcT<T*,T**>(UnRef<T*>),atype<T*>(),
         new E_F1_funcT_Type(atype<T*>(),this,UnRef<T*>),atype<T*>(),
         ::InitializePtr<T*>,dl){}
         
/* class  FuncForEachType : public basicForEachType {public:
  FuncForEachType(const basicForEachType * t);
  const basicForEachType *  rtype;  
 };        
*/
  

inline basicForEachType::basicForEachType(const type_info  & k,
                                          const size_t s,
                                          const E_F1_funcT_Type * p,
                                          basicForEachType *rr,
                                          Function1 iv,Function1 id) 
      : ktype(&k),//ktypefunc(0),
        size(s),
        un_ptr(p),
        un_ptr_type(rr?rr:this), 
        InitExp(iv),
        casting(0), // no casting to 
        //funct_type(0),
        destroy(id) {} 
/*
inline basicForEachType::basicForEachType(const type_info  & k, const type_info  & kf,
                                          const size_t s,
                                          const E_F1_funcT_Type * p,
                                          basicForEachType *rr,
                                          Function1 iv,Function1 id) 
      : ktype(&k),ktypefunc(&kf),
        size(s),
        un_ptr(p),
        un_ptr_type(rr?rr:this), 
        InitExp(iv),        
        destroy(id) ,
        funct_type(new FuncForEachType(this)){} 
        
        
*/

inline C_F0 & operator+=(C_F0 & a,C_F0 &b)
{
   C_F0 r = C_F0(TheOperators,"+",a,b);
   a=r;
   return a;
}


template<class T>
  void Dcl_TypeandPtr (Function1 i,Function1 d,Function1 pi,Function1 pd)
   {
      map_type[typeid(T).name()] = new ForEachType<T>(i,d); 
      map_type[typeid(T*).name()] = new ForEachTypePtr<T>(pi,pd); 
   }


template<class T>
  void Dcl_TypeandPtr (Function1 pi,Function1 pd)
   {
      map_type[typeid(T).name()] = new ForEachType<T>(); 
      map_type[typeid(T*).name()] = new ForEachTypePtr<T>(pi,pd); 
   }
   
template<class T>
  void Dcl_TypeandPtr (Function1 pd)
   {
      map_type[typeid(T).name()] = new ForEachType<T>(); 
      map_type[typeid(T*).name()] = new ForEachTypePtr<T>(pd); 
   }
   
template<class T>
  void Dcl_TypeandPtr ()
   {
      map_type[typeid(T).name()] = new ForEachType<T>(); 
      map_type[typeid(T*).name()] = new ForEachTypePtr<T>(); 
   }
   
template<class T>
  aType Dcl_Type (Function1 iv=0,Function1 id=0)
   {
     if (sizeof(T) >sizeof(AnyData)) {
       cerr << " the type   " << typeid(T).name() << " is too large " << sizeof(AnyData) <<  endl;
       throwassert(sizeof(T) <=sizeof(AnyData));}
     return map_type[typeid(T).name()] = new ForEachType<T>(iv,id); 
    
   }

template<class T>
  void Add(const char * k,const char * op,OneOperator *p0,OneOperator *p1=0,
      OneOperator *p2=0,OneOperator *p3=0,OneOperator *p4=0,
      OneOperator *p5=0,OneOperator *p6=0)  
     {atype<T>()->Add(k,op,p0,p1,p2,p3,p4,p5,p6);}     

inline C_F0 operator *(const C_F0 &a,const C_F0 &b)
{    
  return a==*pOne ? b : ( b ==*pOne ? a : C_F0(TheOperators,"*",a,b)) ;}
  
inline C_F0 &operator +=(C_F0 &a,const C_F0 &b)
{  
   C_F0 r=C_F0(TheOperators,"+",a,b);
   a=r;
   return a;}
   
//inline  bool CC_F0::Empty() const {return !f || f->Empty();}
inline  void CC_F0::operator=(const CListOfInst& c) 
  { C_F0 cc=c;f=cc.f;r=cc.r;}
inline   CListOfInst &  CListOfInst::operator+=(const CC_F0 & a)
  { if( !a.Empty()){ f->Add(a);r=a.left();};return *this;} 
  
inline Type_Expr basicForEachType::SetParam(const C_F0 & c,const ListOfId * ,size_t & ) const
     { cerr << " int basicForEachType " << name() << endl; 
       InternalError("basicForEachType::SetParam non defined");  }//return make_pair<aType,const E_F0  *>(c.left(),c.LeftValue());}
     


/*

//  ---  pour les cast ------
class  OneOpCast: public OneOperator { 
    typedef const E_F1_funcT_Type *  CastFunc;
    CastFunc  f;
    public: 
    E_F0 * code(const basicAC_F0 & args) const   { return  new  E_F0_Func1(f->f,args[0]);} 
    OneOpCast(CastFunc  ff): OneOperator(ff->r,ff->a),f(ff){}
};
*/

// 
inline  bool  basicForEachType::CastingFrom(aType t) const  {
     throwassert(this && t);
     if ( t == this) return true; 
     return casting->FindSameR(ArrayOfaType(t,false));
  }

inline  void CerrCast(const pair<const basicForEachType*,const E_F1_funcT_Type *> & i)
{ 
   cerr << "\t" <<  *i.first << ":" << i.second << endl;
}

inline 	 C_F0 basicForEachType::CastTo(const C_F0 & e) const 
{
 throwassert(this);
 aType t = e.left();
 if (this== t) return e;
  
 
  C_F0 ce=e;
  basicAC_F0 at;
  at=ce;
  OneOperator * opcast =casting->FindSameR(ArrayOfaType(t,false));  
  if ( opcast )  
    if ( *opcast == at ) // left value
      return C_F0(opcast->code(at),this);
    else  { // rigth value 
     aType tr = e.right();
     ce = C_F0(e.RightValue(),tr);
     at = ce;
     return C_F0(opcast->code(at),this); }
  else    
      { cerr << "Impossible to cast " << *e.left() << " in " << *this << endl;
           if (casting)  casting->Show(cerr)  ;
           CompileError();} 
 return C_F0();
}
//  1 variable pour les operation de cast 
class E_F1_funcT_Type: public OneOperator{ public:
//  const basicForEachType *r,*a;
  Function1 f;
    E_F0 * code(const basicAC_F0 & args) const   { return  new  E_F0_Func1(f,args[0]);} 
  
  E_F1_funcT_Type(const basicForEachType *rr,const basicForEachType *aa,Function1 ff)
    : OneOperator(rr,aa), f(ff) {}
  
//: r(rr),a(aa),f(ff) {}
//  friend ostream & operator<<(ostream & f,const E_F1_funcT_Type & e) { f << *e.a << " -> " << *e.r ;return f;}
};

template<class R,class A>
class E_F1_funcT :public  E_F1_funcT_Type{ public:   
   E_F1_funcT(Function1 ff) : E_F1_funcT_Type(map_type[typeid(R).name()],map_type[typeid(A).name()],ff){}
   E_F1_funcT(aType r,aType a,Function1 ff) : E_F1_funcT_Type(r,a,ff){}
};

inline Expression  basicForEachType::RightValueExpr(Expression f) const 
{
  if (un_ptr) return new  E_F0_Func1(un_ptr->f,f);
   else return f;        
}

inline void CompileError(string msg,aType r){ 
 string m= r ? msg + "  type: " + r->name() : msg ;
   lgerror(m.c_str());
 }
 
 inline void ExecError(string msg){ 
   cerr << "Fatal error : ExecError " << msg << endl;
   throw(ErrorExec("exit",1));
 }
 


inline  void CC_F0::operator=(const AC_F0& a) {  f=new E_Array(a); r= atype<E_Array>();};

inline  UnId::UnId(const char * idd,const C_F0 & ee,aType rr=0,bool reff=false) 
  :id(idd),r(rr),e(ee),array(0),re(ee.left()) ,ref(reff){}

inline Stack newStack(size_t l)
 {
  Stack thestack = new char[l];
  for (int i = 0;i< l/sizeof(long);i++) ((long*) thestack)[i]=0;
  ((char **) thestack)[0] = new char [1000]; 
  return thestack;
}

inline void deleteStack(Stack s) 
 {
    delete [] (((char **)  s)[0]);
    delete [] (char *) s;
 }

class E_exception : public exception { public:
  enum CODE_exception { UNKNOWN,e_break,e_continue,e_return} ;
  CODE_exception code;  
  AnyType r; // for return 
  public:
  E_exception(CODE_exception c,AnyType rr=Nothing) : code(c),r(rr)  {}
  const int type() {return code;}
  virtual const char *  what() const throw() { return "E_exception (break,continue ou return) "; }           
};


class E_throw : public E_F0mps { public:
   E_exception::CODE_exception kind;    
   Expression ret; // return value
   E_throw(E_exception::CODE_exception c,Expression e=0) :kind(c),ret(e) {}
   AnyType operator()(Stack s)  const { 
     (ret ? throw(E_exception(kind,(*ret)(s)))
          : throw(E_exception(kind)));
       return Nothing; }
   operator aType () const { return atype<void>();} 
      
 } ;

class E_block :  public E_F0mps { public:
  const int n;
  Expression  * code;
  int * linenumber;
  Expression clean;
   E_block(CListOfInst l,C_F0   c)
     : n(l.size()),code(l.ptr()),linenumber(l.nlines()),clean(c) {}
   E_block( C_F0  l,C_F0  c)
     : n(1),code(new Expression),clean(c) { code[0]=l;}
   AnyType operator()(Stack s)  const {
      if (clean) 
       {
         try { 
          for (int i=0;i<n;i++) {
            TheCurrentLine=linenumber[i];
            (*code[i])(s); }}
         catch( E_exception & e) { 
           (*clean)(s); 
            if(verbosity>50)
             cout << " catch " << e.what() << " clean & throw " << endl;
           throw(e); }
       (*clean)(s); 
       }
      else  // not catch  exception if no clean (optimization} 
       for (int i=0;i<n;i++) 
          (*code[i])(s); 
      return Nothing;
   }
    operator aType () const { return atype<void>();}         
   
};
class Routine;
class E_Routine :  public E_F0mps { public:
  Expression code;
  Expression clean;
  aType rt;
  int nbparam;
  Expression * param;
  const char * name;
  E_Routine(const Routine * routine,const basicAC_F0 & args);
   AnyType operator()(Stack s)  const;
  ~E_Routine() { delete [] param;}
  private:
  E_Routine(const E_Routine &);
  void operator=(const E_Routine &);
  operator aType ()  const{ return rt;}         

};
class Routine: public OneOperator{  public:
   size_t offset;
   aType tfunc,tret;
   const char * name;
   const ListOfId &param;
   Block * currentblock;
   Expression  ins;   
   Expression  clean;
   
    E_F0 * code(const basicAC_F0 & args) const  ;
   Routine(aType tf,aType tr,const char * iden, const ListOfId *l,Block * & cb);
   Block * Set(C_F0   instr) ;
};


class TypeLineFunction: public ForEachType<C_F0> {
  public:
  TypeLineFunction() : ForEachType<C_F0>(0,0) {}
  
  void SetArgs(const ListOfId *lid) const {
     if (lid) CompileError("No Argument in line function");
      } 
     
  Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const 
    {  return Type_Expr(c.left(),c.LeftValue());  } 
    
  C_F0 Initialization(const Type_Expr & e) const 
    {  return C_F0(); }  // nothing to initialize 
    
};

template<class T>
class Transpose{ public:
  T  t;
  Transpose( T  v) : t(v) {}
};

class E_F0_Optimize : public E_F0 { 
 deque<pair<Expression,int> > l;
 MapOfE_F0 m;
 int NBbitem;
 int ret;
 public:
  E_F0_Optimize(deque<pair<Expression,int> > &ll,MapOfE_F0 & mm,int rett) :
  l(ll),m(mm),ret(rett),NBbitem(1) {}
  
    virtual AnyType operator()(Stack s)  const {
      int k= l.size();
      for (int i=0;i<k;i++)
        {  int offset = l[i].second;
           *static_cast<AnyType *>(static_cast<void *>((char*)s+offset))= (*l[i].first)(s);
          // cout << " E_F0_Optimize   " << offset << " " <<  *static_cast<double *>(static_cast<void *>((char*)s+offset)) << endl; ;
        }
      return *static_cast<AnyType *>(static_cast<void *>((char*)s+ret));          
    }
    virtual bool Empty() const {return l.size(); }
   // virtual E_F0 * destroy(Stack ) const {return 0;}
  //  virtual const E_F0 * Parameter(Stack ) const {return this;}
    virtual size_t nbitem() const {  return NBbitem;}
    virtual bool EvaluableWithOutStack() const {return false;} // 
    virtual bool MeshIndependent() const {return true;} // 
    virtual E_F0 * right_E_F0() const { return 0;}
    virtual ~E_F0_Optimize() {}
   // virtual int compare (const E_F0 *t) const { return t-this;} // to give a order in instuction 
    virtual  operator aType ()  const { return  *(l.back().first);}   // the type of the expression  
}; 
 
 
inline    int E_F0::find(const MapOfE_F0 & m) const {  //  exp
       // cout << " ffff :" ;
        MapOfE_F0::const_iterator i= m.find(this); 
        if(i != m.end()) {
            if( (verbosity / 100)% 10 == 1) 
              {
                 cout << "\n    find : ";
                 cout  <<  i->second << " mi=" ;
                 cout << MeshIndependent() << " " ;
                 cout << typeid(*this).name()  ;
                 cout << " cmp = " << compare(i->first) ;
                 cout  << " " << i->first->compare(this) << " ";
                 dump(cout);
               }
             assert( compare(i->first) == 0);
           }     
        return i == m.end() ? 0 : i->second ;
    }
 inline   int E_F0::insert(Expression  opt,deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n) const
    {
     int rr=rr=align8(n);
     pair<Expression,int> p(this,rr);
     if( (verbosity / 100)% 10 == 1) 
       cout << " --  insert opt " << n << " " << *this << endl;     
       n += sizeof(AnyType);         
       l.push_back(make_pair<Expression,int>(opt,rr)); 
       m.insert(p); 
       return rr;
     }

extern queue<pair<const E_Routine*,int> > debugstack;

extern basicForEachType *  typevarreal,  * typevarcomplex;  //  type of real and complex variable
  
#endif


