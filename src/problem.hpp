#ifndef PROBLEM_HPP_
#define PROBLEM_HPP_
extern Block *currentblock;

template<class K> class Matrice_Creuse;
template<class K> class MatriceCreuse;
namespace  Fem2D {
  template<class K> class SolveGCPrecon;
  template<class K> class SolveGMRESPrecon;
  template<class K> class SolveGMRESDiag;  
}

template<class K> class SolveGCDiag; 
class Plot;
class v_fes;

typedef FEbase<double> * pferbase ;
typedef FEbaseArray<double> * pferbasearray ;
typedef pair<pferbase,int> pfer ;
typedef pair<pferbasearray,int> pferarray ;
//typedef pair<pmesh *,int> pmesharray ;

typedef   LinearComb<MGauche,C_F0> Finconnue;
typedef   LinearComb<MDroit,C_F0> Ftest;
typedef   LinearComb<pair<MGauche,MDroit>,C_F0> Foperator;

inline int intOp(const MGauche &i) {return i.second;}
inline int intOp(const MDroit &i) {return i.second;}
inline int intOp(pair<MGauche,MDroit> & p) {return Max(intOp(p.first),intOp(p.second));}

inline void SetOp(KN_<bool> & d,const MGauche &i) 
          { d[i.second% last_operatortype]=true;}
inline void SetOp(KN_<bool> & d,const MDroit &i)
          { d[(int) i.second % last_operatortype]=true;}
inline void SetOp(KN_<bool> & d,const pair<MGauche,MDroit> & p) 
          {SetOp(d,p.first);SetOp(d,p.second);}

typedef  const Finconnue  finconnue;
typedef  const Ftest ftest;
typedef  const Foperator foperator;

Expression IsFebaseArray(Expression f);

void SetArgsFormLinear(const ListOfId *lid,int ordre);

inline C_F0 CCastToR(const C_F0 & f){ return C_F0(atype<double>()->CastTo(f),atype<double>());}
inline bool BCastToR(const C_F0 & f){ return atype<double>()->CastingFrom(f.left());}



inline C_F0 CCastToC(const C_F0 & f){ return C_F0(atype<Complex>()->CastTo(f),atype<Complex>());}
inline bool BCastToC(const C_F0 & f){ return atype<Complex>()->CastingFrom(f.left());}

template<class Result>
inline Expression CastTo(const C_F0 & f) { return atype<Result>()->CastTo(f);}

template<class Result>
inline bool BCastTo(const C_F0 & f) { return atype<Result>()->CastingFrom(f.left());}

inline void Check(bool  v,const char * mess)
{
  if (!v) { cerr << " Erreur " << mess ;
  throw(ErrorExec(mess,1));
  }
}           


class EConstantTypeOfFE :public E_F0
{ 
  TypeOfFE * v;
public:
  AnyType operator()(Stack ) const { /*cout << " ()" << v << endl*/;return SetAny<TypeOfFE*>(v);}
  EConstantTypeOfFE( TypeOfFE * o):v(o) { /*cout << "New constant " << o << endl;*/}
  size_t nbitem() const { assert(v);return v->N ;} 
   operator aType () const { return atype<TypeOfFE*>();} 
};



struct TypeSolveMat {
  enum TSolveMat { NONESQUARE=0, LU=1, CROUT=2, CHOLESKY=3, GC = 4 , GMRES = 5, UMFpack=6 };
  TSolveMat t;
  bool sym;
  bool profile;
  TypeSolveMat(TSolveMat tt=LU) :t(tt),
                                 sym(t == CROUT || t ==CHOLESKY  ||  t==GC ),
                                 profile(t != GC && t != GMRES && t != NONESQUARE && t != UMFpack ) {}
};


class TabFuncArg { public:
  typedef double R;  
  typedef  Fem2D::R2 R2;
  Stack s;
  int nb;
  Expression * e;
  void  eval(R *v) const 
  { for (int i=0;i<nb;i++)
    if (e[i])
      { 
        v[i] = GetAny<R>( (*e[i])(s));
        //    cout <<" (" << i << " " << *v << ") ";
      }
    else 
      v[i] = 0;
  }
  R2 eval_X() const {
    return R2(GetAny<R>( (*e[nb-2])(s)),GetAny<R>( (*e[nb-1])(s)));
  }
  void  eval_2(R *v) const 
  { 
    for (int i=0;i<nb-2;i++)
      if (e[i])
        v[i] = GetAny<R>( (*e[i])(s));
      else 
        v[i] = 0;
  }
  TabFuncArg(Stack ss,int n) : s(ss),nb(n),e(new Expression[n]) {}
  void operator=(int j) { throwassert(j==0); for (int i=0;i<nb;i++) e[i]=0;} // resert
  Expression  & operator[](int i){return e[i];}
  ~TabFuncArg() {delete [] e;}
private: // no copy 
  TabFuncArg(const TabFuncArg & );
  void operator=(const TabFuncArg & );  
};


class C_args: public E_F0mps  {public:
  typedef const C_args *  Result;
  list<C_F0> largs;
  typedef list<C_F0> ::const_iterator const_iterator ;
  // il faut expendre 
  C_args() :largs(){}
  C_args(C_F0 c) : largs() { largs.push_back(c);}
  C_args(  const basicAC_F0 & args) :largs(){ 
    int n=args.size();
    for (int i=0;i< n;i++)
      {
        if (args[i].left() == atype<const C_args *>())
          {
            const C_args * a = dynamic_cast<const C_args *>(args[i].LeftValue());
            for (list<C_F0>::const_iterator i=a->largs.begin();i!=a->largs.end();i++)
              largs.push_back(*i);                   
          } 
        else 
          largs.push_back(args[i]);
      };}
  static ArrayOfaType  typeargs() { return ArrayOfaType(true);}
  AnyType operator()(Stack ) const  { return SetAny<const C_args *>(this);}
  operator aType () const { return atype<const C_args *>();}         

  static  E_F0 * f(const basicAC_F0 & args) { return new C_args(args);}    
  bool IsLinearOperator() const;
  bool IsBilinearOperator() const;
};

class C_args_minus: public C_args  {public:
  C_args_minus(  const basicAC_F0 & args) { 
    int n=args.size();
    throwassert(n==2);    
    if (args[0].left() == atype<const C_args *>())
      {
        const C_args * a = dynamic_cast<const C_args *>(args[0].LeftValue());
        throwassert(a);
        for (list<C_F0>::const_iterator i=a->largs.begin();i!=a->largs.end();i++)
          largs.push_back(*i);                   
      }
    else 
      largs.push_back(args[0]);
    
    largs.push_back(C_F0(TheOperators,"-",args[1]));
     }
  
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<const C_args *>(),true);}
  static  E_F0 * f(const basicAC_F0 & args) { return new C_args_minus(args);}    
};

bool isVF(const list<C_F0> & largs); 

template<typename F>
class Minus_Form: public  E_F0mps    {public:
  typedef const F * Result;
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<const F *>());}
  static  E_F0 * f(const basicAC_F0 & args) { 
    int n=args.size();
    throwassert(n==1);
    aType tF=atype<Result>();
    throwassert(args[0].left() == tF);
    Result f = dynamic_cast<Result>(args[0].LeftValue());
    throwassert(f);
    // F mf = -*f;
    F * rf=new F(-*f);
    return  rf;
  }
    operator aType () const { return atype<Result>();} 
    
    
};

template<class RR=double>
class BC_set : public E_F0mps { public:
  typedef const BC_set* Result;
  vector<Expression> on;
  vector<pair<int,Expression> > bc; //  n¡ de l'inconnue+ valeur
  BC_set(  const basicAC_F0 & args) :on(args.size()){
    int n = args.size();      
    throwassert(args.named_parameter);
    AC_F0::const_iterator ii=args.named_parameter->begin();
    AC_F0::const_iterator ie=args.named_parameter->end();
    bc.resize(args.named_parameter->size());
    
    for (int kk=0;ii!=ie;kk++,ii++)
      { //
        C_F0 x=Find(ii->first);
        if (x.left() != atype<const finconnue *>())
          CompileError("We wait for a unkown  u=... of the problem");
        const finconnue * uu = dynamic_cast<const finconnue *>(x.LeftValue());
        throwassert(uu);
        const MGauche *ui=uu->simple();
        throwassert(ui && ui->second == op_id);
        // cout << ii->first << " n¡" << ui->first <<   " = ? " << endl;
        bc[kk]= make_pair(ui->first,CastTo<double>(ii->second));
        ii->second;
      }     
    for (int i=0;i<n;i++)
      on[i]=CastTo<long>(args[i]);
  }
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<long>(),true);}
  AnyType operator()(Stack ) const  { return SetAny<Result>(this);}
  operator aType () const { return atype<Result>();}         
  
  static  E_F0 * f(const basicAC_F0 & args) { return new BC_set(args);}         
  //    void init(Stack stack) const {}
  
};
class CDomainOfIntegration: public E_F0mps { 
public:
  static const int n_name_param =7;
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  enum typeofkind  { int2d=0, int1d=1, intalledges=2} ; 
  typeofkind  kind; //  0 
  typedef const CDomainOfIntegration* Result;
  Expression Th; 
  vector<Expression> what;
  CDomainOfIntegration( const basicAC_F0 & args,typeofkind b=int2d) :kind(b), what(args.size()-1)
  {
    args.SetNameParam(n_name_param,name_param,nargs);
    Th=CastTo<pmesh>(args[0]);
    int n=args.size();
    for (int i=1;i<n;i++)
      what[i-1]=CastTo<long>(args[i]); 
   // cout << " CDomainOfIntegration " << this << endl;       
  }
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<pmesh>(), true);}// all type
  AnyType operator()(Stack ) const  { return SetAny<const CDomainOfIntegration *>(this);}
       operator aType () const { return atype<const CDomainOfIntegration *>();}         

  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args);} 
  const Fem2D::QuadratureFormular & FIT(Stack) const ;   
  const Fem2D::QuadratureFormular1d & FIE(Stack) const ;  
  bool UseOpt(Stack s) const  {  return nargs[5] ? GetAny<bool>( (*(nargs[5]))(s) )  : 1;}
  double  binside(Stack s) const { return nargs[6] ? GetAny<double>( (*(nargs[6]))(s) )  : 0;} // truc pour FH
};  

class CDomainOfIntegrationBorder: public CDomainOfIntegration { 
public:
  CDomainOfIntegrationBorder( const basicAC_F0 & args) :CDomainOfIntegration(args,int1d) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int1d);}    
};

class CDomainOfIntegrationAllEdges: public CDomainOfIntegration { 
public:
  CDomainOfIntegrationAllEdges( const basicAC_F0 & args) :CDomainOfIntegration(args,intalledges) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,intalledges);}    
};






class Problem  : public Polymorphic { 
  typedef double R;
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =8;
  int Nitem,Mitem;
public:
  struct Data {
    const Mesh * pTh; 
    CountPointer<const FESpace> Uh;
    CountPointer<const FESpace> Vh;
    CountPointer<MatriceCreuse<R> > A;   
    void init()  {pTh=0; Uh.init(),Vh.init();A.init();}
    void destroy() { 
      pTh=0;
      Uh.destroy(); 
      Vh.destroy();
      A.destroy();}        
  } ;
  const OneOperator *precon;
  
  mutable vector<Expression> var; // list des var pour les solutions et test 
  bool complextype,VF;
  const C_args *op; // the list of all operator 
  //   Expression noinit,type,epsilon;
  Expression nargs[n_name_param];
  
  const size_t  offset;         
  Problem(const C_args * ca,const ListOfId &l,size_t & top) ;
  static ArrayOfaType  typeargs() { return ArrayOfaType(true);}// all type
  
  Data * dataptr (Stack stack) const {return   (Data *) (void *) (((char *) stack)+offset);}
  void init(Stack stack) const  {
    // cout << " init  " << (char *) dataptr(stack) - (char*) stack  << " " << offset <<  endl; 
    dataptr(stack)->init();}
  void destroy(Stack stack)  const  {dataptr(stack)->destroy();}
  AnyType operator()(Stack stack) const;
  
  bool Empty() const {return false;}     
  size_t nbitem() const { return Nitem;}     
};


class Solve : public Problem { public:
  // just a problem with implicit solve 
  Solve(const C_args * ca,const ListOfId &l,size_t & top) 
    : Problem(ca,l,top) {} 
};

class FormBilinear : public E_F0mps { public:
  typedef const FormBilinear* Result;
  typedef const CDomainOfIntegration * A;
  typedef const foperator  *  B;
  A  di;
  B b;    
  FormBilinear(const basicAC_F0 & args) {
   di= dynamic_cast<A>(CastTo<A>(args[0]));
   B  bb= dynamic_cast<B>(CastTo<B>(args[1]));
   b = bb->Optimize(currentblock);
  // delete bb; il ne faut pas detruire .. car bb peut etre dans une variable 
    throwassert(di && b); }
  
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
  AnyType operator()(Stack ) const { return SetAny<Result>(this);}
   operator aType () const { return atype<Result>();}         

  static  E_F0 * f(const basicAC_F0 & args) { return new FormBilinear(args);}  
  FormBilinear(A a,Expression bb) : di(a),b(dynamic_cast<B>(bb)->Optimize(currentblock)) {throwassert(b);}   
  FormBilinear operator-() const { return  FormBilinear(di,C_F0(TheOperators,"-",C_F0(b,atype<B>())));}
  bool VF() const { return MaxOp(b) >= last_operatortype;}

  //  void init(Stack stack) const {}
};


class FormLinear : public E_F0mps { public:
  typedef const FormLinear* Result;
  typedef const CDomainOfIntegration * A;
  typedef const ftest  *  B;
  A  di;
  B l;
   
  FormLinear(const basicAC_F0 & args) {
    di= dynamic_cast<A>(CastTo<A>(args[0]));
    assert(di);
    Expression a1=CastTo<B>(args[1]);
    assert(a1);
   // cout << " ---FormLinear: "<< a1 << "  " << typeid(*a1).name() << *a1 <<endl;
    B ll= dynamic_cast<B>(a1);
    assert(ll);
    l = ll->Optimize(currentblock);
    // delete ll; // il ne faut pas detruire car ll peut etre dans une variable
    assert(l);
    throwassert(di && l); 
  }
  bool VF() const { return MaxOp(l) >= last_operatortype;}
 
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
  AnyType operator()(Stack ) const  { return SetAny<Result>(this);}
   operator aType () const { return atype<Result>();}         
  
  static  E_F0 * f(const basicAC_F0 & args) { return new FormLinear(args);}    
  FormLinear(A a,Expression bb) : di(a),l(dynamic_cast<B>(bb)->Optimize(currentblock)) {throwassert(l);}   
  FormLinear operator-() const { return  FormLinear(di,C_F0(TheOperators,"-",C_F0(l,atype<B>())));}
  //  void init(Stack stack) const {}
  
};

class Call_FormLinear: public E_F0mps { public:
  list<C_F0> largs;
  typedef list<C_F0>::const_iterator const_iterator;
  const int N;
  Expression ppfes; 
  Expression *nargs;
  
  Call_FormLinear(Expression * na,Expression  LL, Expression ft) ;
  AnyType operator()(Stack stack) const 
  { InternalError(" bug: no eval of Call_FormLinear ");} 
   operator aType () const { return atype<void>();}         
  
};

class Call_FormBilinear: public E_F0mps { public:
  list<C_F0> largs;
  typedef list<C_F0>::const_iterator const_iterator;
  
  const int N,M;
  Expression *nargs;
  Expression euh,evh; 
  Call_FormBilinear(Expression * na,Expression  LL, Expression fi,Expression fj) ;
  AnyType operator()(Stack stack) const  
  { InternalError(" bug: no eval of Call_FormBilinear ");} 
   operator aType () const { return atype<void>();}         
  
  
};

template<class T>
struct OpCall_FormLinear: public OneOperator {
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =1;

  E_F0 * code(const basicAC_F0 & args) const 
  { 
 Expression * nargs = new Expression[n_name_param];
  args.SetNameParam(n_name_param,name_param,nargs);
  
  return  new Call_FormLinear(nargs,to<const C_args*>(args[0]),to<pfes*>(args[1]));} 
  OpCall_FormLinear() : 
    OneOperator(atype<const Call_FormLinear*>(),atype<const T*>(),atype<pfes*>()) {}
};

template<class T>
struct OpCall_FormLinear2: public OneOperator {
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =1;

  E_F0 * code(const basicAC_F0 & args) const 
  {
  
 Expression * nargs = new Expression[n_name_param];
  args.SetNameParam(n_name_param,name_param,nargs);
    
    Expression p=args[1];
    if ( ! p->EvaluableWithOutStack() ) 
      { CompileError("  a(long,Vh) , The long  must be a constant, and  = 0, sorry");}
    long pv = GetAny<long>((*p)(0));
    if ( pv ) 
      { CompileError("  a(long,Vh) , The long  must be a constant == 0, sorry");}       
    return  new Call_FormLinear(nargs,to<const C_args*>(args[0]),to<pfes*>(args[2]));} 
  OpCall_FormLinear2() : 
    OneOperator(atype<const Call_FormLinear*>(),atype<const T*>(),atype<long>(),atype<pfes*>()) {}
};
template<class T>
struct OpCall_FormBilinear: public OneOperator {
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =9;
  
  E_F0 * code(const basicAC_F0 & args) const 
  { Expression * nargs = new Expression[n_name_param];
  args.SetNameParam(n_name_param,name_param,nargs);
  // cout << " OpCall_FormBilinear " << *args[0].left() << " " << args[0].LeftValue() << endl;
  return  new Call_FormBilinear(nargs,to<const C_args*>(args[0]),to<pfes*>(args[1]),to<pfes*>(args[2]));} 
  OpCall_FormBilinear() : 
    OneOperator(atype<const Call_FormBilinear*>(),atype<const T *>(),atype<pfes*>(),atype<pfes*>()) {}
};

template<class R>  //  to make   x=linearform(x)
struct OpArraytoLinearForm: public OneOperator {
  typedef Call_FormLinear::const_iterator const_iterator;
  
  class Op : public E_F0mps { public:
    const Call_FormLinear *l;
    Expression x;
    AnyType operator()(Stack s)  const ;// {ExecError("Internal Error:to do");} 
    Op(Expression xx,Expression  ll) 
      : l(dynamic_cast<const Call_FormLinear *>(ll)),x(xx) {assert(l);}
     operator aType () const { return atype<KN<R> *>();} 
      
  };
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new Op(to<KN<R>*>(args[0]),args[1]);} 
  OpArraytoLinearForm() : 
    OneOperator(atype<KN<R>*>(),atype<KN<R>*>(),atype<const Call_FormLinear*>()) {}
};

template<class R>  //  to make   A=linearform(x)
struct OpMatrixtoBilinearForm: public OneOperator {
  typedef Call_FormBilinear::const_iterator const_iterator;
  
  class Op : public E_F0mps { public:
    const Call_FormBilinear *b;
    Expression a;
    AnyType operator()(Stack s)  const ;
    Op(Expression aa,Expression  bb) 
      : b(dynamic_cast<const Call_FormBilinear *>(bb)),a(aa) 
    {assert(b && b->nargs);}
     operator aType () const { return atype<Matrice_Creuse<R>  *>();} 

  };
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new Op(to<Matrice_Creuse<R>*>(args[0]),args[1]);} 
  OpMatrixtoBilinearForm() : 
    OneOperator(atype<Matrice_Creuse<R>*>(),atype<Matrice_Creuse<R>*>(),atype<const Call_FormBilinear*>()) {}
};


class IntFunction  : public E_F0mps { public:
  typedef double Result;
  typedef const CDomainOfIntegration * A;
  typedef double  B;
  A  di;
  Expression fonc;    
  IntFunction(const basicAC_F0 & args) {
    di= dynamic_cast<A>(CastTo<A>(args[0]));
    fonc= CastTo<B>(args[1]);
    throwassert(di && fonc); }
  
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
  AnyType operator()(Stack ) const;
  static  E_F0 * f(const basicAC_F0 & args) { return new IntFunction(args);}    
  //  IntFunction(A a,Expression bb) : di(a),fonc(bb) {}   
   operator aType () const { return atype<Result>();}         
  
};


extern Block *currentblock;

class TypeFormOperator: public ForEachType<const C_args*> {
public:
  TypeFormOperator() : ForEachType<const C_args*>(0,0) {}
  void SetArgs(const ListOfId *lid) const {
    SetArgsFormLinear(lid,2);    } 
  
  Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const 
  { return Type_Expr(this,CastTo(c));}
  
  inline  C_F0 Initialization(const Type_Expr & e) const {return C_F0();}
  
};

class TypeFormBilinear: public ForEachType<const FormBilinear*> {
public:
  TypeFormBilinear() : ForEachType<const FormBilinear*>(0,0) {}
  void SetArgs(const ListOfId *lid) const {
    SetArgsFormLinear(lid,2);  
  }
  
  Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const 
  { return Type_Expr(this,CastTo(c));}
  
  
  C_F0 Initialization(const Type_Expr & e) const 
  {  
    // cout << "Initialization " << *e.first << endl;
    return C_F0(); }  // nothing to initialize 
  Type_Expr construct(const Type_Expr & e) const 
  {  
    //cout << "construct " << *e.first << endl;
    return e; } 
  
};

template<bool exec_init,class Problem>
class TypeSolve : public ForEachType<const Problem*>   {
public:
  TypeSolve() : ForEachType<const Problem*>(0,0) {}
  
  void SetArgs(const ListOfId *lid) const {
    SetArgsFormLinear(lid,2);  
    
    
  }  
  Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const 
  {   if (c.left() != atype<const C_args*>())
    CompileError(" Problem  a(...) = invalide type ",c.left());
  const C_args * ca = dynamic_cast<const C_args *>(c.LeftValue());
  Problem * pb=new Problem(ca,*l,top);
  SHOWVERB(cout << "solve:SetParam " << ca << " pb=" << pb << endl);
  
  return Type_Expr(this,pb);        
  
  } 
  
  class SolveInit: public  E_F0 { public:
    const Problem * a;
    AnyType operator()(Stack s)  const {
      a->init(s);       
      return exec_init ? (*a)(s) : Nothing ;
    } 
    SolveInit(const Type_Expr &  te) : a(dynamic_cast<const Problem *>(te.second))
    {  SHOWVERB(cout << "SolveInit " << te.second << endl);
    throwassert(a);}
  };
  
  class SolveDel: public  E_F0 { public:
    const Problem * a;
    SolveDel(const C_F0 & c) : a(dynamic_cast<const Problem *>(c.LeftValue()))
    {   
      SHOWVERB(cout << "SolveDel " << c.left()  << endl);
      throwassert(a);}
    
    AnyType operator()(Stack s)  const { 
      a->destroy(s);
      return Nothing;
    }};
  
  Expression Destroy(const C_F0 & c) const  
  { return new SolveDel(c);}
  
  bool ExistDestroy() const {return true;} 
  
  C_F0 Initialization(const Type_Expr & e) const 
  {  return C_F0( new SolveInit(e) ,atype<void>()); }  
};




class TypeFormLinear: public ForEachType<const FormLinear*> {
public:
  TypeFormLinear() : ForEachType<const FormLinear*>(0,0) {}
  
  void SetArgs(const ListOfId *lid) const {
    SetArgsFormLinear(lid,1);  } 
  
  Type_Expr SetParam(const C_F0 & c,const ListOfId *l,size_t & top) const 
  { return Type_Expr(this,CastTo(c));}
  //  {  return Type_Expr(c.left(),c.LeftValue());  } // 
  
  C_F0 Initialization(const Type_Expr & e) const 
  {  return C_F0(); }  // nothing to initialize 
  
};



template<class K> class Matrice_Creuse  { public:
  CountPointer<FESpace> Uh,Vh;
  pfes  *pUh,*pVh; // pointeur sur la variable stockant FESpace;  
  CountPointer<MatriceCreuse<K> > A;  
  TypeSolveMat typemat;
  void init() {
    A.init(),Uh.init();Vh.init();
    typemat=TypeSolveMat(TypeSolveMat::NONESQUARE);}
  void destroy() {
    A.destroy();
    Uh.destroy();
    Vh.destroy();}   
  Matrice_Creuse( MatriceCreuse<K> * aa,const pfes  *ppUh,const pfes  *ppVh)
    :A(aa),pUh(ppUh),pVh(ppVh),Uh(*ppUh),Vh(*ppVh) {}
};
template<class K>  ostream & operator << (ostream & f,const Matrice_Creuse<K> & A) 
{ if ( !A.A) f << " unset sparce matrix " << endl;
 else f << *A.A ;
 return f;  }

template<class K> class Matrice_Creuse_Transpose  { public:
  Matrice_Creuse<K> * A;
  
  Matrice_Creuse_Transpose(Matrice_Creuse<K> * AA) : A(AA) {assert(A);}
  operator MatriceCreuse<K> & () const {return *A->A;}
  operator Matrice_Creuse<K> * () const {return A;}
};

template<class K> class Matrice_Creuse_inv  { public:
  Matrice_Creuse<K> * A;
  Matrice_Creuse_inv(Matrice_Creuse<K> * AA) : A(AA) {assert(A);}
  operator MatriceCreuse<K> & () const {return *A->A;}
  operator Matrice_Creuse<K> * () const {return A;}
};



template<class K>
class pb2mat : public E_F0 { public:
  typedef Matrice_Creuse<K> *  Result;
  const Problem * pb;
  pb2mat(const basicAC_F0 & args) : pb(dynamic_cast<const Problem *>(args[0].left()))  
  {throwassert(pb);}
  static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<const Problem *>());}
  
  static  E_F0 * f(const basicAC_F0 & args) { return new Plot(args);} 
  
  AnyType operator()(Stack s) const {
    Problem::Data *data= pb->dataptr(stack); 
    throwassert( !!data->A);  
    return  SetAny<Matrice_Creuse<K> * >(&data->A) ;}
  
  
};  



namespace Fem2D {
  
  inline void F_Pi_h(R* v, const R2 & P,const baseFElement & K,int i,const R2 & Phat,void * arg)
  {
    TabFuncArg &tabe(*(TabFuncArg*)arg);
    MeshPoint & mp = *MeshPointStack(tabe.s);
    MeshPointStack(tabe.s)->set(P,Phat,K);
    tabe.eval(v);
    // if (Norme2_2(P-mp.P) > 1e-10) 
    //  cout << " bug?? F_Pi_h " << endl; 
    
  }
  
  inline void FoX_1_Pi_h(R* v, const R2 & P,const baseFElement & K,int i,const R2 & Phat,void * arg)
  {
    TabFuncArg &tabe(*(TabFuncArg*)arg);
    MeshPointStack(tabe.s)->set(P,Phat,K);
    R2 X=tabe.eval_X();
    MeshPointStack(tabe.s)->set(X.x,X.y);
    tabe.eval_2(v);
  }
  
  bool AssembleVarForm(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                       MatriceCreuse<R>  * A,KN<R> * B,const list<C_F0> &largs );
  
  void AssembleLinearForm(Stack stack,const Mesh & Th,const FESpace & Vh,KN<R> * B,const  FormLinear * const l);
  void AssembleBilinearForm(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                            MatriceCreuse<R>  & A, const  FormBilinear * b  );
  void AssembleBC(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                  MatriceCreuse<R>  * A,KN<R> * B,KN<R> * X, const  BC_set<R> * bc , R tgv   );
  void AssembleBC(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                  MatriceCreuse<R>  * A,KN<R> * B,KN<R> * X, const list<C_F0> &largs , R tgv  );
  
  void  Element_rhs(const FElement & Kv,int ie,int label,const LOperaD &Op,R * p,void * stack,KN_<R> & B,bool all);
  void  Element_rhs(const FElement & Kv,const LOperaD &Op,R * p,void * stack,KN_<R> & B);
  void  Element_Op(MatriceElementairePleine<R> & mat,const FElement & Ku,const FElement & Kv,R * p,int ie,int label, void *stack);
  void  Element_Op(MatriceElementaireSymetrique<R> & mat,const FElement & Ku,R * p,int ie,int label, void * stack);
}

template<class R>
AnyType OpArraytoLinearForm<R>::Op::operator()(Stack stack)  const 
{
  
  KN<R> & xx( *GetAny<KN<R> *>((*x)(stack) ));
  pfes  &  pp= *GetAny<pfes * >((*l->ppfes)(stack));
  FESpace * pVh = *pp ;
  FESpace & Vh = *pVh ;
  R tgv= 1e30;
  if (l->nargs[0]) tgv= GetAny<double>((*l->nargs[0])(stack));  
  xx=0; 
  if ( AssembleVarForm(stack,Vh.Th,Vh,Vh,false,0,&xx,l->largs) )
    AssembleBC(stack,Vh.Th,Vh,Vh,false,0,&xx,0,l->largs,tgv);
  // cout << xx.min() << " " << xx.max() << endl;
  /* for (const_iterator i=l->largs.begin();i!=l->largs.end();i++)
     {
     Expression e=i->LeftValue();
     aType r = i->left();
     if (r==atype<const  FormLinear *>()) 
     AssembleLinearForm(s,Vh.Th,Vh,&xx,dynamic_cast<const  FormLinear *>(e));
     else 
     InternalError("OpArraytoLinearForm");
     } */
  return SetAny<KN<R> *>(&xx);
}

template<class R>
void SetSolver(Stack stack,MatriceCreuse<R> & A,const TypeSolveMat *typemat,bool VF,double eps,int NbSpace,int itmax,const OneOperator *precon,int umfpackstrategy, double tgv)
{
  using namespace Fem2D;
  if (typemat->profile)
    {
      MatriceProfile<R> & AA(dynamic_cast<MatriceProfile<R> &>(A));
      throwassert(&AA);
      switch (typemat->t) {
        
      case TypeSolveMat::LU       : AA.typesolver=FactorizationLU; break;
      case TypeSolveMat::CROUT    : AA.typesolver=FactorizationCrout; break;
      case TypeSolveMat::CHOLESKY :  AA.typesolver=FactorizationCholeski; break;
      default:
        cerr << " type resolution " << typemat->t <<" sym=" <<  typemat->profile <<  endl;
        CompileError("type resolution inconnue"); break;       
      }
    }
  else 
    {
      typedef typename MatriceMorse<R>::VirtualSolver VirtualSolver;
      if(verbosity>5) cout << " Matrice morse GC Precond diag" << endl;
      MatriceMorse<R> & AA(dynamic_cast<MatriceMorse<R> &>(A));
      throwassert(&AA);
      //     throwassert(typemat->t==TypeSolveMat::GC);
      switch (typemat->t) {
      case    TypeSolveMat::GC:   
        if (precon)
          AA.SetSolverMaster(static_cast<const VirtualSolver *>(
                                                                new Fem2D::SolveGCPrecon<R>(AA,precon,stack,eps)));
        else 
          AA.SetSolverMaster(static_cast<const VirtualSolver *>(
                                                                new SolveGCDiag<R>(AA,eps)));
        break; 
      case TypeSolveMat::GMRES :
        //        InternalError("GMRES solveur to do");
        if (precon)
          AA.SetSolverMaster(new SolveGMRESPrecon<R>(AA,precon,stack,NbSpace,itmax,eps));
        else 
          AA.SetSolverMaster(new SolveGMRESDiag<R>(AA,NbSpace,itmax,eps));
        break;
#ifdef UMFPACK         
        case TypeSolveMat::UMFpack :
            AA.SetSolverMaster(new SolveUMFPack<R>(AA,umfpackstrategy,tgv));
        break;
           
#endif         
        
      
      default:
      
        if (verbosity >5)
          cout << "  SetSolver:: no  solver by  default " << endl;
        // cerr << " type resolution " << typemat->t << endl;
        //  CompileError("type resolution inconnue"); break;       
      }
      
    }
}
template<class R>
AnyType OpMatrixtoBilinearForm<R>::Op::operator()(Stack stack)  const 
{
  assert(b && b->nargs);// *GetAny<pfes * >
  pfes  * pUh= GetAny<pfes *>((*b->euh)(stack));
  pfes  * pVh= GetAny<pfes *>((*b->evh)(stack));
  const FESpace & Uh =  *(FESpace*) **pUh ;
  const FESpace & Vh =  *(FESpace*) **pVh ;
  MatriceProfile<R> *pmatpf=0;
  bool VF=isVF(b->largs);
//  assert(!VF);
  bool factorize=false;
  long NbSpace = 50; 
  long itermax=0; 
  double eps=1e-6;
  R tgv = 1e30;
  int umfpackstrategy=0;
  TypeSolveMat typemat(  & Uh == & Vh  ? TypeSolveMat::GMRES :TypeSolveMat::NONESQUARE);
  bool initmat=true;
  if (b->nargs[0]) initmat= ! GetAny<bool>((*b->nargs[0])(stack));
  if (b->nargs[1]) typemat= *GetAny<TypeSolveMat *>((*b->nargs[1])(stack));
  if (b->nargs[2]) eps= GetAny<double>((*b->nargs[2])(stack));
  if (b->nargs[4]) NbSpace= GetAny<long>((*b->nargs[4])(stack));
  if (b->nargs[6]) tgv= GetAny<double>((*b->nargs[6])(stack));
  if (b->nargs[7]) factorize= GetAny<bool>((*b->nargs[7])(stack));
  if (b->nargs[8]) umfpackstrategy= GetAny<long>((*b->nargs[8])(stack));
  const OneOperator *precon = 0; //  a changer 
  if ( b->nargs[3])
    {
      const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(b->nargs[3]);
      assert(op);
      precon = op->Find("(",ArrayOfaType(atype<KN<double>* >(),false));
    }
  
  Matrice_Creuse<R> & A( * GetAny<Matrice_Creuse<R>*>((*a)(stack)));
  if  ( (pUh != A.pUh ) || (pVh != A.pVh  || A.typemat.t != typemat.t) )
    { 
      A.Uh.destroy();
      A.Vh.destroy();
    }
  const Mesh & Th = Uh.Th;
  assert( &Uh.Th==&Vh.Th);
  A.typemat = typemat;
  if ( & Uh != A.Uh || & Vh != A.Vh ) 
    { // reconstruct all the matrix
      if (typemat.profile)
        { A.A.master( new MatriceProfile<R>(Vh,VF) ); assert( &Uh == & Vh);}
      else if (typemat.sym )
        { A.A.master( new  MatriceMorse<R>(Vh,typemat.sym,VF) ); 
        assert( &Uh == & Vh);}
      else 
        A.A.master( new  MatriceMorse<R>(Vh,Uh,VF) ); // lines corresponding to test functions 
    }
  *A.A=0; // reset value of the matrix
  if ( AssembleVarForm( stack,Th,Uh,Vh,typemat.sym,A.A,0,b->largs) )
    AssembleBC( stack,Th,Uh,Vh,typemat.sym,A.A,0,0,b->largs,tgv);
  if( factorize ) {
    MatriceProfile<R> * pf = dynamic_cast<MatriceProfile<R> *>((MatriceCreuse<R> *) A.A);
    assert(pf);
    switch (typemat.t) {
    case TypeSolveMat::LU: pf->LU(Abs(eps));break;
    case TypeSolveMat::CROUT: pf->crout(Abs(eps));break;
    case TypeSolveMat::CHOLESKY: pf->cholesky(Abs(eps));break;
    default: ExecError("Sorry no foctorize for this type for matrix"); 
    }
    
  }    
  SetSolver(stack,*A.A,&typemat,VF,eps,NbSpace,itermax,precon,umfpackstrategy,tgv);
  
  return SetAny<Matrice_Creuse<R>  *>(&A);
  
}

class MatrixInterpolation : public OneOperator { public:  

    class Op : public E_F0info { public:
       typedef pfes * A;
       Expression a,b; 
       
       static const int n_name_param =1;
       static basicAC_F0::name_and_type name_param[] ;
        Expression nargs[n_name_param];

       public:
       Op(const basicAC_F0 &  args,Expression aa,Expression bb) : a(aa),b(bb) {
         args.SetNameParam(n_name_param,name_param,nargs);
       } 
    };
    
   MatrixInterpolation() : OneOperator(atype<const MatrixInterpolation::Op *>(),atype<pfes *>(),atype<pfes *>()) {}
  
    E_F0 * code(const basicAC_F0 & args) const 
     { 
       return  new Op(args,t[0]->CastTo(args[0]),
                           t[1]->CastTo(args[1])); 
     }
};


AnyType SetMatrixInterpolation(Stack,Expression ,Expression);

#endif
