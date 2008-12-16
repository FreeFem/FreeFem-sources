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

typedef FEbase<double,v_fes> * pferbase ;
typedef FEbaseArray<double,v_fes> * pferbasearray ;
typedef pair<pferbase,int> pfer ;
typedef pair<pferbasearray,int> pferarray ;

typedef FEbase<Complex,v_fes> * pfecbase ;
typedef FEbaseArray<Complex,v_fes> * pfecbasearray ;
typedef pair<pfecbase,int> pfec ;
typedef pair<pfecbasearray,int> pfecarray ;


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

inline unsigned int GetDiffOp(const MGauche &i, int& lastop) 
   {int op=(i.second% last_operatortype);
     lastop=max(lastop,op) ;
      return 1<<op;}
inline unsigned int GetDiffOp(const MDroit &i, int& lastop)
   {int op=(i.second% last_operatortype);
     lastop=max(lastop,op) ;
      return 1<<op;}
inline unsigned int GetDiffOp(const  pair<MGauche,MDroit> &p, int& lastop) 
{ return GetDiffOp(p.first,lastop)|GetDiffOp(p.second,lastop);}

typedef  const Finconnue  finconnue;
typedef  const Ftest ftest;
typedef  const Foperator foperator;

Expression IsFebaseArray(Expression f);

void SetArgsFormLinear(const ListOfId *lid,int ordre);


inline ostream & operator<<(ostream & f,const  TypeSolveMat & tm)
{
  switch(tm.t) {
   case TypeSolveMat::NONESQUARE:  f << "No Square (Sparse Morse)"; break;
   case TypeSolveMat::LU:  f << "LU (Skyline)"; break;
   case TypeSolveMat::CROUT:  f << "CROUT (Skyline)"; break;
   case TypeSolveMat::CHOLESKY:  f << "CHOLESKY (Skyline)"; break;
   case TypeSolveMat::GC:  f << "CG (Sparse Morse)"; break;
   case TypeSolveMat::GMRES:  f << "GMRES (Sparse Morse)"; break;
   case TypeSolveMat::SparseSolver:  f << "SparseSolver (Sparse Morse)"; break;
   default: f << "Unknown  bug???";
   }
  return f;
}


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
    ffassert(n==2);    
    if (args[0].left() == atype<const C_args *>())
      {
        const C_args * a = dynamic_cast<const C_args *>(args[0].LeftValue());
        ffassert(a);
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
    ffassert(n==1);
    aType tF=atype<Result>();
    ffassert(args[0].left() == tF);
    Result f = dynamic_cast<Result>(args[0].LeftValue());
    ffassert(f);
    // F mf = -*f;
    F * rf=new F(-*f);
    return  rf;
  }
    operator aType () const { return atype<Result>();} 
    
    
};

//template<class RR=double>
class BC_set : public E_F0mps { public:
   bool  complextype;
  typedef const BC_set* Result;
  vector<Expression> on;
  vector<pair<int,Expression> > bc; //  n¡ de l'inconnue+ valeur
  BC_set(  const basicAC_F0 & args) :on(args.size()){
    int n = args.size();      
    ffassert(args.named_parameter);
    AC_F0::const_iterator ii=args.named_parameter->begin();
    AC_F0::const_iterator ie=args.named_parameter->end();
    bc.resize(args.named_parameter->size());
    complextype=false;
    for (int kk=0;ii!=ie;kk++,ii++)
     {
      if( ! BCastTo<double>(ii->second)) 
       complextype = true; 
     } 
    ii=args.named_parameter->begin();
    for (int kk=0;ii!=ie;kk++,ii++)
      { //
        C_F0 x=Find(ii->first);
        if (x.left() != atype<const finconnue *>())
          CompileError("We expected an unkown  u=... of the problem");
        const finconnue * uu = dynamic_cast<const finconnue *>(x.LeftValue());
        ffassert(uu);
        const MGauche *ui=uu->simple();
        ffassert(ui && ui->second == op_id);
        // cout << ii->first << " n¡" << ui->first <<   " = ? " << endl;
        if (complextype)
        bc[kk]= make_pair(ui->first,CastTo<Complex>(ii->second));
        else
        bc[kk]= make_pair(ui->first,CastTo<double>(ii->second));
        ii->second;
      }     
    for (int i=0;i<n;i++)
      on[i]=CastTo<long>(args[i]);
  }

  /*
  // ajout modif FH  mai 2007 XXXXXXXXXXXXX
  void mapping(C_F0 (*f)(const C_F0 &)) {
      for ( vector<pair<int,Expression> >::iterator k=bc.begin();k!=bc.end();k++)
	  k->second=(*f)(k->second) ;}
  // fin ajout
   */ 
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
  enum typeofkind  { int2d=0, int1d=1, intalledges=2,intallVFedges=3, int3d = 4, intallfaces= 5 } ; //3d
  typeofkind  kind; //  0 
  int d; // 3d
  typedef const CDomainOfIntegration* Result;
  Expression Th; 
  vector<Expression> what;
  CDomainOfIntegration( const basicAC_F0 & args,typeofkind b=int2d,int ddim=2) // 3d
    :kind(b),d(ddim), what(args.size()-1)
     
  {
    args.SetNameParam(n_name_param,name_param,nargs);
    if(d==2) // 3d
      Th=CastTo<pmesh>(args[0]);
    else if(d==3)
      Th=CastTo<pmesh3>(args[0]);
    else ffassert(0); // a faire 
    int n=args.size();

    for (int i=1;i<n;i++)
      what[i-1]=CastTo<long>(args[i]); 
    // cout << " CDomainOfIntegration " << this << endl;       
  }
  static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh>(), true);} // all type
  AnyType operator()(Stack ) const  { return SetAny<const CDomainOfIntegration *>(this);}
  operator aType () const { return atype<const CDomainOfIntegration *>();}         
  
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args);} 
  const Fem2D::QuadratureFormular & FIT(Stack) const ;   
  const Fem2D::QuadratureFormular1d & FIE(Stack) const ;  
  const Fem2D::GQuadratureFormular<R3> & FIV(Stack) const ;  // 3d
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

class CDomainOfIntegrationVFEdges: public CDomainOfIntegration { 
public:
  CDomainOfIntegrationVFEdges( const basicAC_F0 & args) :CDomainOfIntegration(args,intallVFedges) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,intallVFedges);}    
};
// add for the 3d case .. // 3d 
class CDomainOfIntegration3d: public CDomainOfIntegration { 
public:
  CDomainOfIntegration3d( const basicAC_F0 & args) :CDomainOfIntegration(args,int3d,3) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int3d,3);}
  static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh3>(), true);} // all type    
};

class CDomainOfIntegrationBorder3d: public CDomainOfIntegration { 
public:
  CDomainOfIntegrationBorder3d( const basicAC_F0 & args) :CDomainOfIntegration(args,int2d,3) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,int2d,3);}    
  static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh3>(), true);} // all type
};

class CDomainOfIntegrationAllFaces: public CDomainOfIntegration { 
public:
  CDomainOfIntegrationAllFaces( const basicAC_F0 & args) :CDomainOfIntegration(args,intallfaces,3) {}
  static  E_F0 * f(const basicAC_F0 & args) { return new CDomainOfIntegration(args,intallfaces,3);}    
  static  ArrayOfaType  typeargs() {  return ArrayOfaType(atype<pmesh3>(), true);} // all type
};

// end add





// hack build template 
template<class T> struct CadnaType{
   typedef T  Scalaire;
 }; 
#ifdef  HAVE_CADNA
#include <cadnafree.h>
// specialisation 
template<> struct CadnaType<complex<double> >{
   typedef complex<double_st> Scalaire;
 };
template<> struct CadnaType<double> {
   typedef double_st Scalaire;
 };
inline double_st conj(const double_st &x){ return x;};
inline complex<double_st> conj(const complex<double_st> &x){ return complex<double_st>(x.real(),-x.imag());}
 
 inline double norm(complex<double_st> x){return x.real()*x.real()+x.imag()*x.imag();} 
 inline double norm(double_st x){return x*x;} 
inline int cestac(const complex<double_st> & z) 
{return min(cestac(z.real()),cestac(z.imag()));}
#endif

class Problem  : public Polymorphic { 
  //  typedef double R;
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =13; // modi FH oct 2005 add tol_pivot 02/ 2007 add nbiter   
  int Nitem,Mitem;
  const int dim; 
public:
 template<class FESpace>   
  struct Data {
    typedef typename FESpace::Mesh Mesh;
    const Mesh * pTh; 
    CountPointer<const FESpace> Uh;
    CountPointer<const FESpace> Vh;
    CountPointer<MatriceCreuse<double> > AR;   
    CountPointer<MatriceCreuse<Complex> > AC; 
    typedef  CadnaType<double>::Scalaire double_st;  
    typedef  CadnaType<complex<double> >::Scalaire cmplx_st;  
    MatriceCreuse<double_st>  * AcadnaR;
    MatriceCreuse<cmplx_st>  * AcadnaC;

    void init()  {pTh=0; AcadnaR=0;AcadnaC=0; Uh.init(),Vh.init();AR.init();AC.init();}
    void destroy() { 
      pTh=0;
      Uh.destroy(); 
      Vh.destroy();
      AR.destroy();       
      AC.destroy();
      
      if(AcadnaR) AcadnaR->destroy();
      if(AcadnaC) AcadnaC->destroy();
      }        
  } ;
  const OneOperator *precon;
  
  C_args *op; // the list of all operator 
  mutable vector<Expression> var; // list des var pour les solutions et test 
  bool complextype,VF;

  //   Expression noinit,type,epsilon;
  Expression nargs[n_name_param];
  
  const size_t  offset;         
  Problem(const C_args * ca,const ListOfId &l,size_t & top) ;
  static ArrayOfaType  typeargs() { return ArrayOfaType(true);}// all type
  
  Data<FESpace> * dataptr (Stack stack) const {return   (Data<FESpace> *) (void *) (((char *) stack)+offset);}
  Data<FESpace3> * dataptr3 (Stack stack) const {return   (Data<FESpace3> *) (void *) (((char *) stack)+offset);}
  void init(Stack stack) const  {
    // cout << " init  " << (char *) dataptr(stack) - (char*) stack  << " " << offset <<  endl; 
    if(dim==2)
	dataptr(stack)->init();
    else 
	dataptr3(stack)->init();
      }
  void destroy(Stack stack)  const  {
      if(dim==2) dataptr(stack)->destroy();
      else dataptr3(stack)->destroy();
  }

  template<class R,class FESpace,class v_fes>
  AnyType eval(Stack stack,Data<FESpace> * data,CountPointer<MatriceCreuse<R> > & dataA,
  MatriceCreuse< typename CadnaType<R>::Scalaire >  * & dataCadna) const;

  AnyType operator()(Stack stack) const
  {
   if(dim==2)
      {
     Data<FESpace> *data= dataptr(stack);
    if (complextype) 
	return eval<Complex,FESpace,v_fes>(stack,data,data->AC,data->AcadnaC);
    else
     return eval<double,FESpace,v_fes>(stack,data,data->AR,data->AcadnaR);
      }
    
    else if(dim==3)
      {
	Data<FESpace3> *data= dataptr3(stack);
      if (complextype) 
	  return eval<Complex,FESpace3,v_fes3>(stack,data,data->AC,data->AcadnaC);
      else
	  return eval<double,FESpace3,v_fes3>(stack,data,data->AR,data->AcadnaR);
      }
    else ffassert(0);
  }
    
  
  
  bool Empty() const {return false;}     
  size_t nbitem() const { return Nitem;}     
};


class Solve : public Problem { public:
  // just a problem with implicit solve 
  Solve(const C_args * ca,const ListOfId &l,size_t & top) 
    : Problem(new C_args(*ca),l,top) {} 
};

class FormBilinear : public E_F0mps { public:
  typedef const FormBilinear* Result;
  typedef const CDomainOfIntegration * A;
  typedef const foperator  *  B;
  A  di;
  Foperator * b;    
  FormBilinear(const basicAC_F0 & args) {
   di= dynamic_cast<A>(CastTo<A>(args[0]));
   B  bb= dynamic_cast<B>(CastTo<B>(args[1]));
  // b = bb->Optimize(currentblock); // FH1004
  b=new Foperator(*bb); // FH1004  no optimisation here because we don't the type of the bilinear form here.
  //  the opimisation is done after in FieldOfForm routine 
  // to find if the form is real or complex 
  
  // delete bb; il ne faut pas detruire .. car bb peut etre dans une variable 
    ffassert(di && b); }
  
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
  AnyType operator()(Stack ) const { return SetAny<Result>(this);}
   operator aType () const { return atype<Result>();}         

  static  E_F0 * f(const basicAC_F0 & args) { return new FormBilinear(args);}  
  FormBilinear(A a,Expression bb) : di(a),b(new Foperator(*dynamic_cast<B>(bb))/*->Optimize(currentblock) FH1004 */) 
  {ffassert(b);}   
  FormBilinear operator-() const { return  FormBilinear(di,C_F0(TheOperators,"-",C_F0(b,atype<B>())));}
  bool VF() const { return MaxOp(b) >= last_operatortype;}
  int dim() const {return di->d;}
  FormBilinear(const FormBilinear & fb) : di(fb.di),b(new Foperator(*fb.b) ) {}
  //  void init(Stack stack) const {}
};


//template<class v_fes>
class FormLinear : public E_F0mps { public:
  typedef const FormLinear* Result;
  typedef const CDomainOfIntegration * A;
  typedef const ftest  *  B;
  A  di;
  Ftest * l;
   
  FormLinear(const basicAC_F0 & args) {
    di= dynamic_cast<A>(CastTo<A>(args[0]));
    assert(di);
    Expression a1=CastTo<B>(args[1]);
    assert(a1);
   // cout << " ---FormLinear: "<< a1 << "  " << typeid(*a1).name() << *a1 <<endl;
    B ll= dynamic_cast<B>(a1);
    assert(ll);
    l = new Ftest(*ll); // FH1004 ->Optimize(currentblock);  same as bilinear
    // delete ll; // il ne faut pas detruire car ll peut etre dans une variable
    assert(l);
    ffassert(di && l); 
  }
  bool VF() const { return MaxOp(l) >= last_operatortype;}
 
  static ArrayOfaType  typeargs() { return ArrayOfaType(atype<A>(),atype<B>());}// all type
  AnyType operator()(Stack ) const  { return SetAny<Result>(this);}
   operator aType () const { return atype<Result>();}         
  int dim() const {return di->d;}
  static  E_F0 * f(const basicAC_F0 & args) { return new FormLinear(args);}    
  FormLinear(A a,Expression bb) : di(a),l(new Ftest(*dynamic_cast<B>(bb))/*->Optimize(currentblock) FH1004 */) {ffassert(l);}   
  FormLinear operator-() const { return  FormLinear(di,C_F0(TheOperators,"-",C_F0(l,atype<B>())));}
  //  void init(Stack stack) const {}
  FormLinear(const FormLinear & fb) : di(fb.di),l(new Ftest(*fb.l) ) {}
  
};

template<class VFES>
class Call_FormLinear: public E_F0mps 
{ 
public:
  const int d;
  list<C_F0> largs;
  Expression *nargs;
  typedef list<C_F0>::const_iterator const_iterator;
  const int N;
  Expression ppfes; 
  
  Call_FormLinear(int dd,Expression * na,Expression  LL, Expression ft) ;
  AnyType operator()(Stack stack) const 
  { InternalError(" bug: no eval of Call_FormLinear ");} 
   operator aType () const { return atype<void>();}         
  
};

template<class VFES>
class Call_FormBilinear: public E_F0mps 
{
public:
  const int d;
  Expression *nargs;
  list<C_F0> largs;
  typedef list<C_F0>::const_iterator const_iterator;
  
  const int N,M;
  Expression euh,evh; 
  Call_FormBilinear(int dd,Expression * na,Expression  LL, Expression fi,Expression fj) ;
  AnyType operator()(Stack stack) const  
  { InternalError(" bug: no eval of Call_FormBilinear ");} 
   operator aType () const { return atype<void>();}         
    
};


struct OpCall_FormLinear_np {
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =1;
  
};

struct OpCall_FormBilinear_np {
  static basicAC_F0::name_and_type name_param[] ;
  static const int n_name_param =12; // 9-> 11 FH 31/10/2005  11->12 nbiter 02/2007
};

template<class T,class v_fes>
struct OpCall_FormLinear
  : public OneOperator,
    public OpCall_FormLinear_np
{
  typedef v_fes *pfes;
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    Expression * nargs = new Expression[n_name_param];
    args.SetNameParam(n_name_param,name_param,nargs);
    
    return  new Call_FormLinear<v_fes>(v_fes::d,nargs,to<const C_args*>(args[0]),to<pfes*>(args[1]));} 
  OpCall_FormLinear() : 
    OneOperator(atype<const Call_FormLinear<v_fes>*>(),atype<const T*>(),atype<pfes*>()) {}
};


template<class T,class v_fes>
struct OpCall_FormLinear2
  : public OneOperator,
    public OpCall_FormLinear_np
{
  static const int d=v_fes::d;
  typedef v_fes *pfes;
  E_F0 * code(const basicAC_F0 & args) const 
  {    
    Expression * nargs = new Expression[this->n_name_param];
    args.SetNameParam(this->n_name_param,this->name_param,nargs);
    
    Expression p=args[1];
    if ( ! p->EvaluableWithOutStack() ) 
      { CompileError("  a(long,Vh) , The long  must be a constant, and  = 0, sorry");}
    long pv = GetAny<long>((*p)(NullStack));
    if ( pv ) 
      { CompileError("  a(long,Vh) , The long  must be a constant == 0, sorry");}       
    return  new Call_FormLinear<v_fes>(v_fes::d,nargs,to<const C_args*>(args[0]),to<pfes*>(args[2]));} 
  OpCall_FormLinear2() : 
    OneOperator(atype<const Call_FormLinear<v_fes>*>(),atype<const T*>(),atype<long>(),atype<pfes*>()) {}
};

template<class T,class v_fes>
struct OpCall_FormBilinear
  : public OneOperator ,
  OpCall_FormBilinear_np
{
  typedef v_fes *pfes;
  static const int d=v_fes::d;
  
  E_F0 * code(const basicAC_F0 & args) const 
  { Expression * nargs = new Expression[n_name_param];
    args.SetNameParam(n_name_param,name_param,nargs);
    // cout << " OpCall_FormBilinear " << *args[0].left() << " " << args[0].LeftValue() << endl;
    return  new Call_FormBilinear<v_fes>(v_fes::d,nargs,to<const C_args*>(args[0]),to<pfes*>(args[1]),to<pfes*>(args[2]));} 
  OpCall_FormBilinear() : 
    OneOperator(atype<const Call_FormBilinear<v_fes>*>(),atype<const T *>(),atype<pfes*>(),atype<pfes*>()) {}
};





bool FieldOfForm( list<C_F0> & largs ,bool complextype);
template<class A>  struct IsComplexType { static const bool value=false;};
template<>  struct IsComplexType<Complex> { static const bool value=true;};

template<class R,class v_fes>  //  to make   x=linearform(x)
struct OpArraytoLinearForm
  : public OneOperator
{
  typedef typename Call_FormLinear<v_fes>::const_iterator const_iterator;
  const bool isptr;
  const bool init;
  const bool zero;
  class Op : public E_F0mps 
  {
  public:
    Call_FormLinear<v_fes> *l;
    Expression x;
    const  bool isptr; 
    const  bool init; 
    const  bool zero;
    AnyType operator()(Stack s)  const ;
    Op(Expression xx,Expression  ll,bool isptrr,bool initt,bool zzero) 
      : l(new Call_FormLinear<v_fes>(*dynamic_cast<const Call_FormLinear<v_fes> *>(ll))),
        x(xx),
        isptr(isptrr),init(initt),zero(zzero)
        {assert(l);FieldOfForm(l->largs,IsComplexType<R>::value); }
     operator aType () const { return atype<KN<R> *>();} 
      
  };
  E_F0 * code(const basicAC_F0 & args) const 
  { if(isptr) return new Op(to<KN<R> *>(args[0]),args[1],isptr,init,zero);
    else      return new Op(to<KN_<R> >(args[0]),args[1],isptr,init,zero);}
   
  
  //  OpArraytoLinearForm(const basicForEachType * tt) : 
  //    OneOperator(atype<KN_<R> >(),tt,atype<const Call_FormLinear*>()),init(false),isptr(false) {}
  
  OpArraytoLinearForm(const basicForEachType * tt,bool isptrr, bool initt,bool zzero=1) : 
    OneOperator(atype<KN_<R> >(),tt,atype<const Call_FormLinear<v_fes>*>()),
    isptr(isptrr), init(initt),zero(zzero) {}
  
};


template<class R,class v_fes>  //  to make   A=linearform(x)
struct OpMatrixtoBilinearForm
  : public OneOperator
{
  typedef typename Call_FormBilinear<v_fes>::const_iterator const_iterator;
  int init; 
  class Op : public E_F0mps {
  public:
    Call_FormBilinear<v_fes> *b;
    Expression a;
    int init;
    AnyType operator()(Stack s)  const ;    

    Op(Expression aa,Expression  bb,int initt) 
      : b(new Call_FormBilinear<v_fes>(* dynamic_cast<const Call_FormBilinear<v_fes> *>(bb))),a(aa),init(initt) 
  {assert(b && b->nargs);FieldOfForm(b->largs,IsComplexType<R>::value)  ;}
  operator aType () const { return atype<Matrice_Creuse<R>  *>();} 
    
  };
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new Op(to<Matrice_Creuse<R>*>(args[0]),args[1],init);} 
  OpMatrixtoBilinearForm(int initt=0) : 
    OneOperator(atype<Matrice_Creuse<R>*>(),atype<Matrice_Creuse<R>*>(),atype<const Call_FormBilinear<v_fes>*>()),
    init(initt)
 {}
};

template<class R>
class IntFunction  : public E_F0mps { public:
  typedef R Result;
  typedef const CDomainOfIntegration * A;
  typedef R  B;
  A  di;
  Expression fonc;    
  IntFunction(const basicAC_F0 & args) {
    di= dynamic_cast<A>(CastTo<A>(args[0]));
    fonc= CastTo<B>(args[1]);
    ffassert(di && fonc); }
  
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
    CompileError(" Problem  a(...) = invalid type ",c.left());
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
    ffassert(a);}
  };
  
  class SolveDel: public  E_F0 { public:
    const Problem * a;
    SolveDel(const C_F0 & c) : a(dynamic_cast<const Problem *>(c.LeftValue()))
    {   
      SHOWVERB(cout << "SolveDel " << c.left()  << endl);
      ffassert(a);}
    
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



template<class K> class Matrice_Creuse  { 
  //  CountPointer<FESpace> Uh,Vh;
  //pfes  *pUh,*pVh; // pointeur sur la variable stockant FESpace;  
public:
  UniqueffId Uh,Vh; // pour la reconstruction 
  //  const void * pUh,pVh; //  pointeur pour la reconstruction 
  CountPointer<MatriceCreuse<K> > A;  
  TypeSolveMat typemat;
  void init() {
    A.init();Uh.init();Vh.init();
    typemat=TypeSolveMat(TypeSolveMat::NONESQUARE);}
  Matrice_Creuse() { init();}
  void destroy() {
    A.destroy();
    //    Uh.destroy();
    //Vh.destroy();
  }   
  Matrice_Creuse( MatriceCreuse<K> * aa)//,const pfes  *ppUh,const pfes  *ppVh)
    :A(aa){}//,pUh(ppUh),pVh(ppVh),Uh(*ppUh),Vh(*ppVh) {}
  Matrice_Creuse( MatriceCreuse<K> * aa,const UniqueffId *pUh,const UniqueffId *pVh)//,const pfes  *ppUh,const pfes  *ppVh)
    :A(aa),Uh(*pUh),Vh(*pVh) {}//,pUh(ppUh),pVh(ppVh),Uh(*ppUh),Vh(*ppVh) {}
  long N() const {return  A ? A->n : 0;}
  long M() const { return A ? A->m : 0;}
  
};

template<class K> class Matrice_Creuse_Transpose;

 template<class KA,class KB>   class Matrix_Prod { public:  
  Matrice_Creuse<KA> *A;
  Matrice_Creuse<KB> *B;
  bool ta,tb;
  Matrix_Prod(Matrice_Creuse<KA> *AA,Matrice_Creuse<KB> *BB) : A(AA),B(BB),ta(false),tb(false) {assert(AA && BB);}
  Matrix_Prod(Matrice_Creuse_Transpose<KA> AA,Matrice_Creuse<KB> *BB)           : A(AA),B(BB),ta(true),tb(false) {assert(AA && BB);}
  Matrix_Prod(Matrice_Creuse<KA> *AA,Matrice_Creuse_Transpose<KB> BB)           : A(AA),B(BB),ta(false),tb(true) {assert(AA && BB);}
  Matrix_Prod(Matrice_Creuse_Transpose<KA> AA,Matrice_Creuse_Transpose<KB> BB) : A(AA),B(BB),ta(true),tb(true) {assert(AA && BB);}
 };
 
template<class K>  ostream & operator << (ostream & f,const Matrice_Creuse<K> & A) 
{ if ( !A.A) f << " unset sparse matrix " << endl;
 else f << *A.A ;
 return f;  }

template<class K>  istream & operator >> (istream & f,Matrice_Creuse<K> & A) 
{ 
    if ( WhichMatrix(f)== 2   )
    {
      //	A.pUh=0;
      //A.pVh=0; 
	A.A.master(new MatriceMorse<K>(f));
	A.typemat=(A.A->n == A.A->m) ? TypeSolveMat(TypeSolveMat::GMRES) : TypeSolveMat(TypeSolveMat::NONESQUARE); //  none square matrice (morse)
	
    }
    else {  
	cerr << " unkwon type of matrix " << endl;
	ExecError("Erreur read matrix ");	
	A.A =0; }
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





namespace Fem2D {
  
  inline void F_Pi_h(R* v, const R2 & P,const baseFElement & K,int i,const R2 & Phat,void * arg)
  {
    TabFuncArg &tabe(*(TabFuncArg*)arg);
    //MeshPoint & mp = *MeshPointStack(tabe.s);
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
  
  template<class R,typename MC,class FESpace >  bool AssembleVarForm(Stack stack,const typename FESpace::Mesh & Th,
								     const FESpace & Uh,const FESpace & Vh,bool sym,
								     MC  * A,KN_<R> * B,const list<C_F0> &largs );

  template<class R,class FESpace>   void AssembleBC(Stack stack,const typename FESpace::Mesh & Th,
						    const FESpace & Uh,const FESpace & Vh,bool sym,
						    MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X,
						    const list<C_F0> &largs , double tgv  );
  
  
  template<class R>   void AssembleLinearForm(Stack stack,const Mesh & Th,const FESpace & Vh,KN_<R> * B,const  FormLinear * const l);
  
  //  template<class R>   void AssembleBC(Stack stack,const Mesh3 & Th,const FESpace3 & Uh,const FESpace3 & Vh,bool sym,
  //			      MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const list<C_F0> &largs , double tgv  );
  
  
  template<class R>   void AssembleLinearForm(Stack stack,const Mesh3 & Th,const FESpace3 & Vh,KN_<R> * B,const  FormLinear * const l);
  template<class R>   void AssembleBC(Stack stack,const Mesh & Th3,const FESpace & Uh3,const FESpace & Vh3,bool sym,
				      MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc , double tgv   );
  
  template<class R>   void  Element_rhs(const FElement3 & Kv,int ie,int label,const LOperaD &Op,double * p,void * stack,KN_<R> & B,bool all);
  template<class R>   void  Element_rhs(const FElement3 & Kv,const LOperaD &Op,double * p,void * stack,KN_<R> & B);
  template<class R>   void  Element_Op(MatriceElementairePleine<R,FESpace3> & mat,const FElement3 & Ku,const FElement3 & Kv,double * p,int ie,int label, void *stack);
  template<class R>   void  Element_Op(MatriceElementaireSymetrique<R,FESpace3> & mat,const FElement3 & Ku,double * p,int ie,int label, void * stack);
  
  template<class R,class FESpace>
  void AssembleBC(Stack stack,const typename FESpace::Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                  MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const list<C_F0> &largs , double tgv  );
  
  // fin 3d 
  
  template<class R>   void  Element_rhs(const FElement & Kv,int ie,int label,const LOperaD &Op,double * p,void * stack,KN_<R> & B,bool all);
  template<class R>   void  Element_rhs(const FElement & Kv,const LOperaD &Op,double * p,void * stack,KN_<R> & B);
  template<class R>   void  Element_Op(MatriceElementairePleine<R,FESpace> & mat,const FElement & Ku,const FElement & Kv,double * p,int ie,int label, void *stack);
  template<class R>   void  Element_Op(MatriceElementaireSymetrique<R,FESpace> & mat,const FElement & Ku,double * p,int ie,int label, void * stack);
  
/*template<class R>   void AssembleBilinearForm(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
                            MatriceCreuse<R>  & A, const  FormBilinear * b  );
*/ // --------- FH 120105                           

  //template<class R>   void AssembleBC(Stack stack,const Mesh & Th,const FESpace & Uh,const FESpace & Vh,bool sym,
  //                MatriceCreuse<R>  * A,KN_<R> * B,KN_<R> * X, const  BC_set * bc , double tgv   );
  
  
  //------



}


template<class R,class v_fes>
AnyType OpArraytoLinearForm<R,v_fes>::Op::operator()(Stack stack)  const 
{
  typedef v_fes *pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;  

  pfes  &  pp= *GetAny<pfes * >((*l->ppfes)(stack));
  FESpace * pVh = *pp ;
  FESpace & Vh = *pVh ;
  double tgv= 1e30;
  if (l->nargs[0]) tgv= GetAny<double>((*l->nargs[0])(stack));  

  KN<R> *px=0;
  if(isptr)
    {
     px = GetAny<KN<R> * >((*x)(stack) );
     if(init ) 
       px->init(Vh.NbOfDF); 
     ffassert(px->N() == Vh.NbOfDF); 
   }
  KN_<R>  xx( px ? *(KN_<R> *) px : GetAny<KN_<R> >((*x)(stack) ));
  if(zero)
  xx=R(); 
  if ( AssembleVarForm<R,MatriceCreuse<R>,FESpace >(stack,Vh.Th,Vh,Vh,false,0,&xx,l->largs) )
    AssembleBC<R,FESpace>(stack,Vh.Th,Vh,Vh,false,0,&xx,0,l->largs,tgv);
  return SetAny<KN_<R> >(xx);
}

template<class R>
void SetSolver(Stack stack,MatriceCreuse<R> & A,const TypeSolveMat *typemat,bool VF,double eps,int NbSpace,int itmax,const OneOperator * const precon,int umfpackstrategy, double tgv,
double tol_pivot,double tol_pivot_sym)
{ 
  using namespace Fem2D;
  if (typemat->profile)
    {
      MatriceProfile<R> & AA(dynamic_cast<MatriceProfile<R> &>(A));
      ffassert(&AA);
      switch (typemat->t) {
        
      case TypeSolveMat::LU       : AA.typesolver=FactorizationLU; break;
      case TypeSolveMat::CROUT    : AA.typesolver=FactorizationCrout; break;
      case TypeSolveMat::CHOLESKY :  AA.typesolver=FactorizationCholeski; break;
      default:
        cerr << " type resolution " << typemat->t <<" sym=" <<  typemat->profile <<  endl;
        CompileError("type resolution unknown"); break;       
      }
    }
  else 
    {
      typedef typename MatriceMorse<R>::VirtualSolver VirtualSolver;
      if(verbosity>5) cout << " Morse matrix GC Precond diag" << endl;
      MatriceMorse<R> & AA(dynamic_cast<MatriceMorse<R> &>(A));
      ffassert(&AA);
      //     ffassert(typemat->t==TypeSolveMat::GC);
      switch (typemat->t) {
      case    TypeSolveMat::GC:   
        if (precon)
          AA.SetSolverMaster(static_cast<const VirtualSolver *>(
                                                                new Fem2D::SolveGCPrecon<R>(AA,precon,stack,itmax,eps)));
        else 
          AA.SetSolverMaster(static_cast<const VirtualSolver *>(
                                                                new SolveGCDiag<R>(AA,itmax,eps)));
        break; 
      case TypeSolveMat::GMRES :
        //        InternalError("GMRES solveur to do");
        if (precon)
          AA.SetSolverMaster(new SolveGMRESPrecon<R>(AA,precon,stack,NbSpace,itmax,eps));
        else 
          AA.SetSolverMaster(new SolveGMRESDiag<R>(AA,NbSpace,itmax,eps));
        break;
//#ifdef HAVE_LIBUMFPACK         
        case TypeSolveMat::SparseSolver :
	    AA.SetSolverMaster(DefSparseSolver<R>::Build(&AA,umfpackstrategy,tgv,eps,tol_pivot,tol_pivot_sym,NbSpace,itmax,(void *)precon,stack )); 
         //   AA.SetSolverMaster(new SolveUMFPack<R>(AA,umfpackstrategy,tgv,eps,tol_pivot,tol_pivot_sym));
        break;
           
//#endif         
        
      
      default:
      
        if (verbosity >5)
          cout << "  SetSolver:: no  default solver " << endl;
        // cerr << " type resolution " << typemat->t << endl;
        //  CompileError("type resolution inconnue"); break;       
      }
      
    }
}

template<class R,class v_fes>
AnyType OpMatrixtoBilinearForm<R,v_fes>::Op::operator()(Stack stack)  const 
{
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;  

  assert(b && b->nargs);// *GetAny<pfes * >
  pfes  * pUh= GetAny<pfes *>((*b->euh)(stack));
  pfes  * pVh= GetAny<pfes *>((*b->evh)(stack));
  const FESpace & Uh =  *(FESpace*) **pUh ;
  const FESpace & Vh =  *(FESpace*) **pVh ;
  bool A_is_square= Uh.NbOfDF == Vh.NbOfDF ;

  // MatriceProfile<R> *pmatpf=0;
  bool VF=isVF(b->largs);
//  assert(!VF);
  bool factorize=false;
  long NbSpace = 50; 
  long itermax=0; 
  double eps=1e-6;
  double tgv = 1e30;
  int umfpackstrategy=0;
  double tol_pivot=-1;
  double tol_pivot_sym=-1;
  TypeSolveMat typemat(  & Uh == & Vh  ? TypeSolveMat::GMRES :TypeSolveMat::NONESQUARE);
  bool initmat=true;
  if (b->nargs[0]) initmat= ! GetAny<bool>((*b->nargs[0])(stack));
  if (b->nargs[1]) typemat= *GetAny<TypeSolveMat *>((*b->nargs[1])(stack));
  if (b->nargs[2]) eps= GetAny<double>((*b->nargs[2])(stack));
  if (b->nargs[4]) NbSpace= GetAny<long>((*b->nargs[4])(stack));
  if (b->nargs[6]) tgv= GetAny<double>((*b->nargs[6])(stack));
  if (b->nargs[7]) factorize= GetAny<bool>((*b->nargs[7])(stack));
  if (b->nargs[8]) umfpackstrategy= GetAny<long>((*b->nargs[8])(stack));
  if (b->nargs[9]) tol_pivot= GetAny<double>((*b->nargs[9])(stack));
  if (b->nargs[10]) tol_pivot_sym= GetAny<double>((*b->nargs[10])(stack));
  if (b->nargs[11]) itermax= GetAny<long>((*b->nargs[11])(stack));
  
  if (! A_is_square && typemat != TypeSolveMat::NONESQUARE) 
   {
     cout << " -- Error the solver << "<< typemat <<"  is set  on rectangular matrix  " << endl;
     ExecError("A solver is set on a none square matrix!");
    typemat= TypeSolveMat::NONESQUARE;
   }
  const OneOperator *precon = 0; //  a changer 
  if ( b->nargs[3])
    {
      const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(b->nargs[3]);
      ffassert(op);
      precon = op->Find("(",ArrayOfaType(atype<KN<double>* >(),false));
    }
  //  for the gestion of the PTR. 
  WhereStackOfPtr2Free(stack)=new StackOfPtr2Free(stack);// FH aout 2007 
  
  Matrice_Creuse<R> & A( * GetAny<Matrice_Creuse<R>*>((*a)(stack)));
  if(init) A.init(); //  
  /*  if  ( (pUh != A.pUh ) || (pVh != A.pVh  || A.typemat.t != typemat.t) )
    { 
      A.Uh.destroy();
      A.Vh.destroy();
      }*/
  const Mesh & Th = Uh.Th;
  bool same=isSameMesh(b->largs,&Uh.Th,&Vh.Th,stack);     
  if ( same)
   {
     A.typemat = typemat;
     if ( A.Uh != Uh  || A.Vh != Vh ) 
       { // reconstruct all the matrix
	 A.A=0; // to delete  old  matrix ADD FH 16112005 
	 A.Uh=Uh;
	 A.Vh=Vh;
	 if (typemat.profile)
	   { A.A.master( new MatriceProfile<R>(Vh,VF) ); ffassert( &Uh == & Vh);}
	 else if (typemat.sym )
	   {  A.A.master( new  MatriceMorse<R>(Vh,typemat.sym,VF) ); 
	     ffassert( &Uh == & Vh);}
	 else 
	   {
	     A.A.master( new  MatriceMorse<R>(Vh,Uh,VF) ); // lines corresponding to test functions 
	   }
       }
     *A.A=R(); // reset value of the matrix
     
     if ( AssembleVarForm<R,MatriceCreuse<R>,FESpace >( stack,Th,Uh,Vh,typemat.sym,A.A,0,b->largs) )
       AssembleBC<R,FESpace>( stack,Th,Uh,Vh,typemat.sym,A.A,0,0,b->largs,tgv);
   }
  else
   { // add FH 17 06 2005  int on different meshes. 
     map<pair<int,int>, R >   AAA;
     bool bc=AssembleVarForm<R,map<pair<int,int>, R >,FESpace  >( stack,Th,Uh,Vh,typemat.sym,&AAA,0,b->largs);
     if (typemat.profile)
        { ExecError(" Sorry, construction of Skyline matrix with different meshes is not implemented! ");}
      else 
        { A.A.master( new  MatriceMorse<R>(Vh.NbOfDF,Uh.NbOfDF,AAA,typemat.sym) ); }
      if (bc)
           AssembleBC<R>( stack,Th,Uh,Vh,typemat.sym,A.A,0,0,b->largs,tgv);
    
   }
  if( A_is_square && factorize ) {
    MatriceProfile<R> * pf = dynamic_cast<MatriceProfile<R> *>((MatriceCreuse<R> *) A.A);
    assert(pf);
    switch (typemat.t) {
    case TypeSolveMat::LU: pf->LU(Abs(eps));break;
    case TypeSolveMat::CROUT: pf->crout(Abs(eps));break;
    case TypeSolveMat::CHOLESKY: pf->cholesky(Abs(eps));break;
    default: ExecError("Sorry no factorize for this type for matrix"); 
    }
    
  }    
  if (A_is_square) 
    SetSolver(stack,*A.A,&typemat,VF,eps,NbSpace,itermax,precon,umfpackstrategy,tgv,tol_pivot,tol_pivot_sym);
  
  return SetAny<Matrice_Creuse<R>  *>(&A);
  
}



bool SetGMRES();
bool SetCG();
#ifdef HAVE_LIBUMFPACK
bool SetUMFPACK();
#endif
/*
template<class R>
AnyType ProdMat(Stack,Expression ,Expression);
template<class R> AnyType DiagMat(Stack,Expression ,Expression);
template<class R> AnyType CopyTrans(Stack stack,Expression emat,Expression eA);
template<class R> AnyType CopyMat(Stack stack,Expression emat,Expression eA);
template<class R> AnyType CombMat(Stack stack,Expression emat,Expression combMat);
template<class R> AnyType MatFull2Sparse(Stack stack,Expression emat,Expression eA);
*/
namespace FreeFempp {

template<class R>
class TypeVarForm { public:
    aType tFB;                       
  //    aType tFB3;                       
    aType tMat;                       
    aType tMat3;                       
    aType tFL;                       
  //aType tFL3;                       
    aType tTab;                       
    aType tMatX;  
    aType tMatTX;  
    aType tDotStar;
    aType tBC ;                    
  // aType tBC3 ;                    
TypeVarForm() :
  tFB( atype<const  FormBilinear *>() ),                      
  //tFB3( atype<const  FormBilinear<v_fes3> *>() ),                      
  tMat( atype<Matrice_Creuse<R>*>() ),                       
  //  tMat3( atype<Matrice_Creuse<R,v_fes3>*>() ),                       
  tFL( atype<const  FormLinear *>() ),                       
  //tFL3( atype<const  FormLinear<v_fes3> *>() ),                       
  tTab( atype<KN<R> *>() ),                       
  tMatX( atype<typename VirtualMatrice<R>::plusAx >() ),  
  tMatTX( atype<typename VirtualMatrice<R>::plusAtx >() ),  
  tDotStar(atype< DotStar_KN_<R> >() ),
  tBC( atype<const  BC_set  *>())                     
  {  }


 static TypeVarForm *Global;
}; 

}
#endif
