#ifdef __MWERKS__
#pragma optimization_level 0
//#pragma inline_depth(0)
#endif

#include  <cmath>
#include  <iostream>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include <cstdio>
#include "MatriceCreuse_tpl.hpp"

//#include "fem3.hpp"
#include "MeshPoint.hpp"
#include <complex>
#include "Operator.hpp" 

#include <set>
#include <vector>

#include "lex.hpp"
#include "lgfem.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include "CGNL.hpp"
namespace bamg { class Triangles; }
namespace Fem2D { void DrawIsoT(const R2 Pt[3],const R ff[3],const RN_ & Viso); }

#include "BamgFreeFem.hpp"

static bool TheWait=false;
bool  NoWait=false;

extern long verbosity;
void init_lgmesh() ;


 class E_F_StackF0F0opt2 :public  E_F0mps { public:
  typedef   AnyType (*func)(Stack,Expression ,Expression ) ; 
  func f;
  Expression a0,a1;
  E_F_StackF0F0opt2(func ff,Expression aa0,Expression aa1) 
    : f(ff),a0(aa0),a1(aa1) {
     deque<pair<Expression,int> > ll;
     MapOfE_F0 m;
     size_t top =  currentblock->OffSet(0), topbb=top; // FH. bofbof ??? 
     int ret =aa1->Optimize(ll, m, top);
     a1 =   new E_F0_Optimize(ll,m,ret); 
    currentblock->OffSet(top-topbb);
  }
  AnyType operator()(Stack s)  const 
    {return  (*f)(s, a0 , a1) ;}  

};



const E_Array * Array(const C_F0 & a) { 
  if (a.left() == atype<E_Array>() )
    return dynamic_cast<const E_Array *>(a.LeftValue());
  else 
    return 0;
}
bool Box2(const C_F0 & bb, Expression * box)
{
  const E_Array * a= Array(bb);
  if(a && a->size() == 2)
   {
    box[0] = to<double>((*a)[0]);
    box[1] = to<double>((*a)[1]);
    return true;
    }
  else 
    return false;
  
}
bool Box2x2(Expression  bb, Expression * box)
{
  const E_Array * a= dynamic_cast<const E_Array *>(bb);
  if(a && a->size() == 2)
     return Box2((*a)[0],box)  &&  Box2((*a)[1],box+2) ;
  else 
    return false;
}
void dump_table()
{
   cout << " dump the table of the language " << endl;
   cout << " ------------------------------ " <<endl <<endl;
	map<const string,basicForEachType *>::const_iterator i; ;
	       
	for (i= map_type.begin();i !=map_type.end();i++)
	  {
	   cout << " type : " << i->first << endl;
	   if( i->second )
	     i->second->ShowTable(cout);
	   else cout << " Null \n";
	   cout << "\n\n";     
	  }

	for (i= map_type.begin();i !=map_type.end();i++)
	  {
	   cout << " type : " << i->first << endl;
	   if( i->second )
	     i->second->ShowTable(cout);
	   else cout << " Null \n";
	   cout << "\n\n";    
	   
	  }
	  
	 cout << "--------------------- " << endl;
	 cout << *TheOperators << endl;
	 cout << "--------------------- " << endl;

}
/*
 class LocalArrayVariable:public E_F0
 { 
  size_t offset;
  aType t; //  type of the variable just for check  
  Expression  n; // expression of the size 
  public:
  AnyType operator()(Stack s) const { 
    SHOWVERB( cout << "\n\tget var " << offset << " " <<  t->name() << endl);  
   return PtrtoAny(static_cast<void *>(static_cast<char *>(s)+offset),t);}
  LocalArrayVariable(size_t o,aType tt,Expression nn):offset(o),t(tt),n(nn) {throwassert(tt);
    SHOWVERB(cout << "\n--------new var " << offset << " " <<  t->name() << endl);
    }
};
*/
template<class A,class B>  A Build(B b) {  return A(b);}

template<class K>
class OneBinaryOperatorA_inv : public OneOperator { public:  
  OneBinaryOperatorA_inv() : OneOperator(atype<Matrice_Creuse_inv<K> >(),atype<Matrice_Creuse<K> *>(),atype<long>()) {}
    E_F0 * code(const basicAC_F0 & args) const 
     { Expression p=args[1];
       if ( ! p->EvaluableWithOutStack() ) 
        { 
          bool bb=p->EvaluableWithOutStack();
          cout << bb << " " <<  * p <<  endl;
          CompileError(" A^p, The p must be a constant == -1, sorry");}
       long pv = GetAny<long>((*p)(0));
        if (pv !=-1)   
         { char buf[100];
           sprintf(buf," A^%d, The pow must be  == -1, sorry",pv);
           CompileError(buf);}     
       return  new E_F_F0<Matrice_Creuse_inv<R>,Matrice_Creuse<K> *>(Build<Matrice_Creuse_inv<R>,Matrice_Creuse<K> *>,t[0]->CastTo(args[0])); 
    }
};



basicAC_F0::name_and_type  MatrixInterpolation::Op::name_param[]= {
     "t=", &typeid(bool) 
};

bool In(long *viso,int n,long v) 
{
  int i=0,j=n,k;
  if  (v <viso[0] || v >viso[j-1]) 
    return false;  
  while (i<j-1)    
   if ( viso[k=(i+j)/2]> v) j=k;
   else i=k;
  return viso[i]=v;
}



class  LinkToInterpreter { public:
 Type_Expr   P,N,x,y,z,label,region,nu_triangle,nu_edge,lenEdge,hTriangle,area,inside;
  LinkToInterpreter() ; 
};

LinkToInterpreter * l2interpreter;


  using namespace Fem2D;

inline pmesh  ReadMesh( string * const & s) {
  Mesh * m=new Mesh(*s);
  R2 Pn,Px;
  m->BoundingBox(Pn,Px);
  m->quadtree=new FQuadTree(m,Pn,Px,m->nv);
  delete s;
  return m;
 }
 
 

template<class Result,class A>
 class E_F_A_Ptr_o_R :public  E_F0 { public:
    typedef Result A::* ptr;
  ptr p; 
  Expression a0;
  E_F_A_Ptr_o_R(Expression aa0,ptr pp) 
    : a0(aa0),p(pp) {}
  AnyType operator()(Stack s)  const {
    return SetAny<Result*>(&(GetAny<A*>((*a0)(s))->*p));}
  bool MeshIndependent() const {return a0->MeshIndependent();} // 
    
};
//  ----
//  remarque pas de template, cela ne marche pas encore ......
 class E_P_Stack_P   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<R3*>(&MeshPointStack(s)->P);} 
    operator aType () const { return atype<R3*>();}         
    
}; 
 class E_P_Stack_Px   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<R*>(&MeshPointStack(s)->P.x);} 
    operator aType () const { return atype<R*>();}         
    
}; 
 class E_P_Stack_Py   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const {throwassert(* (long *) s);
    return SetAny<R*>(&MeshPointStack(s)->P.y);} 
    operator aType () const { return atype<R*>();}         
    
}; 
 class E_P_Stack_Pz   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<R*>(&MeshPointStack(s)->P.z);} 
    operator aType () const { return atype<R*>();}         
    
}; 

 class E_P_Stack_N   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<R3*>(&MeshPointStack(s)->N);} 
        operator aType () const { return atype<R3*>();}         

}; 
 class E_P_Stack_Nx   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<R*>(&MeshPointStack(s)->N.x);} 
    operator aType () const { return atype<R*>();}         
    
}; 
 class E_P_Stack_Ny   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<R*>(&MeshPointStack(s)->N.y);} 
    operator aType () const { return atype<R*>();}             
}; 
 class E_P_Stack_Nz   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<R*>(&MeshPointStack(s)->N.z);} 
    operator aType () const { return atype<R*>();}         
    
}; 

 class E_P_Stack_Region   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<long*>(&MeshPointStack(s)->region);} 
    operator aType () const { return atype<long*>();}         
    
}; 
 class E_P_Stack_Label   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<long*>(&MeshPointStack(s)->label);} 
    operator aType () const { return atype<long *>();}         
    
}; 
 class E_P_Stack_Mesh   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<pmesh >(const_cast<pmesh>(MeshPointStack(s)->Th));} 
  operator aType () const { return atype<pmesh>();}         

}; 
 class E_P_Stack_Nu_Triangle   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<long>(MeshPointStack(s)->t);} 
    operator aType () const { return atype<long>();}         
    
}; 
 class E_P_Stack_Nu_Vertex   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<long>(MeshPointStack(s)->v);} 
    operator aType () const { return atype<long>();}         
    
}; 
 class E_P_Stack_Nu_Face   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<long>(MeshPointStack(s)->f);} 
    operator aType () const { return atype<long>();}         
    
}; 
 class E_P_Stack_Nu_Edge   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<long>(MeshPointStack(s)->e);} 
    operator aType () const { return atype<long>();}         
    
}; 
 class E_P_Stack_inside   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    return SetAny<double>(MeshPointStack(s)->outside? 0.0 : 1.0 );} 
    operator aType () const { return atype<double>();}         
    
}; 

class E_P_Stack_lenEdge   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    MeshPoint * mp=MeshPointStack(s);
    assert(mp->T && mp ->e >=0);
    double l= mp->T->lenEdge(mp->e);
    return SetAny<double>(l);} 
    operator aType () const { return atype<double>();}         
    
}; 

class E_P_Stack_hTriangle   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    MeshPoint * mp=MeshPointStack(s);
    assert(mp->T) ;
    double l= mp->T->h();
    return SetAny<double>(l);} 
    operator aType () const  { return atype<double>();}         
    
};

class E_P_Stack_nTonEdge   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    MeshPoint * mp=MeshPointStack(s);
    assert(mp->T && mp->e > -1  ) ;
    long l=mp->Th->nTonEdge(mp->t,mp->e);
    // cout << " nTonEdge " << l << endl;
    return SetAny<long>( l) ;} 
    operator aType () const  { return atype<long>();}         
    
}; 

class E_P_Stack_areaTriangle   :public  E_F0mps { public: 
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    MeshPoint * mp=MeshPointStack(s);
    assert(mp->T) ;
    double l= mp->T->area;
    return SetAny<double>(l);} 
    operator aType () const { return atype<double>();}         
    
}; 

/*
class New_MeshPoint : public E_F0mps {

};*/


 

 
//inline pfes MakePtr(pmesh * const &  a, TypeOfFE * const & tef){ return pfes(new pfes_tef(a,tef)) ;}
//inline pfes MakePtr(pmesh * const &  a){ return pfes(new pfes_tef(a,&P1Lagrange)) ;}
//inline pfes MakePtr(pfes * const &  a,long const & n){ return pfes(new pfes_fes(a,n)) ;}
/*
class E_pfes : public E_F0 {
  const int N;
  Expression *Th,*tef; 
  E_pfes(Expression *TTh,*ttef) : Th(TTh),tef(ttef),N(ttef?ttef->nbitem():0) {}
   virtual AnyType operator()(Stack)  const {
   return AnyType<pfes*>
  }
  
  virtual size_t nbitem() const {return 1;}  
};
  OneOperator_pfes():OneOperator(atype<void>(),atype<E_Array>(),atype<E_Array>()) {}
    E_F0 * code(const basicAC_F0 & args) const ; 
  
};*/

template<class R>
class LinearCG : public OneOperator 
{ public:
  typedef KN<R> Kn;
  typedef KN_<R> Kn_;
  const int cas;

 class MatF_O: VirtualMatrice<R> { public:
   Stack stack;
   C_F0 c_x;
   mutable  Kn x;
   Expression  mat;
   typedef  typename VirtualMatrice<R>::plusAx plusAx;
   MatF_O(int n,Stack stk,const OneOperator * op) 
     : stack(stk),
       x(n),c_x(CPValue(x)),
       mat(op->code(basicAC_F0_wa(c_x))) {}
   ~MatF_O() { delete mat;delete c_x.LeftValue();}
   void addMatMul(const  Kn_  & xx, Kn_ & Ax) const { 
      assert(xx.N()==Ax.N());
      x =xx;
      Ax  += *GetAny<Kn*>((*mat)(stack)); } 
    plusAx operator*(const Kn &  x) const {return plusAx(this,x);} 
};  
 

  class E_LCG: public E_F0mps { public:
   const int cas;// <0 => Nolinear
   static basicAC_F0::name_and_type name_param[] ;
   static const int n_name_param =4;
   Expression nargs[n_name_param];
  const OneOperator *A, *C; 
  Expression X,B;

  E_LCG(const basicAC_F0 & args,int cc) :cas(cc)
   {
      args.SetNameParam(n_name_param,name_param,nargs);
      {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(args[0].LeftValue());
         assert(op);
         A = op->Find("(",ArrayOfaType(atype<Kn* >(),false)); }
      if (nargs[2]) 
      {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[2]);
         assert(op); 
         C = op->Find("(",ArrayOfaType(atype<Kn* >(),false)); }
       else  C =0;
      X = to<Kn*>(args[1]);
      if (args.size()>2)
        B = to<Kn*>(args[2]);
      else 
        B=0;
   }
     
     virtual AnyType operator()(Stack stack)  const {
      Kn &x = *GetAny<Kn *>((*X)(stack));
      if (B) {
        Kn &b = *GetAny<Kn *>((*B)(stack));
        R p = (b,b);
        if (p) 
           ExecError("Sorry LinearCG work only with nul right hand side, so put the right hand in the function");
      }
      int n=x.N();
      MatF_O AA(n,stack,A);
      double eps = 1.0e-6;
      int nbitermax=  100;
      if (nargs[0]) eps= GetAny<double>((*nargs[0])(stack));
      if (nargs[1]) nbitermax = GetAny<long>((*nargs[1])(stack));
      if (nargs[3]) eps= *GetAny<double*>((*nargs[3])(stack));
      
      int ret;
      if (cas<0) {
       if (C) 
         { MatF_O CC(n,stack,C);
           ret = NLCG(AA,CC,x,nbitermax,eps, 51L-Min(Abs(verbosity),50L) );}
        else 
           ret = NLCG(AA,MatriceIdentite<R>(),x,nbitermax,eps, 51L-Min(Abs(verbosity),50L));
        }
      else 
      if (C) 
       { MatF_O CC(n,stack,C);
         ret = ConjuguedGradient2(AA,CC,x,nbitermax,eps, 51L-Min(Abs(verbosity),50L) );}
      else 
         ret = ConjuguedGradient2(AA,MatriceIdentite<R>(),x,nbitermax,eps, 51L-Min(Abs(verbosity),50L));
      if( nargs[3]) *GetAny<double*>((*nargs[3])(stack)) = -(eps);
      return SetAny<long>(ret);
       
     }  
    operator aType () const { return atype<long>();}         
     
  };
  
  E_F0 * code(const basicAC_F0 & args) const {
    return new E_LCG(args,cas);}
  LinearCG() :   OneOperator(atype<long>(),
                             atype<Polymorphic*>(),
                             atype<KN<R> *>(),atype<KN<R> *>()),cas(2){}
  LinearCG(int cc) :   OneOperator(atype<long>(),
                             atype<Polymorphic*>(),
                             atype<KN<R> *>()),cas(cc){}
   
};


basicAC_F0::name_and_type  LinearCG<R>::E_LCG::name_param[]= {
     "eps", &typeid(double)  ,
     "nbiter",&typeid(long) ,
     "precon",&typeid(Polymorphic*),
     "veps" ,  &typeid(double*) 
};


template<class R>
class LinearGMRES : public OneOperator 
{ public:
  typedef KN<R> Kn;
  typedef KN_<R> Kn_;
  const int cas;

 class MatF_O: VirtualMatrice<R> { public:
   Stack stack;
   C_F0 c_x;
   mutable  Kn x;
   Expression  mat;
   typedef  typename VirtualMatrice<R>::plusAx plusAx;
   MatF_O(int n,Stack stk,const OneOperator * op) 
     : stack(stk),
       x(n),c_x(CPValue(x)),
       mat(op->code(basicAC_F0_wa(c_x))) {}
   ~MatF_O() { delete mat;delete c_x.LeftValue();}
   void addMatMul(const  Kn_  & xx, Kn_ & Ax) const { 
      assert(xx.N()==Ax.N());
      x =xx;
      Ax  += *GetAny<Kn*>((*mat)(stack)); } 
    plusAx operator*(const Kn &  x) const {return plusAx(this,x);} 
};  
 

  class E_LGMRES: public E_F0mps { public:
   const int cas;// <0 => Nolinear
   static basicAC_F0::name_and_type name_param[] ;
   static const int n_name_param =5;
   Expression nargs[n_name_param];
  const OneOperator *A, *C; 
  Expression X,B;

  E_LGMRES(const basicAC_F0 & args,int cc) :cas(cc)
   {
      args.SetNameParam(n_name_param,name_param,nargs);
      {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(args[0].LeftValue());
         assert(op);
         A = op->Find("(",ArrayOfaType(atype<Kn* >(),false)); }
      if (nargs[2]) 
      {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[2]);
         assert(op); 
         C = op->Find("(",ArrayOfaType(atype<Kn* >(),false)); }
       else  C =0;
      X = to<Kn*>(args[1]);
      if (args.size()>2)
        B = to<Kn*>(args[2]);
      else 
        B=0;
   }
     
     virtual AnyType operator()(Stack stack)  const {
      Kn &x = *GetAny<Kn *>((*X)(stack));
      Kn b(x.n);
     
      if (B)   b = *GetAny<Kn *>((*B)(stack));
      else     b=0;
      int n=x.N();
      int dKrylov=50;
      MatF_O AA(n,stack,A);
      double eps = 1.0e-6;
      int nbitermax=  100;
      if (nargs[0]) eps= GetAny<double>((*nargs[0])(stack));
      if (nargs[1]) nbitermax = GetAny<long>((*nargs[1])(stack));
      if (nargs[3]) eps= *GetAny<double*>((*nargs[3])(stack));
      if (nargs[4]) dKrylov= GetAny<long>((*nargs[4])(stack));
      
      int ret;
      if(verbosity>4)
        cout << "  ..GMRES: eps= " << eps << " max iter " << nbitermax 
             << " dim of Krylov space " << dKrylov << endl;
        KNM<R> H(dKrylov+1,dKrylov+1);
       int k=dKrylov,nn=n;
       double epsr=eps;
      // int res=GMRES(a,(KN<R> &)x, (const KN<R> &)b,*this,H,k,nn,epsr);
      if (cas<0) {
        ErrorExec("NL GMRES:  to do! sorry ",1);
/*       if (C) 
         { MatF_O CC(n,stack,C);
           ret = NLGMRES(AA,CC,x,nbitermax,eps, 51L-Min(Abs(verbosity),50L) );}
        else 
           ret = NLGMRES(AA,MatriceIdentite<R>(),x,nbitermax,eps, 51L-Min(Abs(verbosity),50L));
         ConjuguedGradient  */
        }
      else 
       {
       if (C)
        { MatF_O CC(n,stack,C); 
         ret=GMRES(AA,(KN<R> &)x, (const KN<R> &)b,CC,H,k,nbitermax,epsr);}
       else
         ret=GMRES(AA,(KN<R> &)x, (const KN<R> &)b,MatriceIdentite<R>(),H,k,nbitermax,epsr);       
       }
       /*
      if (C) 
       { MatF_O CC(n,stack,C);
         
         ret = ConjuguedGradient2(AA,CC,x,nbitermax,eps, 51L-Min(Abs(verbosity),50L) );}
      else 
         ret = ConjuguedGradient2(AA,MatriceIdentite<R>(),x,nbitermax,eps, 51L-Min(Abs(verbosity),50L));*/
         
     // if( nargs[3]) *GetAny<double*>((*nargs[3])(stack)) = -(eps);
      return SetAny<long>(ret);
       
     }  
    operator aType () const { return atype<long>();}         
     
  };
  
  E_F0 * code(const basicAC_F0 & args) const {
    return new E_LGMRES(args,cas);}
  LinearGMRES() :   OneOperator(atype<long>(),
                             atype<Polymorphic*>(),
                             atype<KN<R> *>(),atype<KN<R> *>()),cas(2){}
  LinearGMRES(int cc) :   OneOperator(atype<long>(),
                             atype<Polymorphic*>(),
                             atype<KN<R> *>()),cas(cc){}
   
};


basicAC_F0::name_and_type  LinearGMRES<R>::E_LGMRES::name_param[]= {
     "eps", &typeid(double)  ,
     "nbiter",&typeid(long) ,
     "precon",&typeid(Polymorphic*),
     "veps" ,  &typeid(double*) ,
     "dimKrylov", &typeid(long) 
};

template<typename int2>
typename map<int,int2>::iterator closeto(map<int,int2> & m, int k)
 {
  typename map<int,int2>::iterator i=  m.find(k);
   if (i==m.end()) 
    {
     i=  m.find(k+1);
     if (i==m.end()) 
      i=  m.find(k-1);
    }
   return i;
 }
 
template<class T,int N>
class Smallvect { public: 
 T v[N];
 T & operator[](int i){return v[i];}
};


template<class T> 
int numeroteclink(KN_<T> & ndfv) 
{
         int nbdfv =0;
         for (int i=0;i<ndfv.N();i++)
           if (ndfv[i]>=i)
           {
             int j=i,ii,kkk=0;
             
             do {
               ii=ndfv[j];
               assert(kkk++<10);
             //  assert(ii>=nbdfv);
               ndfv[j]=nbdfv ;
               j=ii;
               }             
             while (j!=nbdfv);
             if (verbosity > 100) 
             cout << "    ndf: " <<  j << " " <<  ii  << " <- " << nbdfv << " " <<  kkk <<  endl;
              nbdfv++;
           }
       return nbdfv;
}


bool BuildPeriodic( 
  int nbcperiodic,
  Expression *periodic,
  const Mesh &Th,Stack stack,
  int & nbdfv, KN<int> & ndfv,int & nbdfe, KN<int> & ndfe) { 
  
/*
  build numbering of vertex form 0 to nbdfv-1
  and build numbering  of  edge form 0 to nbdfe-1
  we removing common vextex or common edge    
  --  we suppose one df by vertex 
      nbdfv number of df on vertex 
      ndfv[i]  given the numero of the df of the vertex 
  -- we suppose 1 df
*/  
   typedef Smallvect<int,2> int2;    
    if (nbcperiodic ) {
      
   //    KN<int> ndfv(Th.nv);
    //   KN<int> ndfe(Th.neb);
       assert(ndfv.N()==Th.nv);
       assert(ndfe.N()==Th.neb);
        
       MeshPoint *mp=MeshPointStack(stack),smp=*mp;   
       int n= nbcperiodic;
       if (verbosity >2)
         cout << " Nb of pair of periodic conditions: = " << n <<  endl;
       int * link1=0;
       int * link2=0;
       KN<int*> plk1(n),plk2(n);
       KN<int> nlk1(n),nlk2(n);
       KN<int> lab1(n),lab2(n);
#ifndef  HUGE_VAL      
       const double infty= numeric_limits<double>::infinity();
#else
       const double infty= HUGE_VAL;
#endif       
       int nblink1, nblink2;
       int *plink1 , *plink2;
        for (int step=0;step<2;step++)
         {
           nblink1=0,     nblink2=0;
           plink1=link1,  plink2=link2;
           for (int ip=0, k=0;ip<n;ip++,k+=4)
            {
               int label1=GetAny<long>((*periodic[k+0])(stack));
               int label2=GetAny<long>((*periodic[k+2])(stack));
               lab1[ip]=label1;
               lab2[ip]=label2;
               
               int l1=nblink1;
               int l2=nblink2;
               plk1[ip]= plink1;
               plk2[ip]= plink2;
               
               for (int ke=0;ke<Th.neb;ke++)
                {
                 if (Th.bedges[ke].lab==label1)
                  {
                    if (plink1) *plink1++=ke;
                    nblink1++;
                  }
                 else if (Th.bedges[ke].lab==label2)
                  {
                    if (plink2) *plink2++=ke;
                    nblink2++;
                   }
                 }
               nlk1[ip]= nblink1-l1;
               nlk2[ip]= nblink2-l2;              
            }
            if(step) break; // no reallocl 
            if (verbosity >3)
            cout << "  Periodic = " << nblink1 << " " << nblink2 << " step=" << step << endl;
            link1 = new int[nblink1];
            link2 = new int[nblink2];
            if(nblink1 != nblink2)
            {
             ExecError("Periodic:  the both number of edges is not the same ");
            }
         }
        if ( nblink1 >0) 
        {
        for (int ip=0, k=0;ip<n;ip++,k+=4)
          {
            map<int,int2> m;
            const int kk1=1,kk2=3;
            int label1=lab1[ip],label2=lab2[ip];
            int n1=nlk1[ip],n2=nlk2[ip];
            int *pke1=plk1[ip], *pke2=plk2[ip];
            double xmn=infty,xmx=-infty,hmn=infty;
            if (verbosity >1)
            cout << "  --Update: periodic  couple label1= " << label1 << ", n edges= " << n1 << "; "
                                          << ", label2= " << label2<<  ", n edges= " << n2 <<endl; 
            if (n1 != n2) ExecError("periodic BC:  the number of edges is not the same");
            for (int i1=0;i1<n1;i1++)
             {
              const BoundaryEdge & e =Th.bedges[pke1[i1]];
              if (e.lab==label1) 
                {           
                 mp->set(e[0].x,e[0].y);
                 double x0=GetAny<double>((*periodic[k+kk1])(stack));
                 mp->set(e[1].x,e[1].y);
                 double x1=GetAny<double>((*periodic[k+kk1])(stack));
                 xmn=Min(x1,x0,xmn);
                 xmx=Max(x1,x0,xmx); 
                 hmn=Min(hmn,Abs(x1-x0));
               }                                
             }
            assert(hmn>1.0e-20);
            double coef = 8/hmn;
            double x0 = xmn;
            if (verbosity > 2)
            cout << "  --Update: periodic " << xmn << " " << xmx << " " << " h=" << hmn << endl;
            assert(coef>1e-10 && (xmx-xmn)*coef < 1.e7 );
             
           //  map construction ----
           for (int i1=0;i1<n1;i1++)
             {
              int ie=pke1[i1];
              const BoundaryEdge & e =Th.bedges[pke1[i1]];
              if (e.lab==label1)
                 for (int ne=0;ne<2;ne++)
                  {
                   int2 i2;
                   i2[0]=ie;
                   i2[1]=-1;
                   mp->set(e[ne].x,e[ne].y);
                   double xx=GetAny<double>((*periodic[k+kk1])(stack));
                   int i0= (int) ((xx-x0)*coef);
                   map<int,int2>::iterator im=closeto(m,i0);
                   if (im==m.end())
                    {
                     if (verbosity >50)
                     cout << xx << " " << i0 << " " << ie << endl;
                     im=m.insert(make_pair<int,int2>(i0,i2)).first;
                    }
                   else {
                     if (verbosity >50)
                     cout << xx << " " << i0 << " " << ie << " :  " << im->second[0] << " " << im->second[1] << endl;
                    assert( (im->second[1] < 0) && (im->second[0] >=0) );
                    im->second[1]=ie;}
                   
               }                                
             }
            
            for (int i2=0;i2<n2;i2++)
             {
              int ie2=pke2[i2];
              const BoundaryEdge & e =Th.bedges[ie2];
              if (e.lab==label2)
                {
                 if (verbosity >50)
                 cout << Th(e[0]) << " " << Th(e[1]) << ":: ";
                 mp->set(e[0].x,e[0].y);
                 double xx0=GetAny<double>((*periodic[k+kk2])(stack));
                 mp->set(e[1].x,e[1].y);
                 double xx1=GetAny<double>((*periodic[k+kk2])(stack));
                 int i0= int((xx0-x0)*coef);
                 int i1= int((xx1-x0)*coef);
                 map<int,int2>::iterator im0=closeto(m,i0);
                 map<int,int2>::iterator im1=closeto(m,i1);
                 if(im0 == m.end() || im1 == m.end() )
                 	{  ExecError("periodic: Sorry one vertex of edge is losted "); }
                  int ie1=-1;
                 if      ((ie1=im0->second[0])==im1->second[1]) ;
                 else if ((ie1=im0->second[0])==im1->second[1]) ;
                 else if ((ie1=im0->second[1])==im1->second[1]) ;
                 else if ((ie1=im0->second[1])==im1->second[0]) ;
                 else if ((ie1=im0->second[0])==im1->second[0]) ;
                 else
                  {
                   cout << ie2 << " ~ " << im0->second[0] << " " << im0->second[1] << ", " 
                                        << im1->second[0] << " " << im1->second[1] << endl;
                   ExecError("periodic: Sorry one egde is losted "); }
                 const BoundaryEdge & ep =Th.bedges[ie1];
                 mp->set(ep[0].x,ep[0].y);
                 double yy0=GetAny<double>((*periodic[k+kk2])(stack));
                 mp->set(ep[1].x,ep[1].y);
                 double yy1=GetAny<double>((*periodic[k+kk2])(stack));
                 
                 pke1[i2]=ie1*2+ ( ( (yy1-yy0) < 0)  == ( (xx1-xx0) < 0) ) ;
                 if (verbosity >50)
                 cout << "  edge " << ie1 << " <=> " << ie2 << " " 
                      << ( ( (yy1-yy0) < 0)  == ( (xx1-xx0) < 0) ) << "; "
                      << xx0 << " " <<xx1<<" <=> "  << yy0 << " " <<yy1<< 
                      "  ::  " <<  Th(ep[0]) << " " << Th(ep[1]) << endl ;

                }
              }
              
             
          }
          
        *mp = smp;
        for (int i=0;i<Th.neb;i++)
            ndfe[i]=i;// circular link
        for (int i=0;i<Th.nv;i++)
          ndfv[i]=i;// circular link
        for (int i=0;i<nblink1;i++)
         {
           int ie1=link1[i]/2;
           int sens = link1[i]%2;
           int ie2=link2[i];
           assert(ie1!=ie2);
           {
             int k=ie2,kk;
            do  
                  k=ndfe[kk=k];
            while  ( (k!=ie1) && (k!=ie2) );
            
            if (k != ie1) {  // merge of two list 
                  int iee1= ndfe[ie1];
                  ndfe[k] =iee1;
                  ndfe[ie1]=kk;  
              }                  
           }
           for (int ke2=0;ke2<2;ke2++)
             {
                int ke1=ke2;
                if(!sens) ke1=1-ke1; 
               int iv1=Th(Th.bedges[ie1][ke1]);
               int iv2=Th(Th.bedges[ie2][ke2]);
               if (verbosity >50)
               cout << "  vertex:" << iv1 << ": " << Th(iv1) << ", <=> " << iv2 << ": " <<  Th(iv2) << endl;
               if (iv1 > iv2) Exchange(iv1,iv2);
               int j=0;
               int k=iv2,kk;
                do  {
                  k=ndfv[kk=k];
                  if(verbosity>60) cout << k << " " << kk << " " << iv1 << " " << iv2 <<endl;
                  assert(j++<nblink1*4);
                  assert( k != kk || (k==iv2) || (k==iv1)  );
                   }
                 while  ( (k!=iv2) && (k!=iv1) );

               if (k != iv1) {  // merge of two list 
                  int ivv1= ndfv[iv1];
                  ndfv[k] =ivv1;
                  ndfv[iv1]=kk;  
                  if (verbosity >50)
                  cout << "  vertex " << k << " <=> " << iv1 <<  " " << ivv1 << " " << kk << endl;
                 }                  
             }
            
         } 
         // generation de numero de dlt
         
          nbdfv = numeroteclink(ndfv) ; 
          nbdfe = numeroteclink(ndfe) ; 
          if (verbosity>2) 
            cout << " -- nb df on vertices " << nbdfv << endl;
        delete [] link1;
        delete [] link2;
        return true; //new FESpace(**ppTh,*tef,nbdfv,ndfv,nbdfe,ndfe);
      }
   }
   return false;   
}


bool  v_fes::buildperiodic(Stack stack,int & nbdfv, KN<int> & ndfv,int & nbdfe, KN<int> & ndfe) { 
  return BuildPeriodic(nbcperiodic,periodic,**ppTh,stack,nbdfv,ndfv,nbdfe,ndfe);
/*
   typedef Smallvect<int,2> int2;    
    if (nbcperiodic ) {
       const Mesh &Th(**ppTh);
   //    KN<int> ndfv(Th.nv);
    //   KN<int> ndfe(Th.neb);
       assert(ndfv.N()==Th.nv);
       assert(ndfe.N()==Th.neb);
        
       MeshPoint *mp=MeshPointStack(stack),smp=*mp;   
       int n= nbcperiodic;
       if (verbosity >2)
         cout << " Nb of pair of periodic conditions: = " << n <<  endl;
       int * link1=0;
       int * link2=0;
       KN<int*> plk1(n),plk2(n);
       KN<int> nlk1(n),nlk2(n);
       KN<int> lab1(n),lab2(n);
#ifndef  HUGE_VAL      
       const double infty= numeric_limits<double>::infinity();
#else
       const double infty= HUGE_VAL;
#endif       
       int nblink1, nblink2;
       int *plink1 , *plink2;
        for (int step=0;step<2;step++)
         {
           nblink1=0,     nblink2=0;
           plink1=link1,  plink2=link2;
           for (int ip=0, k=0;ip<n;ip++,k+=4)
            {
               int label1=GetAny<long>((*periodic[k+0])(stack));
               int label2=GetAny<long>((*periodic[k+2])(stack));
               lab1[ip]=label1;
               lab2[ip]=label2;
               
               int l1=nblink1;
               int l2=nblink2;
               plk1[ip]= plink1;
               plk2[ip]= plink2;
               
               for (int ke=0;ke<Th.neb;ke++)
                {
                 if (Th.bedges[ke].lab==label1)
                  {
                    if (plink1) *plink1++=ke;
                    nblink1++;
                  }
                 else if (Th.bedges[ke].lab==label2)
                  {
                    if (plink2) *plink2++=ke;
                    nblink2++;
                   }
                 }
               nlk1[ip]= nblink1-l1;
               nlk2[ip]= nblink2-l2;              
            }
            if(step) break; // no reallocl 
            if (verbosity >3)
            cout << "  Periodic = " << nblink1 << " " << nblink2 << " step=" << step << endl;
            link1 = new int[nblink1];
            link2 = new int[nblink2];
            if(nblink1 != nblink2)
            {
             ExecError("Periodic:  the both number of edges is not the same ");
            }
         }
        if ( nblink1 >0) 
        {
        for (int ip=0, k=0;ip<n;ip++,k+=4)
          {
            map<int,int2> m;
            const int kk1=1,kk2=3;
            int label1=lab1[ip],label2=lab2[ip];
            int n1=nlk1[ip],n2=nlk2[ip];
            int *pke1=plk1[ip], *pke2=plk2[ip];
            double xmn=infty,xmx=-infty,hmn=infty;
            if (verbosity >1)
            cout << "  --Update: periodic  couple label1= " << label1 << ", n edges= " << n1 << "; "
                                          << ", label2= " << label2<<  ", n edges= " << n2 <<endl; 
            if (n1 != n2) ExecError("periodic BC:  the number of edges is not the same");
            for (int i1=0;i1<n1;i1++)
             {
              const BoundaryEdge & e =Th.bedges[pke1[i1]];
              if (e.lab==label1) 
                {           
                 mp->set(e[0].x,e[0].y);
                 double x0=GetAny<double>((*periodic[k+kk1])(stack));
                 mp->set(e[1].x,e[1].y);
                 double x1=GetAny<double>((*periodic[k+kk1])(stack));
                 xmn=Min(x1,x0,xmn);
                 xmx=Max(x1,x0,xmx); 
                 hmn=Min(hmn,Abs(x1-x0));
               }                                
             }
            assert(hmn>1.0e-20);
            double coef = 8/hmn;
            double x0 = xmn;
            if (verbosity > 2)
            cout << "  --Update: periodic " << xmn << " " << xmx << " " << " h=" << hmn << endl;
            assert(coef>1e-10 && (xmx-xmn)*coef < 1.e7 );
             
           //  map construction ----
           for (int i1=0;i1<n1;i1++)
             {
              int ie=pke1[i1];
              const BoundaryEdge & e =Th.bedges[pke1[i1]];
              if (e.lab==label1)
                 for (int ne=0;ne<2;ne++)
                  {
                   int2 i2;
                   i2[0]=ie;
                   i2[1]=-1;
                   mp->set(e[ne].x,e[ne].y);
                   double xx=GetAny<double>((*periodic[k+kk1])(stack));
                   int i0= (int) ((xx-x0)*coef);
                   map<int,int2>::iterator im=closeto(m,i0);
                   if (im==m.end())
                    {
                     if (verbosity >50)
                     cout << xx << " " << i0 << " " << ie << endl;
                     im=m.insert(make_pair<int,int2>(i0,i2)).first;
                    }
                   else {
                     if (verbosity >50)
                     cout << xx << " " << i0 << " " << ie << " :  " << im->second[0] << " " << im->second[1] << endl;
                    assert( (im->second[1] < 0) && (im->second[0] >=0) );
                    im->second[1]=ie;}
                   
               }                                
             }
            
            for (int i2=0;i2<n2;i2++)
             {
              int ie2=pke2[i2];
              const BoundaryEdge & e =Th.bedges[ie2];
              if (e.lab==label2)
                {
                 if (verbosity >50)
                 cout << Th(e[0]) << " " << Th(e[1]) << ":: ";
                 mp->set(e[0].x,e[0].y);
                 double xx0=GetAny<double>((*periodic[k+kk2])(stack));
                 mp->set(e[1].x,e[1].y);
                 double xx1=GetAny<double>((*periodic[k+kk2])(stack));
                 int i0= int((xx0-x0)*coef);
                 int i1= int((xx1-x0)*coef);
                 map<int,int2>::iterator im0=closeto(m,i0);
                 map<int,int2>::iterator im1=closeto(m,i1);
                 if(im0 == m.end() || im1 == m.end() )
                 	{  ExecError("periodic: Sorry one vertex of edge is losted "); }
                  int ie1=-1;
                 if      ((ie1=im0->second[0])==im1->second[1]) ;
                 else if ((ie1=im0->second[0])==im1->second[1]) ;
                 else if ((ie1=im0->second[1])==im1->second[1]) ;
                 else if ((ie1=im0->second[1])==im1->second[0]) ;
                 else if ((ie1=im0->second[0])==im1->second[0]) ;
                 else
                  {
                   cout << ie2 << " ~ " << im0->second[0] << " " << im0->second[1] << ", " 
                                        << im1->second[0] << " " << im1->second[1] << endl;
                   ExecError("periodic: Sorry one egde is losted "); }
                 const BoundaryEdge & ep =Th.bedges[ie1];
                 mp->set(ep[0].x,ep[0].y);
                 double yy0=GetAny<double>((*periodic[k+kk2])(stack));
                 mp->set(ep[1].x,ep[1].y);
                 double yy1=GetAny<double>((*periodic[k+kk2])(stack));
                 
                 pke1[i2]=ie1*2+ ( ( (yy1-yy0) < 0)  == ( (xx1-xx0) < 0) ) ;
                 if (verbosity >50)
                 cout << "  edge " << ie1 << " <=> " << ie2 << " " 
                      << ( ( (yy1-yy0) < 0)  == ( (xx1-xx0) < 0) ) << "; "
                      << xx0 << " " <<xx1<<" <=> "  << yy0 << " " <<yy1<< 
                      "  ::  " <<  Th(ep[0]) << " " << Th(ep[1]) << endl ;

                }
              }
              
             
          }
          
        *mp = smp;
        for (int i=0;i<Th.neb;i++)
            ndfe[i]=i;// circular link
        for (int i=0;i<Th.nv;i++)
          ndfv[i]=i;// circular link
        for (int i=0;i<nblink1;i++)
         {
           int ie1=link1[i]/2;
           int sens = link1[i]%2;
           int ie2=link2[i];
           assert(ie1!=ie2);
           {
             int k=ie2,kk;
            do  
                  k=ndfe[kk=k];
            while  ( (k!=ie1) && (k!=ie2) );
            
            if (k != ie1) {  // merge of two list 
                  int iee1= ndfe[ie1];
                  ndfe[k] =iee1;
                  ndfe[ie1]=kk;  
              }                  
           }
           for (int ke2=0;ke2<2;ke2++)
             {
                int ke1=ke2;
                if(!sens) ke1=1-ke1; 
               int iv1=Th(Th.bedges[ie1][ke1]);
               int iv2=Th(Th.bedges[ie2][ke2]);
               if (verbosity >50)
               cout << "  vertex:" << iv1 << ": " << Th(iv1) << ", <=> " << iv2 << ": " <<  Th(iv2) << endl;
               if (iv1 > iv2) Exchange(iv1,iv2);
               int j=0;
               int k=iv2,kk;
                do  {
                  k=ndfv[kk=k];
                  if(verbosity>60) cout << k << " " << kk << " " << iv1 << " " << iv2 <<endl;
                  assert(j++<nblink1*4);
                  assert( k != kk || (k==iv2) || (k==iv1)  );
                   }
                 while  ( (k!=iv2) && (k!=iv1) );

               if (k != iv1) {  // merge of two list 
                  int ivv1= ndfv[iv1];
                  ndfv[k] =ivv1;
                  ndfv[iv1]=kk;  
                  if (verbosity >50)
                  cout << "  vertex " << k << " <=> " << iv1 <<  " " << ivv1 << " " << kk << endl;
                 }                  
             }
            
         } 
         // generation de numero de dlt
         
          nbdfv = numeroteclink(ndfv) ; 
          nbdfe = numeroteclink(ndfe) ; 
          if (verbosity>2) 
            cout << " -- nb df on vertices " << nbdfv << endl;
        delete [] link1;
        delete [] link2;
        return true; //new FESpace(**ppTh,*tef,nbdfv,ndfv,nbdfe,ndfe);
      }
   }
   return false;   
   */
}
#ifdef ZZZZZZZZ
FESpace * pfes_tef::update() { 
   typedef Smallvect<int,2> int2; 
    
    if (nbcperiodic ) {
       const Mesh &Th(**ppTh);
       KN<int> ndfv(Th.nv);
       KN<int> ndfe(Th.neb);
       int nbdfv,nbdfe;
/*       
       MeshPoint *mp=MeshPointStack(stack),smp=*mp;   
       int n= nbcperiodic;
       if (verbosity >2)
         cout << " Nb of pair of periodic conditions: = " << n <<  endl;
       int * link1=0;
       int * link2=0;
       KN<int*> plk1(n),plk2(n);
       KN<int> nlk1(n),nlk2(n);
       KN<int> lab1(n),lab2(n);
#ifndef  HUGE_VAL      
       const double infty= numeric_limits<double>::infinity();
#else
       const double infty= HUGE_VAL;
#endif       
       int nblink1, nblink2;
       int *plink1 , *plink2;
        for (int step=0;step<2;step++)
         {
           nblink1=0,     nblink2=0;
           plink1=link1,  plink2=link2;
           for (int ip=0, k=0;ip<n;ip++,k+=4)
            {
               int label1=GetAny<long>((*periodic[k+0])(stack));
               int label2=GetAny<long>((*periodic[k+2])(stack));
               lab1[ip]=label1;
               lab2[ip]=label2;
               
               int l1=nblink1;
               int l2=nblink2;
               plk1[ip]= plink1;
               plk2[ip]= plink2;
               
               for (int ke=0;ke<Th.neb;ke++)
                {
                 if (Th.bedges[ke].lab==label1)
                  {
                    if (plink1) *plink1++=ke;
                    nblink1++;
                  }
                 else if (Th.bedges[ke].lab==label2)
                  {
                    if (plink2) *plink2++=ke;
                    nblink2++;
                   }
                 }
               nlk1[ip]= nblink1-l1;
               nlk2[ip]= nblink2-l2;              
            }
            if(step) break; // no reallocl 
            if (verbosity >3)
            cout << "  Periodic = " << nblink1 << " " << nblink2 << " step=" << step << endl;
            link1 = new int[nblink1];
            link2 = new int[nblink2];
            if(nblink1 != nblink2)
            {
             ExecError("Periodic:  the both number of edges is not the same ");
            }
         }
        if ( nblink1 >0) 
        {
        for (int ip=0, k=0;ip<n;ip++,k+=4)
          {
            map<int,int2> m;
            const int kk1=1,kk2=3;
            int label1=lab1[ip],label2=lab2[ip];
            int n1=nlk1[ip],n2=nlk2[ip];
            int *pke1=plk1[ip], *pke2=plk2[ip];
            double xmn=infty,xmx=-infty,hmn=infty;
            if (verbosity >1)
            cout << "  --Update: periodic  couple label1= " << label1 << ", n edges= " << n1 << "; "
                                          << ", label2= " << label2<<  ", n edges= " << n2 <<endl; 
            if (n1 != n2) ExecError("periodic BC:  the number of edges is not the same");
            for (int i1=0;i1<n1;i1++)
             {
              const BoundaryEdge & e =Th.bedges[pke1[i1]];
              if (e.lab==label1) 
                {           
                 mp->set(e[0].x,e[0].y);
                 double x0=GetAny<double>((*periodic[k+kk1])(stack));
                 mp->set(e[1].x,e[1].y);
                 double x1=GetAny<double>((*periodic[k+kk1])(stack));
                 xmn=Min(x1,x0,xmn);
                 xmx=Max(x1,x0,xmx); 
                 hmn=Min(hmn,Abs(x1-x0));
               }                                
             }
            assert(hmn>1.0e-20);
            double coef = 8/hmn;
            double x0 = xmn;
            if (verbosity > 2)
            cout << "  --Update: periodic " << xmn << " " << xmx << " " << " h=" << hmn << endl;
            assert(coef>1e-10 && (xmx-xmn)*coef < 1.e7 );
             
           //  map construction ----
           for (int i1=0;i1<n1;i1++)
             {
              int ie=pke1[i1];
              const BoundaryEdge & e =Th.bedges[pke1[i1]];
              if (e.lab==label1)
                 for (int ne=0;ne<2;ne++)
                  {
                   int2 i2;
                   i2[0]=ie;
                   i2[1]=-1;
                   mp->set(e[ne].x,e[ne].y);
                   double xx=GetAny<double>((*periodic[k+kk1])(stack));
                   int i0= (int) ((xx-x0)*coef);
                   map<int,int2>::iterator im=closeto(m,i0);
                   if (im==m.end())
                    {
                     if (verbosity >50)
                     cout << xx << " " << i0 << " " << ie << endl;
                     im=m.insert(make_pair<int,int2>(i0,i2)).first;
                    }
                   else {
                     if (verbosity >50)
                     cout << xx << " " << i0 << " " << ie << " :  " << im->second[0] << " " << im->second[1] << endl;
                    assert( (im->second[1] < 0) && (im->second[0] >=0) );
                    im->second[1]=ie;}
                   
               }                                
             }
            
            for (int i2=0;i2<n2;i2++)
             {
              int ie2=pke2[i2];
              const BoundaryEdge & e =Th.bedges[ie2];
              if (e.lab==label2)
                {
                 if (verbosity >50)
                 cout << Th(e[0]) << " " << Th(e[1]) << ":: ";
                 mp->set(e[0].x,e[0].y);
                 double xx0=GetAny<double>((*periodic[k+kk2])(stack));
                 mp->set(e[1].x,e[1].y);
                 double xx1=GetAny<double>((*periodic[k+kk2])(stack));
                 int i0= int((xx0-x0)*coef);
                 int i1= int((xx1-x0)*coef);
                 map<int,int2>::iterator im0=closeto(m,i0);
                 map<int,int2>::iterator im1=closeto(m,i1);
                 if(im0 == m.end() || im1 == m.end() )
                 	{  ExecError("periodic: Sorry one vertex of edge is losted "); }
                  int ie1=-1;
                 if      ((ie1=im0->second[0])==im1->second[1]) ;
                 else if ((ie1=im0->second[0])==im1->second[1]) ;
                 else if ((ie1=im0->second[1])==im1->second[1]) ;
                 else if ((ie1=im0->second[1])==im1->second[0]) ;
                 else if ((ie1=im0->second[0])==im1->second[0]) ;
                 else
                  {
                   cout << ie2 << " ~ " << im0->second[0] << " " << im0->second[1] << ", " 
                                        << im1->second[0] << " " << im1->second[1] << endl;
                   ExecError("periodic: Sorry one egde is losted "); }
                 const BoundaryEdge & ep =Th.bedges[ie1];
                 mp->set(ep[0].x,ep[0].y);
                 double yy0=GetAny<double>((*periodic[k+kk2])(stack));
                 mp->set(ep[1].x,ep[1].y);
                 double yy1=GetAny<double>((*periodic[k+kk2])(stack));
                 
                 pke1[i2]=ie1*2+ ( ( (yy1-yy0) < 0)  == ( (xx1-xx0) < 0) ) ;
                 if (verbosity >50)
                 cout << "  edge " << ie1 << " <=> " << ie2 << " " 
                      << ( ( (yy1-yy0) < 0)  == ( (xx1-xx0) < 0) ) << "; "
                      << xx0 << " " <<xx1<<" <=> "  << yy0 << " " <<yy1<< 
                      "  ::  " <<  Th(ep[0]) << " " << Th(ep[1]) << endl ;

                }
              }
              
             
          }
          
        *mp = smp;
        for (int i=0;i<Th.neb;i++)
            ndfe[i]=i;// circular link
        for (int i=0;i<Th.nv;i++)
          ndfv[i]=i;// circular link
        for (int i=0;i<nblink1;i++)
         {
           int ie1=link1[i]/2;
           int sens = link1[i]%2;
           int ie2=link2[i];
           assert(ie1!=ie2);
           {
             int k=ie2,kk;
            do  
                  k=ndfe[kk=k];
            while  ( (k!=ie1) && (k!=ie2) );
            
            if (k != ie1) {  // merge of two list 
                  int iee1= ndfe[ie1];
                  ndfe[k] =iee1;
                  ndfe[ie1]=kk;  
              }                  
           }
           for (int ke2=0;ke2<2;ke2++)
             {
                int ke1=ke2;
                if(!sens) ke1=1-ke1; 
               int iv1=Th(Th.bedges[ie1][ke1]);
               int iv2=Th(Th.bedges[ie2][ke2]);
               if (verbosity >50)
               cout << "  vertex:" << iv1 << ": " << Th(iv1) << ", <=> " << iv2 << ": " <<  Th(iv2) << endl;
               if (iv1 > iv2) Exchange(iv1,iv2);
               int j=0;
               int k=iv2,kk;
                do  {
                  k=ndfv[kk=k];
                  if(verbosity>60) cout << k << " " << kk << " " << iv1 << " " << iv2 <<endl;
                  assert(j++<nblink1*4);
                  assert( k != kk || (k==iv2) || (k==iv1)  );
                   }
                 while  ( (k!=iv2) && (k!=iv1) );

               if (k != iv1) {  // merge of two list 
                  int ivv1= ndfv[iv1];
                  ndfv[k] =ivv1;
                  ndfv[iv1]=kk;  
                  if (verbosity >50)
                  cout << "  vertex " << k << " <=> " << iv1 <<  " " << ivv1 << " " << kk << endl;
                 }                  
             }
            
         } 
         // generation de numero de dlt
         
          int nbdfv = numeroteclink(ndfv) ; 
          int nbdfe = numeroteclink(ndfe) ; 
          if (verbosity>2) 
            cout << " -- nb df on vertices " << nbdfv << endl;
        delete [] link1;
        delete [] link2;
*/        
        return  new FESpace(**ppTh,*tef,nbdfv,ndfv,nbdfe,ndfe);
      }
/*     else 
      {
        if (verbosity) 
          cout << " -- Warning no edge for periodic condition." << endl;
       return new FESpace(**ppTh,*tef);
       }
      }*/
     else 
       return  new FESpace(**ppTh,*tef);
}
#endif

struct OpMake_pfes: public OneOperator { 
    static const int n_name_param =1;
    static basicAC_F0::name_and_type name_param[] ;

  class Op: public E_F0mps { public:
  
     Expression eppTh;
     Expression eppfes;
     const E_Array & atef;
     int nb;
     int nbcperiodic;
     Expression *periodic;
     
     Op(Expression ppfes,Expression ppTh, const E_Array & aatef,int nbp,Expression * pr) 
       : eppTh(ppTh),eppfes(ppfes),atef(aatef),nbcperiodic(nbp),periodic(pr) {       
       }
      ~Op() { if(periodic) delete []periodic;}
     AnyType operator()(Stack s)  const {       
       const pmesh * ppTh = GetAny<pmesh *>( (*eppTh)(s) );
       AnyType r = (*eppfes)(s) ;
       const TypeOfFE ** tef= new  const TypeOfFE * [ atef.size()];
       for (int i=0;i<atef.size();i++)
         tef[i]= GetAny<TypeOfFE *>(atef[i].eval(s));
       pfes * ppfes = GetAny<pfes *>(r);
       bool same = true;
       for (int i=1;i<atef.size();i++)
         same &= atef[i].LeftValue() == atef[1].LeftValue();
      /* if ( same)
       *ppfes = new pfes_tefk(ppTh,tef,atef.size(),stack,nbcperiodic,periodic);
       else */
       *ppfes = new pfes_tefk(ppTh,tef,atef.size(),s    ,nbcperiodic,periodic);
       (**ppfes).decrement();

     //  delete [] tef;
       return r;}
  } ;
    E_F0 * code(const basicAC_F0 & args)  const
     { 
    int nbcperiodic=0;
    Expression *periodic=0;
    Expression nargs[n_name_param];
    
      args.SetNameParam(n_name_param,name_param,nargs);
      GetPeriodic(nargs[0],nbcperiodic,periodic);

      aType t_tfe= atype<TypeOfFE*>();
       const E_Array * a2(dynamic_cast<const E_Array *>(args[2].LeftValue()));
       assert(a2);
       int N = a2->size(); ;
       if (!N) CompileError(" We wait an array of Type of Element ");
       for (int i=0;i< N; i++) 
          if ((*a2)[i].left() != t_tfe) CompileError(" We wait an array of  Type of Element ");
        
       return  new OpMake_pfes::Op(args[0],args[1],*a2,nbcperiodic,periodic);} 
   OpMake_pfes() : 
     OneOperator(atype<pfes*>(),atype<pfes*>(),atype<pmesh *>(),atype<E_Array>()) {}
};

inline pfes* MakePtr2(pfes * const &p,pmesh * const &  a, TypeOfFE * const & tef)
   { *p=new pfes_tef(a,tef) ;
      (**p).decrement();
      return p;}
      

class OP_MakePtr2 { public:
  class Op : public E_F0mps  { public:
 //  static int GetPeriodic(Expression  bb, Expression & b,Expression & f);
    static const int n_name_param =1;
    static basicAC_F0::name_and_type name_param[] ;
    Expression nargs[n_name_param];
    typedef pfes * R;
    typedef pfes * A;
    typedef pmesh * B;
    typedef TypeOfFE * C;
    Expression a,b,c;
    int nbcperiodic ;
    Expression *periodic;
    Op(const basicAC_F0 & args);
    
    AnyType operator()(Stack s) const  {
      A p= GetAny<A>( (*a)(s) );
      B th= GetAny<B>( (*b)(s) );
      C tef= GetAny<C>( (*c)(s) );   
     //  cout << " ----------- " << endl;  
      *p=new pfes_tef(th,tef,s,nbcperiodic,periodic) ;
      (**p).decrement();
      return  SetAny<R>(p);
    } 
  }; // end Op class 
  
  typedef Op::R Result;
   static  E_F0 * f(const basicAC_F0 & args) { return  new Op(args);}
   static ArrayOfaType  typeargs() {
     return ArrayOfaType(       
       atype<Op::A>(),
       atype<Op::B>(),
       atype<Op::C>(),false ) ;}
};  

void GetPeriodic(Expression perio,    int & nbcperiodic ,    Expression * &periodic)
{
      if ( perio) 
       {
         if( verbosity>1) 
         cout << "  -- Periodical Condition to do" << endl;
         const E_Array * a= dynamic_cast<const  E_Array *>(perio);
         assert(a);
         int n = a->size();
        nbcperiodic= n/2;
        if( verbosity>1) 
        cout << "    the number of periodicBC " << n << endl;
        if ( 2*nbcperiodic != n ) CompileError(" Sorry the number of periodicBC must by even"); 
        periodic = new Expression[n*2]; 
        for (int i=0,j=0;i<n;i++,j+=2)
          if (GetPeriodic((*a)[i],periodic[j],periodic[j+1])==0)
            CompileError(" a sub array of periodic BC must be [label, realfunction ]");
        }


}

    
OP_MakePtr2::Op::Op(const basicAC_F0 & args)
  : a(to<A>(args[0])),b(to<B>(args[1])),c(to<C>(args[2]))
     {
      nbcperiodic=0;
      periodic=0;
      args.SetNameParam(n_name_param,name_param,nargs);
      GetPeriodic(nargs[0],nbcperiodic,periodic);
     }

int GetPeriodic(Expression  bb, Expression & b,Expression & f)
    {
      const E_Array * a= dynamic_cast<const E_Array *>(bb);
      if(a && a->size() == 2)
       {
         b= to<long>((*a)[0]);
         f= to<double>((*a)[1]);
         return 1;
       }
      else 
       return 0;
   }

basicAC_F0::name_and_type  OP_MakePtr2::Op::name_param[]= {
     "periodic", &typeid(E_Array) 
};
basicAC_F0::name_and_type  OpMake_pfes::name_param[]= {
     "periodic", &typeid(E_Array) 
};

inline pfes* MakePtr2(pfes * const &p,pmesh * const &  a){ 
      *p=new pfes_tef(a,&P1Lagrange);
      (**p).decrement();
       return p ;}
       
inline pfes* MakePtr2(pfes * const &p,pfes * const &  a,long const & n){
       *p= new pfes_fes(a,n);
      (**p).decrement();
       return p ;}
       
 long FindTxy(Stack s,pmesh * const &  ppTh,const double & x,const double & y)
{
   R2 P(x,y),PHat;
   bool outside;
     MeshPoint & mp = *MeshPointStack(s);  
   const Mesh * pTh= *ppTh;
   if(pTh == 0) return 0;  
   const Triangle * K=pTh->Find(mp.P,PHat,outside);
   if (!outside)
     mp.set(*pTh,P,PHat,*K,K->lab);
   else return 0;
   return 1;
}

 


 

template<class K>    
KN<K> * pfer2vect( pair<FEbase<R> *,int> p)
 {  
    KN<K> * x=p.first->x();
    if ( !x) {  // defined 
      FESpace * Vh= p.first->newVh();     
      throwassert( Vh);
      *p.first = x = new KN<K>(Vh->NbOfDF);
      *x=0; 
    }
    return x;}

template<class K>        
long pfer_nbdf(pair<FEbase<K> *,int> p)
 {  
   if (!p.first->Vh) p.first->Vh= p.first->newVh();
   throwassert( !!p.first->Vh);
   return p.first->Vh->NbOfDF;
 }
 
double pmesh_area(pmesh * p)
 { throwassert(p && *p) ;  return (**p).area ;}
long pmesh_nt(pmesh * p)
 { throwassert(p && *p) ;  return (**p).nt ;}
long pmesh_nv(pmesh * p)
 { throwassert(p && *p) ;  return (**p).nv ;}


class Op_CopyArray : public OneOperator { public:
    Op_CopyArray():OneOperator(atype<void>(),atype<E_Array>(),atype<E_Array>()) {}
    E_F0 * code(const basicAC_F0 & args) const ; 
};  
                       
template<class R,int dd>
AnyType pfer2R(Stack s,const AnyType &a)
{
  pfer ppfe=GetAny<pfer>(a);
   FEbase<R> & fe( *ppfe.first);
  int componante=ppfe.second;
  if ( !fe.x()) {
   if ( !fe.x()){
    // CompileError(" Sorry unset fem array ");
     return   SetAny<R>(0.0);
    }
  }
   
  const FESpace & Vh(*fe.Vh);
  const Mesh & Th(Vh.Th);
  MeshPoint & mp = *MeshPointStack(s);
  const Triangle *K;
  R2 PHat;
  bool outside=false;
  bool qnu=true;
  if ( mp.Th == &Th && mp.T) 
   {
    qnu=false;
    K=mp.T;
    PHat=mp.PHat;
   }
  else if ( mp.other.Th == & Th && mp.other.P.x == mp.P.x && mp.other.P.y == mp.P.y )
   {
    K=mp.other.T;
    PHat=mp.other.PHat;
    outside = mp.other.outside;
   } 
  else {
    if (mp.isUnset()) ExecError("Try to get unset x,y, ...");
    K=Th.Find(mp.P,PHat,outside);
    mp.other.set(Th,(R2) mp.P,PHat,*K,0,outside);
    }
  const FElement KK(Vh[Th(K)]);
  if (outside && !KK.tfe->NbDfOnVertex && !KK.tfe->NbDfOnEdge) 
    return   SetAny<R>(0.0); 
/*  if (!outside) 
    {
      if ( Norme2_2( (*K)(PHat) - mp.P ) > 1e-12 )
        cout << "bug ??  " << Norme2_2( (*K)(PHat) - mp.P ) << " " << mp.P << " " << (*K)(PHat) << endl;
    } */
/*  int nbdf=KK.NbDoF();
  
  int N= KK.N;
  KN_<R> U(*fe.x());
  KNMK<R> fb(nbdf,N,3); //  the value for basic fonction
  KN<R> fk(nbdf);
  for (int i=0;i<nbdf;i++) // get the local value
    fk[i] = U[KK(i)];
    //  get value of basic function
  KK.BF(PHat,fb);  
  
  R r = (fb('.',componante,dd),fk);  
*/
  
 const R rr = KK(PHat,*fe.x(),componante,dd);
//  cout << " " << rr << endl;
//  R2 B(mp.P);
/*   if ( r < 0.08001 &&  Norme2_2(mp.P) > 0.05  && Norme2_2(mp.P) < 0.4*0.4  ) 
     {
     int vv=verbosity;
      cout << " f()  triangle  " << Th(K) << " " << mp.P << " " << PHat << " =  " << r << " " <<outside ;
      if (outside) {  verbosity = 200;
         K=Th.Find(mp.P,PHat,outside);
          cout << Th(K) << " " << outside << endl;}
       cout << endl; verbosity=vv;
     } */
//  if ( qnu )   
 //  cout << " f()  triangle       " << Th(K) << " " << mp.P << " " << PHat << " =  " << r <<  endl;
  return SetAny<R>(rr);
}

AnyType set_fe (Stack s,Expression ppfe, Expression e)
  { 
   long kkff = Mesh::kfind,  kkth = Mesh::kthrough;

   
    MeshPoint *mps=MeshPointStack(s),mp=*mps;  
    pair<FEbase<R> *,int>  pp=GetAny<pair<FEbase<R> *,int> >((*ppfe)(s));
    FEbase<R> & fe(*pp.first);
    const  FESpace & Vh(*fe.newVh());
    KN<R> gg(Vh.MaximalNbOfDF()); 
    const  Mesh & Th(Vh.Th);
 //   R F[100]; // buffer 
    TabFuncArg tabexp(s,Vh.N);
    tabexp[0]=e;
    throwassert(Vh.N==1); 
    KN<R> * y=new  KN<R>(Vh.NbOfDF);
    KN<R> & yy(*y);
    KN<R> Viso(100);
    R2 Ptt[3];
    for (int i=0;i<Viso.N();i++)
      Viso[i]=0.01*i; 
      
    FElement::aIPJ ipj(Vh[0].Pi_h_ipj()); 
    FElement::aR2  PtHat(Vh[0].Pi_h_R2()); 
    KN<R>   Aipj(ipj.N());
    KN<R>   Vp(PtHat.N());
    const E_F0 & ff(* (const  E_F0 *) e ) ;

   if (Vh.isFEMesh() )
    {
      
      assert(Vh.NbOfDF == Th.nv && Vh.N == 1 );
      for (int iv=0;iv<Th.nv;iv++)
       {
         const Vertex & v(Th(iv));
         int ik=Th.Contening(&v);
         const Triangle & K(Th[ik]);
         int il=-1;
         if  ( &K[0] == &v) il=0;
         if  ( &K[1] == &v) il=1;
         if  ( &K[2] == &v) il=2;
         assert(il>=0);
         mps->set(Th,v,TriangleHat[il],K,v.lab);
         yy[iv] = GetAny<double>( ff(s) );
       }
      
    }
   else     
    
     for (int t=0;t<Th.nt;t++)
      {
         FElement K(Vh[t]);
         int nbdf=K.NbDoF();
         gg=0;   
#ifdef OLDPih    
// old method          
         K.Pi_h(gg,F_Pi_h,F,&tabexp);
#else               
         K.Pi_h(Aipj);
         for (int p=0;p<PtHat.N();p++)
          { 
            mps->set(K.T(PtHat[p]),PtHat[p],K);
            Vp[p]=GetAny<double>( ff(s) );
           }
         for (int i=0;i<Aipj.N();i++)
          { 
           const FElement::IPJ &ipj_i(ipj[i]);
           assert(ipj_i.j==0); // car  Vh.N=0
           gg[ipj_i.i] += Aipj[i]*Vp[ipj_i.p];            
          }
#endif          
         for (int df=0;df<nbdf;df++)         
            (*y)[K(df)] =  gg[df] ;
                       
      }
    *mps=mp;
    fe=y;
     kkff = Mesh::kfind - kkff;
     kkth = Mesh::kthrough -kkth;

    if(verbosity>1)
      cout << " -- vector function's bound   " << y->min() << " " << y->max() 
           << " " << kkth << "/" << kkff << " =  " << double(kkth)/Max<double>(1.,kkff) << endl;
    return SetAny<FEbase<R>*>(&fe); 
  }
AnyType set_feoX_1 (Stack s,Expression ppfeX_1, Expression e)
  { // inutile 
    // même chose que  v(X1,X2);
    typedef const interpolate_f_X_1<R>::CODE * code;
    MeshPoint mp= *MeshPointStack(s); 
    code ipp = dynamic_cast<code>(ppfeX_1);
    
    pair<FEbase<R> *,int>  pp=GetAny<pair<FEbase<R> *,int> >((*ipp->f)(s));
    FEbase<R> & fe(*pp.first);
    const  FESpace & Vh(*fe.newVh());
    KN<R> gg(Vh.MaximalNbOfDF()); 
    const  Mesh & Th(Vh.Th);
    R F[100]; // buffer 
    TabFuncArg tabexp(s,Vh.N+2);
    tabexp[0]=e;
    tabexp[1]=ipp->x;
    tabexp[2]=ipp->y;
    
    throwassert(Vh.N==1); 
    KN<R> * y=new  KN<R>(Vh.NbOfDF); 
     for (int t=0;t<Th.nt;t++)
      {
         FElement K(Vh[t]);
         int nbdf=K.NbDoF();
        
         gg=0;
         
         K.Pi_h(gg,FoX_1_Pi_h,F,&tabexp);
         for (int df=0;df<nbdf;df++)
          (*y)[K(df)] =  gg[df] ;
      }
    *MeshPointStack(s)=mp;
    fe=y;
    if(verbosity>1)
    cout << " -- interpole f= g*X^-1, function's bound:  " << y->min() << " " << y->max() << endl; 
    return SetAny<FEbase<R>*>(&fe); 
  }
class E_set_fev: public E_F0mps {public:
  E_Array  aa;
  Expression   ppfe;
  bool optimize;
       vector<size_t>  where_in_stack_opt;
       Expression optiexp0,optiexpK;

  E_set_fev(const E_Array * a,Expression pp) ;
  
  AnyType operator()(Stack)  const ;
  operator aType () const { return atype<void>();} 
             
};


E_set_fev::E_set_fev(const E_Array * a,Expression pp) 
  :aa(*a),ppfe(pp),optimize(true),
   where_in_stack_opt(0),optiexp0(0),optiexpK(0) 
   { 
    aa.map(to<R>) ;
    bool kdump=false;
    if(optimize)
			{ // new code Optimez  -------
			     int n=aa.size();
			    deque<pair<Expression,int> > ll;
			    MapOfE_F0 m;
			    where_in_stack_opt.resize(n);
			    size_t top = currentblock->OffSet(0), topbb=top; // FH. bofbof ??? 
			    for (int i=0; i<n; i++)
			      {
			       Expression ee= aa[i].LeftValue();
			       if (kdump)
			       cout << "Optimize OneOperatorMakePtrFE:  type exp: " << typeid(*ee).name() << " "<<endl;
			       where_in_stack_opt[i]=ee->Optimize(ll, m, top);
			       if (kdump)
			       cout  << "\n\t\t"<< i  << ": " << where_in_stack_opt[i] << endl;
			      }
			      
			    currentblock->OffSet(top-topbb);
			    //  
			    int k=ll.size(),k0=0,k1=0;
			    for (int i=0;i<k;i++)
			        if (ll[i].first->MeshIndependent()) k0++;
			    deque<pair<Expression,int> > l0(k0),l1(k-k0);
			    k0=0,k1=0;
			    for (int i=0;i<k;i++)
			       if (ll[i].first->MeshIndependent()) 
			         {
			          if (kdump)
			          cout << " mi " << ll[i].second << " " << *(ll[i].first) << endl;
			          l0[k0++]=ll[i];
			         }
			        else 
			         {
			          if (kdump)
			          cout << " md " << ll[i].second << " " << *(ll[i].first) << endl;
			          l1[k1++]=ll[i];
			         }
			    if (k0)      
			      optiexp0 = new E_F0_Optimize(l0,m,0);  // constant part
			    if (k1) 
			      optiexpK = new E_F0_Optimize(l1,m,0);  // none constant part
			            
			 }
   
   
   }
AnyType E_set_fev::operator()(Stack s)  const
{   
    MeshPoint *mps=MeshPointStack(s), mp=*mps;   
    FEbase<R> ** pp=GetAny< FEbase<R> **>((*ppfe)(s));
    FEbase<R> & fe(**pp);
    const  FESpace & Vh(*fe.newVh());
    KN<R> gg(Vh.MaximalNbOfDF()); 
    
    
    const  Mesh & Th(Vh.Th);
    const int dim=Vh.N;
    double ** copt=0;
    if (optimize)   copt= new double *[dim];
    if(copt) {
      assert(dim== where_in_stack_opt.size());
      for (int i=0;i<dim;i++)
       {
        int offset=where_in_stack_opt[i];
        assert(offset>10);
        copt[i]= static_cast<double *>(static_cast<void *>((char*)s+offset));
        *(copt[i])=0;
       }
    if (optiexp0) (*optiexp0)(s); // init 
    }
    assert(dim<100);
 //   R F[100]; // buffer 
    TabFuncArg tabexp(s,Vh.N);
//   const E_Array * aa = dynamic_cast<const E_Array *>(e);
   throwassert( aa.size() == Vh.N);
   for (int i=0;i<dim;i++)
     tabexp[i]=aa[i]; 
   KN<R> * y=new  KN<R>(Vh.NbOfDF);
   KN<R> & yy(*y);
    FElement::aIPJ ipj(Vh[0].Pi_h_ipj()); 
    FElement::aR2  PtHat(Vh[0].Pi_h_R2()); 
    KN<R>   Aipj(ipj.N());
    KNM<R>  Vp(dim,PtHat.N());
   
   if (Vh.isFEMesh() )
    {
      
      assert(Vh.NbOfDF == Th.nv && dim == 1 );
      for (int iv=0;iv<Th.nv;iv++)
       {
         const E_F0 & ff(* (const  E_F0 *) aa[0]  ) ;
         const Vertex & v(Th(iv));
         int ik=Th.Contening(&v);
         const Triangle & K(Th[ik]);
         int il=-1;
         if  ( &K[0] == &v) il=0;
         if  ( &K[1] == &v) il=1;
         if  ( &K[2] == &v) il=2;
         assert(il>=0);
         mps->set(Th,v,TriangleHat[il],K,v.lab);
         if (copt) {
           if (optiexpK) (*optiexpK)(s); 
            yy[iv] =  *(copt[0]);
          }
          else 
           yy[iv] = GetAny<double>( ff(s) );
       }
      
    }
   else
     for (int t=0;t<Th.nt;t++)
      {
         FElement K(Vh[t]);
         int nbdf=K.NbDoF();
        
         gg=0;
         
#ifdef OLDPih    
// old method          
         K.Pi_h(gg,F_Pi_h,F,&tabexp);
#else               
         K.Pi_h(Aipj);
         for (int p=0;p<PtHat.N();p++)
          { 
            mps->set(K.T(PtHat[p]),PtHat[p],K);

            KN_<double> Vpp(Vp('.',p));
            if (copt) { // optimize  version 
             if (optiexpK) (*optiexpK)(s);
             for (int j=0;j<dim;j++)
               Vpp[j] = *(copt[j]);
            }
            else  // old version 
            for (int j=0;j<dim;j++)
             if (tabexp[j]) 
               Vpp[j]=GetAny<double>( (*tabexp[j])(s) );
              else Vpp[j]=0;
           }
           
         for (int i=0;i<Aipj.N();i++)
          { 
           const FElement::IPJ &ipj_i(ipj[i]);
           gg[ipj_i.i] += Aipj[i]*Vp(ipj_i.j,ipj_i.p);           
          }
#endif          

         for (int df=0;df<nbdf;df++)         
           yy[K(df)] =  gg[df] ;
      }
   //  MeshPointStack(s)->unset();
    fe=y;
    if (copt) delete [] copt;
    *MeshPointStack(s) = mp;
    if(verbosity>1)
    cout << " -- function's bound  (DF)" << y->min() << " " << y->max() << endl; 
    return Nothing;
  }
  



inline pferbase MakePtrFE(pfes * const &  a){ 
  pferbase p=new FEbase<R>(a);
  //cout << "MakePtrFE " << p<< endl; 
  return p ;}
  
inline pferbase* MakePtrFE2(pferbase * const &  p,pfes * const &  a){ 
  *p=new FEbase<R>(a);
  //cout << "MakePtrFE2 " << *p<< endl; 
  return p ;}
  
inline pferbasearray* MakePtrFE3(pferbasearray * const &  p,pfes * const &  a,const long & N){ 
  *p=new FEbaseArray<R>(a,N);
  //cout << "MakePtrFE2 " << *p<< endl; 
  return p ;}
/*
inline pmesharray*  MakePtr(pmesharray*  const &  p,long   const &  a){ 
  p->first=new pmesh [a];
  p->second=a;
  for (int i=0;i<a;i++) 
     p->first[i]=0; // nuset 
  return p ;}
*/  
class  OneOperatorMakePtrFE : public OneOperator { public:
 // il faut Optimize 
    typedef pferbase* R;
    typedef pfes* B;
    typedef double K;
    class CODE : public E_F0mps  { public:
       Expression fer,fes;
       E_set_fev * e_set_fev;
       const E_Array * v;
       CODE(const basicAC_F0 & args) : 
         fer(to<R>(args[0])),
         fes(to<B>(args[1])),
         e_set_fev(0) 
          {
          	if (BCastToR(args[2]) )
              v = new E_Array(basicAC_F0_wa(to<K>(args[2]))); 
            else 
              v = dynamic_cast<const E_Array *>( args[2].LeftValue() );
            if (!v) ErrorCompile(" We wait  a double expression or a array expression",1);
            //v->map(to<K>);
			  e_set_fev=  new   E_set_fev(v,fer);
 
          }
         
        AnyType operator()(Stack stack)  const {
          R  p = GetAny<R>( (*fer)(stack));
          B  a = GetAny<B>( (*fes)(stack));          
          *p=new FEbase<K>(a);
          (*e_set_fev)(stack); 
          return SetAny<R>(p);
        }
          operator aType () const { return atype<R>();}         
        
        
    };
   
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new CODE(args);}
    OneOperatorMakePtrFE(aType tt):  // tt= aType<double>() or aType<E_Array>()  
      OneOperator(map_type[typeid(R).name()],map_type[typeid(R).name()],map_type[typeid(B).name()],tt)
      {}
};


// ---  

template<class Result,class A>
class  OneOperator_Ptr_o_R: public OneOperator {
  //  aType r; //  return type 
    typedef Result A::* ptr;
    ptr p;   
    public: 
    E_F0 * code(const basicAC_F0 & args) const 
     { return  new E_F_A_Ptr_o_R<Result,A>(t[0]->CastTo(args[0]),p);} 
    OneOperator_Ptr_o_R(ptr pp):
       OneOperator(atype<Result*>(),atype<A*>()),p(pp) {}     
};

template<class K>  K  *PAddition(const K * a,const K *  b)  {return  new K(*a+*b);}





// 
class fCLD { public:
  typedef pair<int,Label> Key;
  typedef  map<Key,Expression>::iterator iterator;
  map<Key,Expression> *l;
  fCLD (){ l=new  map<Key,Expression>;}
  void operator=(const fCLD & a){ *l=*a.l;}
  void destroy() { delete l;l=0;}
  ~fCLD(){ delete l;l=0;}
  void Add(finconnue *v,int lab,C_F0 f) { 
    const MGauche *pn= v->simple();
     Check(pn,"Def CL Dirichet ");
    
    Label r(lab);
    Key k(make_pair(pn->first,r));
    iterator i=l->find(k);
    Check( i != l->end() ,"Def CL Dirichet already exist");
    l->insert(make_pair(k,CastTo<double>(f)));        
  }
};

//  pour stocker des expression de compilation 


class Convect : public E_F0mps  { public:
    typedef double  Result; // return type 
    Expression u,v,ff,dt;    
    Convect(const basicAC_F0 & args)  
    {
      args.SetNameParam(); 
      const E_Array * a = dynamic_cast<const E_Array *>(args[0].LeftValue());
            throwassert(a);
       if (a->size() != 2) 
          { CompileError("convect d'un vecteur à qui n'a pas 2 composantes");}
       u= CastTo<double>((*a)[0]);
       v= CastTo<double>((*a)[1]);
       
       dt=CastTo<double>(args[1]);
       ff=CastTo<double>(args[2]);
     }
    
    static ArrayOfaType  typeargs() 
      { return  ArrayOfaType(atype<E_Array>(),atype<double>(),atype<double>());}
      
    static  E_F0 * f(const basicAC_F0 & args) { return new Convect(args);} 
    AnyType operator()(Stack s) const ; 
    operator aType () const { return atype<Result>();}         
    
};    
   



class Plot :  public E_F0mps { public:
    typedef KN_<R>  tab;
    typedef pferbase sol;
    
    typedef long  Result;
    struct Expression2 {
     long what; // 0 mesh, 1 iso, 2 vector, 3 curve , 4 border 
     bool composant;
     Expression e[2];
     Expression2() {e[0]=0;e[1]=0;composant=false;what=0;}
     Expression &operator[](int i){return e[i];}
     sol eval(int i,Stack s,int & cmp) const  {  cmp=-1;
       if (e[i]) {
        if (!composant) {pfer p= GetAny< pfer >((*e[i])(s)); cmp=p.second;return p.first;}
        else {return GetAny< pferbase >((*e[i])(s));}
        }
       else return 0;}
     const Mesh & evalm(int i,Stack s) const  { throwassert(e[i]);return  * GetAny< pmesh >((*e[i])(s)) ;}
     const E_BorderN * evalb(int i,Stack s) const  { throwassert(e[i]);return   GetAny< const E_BorderN *>((*e[i])(s)) ;}
     tab  evalt(int i,Stack s) const  { throwassert(e[i]);return  GetAny<tab>((*e[i])(s)) ;}
    };

   static basicAC_F0::name_and_type name_param[] ;
   static const int n_name_param =15;
   Expression bb[4];
    vector<Expression2> l;
    Expression nargs[n_name_param];
    Plot(const basicAC_F0 & args) : l(args.size()) 
    {

      args.SetNameParam(n_name_param,name_param,nargs);
       if ( nargs[8] )
         Box2x2( nargs[8] , bb);   
         
      for (int i=0;i<l.size();i++)
       
         if (args[i].left()==atype<E_Array>())
          {
            l[i].composant=false;
            const E_Array * a = dynamic_cast<const E_Array *>(args[i].LeftValue());
            throwassert(a);
            if (a->size() >2) { CompileError("plot d'un vecteur à + 2 composantes");}
            if (BCastTo<pfer>((*a)[0]))
             {
              l[i].what=2;
              for (int j=0;j<a->size();j++)             
               l[i][j]= CastTo<pfer>((*a)[j]);
             }
            else 
             {
              l[i].what=3;
              for (int j=0;j<a->size();j++)             
               l[i][j]= CastTo<tab>((*a)[j]);
             }
          }
         else if (BCastTo<pferbase>(args[i])) {
          l[i].composant=true;
          l[i][0]=CastTo<pferbase>(args[i]); }
         else if (BCastTo<pfer>(args[i])) {
          l[i].composant=false;
          l[i].what=1;
          l[i][0]=CastTo<pfer>(args[i]);}
         else if (BCastTo<pmesh>(args[i])){
          l[i].what=0;
          l[i][0]=CastTo<pmesh>(args[i]);}
         else if (BCastTo<const E_BorderN *>(args[i])){
          l[i].what=4;
          l[i][0]=CastTo<const E_BorderN *>(args[i]);}
          
         else {
           CompileError("Sorry no way to plot this kind of data");
         }
    }
    
    static ArrayOfaType  typeargs() { return  ArrayOfaType(true);}// all type
    static  E_F0 * f(const basicAC_F0 & args) { return new Plot(args);} 
    AnyType operator()(Stack s) const ;
}; 


 basicAC_F0::name_and_type Plot::name_param[Plot::n_name_param] = {
     "coef", &typeid(double),
     "cmm", &typeid(string*),
     "ps", &typeid(string*)  ,
     "wait", &typeid(bool) ,
     "fill", &typeid(bool) ,
     "value", &typeid(bool) ,
     "clean", &typeid(bool) ,     
     "aspectratio", &typeid(bool),  
     "bb",&typeid(E_Array) ,
     "nbiso", &typeid(long), 
     "nbarrow", &typeid(long), 
     "viso", &typeid(KN_<double>),       
     "varrow", &typeid(KN_<double>),
     "bw",&typeid(bool),
     "grey", &typeid(bool)
     
   };




//template<class RR,class A>  void PrintP(RR* a, A  b){  *a <<*b;}
LinkToInterpreter::LinkToInterpreter()
{
   //P,N,x,y,z,label,region,nu_triangle;
   P=make_Type_Expr(atype<R3*>(),new E_P_Stack_P);
   x=make_Type_Expr(atype<R*>(),new E_P_Stack_Px);
   y=make_Type_Expr(atype<R*>(),new E_P_Stack_Py);
   z=make_Type_Expr(atype<R*>(),new E_P_Stack_Pz);
   N=make_Type_Expr(atype<R3*>(),new E_P_Stack_N);
   
   region=make_Type_Expr(new E_P_Stack_Region,atype<long*>());
   label=make_Type_Expr(new E_P_Stack_Label,atype<long*>());
   nu_triangle= make_Type_Expr(atype<long>(),new E_P_Stack_Nu_Triangle);
   nu_edge= make_Type_Expr(atype<long>(),new E_P_Stack_Nu_Edge);
   lenEdge    = make_Type_Expr(atype<R>(),new E_P_Stack_lenEdge);
   hTriangle  = make_Type_Expr(atype<R>(),new E_P_Stack_hTriangle);
   area       = make_Type_Expr(atype<R>(),new E_P_Stack_areaTriangle);
   inside     = make_Type_Expr(atype<R>(),new E_P_Stack_inside);
  Global.New("x",x);
  Global.New("y",y);
  Global.New("z",z);
  Global.New("label",label);
  Global.New("region",region);
  Global.New("nuTriangle",nu_triangle);   
  Global.New("nuEdge",nu_edge);   
  Global.New("P",P);   
  Global.New("N",N);   
  
  Global.New("lenEdge",lenEdge);   
  Global.New("area",area);   
  Global.New("hTriangle",hTriangle);
  Global.New("inside",inside);   
  Global.New("nTonEdge",make_Type_Expr(atype<long>(),new E_P_Stack_nTonEdge));   
  
}




template<class K>
struct set_eqmatrice_creuse_fbl: public binary_function<Matrice_Creuse<K>*,const Matrice_Creuse<K>  *,const  C_args * > {
  static Matrice_Creuse<K>* f(Matrice_Creuse<K>* const  & a,const  C_args * const & b)  {
  // 1  verif the FESpace 
  
  // 2 set = or += 
    
    throwassert(0);
    return a;}
};

template<class K>
struct set_eqvect_fl: public binary_function<KN<K>*,const  FormLinear *,KN<K>*> {
  static KN<K>* f(KN<K>* const  & a,const  FormLinear * const & b)  {
    throwassert(0);
    return a;}
};
  
  AnyType IntFunction::operator()(Stack stack) const  { 
  MeshPoint mp=* MeshPointStack(stack);
 double r=0;
 
             SHOWVERB(cout << " int " << endl);
             const vector<Expression>  & what(di->what);
             const Mesh  & Th = * GetAny<pmesh>( (*di->Th)(stack) );
             throwassert(&Th);
             CDomainOfIntegration::typeofkind kind = di->kind;
             set<int> setoflab;
             bool all=true; 
             if ( verbosity>3) 
               if (kind==CDomainOfIntegration::int1d) cout << "  -- boundary int border " ;
               else if (kind==CDomainOfIntegration::intalledges) cout << "  -- boundary int all edges " ;
               else cout << "  -- boundary int  " ;
             for (int i=0;i<what.size();i++)
               {long  lab  = GetAny<long>( (*what[i])(stack));
                setoflab.insert(lab);
                if ( verbosity>3) cout << lab << " ";
                all=false;
               }
               
             
             if (verbosity >3) 
               if (all) cout << " all " << endl ;
               else cout << endl;
               
             if (kind==CDomainOfIntegration::int1d)
               {
                 const QuadratureFormular1d & FI = QF_GaussLegendre2;
               
                 for( int e=0;e<Th.neb;e++)
                  {
                     if (all || setoflab.find(Th.bedges[e].lab) != setoflab.end())   
                       {     
                                 
                        int ie,i =Th.BoundaryTriangle(e,ie);
                        const Triangle & K(Th[i]);   
                        R2 E=K.Edge(ie);
                        R le = sqrt((E,E)); 
                        R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
                        PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
                        
                        for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                          {
                            QuadratureFormular1d::Point pi( FI[npi]);
                            R sa=pi.x,sb=1-sa;
                            R2 Pt(PA*sa+PB*sb ); //  
                            MeshPointStack(stack)->set(Th,K(Pt),Pt,K,Th.bedges[e].lab,R2(E.y,-E.x)/le,ie);
                            r += le*pi.a*GetAny<double>( (*fonc)(stack));
 		          }
                       }
                    }
                 }
             else if (kind==CDomainOfIntegration::int2d) {
             
             const QuadratureFormular & FI = QuadratureFormular_T_2;
             for (int i=0;i< Th.nt; i++) 
              {
                const Triangle & K(Th[i]);
                if (all || setoflab.find(Th[i].lab) != setoflab.end()) 
                  for (int npi=0; npi<FI.n;npi++)
                     {
                       QuadraturePoint pi(FI[npi]);
                       MeshPointStack(stack)->set(Th,K(pi),pi,K,K.lab);                       
                       r += K.area*pi.a*GetAny<double>( (*fonc)(stack)); 
                     }
               }
               }
             else   if (kind==CDomainOfIntegration::intalledges)
               {
                 const QuadratureFormular1d & FI = QF_GaussLegendre2;
                 for (int i=0;i< Th.nt; i++) 
                   if (all || setoflab.find(Th[i].lab) != setoflab.end()) 
                    for( int ie=0;ie<3;ie++)
                       {                                
                        const Triangle & K(Th[i]);   
                        R2 E=K.Edge(ie);
                        R le = sqrt((E,E)); 
                        R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
                        PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
                        
                        for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                          {
                            QuadratureFormular1d::Point pi( FI[npi]);
                            R sa=pi.x,sb=1-sa;
                            R2 Pt(PA*sa+PB*sb ); //  
                            MeshPointStack(stack)->set(Th,K(Pt),Pt,K,Th[ie].lab,R2(E.y,-E.x)/le,ie);
                            r += le*pi.a*GetAny<double>( (*fonc)(stack));
 		                   }
                       }
                    }
                else
              {
               InternalError("CDomainOfIntegration kind unkwon");
             }            
           
  *MeshPointStack(stack)=mp;
  return SetAny<double>(r);
  }
  
void Show(const char * s,int k=1)
{
  if(k) {
  couleur(1);
  float xmin,xmax,ymin,ymax;
       getcadre(xmin,xmax,ymin,ymax);
       rmoveto(xmin+(xmax-xmin)/100,ymax-(k)*(ymax-ymin)/30);
       plotstring(s);
       //  couleur(1);	
  }
}

AnyType Plot::operator()(Stack s) const  { 
   if (!withrgraphique) {initgraphique();withrgraphique=true;}
    viderbuff();
    MeshPoint *mps=MeshPointStack(s),mp=*mps ;

    R boundingbox[4];
    double coeff=1;
    bool wait=TheWait;
    bool value=false;
    bool fill=false;
    bool aspectratio=false;
    bool clean=true;
    bool uaspectratio=false;
    bool pViso=false,pVarrow=false;
    int Niso=20,Narrow=20;
    
    KN<R> Viso,Varrow;
        
    bool bw=false;
    string * psfile=0;
    string * cm=0;
    pferbase  fe=0,fe1=0;
    int cmp0,cmp1;
    bool grey=getgrey();
    bool greyo=grey;
    if (nargs[0]) coeff= GetAny<double>((*nargs[0])(s));
    if (nargs[1]) cm = GetAny<string *>((*nargs[1])(s));
    if (nargs[2]) psfile= GetAny<string*>((*nargs[2])(s));
    if (nargs[3]) wait= GetAny<bool>((*nargs[3])(s));
    if (nargs[4]) fill= GetAny<bool>((*nargs[4])(s));
    if (nargs[5]) value= GetAny<bool>((*nargs[5])(s));
    if (nargs[6]) clean= GetAny<bool>((*nargs[6])(s));
    if (nargs[7]) uaspectratio=true,uaspectratio= GetAny<bool>((*nargs[7])(s));
    if (nargs[8])  
        for (int i=0;i<4;i++)
            boundingbox[i]= GetAny<double>((*bb[i])(s));
    if (nargs[9]) Niso=  GetAny<long>((*nargs[9])(s));
    if (nargs[10]) Narrow=  GetAny<long>((*nargs[10])(s));
    if (nargs[11]) { 
          KN_<double> v =GetAny<KN_<double> >((*nargs[11])(s)) ;
          Niso=v.N();
          Viso.init(Niso);
          Viso=v;
          pViso=true;}
          
    if (nargs[12]) {
          KN_<double> v =GetAny<KN_<double> >((*nargs[12])(s)) ;
          Niso=v.N();
          Varrow.init(Niso);
          Varrow=v;
          pVarrow=true;
          }
                 
    if (nargs[13]) bw= GetAny<bool>((*nargs[13])(s));
    if (nargs[14]) grey= GetAny<bool>((*nargs[14])(s));
    setgrey(grey);
    if (Viso.unset()) Viso.init(Niso);
    if (Varrow.unset()) Varrow.init(Narrow);
    
   // KN<R> Viso(Niso);
   // KN<R> Varrow(Narrow);
   // if (pViso) Viso=*pViso;
   // if (pVarrow) Varrow=*pVarrow;
    const Mesh * cTh=0;
    bool vecvalue=false,isovalue=false;
    bool ops=psfile;
    bool drawmeshes=false;
    if ( clean ) {
         reffecran(); 

    if (psfile) {
      openPS(psfile->c_str());
    }
    if (bw) NoirEtBlanc(1);
    R2 Pmin,Pmax;
    R2 uminmax(1e100,-1e100);
    R2 Vminmax(1e100,-1e100);
    bool first=true;
     for (int i=0;i<l.size();i++)
      { R2  P1,P2;
      if (l[i].what==1 || l[i].what==2) 
      {
         
        if( !uaspectratio)   aspectratio= true;
         
         fe   = l[i].eval(0,s,cmp0);
         fe1  = l[i].eval(1,s,cmp1);
         
         if (!fe->x()) continue; 
         
         int nb=fe->x()->N();

         fe->Vh->cmesh->BoundingBox(P1,P2);
          cTh=fe->Vh->cmesh;
         if (fe1==0)
            uminmax = minmax(uminmax,fe->Vh->MinMax(*fe->x(),cmp0));
         else
          {
         if (fe1) 
          {
          if (fe->Vh == fe1->Vh) 
            {  
              KN_<R> u( *fe->x()),v(*fe1->x());     
              Vminmax = minmax(uminmax,fe->Vh->MinMax(u,v,cmp0,cmp1));
             }  
          else
           cerr << " On ne sait tracer que de vecteur sur un meme interpolation " << endl;
          }      
           
        }
      }
      else if (l[i].what==0) 
       {
        if( !uaspectratio) aspectratio= true;
        const  Mesh & Th= l[i].evalm(0,s);
         Th.BoundingBox(P1,P2);
         cTh=&Th;
       }
      else if (l[i].what==4) 
       {
        if( !uaspectratio) aspectratio= true;
        const  E_BorderN * Bh= l[i].evalb(0,s);
        Bh->BoundingBox(s,P1.x,P2.x,P1.y,P2.y);
          
       }
      else if (l[i].what==3) 
       {
         tab  ttx=l[i].evalt(0,s);
         tab  tty=l[i].evalt(1,s);
         tab *tx=&ttx,*ty=&tty;
         P1=R2(tx->min(),ty->min());
         P2=R2(tx->max(),ty->max());   
         if(verbosity>2)  
          cout << "Plot: bound  Pmin=" <<  P1 << ",  Pmax=" << P2 << endl;  
       }
      else continue;
      
         if (first)
          {first=false; Pmin=P1;Pmax=P2;}
         else {
            Pmin.x = Min(Pmin.x,P1.x);
            Pmin.y = Min(Pmin.y,P1.y);
            Pmax.x = Max(Pmax.x,P2.x);
            Pmax.y = Max(Pmax.y,P2.y);
          }
       
      }
      
      {
        R umx=uminmax.y,umn=uminmax.x;
        if (verbosity>5)
        cout << " u bound " <<  uminmax << endl;

        if (verbosity>1)
         cout << "Plot bound [x,y] " <<  Pmin << " max [x,y] " << Pmax << endl;  
        int N=Viso.N();
        int Na=Varrow.N();
        R2 O((Pmin+Pmax)/2);
        R rx(Pmax.x-Pmin.x),ry(Pmax.y-Pmin.y);
        R r =  (Max(rx,ry)*0.55);
        showgraphic();
        if (aspectratio)
          cadreortho((float)O.x,(float)(O.y+r*0.05),(float) r);
        else 
          cadre( (float)(O.x-rx*.55),(float)(O.x+rx*0.55),(float)(O.y-ry*.55),(float)(O.y+ry*.55));
        R d = fill ? (umx-umn)/(N-1)  : (umx-umn)/(N);       
        R x = fill ? umn-d/2 :umn+d/2;
       if (!pViso) 
        for (int i = 0;i < N;i++)
         {Viso[i]=x;x +=d; }
       if (fill && !pViso) {Viso[0]=umn-d;Viso[N-1]=umx+d;}
        x=0; d= sqrt(Vminmax.y)/Na;
       if (!pVarrow)
        for (int i = 0;i < Na;i++)
          {Varrow[i]=x;x +=d; }
   
        SetColorTable(Max(N,Na)+4) ;           
       }
     }  // clean
   float xx0,xx1,yy0,yy1;
   if (nargs[8])
     {
       xx0=min(boundingbox[0],boundingbox[2]);
       xx1=max(boundingbox[0],boundingbox[2]);
       yy0=min(boundingbox[1],boundingbox[3]);
       yy1=max(boundingbox[1],boundingbox[3]);
       if (verbosity>2) 
         cout << "bb=  xmin =" << xx0 << ", max =" << xx1 << ", ymin = " << yy0 << ", ymax = " << yy1 << endl;
       if (aspectratio)
         cadreortho((xx0+xx1)*0.5,(yy0+yy1)*0.5, max(xx1-xx0,yy1-yy0)*0.5);
       else 
        cadre(xx0,xx1,yy0,yy1);
     }
   getcadre(xx0,xx1,yy0,yy1);
   const R ccoeff=coeff;
   bool plotting = true;   
     //  drawing part  ------------------------------
   while (plotting)
   {
    if(verbosity>99) cout << "plot::operator() Drawing part \n";
    plotting = false; 
    for (int i=0;i<l.size();i++)
    if (l[i].what==0) 
     if (fill)
       l[i].evalm(0,s).Draw(0,fill);
     else 
        l[i].evalm(0,s).Draw(0,fill);
    else  if (l[i].what==1 || l[i].what==2)
     {
      fe=  l[i].eval(0,s,cmp0);
      fe1= l[i].eval(1,s,cmp1);;
     if (!fe->x()) continue;
      
   SHOWVERB(   cout << "   Min = " << fe->x->min() << " max = " << fe->x->max() ;
      if(fe1 && verbosity > 1)
        cout << " Min = " << fe1->x->min() << " max = " << fe1->x->max() ;   
      cout << endl;    );
      if (fe1) 
         {
          if (fe->Vh == fe1->Vh)           
           vecvalue=true,fe->Vh->Draw(*fe->x(),*fe1->x(),Varrow,coeff,cmp0,cmp1);
          else
           cerr << " On ne sait tracer que de vecteur sur un meme interpolation " << endl;
          if (drawmeshes) fe->Vh->Th.Draw(0,fill);
         }      
      else 
        
        if (fill)
         isovalue=true,fe->Vh->Drawfill(*fe->x(),Viso,cmp0);
        else 
         isovalue=true,fe->Vh->Draw(*fe->x(),Viso,cmp0);
         
        if (drawmeshes) fe->Vh->Th.Draw(0,fill);

     }
      else if (l[i].what==4) 
       {
        const  E_BorderN * Bh= l[i].evalb(0,s);
        Bh->Plot(s);          
       }
     
     else if(l[i].what==3)
      {
        
        tab x=l[i].evalt(0,s);
        tab y=l[i].evalt(1,s);
       long k= Min(x.N(),y.N());
       // cout << " à faire " << endl;
       // cout << " plot :\n" << * l[i].evalt(0,s) << endl << * l[i].evalt(1,s) << endl;
       rmoveto(x[0],y[0]);
       couleur(2+i);
       for (int i= 1;i<k;i++)
         rlineto(x[i],y[i]);
      }
      if (value) {
          int k=0; 
           if (isovalue) {PlotValue(Viso,k,"IsoValue");k+= Viso.N()+3;}
           if (vecvalue) {PlotValue(Varrow,k,"Vec Value");k+= Varrow.N()+3;}
        } 
     // value=false;
      
     if (cm) {       
       couleur(1);
       DrawCommentaire(cm->c_str(),0.1,0.97);
     }
     if (ops ) {
       ops=false;
       closePS();
     }
      if (wait && ! NoWait) 
       {  
          next:
          float x,y,x0,y0,x1,y1,dx,dy,coef=1.5;
          getcadre(x0,x1,y0,y1);
          char c=Getxyc(x,y);
          dx=(x1-x0)/2.;dy=(y1-y0)/2.;
          
          switch (c) 
           { 
             case '+' :  plotting=true; 
               cadre(x-dx/coef,x+dx/coef,y-dy/coef,y+dy/coef);reffecran();
             break;
             case '-' :  plotting=true; 
               cadre(x-dx*coef,x+dx*coef,y-dy*coef,y+dy*coef);;reffecran();
             break;
             case '=' :  plotting=true; 
               coeff=ccoeff;
               cadre(xx0,xx1,yy0,yy1);;reffecran();
             break;
             case 'r' :  plotting=true; 
               reffecran();
             break;
             case 'a' : 
             case 'c' : coeff /= 1.5; plotting=true;reffecran(); 
               reffecran();
               break;
             case 'A' : 
             case 'C' : coeff *= 1.5;
                 plotting=true;
               reffecran();
                break;
             case 'b' : bw= !bw;NoirEtBlanc(bw)  ;
                 plotting=true;
                 reffecran();
                 break;
              case 'v' : value = !value ; plotting=true;
                 reffecran();
                 break;
              case 'f' : fill = !fill ; plotting=true;
                 reffecran();
                 break;
              case 'g' : setgrey(grey=!getgrey()); plotting=true;
                 reffecran();
                 break;
                 
              case 'm' : 
                 reffecran();
                 drawmeshes= !drawmeshes;
                  plotting=true;
                 break;
              case 'p' : 
                plotting=true;
                reffecran();
                ops=true;
                openPS(0);
                plotting=true;
                break;
              case 'q' :
               couleur(8);
               if (cTh) cTh->quadtree->Draw();
               couleur(1);
                goto next;
              case 's' :
                if(cTh)
               {
                R2 P(x,y),PF(P),Phat;
                bool outside;
                const Vertex * v=cTh->quadtree->NearestVertexWithNormal(P);
                if (!v)  v=cTh->quadtree->NearestVertex(P);
                else {
                     couleur(2);
                     const Triangle * t= cTh->Find( PF,  Phat,outside,&(*cTh)[cTh->Contening(v)]) ;
                     t->Draw(0.8);
                     couleur(2);
                     PF=(*t)(Phat);
                     DrawMark(PF,0.003);
                
                   }
                 couleur(5);
                  DrawMark(P,0.0015);
                
                 couleur(1);
                 if (v)
                  DrawMark(*v,0.005);
                
                 goto next;}
              case '?':
                int i=2;
                 reffecran();               
              	Show("enter a keyboard character in graphic window to do",i++);
              	i+=2;
              	Show("+:  zoom around mouse point factor 1.5 ",i++);
              	Show("-:  unzoom around mouse point factor 1.5  ",i++);
              	Show("=:  reset zooming and unzoomin ",i++);
              	Show("r:  refrech plot ",i++);
              	Show("ac: decrease the arrow coef  ",i++);
              	Show("AC: indecrease the arrow coef  ",i++);
              	Show("b:  toggle black and white / color plotting ",i++);
              	Show("g:  toggle grey / color plotting ",i++);
              	Show("f:  toggle filling beetween iso or not  ",i++);
              	Show("v:  toggle show the numerical value of iso",i++);
              	Show("p:  save pe plot in postscrit file",i++);
              	Show("m:  show  meshes ",i++);
              	Show("t:  find  Triangles ",i++);
              	Show("?:  show this",i++);
              	Show(" otherwise continue ",i+1);
                goto next;
           }   
         if (!pViso || !pVarrow)
          { //  recompute the iso bound
             R2 uminmax(1e100,-1e100);
             R2 Vminmax(1e100,-1e100);
			 for (int i=0;i<l.size();i++)
			  { R2  P1,P2;
			  if (l[i].what==1 || l[i].what==2) 
			  {
			     fe   = l[i].eval(0,s,cmp0);
			     fe1  = l[i].eval(1,s,cmp1);
			     
			     if (!fe->x()) continue; 
			     
			     int nb=fe->x()->N();

			     
			     if (fe1==0)
			        uminmax = minmax(uminmax,fe->Vh->MinMax(*fe->x(),cmp0,false));
			     else
			      {
			     if (fe1) 
			      if (fe->Vh == fe1->Vh) 
			        {  
			          KN_<R> u( *fe->x()),v(*fe1->x());     
			          Vminmax = minmax(uminmax,fe->Vh->MinMax(u,v,cmp0,cmp1,false));
			         }  
			    }
			  }
			  else continue;
			  
			   
			  }
           if (verbosity>5)   cout << " u bound " <<  uminmax << endl;
	        R umx=uminmax.y,umn=uminmax.x;
	        int N=Viso.N();
	        int Na=Varrow.N();
	        R d = fill ? (umx-umn)/(N-1)  : (umx-umn)/(N);       
	        R x = fill ? umn-d/2 :umn+d/2;
	       if (!pViso) 
	        for (int i = 0;i < N;i++)
	         {Viso[i]=x;x +=d; }
	        if (fill && !pViso) {Viso[0]=umn-d;Viso[N-1]=umx+d;}
	        x=0; d= sqrt(Vminmax.y)/Na;
	        if (!pVarrow)
	        for (int i = 0;i < Na;i++)
	          {Varrow[i]=x;x +=d; }
          
          
          
          }
       }  
      *mps=mp;
     } //  end plotting 
     NoirEtBlanc(0)  ;
     setgrey(greyo);
     if (cm) 
       delete cm;
     if (psfile)
       delete psfile;
     viderbuff();
     
     return 0L;}     


 
AnyType Convect::operator()(Stack s) const 
{
  MeshPoint* mp(MeshPointStack(s));
  static MeshPoint mpp,mps;
  static R ddts;
  R ddt = GetAny<double>((*dt)(s));
  if (ddt) 
   {
    MeshPoint mpc(*mp);
    MeshPointStack(s,&mpc);
    if(*mp==mpp && ddt == ddts) 
      mpc=mps;
    else 
     {

            const Mesh & Th(*mp->Th);
            throwassert(mp->Th && mp->T);
            R l[3];
            l[1]=mpc.PHat.x;
            l[2]=mpc.PHat.y;
            l[0]=1-l[1]-l[2];

            int k=0;
            int j; 
            int it=Th(mpc.T);
            while ( (j=WalkInTriangle(Th,it,l,GetAny<double>((*u)(s)),GetAny<double>((*v)(s)),ddt))>=0) 
                { 
                    throwassert( l[j] == 0);
                    int jj  = j;            
                    R a= l[(j+1)%3], b= l[(j+2)%3];
                    int itt =  Th.TriangleAdj(it,j);
                    if(itt==it || itt <0)  break; // le bord 
                    it = itt;
                    l[j]=0;
                    l[(j+1)%3] = b;
                    l[(j+2)%3] = a;
                     mpc.change(R2(l[1],l[2]),Th[it],0);             
                      throwassert(k++<1000);
                }

            mpc.change(R2(l[1],l[2]),Th[it],0);
            mpp=*mp; 
            mps=mpc;         
     }
   }
  ddts=ddt;
  AnyType r= (*ff)(s);
  MeshPointStack(s,mp);
  
  return r;
}




template<class RR,class AA=RR,class BB=AA> 
struct Op2_mulAv: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  { return (*a->A * *b );} 
};

template<class RR,class AA=RR,class BB=AA> 
struct Op2_mulvirtAv: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  { return RR( (*a).A, *b );} 
};




class Op3_pfer2R : public ternary_function<pfer,R,R,R> { public:


  class Op : public E_F0mps { public:
      Expression a,b,c;
       Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {}       
       AnyType operator()(Stack s)  const 
        { 
           R xx(GetAny<R>((*b)(s)));
           R yy(GetAny<R>((*c)(s)));
           MeshPoint & mp = *MeshPointStack(s),mps=mp;
           mp.set(xx,yy,0.0);
           AnyType ret = pfer2R<R,0>(s,(*a)(s));
           mp=mps;
           return  ret;}
   
  };
};

template<class K>
class Op3_K2R : public ternary_function<K,R,R,R> { public:

  class Op : public E_F0mps { public:
      Expression a,b,c;
       Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {}       
       AnyType operator()(Stack s)  const 
        { 
           R xx(GetAny<R>((*b)(s)));
           R yy(GetAny<R>((*c)(s)));
           MeshPoint & mp = *MeshPointStack(s),mps=mp;
           mp.set(xx,yy,0.0);
           AnyType ret = (*a)(s);
           mp=mps;
           return  ret;}
   
  };
};

template<class RR,class AA=RR>
class JumpOp : public E_F0mps  { public:
  typedef RR  R;
  typedef AA A; 
  typedef RR result_type;
  typedef AA argument_type; 

       Expression a;
       public:
       AnyType operator()(Stack stack)  const 
        { // a faire 
           A rd,rg;       
           MeshPoint *mp=MeshPointStack(stack),smp=*mp;              
           rg = GetAny<A>((*a)(stack));
           rd = rg;
           if ( mp->SetAdj() )
             rd = GetAny<A>((*a)(stack));  
           *mp=smp;         
           return  SetAny<R>(rg-rd);
        }
       JumpOp(Expression aa) : a(aa) {} 
    };

template<class RR,class AA=RR>
class MeanOp : public E_F0mps  { public:
  typedef RR  R;
  typedef AA A; 
  typedef RR result_type;
  typedef AA argument_type; 

       Expression a;
       public:
       AnyType operator()(Stack stack)  const 
        { // a faire 
           A rd,rg;       
           MeshPoint *mp=MeshPointStack(stack),smp=*mp;              
           rg = GetAny<A>((*a)(stack));
           rd = rg;
           if ( mp->SetAdj() )
             rd = GetAny<A>((*a)(stack));  
           *mp=smp;         
           return  SetAny<R>(rg-rd);
        }
       MeanOp(Expression aa) : a(aa) {} 
    };

pferbase* get_element(pferbasearray *const & a, long const & n)
{
  return (**a)[n];
}
pfer get_element(pferarray const & a, long const & n)
{
  return pfer( *(*a.first)[n],a.second);
}

lgElement get_element(pmesh const & a, long const & n)
{
  return lgElement(a,n);
}
lgElement get_element(pmesh *const & a, long const & n)
{
  return lgElement(*a,n);
}
lgVertex get_element(lgElement const & a, long const & n)
{
  return a[n];
}
R getx(lgVertex const & a)
{
  return a.x();
}

R gety(lgVertex const & a)
{
  return a.y();
}
long getlab(lgVertex const & a)
{
  return a.lab();
}

/*
pmesh* get_element(pmesharray *const & a, long const & n)
{
  assert(n>=0 && n <=a->second);
  return a->first +n;
}
*/
template<class A> inline AnyType DestroyKN(Stack,const AnyType &x){
  KN<A> * a=GetAny<KN<A>*>(x);
  for (int i=0;i<a->N(); i++)
    (*a)[i]->destroy();
  a->destroy(); 
  return  Nothing;
}
template<class RR,class A,class B>  
RR * get_elementp_(const A & a,const B & b){ 
  if( b<0 || a->N() <= b) 
   { cerr << " Out of bound  0 <=" << b << " < "  << a->N() << " array type = " << typeid(A).name() << endl;
     ExecError("Out of bound in operator []");}
    return  &((*a)[b]);}
    
template<class R>  R * set_initinit( R* const & a,const long & n){ 
 SHOWVERB( cout << " set_init " << typeid(R).name() << " " << n << endl);
  a->init(n);
  for (int i=0;i<n;i++)
    (*a)[i]=0;
   return a;}
    
void init_mesh_array()
 {
  Dcl_Type<KN<pmesh> *>(0,::DestroyKN<pmesh> );
  atype<KN<pmesh>* >()->Add("[","",new OneOperator2_<pmesh*,KN<pmesh>*,long >(get_elementp_<pmesh,KN<pmesh>*,long>));
    TheOperators->Add("<-", 
       new OneOperator2_<KN<pmesh> *,KN<pmesh> *,long>(&set_initinit));
  map_type_of_map[make_pair(atype<long>(),atype<pmesh>())]=atype<KN<pmesh>*>(); // vector

 }

void  init_lgfem() 
{
  cout <<"lg_fem ";
//Dcl_Type<const C_args*>(); // to store compilation expression
 
 Dcl_Type<MeshPoint*>();
 Dcl_Type<R3*>(::Initialize<R3>);
 Dcl_Type<R2*>(::Initialize<R2>);

 map_type[typeid(R3*).name()] = new ForEachType<R3*>(Initialize<R3>);   
  Dcl_TypeandPtr<pmesh>(); 
  Dcl_Type<lgVertex>(); 
  Dcl_Type<lgElement>( ); 

  atype<long>()->AddCast( 
    new E_F1_funcT<long,lgVertex>(Cast<long,lgVertex>),
    new E_F1_funcT<long,lgElement>(Cast<long,lgElement>)
  );


 Dcl_Type<TypeOfFE*>(); 
 Dcl_Type<TypeSolveMat*>();
 Dcl_Type<Matrice_Creuse<R>* >(InitP<Matrice_Creuse<R> >,Destroy<Matrice_Creuse<R> >);
  Dcl_Type<Matrice_Creuse_Transpose<R> >();      // matrice^t   (A')                             
  Dcl_Type<Matrice_Creuse_inv<R> >();      // matrice^-1   A^{-1}                          
  Dcl_Type<VirtualMatrice<R>::plusAx >();       // A*x (A'*x)
  Dcl_Type<VirtualMatrice<R>::plusAtx >();       // A^t*x (A'*x)
  Dcl_Type<VirtualMatrice<R>::solveAxeqb >();       // A^t*x (A'*x)
  Dcl_Type<const MatrixInterpolation::Op *>(); 
//  Dcl_Type<MatriceCreuseDivKN_<double> >();      //  A^(-1)  x
                                    
                                    
 Dcl_TypeandPtr<pferbase>(); // il faut le 2 pour pourvoir initialiser 
 Dcl_TypeandPtr<pferbasearray>(); // il faut le 2 pour pourvoir initialiser 
 Dcl_Type< pfer >(); 
 Dcl_Type< pferarray >(); 
// Dcl_Type< pmesharray *>(); // il faut le 2 pour pourvoir initialiser 

 map_type[typeid(pfes).name()] = new ForEachType<pfes>(); 
 map_type[typeid(pfes*).name()] = new ForEachTypePtrfspace<pfes>();
 //   
 Dcl_Type<const QuadratureFormular *>();
 Dcl_Type<const QuadratureFormular1d *>();
 Global.New("qf1pT",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_1));
 Global.New("qf1pTlump",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_1lump));
 Global.New("qf2pT",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_2));
 Global.New("qf2pT4P1",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_2_4P1));
 Global.New("qf5pT",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_5));
 
 Global.New("qf1pE",CConstant<const QuadratureFormular1d *>(&QF_GaussLegendre1));
 Global.New("qf2pE",CConstant<const QuadratureFormular1d *>(&QF_GaussLegendre2));
 Global.New("qf3pE",CConstant<const QuadratureFormular1d *>(&QF_GaussLegendre3));
 
 //  juste du code genere 
 
 Global.New("wait",CConstant<bool*>(&TheWait));
 Global.New("NoUseOfWait",CConstant<bool*>(&NoWait));

 Dcl_Type<finconnue *>();
 Dcl_Type<ftest *>();
 Dcl_Type<foperator *>();
 Dcl_Type<foperator *>();
 Dcl_Type<const BC_set<R> *>();  // a set of boundary condition 
 Dcl_Type<const Call_FormLinear *>();    //   to set Vector
 Dcl_Type<const Call_FormBilinear *>();  // to set Matrix
 Dcl_Type<interpolate_f_X_1<double>::type>();  // to make  interpolation x=f o X^1 ;
 
 map_type[typeid(const FormBilinear*).name()] = new TypeFormBilinear;
 map_type[typeid(const FormLinear*).name()] = new TypeFormLinear;
 aType t_C_args = map_type[typeid(const C_args*).name()] = new TypeFormOperator;
 map_type[typeid(const Problem*).name()] = new TypeSolve<false,Problem>;
 map_type[typeid(const Solve*).name()] = new TypeSolve<true,Solve>;
 Dcl_Type<const IntFunction*>(); 
 basicForEachType * t_solve=atype<const  Solve *>();
 basicForEachType * t_problem=atype<const  Problem *>();
 basicForEachType * t_fbilin=atype<const  FormBilinear *>();
 basicForEachType * t_flin=atype<const  FormLinear *>();
 basicForEachType * t_BC=atype<const BC_set<R> *>();
 basicForEachType * t_form=atype<const C_args*>();

  Dcl_Type<const CDomainOfIntegration *>();
  
 

 atype<pmesh >()->AddCast( new E_F1_funcT<pmesh,pmesh*>(UnRef<pmesh >)); 
 atype<pfes >()->AddCast(  new E_F1_funcT<pfes,pfes*>(UnRef<pfes>));
 atype<pferbase>()->AddCast(  new E_F1_funcT<pferbase,pferbase>(UnRef<pferbase>));
 
 Add<pfer>("[]",".",new OneOperator1<KN<double> *,pfer>(pfer2vect<R>));
 Add<pfer>("(","",new OneTernaryOperator<Op3_pfer2R,Op3_pfer2R::Op> );
 Add<double>("(","",new OneTernaryOperator<Op3_K2R<R>,Op3_K2R<R>::Op> );
 
 Add<pfer>("n",".",new OneOperator1<long,pfer>(pfer_nbdf<R>));
 Add<pmesh*>("area",".",new OneOperator1<double,pmesh*>(pmesh_area));
 Add<pmesh*>("nt",".",new OneOperator1<long,pmesh*>(pmesh_nt));
 Add<pmesh*>("nv",".",new OneOperator1<long,pmesh*>(pmesh_nv));
 atype<Matrice_Creuse<R> * >()->AddCast(new OneOperatorCode<pb2mat<R> >);
 
TheOperators->Add("*", 
        new OneBinaryOperator<Op2_mulvirtAv<VirtualMatrice<R>::plusAx,Matrice_Creuse<R>*,KN<R>* > >,
        new OneBinaryOperator<Op2_mulvirtAv<VirtualMatrice<R>::plusAtx,Matrice_Creuse_Transpose<R>,KN<R>* > >,
        new OneBinaryOperator<Op2_mulvirtAv<VirtualMatrice<R>::solveAxeqb,Matrice_Creuse_inv<R>,KN<R>* > >
        );
TheOperators->Add("^", new OneBinaryOperatorA_inv<R>());
//   Add all Finite Element "P0","P1","P2","RT0", ... 
  for (ListOfTFE * i=ListOfTFE::all;i;i=i->next)
    {
     assert(i->tfe); // check 
     Global.New(i->name, Type_Expr(atype<TypeOfFE*>(),new  EConstantTypeOfFE(i->tfe)));
    }
   
// Global.New("P1",CConstant<TypeOfFE*>(&P1Lagrange));
// Global.New("P2",CConstant<TypeOfFE*>(&P2Lagrange));
 
 Global.New("LU",CConstant<TypeSolveMat*>(new TypeSolveMat(TypeSolveMat::LU))); 
 Global.New("CG",CConstant<TypeSolveMat*>(new TypeSolveMat(TypeSolveMat::GC)));
 Global.New("Crout",CConstant<TypeSolveMat*>(new TypeSolveMat(TypeSolveMat::CROUT)));
 Global.New("Cholesky",CConstant<TypeSolveMat*>(new TypeSolveMat(TypeSolveMat::CHOLESKY)));
 Global.New("GMRES",CConstant<TypeSolveMat*>(new TypeSolveMat(TypeSolveMat::GMRES)));
#ifdef UMFPACK
 Global.New("UMFPACK",CConstant<TypeSolveMat*>(new TypeSolveMat(TypeSolveMat::UMFpack)));
#else 
 cout << " --  no UMFPACK => replace UMFPACK  by LU  " << endl;
 Global.New("UMFPACK",CConstant<TypeSolveMat*>(new TypeSolveMat(TypeSolveMat::LU)));
#endif 
//  init pmesh  
 Add<pmesh*>("<-","(",
             new OneOperator1_<pmesh,string*>(ReadMesh),
             new OneOperator3_<long,pmesh*,double,double,
                      E_F_stackF0F0F0_<long,pmesh*,double,double> >(& FindTxy )

 );
 
/*TheOperators->Add("<-",);

 );*/
 Add<const MatrixInterpolation::Op *>("<-","(", new MatrixInterpolation);
 
//old --   
//  init FESpace 
   TheOperators->Add("<-",
       new OneOperator2_<pfes*,pfes*,pmesh* >(& MakePtr2 ),
   //    new OneOperator3_<pfes*,pfes*,pmesh*,long >(& MakePtr2 ),
   //    new OneOperator3_<pfes*,pfes*,pmesh*,TypeOfFE* >(& MakePtr2 ),
       new OneOperatorCode<OP_MakePtr2>,
       new OpMake_pfes
        );
      
 
 Add<MeshPoint*>("P",".", new OneOperator_Ptr_o_R<R3,MeshPoint>(  & MeshPoint::P));
 Add<MeshPoint*>("N",".", new OneOperator_Ptr_o_R<R3,MeshPoint>(  & MeshPoint::N));
 Add<R3*>("x",".", new OneOperator_Ptr_o_R<R,R3>(  & R3::x));
 Add<R3*>("y",".", new OneOperator_Ptr_o_R<R,R3>(  & R3::y));
 Add<R3*>("z",".", new OneOperator_Ptr_o_R<R,R3>(  & R3::z));
 Add<R2*>("x",".", new OneOperator_Ptr_o_R<R,R2>(  & R2::x));
 Add<R2*>("y",".", new OneOperator_Ptr_o_R<R,R2>(  & R2::y));
 
 Add<pmesh>("[","",new OneOperator2_<lgElement,pmesh,long>(get_element));
 Add<pmesh*>("[","",new OneOperator2_<lgElement,pmesh*,long>(get_element));

 Add<lgElement>("[","",new OneOperator2_<lgVertex,lgElement,long>(get_element));
 Add<lgVertex>("x",".",new OneOperator1_<R,lgVertex>(getx));
 Add<lgVertex>("y",".",new OneOperator1_<R,lgVertex>(gety));
 Add<lgVertex>("label",".",new OneOperator1_<long,lgVertex>(getlab));
 
 
//  new type  
 zzzfff->Add("R3",atype<R3*>());
 zzzfff->Add("mesh",atype<pmesh*>());
 zzzfff->Add("element",atype<lgElement>());
 zzzfff->Add("vertex",atype<lgVertex>());
// zzzfff->Add("fespace",atype<pfes*>());
 zzzfff->Add("matrix",atype<Matrice_Creuse<R> *>());

 
 Global.Add("buildmesh","(",new OneOperator1s_<pmesh,const E_BorderN *>(BuildMesh));
 Global.Add("buildmeshborder","(",new OneOperator1s_<pmesh,const E_BorderN *>(BuildMeshBorder));
 Global.Add("LinearCG","(",new LinearCG<R>()); // old form  with rhs (must be zer
 Global.Add("LinearGMRES","(",new LinearGMRES<R>()); // old form  with rhs (must be zer
 Global.Add("LinearGMRES","(",new LinearGMRES<R>(1)); // old form  with rhs (must be zer
 Global.Add("LinearCG","(",new LinearCG<R>(1)); //  without right handsize
 Global.Add("NLCG","(",new LinearCG<R>(-1)); //  without right handsize
  zzzfff->AddF("varf",t_form);    //  var. form ~
 zzzfff->AddF("solve",t_solve);
 zzzfff->AddF("problem",t_problem);
 

 Global.Add("jump","(",new OneOperatorCode<Code_VF<Ftest,Code_Jump> >);
 Global.Add("jump","(",new OneOperatorCode<Code_VF<Finconnue,Code_Jump> >);
 Global.Add("average","(",new OneOperatorCode<Code_VF<Ftest,Code_Mean> >);
 Global.Add("average","(",new OneOperatorCode<Code_VF<Finconnue,Code_Mean> >);
 Global.Add("otherside","(",new OneOperatorCode<Code_VF<Ftest,Code_OtherSide> >);
 Global.Add("otherside","(",new OneOperatorCode<Code_VF<Finconnue,Code_OtherSide> >);
 
 Global.Add("dx","(",new OneOperatorCode<CODE_Diff<Ftest,op_dx> >);
 Global.Add("dy","(",new OneOperatorCode<CODE_Diff<Ftest,op_dy> >);
 Global.Add("dx","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dx> >);
 Global.Add("dy","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dy> >);
 
 Global.Add("dxx","(",new OneOperatorCode<CODE_Diff<Ftest,op_dxx> >);
 Global.Add("dxy","(",new OneOperatorCode<CODE_Diff<Ftest,op_dxy> >);
 Global.Add("dyx","(",new OneOperatorCode<CODE_Diff<Ftest,op_dyx> >);
 Global.Add("dyy","(",new OneOperatorCode<CODE_Diff<Ftest,op_dyy> >);
 
 Global.Add("dxx","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dxx> >);
 Global.Add("dyy","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dyy> >);
 Global.Add("dxy","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dxy> >);
 Global.Add("dyx","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dyx> >);
 
 Global.Add("on","(",new OneOperatorCode<BC_set<double> > );
 Global.Add("plot","(",new OneOperatorCode<Plot> );
 Global.Add("convect","(",new OneOperatorCode<Convect> );


 TheOperators->Add("+",
    new OneOperatorCode<CODE_L_Add<Foperator> > ,
    new OneOperatorCode<CODE_L_Add<Ftest> > ,
    new OneOperatorCode<CODE_L_Add<Finconnue> > ,
    new OneOperatorCode<C_args>(t_C_args,t_C_args,t_C_args)        
  );
 TheOperators->Add("-",
    new OneOperatorCode<CODE_L_Minus<Foperator> > ,
    new OneOperatorCode<CODE_L_Minus<Ftest> > ,
    new OneOperatorCode<CODE_L_Minus<Finconnue> >,
    new OneOperatorCode<CODE_L_Sub<Foperator> > ,    
    new OneOperatorCode<CODE_L_Sub<Ftest> > ,
    new OneOperatorCode<CODE_L_Sub<Finconnue> >,
    new OneOperatorCode<C_args_minus>(t_C_args,t_C_args,t_fbilin) ,       
    new OneOperatorCode<C_args_minus>(t_C_args,t_C_args,t_flin) ,
    new OneOperatorCode<Minus_Form<FormBilinear > >,
    new OneOperatorCode<Minus_Form<FormLinear > > 
           
  );
  
  
  atype<const C_args*>()->AddCast( 
      new OneOperatorCode<C_args>(t_C_args,t_fbilin) ,       
      new OneOperatorCode<C_args>(t_C_args,t_flin)  ,     
      new OneOperatorCode<C_args>(t_C_args,t_BC) ,      
      new OneOperatorCode<C_args>(t_C_args,atype<KN<R>*>())  ,
      new OneOperatorCode<C_args>(t_C_args,atype<DotStar_KN_<R> >())  ,           
      new OneOperatorCode<C_args>(t_C_args,atype<Matrice_Creuse<R>*>()) ,
      new OneOperatorCode<C_args>(t_C_args,atype<VirtualMatrice<R>::plusAx >()),
      new OneOperatorCode<C_args>(t_C_args,atype<VirtualMatrice<R>::plusAtx >())   
    );
    
  atype<const C_args*>()->AddCast( 
      new OneOperatorCode<C_args>(t_C_args,atype<DotSlash_KN_<R> >())  
    );         
    
 TheOperators->Add("*",  
    new OneOperatorCode<CODE_L_Mul<Foperator,Ftest,Finconnue> > ,
    new OneOperatorCode<CODE_L_Mul<Foperator,Finconnue,Ftest> > ,
    new OneOperatorCode<CODE_L_MulLR<Finconnue,double> > ,
    new OneOperatorCode<CODE_L_MulLR<Foperator,double> > ,
    new OneOperatorCode<CODE_L_MulLR<Ftest,double> > ,
    new OneOperatorCode<CODE_L_MulRL<double,Finconnue> > ,
    new OneOperatorCode<CODE_L_MulRL<double,Foperator> > ,
    new OneOperatorCode<CODE_L_MulRL<double,Ftest> > 
  );
 TheOperators->Add("/",  
    new OneOperatorCode<CODE_L_DivLR<Finconnue,double> > ,
    new OneOperatorCode<CODE_L_DivLR<Foperator,double> > ,
    new OneOperatorCode<CODE_L_DivLR<Ftest,double> > ),
  
// TheOperators->Add("=",new OneOperatorCode<BC_set1<double> >);

  
 TheOperators->Add("=",
       new OneOperator2<pmesh*,pmesh*,pmesh >(&set_eqdestroy_incr),
       new OneBinaryOperator<set_eqarray<KN<double> ,VirtualMatrice<R>::plusAx > > ,       
       new OneBinaryOperator<set_eqarray<KN<double> ,VirtualMatrice<R>::plusAtx > >  ,      
       new OneBinaryOperator<set_eqarray<KN<double> ,VirtualMatrice<R>::solveAxeqb > >  ,      
       new OpArraytoLinearForm<double >  ,
       new OpMatrixtoBilinearForm<double >,
       new OneOperator2_<Matrice_Creuse<R>*,Matrice_Creuse<R>*,const MatrixInterpolation::Op*,E_F_StackF0F0>(SetMatrixInterpolation)
       );
       
 TheOperators->Add("+=",
       new OneBinaryOperator<set_eqarray_add<KN<double> ,VirtualMatrice<R>::plusAx > > ,       
       new OneBinaryOperator<set_eqarray_add<KN<double> ,VirtualMatrice<R>::plusAtx > >        
       );


       
 TheOperators->Add("<-",new OpMatrixtoBilinearForm<double > );
       
 Add<const  FormLinear   *>("(","",new OpCall_FormLinear<FormLinear> );
 Add<const  FormBilinear *>("(","",new OpCall_FormBilinear<FormBilinear> );
 Add<const  FormBilinear *>("(","",new OpCall_FormLinear2<FormBilinear> );
 Add<const C_args*>("(","",new OpCall_FormLinear2<C_args>);
 Add<const C_args*>("(","",new OpCall_FormBilinear<C_args> );

//  correction du bug morale 
//  Attention il y a moralement un bug
//  les initialisation   x = y   ( passe par l'operateur binaire <-  dans TheOperators
//   les initialisation   x(y)   ( passe par l'operateur unaire <-  de typedebase de x
//  x(y1,..,yn) est un operator n+1   (x,y1,..,yn)
// on passe toujours par x(y) maintenant.
//   -------
       
       
 TheOperators->Add("<-",
       new OneOperator2_<pferbase*,pferbase*,pfes* >(MakePtrFE2),
       new OneOperator3_<pferbasearray*,pferbasearray*,pfes*,long >(MakePtrFE3) //,
     //  new OneOperator2_<pmesharray*,pmesharray*,long >(MakePtr)
       
       
       );
 TheOperators->Add("<-",
       new OneOperatorMakePtrFE(atype<double>()),  //  scalar case
       new OneOperatorMakePtrFE(atype<E_Array>())  //  vect case
       );
//  interpolation   operator 
 TheOperators->Add("=",
       new OneOperator2_<pfer,pfer,double,E_F_StackF0F0opt2 >(set_fe) ) ;     
     //  new OneOperator2_<pferbase,pferbase,E_Array,E_F_StackF0F0 >(set_fev));
       
//  Attention il y a moralement un bug
//  les initialisation   x = y   ( passe par l'operateur binaire <-  dans TheOperators
//   les initialisation   x(y)   ( passe par l'operateur unaire <-  de typedebase de x
//   -------  corrige 

     
 TheOperators->Add("=",
       new OneOperator2<pfes*,pfes*,pfes>(&set_eqdestroy_incr),
       new Op_CopyArray()
       );
       
 TheOperators->Add("<-",
       new OneOperator2_<pfes*,pfes*,pfes>(&set_copy_incr));
 
 TheOperators->Add("<<",
       new OneBinaryOperator<PrintPnd<R3*>  >,
       new OneBinaryOperator<PrintPnd<Matrice_Creuse<R>*> >

       );   
 
 Global.Add("int2d","(",new OneOperatorCode<CDomainOfIntegration>);
 Global.Add("int1d","(",new OneOperatorCode<CDomainOfIntegrationBorder>);
 Global.Add("intalledges","(",new OneOperatorCode<CDomainOfIntegrationAllEdges>);
 Global.Add("jump","(",new OneUnaryOperator<JumpOp<R>,JumpOp<R> >);
 Global.Add("mean","(",new OneUnaryOperator<MeanOp<R>,MeanOp<R> >);

 Add<const CDomainOfIntegration*>("(","",new OneOperatorCode<FormBilinear> );
 Add<const CDomainOfIntegration *>("(","",new OneOperatorCode<FormLinear> );
 
 Add<const CDomainOfIntegration *>("(","",new OneOperatorCode<IntFunction>);

 map_type[typeid(double).name()]->AddCast(
   new E_F1_funcT<double,pfer>(pfer2R<R,0>)
   );
 
// bof  
 Global.Add("dx","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dx>));
 Global.Add("dy","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dy>));
 Global.Add("dxx","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dxx>));
 Global.Add("dyy","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dyy>));
 Global.Add("dxy","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dxy>));
 Global.Add("dyx","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dyx>));

 
  Add<pferbasearray*>("[","",new OneOperator2_<pferbase*,pferbasearray*,long>(get_element));
  Add<pferarray>("[","",new OneOperator2_<pfer,pferarray,long>(get_element));
 // Add<pmesharray>("[","",new OneOperator2_<pmesh*,pmesharray*,long>(get_element));

  TheOperators->Add("\'",       
       new OneOperator1<Matrice_Creuse_Transpose<R>,Matrice_Creuse<R> *>(&Build<Matrice_Creuse_Transpose<R>,Matrice_Creuse<R> *>)
  );
  
 Add<pfer>("(","",new interpolate_f_X_1<R> ); 
  TheOperators->Add("=", new OneOperator2_<void,interpolate_f_X_1<R>::type,double,E_F_StackF0F0 >(set_feoX_1) ) ; 
// init_lgmesh() ;
  init_mesh_array();

 l2interpreter = new LinkToInterpreter;
        

}   


template<class K>
Expression IsFEcomp(const C_F0 &c,int i)
{
  if(atype<typename E_FEcomp<K>::Result>() == c.left())
   {
     const E_FEcomp<K> * e=dynamic_cast<const E_FEcomp<R>*>(c.LeftValue() );
     throwassert(e);
     if (e->comp !=i) return 0;
     else return e->a0;
   }
  else return 0;
}


 E_F0 * Op_CopyArray::code(const basicAC_F0 & args) const  {
       E_F0 * ret=0;
      const E_Array & a= *dynamic_cast<const E_Array*>(args[0].LeftValue());
      const E_Array & b= *dynamic_cast<const E_Array*>(args[1].LeftValue());
      int na=a.size();
      int nb=b.size();
      if (na != nb ) 
        CompileError("Size Of Array incompatible");
        const  E_F0 * rr=0;

          for (int i=0;i<nb;i++)
           {
            if ( rr ==0 ) {rr=IsFEcomp<R>(a[i],i) ; if (rr==0) goto general;}
             else  if (rr!= IsFEcomp<R>(a[i],i))  goto general;
           }
         return new E_set_fev(&b,rr);
 general:
        CompileError("Internal Error: General Copy of Array : to do ");
      return ret;}

C_F0 NewFEvariable(const char * id,Block *currentblock,C_F0 & fespacetype,CC_F0 init)
{
  ListOfId lid;
  lid.push_back(UnId(id));
  return NewFEvariable(&lid,currentblock,fespacetype,init);
}

C_F0 NewFEvariable(ListOfId * pids,Block *currentblock,C_F0 & fespacetype,CC_F0 init)
{
   typedef FEbase<double> FE;
   typedef E_FEcomp<R> FEi;
   typedef E_FEcomp<R>::Result FEiR;
   
   Expression fes =fespacetype;
    aType dcltype=atype<FE **>();
    aType cf0type=atype<C_F0>();
    aType rtype=atype<FEiR>();
    throwassert(pids);
    ListOfId &ids(*pids);
    
    string str("[");
    
    const int n=ids.size();
     throwassert(n>0);
   if ( fes->nbitem() != n) {
      cerr << " the array size must be " << fes->nbitem()  << " not " <<  n << endl;
      CompileError("Invalide array size  for  vectorial fespace function");
   }
   for (int i=0;i<n;i++)
    { 
     str += ids[i].id;
     if(i<n-1) str +=",";
    }
     str += "]";
     bool binit= !init.Empty(); 
     char * name = strcpy(new char[str.size()+1],str.c_str());
     C_F0 ret= binit ? currentblock->NewVar<LocalVariable>(name,dcltype,basicAC_F0_wa(fespacetype,init))
                     : currentblock->NewVar<LocalVariable>(name,dcltype,basicAC_F0_wa(fespacetype));
     C_F0 base = currentblock->Find(name);
     
      for (int i=0;i<n;i++) 
         currentblock->NewID(cf0type,ids[i].id, C_F0(new FEi(base,i,n), rtype) ); 

      return ret ; 
}

 size_t dimFESpaceImage(const basicAC_F0 &args) 
{
   aType t_tfe= atype<TypeOfFE*>();
   aType t_a= atype<E_Array>();
   size_t dim=0;
   for (int i=0;i<args.size();i++)
     if (args[i].left() == t_tfe)
       dim += args[i].LeftValue()->nbitem();
     else if (args[i].left() == t_a)
       {
         const E_Array & ea= *dynamic_cast<const E_Array *>(args[i].LeftValue());
         assert(&ea);
         for (int i=0;i<ea.size();i++)
           if (ea[i].left() == t_tfe)
            dim += ea[i].nbitem();
       }
    dim = dim ? dim : 1;
    // cout << "dimFESpaceImage:  FESpace in R^"<< dim << endl;
    return dim;
}

C_F0 NewFEarray(const char * id,Block *currentblock,C_F0 & fespacetype,CC_F0 sizeofarray)
{ 
  ListOfId lid;
  lid.push_back(UnId(id));
  return NewFEarray(&lid,currentblock,fespacetype,sizeofarray);
}
   
C_F0 NewFEarray(ListOfId * pids,Block *currentblock,C_F0 & fespacetype,CC_F0 sizeofarray)
{
   typedef FEbaseArray<double>  FE;
   typedef  E_FEcomp<R,FE > FEi;
   typedef FEi::Result FEiR;
   
   Expression fes =fespacetype;
    aType dcltype=atype<FE **>();
    aType cf0type=atype<C_F0>();
    aType rtype=atype<FEiR>();
    throwassert(pids);
    ListOfId &ids(*pids);
    
    string str("[");
    
    const int n=ids.size();
     throwassert(n>0);
   if ( fes->nbitem() != n) {
      cerr << " the array size must be " << fes->nbitem()  << " not " <<  n << endl;
      CompileError("Invalide array size  for  vectorial fespace function");
   }
   for (int i=0;i<n;i++)
    { 
     str += ids[i].id;
     if(i<n-1) str +=",";
    }
     str += "]";
     char * name = strcpy(new char[str.size()+1],str.c_str());
     C_F0 ret=  currentblock->NewVar<LocalVariable>(name,dcltype,basicAC_F0_wa(fespacetype,sizeofarray)); 
     C_F0 base = currentblock->Find(name);
     
      for (int i=0;i<n;i++) 
         currentblock->NewID(cf0type,ids[i].id, C_F0(new FEi(base,i,n), rtype) ); 

      return ret ; 

}
#include "InitFunct.hpp"
static addingInitFunct TheaddingInitFunct(-20,init_lgfem);

