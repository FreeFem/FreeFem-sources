/// \file
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
#ifdef __MWERKS__
#pragma optimization_level 0
#endif

#include  <cmath>
#include  <iostream>
#include <cfloat>

using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include <cstdio>
#include "fem.hpp"
#include "Mesh3dn.hpp"

#include "HashMatrix.hpp"

#include "SparseLinearSolver.hpp"

#include "MeshPoint.hpp"
#include <complex>
#include "Operator.hpp"

#include <set>
#include <map>
#include <vector>

#include "lex.hpp"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include "CGNL.hpp"
#include "AddNewFE.h"
#include "array_resize.hpp"
#include "PlotStream.hpp"

// add for the gestion of the endianness of the file.
//PlotStream::fBytes PlotStream::zott; //0123;
//PlotStream::hBytes PlotStream::zottffss; //012345678;
// ---- FH
namespace bamg { class Triangles; }
namespace Fem2D { void DrawIsoT(const R2 Pt[3],const R ff[3],const RN_ & Viso);
   extern GTypeOfFE<Mesh3> &P1bLagrange3d;
   extern GTypeOfFE<Mesh3> &RT03d;
   extern GTypeOfFE<Mesh3> &Edge03d;
   extern GTypeOfFE<MeshS> &P1bLagrange_surf;
void  Expandsetoflab(Stack stack,const CDomainOfIntegration & di,set<int> & setoflab,bool &all);
}

#include "BamgFreeFem.hpp"

static bool TheWait=false;
bool  NoWait=false;
extern bool  NoGraphicWindow;

extern long verbosity;
extern FILE *ThePlotStream; //  Add for new plot. FH oct 2008
void init_lgmesh() ;

namespace FreeFempp {
template<class R>
 TypeVarForm<R> * TypeVarForm<R>::Global;
}

basicAC_F0::name_and_type  OpCall_FormBilinear_np::name_param[]= {
{   "bmat",&typeid(Matrice_Creuse<R>* )},
     LIST_NAME_PARM_MAT
};


basicAC_F0::name_and_type  OpCall_FormLinear_np::name_param[]= {
  "tgv",&typeid(double )
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


bool In(long *viso,int n,long v)
{
  int i=0,j=n,k;
  if  (v <viso[0] || v >viso[j-1])
    return false;
  while (i<j-1)
   if ( viso[k=(i+j)/2]> v) j=k;
   else i=k;
  return (viso[i]=v);
}


class  LinkToInterpreter { public:
 Type_Expr   P,N,x,y,z,label,region,nu_triangle,nu_face,nu_edge,lenEdge,hTriangle,area,inside,volume;
  LinkToInterpreter() ;
};

LinkToInterpreter * l2interpreter;

  using namespace Fem2D;
  using namespace EF23;

template<class Result,class A>
class E_F_A_Ptr_o_R :public  E_F0 { public:
  typedef Result A::* ptr;
  Expression a0;
  ptr p;
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
    ffassert(mp->T && mp ->e >=0 && mp->d==2);
    double l= mp->T->lenEdge(mp->e);
    return SetAny<double>(l);}
    operator aType () const { return atype<double>();}

};

class E_P_Stack_hTriangle   :public  E_F0mps { public:
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    MeshPoint * mp=MeshPointStack(s);
    assert(mp->T) ;
    double l=1e100;
    if( mp->d==2) l=mp->T->h();
    else if ( mp->d==3 && mp->dHat==3 ) l=mp->T3->lenEdgesmax();
    else if ( mp->d==3 && mp->dHat==2 ) l=mp->TS->lenEdgesmax();
    return SetAny<double>(l);}
    operator aType () const  { return atype<double>();}

};

class E_P_Stack_nTonEdge   :public  E_F0mps { public:
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    MeshPoint * mp=MeshPointStack(s);
    assert(mp->T && mp->e > -1 && mp->d==2 ) ;
    long l=mp->Th->nTonEdge(mp->t,mp->e);
    // cout << " nTonEdge " << l << endl;
    return SetAny<long>( l) ;}
    operator aType () const  { return atype<long>();}

};

class E_P_Stack_nElementonB   :public  E_F0mps { public:
    AnyType operator()(Stack s)  const { throwassert(* (long *) s);
      MeshPoint * mp=MeshPointStack(s);
      long l=0;
      if((mp->T) && (mp->e > -1) && (mp->d==2 ))
         l=mp->Th->nTonEdge(mp->t,mp->e);
      else if (mp->d==3 && mp->dHat==3 && mp->T3 && ( mp->f>=0) )
         l= mp->Th3->nElementonB(mp->t,mp->f);
      else if (mp->d==3 && mp->dHat==2 && mp->TS && ( mp->e>=0) )
          l= mp->ThS->nElementonB(mp->t,mp->e);
      else if (mp->d==3 && mp->dHat==1 && mp->TL && ( mp->e>=0) )
          l= mp->ThL->nElementonB(mp->t,mp->e);
      // cout << " nTonEdge " << l << endl;
      return SetAny<long>( l) ;}
    operator aType () const  { return atype<long>();}

};

template<int NBT>
class E_P_Stack_TypeEdge   :public  E_F0mps { public:
    AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    MeshPoint * mp=MeshPointStack(s);
    assert(mp->T && mp->e > -1 && mp->d==2 ) ;
    long l=mp->Th->nTonEdge(mp->t,mp->e)==NBT;
    // cout << " nTonEdge " << l << endl;
    return SetAny<long>( l) ;}
    operator aType () const  { return atype<long>();}
};

class E_P_Stack_areaTriangle   :public  E_F0mps { public:
  AnyType operator()(Stack s)  const { throwassert(* (long *) s);
    MeshPoint * mp=MeshPointStack(s);
    assert(mp->T) ;
      double l=-1; // unset ...
    if(mp->d==2)	
      l= mp->T->area;	
    else if (mp->d==3 && mp->dHat==3 && mp->f >=0) {
	  R3 NN = mp->T3->N(mp->f);
	  l= NN.norme()/2.;
      }
     else if (mp->d==3 && mp->dHat==2 )
       l= mp->TS->mesure();
    else {
	  cout << "erreur : E_P_Stack_areaTriangle" << mp->d << " " << mp->f << endl;
	  ffassert(0); // undef
      }
    return SetAny<double>(l);}
    operator aType () const { return atype<double>();}

};

// Add jan 2018
class E_P_Stack_EdgeOrient  :public  E_F0mps { public:
    AnyType operator()(Stack s)  const { throwassert(* (long *) s);
        MeshPoint &mp= *MeshPointStack(s); // the struct to get x,y, normal , value
        double r=1;
        if(mp.d==2) {
            if( mp.T && mp.e >=0 )
            r=mp.T->EdgeOrientation(mp.e);
        }
        else if(mp.d==3 && mp.dHat==3) {
             if( mp.T3 && mp.f >=0 )
             r = mp.T3->faceOrient(mp.f);
        }
        else if(mp.d==3 && mp.dHat==2) {
            if( mp.TS && mp.e >=0 )
            r = mp.T3->EdgeOrientation(mp.e);
        }
        return r ;
    }
    operator aType () const  { return atype<double>();}
};


// add FH
class E_P_Stack_VolumeTet   :public  E_F0mps { public:
    AnyType operator()(Stack s)  const { throwassert(* (long *) s);
	MeshPoint * mp=MeshPointStack(s);
	assert(mp->T) ;
	double l=-1; // unset ...
	if (mp->d==3 && mp->dHat==3 && mp->T3 )
	  	l= mp->T3->mesure();
    else {
	    cout << "erreur : E_P_Stack_VolumeTet" << mp->d << " " << mp->f << endl;
	    ffassert(0); // undef 
	}
	return SetAny<double>(l);} 
    operator aType () const { return atype<double>();}

};

template<class R>
class  E_StopGC: public StopGC<R> {
public:
    typedef KN<R> Kn;
    typedef KN_<R> Kn_;

    Stack s;
    long n;
    long iter;
    Kn_ xx,gg;
    C_F0 citer,cxx,cgg;
    C_F0 stop;

    E_StopGC(Stack ss,long nn,const  Polymorphic * op): s(ss),n(nn),iter(-1),
    xx(0,0),gg(0,0),
    citer(CConstant<long*>(&iter)),
    cxx(dCPValue(&xx)),
    cgg(dCPValue(&gg)),
    stop(op,"(",citer,cxx,cgg)
    {

    }
    ~E_StopGC()
    {//  a verifier ???? FH....
        delete (E_F0 *) cxx; // ???
        delete (E_F0 *) cgg; // ???
        delete (E_F0 *) citer; // ???
        delete (E_F0 *) stop; // ???
    }
    bool Stop(int iterr, R *x, R * g)
    {
        iter=iterr;
        xx.set(x,n);
        gg.set(g,n);
        return GetAny<bool>(stop.eval(s));
    }
};


template<class R>
class LinearCG : public OneOperator
{ public:
  typedef KN<R> Kn;
  typedef KN_<R> Kn_;
  const int cas;

 class MatF_O: RNM_VirtualMatrix<R> { public:
   Stack stack;
   mutable  Kn x;
   C_F0 c_x;

   Expression  mat1,mat;
   typedef  typename RNM_VirtualMatrix<R>::plusAx plusAx;
   MatF_O(int n,Stack stk,const OneOperator * op)
     : RNM_VirtualMatrix<R>(n),stack(stk),
       x(n),c_x(CPValue(x)),
       mat1(op->code(basicAC_F0_wa(c_x))),
       mat( CastTo<Kn_>(C_F0(mat1,(aType)*op))) {
         }
   ~MatF_O() {
     if(mat1 != mat)
       delete mat;
      delete mat1;
      Expression zzz = c_x;
     delete zzz;

     }
   void addMatMul(const  Kn_  & xx, Kn_ & Ax) const {
      ffassert(xx.N()==Ax.N());
      x =xx;
      Ax  += GetAny<Kn_>((*mat)(stack));
      WhereStackOfPtr2Free(stack)->clean();
       }
    plusAx operator*(const Kn &  x) const {return plusAx(this,x);}
  virtual bool ChecknbLine(int n) const { return true;}
  virtual bool ChecknbColumn(int m) const { return true;}
    
};

  class E_LCG: public E_F0mps { public:
      
      
   const int cas;// <0 => Nolinear
   static const int n_name_param=6;

   static basicAC_F0::name_and_type name_param[] ;


  Expression nargs[n_name_param];
   
  const OneOperator *A, *C;
  Expression X,B;

      
  E_LCG(const basicAC_F0 & args,int cc) :cas(cc)
   {
      args.SetNameParam(n_name_param,name_param,nargs);
      {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(args[0].LeftValue());
         ffassert(op);
         A = op->Find("(",ArrayOfaType(atype<Kn* >(),false));
         ffassert(A);
      }
      if (nargs[2])
      {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[2]);
         ffassert(op);
         C = op->Find("(",ArrayOfaType(atype<Kn* >(),false));
         ffassert(C);
      }
       else  C =0;
      X = to<Kn*>(args[1]);
      if (args.size()>2)
        B = to<Kn*>(args[2]);
      else
        B=0;
   }
     
     virtual AnyType operator()(Stack stack)  const {
       int ret=-1;
       E_StopGC<R> *stop=0;
      try {
      Kn &x = *GetAny<Kn *>((*X)(stack));
      int n=x.N();
      MatF_O AA(n,stack,A);
      double eps = 1.0e-6;
	  double *veps=0;
      int nbitermax=  100;
      long verb = verbosity;

      if (nargs[0]) eps= GetAny<double>((*nargs[0])(stack));
      if (nargs[1]) nbitermax = GetAny<long>((*nargs[1])(stack));
      if (nargs[3]) veps=GetAny<double*>((*nargs[3])(stack));
      if (nargs[4]) verb=Abs(GetAny<long>((*nargs[4])(stack)));
      if (nargs[5]) stop= new E_StopGC<R>(stack,n,dynamic_cast<const  Polymorphic *>(nargs[5]));
      long gcverb=51L-Min(Abs(verb),50L);
      if(verb==0) gcverb = 1000000000;// no print
      if(veps) eps= *veps;
      KN<R>  bzero(B?1:n); // const array zero
      bzero=R();
      KN<R> *bb=&bzero;
      if (B) {
        Kn &b = *GetAny<Kn *>((*B)(stack));
        R p = (b,b);
       if (p== R())
         {
          // ExecError("Sorry LinearCG work only with nul right hand side, so put the right hand in the function");
          }
         bb = &b;
      }
      if (cas<0) {
       if (C)
         { MatF_O CC(n,stack,C);
           ret = NLCG(AA,CC,x,nbitermax,eps, gcverb,stop );}
        else
           ret = NLCG(AA,MatriceIdentite<R>(n),x,nbitermax,eps, gcverb,stop);
        }
      else
      if (C)
       { MatF_O CC(n,stack,C);
         ret = ConjuguedGradient2(AA,CC,x,*bb,nbitermax,eps, gcverb, stop );}
      else
         ret = ConjuguedGradient2(AA,MatriceIdentite<R>(n),x,*bb,nbitermax,eps, gcverb, stop );
      if(veps) *veps = -(eps);
      }
      catch(...)
      {
       if( stop) delete stop;
        throw;
      }
     if( stop) delete stop;

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


template<class R>
basicAC_F0::name_and_type  LinearCG<R>::E_LCG::name_param[]= {
  {   "eps", &typeid(double)  },
  {   "nbiter",&typeid(long) },
  {   "precon",&typeid(Polymorphic*)},
  {   "veps" ,  &typeid(double*) },
  {   "verbosity" ,  &typeid(long)},
  {   "stop" ,  &typeid(Polymorphic*)}
};


template<class R>
class LinearGMRES : public OneOperator
{ public:
  typedef KN<R> Kn;
  typedef KN_<R> Kn_;
  const int cas;

 class MatF_O: RNM_VirtualMatrix<R> { public:
   Stack stack;
   mutable  Kn x;
   C_F0 c_x;
   Kn *b;
   Expression  mat1,mat;
   typedef  typename RNM_VirtualMatrix<R>::plusAx plusAx;
   MatF_O(int n,Stack stk,const OneOperator * op,Kn *bb)
     : RNM_VirtualMatrix<R>(n),
       stack(stk),
       x(n),c_x(CPValue(x)),b(bb),
       mat1(op->code(basicAC_F0_wa(c_x))),
       mat( CastTo<Kn_>(C_F0(mat1,(aType)*op))  /*op->code(basicAC_F0_wa(c_x))*/) {
       }
     ~MatF_O() { if(mat1!=mat) delete mat; delete mat1; delete c_x.LeftValue();}
   void addMatMul(const  Kn_  & xx, Kn_ & Ax) const {

     ffassert(xx.N()==Ax.N());
       x =xx;
      Ax    += GetAny<Kn_>((*mat)(stack));
      if(b && &Ax!=b) Ax += *b; // Ax -b => add b (not in cas of init. b c.a.d  &Ax == b
      WhereStackOfPtr2Free(stack)->clean(); //  add dec 2008
   }
    plusAx operator*(const Kn &  x) const {return plusAx(this,x);}
  virtual bool ChecknbLine(int n) const { return true;}
  virtual bool ChecknbColumn(int m) const { return true;}

};


  class E_LGMRES: public E_F0mps { public:
   const int cas;// <0 => Nolinear
   static basicAC_F0::name_and_type name_param[] ;
   static const int n_name_param =7;
   Expression nargs[n_name_param];
  const OneOperator *A, *C;
  Expression X,B;

  E_LGMRES(const basicAC_F0 & args,int cc) :cas(cc)
   {
      args.SetNameParam(n_name_param,name_param,nargs);
      {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(args[0].LeftValue());
         ffassert(op);
         A = op->Find("(",ArrayOfaType(atype<Kn* >(),false));
          ffassert(A);
      }
      if (nargs[2])
      {  const  Polymorphic * op=  dynamic_cast<const  Polymorphic *>(nargs[2]);
         ffassert(op);
         C = op->Find("(",ArrayOfaType(atype<Kn* >(),false));
          ffassert(C);
      }
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
       E_StopGC<R> *stop=0;
      if (B)   b = *GetAny<Kn *>((*B)(stack));
      else     b= R();
      int n=x.N();
      int dKrylov=50;
      double eps = 1.0e-6;
      int nbitermax=  100;
      long verb = verbosity;
      if (nargs[0]) eps= GetAny<double>((*nargs[0])(stack));
      if (nargs[1]) nbitermax = GetAny<long>((*nargs[1])(stack));
      if (nargs[3]) eps= *GetAny<double*>((*nargs[3])(stack));
      if (nargs[4]) dKrylov= GetAny<long>((*nargs[4])(stack));
      if (nargs[5]) verb=Abs(GetAny<long>((*nargs[5])(stack)));
      if (nargs[6]) stop= new E_StopGC<R>(stack,n,dynamic_cast<const  Polymorphic *>(nargs[6]));

	 long gcverb=51L-Min(Abs(verb),50L);


      int ret=-1;
      if(verbosity>4)
        cout << "  ..GMRES: eps= " << eps << " max iter " << nbitermax
             << " dim of Krylov space " << dKrylov << endl;
        KNM<R> H(dKrylov+1,dKrylov+1);
	int k=dKrylov;//,nn=n;
       double epsr=eps;

         KN<R>  bzero(B?1:n); // const array zero
         bzero=R();
         KN<R> *bb=&bzero;
         if (B) {
             Kn &b = *GetAny<Kn *>((*B)(stack));
             R p = (b,b);
             if (p)
             {
                 // ExecError("Sorry MPILinearCG work only with nul right hand side, so put the right hand in the function");
             }
             bb = &b;
         }
         KN<R> * bbgmres =0;
         if ( !B ) bbgmres=bb; // none zero if gmres without B
         MatF_O AA(n,stack,A,bbgmres);
         if(bbgmres ){
             AA.addMatMul(*bbgmres,*bbgmres); // Ok Ax == b -> not translation of b .
             *bbgmres = - *bbgmres;
             if(verbosity>1) cout << "  ** GMRES set b =  -A(0);  : max=" << bbgmres->max() << " " << bbgmres->min()<<endl;
         }

      if (cas<0) {
        ErrorExec("NL GMRES:  to do! sorry ",1);
        }
      else
       {
       if (C)
        { MatF_O CC(n,stack,C,0);
         ret=GMRES(AA,(KN<R> &)x, *bb,CC,H,k,nbitermax,epsr,verb,stop);}
       else
         ret=GMRES(AA,(KN<R> &)x, *bb,MatriceIdentite<R>(n),H,k,nbitermax,epsr,verb,stop);
       }
      if(verbosity>99)    cout << " Sol GMRES :" << x << endl;
         if(stop) delete stop;
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


template<class R>
basicAC_F0::name_and_type  LinearGMRES<R>::E_LGMRES::name_param[]= {
  {   "eps", &typeid(double)  },
  {   "nbiter",&typeid(long) },
  {   "precon",&typeid(Polymorphic*)},
  {   "veps" ,  &typeid(double*) },
  {   "dimKrylov", &typeid(long) },
  {   "verbosity", &typeid(long) },
  {   "stop" ,  &typeid(Polymorphic*)}
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
 const T & operator[](int i) const {return v[i];}
};
template<class T,int N>
ostream & operator<<(ostream & f,const Smallvect<T,N> & v)
{
    for(int i=0;i<N;++i)  f << v[i] << ' ';
    return f;
}

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
		 ffassert(kkk++<10);
	       assert(nbdfv <= j);
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

bool  InCircularList(const int *p,int i,int k)
//  find k in circular list:  i , p[i], p[p[i]], ...
{
    int j=i,l=0;
    do {
	if (j==k) return true;
	ffassert(l++<10);
        j=p[j];
    } while (j!=i);
    return false;
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

       ffassert(ndfv.N()==Th.nv);
       ffassert(ndfe.N()==Th.neb);

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
		 if(verbosity>5)
			cout << "lab1:  e[" << pke1[i1] << "]  v0:   " <<  e[0].x << " " << e[0].y << "  s = " << x0
			<< "\t v1 " <<  e[1].x << " " << e[1].y << "  s = " << x1 << endl;
		    xmn=Min(x1,x0,xmn);
                 xmx=Max(x1,x0,xmx);
                 hmn=Min(hmn,Abs(x1-x0));
               }
             }
            ffassert(hmn>1.0e-20);
            double coef = 8/hmn;
            double x0 = xmn;
            if (verbosity > 2)
            cout << "  --Update: periodic " << xmn << " " << xmx << " " << " h=" << hmn << endl;
            ffassert(!n1 || (coef>1e-10 && (xmx-xmn)*coef < 1.e7 ));

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
                     im=m.insert(pair<int,int2>(i0,i2)).first;
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
                    cout << i2 << " : " <<Th(e[0]) << " " << Th(e[1]) << ":: ";
                 mp->set(e[0].x,e[0].y);
                 double xx0=GetAny<double>((*periodic[k+kk2])(stack));
                 mp->set(e[1].x,e[1].y);
                 double xx1=GetAny<double>((*periodic[k+kk2])(stack));
		 if(verbosity>5 )
		      cout << "lab2:  e[" << pke2[i2] << "]  v0:   " <<  e[0].x << " " << e[0].y << "  s = " << xx0
		      << "\t v1 " <<  e[1].x << " " << e[1].y << "  s = " << xx1 << endl;

                 int i0= int((xx0-x0)*coef);
                 int i1= int((xx1-x0)*coef);
                 map<int,int2>::iterator im0=closeto(m,i0);
                 map<int,int2>::iterator im1=closeto(m,i1);
                 if(im0 == m.end() || im1 == m.end() )
                 	{

			  cout << "Abscisse: s0 = "<< xx0 << " <==> s1 " << xx1  <<endl;
			  ExecError("periodic: Sorry one vertex of edge is losted "); }
                  int ie1=-1;
                 if      (((ie1=im0->second[0])==im1->second[1]) && (ie1>=0)) ;
                 else if (((ie1=im0->second[0])==im1->second[1]) && (ie1>=0)) ;
                 else if (((ie1=im0->second[1])==im1->second[1]) && (ie1>=0)) ;
                 else if (((ie1=im0->second[1])==im1->second[0]) && (ie1>=0)) ;
                 else if (((ie1=im0->second[0])==im1->second[0]) && (ie1>=0)) ;
                 else
                  {
                   cout << ie2 << " ~ " << im0->second[0] << " " << im0->second[1] << ", "
                                        << im1->second[0] << " " << im1->second[1] << endl;
                   ExecError("periodic: Sorry one egde is losted "); }
		 if(verbosity>50)
		   cout << " ( " << im0->second << " , " << im1->second << " ) .. ";
		 ffassert(ie1>=0 && ie1 < Th.neb );
                 const BoundaryEdge & ep =Th.bedges[ie1];
                 mp->set(ep[0].x,ep[0].y);
                 double yy0=GetAny<double>((*periodic[k+kk1])(stack));
                 mp->set(ep[1].x,ep[1].y);
                 double yy1=GetAny<double>((*periodic[k+kk1])(stack));
		if(verbosity>50)
			cout << " e0: s  "<< xx0 << " " << xx1 << "e1 s "<< yy0 << " " << yy1  ;

                 pke1[i2]=ie1*2+ ( ( (yy1-yy0) < 0)  == ( (xx1-xx0) < 0) ) ;

                 if (verbosity >50)
                 cout << " \t  edge " << ie1 << " <=> " << ie2 << " "
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
	   if(!InCircularList(ndfe,ie1,ie2))   // merge of two list
	      Exchange(ndfe[ie1],ndfe[ie2]);
	   for (int ke2=0;ke2<2;ke2++)
             {
                int ke1=ke2;
                if(!sens) ke1=1-ke1;
                int iv1=Th(Th.bedges[ie1][ke1]);
                int iv2=Th(Th.bedges[ie2][ke2]);
                if (!InCircularList(ndfv,iv1,iv2)) {  // merge of two list
		   Exchange(ndfv[iv2],ndfv[iv1]);
                   if (verbosity >50)
		   {
		     cout << "  vertex " << iv1 <<  "<==> " << iv2 << " list : " << iv1;
		     int i=iv1,k=0;
		     while ( (i=ndfv[i]) != iv1 && k++<10)
			 cout << ", "<< i ;
		       cout << endl;
		   }}
             }

         }
         // generation de numero de dlt

          nbdfv = numeroteclink(ndfv) ;
          nbdfe = numeroteclink(ndfe) ;
          if (verbosity>2)
            cout << "  -- nb df on vertices " << nbdfv << endl;
        delete [] link1;
        delete [] link2;
        return true; //new FESpace(**ppTh,*tef,nbdfv,ndfv,nbdfe,ndfe);
      }
        else {
	      delete [] link1;
	      delete [] link2;
        }
   }
   return false;
}


bool  v_fes::buildperiodic(Stack stack,int & nbdfv, KN<int> & ndfv,int & nbdfe, KN<int> & ndfe) {
  return BuildPeriodic(nbcperiodic,periodic,**ppTh,stack,nbdfv,ndfv,nbdfe,ndfe);

}
#ifdef ZZZZZZZZ
FESpace * pfes_tef::update() {
   typedef Smallvect<int,2> int2;

    if (nbcperiodic ) {
       const Mesh &Th(**ppTh);
       KN<int> ndfv(Th.nv);
       KN<int> ndfe(Th.neb);
       int nbdfv,nbdfe;
        return  new FESpace(**ppTh,*tef,nbdfv,ndfv,nbdfe,ndfe);
      }
     else
          return  new FESpace(**ppTh,*tef);
}
#endif

struct OpMake_pfes_np {
  static const int n_name_param =1;
  static basicAC_F0::name_and_type name_param[] ;
};

basicAC_F0::name_and_type  OpMake_pfes_np::name_param[]= {
  "periodic", &typeid(E_Array)
};

// by default, in DSL a FE is 2D, in the building of fespace, if mesh3-S used them associate the corresponding FE
// mapping between TypeOfFE2 and TypeOfFES
map<TypeOfFE *,TypeOfFE3 *> TEF2dto3d;
AnyType TypeOfFE3to2(Stack,const AnyType &b) {
    TypeOfFE3 *t3=0;
    TypeOfFE  *t2=GetAny<TypeOfFE *>(b);
    map<TypeOfFE *,TypeOfFE3 *>::const_iterator i=TEF2dto3d.find(t2);
    if(i != TEF2dto3d.end())
	t3=i->second;

    if(t3==0)
      {
	  cerr << " sorry no cast to this 3d finite element " <<endl;
	  ExecError( " sorry no cast to this 3d finite element ");
      }
    return t3;
}

// mapping between TypeOfFE2 and TypeOfFES
map<TypeOfFE *,TypeOfFES *> TEF2dtoS;
AnyType TypeOfFESto2(Stack,const AnyType &b) {
    TypeOfFES *tS=0;
    TypeOfFE  *t2=GetAny<TypeOfFE *>(b);
    map<TypeOfFE *,TypeOfFES *>::const_iterator i=TEF2dtoS.find(t2);
    if(i != TEF2dtoS.end())
        tS=i->second;

    if(tS==0)
    {
        cerr << " sorry no cast to this surface finite element " <<endl;
        ExecError( " sorry no cast to this surface finite element ");
    }
    return tS;
}

// mapping between TypeOfFE2 and TypeOfFEL
map<TypeOfFE *,TypeOfFEL *> TEF2dtoL;
AnyType TypeOfFELto2(Stack,const AnyType &b) {
    TypeOfFEL *tL=0;
    TypeOfFE  *t2=GetAny<TypeOfFE *>(b);
    map<TypeOfFE *,TypeOfFEL *>::const_iterator i=TEF2dtoL.find(t2);
    if(i != TEF2dtoL.end())
        tL=i->second;
    
    if(tL==0)
    {
        cerr << " sorry no cast to this surface finite element " <<endl;
        ExecError( " sorry no cast to this surface finite element ");
    }
    return tL;
}

TypeOfFE * FindFE2(const char * s)
{
    for (ListOfTFE * i=ListOfTFE::all;i;i=i->next)
	if(strcmp(i->name,s)==0)
	    return i->tfe;
    cout << " s =" << s << endl;
    lgerror("FindFE2 ");
    return 0;
}


typedef TypeOfFE TypeOfFE2;
template<class pfes,class Mesh,class TypeOfFE,class pfes_tefk>
struct OpMake_pfes: public OneOperator , public OpMake_pfes_np {

  struct Op: public E_F0mps {
  public:

    Expression eppTh;
    Expression eppfes;
    const E_Array & atef;
    int nb;
    int nbcperiodic;
    Expression *periodic;
    KN<int>  tedim;
    Op(Expression ppfes,Expression ppTh, const E_Array & aatef,int nbp,Expression * pr,KN<int> &ttedim)
      : eppTh(ppTh),eppfes(ppfes),atef(aatef),nbcperiodic(nbp),periodic(pr),tedim(ttedim) {
    }
    ~Op() { if(periodic) delete []periodic;}
    AnyType operator()(Stack s)  const {
      const int d = Mesh::Rd::d;
      const Mesh ** ppTh = GetAny<const Mesh  **>( (*eppTh)(s) );
      AnyType r = (*eppfes)(s) ;
      const TypeOfFE ** tef= new  const TypeOfFE * [ atef.size()];
      for (int i=0;i<atef.size();i++)
	if(tedim[i]==d)
	  tef[i]= GetAny<TypeOfFE *>(atef[i].eval(s));
	else if(tedim[i]==2 && d ==3)
	  tef[i]= GetAny<TypeOfFE *>(TypeOfFE3to2(s,atef[i].eval(s)));
	else ffassert(0);

      pfes * ppfes = GetAny<pfes *>(r);
      bool same = true;
      for (int i=1;i<atef.size();i++)
      same &= atef[i].LeftValue() == atef[1].LeftValue();
      *ppfes = new pfes_tefk(ppTh,tef,atef.size(),s,nbcperiodic,periodic);
      return r;}
  } ;

  E_F0 * code(const basicAC_F0 & args)  const
  {
    int nbcperiodic=0;
    Expression *periodic=0;
    Expression nargs[n_name_param];

    args.SetNameParam(n_name_param,name_param,nargs);
      GetPeriodic(Mesh::Rd::d,nargs[0],nbcperiodic,periodic);
    aType t_tfe= atype<TypeOfFE*>();
    aType t_tfe2= atype<TypeOfFE2*>();
    int d=  TypeOfFE::Rd::d;
    string sdim= d ?  " 2d : " : " 3d : " ;
    const E_Array * a2(dynamic_cast<const E_Array *>(args[2].LeftValue()));
    ffassert(a2);
    int N = a2->size(); ;
    if (!N) CompileError(sdim+" We wait an array of Type of Element ");
      KN<int> tedim(N);
    for (int i=0;i< N; i++)
      if ((*a2)[i].left() == t_tfe)
	  tedim[i]=d;
      else if ((*a2)[i].left() ==t_tfe2)
	  tedim[i]=2;
      else
	    CompileError(sdim+" We wait an array of  Type of Element ");
    //    ffassert(0);
    return  new Op(args[0],args[1],*a2,nbcperiodic,periodic,tedim);
  }
  OpMake_pfes() :
    OneOperator(atype<pfes*>(),atype<pfes*>(),atype<const Mesh **>(),atype<E_Array>()) {}
};

inline pfes* MakePtr2(pfes * const &p,pmesh * const &  a, TypeOfFE * const & tef)
{ *p=new pfes_tef(a,tef) ;
  return p;}

inline pfes3* MakePtr3(pfes3 * const &p,pmesh3 * const &  a, TypeOfFE3 * const & tef)
{ *p=new pfes3_tef(a,tef) ;
  return p;}

inline pfesS* MakePtrS(pfesS * const &p,pmeshS * const &  a, TypeOfFES * const & tef)
{ *p=new pfesS_tef(a,tef) ;
    return p;}

inline pfesL* MakePtrL(pfesL * const &p,pmeshL * const &  a, TypeOfFEL * const & tef)
{ *p=new pfesL_tef(a,tef) ;
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
	*p=new pfes_tef(th,tef,s,nbcperiodic,periodic) ;
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


class OP_MakePtr3 { public:
    class Op : public E_F0mps  { public:
	static const int n_name_param =1;
      static basicAC_F0::name_and_type name_param[] ;
      Expression nargs[n_name_param];
      typedef pfes3 * R;
      typedef pfes3 * A;
      typedef pmesh3 * B;
      typedef TypeOfFE3 * C;
      Expression a,b,c;
      int nbcperiodic ;
      Expression *periodic;
      Op(const basicAC_F0 & args);

      AnyType operator()(Stack s) const  {
	A p= GetAny<A>( (*a)(s) );
	B th= GetAny<B>( (*b)(s) );
	C tef= GetAny<C>( (*c)(s) );
	*p=new pfes3_tef(th,tef,s,nbcperiodic,periodic) ;
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

class OP_MakePtrS { public:
    class Op : public E_F0mps  { public:
        //  static int GetPeriodic(Expression  bb, Expression & b,Expression & f);
        static const int n_name_param =1;
        static basicAC_F0::name_and_type name_param[] ;
        Expression nargs[n_name_param];
        typedef pfesS * R;
        typedef pfesS * A;
        typedef pmeshS * B;
        typedef TypeOfFES * C;
        Expression a,b,c;
        int nbcperiodic ;
        Expression *periodic;
        Op(const basicAC_F0 & args);

        AnyType operator()(Stack s) const  {
            A p= GetAny<A>( (*a)(s) );
            B th= GetAny<B>( (*b)(s) );
            C tef= GetAny<C>( (*c)(s) );
            *p=new pfesS_tef(th,tef,s,nbcperiodic,periodic) ;
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

class OP_MakePtrL { public:
    class Op : public E_F0mps  { public:
        //  static int GetPeriodic(Expression  bb, Expression & b,Expression & f);
        static const int n_name_param =1;
        static basicAC_F0::name_and_type name_param[] ;
        Expression nargs[n_name_param];
        typedef pfesL * R;
        typedef pfesL * A;
        typedef pmeshL * B;
        typedef TypeOfFEL * C;
        Expression a,b,c;
        int nbcperiodic ;
        Expression *periodic;
        Op(const basicAC_F0 & args);
        
        AnyType operator()(Stack s) const  {
            A p= GetAny<A>( (*a)(s) );
            B th= GetAny<B>( (*b)(s) );
            C tef= GetAny<C>( (*c)(s) );
            *p=new pfesL_tef(th,tef,s,nbcperiodic,periodic) ;
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


void GetPeriodic(const int d,Expression perio,    int & nbcperiodic ,    Expression * &periodic)
{
    ffassert(d==2 || d ==3);
      if ( perio)
       {
         if( verbosity>1)
         cout << "  -- Periodical Condition to do" << endl;
         const E_Array * a= dynamic_cast<const  E_Array *>(perio);
         ffassert(a);
         int n = a->size();
        nbcperiodic= n/2;
        if( verbosity>1)
        cout << "    the number of periodicBC " << n << endl;
        if ( 2*nbcperiodic != n ) CompileError(" Sorry the number of periodicBC must by even");
        periodic = new Expression[n*d];
        for (int i=0,j=0;i<n;i++,j+=d)
	  if(d==2)
	    { if (GetPeriodic((*a)[i],periodic[j],periodic[j+1])==0)
            CompileError(" a sub array of periodic BC must be [label, realfunction ]");
	    }
	  else if (d==3)
	    { if (GetPeriodic((*a)[i],periodic[j],periodic[j+1],periodic[j+2])==0)
		CompileError(" a sub array of periodic BC must be [label, realfunction , realfunction]");
	    }
	  else ffassert(0);
        }
}


OP_MakePtr2::Op::Op(const basicAC_F0 & args)
  : a(to<A>(args[0])),b(to<B>(args[1])),c(to<C>(args[2]))
     {
      nbcperiodic=0;
      periodic=0;
      args.SetNameParam(n_name_param,name_param,nargs);
      GetPeriodic(2,nargs[0],nbcperiodic,periodic);
     }
//3D volume
OP_MakePtr3::Op::Op(const basicAC_F0 & args)
  : a(to<A>(args[0])),b(to<B>(args[1])),c(to<C>(args[2]))
{
  nbcperiodic=0;
  periodic=0;
  args.SetNameParam(n_name_param,name_param,nargs);
  GetPeriodic(3,nargs[0],nbcperiodic,periodic);
}
//3D surface
OP_MakePtrS::Op::Op(const basicAC_F0 & args)
: a(to<A>(args[0])),b(to<B>(args[1])),c(to<C>(args[2]))
{
    nbcperiodic=0;
    periodic=0;
    args.SetNameParam(n_name_param,name_param,nargs);
    GetPeriodic(3,nargs[0],nbcperiodic,periodic);
}
//3D curve
OP_MakePtrL::Op::Op(const basicAC_F0 & args)
: a(to<A>(args[0])),b(to<B>(args[1])),c(to<C>(args[2]))
{
    nbcperiodic=0;
    periodic=0;
    args.SetNameParam(n_name_param,name_param,nargs);
    GetPeriodic(3,nargs[0],nbcperiodic,periodic);    
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
int GetPeriodic(Expression  bb, Expression & b,Expression & f1,Expression & f2)
{
    const E_Array * a= dynamic_cast<const E_Array *>(bb);
    if(a && a->size() == 3)
      {
	  b= to<long>((*a)[0]);
	  f1= to<double>((*a)[1]);
	  f2= to<double>((*a)[2]);
	  return 1;
      }
    else
	return 0;
}

basicAC_F0::name_and_type  OP_MakePtr2::Op::name_param[]= {
  "periodic", &typeid(E_Array)
};

basicAC_F0::name_and_type  OP_MakePtr3::Op::name_param[]= {
  "periodic", &typeid(E_Array)
};

basicAC_F0::name_and_type  OP_MakePtrS::Op::name_param[]= {
    "periodic", &typeid(E_Array)
};

basicAC_F0::name_and_type  OP_MakePtrL::Op::name_param[]= {
    "periodic", &typeid(E_Array)
};

inline pfes* MakePtr2(pfes * const &p,pmesh * const &  a){
      *p=new pfes_tef(a,&P1Lagrange);
       return p ;}

inline pfes* MakePtr2(pfes * const &p,pfes * const &  a,long const & n){
       *p= new pfes_fes(a,n);
       return p ;}

 long FindTxy(Stack s,pmesh * const &  ppTh,const double & x,const double & y)
{
   R2 P(x,y),PHat;
   bool outside;
     MeshPoint & mp = *MeshPointStack(s);
   const Mesh * pTh= *ppTh;
   if(pTh == 0) return 0;
   const Triangle * K=pTh->Find(mp.P.p2(),PHat,outside);
   if (!outside)
     mp.set(*pTh,P,PHat,*K,K->lab);
   else return 0;
   return 1;
}


template<class K>
KN<K> * pfer2vect( pair<FEbase<K,v_fes> *,int> p)
 {
    KN<K> * x=p.first->x();
    if ( !x) {  // defined
      FESpace * Vh= p.first->newVh(); // cout << "test 1" << endl;
      throwassert( Vh);
      *p.first = x = new KN<K>(Vh->NbOfDF);
      *x=K();
    }
    return x;}

template<class K>
pmesh pfer_Th(pair<FEbase<K,v_fes> *,int> p)
{
    if (!p.first->Vh) p.first->Vh= p.first->newVh();
    throwassert( !!p.first->Vh);
    return &p.first->Vh->Th;
}
//  add to refesh Th in fespace oct 2018
template<class K,class v_fes>
bool pfer_refresh(pair<FEbase<K,v_fes> *,int> p)
{
    int n = 0;
    if( p.first->x() ) n = p.first->x()->N();
    const FESpace * Vho = p.first->Vh;
    p.first->Vh= p.first->newVh();
    if ( n && n !=  p.first->Vh->NbOfDF)
        * p.first  = new KN<K>(p.first->Vh->NbOfDF,K());
    return Vho != (const FESpace * ) p.first->Vh;// FE change !!!
}
template<class K>
long pfer_nbdf(pair<FEbase<K,v_fes> *,int> p)
 {
   if (!p.first->Vh) p.first->Vh= p.first->newVh();
   throwassert( !!p.first->Vh);
   return p.first->Vh->NbOfDF;
 }

double pmesh_area(pmesh * p)
 { throwassert(p) ;  return *p ? (**p).area : 0.0;}
double pmesh_bordermeasure(pmesh * p)
{ throwassert(p) ;  return *p ? (**p).lenbord : 0.0;}

long pmesh_nt(pmesh * p)
 { throwassert(p) ;  return *p ? (**p).nt : 0;}
long pmesh_nbe(pmesh * p)
{ throwassert(p) ;  return *p ? (**p).neb : 0;}
long pmesh_nv(pmesh * p)
 { throwassert(p) ;  return *p ? (**p).nv : 0;}

double pmesh_hmax(pmesh * p)
{ throwassert(p && *p) ;
    double hmax2 =0;
    const Mesh & Th = **p;
    for(int k=0; k< Th.nt; ++k)
        for(int e=0; e<3; ++e)
            hmax2=max(hmax2,Th[k].lenEdge2(e));
    return sqrt(hmax2);}

double pmesh_hmin(pmesh * p)
{ throwassert(p && *p) ;
    double hmin2 =1e100;
    const Mesh & Th = **p;
    for(int k=0; k< Th.nt; ++k)
        for(int e=0; e<3; ++e)
            hmin2=min(hmin2,Th[k].lenEdge2(e));
    return sqrt(hmin2);}

long pVh_ndof(pfes * p)
 { throwassert(p && *p);
   FESpace *fes=**p; ;  return fes->NbOfDF ;}
pmesh pVh_Th(pfes * p)
{ throwassert(p && *p);
    FESpace *fes=**p; ;  return &fes->Th ;}

long pVh_nt(pfes * p)
 { throwassert(p && *p);
   FESpace *fes=**p; ;  return fes->NbOfElements ;}
long pVh_ndofK(pfes * p)
 { throwassert(p && *p);
   FESpace *fes=**p;   return (*fes)[0].NbDoF() ;}

long mp_nuTriangle(MeshPoint * p)
 { throwassert(p  && p->Th && p->T);
   long nu=0;
   if(p->d==2)
     nu=(*p->Th)(p->T);
   else if  (p->d==3 && p->dHat==3)
     nu=(*p->Th3)(p->T3);
   else if  (p->d==3 && p->dHat==2)
       nu=(*p->ThS)(p->TS);
   else if  (p->d==3 && p->dHat==1)
       nu=(*p->ThL)(p->TL);
   else ffassert(0);
   delete p;
   return nu ;}

long mp_region(MeshPoint * p)
 {
   long  nu(p->region);
   delete p;
   return nu ;}


class pVh_ndf : public ternary_function<pfes *,long,long,long>
{ public:

  class Op : public E_F0mps { public:
      Expression a,b,c;
       Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {}
       AnyType operator()(Stack s)  const
        {
           pfes * p(GetAny<pfes *>((*a)(s)));
           long  k(GetAny<long>((*b)(s)));
           long  i(GetAny<long>((*c)(s)));
           throwassert(p && *p);
           FESpace *fes=**p;
           throwassert(fes && k >=0 && k < fes->NbOfElements );
           FElement K=(*fes)[k];
           throwassert(i>=0 && i <K.NbDoF() );
           long ret(K(i));
           return  ret;
         }
  };
};

class Op_CopyArray : public OneOperator { public:
    Op_CopyArray():OneOperator(atype<void>(),atype<E_Array>(),atype<E_Array>()) {}
    E_F0 * code(const basicAC_F0 & args) const ;
};

template<class R,int dd>
AnyType pfer2R(Stack s,const AnyType &a)
{
  pair< FEbase<R,v_fes> *  ,int> ppfe=GetAny<pair< FEbase<R,v_fes> *,int> >(a);
  FEbase<R,v_fes> & fe( *ppfe.first);
  int componante=ppfe.second;
  if ( !fe.x()) {
   if ( !fe.x()){
     return   SetAny<R>(0.0);
    }
  }

  const FESpace & Vh(*fe.Vh);
  const Mesh & Th(Vh.Th);
  assert(Th.ntet==0 && Th.volume==0 && Th.triangles != 0);
  MeshPoint & mp = *MeshPointStack(s);
  const Triangle *K;
  R2 PHat;
  bool outside=false;
  bool qnu=true;
  if ( mp.Th == &Th && mp.T)
   {
    qnu=false;
    K=mp.T;
    PHat=mp.PHat.p2();
   }
  else if ( mp.other.Th == & Th && mp.other.P.x == mp.P.x && mp.other.P.y == mp.P.y )
   {
    K=mp.other.T;
    PHat=mp.other.PHat.p2();
    outside = mp.other.outside;
   }
  else {
    if (mp.isUnset()) ExecError("Try to get unset x,y, ...");
    K=Th.Find(mp.P.p2(),PHat,outside);
    mp.other.set(Th,mp.P.p2(),PHat,*K,0,outside);
    }
  const FElement KK(Vh[Th(K)]);
  if (outside && !KK.tfe->NbDfOnVertex && !KK.tfe->NbDfOnEdge)
    return   SetAny<R>(0.0);

 const R rr = KK(PHat,*fe.x(),componante,dd);
  return SetAny<R>(rr);
}


template<class R>
AnyType set_fe (Stack s,Expression ppfe, Expression e)
  {
   long kkff = Mesh::kfind,  kkth = Mesh::kthrough;
      StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);

    MeshPoint *mps=MeshPointStack(s),mp=*mps;
    pair<FEbase<R,v_fes> *,int>  pp=GetAny<pair<FEbase<R,v_fes> *,int> >((*ppfe)(s));
    FEbase<R,v_fes> & fe(*pp.first);
    const  FESpace * pVh(fe.newVh());

    if(!pVh ) ExecError("Unset FEspace (Null mesh ? ) on  uh= ");
    const  FESpace & Vh= *pVh;
    KN<R> gg(Vh.MaximalNbOfDF());
    const  Mesh & Th(Vh.Th);
    TabFuncArg tabexp(s,Vh.N);
    tabexp[0]=e;

    if(Vh.N!=1)
    {  cerr << " Try to set a  vectorial  FE function  (nb  componant=" <<  Vh.N << ") with one scalar " << endl;
       ExecError(" Error interploation (set)  FE function (vectorial) with a scalar");
    }
    KN<R> * y=new  KN<R>(Vh.NbOfDF);
    KN<R> & yy(*y);
    KN<R> Viso(100);
    for (int i=0;i<Viso.N();i++)
      Viso[i]=0.01*i;

    FElement::aIPJ ipj(Vh[0].Pi_h_ipj());
    FElement::aR2  PtHat(Vh[0].Pi_h_R2());
    KN<double>   Aipj(ipj.N());
    KN<R>   Vp(PtHat.N());
    const E_F0 & ff(* (const  E_F0 *) e ) ;

   if (Vh.isFEMesh() )
    {
      ffassert(Vh.NbOfDF == Th.nv && Vh.N == 1 );
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
         yy[iv] = GetAny<R>( ff(s) );
          sptr->clean(); // modif FH mars 2006  clean Ptr
       }
    }
   else

     for (int t=0;t<Th.nt;t++)
      {
         FElement K(Vh[t]);
         int nbdf=K.NbDoF();
         gg=R();
#ifdef OLDPih
// old method
         K.Pi_h(gg,F_Pi_h,F,&tabexp);
#else
         K.Pi_h(Aipj);
         for (int p=0;p<PtHat.N();p++)
          {
            mps->set(K.T(PtHat[p]),PtHat[p],K);
            Vp[p]=GetAny<R>( ff(s) );
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
           sptr->clean(); // modif FH mars 2006  clean Ptr

      }
    *mps=mp;
    fe=y;
     kkff = Mesh::kfind - kkff;
     kkth = Mesh::kthrough -kkth;

    if(verbosity>1)
      ShowBound(*y,cout)
           << " " << kkth << "/" << kkff << " =  " << double(kkth)/Max<double>(1.,kkff) << endl;
    return SetAny<FEbase<R,v_fes>*>(&fe);
  }
AnyType set_feoX_1 (Stack s,Expression ppfeX_1, Expression e)
  { // inutile
    // meme chose que  v(X1,X2);
      StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);
    typedef const interpolate_f_X_1<R>::CODE * code;
    MeshPoint mp= *MeshPointStack(s);
    code ipp = dynamic_cast<code>(ppfeX_1);

    pair<FEbase<R,v_fes> *,int>  pp=GetAny<pair<FEbase<R,v_fes> *,int> >((*ipp->f)(s));
    FEbase<R,v_fes> & fe(*pp.first);
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

         gg=R();

         K.Pi_h(gg,FoX_1_Pi_h,F,&tabexp);
         for (int df=0;df<nbdf;df++)
          (*y)[K(df)] =  gg[df] ;
          sptr->clean(); // modif FH mars 2006  clean Ptr
      }
    *MeshPointStack(s)=mp;
    fe=y;
    if(verbosity>1)
    cout << "  -- interpole f= g*X^-1, function's bound:  " << y->min() << " " << y->max() << endl;
    return SetAny<FEbase<R,v_fes>*>(&fe);
  }

template<class K>
E_set_fev<K>::E_set_fev(const E_Array * a,Expression pp,int ddim)
  :dim(ddim), aa(*a),ppfe(pp),optimize(true),
   where_in_stack_opt(),optiexp0(),optiexpK()
{
  aa.map(to<K>) ;
  bool kdump=false;
  if(optimize)
    { // new code Optimized  -------
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


template<class K>
AnyType E_set_fev<K>::operator()(Stack s)  const
{
  if(dim== 2)  return  Op2d(s);
  else if(dim == 3)  return Op3d(s);
  else if(dim == 4)  return OpS(s);
  return Nothing;
}


template<class K>
AnyType E_set_fev<K>::Op3d(Stack s)  const
{
  //  voir E_set_fev3  ( pb de consitance a revoir FH)
  ffassert(0); // a faire
}

template<class K>
AnyType E_set_fev<K>::OpS(Stack s)  const
{
    //  voir E_set_fevS  ( pb de consitance a revoir FH)
    ffassert(0); // a faire
}

template<class K>
AnyType E_set_fev<K>::Op2d(Stack s)  const
{
  StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);
  MeshPoint *mps=MeshPointStack(s), mp=*mps;
  FEbase<K,v_fes> ** pp=GetAny< FEbase<K,v_fes> **>((*ppfe)(s));
  FEbase<K,v_fes> & fe(**pp);
  const  FESpace & Vh(*fe.newVh());
  KN<K> gg(Vh.MaximalNbOfDF());


  const  Mesh & Th(Vh.Th);
  const int dim=Vh.N;
  K ** copt=0;
  if (optimize)   copt= new K *[dim];
  if(copt) {
    assert((size_t) dim== where_in_stack_opt.size());
    for (int i=0;i<dim;i++)
      {
        int offset=where_in_stack_opt[i];
        assert(offset>10);
        copt[i]= static_cast<K *>(static_cast<void *>((char*)s+offset));
        *(copt[i])=0;
       }
    if (optiexp0) (*optiexp0)(s); // init
  }

#ifdef OLDPih
  ffassert(dim<100);
#endif

  TabFuncArg tabexp(s,Vh.N);
  ffassert( aa.size() == Vh.N);
  for (int i=0;i<dim;i++)
    tabexp[i]=aa[i];

  KN<K> * y=new  KN<K>(Vh.NbOfDF);
  KN<K> & yy(*y);

  FElement::aIPJ ipj(Vh[0].Pi_h_ipj());
  FElement::aR2  PtHat(Vh[0].Pi_h_R2());

  KN<double>   Aipj(ipj.N());

  KN<K>  Vp1(dim*PtHat.N());


  if (Vh.isFEMesh() )
    {

      ffassert(Vh.NbOfDF == Th.nv && dim == 1 );
      for (int iv=0;iv<Th.nv;iv++)
	{
	  const E_F0 & ff(* (const  E_F0 *) aa[0]  ) ;
	  const Vertex & v(Th(iv));
	  int ik=Th.Contening(&v);
	  const Triangle & Kt(Th[ik]);
	  int il=-1;
	  if  ( &Kt[0] == &v) il=0;
	  if  ( &Kt[1] == &v) il=1;
	  if  ( &Kt[2] == &v) il=2;
	  assert(il>=0);
	  mps->set(Th,v,TriangleHat[il],Kt,v.lab);
	  if (copt) {
	    if (optiexpK) (*optiexpK)(s);
            yy[iv] =  *(copt[0]);
          }
          else
	    yy[iv] = GetAny<K>( ff(s) );
	  sptr->clean(); // modif FH mars 2006  clean Ptr
	}

    }
  else
    for (int t=0;t<Th.nt;t++)
      {
	FElement Kt(Vh[t]);
	int nbdf=Kt.NbDoF();

	gg=K();

#ifdef OLDPih
	// old method
	Kt.Pi_h(gg,F_Pi_h,F,&tabexp);
#else
	Kt.Pi_h(Aipj);

	for (int p=0;p<PtHat.N();p++)
          {
            mps->set(Kt.T(PtHat[p]),PtHat[p],Kt);

            KN_<K> Vpp(Vp1,SubArray(dim,p*dim)); // a Change FHHHHHHHH
            if (copt) { // optimize  version
             if (optiexpK) (*optiexpK)(s);
             for (int j=0;j<dim;j++)
               Vpp[j] = *(copt[j]);
            }
            else  // old version
            for (int j=0;j<dim;j++)
             if (tabexp[j])
               Vpp[j]=GetAny<K>( (*tabexp[j])(s) );
              else Vpp[j]=0;

           }

         for (int i=0;i<Aipj.N();i++)
          {
           const FElement::IPJ &ipj_i(ipj[i]);
             gg[ipj_i.i] += Aipj[i]*Vp1(ipj_i.j+ipj_i.p*dim); // index a la main
            sptr->clean(); // modif FH mars 2006  clean Ptr
          }
#endif

         for (int df=0;df<nbdf;df++)
           yy[Kt(df)] =  gg[df] ;
      }
    fe=y;
    if (copt) delete [] copt;
    *MeshPointStack(s) = mp;
    if(verbosity>1)
        ShowBound(*y,cout) << endl ;
    return Nothing;
  }


template<class K>
inline FEbase<K,v_fes> * MakePtrFE(pfes * const &  a){
  FEbase<K,v_fes> * p=new FEbase<K,v_fes>(a);
  return p ;}

template<class K>
inline FEbase<K,v_fes> ** MakePtrFE2(FEbase<K,v_fes> * * const &  p,pfes * const &  a){
  *p=new FEbase<K,v_fes>(a);
  return p ;}

template<class K>
inline FEbaseArray<K,v_fes> ** MakePtrFE3(FEbaseArray<K,v_fes> * * const &  p,pfes * const &  a,const long & N){
  *p=new FEbaseArray<K,v_fes>(a,N);
  return p ;}

template<class K>
class  OneOperatorMakePtrFE : public OneOperator
{
public:
  typedef  FEbase<K,v_fes> ** R;
  typedef pfes* B;
  class CODE : public E_F0mps
  {
  public:
    Expression fer,fes;
    E_set_fev<K> * e_set_fev;
    const E_Array * v;
    CODE(const basicAC_F0 & args)
      :
      fer(to<R>(args[0])),
      fes(to<B>(args[1])),
      e_set_fev(0)
    {
      if (BCastTo<K>(args[2]) )
	v = new E_Array(basicAC_F0_wa(to<K>(args[2])));
      else
	v = dynamic_cast<const E_Array *>( args[2].LeftValue() );
      if (!v) {
	cout << "Error: type of arg :" << *args[2].left()  << " in " << typeid(K).name() << " case " << endl;
	ErrorCompile(" We wait  a double/complex expression or a array expression",1);
      }
        e_set_fev=  new   E_set_fev<K>(v,fer,v_fes::dHat);

    }

    AnyType operator()(Stack stack)  const {
      R  p = GetAny<R>( (*fer)(stack));
      B  a = GetAny<B>( (*fes)(stack));
      *p=new FEbase<K,v_fes>(a);
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

template<class Result,class A>
class  OneOperator_Ptr_o_R: public OneOperator {
    typedef Result A::* ptr;
    ptr p;
    public:
    E_F0 * code(const basicAC_F0 & args) const
     { return  new E_F_A_Ptr_o_R<Result,A>(t[0]->CastTo(args[0]),p);}
    OneOperator_Ptr_o_R(ptr pp):
       OneOperator(atype<Result*>(),atype<A*>()),p(pp) {}
};

template<class K>  K  *PAddition(const K * a,const K *  b)  {return  new K(*a+*b);}


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
    Check( i != l->end() ,"Def CL Dirichet already exists");
    l->insert(make_pair(k,CastTo<double>(f)));
  }
};

//  pour stocker des expression de compilation


class Convect : public E_F0mps  { public:
    typedef double  Result; // return type
    Expression u,v,w,ff,dt;
    long state;
    static Expression ou,ov,ow,odt; // previous def
    static long  count;
    int d;
    Convect(const basicAC_F0 & args)  : u(0),v(0),w(0),ff(0),dt(0)
    {
      args.SetNameParam();
      const E_Array * a = dynamic_cast<const E_Array *>(args[0].LeftValue());
            ffassert(a);
	d= a->size();
       if (d == 3)
	   w= CastTo<double>((*a)[2]);
       else if (d != 2)
          { CompileError("convect vector have only 2 or 3 componant");}
       u= CastTo<double>((*a)[0]);
       v= CastTo<double>((*a)[1]);

       dt=CastTo<double>(args[1]);
       ff=CastTo<double>(args[2]);
       // save previous  state
        if( !(( ou && u->compare(ou)==0)  && (v->compare(ov)==0) && ( (w==0) || (w->compare(ov)==0))
              && (dt->compare(odt)==0) )) {
             count++;

        }
        state= count;// use of optim of convect
        if(verbosity>3)
          cout << "\n  -- Convert number of Convect case:  .... "<< state << endl;
        ou=u;
        ov=v;
        ow=w;
        odt=dt;
     }

    static ArrayOfaType  typeargs()
      { return  ArrayOfaType(atype<E_Array>(),atype<double>(),atype<double>());}

    static  E_F0 * f(const basicAC_F0 & args) { return new Convect(args);}
    AnyType operator()(Stack s) const ;
    AnyType eval2(Stack s) const ;
    AnyType eval3(Stack s) const ;
    AnyType eval3old(Stack s) const ;  // old version a supprime en  2017
    operator aType () const { return atype<Result>();}

};
Expression Convect::ou=0;
Expression Convect::ov=0;
Expression Convect::ow=0;
Expression Convect::odt=0;
long  Convect::count=0;
/// <<Plot>> used for the [[plot_keyword]]

class Plot :  public E_F0mps /* [[file:AFunction.hpp::E_F0mps]] */ {
public:
    typedef KN_<R>  tab;
    typedef KN<KN<R> >*  pttab;
    typedef pferbase sol;
    typedef pferbasearray asol;
    typedef pf3rbase sol3;
    typedef pf3rbasearray asol3;
    typedef pfSrbase solS;
    typedef pfSrbasearray asolS;
    typedef pfLrbase solL;
    typedef pfLrbasearray asolL;
    
    typedef pfecbase solc;
    typedef pfecbasearray asolc;
    typedef pf3cbase solc3;
    typedef pf3cbasearray asolc3;
    typedef pfScbase solcS;
    typedef pfScbasearray asolcS;
    typedef pfLcbase solcL;
    typedef pfLcbasearray asolcL;
    
    typedef long  Result;
    struct ListWhat {
	int what,i;
	int cmp[4];
	int n;
       union {
        long l[4];
	void * v[4];//  for
        const void * cv[4];//  for
       };
	pmesh th() { assert(v[0] && what==0); return static_cast<pmesh>(cv[0]);}
	pmesh3 th3() { assert(v[0] && what==5); return static_cast<pmesh3>(cv[0]);}
 	pmeshS thS() { assert(v[0] && what==50); return static_cast<pmeshS>(cv[0]);}
    pmeshL thL() { assert(v[0] && what==55); return static_cast<pmeshL>(cv[0]);}
        
	void Set(int nn=0,void **vv=0,int *c=0) {
	    cmp[0]=cmp[1]=cmp[2]=cmp[3]=-1;
	    v[0]=v[1]=v[2]=v[3]=0;
        n=nn;
	    for(int i=0;i<nn;++i)
	      {
		  if(c) cmp[i]=c[i];
		  if(vv) v[i]=vv[i];
	      }
	}
        void Set(int nn,const void **vv,int *c=0) {
            cmp[0]=cmp[1]=cmp[2]=cmp[3]=-1;
            v[0]=v[1]=v[2]=v[3]=0;
            cv[0]=cv[1]=cv[2]=cv[3]=0;
            n=nn;
            for(int i=0;i<nn;++i)
            {
                if(c) cmp[i]=c[i];
                if(vv) cv[i]=vv[i];
            }
        }

	ListWhat(int w=-1,int ii=-1 )
	: what(w),i(ii) {Set();}
	ListWhat(int what,int ii,int n,void ** f0,int *c)
	: what(what),i(ii){ Set(n,f0,c);}
        ListWhat(int what,int ii,int n,const void ** f0,int *c)
        : what(what),i(ii){ Set(n,f0,c);}

	ListWhat(int what,int ii,void * f0)
	: what(what),i(ii){ Set(1,&f0,0);}
        ListWhat(int what,int ii,const void * f0)
        : what(what),i(ii){ Set(1,&f0,0);}

        ListWhat(int what,int ii,int nn,long *f0)
        : what(what),i(ii),n(nn){
            cmp[0]=cmp[1]=cmp[2]=cmp[3]=-1;
            v[0]=v[1]=v[2]=v[3]=0;
            for(int i=0; i<n; ++i)
                l[i]=f0[i];
        }

	template<typename S>
	void eval(S *f,int *c)
	{ for(int i=0;i<3;++i) {
        f[i]= static_cast<S>(v[i]);
	    c[i]= cmp[i];  }
	}

	template<typename M>
	M eval()
	{ assert(v[0]);
	    return static_cast<M>(v[0]);
	}
	void eval(sol & f0,int & cmp0, sol &f1,int &cmp1)
	{
	  f0=static_cast<sol>(v[0]);
	  f1=static_cast<sol>(v[1]);
	  cmp0=cmp[0];
	  cmp1=cmp[1];
	}
    void eval(solS & f0,int & cmp0, solS &f1,int &cmp1)    // TODO a template func
    {
       f0=static_cast<solS>(v[0]);
       f1=static_cast<solS>(v[1]);
       cmp0=cmp[0];
       cmp1=cmp[1];
    }

    };

  /// <<Expression2>>
    struct Expression2
     {//  FH. change nov 2016  add one  expression  for  colored curve ...
	long what; // 0 mesh, 1 iso, 2 vector, 3 curve , 4 border , 5  mesh3, 6 iso 3d,
	// 7: vector 3d  ( +10 -> complex visu ???? )
	// 101 array of iso 2d  , 106 array of iso 3d  , 100  array of meshes
        // 103  array of curves ...
	bool composant;
	Expression e[4];
	Expression2() {e[0]=0;e[1]=0;e[2]=0;e[3]=0;composant=false;what=0;}
	Expression &operator[](int i){return e[i];}

        int EvalandPush(Stack s,int ii,vector<ListWhat> & ll,vector<AnyType> & lat ) const
         {  // add for curve ... and multi curve ...
             //  store date in lat..
             long  f[4];
             for(int i=0; i< 4;++i)
             {
                 f[i]=-1;
                 if( e[i] ) { // eval ..
                     f[i]= lat.size();
                     lat.push_back((*e[i])(s));
                 }
             }
             ll.push_back(ListWhat(what,ii,4,f));
             return 4;
         }

	template<class S>
	int EvalandPush(Stack s,int ii,vector<ListWhat> & ll ) const
	{
	    int n=-1;
	    S f[3]={0,0,0};
	    int cmp[3]={-1,-1,-1};

	    for(int i=0;i<3;++i)
		if (e[i]) {
		    if (!composant)
		      { pair<S,int> p= GetAny< pair<S,int> >((*e[i])(s));
			  n=i;cmp[i]=p.second;
			  f[i]= p.first;}
		    else { cmp[i]=0;
			f[i]=GetAny< S >((*e[i])(s));
			n=i;}
		}
	    ll.push_back(ListWhat(what,ii,n+1,f,cmp));
	    return n;}

         int AEvalandPush(Stack s,int ii,vector<ListWhat> & ll,vector<AnyType> & lat) const
         {

          pttab pt[4]={0,0,0,0};
          pt[0] =evalptt(0,s);
          pt[1]=evalptt(1,s);
          if (e[2]) pt[2]=evalptt(2,s);
         if (e[3]) pt[3]=evalptt(3,s);
         int kt = min(pt[0]->N() , pt[1]->N());

         for( int j=0; j<kt; ++j)
         {
          int what = 13;
          long  f[4];
         if(verbosity>99)
             cout << " plot : A curve "<< j << " " ;
          for(int k=0; k<4; ++k)
          {
              pttab ptk=pt[k];
              f[k]=-1;
              if(ptk)
              {
                  KN_<double> t=(*ptk)[j];
                  f[k]=lat.size();
                  if(verbosity>99) cout << " ("<<k << ") " << t.N() << " "<< f[k] ;

                  lat.push_back(SetAny<KN_<double> >(t));
              }
          }
         if(verbosity>99) cout << endl;

          ll.push_back(ListWhat(what,ii,4,f));
         }
         return 4;
         }


	template<class A,class S> // ok of mesh too because composant=true;
 	int AEvalandPush(Stack s,int ii,vector<ListWhat> & ll ) const
	{  typedef pair<A,int> PA;
	    int nn=-1;
	    A f[3];
	    union {
		S fj[3];
		void *fv[3];
	    };
	    f[0]=f[1]=f[2]=0;
	    int cmp[3]={-1,-1,-1};

	    for(int i=0;i<3;++i)
		if (e[i])
		  {
		    if (!composant)
		     { PA p= GetAny< PA >((*e[i])(s)); cmp[i]=p.second;f[i]=p.first; nn=i;}
	            else
		     { f[i]= GetAny< A >((*e[i])(s)); cmp[i]=0; nn=i;}
		  }
	    	else break;
	    nn++;
	    int n = f[0]->N;
	    if(verbosity>50) // add 01/2011 FH ????
	      cout << "add  N = " << n << " " << nn  << " "<< what << endl;
	    for(int j=0;j<n;++j)
	      {

		int m=-1;
		fj[0]=fj[1]=fj[2]=0;// clean
		for (int i=0;i<nn;++i)
		  {
		    fj[i]=  *f[i]->operator[](j);
		    if(fj[i] && fj[i]->x()) m=i;
		    else break;
		}
		if(m>=0)  {
		    ll.push_back(ListWhat(what%100,ii,m+1,fv,cmp));
		    if(verbosity>100)
		    cout << ".";
		}

	      }
	    if(verbosity>100)
	    cout << endl;
	    return nn;
	}
     template<class S>
       int MEvalandPush(Stack s,int ii,vector<ListWhat> & ll ) const
       {  typedef KN<S> * A;

	   A ath;

	   ath= GetAny< A >((*e[0])(s));
	   int n=0;
	   if(ath) n = ath->N();
	    S th;

	   for(int j=0;j<n;++j)
	     {
	       th= ath->operator[](j);
	       if(th)
		ll.push_back(ListWhat(what%100,ii,static_cast<const void *>(th)));

	     }
	   return n;
       }

	sol eval(int i,Stack s,int & cmp) const  {  cmp=-1;
	    if (e[i]) {
		if (!composant) {pfer p= GetAny< pfer >((*e[i])(s)); cmp=p.second;return p.first;}
		else {return GetAny< pferbase >((*e[i])(s));}
	    }
	    else return 0;}
	sol3 eval3(int i,Stack s,int & cmp) const  {  cmp=-1;
	    if (e[i]) {
		if (!composant) {pf3r p= GetAny< pf3r >((*e[i])(s)); cmp=p.second;return p.first;}
		else {return GetAny< pf3rbase >((*e[i])(s));}
	    }
	    else return 0;}
	// add FH Japon 2010 ..	for complex visu ...  to complex ....  try to uniformize ...
	solc evalc(int i,Stack s,int & cmp) const  {  cmp=-1;
	    if (e[i]) {
		if (!composant) {pfec p= GetAny< pfec >((*e[i])(s)); cmp=p.second;return p.first;}
		else {return GetAny< pfecbase >((*e[i])(s));}
	    }
	    else return 0;}
	solc3 evalc3(int i,Stack s,int & cmp) const  {  cmp=-1;
	    if (e[i]) {
		if (!composant) {pf3c p= GetAny< pf3c >((*e[i])(s)); cmp=p.second;return p.first;}
		else {return GetAny< pf3cbase >((*e[i])(s));}
	    }
	    else return 0;}


	asol evala(int i, Stack s,int & cmp) const  {  cmp=-1;
	    if (e[i])
	      {pferarray p= GetAny< pferarray >((*e[i])(s)); cmp=p.second;return p.first;}
	    else return 0;}
	asol3 evala3(int i, Stack s,int & cmp) const  {  cmp=-1;
	    if (e[i])
	      {pf3rarray p= GetAny< pf3rarray >((*e[i])(s)); cmp=p.second;return p.first;}
	    else return 0;}
   asolS evalaS(int i, Stack s,int & cmp) const  {  cmp=-1;
         if (e[i])
          {pfSrarray p= GetAny< pfSrarray >((*e[i])(s)); cmp=p.second;return p.first;}
        else return 0;}

	asolc evalca(int i, Stack s,int & cmp) const  {  cmp=-1;
	    if (e[i])
	      {pfecarray p= GetAny< pfecarray >((*e[i])(s)); cmp=p.second;return p.first;}
	    else return 0;}
	asolc3 evalca3(int i, Stack s,int & cmp) const  {  cmp=-1;
	    if (e[i])
	      {pf3carray p= GetAny< pf3carray >((*e[i])(s)); cmp=p.second;return p.first;}
	    else return 0;}
    asolcS evalcaS(int i, Stack s,int & cmp) const  {  cmp=-1;
         if (e[i])
          {pfScarray p= GetAny< pfScarray >((*e[i])(s)); cmp=p.second;return p.first;}
        else return 0;}

	const Mesh & evalm(int i,Stack s) const  { throwassert(e[i]);return  * GetAny< pmesh >((*e[i])(s)) ;}
	KN<pmesh> * evalma(int i,Stack s) const  { throwassert(e[i]);return   GetAny< KN<pmesh> * >((*e[i])(s)) ;}
	const Mesh3 & evalm3(int i,Stack s) const  { throwassert(e[i]);return  * GetAny< pmesh3 >((*e[i])(s)) ;}
    const MeshS & evalmS(int i,Stack s) const  { throwassert(e[i]);return  * GetAny< pmeshS >((*e[i])(s)) ;}
    const MeshL & evalmL(int i,Stack s) const  { throwassert(e[i]);return  * GetAny< pmeshL >((*e[i])(s)) ;}
	const E_BorderN * evalb(int i,Stack s) const  { throwassert(e[i]);return   GetAny< const E_BorderN *>((*e[i])(s)) ;}
    const E_Curve3N * evalc(int i,Stack s) const  { throwassert(e[i]);return   GetAny< const E_Curve3N *>((*e[i])(s)) ;}
	tab  evalt(int i,Stack s) const  { throwassert(e[i]);return  GetAny<tab>((*e[i])(s)) ;}
         pttab  evalptt(int i,Stack s) const  { throwassert(e[i]);return  GetAny<pttab>((*e[i])(s)) ;}
    };

  // see [[Plot_name_param]]
  static basicAC_F0::name_and_type name_param[] ;

  /// <<number_of_distinct_named_parameters_for_plot>> FFCS: added new parameters for VTK graphics. See
  /// [[Plot_name_param]] for new parameter names
  static const int n_name_param=43;

  Expression bb[4];

  /// [[Expression2]] is a description of an object to plot
  vector<Expression2> l;
    typedef KN<KN<double> > * ptaboftab;
    Expression nargs[n_name_param];

    Plot(const basicAC_F0 & args) : l(args.size())
    {
        args.SetNameParam(n_name_param,name_param,nargs);
        if ( nargs[8] )
            Box2x2( nargs[8] , bb);

        // scan all the parameters of the plot() call
        for (size_t i=0;i<l.size();i++)

            // argument is an [[file:AFunction.hpp::E_Array]] (= array of E_F0)
            if (args[i].left()==atype<E_Array>())
            {
                l[i].composant=false;
                const E_Array * a = dynamic_cast<const E_Array *>(args[i].LeftValue());
                ffassert(a);
                int asizea=a->size();
                if(asizea==0) CompileError("plot of vector with 0 of components(!= 2 or 3) ");
                bool bpfer=  BCastTo<pfer>((*a)[0]);
                bool bpf3r=  BCastTo<pf3r>((*a)[0]);
                bool bpfSr=  BCastTo<pfSr>((*a)[0]); // for 3D surface with reals
                bool bpfLr=  BCastTo<pfLr>((*a)[0]); // for 3D curve with reals
                bool bpfec=  BCastTo<pfec>((*a)[0]);
                bool bpf3c=  BCastTo<pf3c>((*a)[0]);
                bool bpfSc=  BCastTo<pfSc>((*a)[0]); // for 3D surface with complex
                bool bpfLc=  BCastTo<pfLc>((*a)[0]); // for 3D curve with complex
                bool bptab=  BCastTo<tab> ((*a)[0]);
                bool bpttab=  BCastTo<pttab>((*a)[0]);
                bool bpferarray=  BCastTo<pferarray>((*a)[0]);
                bool bpfecarray=  BCastTo<pfecarray>((*a)[0]);
                if ( bpfer && asizea <3)
                {
                    l[i].what=asizea;
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pfer>((*a)[j]);
                }
                else if ( bpfec && asizea <3)
                {
                    l[i].what=10+asizea;
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pfec>((*a)[j]);
                }
                else if (bptab && asizea>=2 && asizea<=4)// change nov 2016 FH  for color corve
                {
                    l[i].what=3;
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<tab>((*a)[j]);
                }
                else if (bpttab && asizea>=2 && asizea<=4)// change nov 2016 FH  for  arry of curve
                {
                    l[i].what=103;
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pttab>((*a)[j]);
                }
                else if (asizea == 3 && bpf3r ) // 3d vector ...
                {
                    l[i].what=7; // new 3d vector
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pf3r>((*a)[j]);

                }
                else if (asizea == 3 && bpf3c ) // 3d vector ...
                {
                    l[i].what=17; // new 3d vector
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pf3c>((*a)[j]);

                }
                else if (asizea == 3 && bpfSr ) // 3D vector for surface...
                {
                    l[i].what=9; // new 3d vector
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pfSr>((*a)[j]);

                }
                else if (asizea == 3 && bpfSc ) // 3D vector for surface...
                {
                    l[i].what=19; // new 3d vector
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pfSc>((*a)[j]);

                }
                else if (asizea == 3 && bpfLr ) // 3D vector for curve...
                {
                    l[i].what=15; // new 3d vector
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pfLr>((*a)[j]);
                    
                }
                else if (asizea == 3 && bpfLc ) // 3D vector for curve...
                {
                    l[i].what=21; // new 3d vector
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pfLc>((*a)[j]);
                    
                }
                else if (bpferarray && asizea == 2)
                {
                    l[i].what=102;
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pferarray>((*a)[j]);
                }
                else if (bpfecarray && asizea == 2)
                {
                    l[i].what=112;
                    for (int j=0;j<a->size();j++)
                        l[i][j]= CastTo<pfecarray>((*a)[j]);
                }
                else { CompileError("plot of array with wrong  number of components (!= 2 or 3) ");}
            }
            else if (BCastTo<pferbase>(args[i])) { // [[file:problem.hpp::pferbase]] [[file:lgmesh3.hpp::BCastTo]]
                l[i].what=1; //  iso value 2d
                l[i].composant=true;
                l[i][0]=CastTo<pferbase>(args[i]); }
            else if (BCastTo<pfer>(args[i])) { // [[file:problem.hpp::pfer]]
                l[i].composant=false;
                l[i].what=1; //  iso value 2d
                l[i][0]=CastTo<pfer>(args[i]);}
            else if (BCastTo<pfecbase>(args[i])) { // [[file:problem.hpp::pfecbase]]
                l[i].what=11; //  iso value 2d
                l[i].composant=true;
                l[i][0]=CastTo<pfecbase>(args[i]); }
            else if (BCastTo<pfec>(args[i])) { // [[file:problem.hpp::pfec]]
                l[i].composant=false;
                l[i].what=11; //  iso value 2d
                l[i][0]=CastTo<pfec>(args[i]);}


            else if (BCastTo<pf3r>(args[i])) { // [[file:lgmesh3.hpp::pf3r]]
                l[i].composant=false;
                l[i].what=6; //  real iso value 3D volume
                l[i][0]=CastTo<pf3r>(args[i]);}
            else if (BCastTo<pf3c>(args[i])) { // [[file:lgmesh3.hpp::pf3c]]
                l[i].composant=false;
                l[i].what=16; //  complex iso value 3D volume
                l[i][0]=CastTo<pf3c>(args[i]);}

            else if (BCastTo<pfSr>(args[i])) { // [[file:lgmesh3.hpp::pfSr]]
                l[i].composant=false;
                l[i].what=8; //  real iso value 3D surface
                l[i][0]=CastTo<pfSr>(args[i]);}
            else if (BCastTo<pfSc>(args[i])) { // [[file:lgmesh3.hpp::pfSc]]
                l[i].composant=false;
                l[i].what=18;   // complex iso value 3D surface
                l[i][0]=CastTo<pfSc>(args[i]);}
            else if (BCastTo<pfLr>(args[i])) { // [[file:lgmesh3.hpp::pfSr]]
                l[i].composant=false;
                l[i].what=14; //  real iso value 3D curve
                l[i][0]=CastTo<pfLr>(args[i]);}
            else if (BCastTo<pfLc>(args[i])) { // [[file:lgmesh3.hpp::pfSc]]
                l[i].composant=false;
                l[i].what=20;   // complex iso value 3D curve
                l[i][0]=CastTo<pfLc>(args[i]);}

            else if (BCastTo<pferarray>(args[i])) {
                l[i].composant=false;
                l[i].what=101; //  iso value array iso value 2d
                l[i][0]=CastTo<pferarray>(args[i]);}
            else if (BCastTo<pfecarray>(args[i])) {
                l[i].composant=false;
                l[i].what=111; //  iso value array iso value 2d
                l[i][0]=CastTo<pfecarray>(args[i]);}

            else if (BCastTo<pf3rarray>(args[i])) { // [[file:lgmesh3.hpp::pf3rarray]]
                l[i].composant=false;
                l[i].what=106; // iso value array iso value 3d
                l[i][0]=CastTo<pf3rarray>(args[i]);}
            else if (BCastTo<pf3carray>(args[i])) { // [[file:lgmesh3.hpp::pf3carray]]
                l[i].composant=false;
                l[i].what=116; // iso value array iso value 3d
                l[i][0]=CastTo<pf3carray>(args[i]);}

            else if (BCastTo<pfSrarray>(args[i])) { // [[file:lgmesh3.hpp::pfSrarray]]
                l[i].composant=false;
                l[i].what=108; //arry iso value array iso value 3d
                l[i][0]=CastTo<pfSrarray>(args[i]);}
            else if (BCastTo<pfScarray>(args[i])) { // [[file:lgmesh3.hpp::pfScarray]]
                l[i].composant=false;
                l[i].what=118;  //arry iso value array iso value 3d
                l[i][0]=CastTo<pfScarray>(args[i]);}
            else if (BCastTo<pfLrarray>(args[i])) { // [[file:lgmesh3.hpp::pfSrarray]]
                l[i].composant=false;
                l[i].what=109; //arry iso value array iso value 3d
                l[i][0]=CastTo<pfLrarray>(args[i]);}
            else if (BCastTo<pfLcarray>(args[i])) { // [[file:lgmesh3.hpp::pfScarray]]
                l[i].composant=false;
                l[i].what=119;  //arry iso value array iso value 3d
                l[i][0]=CastTo<pfLcarray>(args[i]);}
            else if (BCastTo<pmesh>(args[i])){
                l[i].composant=true;
                l[i].what=0; // mesh ...
                l[i][0]=CastTo<pmesh>(args[i]);}
            else if (BCastTo<pmesh3>(args[i])){
                l[i].composant=true;
                l[i].what=5;// 3d mesh (real volume or if contains a meshS pointer...
                l[i][0]=CastTo<pmesh3>(args[i]);}
            else if (BCastTo<pmeshS>(args[i])){
                l[i].composant=true;
                l[i].what=50;// 3d surface mesh
                l[i][0]=CastTo<pmeshS>(args[i]);}
            else if (BCastTo<pmeshL>(args[i])){
                l[i].composant=true;
                l[i].what=55;// 3d line mesh
                l[i][0]=CastTo<pmeshL>(args[i]);}
            else if (BCastTo<const E_BorderN *>(args[i])){
                l[i].what=4; // border 2d
                l[i].composant=true;
                l[i][0]=CastTo<const E_BorderN *>(args[i]);}
            else if (BCastTo<const E_Curve3N *>(args[i])){
                l[i].what=40; // border 3d
                l[i].composant=true;
                l[i][0]=CastTo<const E_Curve3N *>(args[i]);}
            else if (BCastTo<KN<pmesh> *>(args[i])){
                l[i].composant=true;
                l[i].what=100; //  mesh 2d array
                l[i][0]=CastTo<KN<pmesh> *>(args[i]);}

            else {
                CompileError("Sorry no way to plot this kind of data");
            }

    }


  static ArrayOfaType  typeargs() { return  ArrayOfaType(true);}// all type

  /// <<Plot_f>> Creates a Plot object with the list of arguments obtained from the script during the grammatical
  /// analysis of the script (in lg.ypp)

    static  E_F0 * f(const basicAC_F0 & args) {;return new Plot(args);}

  /// Evaluates the contents of the Plot object during script evaluation. Implemented at [[Plot_operator_brackets]]

  AnyType operator()(Stack s) const ;
};

/// <<Plot_name_param>>

basicAC_F0::name_and_type Plot::name_param[Plot::n_name_param] = {
  {   "coef", &typeid(double)},
  {   "cmm", &typeid(string*)},
  {   "ps", &typeid(string*)  },
  {   "wait", &typeid(bool) },
  {   "fill", &typeid(bool) },
  {   "value", &typeid(bool) },
  {   "clean", &typeid(bool) },
  {   "aspectratio", &typeid(bool)},
  {   "bb",&typeid(E_Array) },
  {   "nbiso", &typeid(long)},
  {   "nbarrow", &typeid(long)},
  {   "viso", &typeid(KN_<double>)},
  {   "varrow", &typeid(KN_<double>)},
  {   "bw",&typeid(bool)},
  {   "grey", &typeid(bool)},
  {   "hsv", &typeid(KN_<double>)},
  {   "boundary", &typeid(bool)}, // 16
  {   "dim", &typeid(long)}, // 2 or 3
  {   "add", &typeid(bool)}, // add to previous plot
  {   "prev", &typeid(bool)}, // keep previou  view point
  {   "ech", &typeid(double)}, // keep previou  view point

  // FFCS: more options for VTK graphics (numbers are required for processing)

  {"ZScale",&typeid(double)}, // #1
  {"WhiteBackground",&typeid(bool)}, // #2
  {"OpaqueBorders",&typeid(bool)}, // #3
  {"BorderAsMesh",&typeid(bool)}, // #4
  {"ShowMeshes",&typeid(bool)}, // #5
  {"ColorScheme",&typeid(long)}, // #6
  {"ArrowShape",&typeid(long)}, // #7
  {"ArrowSize",&typeid(double)}, // #8
  {"ComplexDisplay",&typeid(long)}, // #9
  {"LabelColors",&typeid(bool)}, // #10
  {"ShowAxes",&typeid(bool)}, // #11
  {"CutPlane",&typeid(bool)}, // #12
  {"CameraPosition",&typeid(KN_<double>)}, // #13
  {"CameraFocalPoint",&typeid(KN_<double>)}, // #14
  {"CameraViewUp",&typeid(KN_<double>)}, // #15
  {"CameraViewAngle",&typeid(double)}, // #16
  {"CameraClippingRange",&typeid(KN_<double>)}, // #17
  {"CutPlaneOrigin",&typeid(KN_<double>)}, // #18
  {"CutPlaneNormal",&typeid(KN_<double>)}, // #19
  {"WindowIndex",&typeid(long)}, // #20
  {"NbColorTicks",&typeid(long)}, // #21
  {"NbColors",&typeid(long)} // #22
};



template<class K>
class pb2mat : public E_F0 { public:
  typedef Matrice_Creuse<K> *  Result;
  const Problem * pb;
  pb2mat(const basicAC_F0 & args) : pb(dynamic_cast<const Problem *>(args[0].left()))
  {ffassert(pb);}
  static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<const Problem *>());}

  static  E_F0 * f(const basicAC_F0 & args) { return new Plot(args);}

  AnyType operator()(Stack s) const
  {
    Problem::Data<FESpace> *data= pb->dataptr(this->stack);
    if ( SameType<K,double>::OK )
      {
	ffassert( !!data->AR);
	return  SetAny<Matrice_Creuse<K> * >(&data->AR) ;
      }
    else
      {
	ffassert( !!data->AC);
	return SetAny<Matrice_Creuse<K> * >(&data->AC) ;
      }
  }


};

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
   nu_face= make_Type_Expr(atype<long>(),new E_P_Stack_Nu_Face);
   lenEdge    = make_Type_Expr(atype<R>(),new E_P_Stack_lenEdge);
   hTriangle  = make_Type_Expr(atype<R>(),new E_P_Stack_hTriangle);
   area       = make_Type_Expr(atype<R>(),new E_P_Stack_areaTriangle);
   volume       = make_Type_Expr(atype<R>(),new E_P_Stack_VolumeTet);
   inside     = make_Type_Expr(atype<R>(),new E_P_Stack_inside);
  Global.New("x",x);
  Global.New("y",y);
  Global.New("z",z);
  Global.New("label",label);
  Global.New("region",region);
  Global.New("notaregion",CConstant<long>(lnotaregion));
  Global.New("nuTriangle",nu_triangle);
  Global.New("nuTet",nu_triangle);
  Global.New("nuEdge",nu_edge);
  Global.New("nuFace",nu_face);
  Global.New("P",P);
  Global.New("N",N);

  Global.New("lenEdge",lenEdge);
  Global.New("area",area);
  Global.New("volume",volume);
  Global.New("hTriangle",hTriangle);
  Global.New("inside",inside);
  Global.New("nTonEdge",make_Type_Expr(atype<long>(),new E_P_Stack_nTonEdge));
  Global.New("nElementonB",make_Type_Expr(atype<long>(),new E_P_Stack_nElementonB));
  Global.New("edgeOrientation",make_Type_Expr(atype<R>(),new E_P_Stack_EdgeOrient)); // Add FH jan 2018
  Global.New("BoundaryEdge",make_Type_Expr(atype<long>(),new E_P_Stack_TypeEdge<1>));// Add FH jan 2018
  Global.New("InternalEdge",make_Type_Expr(atype<long>(),new E_P_Stack_TypeEdge<2>));// Add FH jan 2018


}


template<class K>
struct set_eqmatrice_creuse_fbl: public binary_function<Matrice_Creuse<K>*,const Matrice_Creuse<K>  *,const  C_args * > {
  static Matrice_Creuse<K>* f(Matrice_Creuse<K>* const  & a,const  C_args * const & b)  {
  // 1  verif the FESpace

  // 2 set = or +=

    ffassert(0);
    return a;}
};

template<class K>
struct set_eqvect_fl: public binary_function<KN<K>*,const  FormLinear *,KN<K>*> {
  static KN<K>* f(KN<K>* const  & a,const  FormLinear * const & b)  {
    ffassert(0);
    return a;}
};

 template<class R>
  AnyType IntFunction<R>::operator()(Stack stack) const  {
  MeshPoint mp=* MeshPointStack(stack);
  StackOfPtr2Free  * wsptr2free=WhereStackOfPtr2Free(stack);
  size_t  swsptr2free =wsptr2free->size();
 R r=0;

 SHOWVERB(cout << " int " << endl);
 const vector<Expression>  & what(di->what);
 const int dim =di->d;
 const bool surface=di->isMeshS;
 const bool curve=di->isMeshL;
 const GQuadratureFormular<R1>& FIE = di->FIE(stack);
 const GQuadratureFormular<R2> & FIT = di->FIT(stack);
 const GQuadratureFormular<R3> & FIV = di->FIV(stack);

 CDomainOfIntegration::typeofkind kind = di->kind;
 set<int> setoflab;
 bool all=true;
  if( di->withmap()) { ExecError(" no map  in the case (1)??");}
  if (verbosity>3)
  {
      if(dim==2)
      {
          if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") levelset: "<< di->islevelset() << " ,"  ;
          else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
          else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
          else cout << "  --  int 2d   (nQP: "<< FIT.n << " ) in "  ;
      }
      else if(dim==3 && !surface && !curve)
      {
          if (CDomainOfIntegration::int2d==kind) cout << "  -- boundary int border ( nQP: "<< FIT.n << ") ,"  ;
          else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all faces ( nQP: "<< FIT.n << "),"  ;
          else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF face nQP: ("<< FIT.n << ")," ;
          else cout << "  --  int 3d volume  (nQP: "<< FIV.n << " ) in "  ;
      }
      else if(dim==3 && surface && !curve)
      {
          if (CDomainOfIntegration::int1d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") levelset: "<< di->islevelset() << " ,"  ;
          else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
          else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
          else cout << "  --  int 3d surface  (nQP: "<< FIT.n << " ) in "  ;
      }
      else if(dim==3 && !surface  && curve)
      {
          //if (CDomainOfIntegration::int0d==kind) cout << "  -- boundary int border ( nQP: "<< FIE.n << ") levelset: "<< di->islevelset() << " ,"  ;
          //else  if (CDomainOfIntegration::intalledges==kind) cout << "  -- boundary int all edges ( nQP: "<< FIE.n << "),"  ;
          //else  if (CDomainOfIntegration::intallVFedges==kind) cout << "  -- boundary int all VF edges nQP: ("<< FIE.n << ")," ;
          cout << "  --  int 3d curve  (nQP: "<< FIT.n << " ) in "  ;
      }
  }

 Expandsetoflab(stack,*di, setoflab,all);

 if(dim==2)
   {
     if(di->islevelset() && (CDomainOfIntegration::int1d!=kind) &&  (CDomainOfIntegration::int2d!=kind) ) InternalError("So no levelset integration type  case (10 2d)");
     const Mesh  & Th = * GetAny<pmesh>( (*di->Th)(stack) );
     ffassert(&Th);

     if (verbosity >3)
     {
       if (all) cout << " all " << endl ;
       else cout << endl;
     }
     if (kind==CDomainOfIntegration::int1d)
       {
	 const QuadratureFormular1d & FI = FIE;
           if(di->islevelset())
           {
               double llevelset = 0;
               double uset = HUGE_VAL;
               R2 Q[3];
               KN<double> phi(Th.nv);phi=uset;
               double f[3];
               for(int t=0; t< Th.nt;++t)
               {
                   double umx=-HUGE_VAL,umn=HUGE_VAL;
                   for(int i=0;i<3;++i)
                   {
                       int j= Th(t,i);
                       if( phi[j]==uset)
                       {
                           MeshPointStack(stack)->setP(&Th,t,i);
                           phi[j]= di->levelset(stack);//zzzz
                       }
                       f[i]=phi[j];
                       umx = std::max(umx,phi[j]);
                       umn = std::min(umn,phi[j]);

                   }
                   if( umn <=0 && umx >= 0)
                   {

                       int np= IsoLineK(f,Q,1e-10);
                       if(np==2)
                       {
                           const Triangle & K(Th[t]);
                           R2 PA(K(Q[0])),PB(K(Q[1]));
                           R2 NAB(PA,PB);
                           double  lAB=sqrt((NAB,NAB));
                           NAB = NAB.perp()/lAB;
                           llevelset += lAB;
                           for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                           {
                               QuadratureFormular1dPoint pi( FI[npi]);
                               double sa=pi.x,sb=1.-sa;
                               R2 Pt(Q[0]*sa+Q[1]*sb ); //
                               MeshPointStack(stack)->set(Th,K(Pt),Pt,K,-1,NAB,-1);
                               r += lAB*pi.a*GetAny<R>( (*fonc)(stack));
                           }
                       }

                   }
                   wsptr2free->clean(swsptr2free);// ADD FH 11/2017
               }
               if(verbosity > 5) cout << " Lenght level set = " << llevelset << endl;

           }

        else
	 for( int e=0;e<Th.neb;e++)
	   {
	     if (all || setoflab.find(Th.bedges[e].lab) != setoflab.end())
	       {

		 int ie,i =Th.BoundaryElement(e,ie);
		 const Triangle & K(Th[i]);
		 R2 E=K.Edge(ie);
		 double le = sqrt((E,E));
		 R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
		   PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);

		 for (int npi=0;npi<FI.n;npi++) // loop on the integration point
		   {
		     QuadratureFormular1dPoint pi( FI[npi]);
		     double sa=pi.x,sb=1.-sa;
		     R2 Pt(PA*sa+PB*sb ); //
		     MeshPointStack(stack)->set(Th,K(Pt),Pt,K,Th.bedges[e].lab,R2(E.y,-E.x)/le,ie);
		     r += le*pi.a*GetAny<R>( (*fonc)(stack));
		   }
	       }
               wsptr2free->clean(swsptr2free);// ADD FH 11/2017
	   }
       }
     else if (kind==CDomainOfIntegration::int2d)
     {
         const QuadratureFormular & FI =FIT;

         if(di->islevelset())
         { // add FH mars 2014 compute int2d  on phi < 0 ..
             double llevelset = 0;
             double uset = HUGE_VAL;
             R2 Q[3];
             KN<double> phi(Th.nv);phi=uset;
             double f[3],umx,umn;
             for(int t=0; t< Th.nt;++t)
             {
                if (all || setoflab.find(Th[t].lab) != setoflab.end())
                {
                const Triangle & K(Th[t]);

                 double umx=-HUGE_VAL,umn=HUGE_VAL;
                 for(int i=0;i<3;++i)
                 {
                     int j= Th(t,i);
                     if( phi[j]==uset)
                     {
                         MeshPointStack(stack)->setP(&Th,t,i);
                         phi[j]= di->levelset(stack);//zzzz
                     }
                     f[i]=phi[j];
                     umx = std::max(umx,phi[j]);
                     umn = std::min(umn,phi[j]);

                 }
                 double area =K.area;
                 if( umn >=0 ) continue; //  all positif => nothing
                 if( umx >0 )
                    { // coupe ..
                        int i0 = 0, i1 = 1, i2 =2;

                        if( f[i0] > f[i1] ) swap(i0,i1) ;
                        if( f[i0] > f[i2] ) swap(i0,i2) ;
                        if( f[i1] > f[i2] ) swap(i1,i2) ;

                        double c = (f[i2]-f[i1])/(f[i2]-f[i0]); // coef Up Traing
                        if( f[i1] < 0 ) {double y=f[i2]/(f[i2]-f[i1]); c *=y*y; }
                        else {double y=f[i0]/(f[i0]-f[i1]) ; c = 1.- (1.-c)*y*y; };
                        assert( c > 0 && c < 1);
                        area *= 1-c;
                    }
                 //  warning  quadrature  wrong just ok for constante FH, we must also change the quadaturer points ..
                 // just order 1  here ???
                 for (int npi=0; npi<FI.n;npi++)
                    {
                        QuadraturePoint pi(FI[npi]);
                        MeshPointStack(stack)->set(Th,K(pi),pi,K,K.lab);
                        r += area*pi.a*GetAny<R>( (*fonc)(stack));
                    }

                }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }

         }
      else
       for (int i=0;i< Th.nt; i++)
	 {
	   const Triangle & K(Th[i]);
	   if (all || setoflab.find(Th[i].lab) != setoflab.end())
	     for (int npi=0; npi<FI.n;npi++)
	       {
		 QuadraturePoint pi(FI[npi]);
		 MeshPointStack(stack)->set(Th,K(pi),pi,K,K.lab);
		 r += K.area*pi.a*GetAny<R>( (*fonc)(stack));
	       }
             wsptr2free->clean(swsptr2free);// ADD FH 11/2017
	 }
     }
     else   if (kind==CDomainOfIntegration::intalledges)
       {
	 const QuadratureFormular1d & FI = FIE;
	 for (int i=0;i< Th.nt; i++)
	   if (all || setoflab.find(Th[i].lab) != setoflab.end())
           {
	     for( int ie=0;ie<3;ie++)
	       {
		 const Triangle & K(Th[i]);
                   int e0=VerticesOfTriangularEdge[ie][0];
                   int e1=VerticesOfTriangularEdge[ie][1];
                   int i1 = Th(K[e0]),i2 = Th(K[e1]);
                   BoundaryEdge * be = Th.TheBoundaryEdge(i1,i2);
                   int lab = be ? be->lab :  notaregion;
		 R2 E=K.Edge(ie);
		 double le = sqrt((E,E));
		 R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
		   PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);

		 for (int npi=0;npi<FI.n;npi++) // loop on the integration point
		   {
		     QuadratureFormular1dPoint pi( FI[npi]);
		     double sa=pi.x,sb=1-sa;
		     R2 Pt(PA*sa+PB*sb ); //
		     MeshPointStack(stack)->set(Th,K(Pt),Pt,K,lab,R2(E.y,-E.x)/le,ie);// correction FH 6/2/2014
		     r += le*pi.a*GetAny<R>( (*fonc)(stack));
		   }
	       }
               wsptr2free->clean(swsptr2free);// ADD FH 11/2017
           }
       }
     else   if (kind==CDomainOfIntegration::intallVFedges)
       {
	 double untier(1./3.);
	 cerr << " a faire CDomainOfIntegration::intallVFedges " << endl; //%%%%%%%%%
	 ffassert(0);
	 const QuadratureFormular1d & FI = FIE;
	 for (int i=0;i< Th.nt; i++)
	   if (all || setoflab.find(Th[i].lab) != setoflab.end())
	     {
	       const Triangle & K(Th[i]);
	       const R2 GH(untier,untier);
	       const R2 G=K(GH);
	       for( int ie=0;ie<3;ie++)
		 {
		   int ie0=VerticesOfTriangularEdge[ie][0] ;
		   int ie1=VerticesOfTriangularEdge[ie][1] ;
		   const R2 MH=(TriangleHat[ie0]+TriangleHat[ie1])*0.5;
		   const R2 M(K(MH));
		   R2 E(G,M);
		   double le = sqrt((E,E));

		   for (int npi=0;npi<FI.n;npi++) // loop on the integration point
		     {
		       QuadratureFormular1dPoint pi( FI[npi]);
		       double sa=pi.x,sb=1-sa;
		       R2 Pt(GH*sa+MH*sb ); //
		       MeshPointStack(stack)->set(Th,K(Pt),Pt,K,Th[ie].lab,R2(E.y,-E.x)/le,ie,1);
		       r += le*pi.a*GetAny<R>( (*fonc)(stack));
		     }
		 }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
	     }
       }
     else
       {
	 InternalError("CDomainOfIntegration kind unkown");
       }
   }

 // volume 3d
 else if(dim==3 && !surface  && !curve)
   {
     if(di->islevelset() &&  (CDomainOfIntegration::int2d!=kind) && (CDomainOfIntegration::int3d!=kind) )
         InternalError("So no levelset integration type on no int2d / int3d case (10 3d)");

     const Mesh3  & Th = * GetAny<pmesh3>( (*di->Th)(stack) );
     ffassert(&Th);

     if (verbosity >3)
       {
	 if (all) cout << " all " << endl ;
	 else cout << endl;
       }
       if (kind==CDomainOfIntegration::int2d)
           if(di->islevelset())
           {
               const GQuadratureFormular<R2> & FI = FIT;
               double llevelset = 0;
               const double uset = std::numeric_limits<double>::max();
               R3 Q[4];
               KN<double> phi(Th.nv);
               phi=uset;
               double f[4];

               for(int t=0; t< Th.nt;++t)
               {
                   double umx=std::numeric_limits<double>::min(),umn=std::numeric_limits<double>::max();
                   for(int i=0;i<4;++i)
                   {
                       int j= Th(t,i);
                       if( phi[j]==uset)
                       {
                           MeshPointStack(stack)->setP(&Th,t,i);
                           phi[j]= di->levelset(stack);//zzzz
                       }
                       f[i]=phi[j];
                       umx = std::max(umx,f[i]);
                       umn = std::min(umn,f[i]);

                   }
                   if( umn <=0 && umx >= 0)
                   {

                       int np= IsoLineK(f,Q,1e-10);

                       double l[3];
                       if(np>2)
                       {
                        if(verbosity>999)  cout << t << " int levelset : " << umn << " .. " << umx << " np " << np <<" "
                         << f[0] << " " << f[1] << " "<< f[2] << " "<< f[3] << " "<<endl;

                           const Mesh3::Element & K(Th[t]);
                           double epsmes3=K.mesure()*K.mesure()*1e-18;
                           R3 PP[4];
                           for(int i=0; i< np; ++i)
                               PP[i]= K(Q[i]);
                           for( int i =0; i+1 < np; i+=2)
                           { // 0,1,, a and 2,3,0.
                               int i0=i,i1=i+1,i2=(i+2)%np;
                               R3 NN= R3(PP[i0],PP[i1])^R3(PP[i0],PP[i2]);
                               double mes2 = (NN,NN);
                               double mes = sqrt(mes2);
                               if(mes2*mes <epsmes3) continue; //  too small
                               NN /= mes;
                               mes *= 0.5; //   warning correct FH 050109
                               llevelset += mes;
                               for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                               {
                                   GQuadraturePoint<R2>  pi( FI[npi]);
                                   pi.toBary(l);
                                   R3 Pt( l[0]*Q[i0]+l[1]*Q[i1]+l[2]*Q[i2]); //
                                   MeshPointStack(stack)->set(Th,K(Pt),Pt,K,-1,NN,-1);
                                   r += mes*pi.a*GetAny<R>( (*fonc)(stack));
                               }
                           }
                       }

                   }
                   wsptr2free->clean(swsptr2free);// ADD FH 11/2017
               }
               if(verbosity > 5) cout << " Area level set = " << llevelset << endl;

           }

         else

       {
	 const GQuadratureFormular<R2> & FI = FIT;
         int lab;
	 for( int e=0;e<Th.nbe;e++)
	   {
	     if (all || setoflab.find(lab=Th.be(e).lab) != setoflab.end())
	       {

		 int ie,i =Th.BoundaryElement(e,ie);
		 const Mesh3::Element & K(Th[i]);
		 R3 NN=K.N(ie);
		 double mes = sqrt((NN,NN));
		 NN /= mes;
		 mes *= 0.5; //   warning correct FH 050109
		 for (int npi=0;npi<FI.n;npi++) // loop on the integration point
		   {
		     GQuadraturePoint<R2> pi( FI[npi]);
		     R3 Pt(K.PBord(ie,pi)); //
		     MeshPointStack(stack)->set(Th,K(Pt),Pt,K,lab,NN,ie);
		     r += mes*pi.a*GetAny<R>( (*fonc)(stack));
		   }
                   wsptr2free->clean(swsptr2free);// ADD FH 11/2017
	       }
	   }
       }
     else if (kind==CDomainOfIntegration::int3d)
     {

         if(di->islevelset())
         {
             GQuadratureFormular<R3> FI(FIV.n*3);
             double llevelset = 0;
             const double uset = std::numeric_limits<double>::max();
             R3 Q[3][4];
             double vol6[3];
             KN<double> phi(Th.nv);
             phi=uset;
             double f[4];

             for (int t=0;t< Th.nt; t++)
             {

		 const Mesh3::Element & K(Th[t]);
		 if (all || setoflab.find(K.lab) != setoflab.end())
                 {
                     double umx=std::numeric_limits<double>::min(),umn=std::numeric_limits<double>::max();
                     for(int i=0;i<4;++i)
                     {
                         int j= Th(t,i);
                         if( phi[j]==uset)
                         {
                             MeshPointStack(stack)->setP(&Th,t,i);
                             phi[j]= di->levelset(stack);//zzzz
                         }
                         f[i]=phi[j];
                     }
                     int ntets= UnderIso(f,Q, vol6,1e-14);
                     setQF<R3>(FI,FIV,QuadratureFormular_Tet_1, Q,vol6,ntets);
                     for (int npi=0; npi<FI.n;npi++)
                     {
                         GQuadraturePoint<R3> pi(FI[npi]);
                         MeshPointStack(stack)->set(Th,K(pi),pi,K,K.lab);
                         r += K.mesure()*pi.a*GetAny<R>( (*fonc)(stack));
                     }
                 }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }


         }
        else
         {
       const GQuadratureFormular<R3> & FI =FIV;
             for (int i=0;i< Th.nt; i++)
	       {
		 const Mesh3::Element & K(Th[i]);
		 if (all || setoflab.find(K.lab) != setoflab.end())
		   for (int npi=0; npi<FI.n;npi++)
                     {
                       GQuadraturePoint<R3> pi(FI[npi]);
                       MeshPointStack(stack)->set(Th,K(pi),pi,K,K.lab);
                       r += K.mesure()*pi.a*GetAny<R>( (*fonc)(stack));
                     }
                   wsptr2free->clean(swsptr2free);// ADD FH 11/2017
               }
         }
     }
     else   if (kind==CDomainOfIntegration::intalledges)
       {
         const GQuadratureFormular<R2> & FI = FIT;
         int lab;
	 for (int i=0;i< Th.nt; i++)
	   if (all || setoflab.find(Th[i].lab) != setoflab.end())
	     for( int ie=0;ie<4;ie++) // Coorection mai 2018 FH ???????????????? never tested
	       {

		 const Mesh3::Element  & K(Th[i]);
		 R3 NN=K.N(ie);
		 double mes = NN.norme();
		   NN /= mes;
		   mes*=0.5;//  correction 05/01/09 FH
		 for (int npi=0;npi<FI.n;npi++) // loop on the integration point
		   {
		     GQuadraturePoint<R2> pi( FI[npi]);
		     R3 Pt(K.PBord(ie,pi)); //
		     MeshPointStack(stack)->set(Th,K(Pt),Pt,K,lab,NN,ie);
		     r += mes*pi.a*GetAny<R>( (*fonc)(stack));
		   }
	       }
           wsptr2free->clean(swsptr2free);// ADD FH 11/2017
       }

   }

 else if(dim==3 && surface  && !curve)
 {
     if(di->islevelset() && (CDomainOfIntegration::int1d!=kind) &&  (CDomainOfIntegration::int2d!=kind) ) InternalError("So no levelset integration type  case (10 3d)");
     const MeshS  & Th = * GetAny<pmeshS>( (*di->Th)(stack) );
     ffassert(&Th);

     if (verbosity >3)
     {
         if (all) cout << " all " << endl ;
         else cout << endl;
     }
     if (kind==CDomainOfIntegration::int1d)
     {
         const QuadratureFormular1d & FI = FIE;
         if(di->islevelset())
         {
             double llevelset = 0;
             double uset = HUGE_VAL;
             R2 Q[3];
             KN<double> phi(Th.nv);phi=uset;
             double f[3];
             for(int t=0; t< Th.nt;++t)
             {
                 double umx=-HUGE_VAL,umn=HUGE_VAL;
                 for(int i=0;i<3;++i)
                 {
                     int j= Th(t,i);
                     if( phi[j]==uset)
                     {
                         MeshPointStack(stack)->setP(&Th,t,i);
                         phi[j]= di->levelset(stack);
                     }
                     f[i]=phi[j];
                     umx = std::max(umx,phi[j]);
                     umn = std::min(umn,phi[j]);

                 }
                 if( umn <=0 && umx >= 0)
                 {

                     int np= IsoLineK(f,Q,1e-10);
                     if(np==2)
                     {
                         const TriangleS & K(Th[t]);
                         double epsmes3=K.mesure()*K.mesure()*1e-18;
                         R3 PA(K(Q[0])),PB(K(Q[1])), PC(K(Q[2]));
                         R3 NN= R3(PB,PA)^R3(PC,PA); //R3 NAB(PA,PB);

                         double mes2 = (NN,NN);
                         double mes = sqrt(mes2);

                         if(mes2*mes <epsmes3) continue; //  too small
                         NN /= mes;
                         llevelset += mes;

                         for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                         {
                             QuadratureFormular1dPoint pi( FI[npi]);
                             double sa=pi.x,sb=1.-sa;
                             R2 Pt(Q[0]*sa+Q[1]*sb ); //
                             MeshPointStack(stack)->set(Th,K(Pt),Pt,K,-1,NN,-1);
                             r += mes*pi.a*GetAny<R>( (*fonc)(stack));

                         }
                     }

                 }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }
             if(verbosity > 5) cout << " Lenght level set = " << llevelset << endl;

         }

         else
             for( int e=0;e<Th.nbe;e++)
             {
                 if (all || setoflab.find(Th.be(e).lab) != setoflab.end())
                 {

                     int ie,i =Th.BoundaryElement(e,ie);
                     const TriangleS & K(Th[i]);
                     R3 E=K.Edge(ie);
                     double le = sqrt((E,E));
                     R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
                     PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);

                     for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                     {
                         QuadratureFormular1dPoint pi( FI[npi]);
                         double sa=pi.x,sb=1.-sa;
                         R2 Pt(PA*sa+PB*sb ); //
                         MeshPointStack(stack)->set(Th,K(Pt),Pt,K,Th.be(e).lab,R2(E.y,-E.x)/le,ie);
                         r += le*pi.a*GetAny<R>( (*fonc)(stack));
                     }
                 }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }
     }
     else if (kind==CDomainOfIntegration::int2d)
     {
         const QuadratureFormular & FI =FIT;

         if(di->islevelset())
         { // add FH mars 2014 compute int2d  on phi < 0 ..
             double llevelset = 0;
             double uset = HUGE_VAL;
             R2 Q[3];
             KN<double> phi(Th.nv);phi=uset;
             double f[3],umx,umn;
             for(int t=0; t< Th.nt;++t)
             {
                 if (all || setoflab.find(Th[t].lab) != setoflab.end())
                 {
                     const TriangleS & K(Th[t]);

                     double umx=-HUGE_VAL,umn=HUGE_VAL;
                     for(int i=0;i<3;++i)
                     {
                         int j= Th(t,i);
                         if( phi[j]==uset)
                         {
                             MeshPointStack(stack)->setP(&Th,t,i);
                             phi[j]= di->levelset(stack);
                         }
                         f[i]=phi[j];
                         umx = std::max(umx,phi[j]);
                         umn = std::min(umn,phi[j]);

                     }
                     double area =K.mesure();
                     if( umn >=0 ) continue; //  all positif => nothing
                     if( umx >0 )
                     { // coupe ..
                         int i0 = 0, i1 = 1, i2 =2;

                         if( f[i0] > f[i1] ) swap(i0,i1) ;
                         if( f[i0] > f[i2] ) swap(i0,i2) ;
                         if( f[i1] > f[i2] ) swap(i1,i2) ;

                         double c = (f[i2]-f[i1])/(f[i2]-f[i0]); // coef Up Traing
                         if( f[i1] < 0 ) {double y=f[i2]/(f[i2]-f[i1]); c *=y*y; }
                         else {double y=f[i0]/(f[i0]-f[i1]) ; c = 1.- (1.-c)*y*y; };
                         assert( c > 0 && c < 1);
                         area *= 1-c;
                     }
                     //  warning  quadrature  wrong just ok for constante FH, we must also change the quadaturer points ..
                     // just order 1  here ???
                     for (int npi=0; npi<FI.n;npi++)
                     {
                         QuadraturePoint pi(FI[npi]);
                         MeshPointStack(stack)->set(Th,K(pi),pi,K,K.lab);
                         r += area*pi.a*GetAny<R>( (*fonc)(stack));
                     }

                 }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }

         }
         else
             for (int i=0;i< Th.nt; i++)
             {
                 const TriangleS & K(Th[i]);
                 if (all || setoflab.find(Th[i].lab) != setoflab.end())
                     for (int npi=0; npi<FI.n;npi++)
                     {
                         QuadraturePoint pi(FI[npi]);
                         MeshPointStack(stack)->set(Th,K(pi),pi,K,K.lab);
                         r += K.mesure()*pi.a*GetAny<R>( (*fonc)(stack));
                     }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }
     }
     else   if (kind==CDomainOfIntegration::intalledges)
     {
         const QuadratureFormular1d & FI = FIE;
         int lab;
         for (int i=0;i< Th.nt; i++)
             if (all || setoflab.find(Th[i].lab) != setoflab.end())
                 for( int ie=0;ie<3;ie++)
                 {

                     const MeshS::Element  & K(Th[i]);
                     R3 NN=K.N(ie);
                     double mes = NN.norme();
                     NN /= mes;
                     for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                    {
                         QuadratureFormular1dPoint pi( FI[npi]);
                         R2 Pt(K.PBord(ie,pi)); // cout
                         MeshPointStack(stack)->set(Th,K(Pt),Pt,K,lab,NN,ie);
                         r += mes*pi.a*GetAny<R>( (*fonc)(stack));
                     }
                 }
         wsptr2free->clean(swsptr2free);// ADD FH 11/2017
     }
     else   if (kind==CDomainOfIntegration::intallVFedges)
     {
         ffassert(0); //TODO AXEL
     }
     else
     {
         InternalError("CDomainOfIntegration kind unkown");
     }
 }
 else if(dim==3 && !surface  && curve)
 {
     ffassert(0);
     /*if(di->islevelset() && (CDomainOfIntegration::int1d!=kind) &&  (CDomainOfIntegration::int2d!=kind) ) InternalError("So no levelset integration type  case (10 3d)");
     const MeshS  & Th = * GetAny<pmeshS>( (*di->Th)(stack) );
     ffassert(&Th);
     
     if (verbosity >3)
     {
         if (all) cout << " all " << endl ;
         else cout << endl;
     }
     if (kind==CDomainOfIntegration::int1d)
     {
         const QuadratureFormular1d & FI = FIE;
         if(di->islevelset())
         {
             double llevelset = 0;
             double uset = HUGE_VAL;
             R2 Q[3];
             KN<double> phi(Th.nv);phi=uset;
             double f[3];
             for(int t=0; t< Th.nt;++t)
             {
                 double umx=-HUGE_VAL,umn=HUGE_VAL;
                 for(int i=0;i<3;++i)
                 {
                     int j= Th(t,i);
                     if( phi[j]==uset)
                     {
                         MeshPointStack(stack)->setP(&Th,t,i);
                         phi[j]= di->levelset(stack);
                     }
                     f[i]=phi[j];
                     umx = std::max(umx,phi[j]);
                     umn = std::min(umn,phi[j]);
                     
                 }
                 if( umn <=0 && umx >= 0)
                 {
                     
                     int np= IsoLineK(f,Q,1e-10);
                     if(np==2)
                     {
                         const TriangleS & K(Th[t]);
                         double epsmes3=K.mesure()*K.mesure()*1e-18;
                         R3 PA(K(Q[0])),PB(K(Q[1])), PC(K(Q[2]));
                         R3 NN= R3(PB,PA)^R3(PC,PA); //R3 NAB(PA,PB);
                         
                         double mes2 = (NN,NN);
                         double mes = sqrt(mes2);
                         
                         if(mes2*mes <epsmes3) continue; //  too small
                         NN /= mes;
                         llevelset += mes;
                         
                         for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                         {
                             QuadratureFormular1dPoint pi( FI[npi]);
                             double sa=pi.x,sb=1.-sa;
                             R2 Pt(Q[0]*sa+Q[1]*sb ); //
                             MeshPointStack(stack)->set(Th,K(Pt),Pt,K,-1,NN,-1);
                             r += mes*pi.a*GetAny<R>( (*fonc)(stack));
                             
                         }
                     }
                     
                 }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }
             if(verbosity > 5) cout << " Lenght level set = " << llevelset << endl;
             
         }
         
         else
             for( int e=0;e<Th.nbe;e++)
             {
                 if (all || setoflab.find(Th.be(e).lab) != setoflab.end())
                 {
                     
                     int ie,i =Th.BoundaryElement(e,ie);
                     const TriangleS & K(Th[i]);
                     R3 E=K.Edge(ie);
                     double le = sqrt((E,E));
                     R2 PA(TriangleHat[VerticesOfTriangularEdge[ie][0]]),
                     PB(TriangleHat[VerticesOfTriangularEdge[ie][1]]);
                     
                     for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                     {
                         QuadratureFormular1dPoint pi( FI[npi]);
                         double sa=pi.x,sb=1.-sa;
                         R2 Pt(PA*sa+PB*sb ); //
                         MeshPointStack(stack)->set(Th,K(Pt),Pt,K,Th.be(e).lab,R2(E.y,-E.x)/le,ie);
                         r += le*pi.a*GetAny<R>( (*fonc)(stack));
                     }
                 }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }
     }
     else if (kind==CDomainOfIntegration::int2d)
     {
         const QuadratureFormular & FI =FIT;
         
         if(di->islevelset())
         { // add FH mars 2014 compute int2d  on phi < 0 ..
             double llevelset = 0;
             double uset = HUGE_VAL;
             R2 Q[3];
             KN<double> phi(Th.nv);phi=uset;
             double f[3],umx,umn;
             for(int t=0; t< Th.nt;++t)
             {
                 if (all || setoflab.find(Th[t].lab) != setoflab.end())
                 {
                     const TriangleS & K(Th[t]);
                     
                     double umx=-HUGE_VAL,umn=HUGE_VAL;
                     for(int i=0;i<3;++i)
                     {
                         int j= Th(t,i);
                         if( phi[j]==uset)
                         {
                             MeshPointStack(stack)->setP(&Th,t,i);
                             phi[j]= di->levelset(stack);
                         }
                         f[i]=phi[j];
                         umx = std::max(umx,phi[j]);
                         umn = std::min(umn,phi[j]);
                         
                     }
                     double area =K.mesure();
                     if( umn >=0 ) continue; //  all positif => nothing
                     if( umx >0 )
                     { // coupe ..
                         int i0 = 0, i1 = 1, i2 =2;
                         
                         if( f[i0] > f[i1] ) swap(i0,i1) ;
                         if( f[i0] > f[i2] ) swap(i0,i2) ;
                         if( f[i1] > f[i2] ) swap(i1,i2) ;
                         
                         double c = (f[i2]-f[i1])/(f[i2]-f[i0]); // coef Up Traing
                         if( f[i1] < 0 ) {double y=f[i2]/(f[i2]-f[i1]); c *=y*y; }
                         else {double y=f[i0]/(f[i0]-f[i1]) ; c = 1.- (1.-c)*y*y; };
                         assert( c > 0 && c < 1);
                         area *= 1-c;
                     }
                     //  warning  quadrature  wrong just ok for constante FH, we must also change the quadaturer points ..
                     // just order 1  here ???
                     for (int npi=0; npi<FI.n;npi++)
                     {
                         QuadraturePoint pi(FI[npi]);
                         MeshPointStack(stack)->set(Th,K(pi),pi,K,K.lab);
                         r += area*pi.a*GetAny<R>( (*fonc)(stack));
                     }
                     
                 }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }
             
         }
         else
             for (int i=0;i< Th.nt; i++)
             {
                 const TriangleS & K(Th[i]);
                 if (all || setoflab.find(Th[i].lab) != setoflab.end())
                     for (int npi=0; npi<FI.n;npi++)
                     {
                         QuadraturePoint pi(FI[npi]);
                         MeshPointStack(stack)->set(Th,K(pi),pi,K,K.lab);
                         r += K.mesure()*pi.a*GetAny<R>( (*fonc)(stack));
                     }
                 wsptr2free->clean(swsptr2free);// ADD FH 11/2017
             }
     }
     else   if (kind==CDomainOfIntegration::intalledges)
     {
         const QuadratureFormular1d & FI = FIE;
         int lab;
         for (int i=0;i< Th.nt; i++)
             if (all || setoflab.find(Th[i].lab) != setoflab.end())
                 for( int ie=0;ie<3;ie++)
                 {
                     
                     const MeshS::Element  & K(Th[i]);
                     R3 NN=K.N(ie);
                     double mes = NN.norme();
                     NN /= mes;
                     for (int npi=0;npi<FI.n;npi++) // loop on the integration point
                     {
                         QuadratureFormular1dPoint pi( FI[npi]);
                         R2 Pt(K.PBord(ie,pi)); // cout
                         MeshPointStack(stack)->set(Th,K(Pt),Pt,K,lab,NN,ie);
                         r += mes*pi.a*GetAny<R>( (*fonc)(stack));
                     }
                 }
         wsptr2free->clean(swsptr2free);// ADD FH 11/2017
     }
     else   if (kind==CDomainOfIntegration::intallVFedges)
     {
         ffassert(0); //TODO AXEL
     }
     else
     {
         InternalError("CDomainOfIntegration kind unkown");
     }*/
 }
 else {  InternalError("CDomainOfIntegration dim unkown");}

 *MeshPointStack(stack)=mp;
     return SetAny<R>(r);
  }


void Show(const char * s,int k=1)
{
  if(k) {
  couleur(1);
  float xmin,xmax,ymin,ymax;
       getcadre(xmin,xmax,ymin,ymax);
       rmoveto(xmin+(xmax-xmin)/100,ymax-(k)*(ymax-ymin)/30);
       plotstring(s);
  }
}
template<class K,class v_fes>
int Send2d(PlotStream & theplot,Plot::ListWhat & lli,map<const typename v_fes::FESpace::Mesh *,long> & mapth)
{
    typedef FEbase<K,v_fes> * pfek ;
    pfek fe[3]={0,0,0};
    int cmp[3]={-1,-1,-1};
    int err=1;
    long what=lli.what;
    int lg,nsb;
    lli.eval(fe,cmp);
    if (fe[0]->x() && what %10 ==1)
	{
	      err=0;
	      theplot << what ;
	      theplot <<mapth[ &(fe[0]->Vh->Th)];// numero du maillage

	      KN<K> V1=fe[0]->Vh->newSaveDraw(*fe[0]->x(),cmp[0],lg,nsb);

	      // construction of the sub division ...
	      int nsubT=NbOfSubTriangle(nsb);
	      int nsubV=NbOfSubInternalVertices(nsb);
	      KN<R2> Psub(nsubV);
	      KN<int> Ksub(nsubT*3);
	      for(int i=0;i<nsubV;++i)
	      Psub[i]=SubInternalVertex(nsb,i);
	      for(int sk=0,p=0;sk<nsubT;++sk)
	      for(int i=0;i<3;++i,++p)
	      Ksub[p]=numSubTriangle(nsb,sk,i);

	      if(verbosity>9)
	      cout << " Send plot:what: " << what << " " << nsb << " "<< V1.N()
	      << " Max "  << V1.max() << " min " << V1.min() << endl;
	      theplot << Psub ;
	      theplot << Ksub ;
	      theplot << V1;

	      }
    else if (fe[0]->x() && fe[1]->x() &&what %10 ==2)
      {

	{
	  err=0;
	  theplot << what ;


	  KN<K> V1=fe[0]->Vh->newSaveDraw(*fe[0]->x(),*fe[1]->x(),cmp[0],cmp[1],lg,nsb);
	  // construction of the sub division ...
	  int nsubT=NbOfSubTriangle(nsb);
	  int nsubV=NbOfSubInternalVertices(nsb);
	  KN<R2> Psub(nsubV);
	  KN<int> Ksub(nsubT*3);
	  for(int i=0;i<nsubV;++i)
	      Psub[i]=SubInternalVertex(nsb,i);
	  for(int sk=0,p=0;sk<nsubT;++sk)
	      for(int i=0;i<3;++i,++p)
		  Ksub[p]=numSubTriangle(nsb,sk,i);

	  theplot <<mapth[ &(fe[0]->Vh->Th)];// numero du maillage
	  theplot << Psub ;
	  theplot << Ksub ;
	  theplot <<  V1;
	}

      }
    return err;

}

template<class K,class v_fes>
int Send3d(PlotStream & theplot,Plot::ListWhat &lli,map<const typename v_fes::FESpace::Mesh *,long> &mapth3)
{
    typedef FEbase<K,v_fes> * pfek3 ;
    pfek3 fe3[3]={0,0,0};
    int cmp[3]={-1,-1,-1};
    int err=1;
    long what=lli.what;
    int lg,nsb;
    if (what%10==6)
    {
        int lg,nsb;
        lli.eval(fe3,cmp);
        // FFCS is able to display 3d complex data
        {
            if (fe3[0]->x())
            {
                err=0;
                theplot << what ;
                theplot <<mapth3[ &(fe3[0]->Vh->Th)];// numero du maillage
                KN<R3> Psub;
                KN<int> Ksub;
                KN<K> V1=fe3[0]->Vh->newSaveDraw(*fe3[0]->x(),cmp[0],lg,Psub,Ksub,0);
                if(verbosity>9)
                    cout << " Send plot:what: " << what << " " << nsb << " "<< V1.N()
                    << " "  << V1.max() << " " << V1.min() << endl;
                theplot << Psub ;
                theplot << Ksub ;
                theplot << V1;
            }
        }
    }
    else  if (what%10==7)
    {
        int lg,nsb;
        lli.eval(fe3,cmp);
        // FFCS is able to display 3d complex data
        {
            if (fe3[0]->x()&& fe3[1]->x() && fe3[2]->x())
            {
                err=0;
                theplot << what ;
                theplot <<mapth3[ &(fe3[0]->Vh->Th)];// numero du maillage
                KN<R3> Psub1,Psub2,Psub3; // bf Bof ...
                KN<int> Ksub1,Ksub2,Ksub3;
                KN<K> V1=fe3[0]->Vh->newSaveDraw(*fe3[0]->x(),cmp[0],lg,Psub1,Ksub1,0);
                KN<K> V2=fe3[1]->Vh->newSaveDraw(*fe3[1]->x(),cmp[1],lg,Psub2,Ksub2,0);
                KN<K> V3=fe3[2]->Vh->newSaveDraw(*fe3[2]->x(),cmp[2],lg,Psub3,Ksub3,0);
                if(verbosity>9)
                    cout << " Send plot:what: " << what << " " << nsb << " "<< V1.N()
                    << " "  << V1.max() << " " << V1.min() << endl;
                theplot << Psub1 ;
                theplot << Ksub1 ;
                ffassert( V1.N() == V2.N()  &&V1.N() == V3.N());
                KNM<K> V123(3,V1.N()); // warning fortran numbering ...
                V123(0,'.')=V1;
                V123(1,'.')=V2;
                V123(2,'.')=V3;
                // FFCS: should be able to deal with complex as well
                theplot << (KN_<K>&) V123;

            }
        }
    }
    return err;
}



template<class K,class v_fes>
int SendS(PlotStream & theplot,Plot::ListWhat &lli,map<const MeshS *,long> &mapthS)
{
    typedef FEbase<K,v_fes> * pfekS ;
    pfekS feS[3]={0,0,0};
    int cmp[3]={-1,-1,-1};
    int err=1;
    long what=lli.what;
    int lg,nsb=0;
    lli.eval(feS,cmp);


    if (what%10==8)
    {
        err=0;
        theplot << what ;
        theplot <<mapthS[ &(feS[0]->Vh->Th)];// numero du maillage
        KN<R2> Psub;
        KN<int> Ksub;
        KN<K> V1=feS[0]->Vh->newSaveDraw(*feS[0]->x(),cmp[0],lg,Psub,Ksub,0);

        if(verbosity>9)
            cout << " Send plot:what: " << what << " " << nsb << " "<<  V1.N()
            << " Max "  << V1.max() << " min " << V1.min() << endl;
        theplot << Psub ;
        theplot << Ksub ;
        theplot << V1;

    }
     else  if (what%10==9)
     { ffassert(0);
     }
    return err;
}

template<class K,class v_fes>
int SendL(PlotStream & theplot,Plot::ListWhat &lli,map<const MeshL *,long> &mapthS)
{
    typedef FEbase<K,v_fes> * pfekS ;
    pfekS feL[3]={0,0,0};
    int cmp[3]={-1,-1,-1};
    int err=1;
    long what=lli.what;
    int lg,nsb=0;
    lli.eval(feL,cmp);
    
    
    if (what==14 || what==20)
    {
        err=0;
        theplot << what ;
        theplot <<mapthS[ &(feL[0]->Vh->Th)];// numero du maillage
        KN<R1> Psub;
        KN<int> Ksub;
        KN<K> V1=feL[0]->Vh->newSaveDraw(*feL[0]->x(),cmp[0],lg,Psub,Ksub,0);
        
        if(verbosity>9)
            cout << " Send plot:what: " << what << " " << nsb << " "<<  V1.N()
            << " Max "  << V1.max() << " min " << V1.min() << endl;
        theplot << Psub ;
        theplot << Ksub ;
        theplot << V1;
        
    }
    else  if (what==15 || what==21)
    { ffassert(0);
    }
    return err;
}
//  missing function
inline  void NewSetColorTable(int nb,float *colors=0,int nbcolors=0,bool hsv=true)
{
    if(colors && nbcolors)
        SetColorTable1(nb,hsv,nbcolors,colors);
    else
        SetColorTable(nb);
}


/// <<Plot_operator_brackets>> from class [[Plot]]
AnyType Plot::operator()(Stack s) const{

    // remap  case 107 and 108 , 109  for array of FE.
    vector<ListWhat> ll;
    vector<AnyType> lat;
    ll.reserve(l.size());
    // generation de la list de plot ...
    for (size_t i=0;i<l.size();i++)
    {

        switch (l[i].what) {
            case 0 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 1 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 2 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 3 : l[i].EvalandPush(s,i,ll,lat); break;
            case 5 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 6 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 7 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 8 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 9 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 11 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 12 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 14 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 15 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 16 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 17 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 18 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 19 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 20 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 21 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 50 : l[i].EvalandPush<void *>(s,i,ll); break;
            case 55 : l[i].EvalandPush<void *>(s,i,ll); break;
            case  100 : l[i].MEvalandPush< pmesh>(s,i,ll); break;
            case  101 : l[i].AEvalandPush<asol, sol>(s,i,ll); break;
            case  102 : l[i].AEvalandPush<asol, sol>(s,i,ll); break;
            case  103 : l[i].AEvalandPush(s,i,ll,lat); break;
            case  105 : l[i].MEvalandPush< pmesh3>(s,i,ll); break;
            case  106 : l[i].AEvalandPush<asol3, sol3>(s,i,ll); break;
            case  108 : l[i].AEvalandPush<asolS, solS>(s,i,ll); break;
            case  109 : l[i].AEvalandPush<asolL, solL>(s,i,ll); break;
            case  111 : l[i].AEvalandPush<asolc, solc>(s,i,ll); break;
            case  112 : l[i].AEvalandPush<asolc, solc>(s,i,ll); break;
            case  116 : l[i].AEvalandPush<asolc3, solc3>(s,i,ll); break;
            case  118 : l[i].AEvalandPush<asolcS, solcS>(s,i,ll); break;
            case  119 : l[i].AEvalandPush<asolcL, solcL>(s,i,ll); break;

            default:
                ffassert(l[i].what<100) ; // missing piece of code FH (jan 2010) ...
                ll.push_back(ListWhat(l[i].what,i));
                break;
        }

    }

    if(ThePlotStream)
    {
        /*
         les different item of the plot are given by the number what:
         what = 0 -> mesh
         what = 1 -> real scalar field (FE function  2d)
         what = 2 -> 2d vector field (two FE function  2d)
         what = 3 -> curve def by 2,.., 4 array
         what = 4 -> border 2d
         what = 5 -> 3d mesh, real volume or if contains a surface mesh
         what = 6 -> real FE function 3D volume
         what = 7 -> real 3d vector field (tree FE function  3d) volume
         what = 8 -> real FE function 3D surface
         what = 9 -> real 3d vector field (tree FE function  3d) surface
         what = 11 -> complex scalar field (FE function  2d) real
         what = 13 -> curve def by 4  array : x,y,z,  value for color ...
         what = 14 -> real FE function 3D curve
         what = 15 -> real 3d vector field (tree FE function  3d) curve
         what = 16 -> complex FE function 3D volume
         what = 17 -> complex 3d vector field (tree FE function  3d) volume
         what = 18 -> complex FE function 3D surface
         what = 19 -> complex 3d vector field (tree FE function  3d) surface
         what = 20 -> complex FE function 3D curve
         what = 21 -> complex 3d vector field (tree FE function  3d) curve
         what = 40 -> 3D border mesh
         what = 50 -> 3D surface mesh
         what = 55 -> 3D line mesh
         what = 100,101,106,109 ->   remap with real ... 2d, 3D volume, 3D surface, 3D curve
         what = 111, 116, 117 ,118,119 ->  remap with complex ... 2d, 3D volume, 3D surface, 3D curve
         what = -1 -> error, item empty
         */
        PlotStream theplot(ThePlotStream);
        pferbase  fe[3]={0,0,0};
        pf3rbase  fe3[3]={0,0,0};
        pfSrbase  feS[3]={0,0,0};
        pfLrbase  feL[3]={0,0,0};
        double echelle=1;
        int cmp[3]={-1,-1,-1};
        theplot.SendNewPlot();
        if (nargs[0]) (theplot<< 0L) <=  GetAny<double>((*nargs[0])(s));
        if (nargs[1]) (theplot<< 1L) <=  GetAny<string *>((*nargs[1])(s));
        if (nargs[2]) (theplot<< 2L) <=  GetAny<string*>((*nargs[2])(s));
        if (nargs[3]) (theplot<< 3L)  <= (bool) (!NoWait &&  GetAny<bool>((*nargs[3])(s)));
        else (theplot<< 3L)  <=  (bool)  (TheWait&& !NoWait);
        if (nargs[4]) (theplot<< 4L)  <= GetAny<bool>((*nargs[4])(s));
        if (nargs[5]) (theplot<< 5L) <=  GetAny<bool>((*nargs[5])(s));
        if (nargs[6]) (theplot<< 6L) <=  GetAny<bool>((*nargs[6])(s));
        if (nargs[7]) (theplot<< 7L)  <= GetAny<bool>((*nargs[7])(s));
        if (nargs[8])
        {  KN<double> bbox(4);
            for (int i=0;i<4;i++)
                bbox[i]= GetAny<double>((*bb[i])(s));

            (theplot<< 8L) <= bbox ;
        }
        if (nargs[9])  (theplot<< 9L)   <=  GetAny<long>((*nargs[9])(s));
        if (nargs[10])  (theplot<< 10L) <= GetAny<long>((*nargs[10])(s));
        if (nargs[11]) {
            KN_<double> v =GetAny<KN_<double> >((*nargs[11])(s)) ;
            (theplot<< 11L)  <= v   ;}

        if (nargs[12])
            (theplot<< 12L) <=  GetAny<KN_<double> >((*nargs[12])(s)) ;



        if (nargs[13]) (theplot<< 13L)  <= GetAny<bool>((*nargs[13])(s));
        if (nargs[14]) (theplot<< 14L) <= GetAny<bool>((*nargs[14])(s));
        if (nargs[15])
            (theplot<< 15L)  <= GetAny<KN_<double> >((*nargs[15])(s));
        if (nargs[16]) (theplot<< 16L)  <= GetAny<bool>((*nargs[16])(s));
        // add frev 2008 FH for 3d plot ...
        if (nargs[17]) (theplot<< 17L)  <= GetAny<long>((*nargs[17])(s));
        if (nargs[18]) (theplot<< 18L)  <= GetAny<bool>((*nargs[18])(s));
        if (nargs[19]) (theplot<< 19L)  <= GetAny<bool>((*nargs[19])(s));
        if (nargs[20]) (theplot<< 20L)  <= (echelle=GetAny<double>((*nargs[20])(s)));

        // FFCS: extra plot options for VTK (indexed from 1 to keep these lines unchanged even if the number of standard
        // FF parameters above changes) received in [[file:../ffcs/src/visudata.cpp::receiving_plot_parameters]]. When
        // adding a parameter here, do _NOT_ forget to change the size of the array at
        // [[number_of_distinct_named_parameters_for_plot]] and to name the new parameters at [[Plot_name_param]]. Also
        // update the list of displayed values at [[file:../ffcs/src/plot.cpp::Plotparam_listvalues]] and read the
        // parameter value from the pipe at [[file:../ffcs/src/visudata.cpp::receiving_plot_parameters]].

#define VTK_START 20
#define SEND_VTK_PARAM(index,type)                    \
if(nargs[VTK_START+index])                    \
(theplot<<(long)(VTK_START+index))                \
<=GetAny<type>((*nargs[VTK_START+index])(s));

        SEND_VTK_PARAM(1,double); // ZScale
        SEND_VTK_PARAM(2,bool); // WhiteBackground
        SEND_VTK_PARAM(3,bool); // OpaqueBorders
        SEND_VTK_PARAM(4,bool); // BorderAsMesh
        SEND_VTK_PARAM(5,bool); // ShowMeshes
        SEND_VTK_PARAM(6,long); // ColorScheme
        SEND_VTK_PARAM(7,long); // ArrowShape
        SEND_VTK_PARAM(8,double); // ArrowSize
        SEND_VTK_PARAM(9,long); // ComplexDisplay
        SEND_VTK_PARAM(10,bool); // LabelColors
        SEND_VTK_PARAM(11,bool); // ShowAxes
        SEND_VTK_PARAM(12,bool); // CutPlane
        SEND_VTK_PARAM(13,KN_<double>); // CameraPosition
        SEND_VTK_PARAM(14,KN_<double>); // CameraFocalPoint
        SEND_VTK_PARAM(15,KN_<double>); // CameraViewUp
        SEND_VTK_PARAM(16,double); // CameraViewAngle
        SEND_VTK_PARAM(17,KN_<double>); // CameraClippingRange
        SEND_VTK_PARAM(18,KN_<double>); // CutPlaneOrigin
        SEND_VTK_PARAM(19,KN_<double>); // CutPlaneNormal
        SEND_VTK_PARAM(20,long); // WindowIndex
        SEND_VTK_PARAM(21,long); // NbColorTicks
        SEND_VTK_PARAM(22,long); // NbColors

        theplot.SendEndArgPlot();
        map<const Mesh *,long> mapth;
        map<const Mesh3 *,long> mapth3;
        map<const MeshS *,long> mapthS;
        map<const MeshL *,long> mapthL;
        
        long kth=0,kth3=0,kthS=0,kthL=0;
        //  send all the mesh:
        for (size_t ii=0;ii<ll.size();ii++) {
          int i=ll[ii].i;
          long what = ll[i].what;
          const Mesh *th=0;
          const Mesh3 *th3=0;
          const MeshS *thS=0;
          const MeshL *thL=0;
          // Prepare for the sending mesh 2d
          if(what==0)
            th= ll[ii].th();
          // Prepare for the sending mesh 3d with differenciation 3D line / surface / volumuic / volum+surfac / line+volum+surfac
          if( what==5 )
              th3= & (l[i].evalm3(0,s)); // 3d mesh3 -> if contains a meshS or meshL pointer, send only the principal mesh: the mesh3
          if( what==50 )
            thS= (&(l[i].evalmS(0,s))); // 3d meshS
          if( what ==55 )
            thL= (&(l[i].evalmL(0,s))); // 3d  meshL
          // Prepare for the sending 2d iso values for ffglut
          else if (what==1 || what==2|| what==11 || what==12) {
            ll[ii].eval(fe,cmp);
            if (fe[0]->x()) th=&fe[0]->Vh->Th;
              if(fe[1] && fe[1]->x())
                    ffassert(th == &fe[1]->Vh->Th);
          }
          // Prepare for the sending 3D volume iso values for ffglut
          else if (what==6 || what==7 || what==16 || what==17) {
            ll[ii].eval(fe3,cmp);
            if (fe3[0]->x()) th3=&fe3[0]->Vh->Th;
            if (fe3[1]) ffassert(th3 == &fe3[1]->Vh->Th);
            if (fe3[2]) ffassert(th3 == &fe3[2]->Vh->Th);

          }
          // Prepare for the sending 3D surface iso values for ffglut
          else if (what==8 || what==9 || what==18 || what==19) {
            ll[ii].eval(feS,cmp);
            if (feS[0]->x()) thS=&feS[0]->Vh->Th;
            if (feS[1] && feS[1]->x()) ffassert(thS == &feS[1]->Vh->Th);
            if (feS[2] && feS[2]->x()) ffassert(thS == &feS[2]->Vh->Th);
          }
          // Prepare for the sending 3D curve iso values for ffglut
          else if (what==14 || what==15 || what==20 || what==21) {
            ll[ii].eval(feL,cmp);
            if (feL[0]->x()) thL=&feL[0]->Vh->Th;
            if (feL[1] && feL[1]->x()) ffassert(thL == &feL[1]->Vh->Th);
            if (feL[2] && feL[2]->x()) ffassert(thL == &feL[2]->Vh->Th);
          }
          // test on the type meshes --- 2D, 3D volume and 3D surface
          if(th && mapth.find(th)==mapth.end())
            mapth[th]=++kth;
          if(th3 && (mapth3.find(th3)==mapth3.end()))
            mapth3[th3]=++kth3;
          if(thS && (mapthS.find(thS)==mapthS.end()))
            mapthS[thS]=++kthS;
          if(thL && (mapthL.find(thL)==mapthL.end()))
            mapthL[thL]=++kthL;
        }
        // send of meshes 2d
        theplot.SendMeshes();
        theplot << kth ;

        for (map<const Mesh *,long>::const_iterator i=mapth.begin();i != mapth.end(); ++i)
            theplot << i->second << *  i->first ;

       // only send of volume meshes 3D if completed mesh3 (meshS!=NULL or/and !=meshL!=NULL)
        if(kth3) {
            theplot.SendMeshes3();
            theplot << kth3 ;
            for (map<const Mesh3 *,long>::const_iterator i=mapth3.begin();i != mapth3.end(); ++i)
                theplot << i->second << *  i->first ;
           }
       if(kthS) {
            theplot.SendMeshesS();
            theplot << kthS ;
            for (map<const MeshS *,long>::const_iterator i=mapthS.begin();i != mapthS.end(); ++i)
                theplot << i->second << *  i->first ;
       }
        if(kthL) {
          theplot.SendMeshesL();
            theplot << kthL ;
          for (map<const MeshL *,long>::const_iterator i=mapthL.begin();i != mapthL.end(); ++i)
            theplot << i->second << *  i->first ;
        }

        // end of what ploting for meshes
        theplot.SendPlots();


        theplot <<(long) ll.size();
        for (size_t ii=0;ii<ll.size();ii++)
        {
            int i=ll[ii].i;
            long what = ll[ii].what;
            int err =1;//  by default we are in error
            const Mesh *pTh=0;
            const Mesh3 *pTh3=0;
            const MeshS *pThS=0;
            const MeshL *pThL=0;

            // send 2d meshes for ffglut
            if(what ==0)
            {
                pTh=ll[ii].th();
                if(pTh) {
                    err=0;
                    theplot << what ;
                    theplot <<mapth[ ll[ii].th() ];// numero du maillage 2d
                }
            }
            // send 2d iso values for ffglut
            else if (what==1 || what==2 )
                err = Send2d<R,v_fes>( theplot,ll[ii] ,mapth);
            else if (what==11 || what==12 )
                err = Send2d<Complex,v_fes>( theplot,ll[ii] ,mapth);

            else if (what==3 || what==13 )
            {

                what=13;
                theplot << what  ; //
                KN<double> z0;
                if(verbosity>99)
                    cout << " sendplot curve "<< what << " " << ii ;

                for(int k=0; k<4; ++k)
                {
                    int ilat=ll[ii].l[k];
                    if (ilat>=0)
                    {
                        KN_<double> t = GetAny<KN_<double> >(lat[ilat]);
                        theplot << t;
                        if(verbosity>99) cout << " (" << k <<" " << ilat << ") " << t.N();
                    }
                    else theplot << z0;// empty arry ...
                }
                if(verbosity>99)
                    cout <<endl;
                err=0;
            }

            else if (l[i].what==4 ) {
                err=0;
                theplot << what ;
                const  E_BorderN * Bh= l[i].evalb(0,s);
                Bh->SavePlot(s,theplot);
            }
            // send volume 3d meshes for ffglut
            else if(what ==5) {
                pTh3=&l[i].evalm3(0,s);
                if(pTh3) {
                    err=0;
                    theplot << what ;
                    theplot <<mapth3[ &l[i].evalm3(0,s)];// numero du maillage 3D volume
                }
            }
            else if (l[i].what==40 ) {
                err=0;
                theplot << what ;
                const  E_Curve3N * Bh3= l[i].evalc(0,s);
                Bh3->SavePlot(s,theplot);
            }
            else if(what ==50) {
                pThS=&(l[i].evalmS(0,s));
                    if(pThS) {
                        err=0;
                        theplot << what ;
                        theplot <<mapthS[ &(l[i].evalmS(0,s))];// numero du maillage 3D surface
                 }
            }
            else if(what ==55) {
                pThL=&(l[i].evalmL(0,s));
                if(pThL) {
                    err=0;
                    theplot << what ;
                    theplot <<mapthL[ &(l[i].evalmL(0,s))];// numero du maillage 3D line
                }
            }

            // send 3D volume iso values for ffglut
            else  if (what==6 ||what==7 )
                err = Send3d<R,v_fes3>( theplot,ll[ii] ,mapth3);
            else if (what==16 || what==17 )
                err = Send3d<Complex,v_fes3>( theplot,ll[ii] ,mapth3);
            // send 3D surface iso values for ffglut
            else  if (what==8 ||what==9 )
                err = SendS<R,v_fesS>( theplot,ll[ii] ,mapthS);
            else if (what==18 || what==19 )
                err = SendS<Complex,v_fesS>( theplot,ll[ii] ,mapthS);
            // send 3D curve iso values for ffglut
            else  if (what==14 ||what==15 )
                err = SendL<R,v_fesL>( theplot,ll[ii] ,mapthL);
            else if (what==20 || what==21 )
                err = SendL<Complex,v_fesL>( theplot,ll[ii] ,mapthL);
            else
                ffassert(0);// erreur type theplot inconnue
            if(err==1)
            { if(verbosity)
                cerr << "Warning: May be a bug in your script, \n"
                << " a part of the plot is wrong t (mesh or FE function, curve)  => skip the item  " << i+1
                << " in plot command " << endl;
                theplot << -1L << (long) i ;
            }


        }
        theplot.SendEndPlot();
    }

    // begin to the post scrit procedure
    if (!withrgraphique) {initgraphique();withrgraphique=true;}
    viderbuff();
    MeshPoint *mps=MeshPointStack(s),mp=*mps ;
    int nbcolors=0;
    float *colors=0;
    bool hsv=true; // hsv  type
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
    double ArrowSize=-1;
    // PPPPP
    KN<R> Viso,Varrow;

    bool bw=false;
    string * psfile=0;
    string * cm=0;
    pferbase  fe=0,fe1=0;
    int cmp0,cmp1;
    bool grey=getgrey();
    bool greyo=grey;
    bool drawborder=true;
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
    if (nargs[15]) {
        KN_<double> cc= GetAny<KN_<double> >((*nargs[15])(s));
        nbcolors= cc.N()/3;
        if ( nbcolors > 1&& nbcolors < 100)
        {
            colors = new float [nbcolors*3];
            for (int i=0; i<3*nbcolors; i++) colors[i]=cc[i];
        }
        else nbcolors = 0;
    }
    if (nargs[16]) drawborder= GetAny<bool>((*nargs[16])(s));
    int   dimplot=2;
    if (nargs[17]) dimplot= GetAny<long>((*nargs[17])(s));
    bool addtoplot=false, keepPV=false;
    if (nargs[18]) addtoplot= GetAny<bool>((*nargs[18])(s));
    if (nargs[19]) keepPV= GetAny<bool>((*nargs[19])(s));
    if (nargs[VTK_START+8]) ArrowSize = GetAny<double>((*nargs[VTK_START+8])(s));
    //  for the gestion of the PTR.
    WhereStackOfPtr2Free(s)=new StackOfPtr2Free(s);// FH aout 2007

    setgrey(grey);
    if (Viso.unset()) Viso.init(Niso);
    if (Varrow.unset()) Varrow.init(Narrow);

    const Mesh * cTh=0;
    bool vecvalue=false,isovalue=false;
    bool ops=psfile;
    bool drawmeshes=false;
    if ( clean ) {

        // ALH - 28/3/15 - Open PS file before blanking the current picture because Javascript needs to know any "ps="
        // parameter to send the graphical commands to the right canvas.

        if (psfile) {
            // [[file:../Graphics/sansrgraph.cpp::openPS]]
            openPS(psfile->c_str());
        }

        reffecran();

        if (bw) NoirEtBlanc(1);
        R2 Pmin,Pmax;
        R2 uminmax(1e100,-1e100);
        R2 Vminmax(1e100,-1e100);
        bool first=true;
        for (size_t ii=0;ii<ll.size();ii++)
        {
            int i=ll[ii].i;
            long what = ll[ii].what;
            R2  P1,P2;
            R3  P11,P22;
            if (what==1 || what==2)
            {

                if( !uaspectratio)   aspectratio= true;
                ll[ii].eval(fe,cmp0,fe1,cmp1);

                if (!fe->x()) continue;

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
                            Vminmax = minmax(Vminmax,fe->Vh->MinMax(u,v,cmp0,cmp1));
                        }
                        else
                            cerr << " On ne sait tracer que de vecteur sur un meme type element finite. " << endl;
                    }

                }
            }
            else if (l[i].what==0)
            {
                if( !uaspectratio) aspectratio= true;
                const  Mesh & Th= *ll[ii].th();
                Th.BoundingBox(P1,P2);
                cTh=&Th;
            }
            else if (l[i].what==4)
            {
                if( !uaspectratio) aspectratio= true;
                const  E_BorderN * Bh= l[i].evalb(0,s);
                Bh->BoundingBox(s,P1.x,P2.x,P1.y,P2.y);

            }
            else if (l[i].what==40)
            {
                if( !uaspectratio) aspectratio= true;
                const  E_Curve3N * Bh3= l[i].evalc(0,s);
                Bh3->BoundingBox(s,P11.x,P22.x,P11.y,P22.y,P11.z,P22.z);
                
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
                cout << " u bound " <<  uminmax << "  V : " << Vminmax <<endl;

            if (verbosity>1)
                cout << "Plot bound [x,y] " <<  Pmin << " max [x,y] " << Pmax << endl;
            int N=Viso.N();
            int Na=Varrow.N();
            R2 O((Pmin+Pmax)/2);
            R rx(Pmax.x-Pmin.x),ry(Pmax.y-Pmin.y);
            // bug   version 1.41 correct FH to remove div by zero.
            rx = Max(rx,1e-30);
            ry = Max(ry,1e-30);
            // -- end correction
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
            x=0; d= sqrt(Vminmax.y)/(Na-1.001);
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
        bool thfill=fill;
        for (size_t ii=0;ii<ll.size();ii++)
        {
            int i=ll[ii].i;
            long what = ll[i].what;

            if (l[i].what==0)
                if (fill)
                    ll[ii].th()->Draw(0,thfill);
                else
                    ll[ii].th()->Draw(0,thfill);
                else  if (what==1 || what==2)
                {

                    ll[ii].eval(fe,cmp0,fe1,cmp1);
                    if (!fe->x()) continue;
#ifdef VVVVVVV
                    cout << "   Min = " << fe->x->min() << " max = " << fe->x->max() ;
                    if(fe1 && verbosity > 1)
                        cout << " Min = " << fe1->x->min() << " max = " << fe1->x->max() ;
                    cout << endl;
#endif
                    if (fe1)
                    {
                        if (fe->Vh == fe1->Vh)
                            vecvalue=true,fe->Vh->Draw(*fe->x(),*fe1->x(),Varrow,coeff,cmp0,cmp1,colors,nbcolors,hsv,drawborder,ArrowSize);
                        else
                            cerr << " Draw only vector field on same Finites Element , Sorry. " << endl;
                        if (drawmeshes) fe->Vh->Th.Draw(0,fill);
                    }
                    else


                        if (fill)
                            isovalue=true,fe->Vh->Drawfill(*fe->x(),Viso,cmp0,1.,colors,nbcolors,hsv,drawborder);
                        else
                            isovalue=true,fe->Vh->Draw(*fe->x(),Viso,cmp0,colors,nbcolors,hsv,drawborder);

                    if (drawmeshes) fe->Vh->Th.Draw(0,fill);

                }
                else if (l[i].what==4) {
                    const  E_BorderN * Bh= l[i].evalb(0,s);
                    Bh->Plot(s);
                }
                else if(l[i].what==3)
                {
                    penthickness(6);
                    tab x=l[i].evalt(0,s);
                    tab y=l[i].evalt(1,s);
                    KN<double> pz0;
                    KN_<double> z(pz0), v(pz0);
                    if (l[i].e[2]) { z.set(l[i].evalt(2,s));}
                    if (l[i].e[3]) { v.set(l[i].evalt(3,s));}
                    long k= Min(x.N(),y.N());
                    NewSetColorTable(Viso.N()+4,colors,nbcolors,hsv);
                    rmoveto(x[0],y[0]);
                    couleur(2+i);
                    for (int i= 1;i<k;i++)
                        rlineto(x[i],y[i]);
                }
                else {
                    if(verbosity) cout << "  Plot::  Sorry no ps version for this type of plot "
                         << l[i].what <<endl;
                }
            thfill=false;
        }
        if (value) {
            int k=0;
            if (isovalue) {PlotValue(Viso,k,"IsoValue");k+= Viso.N()+3;}
            if (vecvalue) {PlotValue(Varrow,k,"Vec Value");k+= Varrow.N()+3;}
        }

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
                case 'g' : setgrey(grey= !getgrey()); plotting=true;
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
#ifdef DRAWING
                    if (cTh) cTh->quadtree->Draw();
#endif
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
                    Show("Enter a keyboard character in the FreeFem Graphics window in order to:",i++);

                    i+=2;
                    Show("+)  zoom in around the cursor 3/2 times ",i++);
                    Show("-)  zoom out around the cursor 3/2 times  ",i++);
                    Show("=)  reset zooming  ",i++);
                    Show("r)  refresh plot ",i++);
                    Show("ac) increase   the size arrow ",i++);
                    Show("AC) decrease the size arrow  ",i++);
                    Show("b)  switch between black and white or color plotting ",i++);
                    Show("g)  switch between grey or color plotting ",i++);
                    Show("f)  switch between filling iso or not  ",i++);
                    Show("v)  switch between show  the numerical value of iso or not",i++);
                    Show("p)   save  plot in a Postscprit file",i++);
                    Show("m)  switch between show  meshes or not",i++);
                    Show("p)  switch between show  quadtree or not (for debuging)",i++);
                    Show("t)  find  Triangle ",i++);
                    Show("?)  show this help window",i++);
                    Show("any other key : continue ",++i);
                    goto next;
            }
            if (!pViso || !pVarrow)
            { //  recompute the iso bound
                R2 uminmax(1e100,-1e100);
                R2 Vminmax(1e100,-1e100);
                for (size_t i=0;i<l.size();i++)
                { R2  P1,P2;
                    if (l[i].what==1 || l[i].what==2)
                    {
                        fe   = l[i].eval(0,s,cmp0);
                        fe1  = l[i].eval(1,s,cmp1);

                        if (!fe->x()) continue;


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
    if (colors) delete[] colors;
    viderbuff();

    return 0L;
}

AnyType Convect::operator()(Stack s) const
{
    if(d==2)  return eval2(s);
    else  return eval3(s);
}

AnyType Convect::eval2(Stack s) const
{
  MeshPoint* mp(MeshPointStack(s)),mpc(*mp);
  MeshPointStack(s,&mpc);// P  ptr on  variable mpc ...

  static MeshPoint mpp,mps;
  static int stateold=0;
  static int count =0;
  static R ddtp=0;
  R ddt = GetAny<double>((*dt)(s));
  if (ddt)
    {
      if( (stateold == state) && (ddt==ddtp) && (*mp==mpp  ) )// optim same convect at same point nov/2015
      {
          if( verbosity > 3 && count++ < 10)
              cout <<" -- optim convect " <<stateold << "  P= " << mp->P  << ", "<< mp->T
              << " chi(P) = " << mps.P <<  ","<< mps.T <<endl;
	  mpc=mps;
      }
      else
	{

          if( verbosity > 3 && count++ < 10*10) cout <<" -- no optim convect " <<stateold << " " << state
                    << " P= " << mp->P  << " PP= " << mpp.P <<  endl;
          stateold=state;// correction FH..
          ddtp=ddt;
	  const Mesh & Th(*mp->Th);
	  ffassert(mp->Th && mp->T);
	  R l[3];
	  l[1]=mpc.PHat.x;
	  l[2]=mpc.PHat.y;
	  l[0]=1-l[1]-l[2];

	  int k=0;
	  int j;
	  int it=Th(mpc.T);
	  while ( (j=WalkInTriangle(Th,it,l,GetAny<double>((*u)(s)),GetAny<double>((*v)(s)),ddt))>=0)
	    {
	      ffassert( l[j] == 0);
	      //int jj  = j;
	      R a= l[(j+1)%3], b= l[(j+2)%3];
	      int itt =  Th.ElementAdj(it,j);
	      if(itt==it || itt <0)  break; // le bord
	      it = itt;
	      l[j]=0;
	      l[(j+1)%3] = b;
	      l[(j+2)%3] = a;
	      mpc.change(R2(l[1],l[2]),Th[it],0);
	      if(k++>1000)
		{
		  cerr << "Fatal  error  in Convect (R2) operator: loop  => velocity too high ???? or NaN F. Hecht  " << endl;
		  ffassert(0);
		}
	    }

	  mpc.change(R2(l[1],l[2]),Th[it],0);
	  mpp=*mp;// previous value
	  mps=mpc;// convect value
	}
    }
  // warning use poit on &mpc .. bug correct in dec 2015 F.H.
  AnyType r= (*ff)(s);
  if( verbosity > 3 && count++ < 10*10) cout << "  %%%r= "<< GetAny<double>(r) << "  P= " << mp->P  << ", "<< mp->T << endl;
   MeshPointStack(s,mp);// restor old pointeur ..

  return r;
}

inline int FindandAdd(set<int> *st,vector<int> & lst,int k)
{
    if(lst.size() < 10)
    {

    }
    return 0;
}



AnyType Convect::eval3old(Stack s) const
{
    extern long newconvect3;
    if(newconvect3) return eval3(s);//  New Convect in test
    MeshPoint* mp(MeshPointStack(s)),mpc(*mp);
    MeshPointStack(s,&mpc);// P  ptr on  variable mpc ...

     static MeshPoint mpp,mps;// previous state ..
    static int stateold=0;
    static int count =0;
    static R ddtp=0;

    randwalk(-1); // init randwalk

    R ddt = GetAny<double>((*dt)(s));
    if (ddt)
    {
        bool ddd=verbosity>1000;
        if( (stateold == state) && (ddt==ddtp) && (*mp==mpp  ) )// optim same convect at same point nov/2015
        {
            if( verbosity > 3 && count++ < 10)
                cout <<" -- optim convect3 " <<stateold << endl;
            mpc=mps;
        }
        else
        {

           if( verbosity > 3 && count++ < 10*10) cout <<" -- no optim3 convect " <<stateold << " " << state
                << " P= " << mp->P  << " PP= " << mpp.P <<  endl;
           const Mesh3 & Th3(*mp->Th3);
            ffassert(mp->Th3 && mp->T3);
            R3 PHat=mpc.PHat;

            int k=0;
            int j;
            int it=Th3(mpc.T3);
            if(ddd) cout << " IN: " <<  (*mpc.T3)(PHat) << " ; " << mpc.P <<" : " << ddt << endl;
            while ( (j=WalkInTet(Th3,it,PHat,R3(GetAny<double>((*u)(s)),GetAny<double>((*v)(s)),GetAny<double>((*w)(s))),ddt))>=0)
                if(j>3)  {
                    it=j-4;
                    mpc.change(PHat,Th3[it],0);
                    if(ddd) cout << "   **P= "<< (*mpc.T3)(PHat) << " ,  Ph " << PHat << " : j = " << j  <<  " it:  " << it << "ddt=" << ddt ;
                    if(ddt==0) break; // finish ...

                }
                else
                {
                    if(ddd) cout << "P= "<< (*mpc.T3)(PHat) << " ,  Ph " << PHat << " : j = " << j  <<  " it:  " << it ;
#ifdef DEBUG
                    R3  Po=(*mpc.T3)(PHat),Pho=PHat; int ito=it;
#endif
                    int itt =  Th3.ElementAdj(it,j,PHat);
                    if(ddd && itt>=0) cout << "  -> " << itt << " " << j  << "  : Pn " <<  Th3[itt](PHat) << " PHn " << PHat << " , " << ddt << endl;
                    if(itt<0) break;
                    it=itt;
                    mpc.change(PHat,Th3[it],0);
#ifdef DEBUG
                    if(((Po-mpc.P).norme2() > 1e-10))
                    {   cout << ito << " " << &Th3[ito][0] << " " << &Th3[ito][1] << " " << &Th3[ito][2] << " " << &Th3[ito][3] << " "  << endl;
                        cout << it << " " <<  &Th3[it][0] << " " <<  &Th3[it][1] << " " <<  &Th3[it][2] << " " <<  &Th3[it][3] << endl;
                        cout << Pho  << "o Hat " << PHat << endl;
                        cout << Po << " o != " << mpc.P << " diff= "<< (Po-mpc.P).norme2() <<endl;
                        assert(0);

                    }
#endif
                    ffassert(k++<2000);
                }

            mpc.change(PHat,Th3[it],0);
            mpp=*mp;
            mps=mpc;
        }
    }
    AnyType r= (*ff)(s);
    MeshPointStack(s,mp);
    return r;
}

AnyType Convect::eval3(Stack s) const
{  // nouvelle version de convect 3d Feb 2015  version 3.44-01
    MeshPoint* mp(MeshPointStack(s)),mpc(*mp);
    MeshPointStack(s,&mpc);// P  ptr on  variable mpc ...
    static MeshPoint mpp,mps;
    static R ddts;
    randwalk(-1); // init randwalk
    R3 offset;
    R ddt = GetAny<double>((*dt)(s));
    if (ddt)
    {
        bool ddd=verbosity>1000;

        if(*mp==mpp && ddt == ddts)
            mpc=mps;
        else
        {

            const Mesh3 & Th3(*mp->Th3);
            ffassert(mp->Th3 && mp->T3);
            R3 PHat=mpc.PHat;

            int k=0;
            int j;
            int it=Th3(mpc.T3);
            if(ddd) cout << " IN: " <<  (*mpc.T3)(PHat) << " ; " << mpc.P <<" : " << ddt << endl;
            while ( (j=WalkInTetn(Th3,it,PHat,R3(GetAny<double>((*u)(s)),GetAny<double>((*v)(s)),GetAny<double>((*w)(s))),ddt,offset))>=0)
                if(j>3)  {
                    it=j-4;
                    mpc.change(PHat,Th3[it],0);
                    if(ddd) cout << "   **P= "<< (*mpc.T3)(PHat) << " ,  Ph " << PHat << " : j = " << j  <<  " it:  " << it << "ddt=" << ddt ;
                    if(ddt==0) break; // finish ...

                }
                else
                {
                    if(ddd) cout << "P= "<< (*mpc.T3)(PHat) << " ,  Ph " << PHat << " : j = " << j  <<  " it:  " << it ;
#ifdef DEBUG
                    R3  Po=(*mpc.T3)(PHat),Pho=PHat; int ito=it;
#endif
                    int itt =  Th3.ElementAdj(it,j,PHat);
                    if(ddd && itt>=0) cout << "  -> " << itt << " " << j  << "  : Pn " <<  Th3[itt](PHat) << " PHn " << PHat << " , " << ddt << endl;
                    if(itt<0) break;
                    it=itt;
                    mpc.change(PHat,Th3[it],0);
#ifdef DEBUG
                    if(((Po-mpc.P).norme2() > 1e-10))
                    {   cout << ito << " " << &Th3[ito][0] << " " << &Th3[ito][1] << " " << &Th3[ito][2] << " " << &Th3[ito][3] << " "  << endl;
                        cout << it << " " <<  &Th3[it][0] << " " <<  &Th3[it][1] << " " <<  &Th3[it][2] << " " <<  &Th3[it][3] << endl;
                        cout << Pho  << "o Hat " << PHat << endl;
                        cout << Po << " o != " << mpc.P << " diff= "<< (Po-mpc.P).norme2() <<endl;
                        assert(0);

                    }
#endif
                    ffassert(k++<2000);
                }

            mpc.change(PHat,Th3[it],0);
            mpp=*mp;
            mps=mpc;
        }
    }
    ddts=ddt;
    AnyType r= (*ff)(s);
    MeshPointStack(s,mp);
    return r;
}


template<class K>
class Op3_pfe2K : public ternary_function<pair<FEbase<K,v_fes> *,int>,R,R,K> { public:


  class Op : public E_F0mps { public:
    Expression a,b,c;
    Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {/*cout << "Op3_pfe2K" << endl;*/}
    AnyType operator()(Stack s)  const
    {
      R xx(GetAny<R>((*b)(s)));
      R yy(GetAny<R>((*c)(s)));
      MeshPoint & mp = *MeshPointStack(s),mps=mp;
      mp.set(xx,yy,0.0);
      AnyType ret = pfer2R<K,0>(s,(*a)(s));
      mp=mps;
      return  ret;
    }
  };
};


//Add  FH 16032005
class Op3_Mesh2mp : public ternary_function<pmesh*,R,R,MeshPoint *> { public:
  class Op : public E_F0mps { public:
      Expression a,b,c;
       Op(Expression aa,Expression bb,Expression cc) : a(aa),b(bb),c(cc) {}
       AnyType operator()(Stack s)  const
        {
           R xx(GetAny<R>((*b)(s)));
           R yy(GetAny<R>((*c)(s)));
           pmesh *ppTh(GetAny<pmesh*>((*a)(s)));
           if( !ppTh || !*ppTh) ExecError("Op3_Mesh2mp unset mesh ??");
           pmesh pTh(*ppTh);
           MeshPoint * mp = new MeshPoint();
           mp->set(xx,yy,0.0);
           R2 PHat;
           bool outside;
           const Triangle * K=pTh->Find(mp->P.p2(),PHat,outside);
           mp->set(*pTh,(R2) mp->P.p2(),PHat,*K,0,outside);
           return mp;}

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
           rd = 0. ;// to be compatible with varf def... FH april 2014 ... v 3.31.
           if ( mp->SetAdj() )
             rd = GetAny<A>((*a)(stack));
           *mp=smp;
           return  SetAny<R>(rd-rg);//  external - internal
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
           return  SetAny<R>((rg+rd)*0.5);
        }
       MeanOp(Expression aa) : a(aa) {}
    };

template<class RR,class AA=RR>
class OthersideOp: public E_F0mps  { public:
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
        rd = 0.;
        if ( mp->SetAdj() )
            rd = GetAny<A>((*a)(stack));
        *mp=smp;
        return  SetAny<R>(rd);
    }
    OthersideOp(Expression aa) : a(aa) {}
};

long get_size(pferarray const & a)
{
  return a.first->N;
}
long get_size(pfecarray const & a)
{
    return a.first->N;
}
long get_size(pferbasearray *const & a)
{
    return (**a).N;
}
long get_size(pfecbasearray *const & a)
{
    return (**a).N;
}
long get_size(pf3rarray const & a)
{
  return a.first->N;
}
long get_size(pf3rbasearray *const & a)
{
    return (**a).N;
}
long get_size(pf3carray const & a)
{
    return a.first->N;
}
long get_size(pf3cbasearray *const & a)
{
    return (**a).N;
}
long resize(pferbasearray *const & a, long const & n)
{
    (**a).resize(n);
    return n;
}

long resize(pfecbasearray *const & a, long const & n)
{
    (**a).resize(n);
    return n;
}

pferbase* get_element(pferbasearray *const & a, long const & n)
{
  return (**a)[n];
}

pfer get_element(pferarray const & a, long const & n)
{
  return pfer( *(*a.first)[n],a.second);
}

//  complex case
pfecbase* get_element(pfecbasearray *const & a, long const & n)
{
  return (**a)[n];
}
pfec get_element(pfecarray const & a, long const & n)
{
  return pfec( *(*a.first)[n],a.second);
}
//  end complex case

lgElement get_element(pmesh const & a, long const & n)
{
  return lgElement(a,n);
}
lgElement get_element(pmesh *const & a, long const & n)
{
  return lgElement(*a,n);
}



lgBoundaryEdge get_belement(lgBoundaryEdge::BE const & a, long const & n)
{
    return lgBoundaryEdge(a,n);
}

lgElement get_adj(lgElement::Adj const & a, long  * const & n)
{
    return  a.adj(*n);
}

lgVertex get_vertex(pmesh const & a, long const & n)
{
    return lgVertex(a,n);
}
lgVertex get_vertex(pmesh *const & a, long const & n)
{
    return lgVertex(*a,n);
}
lgVertex get_element(lgElement const & a, long const & n)
{
  return a[n];
}
lgVertex get_belement(lgBoundaryEdge const & a, long const & n)
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
long getlab(lgElement const & a)
{
  return a.lab();
}
long getlab(lgBoundaryEdge const & a)
{
    return a.lab();
}

double getarea(lgElement const & a)
{
    return a.area();
}
double getlength(lgBoundaryEdge const & a)
{
    return a.length();
}
lgElement getElement(lgBoundaryEdge const & a)
{
    return a.Element();
}
long EdgeElement(lgBoundaryEdge const & a)
{
    return a.EdgeElement();
}

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

template<class A> inline AnyType DestroyKNmat(Stack,const AnyType &x){
  KN<A> * a=GetAny<KN<A>*>(x);
  for (int i=0;i<a->N(); i++)
    (*a)[i].destroy();
  a->destroy();
  return  Nothing;
}

template<class R>  R * set_initmat( R* const & a,const long & n){
 SHOWVERB( cout << " set_init " << typeid(R).name() << " " << n << endl);
  a->init(n);
  for (int i=0;i<n;i++)
    (*a)[i].init();
   return a;}

void init_mesh_array()
 {
  Dcl_Type<KN<pmesh> *>(0,::DestroyKN<pmesh> );
  atype<KN<pmesh>* >()->Add("[","",new OneOperator2_<pmesh*,KN<pmesh>*,long >(get_elementp_<pmesh,KN<pmesh>*,long>));
    TheOperators->Add("<-",
       new OneOperator2_<KN<pmesh> *,KN<pmesh> *,long>(&set_initinit));
  map_type_of_map[make_pair(atype<long>(),atype<pmesh>())]=atype<KN<pmesh>*>(); // vector

  // resize mars 2006 v2.4-1
  Dcl_Type< Resize<KN<pmesh> > > ();
  Add<KN<pmesh>*>("resize",".",new OneOperator1< Resize<KN<pmesh> >,KN<pmesh>*>(to_Resize));
  Add<Resize<KN<pmesh> > >("(","",new OneOperator2_<KN<pmesh> *,Resize<KN<pmesh> > , long   >(resizeandclean1));


 }
template<class RR,class A,class B>
RR  get_elementp(const A & a,const B & b){
  if( b<0 || a->N() <= b)
   { cerr << " Out of bound  0 <=" << b << " < "  << a->N() << " array type = " << typeid(A).name() << endl;
     ExecError("Out of bound in operator []");}
    return  ((*a)[b]);}

template<class T> T *resizeandclean2(const Resize<T> & t,const long &n)
 {  // resizeandclean1
  int nn= t.v->N(); // old size

  for (int i=n;i<nn;i++)  {(*t.v)[i].destroy();;} // clean
  t.v->resize(n);
  for (int i=nn;i<n;i++)  {(*t.v)[i].init();}
  return  t.v;
 }

template<class PMat>
AnyType ClearReturn(Stack stack, const AnyType & a)
{
    // a ne faire que pour les variables local au return...
    //  pour l'instant on copie pour fqire mqrche
    // a repense  FH  mqi 20014....
    PMat * m = GetAny<PMat * >(a);
    m->increment();
    Add2StackOfPtr2FreeRC(stack,m);
    return m;
}


template <class R>
void DclTypeMatrix()
{
   Dcl_Type<RNM_VirtualMatrix<R>*>( ); // ???????  ZZZZZZ

  Dcl_Type<Matrice_Creuse<R>* >(InitP<Matrice_Creuse<R> >,Destroy<Matrice_Creuse<R> >, ClearReturn<Matrice_Creuse<R> >);
    // newpMatrice_Creuse
    Dcl_Type<newpMatrice_Creuse<R> >();      // to def new Matrice_Creuse
    Dcl_Type<Matrice_Creuse_Transpose<R> >();      // matrice^t   (A')

  Dcl_Type<Matrice_Creuse_Transpose<R> >();      // matrice^t   (A')
  Dcl_Type<Matrice_Creuse_inv<R> >();      // matrice^-1   A^{-1}
  Dcl_Type<Matrice_Creuse_inv_trans<R> >();      // matrice^-1   A'^{-1}
  Dcl_Type<typename RNM_VirtualMatrix<R>::plusAx >();       // A*x (A'*x)
  Dcl_Type<typename RNM_VirtualMatrix<R>::plusAtx >();       // A^t*x (A'*x)
  Dcl_Type<typename RNM_VirtualMatrix<R>::solveAxeqb >();       // A^-1*x (
    Dcl_Type<typename RNM_VirtualMatrix<R>::solveAtxeqb >();       // A^t^-1*x
  Dcl_Type<Matrix_Prod<R,R> >();
  Dcl_Type<list<tuple<R,MatriceCreuse<R> *,bool> >*>();

   // resize mars 2006 v2.4-1

  // array of matrix   matrix[int] A(n);
  typedef  Matrice_Creuse<R> Mat;
  typedef  Mat * PMat;
  typedef  KN<Mat> AMat;
  Dcl_Type<AMat *>(0,::DestroyKNmat<Mat> );

  // init array
    TheOperators->Add("<-",
       new OneOperator2_<AMat *,AMat *,long>(&set_initmat)
                      );
  //  A[i]
  atype<AMat* >()->Add("[","",new OneOperator2_<PMat,AMat*,long >(get_elementp_<Mat,AMat*,long>));

  // resize
  Dcl_Type< Resize<AMat> > ();
  Add<AMat*>("resize",".",new OneOperator1< Resize<AMat >,AMat*>(to_Resize));
  Add<Resize<AMat> >("(","",new OneOperator2_<AMat *,Resize<AMat> , long   >(resizeandclean2));

  // to declare matrix[int]
  map_type_of_map[make_pair(atype<long>(),atype<PMat>())]=atype<AMat*>();

}


template<class A,class B>
AnyType First(Stack,const AnyType &b) {
return   SetAny<A>(GetAny<B>(b).first);}

template<class K>
AnyType AddIncrement(Stack stack, const AnyType & a)
{
    K m = GetAny<K>(a);
    m->increment();
    Add2StackOfPtr2FreeRC(stack,m);
    if(verbosity>1)
    cout << "AddIncrement:: increment + Add2StackOfPtr2FreeRC " << endl;
    return a;
}

// FE 3D volume
Type_Expr CConstantTFE3(const EConstantTypeOfFE3::T & v)
{
    throwassert(map_type[typeid( EConstantTypeOfFE3::T).name()]);
    return make_pair(map_type[typeid( EConstantTypeOfFE3::T).name()],new EConstantTypeOfFE3(v));
}

// FE 3D surface
Type_Expr CConstantTFES(const EConstantTypeOfFES::T & v)
{
    throwassert(map_type[typeid( EConstantTypeOfFES::T).name()] !=0);
    return make_pair(map_type[typeid( EConstantTypeOfFES::T).name()],new EConstantTypeOfFES(v));
}

// FE 3D curve
Type_Expr CConstantTFEL(const EConstantTypeOfFEL::T & v)
{
    throwassert(map_type[typeid( EConstantTypeOfFEL::T).name()] !=0);
    return make_pair(map_type[typeid( EConstantTypeOfFEL::T).name()],new EConstantTypeOfFEL(v));
}

//  end --- call meth be ..
// 2013 resize of array of fe function..
template<typename  T> T fepresize(const Resize1<T> & rt,const long &n) {
    (**(rt.v)).resize(n);
    return rt.v;}
template<typename  T> T feresize(const Resize1<T> & rt,const long &n) {
    rt.v.first->resize(n);
    return rt.v;}
double get_R3(R3 *p,long i){return (*p)[i];}
R3 * set_eqp(R3 *a,R3 *b) { *a=*b; return a;}


void  init_lgfem()
{
  if(verbosity&& (mpirank==0)) cout <<"lg_fem ";
#ifdef HAVE_CADNA
  cout << "cadna ";
  cadna_init(-1); // pas de fichier
#endif

 Dcl_Type<MeshPoint*>();
 Dcl_Type<R3*>(::Initialize<R3>);
 Dcl_Type<R2*>(::Initialize<R2>);

 map_type[typeid(R3*).name()] = new ForEachType<R3*>(Initialize<R3>);
  Dcl_TypeandPtr<pmesh>(0,0, ::InitializePtr<pmesh>,::DestroyPtr<pmesh>,AddIncrement<pmesh>,NotReturnOfthisType);
  Dcl_TypeandPtr<pmesh3>(0,0,::InitializePtr<pmesh3>,::DestroyPtr<pmesh3>,AddIncrement<pmesh3>,NotReturnOfthisType);
  Dcl_TypeandPtr<pmeshS>(0,0,::InitializePtr<pmeshS>,::DestroyPtr<pmeshS>,AddIncrement<pmeshS>,NotReturnOfthisType);
  Dcl_TypeandPtr<pmeshL>(0,0,::InitializePtr<pmeshL>,::DestroyPtr<pmeshL>,AddIncrement<pmeshL>,NotReturnOfthisType);
    
  Dcl_Type<lgVertex>();
  Dcl_Type<lgElement>( );
  Dcl_Type<lgElement::Adj>( );

  Dcl_Type<lgBoundaryEdge::BE >( );
  Dcl_Type<lgBoundaryEdge>( );

  atype<long>()->AddCast(
    new E_F1_funcT<long,lgVertex>(Cast<long,lgVertex>),
    new E_F1_funcT<long,lgElement>(Cast<long,lgElement>),
    new E_F1_funcT<long,lgBoundaryEdge>(Cast<long,lgBoundaryEdge>)
  );

 Dcl_Type<TypeOfFE*>();
 Dcl_Type<TypeOfFE3*>(); // 3D volume
 Dcl_Type<TypeOfFES*>(); // 3D surface
 Dcl_Type<TypeOfFEL*>(); // 3D curve
 map_type[typeid(TypeOfFE3*).name()]->AddCast(
				 new E_F1_funcT<TypeOfFE3*,TypeOfFE*>(TypeOfFE3to2)	);
 map_type[typeid(TypeOfFES*).name()]->AddCast(
                 new E_F1_funcT<TypeOfFES*,TypeOfFE*>(TypeOfFESto2)    );
 map_type[typeid(TypeOfFEL*).name()]->AddCast(
                 new E_F1_funcT<TypeOfFEL*,TypeOfFE*>(TypeOfFELto2)    );
    
 DclTypeMatrix<R>();
 DclTypeMatrix<Complex>();

 Dcl_TypeandPtr<pferbase>(); // il faut le 2 pour pourvoir initialiser
 Dcl_TypeandPtr<pferbasearray>(); // il faut le 2 pour pourvoir initialiser
 Dcl_Type< pfer >();
 Dcl_Type< pferarray >();
 Dcl_Type< pferarray >();

//  pour des Func FE complex   // FH  v 1.43
    Dcl_TypeandPtr<pfecbase>(); // il faut le 2 pour pourvoir initialiser
    Dcl_TypeandPtr<pfecbasearray>(); // il faut le 2 pour pourvoir initialiser
    Dcl_Type< pfec >();
    Dcl_Type< pfecarray >();
    //  FH v 1.43
    // add  mai 2009 FH for 3d eigen value.
    Dcl_Type<FEbaseArrayKn<double> *>();
    Dcl_Type<FEbaseArrayKn<Complex> *>();



    // Dcl_Type< pmesharray *>(); // il faut le 2 pour pourvoir initialiser

    map_type[typeid(pfes).name()] = new ForEachType<pfes>();
    map_type[typeid(pfes*).name()] = new ForEachTypePtrfspace<pfes,2>();

// Dcl type for 3D volume FE
 Dcl_TypeandPtr<pf3rbase>(); // il faut le 2 pour pourvoir initialiser
 Dcl_TypeandPtr<pf3rbasearray>(); // il faut le 2 pour pourvoir initialiser
 Dcl_Type< pf3r >();
 Dcl_Type< pf3rarray >();

//  pour des Func FE complex   // FH  v 1.43
 Dcl_TypeandPtr<pf3cbase>(); // il faut le 2 pour pourvoir initialiser
 Dcl_TypeandPtr<pf3cbasearray>(); // il faut le 2 pour pourvoir initialiser
 Dcl_Type< pf3c >();
 Dcl_Type< pf3carray >();

 // Dcl type for 3D surface FE v 4.00
 Dcl_TypeandPtr<pfSrbase>(); // il faut le 2 pour pourvoir initialiser
 Dcl_TypeandPtr<pfSrbasearray>(); // il faut le 2 pour pourvoir initialiser
 Dcl_Type< pfSr >();
 Dcl_Type< pfSrarray >();

 //  pour des Func FE complex   // FH  v 4.00
 Dcl_TypeandPtr<pfScbase>(); // il faut le 2 pour pourvoir initialiser
 Dcl_TypeandPtr<pfScbasearray>(); // il faut le 2 pour pourvoir initialiser
 Dcl_Type< pfSc >();
 Dcl_Type< pfScarray >();

 // Dcl type for 3D curve FE v 4.5
 Dcl_TypeandPtr<pfLrbase>(); // il faut le 2 pour pourvoir initialiser
 Dcl_TypeandPtr<pfLrbasearray>(); // il faut le 2 pour pourvoir initialiser
 Dcl_Type< pfLr >();
 Dcl_Type< pfLrarray >();
    
 //  pour des Func FE complex
 Dcl_TypeandPtr<pfLcbase>(); // il faut le 2 pour pourvoir initialiser
 Dcl_TypeandPtr<pfLcbasearray>(); // il faut le 2 pour pourvoir initialiser
 Dcl_Type< pfLc >();
 Dcl_Type< pfLcarray >();

    //  cast of eigen value  mai 2009 ...
  map_type[typeid(FEbaseArrayKn<double> *).name()]->AddCast(
							      new E_F1_funcT<FEbaseArrayKn<double> *,pferbasearray>(Cast<FEbaseArrayKn<double> *,pferbasearray> ),
							      new E_F1_funcT<FEbaseArrayKn<double> *,pferarray>(First<FEbaseArrayKn<double> *,pferarray> ),
							      new E_F1_funcT<FEbaseArrayKn<double> *,pf3rbasearray>(Cast<FEbaseArrayKn<double> *,pf3rbasearray> ),
							      new E_F1_funcT<FEbaseArrayKn<double> *,pf3rarray>(First<FEbaseArrayKn<double> *,pf3rarray> )

							      );
  map_type[typeid(FEbaseArrayKn<Complex> *).name()]->AddCast(
							       new E_F1_funcT<FEbaseArrayKn<Complex> *,pfecbasearray>(Cast<FEbaseArrayKn<Complex> *,pfecbasearray> ),
							       new E_F1_funcT<FEbaseArrayKn<Complex> *,pfecarray>(First<FEbaseArrayKn<Complex> *,pfecarray> ),
							       new E_F1_funcT<FEbaseArrayKn<Complex> *,pf3cbasearray>(Cast<FEbaseArrayKn<Complex> *,pf3cbasearray> ),
							       new E_F1_funcT<FEbaseArrayKn<Complex> *,pf3carray>(First<FEbaseArrayKn<Complex> *,pf3carray> )
							       );

 map_type[typeid(pfes3).name()] = new ForEachType<pfes3>();  // 3D volume
 map_type[typeid(pfes3*).name()] = new ForEachTypePtrfspace<pfes3,3>(); // // 3D volume

 map_type[typeid(pfesS).name()] = new ForEachType<pfesS>();  // 3D surface
 map_type[typeid(pfesS*).name()] = new ForEachTypePtrfspace<pfesS,4>(); // 3D surface
    
 map_type[typeid(pfesL).name()] = new ForEachType<pfesL>();  // 3D curve
 map_type[typeid(pfesL*).name()] = new ForEachTypePtrfspace<pfesL,5>(); // 3D curve

 //
 Dcl_Type<const QuadratureFormular *>();
 Dcl_Type<const QuadratureFormular1d *>();
 Dcl_Type<const GQuadratureFormular<R3> *>();

 Global.New("qf1pT",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_1));
 Global.New("qf1pTlump",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_1lump));
 Global.New("qf2pT",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_2));
 Global.New("qf2pT4P1",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_2_4P1));
 Global.New("qf5pT",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_5));

 Global.New("qf7pT",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_7));
 Global.New("qf9pT",CConstant<const QuadratureFormular *>(&QuadratureFormular_T_9));

 Global.New("qf1pE",CConstant<const QuadratureFormular1d *>(&QF_GaussLegendre1));
 Global.New("qf2pE",CConstant<const QuadratureFormular1d *>(&QF_GaussLegendre2));
 Global.New("qf3pE",CConstant<const QuadratureFormular1d *>(&QF_GaussLegendre3));
 Global.New("qf4pE",CConstant<const QuadratureFormular1d *>(&QF_GaussLegendre4));
 Global.New("qf5pE",CConstant<const QuadratureFormular1d *>(&QF_GaussLegendre5));
 Global.New("qf1pElump",CConstant<const QuadratureFormular1d *>(&QF_LumpP1_1D));

  Global.New("qfV1",CConstant<const GQuadratureFormular<R3>  *>(&QuadratureFormular_Tet_1));
  Global.New("qfV2",CConstant<const GQuadratureFormular<R3>  *>(&QuadratureFormular_Tet_2));
  Global.New("qfV5",CConstant<const GQuadratureFormular<R3>  *>(&QuadratureFormular_Tet_5));
  Global.New("qfV1lump",CConstant<const GQuadratureFormular<R3>  *>(&QuadratureFormular_Tet_1lump));


 //  juste du code genere

 Global.New("wait",CConstant<bool*>(&TheWait));
 Global.New("NoUseOfWait",CConstant<bool*>(&NoWait));
 Global.New("NoGraphicWindow",CConstant<bool*>(&NoGraphicWindow));

 Dcl_Type<MeshPoint *>();
 Dcl_Type<finconnue *>();
 Dcl_Type<ftest *>();
 Dcl_Type<foperator *>();
 Dcl_Type<foperator *>();
 Dcl_Type<const BC_set *>();  // a set of boundary condition
 Dcl_Type<const Call_FormLinear<v_fes> *>();    //   to set Vector
 Dcl_Type<const Call_FormBilinear<v_fes> *>();  // to set Matrix
 Dcl_Type<const Call_FormLinear<v_fes3> *>();    //   to set Vector 3D volume
 Dcl_Type<const Call_FormBilinear<v_fes3> *>();  // to set Matrix 3D volume
 Dcl_Type<const Call_FormLinear<v_fesS> *>();    //   to set Vector 3D surface
 Dcl_Type<const Call_FormBilinear<v_fesS> *>();  // to set Matrix 3D surface
 Dcl_Type<const Call_FormLinear<v_fesL> *>();    //   to set Vector 3D curve
 Dcl_Type<const Call_FormBilinear<v_fesL> *>();  // to set Matrix 3D curve
 Dcl_Type<interpolate_f_X_1<double>::type>();  // to make  interpolation x=f o X^1 ;

 map_type[typeid(const FormBilinear*).name()] = new TypeFormBilinear;
 map_type[typeid(const FormLinear*).name()] = new TypeFormLinear;
 aType t_C_args = map_type[typeid(const C_args*).name()] = new TypeFormOperator;
 map_type[typeid(const Problem*).name()] = new TypeSolve<false,Problem>;
 map_type[typeid(const Solve*).name()] = new TypeSolve<true,Solve>;
 Dcl_Type<const IntFunction<double>*>();
 Dcl_Type<const IntFunction<complex<double> >*>();
 basicForEachType * t_solve=atype<const  Solve *>();
 basicForEachType * t_problem=atype<const  Problem *>();
 basicForEachType * t_fbilin=atype<const  FormBilinear *>();
 basicForEachType * t_flin=atype<const  FormLinear *>();
 basicForEachType * t_BC=atype<const BC_set *>();

 /// Doxygen doc
 basicForEachType * t_form=atype<const C_args*>();

  Dcl_Type<const CDomainOfIntegration *>();



 atype<pmesh >()->AddCast( new E_F1_funcT<pmesh,pmesh*>(UnRef<pmesh >));
 atype<pfes >()->AddCast(  new E_F1_funcT<pfes,pfes*>(UnRef<pfes>));

 atype<pferbase>()->AddCast(  new E_F1_funcT<pferbase,pferbase>(UnRef<pferbase>));
 atype<pfecbase>()->AddCast(  new E_F1_funcT<pfecbase,pfecbase>(UnRef<pfecbase>));

 Add<pfer>("[]",".",new OneOperator1<KN<double> *,pfer>(pfer2vect<R>));
 Add<pfec>("[]",".",new OneOperator1<KN<Complex> *,pfec>(pfer2vect<Complex>));


 Add<pfer>("(","",new OneTernaryOperator<Op3_pfe2K<R>,Op3_pfe2K<R>::Op> );
 Add<pfec>("(","",new OneTernaryOperator<Op3_pfe2K<Complex>,Op3_pfe2K<Complex>::Op> );
 Add<double>("(","",new OneTernaryOperator<Op3_K2R<R>,Op3_K2R<R>::Op> );
 Add<Complex>("(","",new OneTernaryOperator<Op3_K2R<Complex>,Op3_K2R<Complex>::Op> );
 Add<pmesh *>("(","",new OneTernaryOperator<Op3_Mesh2mp,Op3_Mesh2mp::Op> );


 Add<MeshPoint *>("nuTriangle",".",new OneOperator1<long,MeshPoint *>(mp_nuTriangle));
 Add<MeshPoint *>("region",".",new OneOperator1<long,MeshPoint *>(mp_region));

    Add<pfer>("refresh",".",new OneOperator1<bool,pfer>(pfer_refresh<R,v_fes>));
    Add<pfec>("refresh",".",new OneOperator1<bool,pfec>(pfer_refresh<Complex,v_fes>));

 Add<pfer>("n",".",new OneOperator1<long,pfer>(pfer_nbdf<R>));
 Add<pfec>("n",".",new OneOperator1<long,pfec>(pfer_nbdf<Complex>));
 Add<pfer>("Th",".",new OneOperator1<pmesh ,pfer>(pfer_Th<R>));
 Add<pfec>("Th",".",new OneOperator1<pmesh,pfec>(pfer_Th<Complex>));

 Add<pmesh*>("area",".",new OneOperator1<double,pmesh*>(pmesh_area));
 Add<pmesh*>("mesure",".",new OneOperator1<double,pmesh*>(pmesh_area));
 Add<pmesh*>("measure",".",new OneOperator1<double,pmesh*>(pmesh_area));
 Add<pmesh*>("bordermeasure",".",new OneOperator1<double,pmesh*>(pmesh_bordermeasure));// add june 2017 F.H
 Add<pmesh*>("nt",".",new OneOperator1<long,pmesh*>(pmesh_nt));
 Add<pmesh*>("nbe",".",new OneOperator1<long,pmesh*>(pmesh_nbe));

 Add<pmesh*>("nv",".",new OneOperator1<long,pmesh*>(pmesh_nv));

  Add<pmesh*>("hmax",".",new OneOperator1<double,pmesh*>(pmesh_hmax));
  Add<pmesh*>("hmin",".",new OneOperator1<double,pmesh*>(pmesh_hmin));

 Add<pfes*>("ndof",".",new OneOperator1<long,pfes*>(pVh_ndof));
 Add<pfes*>("Th",".",new OneOperator1<pmesh,pfes*>(pVh_Th));
 Add<pfes*>("nt",".",new OneOperator1<long,pfes*>(pVh_nt));
 Add<pfes*>("ndofK",".",new OneOperator1<long,pfes*>(pVh_ndofK));
 Add<pfes*>("(","", new OneTernaryOperator<pVh_ndf,pVh_ndf::Op>  );
/* FH: ne peux pas marcher, il faut passer aussi le nouveau Vh
 Add<pfes*>("(","", new OneBinaryOperator_st<pVh_renumber>  );
 */


 atype<Matrice_Creuse<R> * >()->AddCast(new OneOperatorCode<pb2mat<R> >);
 atype<Matrice_Creuse<Complex> * >()->AddCast(new OneOperatorCode<pb2mat<Complex> >);

//   Add all Finite Element "P0","P1","P2","RT0", ...
  for (ListOfTFE * i=ListOfTFE::all;i;i=i->next)
    {
     ffassert(i->tfe); // check
     AddNewFE(i->name,i->tfe);
     }
    static  string LU="LU";
    static  string CG="CG";
    static  string GMRES="GMRES";
    static  string Crout="CROUT";
    static  string Cholesky="Cholesky";
    static  string UMFPACK="UMFPACK";
    static  string sparsesolver="sparsesolver";
    static  string sparsesolverSym="sparsesolverSym";
    Global.New("LU",CConstant<string*>(&LU));
    Global.New(CG.c_str(),CConstant<string*>(&CG));
    Global.New(GMRES.c_str(),CConstant<string*>(&GMRES));
    Global.New("Crout",CConstant<string*>(&Crout));
    Global.New("Cholesky",CConstant<string*>(&Cholesky));
    Global.New(UMFPACK.c_str(),CConstant<string*>(&UMFPACK));
    Global.New(sparsesolver.c_str(),CConstant<string*>(&sparsesolver));
    Global.New(sparsesolverSym.c_str(),CConstant<string*>(&sparsesolverSym));

//old --
//  init FESpace
 TheOperators->Add("<-",
		   new OneOperator2_<pfes*,pfes*,pmesh* >(& MakePtr2 ),
		   new OneOperatorCode<OP_MakePtr2>,
		   new OneOperatorCode<OP_MakePtr3>,
           new OneOperatorCode<OP_MakePtrS>,
           new OneOperatorCode<OP_MakePtrL>,
		   new OpMake_pfes<pfes,Mesh,TypeOfFE,pfes_tefk>,
		   new OpMake_pfes<pfes3,Mesh3,TypeOfFE3,pfes3_tefk>,
           new OpMake_pfes<pfesS,MeshS,TypeOfFES,pfesS_tefk>,      // add for 3D surface  FEspace
           new OpMake_pfes<pfesL,MeshL,TypeOfFEL,pfesL_tefk>
        );
    TheOperators->Add("=",new OneOperator2<R3*,R3*,R3* >(&set_eqp));

 Add<MeshPoint*>("P",".", new OneOperator_Ptr_o_R<R3,MeshPoint>(  & MeshPoint::P));
 Add<MeshPoint*>("N",".", new OneOperator_Ptr_o_R<R3,MeshPoint>(  & MeshPoint::N));
 Add<R3*>("x",".", new OneOperator_Ptr_o_R<R,R3>(  & R3::x));
 Add<R3*>("y",".", new OneOperator_Ptr_o_R<R,R3>(  & R3::y));
 Add<R3*>("z",".", new OneOperator_Ptr_o_R<R,R3>(  & R3::z));
 Add<R2*>("x",".", new OneOperator_Ptr_o_R<R,R2>(  & R2::x));
 Add<R2*>("y",".", new OneOperator_Ptr_o_R<R,R2>(  & R2::y));

 Add<R3*>("[","",new OneOperator2<double,R3*,long>(get_R3));

 Add<pmesh>("[","",new OneOperator2_<lgElement,pmesh,long>(get_element));
 Add<pmesh*>("be",".",new OneOperator1_<lgBoundaryEdge::BE,pmesh*>(Build));
 Add<lgElement>("adj",".",new OneOperator1_<lgElement::Adj,lgElement>(Build));
 Add<lgBoundaryEdge::BE>("(","",new OneOperator2_<lgBoundaryEdge,lgBoundaryEdge::BE,long>(get_belement));
 Add<lgElement::Adj>("(","",new OneOperator2_<lgElement,lgElement::Adj,long*>(get_adj));
 TheOperators->Add("==", new OneBinaryOperator<Op2_eq<lgElement,lgElement> >);
 TheOperators->Add("!=", new OneBinaryOperator<Op2_ne<lgElement,lgElement> >);
 TheOperators->Add("<", new OneBinaryOperator<Op2_lt<lgElement,lgElement> >);
 TheOperators->Add("<=", new OneBinaryOperator<Op2_le<lgElement,lgElement> >);


 Add<pmesh*>("[","",new OneOperator2_<lgElement,pmesh*,long>(get_element));
 Add<pmesh>("(","",new OneOperator2_<lgVertex,pmesh,long>(get_vertex));
 Add<pmesh*>("(","",new OneOperator2_<lgVertex,pmesh*,long>(get_vertex));

 Add<lgElement>("[","",new OneOperator2_<lgVertex,lgElement,long>(get_element));
 Add<lgBoundaryEdge>("[","",new OneOperator2_<lgVertex,lgBoundaryEdge,long>(get_belement));

 Add<lgVertex>("x",".",new OneOperator1_<R,lgVertex>(getx));
 Add<lgVertex>("y",".",new OneOperator1_<R,lgVertex>(gety));
 Add<lgVertex>("label",".",new OneOperator1_<long,lgVertex>(getlab));
 Add<lgElement>("label",".",new OneOperator1_<long,lgElement>(getlab));
 Add<lgElement>("region",".",new OneOperator1_<long,lgElement>(getlab));
 Add<lgElement>("area",".",new OneOperator1_<double,lgElement>(getarea));
 Add<lgElement>("mesure",".",new OneOperator1_<double,lgElement>(getarea));
 Add<lgElement>("measure",".",new OneOperator1_<double,lgElement>(getarea));
 Add<lgBoundaryEdge>("length",".",new OneOperator1_<double,lgBoundaryEdge>(getlength));
 Add<lgBoundaryEdge>("label",".",new OneOperator1_<long,lgBoundaryEdge>(getlab));
 Add<lgBoundaryEdge>("Element",".",new OneOperator1_<lgElement,lgBoundaryEdge>(getElement));
 Add<lgBoundaryEdge>("whoinElement",".",new OneOperator1_<long,lgBoundaryEdge>(EdgeElement));


 // New FF language types. zzzfff is defined at [[file:lex.hpp::zzzfff]] as a pointer to an object of class mylex
 // [[file:lex.hpp::class mylex]]. zzzfff->Add() is at [[file:lex.cpp::void mylex Add Key k aType t]]. The lexer will
 // then be called from the parser via [[file:../lglib/lg.ypp::yylex]]

 zzzfff->Add("R3",atype<R3*>());

 // <<mesh_keyword>> pmesh is a pointer to Mesh [[file:../femlib/fem.hpp::class Mesh]] defined at
 // [[file:lgfem.hpp::typedef Mesh pmesh]]
 // pmesh is a pointer to Mesh
 zzzfff->Add("mesh",atype<pmesh*>());
 // pmesh3 is a pointer to Mesh3 defined at [[file:lgfem.hpp::typedef Mesh3 pmesh3]]
 zzzfff->Add("mesh3",atype<pmesh3*>());
 // pmeshS is a pointer to MeshS defined at [[file:lgfem.hpp::typedef MeshS pmeshS]]
 zzzfff->Add("meshS",atype<pmeshS*>());
 // pmeshL is a pointer to MeshL defined at [[file:lgfem.hpp::typedef MeshL pmeshL]]
 zzzfff->Add("meshL",atype<pmeshL*>());

 zzzfff->Add("element",atype<lgElement>());
 zzzfff->Add("vertex",atype<lgVertex>());
 zzzfff->Add("matrix",atype<Matrice_Creuse<R> *>());
 zzzfff->Add("Cmatrix",atype<Matrice_Creuse<Complex> *>()); // a voir



 Global.Add("LinearCG","(",new LinearCG<R>()); // old form  with rhs (must be zer
 Global.Add("LinearGMRES","(",new LinearGMRES<R>()); // old form
 Global.Add("LinearGMRES","(",new LinearGMRES<R>(1)); // old form  without rhs
 Global.Add("AffineGMRES","(",new LinearGMRES<R>(1)); // New  better
 Global.Add("LinearCG","(",new LinearCG<R>(1)); //  without right handsize
 Global.Add("AffineCG","(",new LinearCG<R>(1)); //  without right handsize
 Global.Add("NLCG","(",new LinearCG<R>(-1)); //  without right handsize

 zzzfff->AddF("varf",t_form);    //  var. form ~  <<varf>>
 zzzfff->AddF("solve",t_solve);
 zzzfff->AddF("problem",t_problem);


 Global.Add("jump","(",new OneOperatorCode<Code_VF<Ftest,Code_Jump> >);
 Global.Add("jump","(",new OneOperatorCode<Code_VF<Finconnue,Code_Jump> >);
 Global.Add("average","(",new OneOperatorCode<Code_VF<Ftest,Code_Mean> >);
 Global.Add("average","(",new OneOperatorCode<Code_VF<Finconnue,Code_Mean> >);
 Global.Add("mean","(",new OneOperatorCode<Code_VF<Ftest,Code_Mean> >);
 Global.Add("mean","(",new OneOperatorCode<Code_VF<Finconnue,Code_Mean> >);
 Global.Add("otherside","(",new OneOperatorCode<Code_VF<Ftest,Code_OtherSide> >);
 Global.Add("otherside","(",new OneOperatorCode<Code_VF<Finconnue,Code_OtherSide> >);

 Global.Add("conj","(",new OneOperatorCode<CODE_conj<Finconnue> >);
 Global.Add("conj","(",new OneOperatorCode<CODE_conj<Ftest> >);
 Global.Add("conj","(",new OneOperatorCode<CODE_conj<Foperator> >);
 TheOperators->Add("\'", new OneOperatorCode<CODE_conj<Finconnue> >);
 TheOperators->Add("\'", new OneOperatorCode<CODE_conj<Ftest> >);
 TheOperators->Add("\'", new OneOperatorCode<CODE_conj<Foperator> >);

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

 Global.Add("on","(",new OneOperatorCode<BC_set > );

 /// <<plot_keyword>> uses [[Plot]] and [[file:AFunction.hpp::OneOperatorCode]] and [[file:AFunction.hpp::Global]]
 Global.Add("plot","(",new OneOperatorCode<Plot> );
 Global.Add("convect","(",new OneOperatorCode<Convect> );


 TheOperators->Add("+",
    new OneOperatorCode<CODE_L_Add<Foperator> > ,
    new OneOperatorCode<CODE_L_Add<Ftest> > ,
    new OneOperatorCode<CODE_L_Add<Finconnue> > ,
    new OneOperatorCode<C_args>(t_C_args,t_C_args,t_C_args) // ,
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
      new OneOperatorCode<C_args>(t_C_args,t_BC)
    );

  atype<const C_args*>()->AddCast(
      new OneOperatorCode<C_args>(t_C_args,atype<DotSlash_KN_<R> >()) ,
      new OneOperatorCode<C_args>(t_C_args,atype<KN<R>*>())  ,
      new OneOperatorCode<C_args>(t_C_args,atype<DotStar_KN_<R> >())  ,
      new OneOperatorCode<C_args>(t_C_args,atype<Matrice_Creuse<R>*>()) ,
      new OneOperatorCode<C_args>(t_C_args,atype<RNM_VirtualMatrix<R>::plusAx >()),
      new OneOperatorCode<C_args>(t_C_args,atype<RNM_VirtualMatrix<R>::plusAtx >())

    );

  atype<const C_args*>()->AddCast(
      new OneOperatorCode<C_args>(t_C_args,atype<DotSlash_KN_<Complex> >()) ,
      new OneOperatorCode<C_args>(t_C_args,atype<KN<Complex>*>())  ,
      new OneOperatorCode<C_args>(t_C_args,atype<DotStar_KN_<Complex> >())  ,
      new OneOperatorCode<C_args>(t_C_args,atype<Matrice_Creuse<Complex>*>()) ,
      new OneOperatorCode<C_args>(t_C_args,atype<RNM_VirtualMatrix<Complex>::plusAx >()),
      new OneOperatorCode<C_args>(t_C_args,atype<RNM_VirtualMatrix<Complex>::plusAtx >())

    );

 TheOperators->Add("*",
    new OneOperatorCode<CODE_L_Mul<Foperator,Ftest,Finconnue> > ,
    new OneOperatorCode<CODE_L_Mul<Foperator,Finconnue,Ftest> > );

// Warning just double or complex in following operator
// ----------------------------------------------------
//   in case of  ambiguity we take the double version
//   case    long -> double
//           long -> complex
 TheOperators->Add("*",
    new OneOperatorCode<CODE_L_MulLR<Finconnue,double>, 20> ,
    new OneOperatorCode<CODE_L_MulLR<Foperator,double>, 20 > ,
    new OneOperatorCode<CODE_L_MulLR<Ftest,double>, 20 > ,
    new OneOperatorCode<CODE_L_MulRL<double,Finconnue>, 20 > ,
    new OneOperatorCode<CODE_L_MulRL<double,Foperator>, 20 > ,
    new OneOperatorCode<CODE_L_MulRL<double,Ftest>, 20 >
  );
 TheOperators->Add("*",
    new OneOperatorCode<CODE_L_MulLR<Finconnue,Complex>, 10 > ,
    new OneOperatorCode<CODE_L_MulLR<Foperator,Complex>, 10 > ,
    new OneOperatorCode<CODE_L_MulLR<Ftest,Complex>, 10 > ,
    new OneOperatorCode<CODE_L_MulRL<Complex,Finconnue>, 10 > ,
    new OneOperatorCode<CODE_L_MulRL<Complex,Foperator>, 10 > ,
    new OneOperatorCode<CODE_L_MulRL<Complex,Ftest>, 10 >
  );

 TheOperators->Add("/",
    new OneOperatorCode<CODE_L_DivLR<Finconnue,double>, 20 > ,
    new OneOperatorCode<CODE_L_DivLR<Foperator,double>, 20 > ,
    new OneOperatorCode<CODE_L_DivLR<Ftest,double>, 20 > );

 TheOperators->Add("/",
    new OneOperatorCode<CODE_L_DivLR<Finconnue,Complex>, 10 > ,
    new OneOperatorCode<CODE_L_DivLR<Foperator,Complex>, 10 > ,
    new OneOperatorCode<CODE_L_DivLR<Ftest,Complex>, 10 > );

// Warning just double or complex in previous operator
// ----------------------------------------------------

// TheOperators->Add("=",new OneOperatorCode<BC_set1<double> >);


 TheOperators->Add("=",
		   new OneOperator2<pmesh*,pmesh*,pmesh >(&set_eqdestroy_incr),

		   new OneBinaryOperator<set_eq_array<KN_<double> ,RNM_VirtualMatrix<R>::plusAx > > ,
		   new OneBinaryOperator<set_eq_array<KN_<double> ,RNM_VirtualMatrix<R>::plusAtx > >  , //ZZ set_eq_array
		   new OneBinaryOperator<set_eq_array<KN_<double> ,RNM_VirtualMatrix<R>::solveAxeqb > >  ,
                   new OneBinaryOperator<set_eq_array<KN_<double> ,RNM_VirtualMatrix<R>::solveAtxeqb > >  ,

		   new OpArraytoLinearForm<double,v_fes>(atype< KN_<double> >(),false,false)  ,
		   new OpMatrixtoBilinearForm<double,v_fes >);


 TheOperators->Add("=",
		   new OpArraytoLinearForm<double,v_fes3>(atype< KN_<double> >(),false,false)  ,//3D volume
		   new OpMatrixtoBilinearForm<double,v_fes3 > , // 3D volume
           new OpArraytoLinearForm<double,v_fesS>(atype< KN_<double> >(),false,false)  , // 3D surface
           new OpMatrixtoBilinearForm<double,v_fesS > , // 3D surface
           new OpArraytoLinearForm<double,v_fesL>(atype< KN_<double> >(),false,false)  , // 3D curve
           new OpMatrixtoBilinearForm<double,v_fesL >); // 3D curve


 TheOperators->Add("<-",
		   new OpArraytoLinearForm<double,v_fes>(atype< KN<double>* >(),true,true) ,
		   new OpArraytoLinearForm<Complex,v_fes>(atype< KN<Complex>* >(),true,true) ,
		   new OpArraytoLinearForm<double,v_fes3>(atype< KN<double>* >(),true,true) , //3D volume
		   new OpArraytoLinearForm<Complex,v_fes3>(atype< KN<Complex>* >(),true,true), //3D volume
           new OpArraytoLinearForm<double,v_fesS>(atype< KN<double>* >(),true,true) , //3D surface
           new OpArraytoLinearForm<Complex,v_fesS>(atype< KN<Complex>* >(),true,true),  //3D surface
           new OpArraytoLinearForm<double,v_fesL>(atype< KN<double>* >(),true,true) , //3D curve
           new OpArraytoLinearForm<Complex,v_fesL>(atype< KN<Complex>* >(),true,true) //3D curve
        );


 TheOperators->Add("=",
		    new OneBinaryOperator<set_eq_array<KN_<Complex> ,RNM_VirtualMatrix<Complex>::plusAx > > ,
		    new OneBinaryOperator<set_eq_array<KN_<Complex> ,RNM_VirtualMatrix<Complex>::plusAtx > >  ,
		    new OneBinaryOperator<set_eq_array<KN_<Complex> ,RNM_VirtualMatrix<Complex>::solveAxeqb > >  ,
                   new OneBinaryOperator<set_eq_array<KN_<Complex> ,RNM_VirtualMatrix<Complex>::solveAtxeqb > >  ,

		   new OpArraytoLinearForm<Complex,v_fes>(atype< KN<Complex>* >(),true,false)  ,
		   new OpMatrixtoBilinearForm<Complex,v_fes >);

 TheOperators->Add("=",
		   new OpArraytoLinearForm<Complex,v_fes3>(atype< KN_<Complex> >(),false,false)   , //3D volume
		   new OpMatrixtoBilinearForm<Complex,v_fes3 > , //3D volume
           new OpArraytoLinearForm<Complex,v_fesS>(atype< KN_<Complex> >(),false,false)   , //3D surface
           new OpMatrixtoBilinearForm<Complex,v_fesS > , //3D surface
           new OpArraytoLinearForm<Complex,v_fesL>(atype< KN_<Complex> >(),false,false)   , //3D curve
           new OpMatrixtoBilinearForm<Complex,v_fesL >) ; //3D surface

 // add august 2007
 TheOperators->Add("<-",
		   new OneBinaryOperator<init_eqarray<KN<double> ,RNM_VirtualMatrix<double>::plusAx > > ,
		   new OneBinaryOperator<init_eqarray<KN<double> ,RNM_VirtualMatrix<double>::plusAtx > >  ,
		   new OneBinaryOperator<init_eqarray<KN<double> ,RNM_VirtualMatrix<double>::solveAxeqb > >  ,
           new OneBinaryOperator<init_eqarray<KN<double> ,RNM_VirtualMatrix<double>::solveAtxeqb > >  ,

		   new OneBinaryOperator<init_eqarray<KN<Complex> ,RNM_VirtualMatrix<Complex>::plusAx > > ,
		   new OneBinaryOperator<init_eqarray<KN<Complex> ,RNM_VirtualMatrix<Complex>::plusAtx > >  ,
		   new OneBinaryOperator<init_eqarray<KN<Complex> ,RNM_VirtualMatrix<Complex>::solveAxeqb > >,
                   new OneBinaryOperator<init_eqarray<KN<Complex> ,RNM_VirtualMatrix<Complex>::solveAtxeqb > >

		);
    // Jan 2018  FH  Vh uh=yu[]; //  A FAIRE FH POUR F NATAF FFFFFFFF
    TheOperators->Add("<-"
                      , new init_FE_eqarray<FFset3<R,v_fes,RNM_VirtualMatrix<R>::plusAx> >(10)
                      , new init_FE_eqarray<FFset3<R,v_fes,RNM_VirtualMatrix<R>::solveAxeqb> >(10)
                      , new init_FE_eqarray<FFset3<R,v_fes,RNM_VirtualMatrix<R>::solveAtxeqb> >(10)
                      , new init_FE_eqarray<FFset3<R,v_fes,RNM_VirtualMatrix<R>::plusAtx> >(10)
                      , new init_FE_eqarray<FFset3<R,v_fes,KN_<R> > >(10)
                      , new init_FE_eqarray<FF_L_args<R,v_fes,Call_FormLinear<v_fes> > >(10)

                      );
    TheOperators->Add("<-"
                      , new init_FE_eqarray<FFset3<Complex,v_fes,RNM_VirtualMatrix<Complex>::plusAx> >(10)
                      , new init_FE_eqarray<FFset3<Complex,v_fes,RNM_VirtualMatrix<Complex>::solveAxeqb> >(10)
                      , new init_FE_eqarray<FFset3<Complex,v_fes,RNM_VirtualMatrix<Complex>::solveAtxeqb> >(10)
                      , new init_FE_eqarray<FFset3<Complex,v_fes,RNM_VirtualMatrix<Complex>::plusAtx> >(10)
                      , new init_FE_eqarray<FFset3<Complex,v_fes,KN_<Complex> > >(10)
                      , new init_FE_eqarray<FF_L_args<Complex,v_fes,Call_FormLinear<v_fes> > >(10)
                      );
    TheOperators->Add("<-"
                      , new init_FE_eqarray<FFset3<R,v_fes3,RNM_VirtualMatrix<R>::plusAx> >(10)  // 3D volume
                      , new init_FE_eqarray<FFset3<R,v_fes3,RNM_VirtualMatrix<R>::solveAxeqb> >(10)  // 3D volume
                      , new init_FE_eqarray<FFset3<R,v_fes3,RNM_VirtualMatrix<R>::plusAtx> >(10)  // 3D volume
                      , new init_FE_eqarray<FFset3<R,v_fes3,KN_<R> > >(10)  // 3D volume
                      , new init_FE_eqarray<FF_L_args<R,v_fes3,Call_FormLinear<v_fes3> > >(10)  // 3D volume

                      );
    TheOperators->Add("<-"
                      , new init_FE_eqarray<FFset3<Complex,v_fes3,RNM_VirtualMatrix<Complex>::plusAx> >(10)  // 3D volume
                      , new init_FE_eqarray<FFset3<Complex,v_fes3,RNM_VirtualMatrix<Complex>::solveAxeqb> >(10)  // 3D volume
                      , new init_FE_eqarray<FFset3<Complex,v_fes3,RNM_VirtualMatrix<Complex>::plusAtx> >(10)  // 3D volume
                      , new init_FE_eqarray<FFset3<Complex,v_fes3,KN_<Complex> > >(10)  // 3D volume
                      , new init_FE_eqarray<FF_L_args<Complex,v_fes3,Call_FormLinear<v_fes3> > >(10)  // 3D volume
                      );

    // Fin
    // Fin
 TheOperators->Add("+=",

		   new OneBinaryOperator<set_eq_array_add<KN_<double> ,RNM_VirtualMatrix<double>::plusAx > > ,
		   new OneBinaryOperator<set_eq_array_add<KN_<double> ,RNM_VirtualMatrix<double>::plusAtx > >,

		   new OpArraytoLinearForm<double,v_fes>(atype< KN_<double> >(),false,false,false)  ,
		   new OpArraytoLinearForm<double,v_fes3>(atype< KN_<double> >(),false,false,false) ,  // 3D volume
		   new OpArraytoLinearForm<double,v_fesS>(atype< KN_<double> >(),false,false,false) ,  // 3D surface
           new OpArraytoLinearForm<double,v_fesL>(atype< KN_<double> >(),false,false,false)    // 3D surface
       );

 TheOperators->Add("+=",

		   new OneBinaryOperator<set_eq_array_add<KN_<Complex> ,RNM_VirtualMatrix<Complex>::plusAx > > ,
		   new OneBinaryOperator<set_eq_array_add<KN_<Complex> ,RNM_VirtualMatrix<Complex>::plusAtx > >,

		   new OpArraytoLinearForm<Complex,v_fes>(atype< KN_<Complex> >(),false,false,false)  ,

		   new OpArraytoLinearForm<Complex,v_fes3>(atype< KN_<Complex> >(),false,false,false) , // 3D volume
           new OpArraytoLinearForm<Complex,v_fesS>(atype< KN_<Complex> >(),false,false,false) , // 3D surface
           new OpArraytoLinearForm<Complex,v_fesL>(atype< KN_<Complex> >(),false,false,false)   // 3D curve
       );



 TheOperators->Add("<-",new OpMatrixtoBilinearForm<double,v_fes >(1) );
 TheOperators->Add("<-",new OpMatrixtoBilinearForm<Complex,v_fes >(1) );

 TheOperators->Add("<-",new OpMatrixtoBilinearForm<double,v_fes3 >(1) ); // 3D volume
 TheOperators->Add("<-",new OpMatrixtoBilinearForm<Complex,v_fes3 >(1) );// 3D volume

 TheOperators->Add("<-",new OpMatrixtoBilinearForm<double,v_fesS >(1) ); // 3D surface
 TheOperators->Add("<-",new OpMatrixtoBilinearForm<Complex,v_fesS >(1) ); // 3D surface
    
 TheOperators->Add("<-",new OpMatrixtoBilinearForm<double,v_fesL >(1) ); // 3D curve
 TheOperators->Add("<-",new OpMatrixtoBilinearForm<Complex,v_fesL >(1) ); // 3D curve

 Add<const  FormLinear   *>("(","",new OpCall_FormLinear<FormLinear,v_fes> );
 Add<const  FormBilinear *>("(","",new OpCall_FormBilinear<FormBilinear,v_fes> );
 Add<const  FormBilinear *>("(","",new OpCall_FormLinear2<FormBilinear,v_fes> );
 Add<const C_args*>("(","",new OpCall_FormLinear2<C_args,v_fes>);
 Add<const C_args*>("(","",new OpCall_FormBilinear<C_args,v_fes> );

 Add<const  FormLinear   *>("(","",new OpCall_FormLinear<FormLinear,v_fes3> );  // 3D volume
 Add<const  FormBilinear *>("(","",new OpCall_FormBilinear<FormBilinear,v_fes3> );  // 3D volume
 Add<const  FormBilinear *>("(","",new OpCall_FormLinear2<FormBilinear,v_fes3> );  // 3D volume
 Add<const C_args*>("(","",new OpCall_FormLinear2<C_args,v_fes3>);  // 3D volume
 Add<const C_args*>("(","",new OpCall_FormBilinear<C_args,v_fes3> );  // 3D volume

 Add<const  FormLinear   *>("(","",new OpCall_FormLinear<FormLinear,v_fesS> );  // 3D surface
 Add<const  FormBilinear *>("(","",new OpCall_FormBilinear<FormBilinear,v_fesS> );  // 3D surface
 Add<const  FormBilinear *>("(","",new OpCall_FormLinear2<FormBilinear,v_fesS> );  // 3D surface
 Add<const C_args*>("(","",new OpCall_FormLinear2<C_args,v_fesS>);  // 3D surface
 Add<const C_args*>("(","",new OpCall_FormBilinear<C_args,v_fesS> );  // 3D surface
    
 Add<const  FormLinear   *>("(","",new OpCall_FormLinear<FormLinear,v_fesL> );  // 3D curve
 Add<const  FormBilinear *>("(","",new OpCall_FormBilinear<FormBilinear,v_fesL> );  // 3D curve
 Add<const  FormBilinear *>("(","",new OpCall_FormLinear2<FormBilinear,v_fesL> );  // 3D curve
 Add<const C_args*>("(","",new OpCall_FormLinear2<C_args,v_fesL>);  // 3D curve
 Add<const C_args*>("(","",new OpCall_FormBilinear<C_args,v_fesL> );  // 3D curve

//  correction du bug morale
//  Attention il y a moralement un bug
//  les initialisation   x = y   ( passe par l'operateur binaire <-  dans TheOperators
//   les initialisation   x(y)   ( passe par l'operateur unaire <-  de typedebase de x (inutile 2007).
//  x(y1,..,yn) est un operator n+1   (x,y1,..,yn)
// on passe toujours par x(y) maintenant.
//   -------


 TheOperators->Add("<-",
       new OneOperator2_<pferbase*,pferbase*,pfes* >(MakePtrFE2),
       new OneOperator3_<pferbasearray*,pferbasearray*,pfes*,long >(MakePtrFE3),

       new OneOperator2_<pfecbase*,pfecbase*,pfes* >(MakePtrFE2),
       new OneOperator3_<pfecbasearray*,pfecbasearray*,pfes*,long >(MakePtrFE3)


       );
 TheOperators->Add("<-",
       new OneOperatorMakePtrFE<double>(atype<double>()),  //  scalar case
       new OneOperatorMakePtrFE<double>(atype<E_Array>()),  //  vect case
       new OneOperatorMakePtrFE<Complex>(atype<Complex>()),  //  scalar complex  case
       new OneOperatorMakePtrFE<Complex>(atype<E_Array>())  //  vect complex case
       );
//  interpolation   operator
 TheOperators->Add("=",
       new OneOperator2_<pfer,pfer,double,E_F_StackF0F0opt2<double> >(set_fe<double>) ,
       new OneOperator2_<pfec,pfec,Complex,E_F_StackF0F0opt2<Complex> >(set_fe<Complex>)

       ) ;

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
       new OneBinaryOperator<PrintPnd<Matrice_Creuse<R>*> >,
       new OneBinaryOperator<PrintPnd<Matrice_Creuse<Complex>*> >

       );

 TheOperators->Add(">>",
		   new OneBinaryOperator<Op_Read<Matrice_Creuse<R> > >,
		   new OneBinaryOperator<Op_Read<Matrice_Creuse<Complex> > >

		   );


 Global.Add("int2d","(",new OneOperatorCode<CDomainOfIntegration>);
 Global.Add("int1d","(",new OneOperatorCode<CDomainOfIntegrationBorder>);
 Global.Add("intalledges","(",new OneOperatorCode<CDomainOfIntegrationAllEdges>);
 Global.Add("intallVFedges","(",new OneOperatorCode<CDomainOfIntegrationVFEdges>);
 Global.Add("jump","(",new OneUnaryOperator<JumpOp<R>,JumpOp<R> >);
 Global.Add("mean","(",new OneUnaryOperator<MeanOp<R>,MeanOp<R> >);
 Global.Add("average","(",new OneUnaryOperator<MeanOp<R>,MeanOp<R> >);
 Global.Add("otherside","(",new OneUnaryOperator<OthersideOp<R>,OthersideOp<R> >);

 Global.Add("jump","(",new OneUnaryOperator<JumpOp<Complex>,JumpOp<Complex> >);
 Global.Add("mean","(",new OneUnaryOperator<MeanOp<Complex>,MeanOp<Complex> >);
 Global.Add("average","(",new OneUnaryOperator<MeanOp<Complex>,MeanOp<Complex> >);
 Global.Add("otherside","(",new OneUnaryOperator<OthersideOp<Complex>,OthersideOp<Complex> >);


 Add<const CDomainOfIntegration*>("(","",new OneOperatorCode<FormBilinear> );
 Add<const CDomainOfIntegration *>("(","",new OneOperatorCode<FormLinear> );

 Add<const CDomainOfIntegration *>("(","",new OneOperatorCode<IntFunction<double>,1 >);
 Add<const CDomainOfIntegration *>("(","",new OneOperatorCode<IntFunction<complex<double> >,0 >);



 map_type[typeid(double).name()]->AddCast(
   new E_F1_funcT<double,pfer>(pfer2R<R,0>)
   );

 map_type[typeid(Complex).name()]->AddCast(
   new E_F1_funcT<Complex,pfec>(pfer2R<Complex,0>)
   );

// bof
 Global.Add("dx","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dx>));
 Global.Add("dy","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dy>));
 Global.Add("dxx","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dxx>));
 Global.Add("dyy","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dyy>));
 Global.Add("dxy","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dxy>));
 Global.Add("dyx","(",new E_F1_funcT<double,pfer>(pfer2R<R,op_dyx>));

 Global.Add("dx","(",new E_F1_funcT<Complex,pfec>(pfer2R<Complex,op_dx>));
 Global.Add("dy","(",new E_F1_funcT<Complex,pfec>(pfer2R<Complex,op_dy>));
 Global.Add("dxx","(",new E_F1_funcT<Complex,pfec>(pfer2R<Complex,op_dxx>));
 Global.Add("dyy","(",new E_F1_funcT<Complex,pfec>(pfer2R<Complex,op_dyy>));
 Global.Add("dxy","(",new E_F1_funcT<Complex,pfec>(pfer2R<Complex,op_dxy>));
 Global.Add("dyx","(",new E_F1_funcT<Complex,pfec>(pfer2R<Complex,op_dyx>));

  Add<pfecbasearray*>("[","",new OneOperator2_<pfecbase*,pfecbasearray*,long>(get_element));  // use FH sep. 2009
  Add<pferbasearray*>("[","",new OneOperator2_<pferbase*,pferbasearray*,long>(get_element));  //  use ???? FH sep. 2009
    // bof bof ..
    // resize of array of Finite element ..  a little hard 2013 FH
    Dcl_Type< Resize1<pfecbasearray* > > ();
    Dcl_Type< Resize1<pferbasearray* > > ();
    Dcl_Type< Resize1<pfecarray > > ();
    Dcl_Type< Resize1<pferarray > > ();
    Add<pfecbasearray*>("resize",".",new OneOperator1<Resize1<pfecbasearray* >,pfecbasearray*>(to_Resize1));  //  FH fev 2013
    Add<pferbasearray*>("resize",".",new OneOperator1<Resize1<pferbasearray* >,pferbasearray*>(to_Resize1));  //   FH fev. 2013
    Add<pferarray>("resize",".",new OneOperator1<Resize1<pferarray >,pferarray>(to_Resize1));  //  FH fev 2013
    Add<pfecarray>("resize",".",new OneOperator1<Resize1<pfecarray >,pfecarray>(to_Resize1));  //   FH fev. 2013
    new OneOperator2_<pferbasearray*,Resize1<pferbasearray* > , long  >(fepresize<pferbasearray*>);
    Add<Resize1<pferbasearray* > >("(","",new  OneOperator2_<pferbasearray*,Resize1<pferbasearray* > , long  >(fepresize));
    Add<Resize1<pfecbasearray* > >("(","",new OneOperator2_<pfecbasearray*,Resize1<pfecbasearray* > , long  >(fepresize));
    Add<Resize1<pferarray > >("(","",new OneOperator2_<pferarray,Resize1<pferarray > , long  >(feresize));
    Add<Resize1<pfecarray > >("(","",new OneOperator2_<pfecarray,Resize1<pfecarray > , long  >(feresize));

    Dcl_Type< Resize1<pf3rbasearray* > > ();
    Dcl_Type< Resize1<pf3rarray > > ();
    Dcl_Type< Resize1<pf3cbasearray* > > (); //   FH oct. 2016
    Dcl_Type< Resize1<pf3carray > > (); //   FH oct. 2016
    Add<pf3rbasearray*>("resize",".",new OneOperator1<Resize1<pf3rbasearray* >,pf3rbasearray*>(to_Resize1));  //   FH fev. 2013
    Add<pf3rarray>("resize",".",new OneOperator1<Resize1<pf3rarray >,pf3rarray>(to_Resize1));  //  FH fev 2013
    Add<pf3cbasearray*>("resize",".",new OneOperator1<Resize1<pf3cbasearray* >,pf3cbasearray*>(to_Resize1));  //   FH oct. 2016
    Add<pf3carray>("resize",".",new OneOperator1<Resize1<pf3carray >,pf3carray>(to_Resize1));  //  FH Oct 2016

    Add<Resize1<pf3rbasearray* > >("(","",new  OneOperator2_<pf3rbasearray*,Resize1<pf3rbasearray* > , long  >(fepresize));
    Add<Resize1<pf3rarray > >("(","",new OneOperator2_<pf3rarray,Resize1<pf3rarray > , long  >(feresize));
    Add<Resize1<pf3cbasearray* > >("(","",new  OneOperator2_<pf3cbasearray*,Resize1<pf3cbasearray* > , long  >(fepresize));//  FH Oct 2016
    Add<Resize1<pf3carray > >("(","",new OneOperator2_<pf3carray,Resize1<pf3carray > , long  >(feresize));//  FH Oct 2016

    Dcl_Type< Resize1<pfSrbasearray* > > ();
    Dcl_Type< Resize1<pfSrarray > > ();
    Dcl_Type< Resize1<pfScbasearray* > > ();
    Dcl_Type< Resize1<pfScarray > > ();
    Add<pfSrbasearray*>("resize",".",new OneOperator1<Resize1<pfSrbasearray* >,pfSrbasearray*>(to_Resize1));
    Add<pfSrarray>("resize",".",new OneOperator1<Resize1<pfSrarray >,pfSrarray>(to_Resize1));
    Add<pfScbasearray*>("resize",".",new OneOperator1<Resize1<pfScbasearray* >,pfScbasearray*>(to_Resize1));
    Add<pfScarray>("resize",".",new OneOperator1<Resize1<pfScarray >,pfScarray>(to_Resize1));

    Add<Resize1<pfSrbasearray* > >("(","",new  OneOperator2_<pfSrbasearray*,Resize1<pfSrbasearray* > , long  >(fepresize));
    Add<Resize1<pfSrarray > >("(","",new OneOperator2_<pfSrarray,Resize1<pfSrarray > , long  >(feresize));
    Add<Resize1<pfScbasearray* > >("(","",new  OneOperator2_<pfScbasearray*,Resize1<pfScbasearray* > , long  >(fepresize));
    Add<Resize1<pfScarray > >("(","",new OneOperator2_<pfScarray,Resize1<pfScarray > , long  >(feresize));

// end of resize ...

  Add<pfecbasearray*>("n",".",new OneOperator1_<long,pfecbasearray*>(get_size));  //  FH fev 2013
  Add<pferbasearray*>("n",".",new OneOperator1_<long,pferbasearray*>(get_size));  //   FH fev. 2013
  Add<pferarray>("n",".",new OneOperator1_<long,pferarray>(get_size));  //  FH fev 2013
  Add<pfecarray>("n",".",new OneOperator1_<long,pfecarray>(get_size));  //   FH fev. 2013

  Add<pf3rbasearray*>("n",".",new OneOperator1_<long,pf3rbasearray*>(get_size));  //   FH fev. 2013
  Add<pf3rarray>("n",".",new OneOperator1_<long,pf3rarray>(get_size));  //  FH fev 2013
  Add<pf3cbasearray*>("n",".",new OneOperator1_<long,pf3cbasearray*>(get_size));  //   FH  Oct 2016
  Add<pf3carray>("n",".",new OneOperator1_<long,pf3carray>(get_size));  //  FH Oct 2016


  Add<pferarray>("[","",new OneOperator2_FE_get_elmnt<double,v_fes>());// new version FH sep 2009
  Add<pfecarray>("[","",new OneOperator2_FE_get_elmnt<Complex,v_fes>());

  Add<pf3rarray>("[","",new OneOperator2_FE_get_elmnt<double,v_fes>());// new version FH sep 2009
  Add<pf3carray>("[","",new OneOperator2_FE_get_elmnt<Complex,v_fes>());// FH Oct 2016

  TheOperators->Add("\'",
       new OneOperator1<Matrice_Creuse_Transpose<R>,Matrice_Creuse<R> *>(&Build<Matrice_Creuse_Transpose<R>,Matrice_Creuse<R> *>),
       new OneOperator1<Matrice_Creuse_Transpose<Complex>,Matrice_Creuse<Complex> *>(&Build<Matrice_Creuse_Transpose<Complex>,Matrice_Creuse<Complex> *>)
  );

  Add<pfer>("(","",new interpolate_f_X_1<R> );
  TheOperators->Add("=", new OneOperator2_<void,interpolate_f_X_1<R>::type,double,E_F_StackF0F0 >(set_feoX_1) ) ;
  init_lgmat();
  init_mesh_array();

  l2interpreter = new LinkToInterpreter;
  using namespace FreeFempp;
  FreeFempp::TypeVarForm<double>::Global = new TypeVarForm<double>();
  FreeFempp::TypeVarForm<Complex>::Global = new TypeVarForm<Complex>();


  Global.New("P13d",CConstantTFE3(&DataFE<Mesh3>::P1));
  Global.New("P23d",CConstantTFE3(&DataFE<Mesh3>::P2));
  Global.New("P03d",CConstantTFE3(&DataFE<Mesh3>::P0));
  Global.New("RT03d",CConstantTFE3(&RT03d));
  Global.New("Edge03d",CConstantTFE3(&Edge03d));
  Global.New("P1b3d",CConstantTFE3(&P1bLagrange3d));

  Global.New("RT0S",CConstantTFES(&DataFE<MeshS>::RT0));
  Global.New("P2S",CConstantTFES(&DataFE<MeshS>::P2));
  Global.New("P1S",CConstantTFES(&DataFE<MeshS>::P1));
  Global.New("P0S",CConstantTFES(&DataFE<MeshS>::P0));
  Global.New("P1bS",CConstantTFES(&P1bLagrange_surf));
    
  Global.New("P2L",CConstantTFEL(&DataFE<MeshL>::P2));
  Global.New("P1L",CConstantTFEL(&DataFE<MeshL>::P1));
  Global.New("P0L",CConstantTFEL(&DataFE<MeshL>::P0));

  TEF2dtoS[FindFE2("P0")]=&DataFE<MeshS>::P0;
  TEF2dtoS[FindFE2("P1")]=&DataFE<MeshS>::P1;
  TEF2dtoS[FindFE2("P2")]=&DataFE<MeshS>::P2;
  TEF2dtoS[FindFE2("P1b")]=&P1bLagrange_surf;
  TEF2dtoS[FindFE2("RT0")]=&DataFE<MeshS>::RT0;
    
  TEF2dtoL[FindFE2("P0")]=&DataFE<MeshL>::P0;
  TEF2dtoL[FindFE2("P1")]=&DataFE<MeshL>::P1;
  TEF2dtoL[FindFE2("P2")]=&DataFE<MeshL>::P2;

  TEF2dto3d[FindFE2("P1")]=&DataFE<Mesh3>::P1;
  TEF2dto3d[FindFE2("P2")]=&DataFE<Mesh3>::P2;
  TEF2dto3d[FindFE2("P0")]=&DataFE<Mesh3>::P0;
  TEF2dto3d[FindFE2("P1b")]=&P1bLagrange3d;
  TEF2dto3d[FindFE2("RT0")]=&RT03d;
}

void clean_lgfem()
{

  delete l2interpreter;
  delete  FreeFempp::TypeVarForm<double>::Global;
  delete  FreeFempp::TypeVarForm<Complex>::Global;
}
template<class K,class v_fes>
Expression IsFEcomp(const C_F0 &c,int i, Expression & rrr,Expression &iii)
{
    Expression r=0;
    if(!i) rrr=0,iii=0;
    if(atype<typename E_FEcomp<K,v_fes>::Result>() == c.left())
      {
	  const E_FEcomp<K,v_fes> * e= dynamic_cast<const E_FEcomp<K,v_fes>*>(c.LeftValue() );
	  if( !e)
	    {
		const E_FEcomp_get_elmnt_array<K,v_fes> * ee= dynamic_cast<const E_FEcomp_get_elmnt_array<K,v_fes>*>(c.LeftValue() );

		if( !ee)
		  {
		      cerr <<" Fatal error Copy array .." << *c.left()<< " composante : " << i << endl;
		      ffassert(ee);
		  }
		else
		  {
		      if (ee->comp ==i) {
			  if (i && ee->a00->a0 != rrr) cerr << " error composante arry vect. incompatible " << ee->comp << " "<< ee->a00->a0 << " != " << rrr << endl;
			  else {
			   r= ee->a0;
			    rrr= ee->a00->a0;
			    iii= ee->a1;
			  }

		      }
		      else cerr << " erreur composante " << ee->comp << " != " << i << endl;


		  }
	    }
	  else
	    {

		if (e->comp ==i) {
		      if (i && e->a0 != rrr) cerr << " error composante incompatible " << e->comp  << endl;
		      else  rrr=r=e->a0;
		}
	      else cerr << " erreur composante " << e->comp << " != " << endl;
	  }
      }
    return r;
}

template<class K,class v_fes>
Expression Op_CopyArrayT(const E_Array & a,const E_Array & b)
{
    typedef FEbaseArray<K,v_fes>  FEba;

    int nb=b.size();
    Expression r=0,rr=0,rrr,iii;
    //  try real voctor value FE interpolation
    rr=IsFEcomp<K,v_fes>(a[0],0,rrr,iii) ;
    if (rr !=0)
      {
          for (int i=1;i<nb;i++)
	      if (!IsFEcomp<K,v_fes>(a[i],i,rrr,iii))
		  CompileError("Copy of Array with incompatible K  vector value FE function () !");;
	  if(iii) {
	      C_F0 aa(rrr,atype<FEba**>()),ii(iii,atype<long>());
	      C_F0 aa_ii(aa,"[",ii);
	      rr=aa_ii.LeftValue();
	  }
	  if(v_fes::dHat==2 && v_fes::d==2) r=new E_set_fev<K>(&b,rr,2);
	  else if(v_fes::dHat==3 && v_fes::d==3) r=new  E_set_fev3<K,v_fes3>(&b,rr);
      else if(v_fes::dHat==2 && v_fes::d==3) r=new  E_set_fevS<K,v_fesS>(&b,rr);
      else if(v_fes::dHat==1 && v_fes::d==3) r=new  E_set_fevL<K,v_fesL>(&b,rr);
      }
    //  try complex vector value FE interpolation
    return r;
}

 E_F0 * Op_CopyArray::code(const basicAC_F0 & args) const  {
       E_F0 * ret=0;
      const E_Array & a= *dynamic_cast<const E_Array*>(args[0].LeftValue());
      const E_Array & b= *dynamic_cast<const E_Array*>(args[1].LeftValue());
      int na=a.size();
      int nb=b.size();
      if (na != nb )
        CompileError("Copy of Array with incompatible size!");
      if(0)
	{ // old code !!!!!!! before removing FH sept. 2009
       Expression rr=0,rrr,iii;
       //  try real voctor value FE interpolation
       rr=IsFEcomp<double,v_fes>(a[0],0,rrr,iii) ;
       if (rr !=0)
        {
          for (int i=1;i<nb;i++)
	    if (!IsFEcomp<double,v_fes>(a[i],i,rrr,iii))
                CompileError("Copy of Array with incompatible real vector value FE function () !");;
	  return  new E_set_fev<double>(&b,rr,2);
         }
       //  try complex vector value FE interpolation

       rr=IsFEcomp<Complex,v_fes>(a[0],0,rrr,iii) ;
       if (rr !=0)
        {
          for (int i=1;i<nb;i++)
	    if (!IsFEcomp<Complex,v_fes>(a[i],i,rrr,iii))
	      CompileError("Copy of Array with incompatible complex vector value FE function () !");;
	  return  new E_set_fev<Complex>(&b,rr,2);
	}

       rr=IsFEcomp<double,v_fes3>(a[0],0,rrr,iii) ;
       if (rr !=0)
        {
          for (int i=1;i<nb;i++)
	    if (!IsFEcomp<double,v_fes3>(a[i],i,rrr,iii))
                CompileError("Copy of Array with incompatible real vector value FE function () !");;
	  return  new E_set_fev3<double,v_fes3>(&b,rr);
         }
       //  try complex vector value FE interpolation

       rr=IsFEcomp<Complex,v_fes3>(a[0],0,rrr,iii) ;
       if (rr !=0)
        {
          for (int i=1;i<nb;i++)
	    if (!IsFEcomp<Complex,v_fes3>(a[i],i,rrr,iii))
                CompileError("Copy of Array with incompatible complex vector value FE function () !");;
	  return  new E_set_fev3<Complex,v_fes3>(&b,rr);
         }
	}
       else
	 {  Expression r=0; // new code FH sep 2009.
	     if(!r)  r=Op_CopyArrayT<double,v_fes>(a,b);
	     if(!r)  r=Op_CopyArrayT<Complex,v_fes>(a,b);
	     if(!r)  r=Op_CopyArrayT<double,v_fes3>(a,b);
	     if(!r)  r=Op_CopyArrayT<Complex,v_fes3>(a,b);
	     if(r) return r;
	 }
        CompileError("Internal Error: General Copy of Array : to do ");
      return ret;}

template<class v_fes,int DIM>
C_F0 NewFEvariableT(ListOfId * pids,Block *currentblock,C_F0 & fespacetype,CC_F0 init,bool cplx,int dim)
{
  ffassert(dim==DIM);
  typedef FEbase<double,v_fes> FE;
  typedef E_FEcomp<R,v_fes> FEi;
  typedef typename FEi::Result FEiR;

  typedef FEbase<Complex,v_fes> CFE;
  typedef E_FEcomp<Complex,v_fes> CFEi;
  typedef typename  CFEi::Result CFEiR;

  Expression fes =fespacetype;

  aType dcltype=atype<FE **>();
  aType cf0type=atype<C_F0>();
  aType rtype=atype<FEiR>();

  if(cplx)
    {
      dcltype=atype<CFE **>();
      rtype=atype<CFEiR>();

    }
  ffassert(pids);
  ListOfId &ids(*pids);

  string str("[");

  const int n=ids.size();
  ffassert(n>0);
  if ( fes->nbitem() != (size_t) n) {
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
     char * name = strcpy(CodeAllocT<char>::New(str.size()+1),str.c_str());
     C_F0 ret;
     // modif  100109 (add Block::  before NewVar for g++ 3.3.3 on Suse 9)
    ret= binit ? currentblock->Block::NewVar<LocalVariable>(name,dcltype,basicAC_F0_wa(fespacetype,init))
      : currentblock->Block::NewVar<LocalVariable>(name,dcltype,basicAC_F0_wa(fespacetype));
     C_F0 base = currentblock->Find(name);
     if (cplx)
       for (int i=0;i<n;i++)
         currentblock->NewID(cf0type,ids[i].id, C_F0(new CFEi(base,i,n), rtype) );
     else
       for (int i=0;i<n;i++)
         currentblock->NewID(cf0type,ids[i].id, C_F0(new FEi(base,i,n), rtype) );
     delete pids; // add FH 25032005

      return ret ;
}


C_F0 NewFEvariable(ListOfId * pids,Block *currentblock,C_F0 & fespacetype,CC_F0 init,bool cplx,int dim)
{
  if(dim==2)
    return NewFEvariableT<v_fes,2>(pids,currentblock,fespacetype,init,cplx,dim);
  else if  (dim==3)
      return NewFEvariableT<v_fes3,3>(pids,currentblock,fespacetype,init,cplx,dim);
  else if  (dim==4)
      return NewFEvariableT<v_fesS,4>(pids,currentblock,fespacetype,init,cplx,dim);
  else if  (dim==5)
      return NewFEvariableT<v_fesL,5>(pids,currentblock,fespacetype,init,cplx,dim);
  else
    CompileError("Invalide fespace on Rd  ( d != 2 or 3) ");
    return C_F0();
}
C_F0 NewFEvariable(const char * id,Block *currentblock,C_F0 & fespacetype,CC_F0 init,bool cplx,int dim)
{
  ListOfId * lid =new ListOfId;
  lid->push_back(UnId(id));
  return NewFEvariable(lid,currentblock,fespacetype,init,cplx,dim);
}


size_t dimFESpaceImage(const basicAC_F0 &args)
{
  aType t_tfe= atype<TypeOfFE*>();
  aType t_tfe3= atype<TypeOfFE3*>();
  aType t_tfeS= atype<TypeOfFES*>();
  aType t_tfeL= atype<TypeOfFEL*>();
  aType t_a= atype<E_Array>();
  size_t dim23=0;

  for (int i=0;i<args.size();i++)
    if (args[i].left() == t_tfe || args[i].left() == t_tfe3 || args[i].left() == t_tfeS || args[i].left() == t_tfeL)
      dim23 += args[i].LeftValue()->nbitem();
    else if (args[i].left() == t_a)
      {
	const E_Array & ea= *dynamic_cast<const E_Array *>(args[i].LeftValue());
	ffassert(&ea);
	for (int i=0;i<ea.size();i++)
	  if (ea[i].left() == t_tfe || ea[i].left() == t_tfe3|| ea[i].left() == t_tfeS || ea[i].left() == t_tfeL)
            dim23 += ea[i].nbitem();
	  else ffassert(0); // bug
      }
  dim23 = dim23 ? dim23 : 1;
  return dim23;
}

aType  typeFESpace(const basicAC_F0 &args)
{

  aType t_m2= atype<pmesh*>();
  aType t_m3= atype<pmesh3*>();
  aType t_mS= atype<pmeshS*>();
  aType t_mL= atype<pmeshL*>();

  aType t_tfe= atype<TypeOfFE*>();
  aType t_tfe3= atype<TypeOfFE3*>();
  aType t_tfeS= atype<TypeOfFES*>();
  aType t_tfeL= atype<TypeOfFEL*>();

  aType atfe[]={t_tfe,t_tfe3,t_tfeS,t_tfeL};
  aType atfs[]={atype<pfes *>(),atype<pfes3 *>(),atype<pfesS *>(),atype<pfesL *>()};
  aType t_a= atype<E_Array>();
  aType ret =0,tl=0;
  aType tMesh=0;
  int dm =-1,id=-2;

  for (int i=0;i<args.size();i++) {
    tl=args[i].left();
    if ( tl == t_m2) {ffassert(dm==2 || dm<0); dm=2;}
    else if( tl  == t_m3 ) {ffassert(dm==3 || dm<0); dm=3;}
    else if( tl  == t_mS ) {ffassert(dm==4 || dm<0); dm=4;}
    else if( tl  == t_mL ) {ffassert(dm==5 || dm<0); dm=5;}
    // array
    else if(tl==t_a) {
      const E_Array & ea= *dynamic_cast<const E_Array *>(args[i].LeftValue());
      ffassert(&ea);
      for (int i=0;i<ea.size();i++) {
        tl=ea[i].left();
        for( int it=0; it<4; ++it)
          if (atfe[it]->CastingFrom(tl)) // Warning  P1 can be cast in 2d or 3d FE ...
            id = it;
      }
    }
    else
      for( int it=0; it<4; ++it)
        if (atfe[it]->CastingFrom(tl))
          id = it;
  }

  if (dm==2) ret =atfs[0]; // 2D fespace
  else if (dm==3) ret = atfs[1]; // 3D fespace (Volume)
  else if (dm==4) ret =atfs[2]; // 3D fespace (surface)
  else if (dm==5) ret =atfs[3]; // 3D fespace (curve)
  else {
    cerr << " typeFESpace:: bug dim: maes/EZFv mesh dim :"  << dm << " type FE "<< id +2 <<endl;
    ffassert(0);
  }
  if(verbosity>99) cout << " typeFESpace " << id << " " << ret->name() << endl;
  return ret;
}


template<class v_fes,int DIM>
C_F0 NewFEarrayT(ListOfId * pids,Block *currentblock,C_F0 & fespacetype,CC_F0 sizeofarray,bool cplx,int dim)
{
  ffassert(dim==DIM || dim!=4);  // TODO
  typedef FEbaseArray<double,v_fes>  FE;
  typedef  E_FEcomp<R,v_fes,FE > FEi;
   typedef typename FEi::Result FEiR;

   typedef FEbaseArray<Complex,v_fes>  CFE;
   typedef  E_FEcomp<Complex,v_fes,CFE > CFEi;
   typedef typename CFEi::Result CFEiR;

   Expression fes =fespacetype;
    aType dcltype=atype<FE **>();
    aType cf0type=atype<C_F0>();
    aType rtype=atype<FEiR>();
    ffassert(pids);
    ListOfId &ids(*pids);
    if(cplx)
      {
        dcltype=atype<CFE **>();
        rtype=atype<CFEiR>();

      }

    string str("[");

    const int n=ids.size();
     ffassert(n>0);
   if ( fes->nbitem() != (size_t) n) {
      cerr << " the array size must be " << fes->nbitem()  << " not " <<  n << endl;
      CompileError("Invalid array size  for  vectorial fespace function");
   }
   for (int i=0;i<n;i++)
    {
     str += ids[i].id;
     if(i<n-1) str +=",";
    }
     str += "]";
     char * name = strcpy(CodeAllocT<char>::New(str.size()+1),str.c_str());
     C_F0 ret=  currentblock->Block::NewVar<LocalVariable>(name,dcltype,basicAC_F0_wa(fespacetype,sizeofarray));
     C_F0 base = currentblock->Find(name);
     if(cplx)
      for (int i=0;i<n;i++)
         currentblock->NewID(cf0type,ids[i].id, C_F0(new CFEi(base,i,n), rtype) );
     else
      for (int i=0;i<n;i++)
         currentblock->NewID(cf0type,ids[i].id, C_F0(new FEi(base,i,n), rtype) );

     delete pids ; // add FH 25032005
      return ret ;

}


C_F0 NewFEarray(ListOfId * pids,Block *currentblock,C_F0 & fespacetype,CC_F0 sizeofarray,bool cplx,int dim)
{
  if(dim==2)
    return NewFEarrayT<v_fes,2>(pids,currentblock,fespacetype,sizeofarray,cplx,dim);
  else if  (dim==3)
    return NewFEarrayT<v_fes3,3>(pids,currentblock,fespacetype,sizeofarray,cplx,dim);
  else if  (dim==4)
      return NewFEarrayT<v_fesS,4>(pids,currentblock,fespacetype,sizeofarray,cplx,dim);
 // else if  (dim==5)
 //     return NewFEarrayT<v_fesL,5>(pids,currentblock,fespacetype,sizeofarray,cplx,dim);
  else
    CompileError("Invalid vectorial fespace on Rd  ( d != 2 or 3) ");
    return C_F0();
}
C_F0 NewFEarray(const char * id,Block *currentblock,C_F0 & fespacetype,CC_F0 sizeofarray,bool cplx,int dim)
{
  ListOfId *lid= new ListOfId;
  lid->push_back(UnId(id));
  return NewFEarray(lid,currentblock,fespacetype,sizeofarray,cplx,dim);
}

namespace Fem2D {
void  Expandsetoflab(Stack stack,const BC_set & bc,set<long> & setoflab)
{
    for (size_t i=0;i<bc.on.size();i++)
        if(bc.onis[i] ==0)
	    {
            long  lab  = GetAny<long>( (*bc.on[i])(stack));
            setoflab.insert(lab);
            if ( verbosity>99) cout << lab << " ";

	    }
        else
	    {
            KN<long>  labs( GetAny<KN_<long> >( (*bc.on[i])(stack)));
            for (long j=0; j<labs.N(); ++j) {
                setoflab.insert(labs[j]);
                if ( verbosity>99) cout << labs[j] << " ";
            }

	    }
    if(verbosity>99)
        cout << endl;

}

void  Expandsetoflab(Stack stack,const CDomainOfIntegration & di,set<int> & setoflab,bool &all)
{
    for (size_t i=0;i<di.what.size();i++)
        if(di.whatis[i] ==0)
	    {
            long  lab  = GetAny<long>( (*di.what[i])(stack));
            setoflab.insert(lab);
            if ( verbosity>3) cout << lab << " ";
            all=false;
	    }
        else
	    {
            KN<long>  labs( GetAny<KN_<long> >( (*di.what[i])(stack)));
            for (long j=0; j<labs.N(); ++j) {
                setoflab.insert(labs[j]);
                if ( verbosity>3) cout << labs[j] << " ";
            }
            all=false;
	    }

}
}


#include "InitFunct.hpp"


static addingInitFunct TheaddingInitFunct(-20,init_lgfem);
