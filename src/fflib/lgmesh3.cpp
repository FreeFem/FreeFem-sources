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
#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"

#include "RNM.hpp"
#include "fem.hpp"



#include "FESpacen.hpp" 
#include "FESpace.hpp" 

//#include "fem3.hpp"
#include "MeshPoint.hpp"
#include <complex>
#include "Operator.hpp" 

#include <set>
#include <vector>
#include <fstream>

#include "lex.hpp"

#include "lgfem.hpp"
using Fem2D::Mesh;
using Fem2D::MeshPoint;

extern bool NoWait; 

typedef Mesh * pmesh;
typedef Mesh3 * pmesh3;

 
template<class Mesh> 
class GlgVertex {
public:
  typedef double R;
  typedef typename Mesh::Rd Rd;
  typedef typename Mesh::Vertex Vertex;
  CountPointer<Mesh> pTh;
  const Vertex *v;
  void Check() const {   if (!v || !pTh) { ExecError("Too bad! Unset Vertex!"); } }
  void init() { v=0;pTh.init();}
  GlgVertex(Mesh * Th,long kk): pTh(Th),v( &(*pTh)(kk)) {}
  GlgVertex(Mesh * Th,const Vertex * kk): pTh(Th),v(kk) {}
  operator int() const { Check(); return (* pTh)(v);} 
  operator Rd*(){ Check(); return v;} 
  R x() const {Check() ; return v->X();}
  R y() const {Check() ; return v->Y();}
  R z() const {Check() ; return v->Z();}
  long lab() const {Check() ; return v->lab;}
  void destroy()  {pTh.destroy();}
};

template<class Mesh>
class GlgElement { public:
    CountPointer<Mesh> pTh;
  typedef typename Mesh::Element Element;
  const Element *k;
  
  GlgElement():  k(0) {}
  void  Check() const  {   if (!k || !pTh) { ExecError("Unset Triangle,Sorry!"); } }
  void init() { k=0;pTh.init();}
  void destroy() {pTh.destroy();}
  GlgElement(Mesh * Th,long kk): pTh(Th),k( &(*pTh)[kk]) {}
  GlgElement(Mesh * Th,Element * kk): pTh(Th),k(kk) {}
  operator int() const { Check(); return (* pTh)(k);} 
  GlgVertex<Mesh> operator [](const long & i) const { Check(); return GlgVertex<Mesh>(pTh,&(*k)[i]);}   
  long lab() const {Check() ; return k ? k->lab : 0;}
  double mes() const {Check() ; return k->mesure() ;}
  long n() const { return k ? Element::nv: 0 ;}

};


GlgElement<Mesh3> get_element(pmesh3 const & a, long const & n) {  return GlgElement<Mesh3>(a,n);}
GlgElement<Mesh3> get_element(pmesh3 *const & a, long const & n) {  return GlgElement<Mesh3>(*a,n);}

GlgVertex<Mesh3> get_vertex(pmesh3 const & a, long const & n){ return GlgVertex<Mesh3>(a,n);}
GlgVertex<Mesh3> get_vertex(pmesh3 *const & a, long const & n){ return GlgVertex<Mesh3>(*a,n);}
GlgVertex<Mesh3> get_element(GlgElement<Mesh3> const & a, long const & n) {  return a[n];}

R getx(GlgVertex<Mesh3> const & a){  return a.x();}
R gety(GlgVertex<Mesh3> const & a){  return a.y();}
R getz(GlgVertex<Mesh3> const & a){  return a.z();}
long  getlab(GlgVertex<Mesh3> const & a){  return a.lab();}
long getlab(GlgElement<Mesh3> const & a){  return a.lab();}
R getmes(GlgElement<Mesh3> const & a){  return a.mes();}

double pmesh_mes(pmesh3 * p) { ffassert(p && *p) ;  return (**p).mes ;}
double pmesh_mesb(pmesh3 * p) { ffassert(p && *p) ;  return (**p).mesb;}
long pmesh_nt(pmesh3 * p) { ffassert(p && *p) ;  return (**p).nt ;}
long pmesh_nv(pmesh3 * p) { ffassert(p && *p) ;  return (**p).nv ;}
long pmesh_nbe(pmesh3 * p) { ffassert(p && *p) ;  return (**p).nbe ;}



class MoveMesh3 :  public E_F0mps { public:  
 
   typedef pmesh  Result;
   Expression getmesh;
   Expression U,V;
   int nbsol;    
    vector<Expression> sol;
   
    MoveMesh3(const basicAC_F0 & args) :nbsol(args.size()-2),sol(args.size()-2)
    {   
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]); 
      const E_Array * a = dynamic_cast<const E_Array *>(args[1].LeftValue());
      
      ffassert(a);
      if (a->size() !=2) CompileError("movemesh(Th,[u,v],...) need 2 componate in array ",atype<pmesh>());
      U=to<double>( (*a)[0]);
      V=to<double>( (*a)[1]);
      
      for (int i=2;i<args.size();i++)
        sol[i-2]=to<double>(args[i]);      
    }   
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<E_Array>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new MoveMesh3(args);} 
    AnyType operator()(Stack s) const ;
      operator aType () const { return atype<Result>();} 

};

class SaveMesh3 :  public E_F0 { public:  
 
   typedef pmesh  Result;
   Expression getmesh;
   Expression filename; 
   Expression xx,yy,zz;  
   SaveMesh3(const basicAC_F0 & args) 
    {   
      xx=0;
      yy=0;
      zz=0;
      args.SetNameParam();
      getmesh=to<pmesh3>(args[0]); 
      filename=to<string*>(args[1]); 
      if (args.size() >2) 
        {
          const E_Array * a = dynamic_cast<const E_Array *>(args[2].LeftValue());
          if (!a) CompileError("savemesh(Th,\"filename\",[u,v,w],...");
          int k=a->size() ;
         // cout << k << endl;
          if ( k!=2 && k !=3) CompileError("savemesh(Th,\"filename\",[u,v,w]) need 2 or 3  componate in array ",atype<pmesh>());
          xx=to<double>( (*a)[0]);
          yy=to<double>( (*a)[1]);
          if(k==3)
           zz=to<double>( (*a)[2]);
         }
      
   }   
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh3>(),atype<string*>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new SaveMesh3(args);} 
    AnyType operator()(Stack s) const ;
  
};


AnyType SaveMesh3::operator()(Stack stack) const 
{
  using  Fem2D::MeshPointStack;
  
  
   pmesh3 Thh = GetAny<pmesh3>((*getmesh)(stack));
   string * fn =  GetAny<string*>((*filename)(stack));
   cout << "SaveMesh3 " << *fn << " " << Thh << endl;
   Thh->Save(*fn);
   return SetAny<pmesh3>(Thh);

}

AnyType MoveMesh3::operator()(Stack stack) const 
{
  ffassert(0);
    /*
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
   MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
   Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   ffassert(Thh);
   long nbv=Thh->nv;
   long nbt=Thh->nt;
   KN<double> u(nbv),v(nbv);
   double infini=DBL_MAX;
   u=infini;
   for (int it=0;it<nbt;it++)
    for (int iv=0;iv<3;iv++)
    {
      int i=(*Thh)(it,iv);
      if ( u[i]==infini) { // if nuset the set 
        mp->setP(Thh,it,iv);
        u[i]=GetAny<double>((*U)(stack));
        v[i]=GetAny<double>((*V)(stack));
      }
    }
    
   Mesh * pth= MoveTheMesh(*Thh,u,v);
   if (pth)
     for (size_t i=0;i<sol.size();i++)
       { //  ale 
          pair<FEbase<double>,int> * s = GetAny<pair<FEbase<double>,int>*>( (*sol[i])(stack));
          ffassert(s->first.Vh);
          ffassert( &s->first.Vh->Th == Thh); // same old mesh
          ffassert(0); // a faire ????
       }
   *mp=mps;
    pth->decrement();   
    return SetAny<pmesh>(pth);
    */
}


inline pmesh3 *  initMesh(pmesh3 * const & p, string * const & s) {
  Mesh3 * m;
  cout << " initMesh " << *s << endl;
  *p= m =new Mesh3(*s); 
  m->BuildGTree();
 //  delete s;  modif mars 2006 auto del ptr 
  return p;
 }
/*
class CheckMoveMesh :  public E_F0mps { public:  
 
   typedef double  Result;
   Expression getmesh;
   Expression U,V;
   int nbsol;    
    vector<Expression> sol;
   
    CheckMoveMesh(const basicAC_F0 & args) :nbsol(args.size()-2),sol(args.size()-2)
    {   
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]); 
      const E_Array * a = dynamic_cast<const E_Array *>(args[1].LeftValue());
      
      ffassert(a);
      if (a->size() !=2) CompileError("CheckMoveMesh(Th,[u,v]) need 2 componate in array ",atype<pmesh>());
      U=to<double>( (*a)[0]);
      V=to<double>( (*a)[1]);
      
      for (int i=2;i<args.size();i++)
        sol[i-2]=to<double>(args[i]);      
    }   
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<E_Array>(),false);}
    static  E_F0 * f(const basicAC_F0 & args){ return new CheckMoveMesh(args);} 
    AnyType operator()(Stack s) const ;
    operator aType () const { return atype<double>();}         
  
};
AnyType CheckMoveMesh::operator()(Stack stack) const 
{
 
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
   MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
   Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   Mesh & Th(*Thh);
   ffassert(Thh);
   long nbv=Thh->nv;
   long nbt=Thh->nt;
   KN<double> u(nbv),v(nbv);
   double infini=DBL_MAX;
   u=infini;
   for (int it=0;it<nbt;it++)
    for (int iv=0;iv<3;iv++)
    {
      int i=(*Thh)(it,iv);
      if ( u[i]==infini) { // if nuset the set 
        mp->setP(Thh,it,iv);
        u[i]=GetAny<double>((*U)(stack));
        v[i]=GetAny<double>((*V)(stack));
      }
    }
     double minarea=DBL_MAX;
    for (int t=0;t<Th.nt;t++)
     {
      int i0=Th(t,0),i1=Th(t,1),i2=Th(t,2);
      minarea=Min(minarea,Area2(R2(u[i0],v[i0]), R2(u[i1],v[i1]),R2(u[i2],v[i2])));
     }
    *mp=mps;
    return SetAny<double>(minarea/2.);

}
*/
void init_lgmesh3() {
   if(verbosity)  cout <<"lg_mesh3 ";

    //   Global.Add("buildmesh","(",new OneOperatorCode<classBuildMesh3>);
    // Global.Add("buildmesh","(",new OneOperatorCode<BuildMeshFile3>);

   TheOperators->Add("<-",
       new OneOperator2_<pmesh3*,pmesh3*,string* >(&initMesh));
       
// use for :   mesh Th = readmesh ( ...);       
  TheOperators->Add("<-",
       new OneOperator2_<pmesh3*,pmesh3*,pmesh3 >(&set_copy_incr));

   Global.Add("savemesh","(",new OneOperatorCode<SaveMesh3>);

   Dcl_Type<GlgVertex<Mesh3> >(); 
   Dcl_Type<GlgElement<Mesh3> >( ); 

   atype<long>()->AddCast( 
			  new E_F1_funcT<long,GlgVertex<Mesh3> >(Cast<long,GlgVertex<Mesh3> >),
			  new E_F1_funcT<long,GlgElement<Mesh3> >(Cast<long,GlgElement<Mesh3> >)
  );

   Add<pmesh3>("[","",new OneOperator2_<GlgElement<Mesh3>,pmesh3,long>(get_element));
   Add<pmesh3*>("[","",new OneOperator2_<GlgElement<Mesh3>,pmesh3*,long>(get_element));
   Add<pmesh3>("(","",new OneOperator2_<GlgVertex<Mesh3>,pmesh3,long>(get_vertex));
   Add<pmesh3*>("(","",new OneOperator2_<GlgVertex<Mesh3>,pmesh3*,long>(get_vertex));

   Add<GlgElement<Mesh3> >("[","",new OneOperator2_<GlgVertex<Mesh3> ,GlgElement<Mesh3> ,long>(get_element));
   Add<GlgVertex<Mesh3> >("x",".",new OneOperator1_<R,GlgVertex<Mesh3> >(getx));
   Add<GlgVertex<Mesh3> >("y",".",new OneOperator1_<R,GlgVertex<Mesh3> >(gety));
   Add<GlgVertex<Mesh3> >("z",".",new OneOperator1_<R,GlgVertex<Mesh3> >(getz));
   Add<GlgVertex<Mesh3> >("label",".",new OneOperator1_<long,GlgVertex<Mesh3> >(getlab));
   Add<GlgElement<Mesh3> >("label",".",new OneOperator1_<long,GlgElement<Mesh3> >(getlab));
   Add<GlgElement<Mesh3> >("region",".",new OneOperator1_<long,GlgElement<Mesh3> >(getlab));
   Add<GlgElement<Mesh3> >("mes",".",new OneOperator1_<double,GlgElement<Mesh3> >(getmes));
   Add<pmesh3*>("mes",".",new OneOperator1<double,pmesh3*>(pmesh_mes));
   Add<pmesh3*>("bordermes",".",new OneOperator1<double,pmesh3*>(pmesh_mesb));
   Add<pmesh3*>("nt",".",new OneOperator1<long,pmesh3*>(pmesh_nt));
   Add<pmesh3*>("nv",".",new OneOperator1<long,pmesh3*>(pmesh_nv));
   Add<pmesh3*>("nbe",".",new OneOperator1<long,pmesh3*>(pmesh_nbe));
   
}
//#include "InitFunct.hpp"
//static addingInitFunct TheaddingInitFunct(-10,init_lgmesh);
