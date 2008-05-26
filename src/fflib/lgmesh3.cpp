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
#include  <cmath>
#include <cstdio>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"

#include "RNM.hpp"
#include "fem.hpp"



#include "FESpacen.hpp" 
#include "FESpace.hpp" 

#include "MatriceCreuse_tpl.hpp"

//#include "fem3.hpp"
#include "MeshPoint.hpp"
#include <complex>
#include "Operator.hpp" 

#include <set>
#include <vector>
#include <fstream>

#include "lex.hpp"

#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"


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

template<class R>
AnyType set_fe3 (Stack s,Expression ppfe, Expression e)
{ 
  typedef v_fes3 v_fes;
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;  
  long kkff = Mesh::kfind,  kkth = Mesh::kthrough;
  StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);
  
  
  MeshPoint *mps=MeshPointStack(s),mp=*mps;  
  pair<FEbase<R,v_fes> *,int>  pp=GetAny<pair<FEbase<R,v_fes> *,int> >((*ppfe)(s));
  FEbase<R,v_fes> & fe(*pp.first);
  const  FESpace & Vh(*fe.newVh());
  KN<R> gg(Vh.MaximalNbOfDF()); 
  const  Mesh & Th(Vh.Th);
  //   R F[100]; // buffer 
  TabFuncArg tabexp(s,Vh.N);
  tabexp[0]=e;
  
  if(Vh.N!=1)
    {  cerr << " Try to set a  vectorial  FE function  (nb  componant=" <<  Vh.N << ") with one scalar " << endl;
       ExecError(" Error interploation (set)  FE function (vectorial) with a scalar");
    }
  KN<R> * y=new  KN<R>(Vh.NbOfDF);
  KN<R> & yy(*y);
  // interpoler
  int npPh = Vh.maxNbPtforInterpolation;
  KNM<R>   Vp(npPh,1);
  KN<R>  Vdf(Vh.MaxNbDFPerElement);
  const E_F0 & ff(* (const  E_F0 *) e ) ;
  
  if (Vh.isFEMesh() )
    {
      
      ffassert(Vh.NbOfDF == Th.nv && Vh.N == 1 );
      for (int iv=0;iv<Th.nv;iv++)
	{
	  const Vertex & v(Th(iv));
	  int ik=Th.Contening(&v);
	  const Element & K(Th[ik]);
	  int il=-1;
	  for(int k=0;k<Element::nv;++k)
	    if  ( &K[k] == &v) il=k;
	  assert(il>=0);
	  mps->setP(&Th,ik,il);
	  yy[iv] = GetAny<R>( ff(s) );
	  sptr->clean(); // modif FH mars 2006  clean Ptr
	}
      
    }
  else     
    {
      InterpolationMatrix<RdHat> ipmat(Vh);    
      for (int t=0;t<Th.nt;t++)
	{
	  FElement K(Vh[t]);
	  int nbdf=K.NbDoF() ;	  
	  ipmat.set(t);
	  for (int p=0;p<ipmat.np;p++)
	    { 
	      const RdHat & PtHat(ipmat.P[p]);
	      mps->set(K.T(PtHat),PtHat,K);
	      Vp[p]=GetAny<R>( ff(s) );
	    }
	  K.Pi_h(Vp,Vdf,ipmat);  
	  for (int df=0;df<nbdf;df++)         
	    (*y)[K(df)] =  Vdf[df] ;
	  sptr->clean(); // modif FH mars 2006  clean Ptr	
	}
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


template<class K,class v_fes>
E_set_fev3<K,v_fes>::E_set_fev3(const E_Array * a,Expression pp) 
  :aa(*a),ppfe(pp),optimize(true),
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


template<class K,class v_fes>   
AnyType E_set_fev3<K,v_fes>::operator()(Stack s)  const
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
  
  ffassert(dim<100);
  //   R F[100]; // buffer 
  
  TabFuncArg tabexp(s,Vh.N);
  //   const E_Array * aa = dynamic_cast<const E_Array *>(e);
  ffassert( aa.size() == Vh.N);
  for (int i=0;i<dim;i++)
    tabexp[i]=aa[i]; 
  
  KN<K> * y=new  KN<K>(Vh.NbOfDF);
  KN<K> & yy(*y);
  int npPh = Vh.maxNbPtforInterpolation;
  KNM<K>   Vp(npPh,dim);
  KN<K>  Vdf(Vh.MaxNbDFPerElement);
  
  
  if (Vh.isFEMesh() )
    {
      RdHat KHat[Element::nv];
      for (int i=1 ; i< Element::nv;++i)
	KHat[i+1][i] = 1;
      
      ffassert(Vh.NbOfDF == Th.nv && dim == 1 );
      for (int iv=0;iv<Th.nv;iv++)
	{
	  const E_F0 & ff(* (const  E_F0 *) aa[0]  ) ;
	  const Vertex & v(Th(iv));
	  int ik=Th.Contening(&v);
	  const Element & Kv(Th[ik]);
          int il=-1;
          for(int k=0;k<Element::nv;++k)
            if  ( &Kv[k] == &v) il=k;
          assert(il>=0);
          mps->set(Th,v,KHat[il],Kv,v.lab);
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
    {
       InterpolationMatrix<RdHat> ipmat(Vh);    
       
       for (int t=0;t<Th.nt;t++)
	 {
	   FElement Kt(Vh[t]);
	   int nbdf=Kt.NbDoF();
	   
	   gg=K();
	   
	   for (int p=0;p<ipmat.np;p++)
	     {
	       const RdHat & PtHat(ipmat.P[p]);
	       mps->set(Kt.T(PtHat),PtHat,Kt);
	       
	       if (copt) { // optimize  version 
		 if (optiexpK) (*optiexpK)(s);
		 for (int j=0;j<dim;j++)
		   Vp(p,j) = *(copt[j]);}
	       else  // old version 
		 for (int j=0;j<dim;j++)
		   if (tabexp[j]) 
		     Vp(p,j)=GetAny<K>( (*tabexp[j])(s) );
		   else Vp(p,j)=0;
	       
	     }
           
	   Kt.Pi_h(Vp,Vdf,ipmat);  
	   for (int df=0;df<nbdf;df++)         
	     yy[Kt(df)] =  gg[df] ;
	   
	   sptr->clean(); // modif FH mars 2006  clean Ptr          
	 } 
    }
  fe=y;
  if (copt) delete [] copt;
  *MeshPointStack(s) = mp;
  if(verbosity>1)
    ShowBound(*y,cout) << endl ;
  return Nothing;
}


template<class K>
inline FEbase<K,v_fes> * MakePtrFE3_(pfes3 * const &  a){ 
  FEbase<K,v_fes3> * p=new FEbase<K,v_fes3>(a);
  //cout << "MakePtrFE " << p<< endl; 
  return p ;}
  
template<class K>
inline FEbase<K,v_fes3> ** MakePtrFE3_2(FEbase<K,v_fes3> * * const &  p,pfes3 * const &  a){ 
  *p=new FEbase<K,v_fes3>(a);
  //cout << "MakePtrFE2 " << *p<< endl; 
  return p ;}

template<class K>  
inline FEbaseArray<K,v_fes3> ** MakePtrFE3_3(FEbaseArray<K,v_fes3> * * const &  p,pfes3 * const &  a,const long & N){ 
  *p=new FEbaseArray<K,v_fes3>(a,N);
  //cout << "MakePtrFE2 " << *p<< endl; 
  return p ;}

template<class K,class v_fes>
class  OneOperatorMakePtrFE3 : public OneOperator 
{
public:
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;  

  // il faut Optimize 
  // typedef double K;
  typedef  FEbase<K,v_fes> ** R;
  typedef pfes* B;
  class CODE : public E_F0mps  
  {
  public:
    Expression fer,fes;
    E_set_fev3<K,v_fes> * e_set_fev3;
    const E_Array * v;
    CODE(const basicAC_F0 & args) 
      : 
      fer(to<R>(args[0])),
      fes(to<B>(args[1])),
      e_set_fev3(0) 
    {
      if (BCastTo<K>(args[2]) )
	v = new E_Array(basicAC_F0_wa(to<K>(args[2]))); 
      else 
	v = dynamic_cast<const E_Array *>( args[2].LeftValue() );
      if (!v) {
	cout << "Error: type of arg :" << *args[2].left()  << " in " << typeid(K).name() << " case " << endl;
	ErrorCompile(" We wait  a double/complex expression or a array expression",1);
      }
      e_set_fev3=  new   E_set_fev3<K,v_fes>(v,fer);
      
    }
    
    AnyType operator()(Stack stack)  const {
      R  p = GetAny<R>( (*fer)(stack));
      B  a = GetAny<B>( (*fes)(stack)); 
      *p=new FEbase<K,v_fes>(a);
      (*e_set_fev3)(stack); 
      return SetAny<R>(p);
    }
    operator aType () const { return atype<R>();}         
    
    
  };
  
  E_F0 * code(const basicAC_F0 & args) const 
  { return  new CODE(args);}
  OneOperatorMakePtrFE3(aType tt):  // tt= aType<double>() or aType<E_Array>()  
    OneOperator(map_type[typeid(R).name()],map_type[typeid(R).name()],map_type[typeid(B).name()],tt)
  {}
};


/*
template<class K,class v_fes>    
KN<K> * pf3r2vect( pair<FEbase<K,v_fes> *,int> p)
 {  
    KN<K> * x=p.first->x();
    if ( !x) {  // defined 
      FESpace3 * Vh= p.first->newVh();     
      throwassert( Vh);
      *p.first = x = new KN<K>(Vh->NbOfDF);
      *x=K(); 
    }
    return x;}
*/
template<class K>        
long pf3r_nbdf(pair<FEbase<K,v_fes3> *,int> p)
 {  
   if (!p.first->Vh) p.first->Vh= p.first->newVh();
   throwassert( !!p.first->Vh);
   return p.first->Vh->NbOfDF;
 }

long pVh3_ndof(pfes3 * p)
 { throwassert(p && *p);
   FESpace3 *fes=**p; ;  return fes->NbOfDF ;}
long pVh3_nt(pfes3 * p)
 { throwassert(p && *p);
   FESpace3 *fes=**p; ;  return fes->NbOfElements ;}
long pVh3_ndofK(pfes3 * p)
 { throwassert(p && *p);
   FESpace3 *fes=**p;   return (*fes)[0].NbDoF() ;}

template<class R,int dd,class v_fes>
AnyType pf3r2R(Stack s,const AnyType &a)
{
  typedef typename  v_fes::pfes pfes;
  typedef typename  v_fes::FESpace FESpace;
  typedef typename  FESpace::Mesh Mesh;
  typedef typename  FESpace::FElement FElement;
  typedef typename  Mesh::Element Element;
  typedef typename  Mesh::Vertex Vertex;  
  typedef typename  Mesh::RdHat RdHat;  
  typedef typename  Mesh::Rd Rd;  
  
  pair< FEbase<R,v_fes> *  ,int> ppfe=GetAny<pair< FEbase<R,v_fes> *,int> >(a);
  FEbase<R,v_fes> & fe( *ppfe.first);
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
  const Element *K;
  RdHat PHat;
  bool outside=false;
  bool qnu=true;
  if ( mp.Th3 == &Th && mp.T) 
   {
     qnu=false;
     K=mp.T3;
     PHat=mp.PHat;
   }
  else if ( mp.other.Th3 == & Th && mp.other.P.x == mp.P.x && mp.other.P.y == mp.P.y )
    {
      K=mp.other.T3;
      PHat=mp.other.PHat.p2();
      outside = mp.other.outside;
    } 
  else {
    if (mp.isUnset()) ExecError("Try to get unset x,y, ...");
    K=Th.Find(mp.P,PHat,outside);
    mp.other.set(Th,mp.P.p2(),PHat,*K,0,outside);
  }
  // cout << "  ---  " << qnu << "  " << mp.P << " " << mp.outside <<  " " << outside << endl;
  const FElement KK(Vh[Th(K)]);
  if (outside && KK.tfe->discontinue) 
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


template<class K,class v_fes>
class Op4_pf32K : public quad_function<pair<FEbase<K,v_fes> *,int>,R,R,R,K> { public:
    
    
    class Op : public E_F0mps { public:
	Expression a,b,c,d;
      Op(Expression aa,Expression bb,Expression cc,Expression dd) 
	: a(aa),b(bb),c(cc),d(dd) {}       
      AnyType operator()(Stack s)  const 
      { 
	R xx(GetAny<R>((*b)(s)));
	R yy(GetAny<R>((*c)(s)));
	R zz(GetAny<R>((*d)(s)));
	MeshPoint & mp = *MeshPointStack(s),mps=mp;
	mp.set(xx,yy,zz);
	AnyType ret = pf3r2R<K,0,v_fes>(s,(*a)(s));
	mp=mps;
	return  ret;}
      
    };
};

template<class K,class v_fes>    
KN<K> * pf3r2vect( pair<FEbase<K,v_fes> *,int> p)
{  
  typedef typename  v_fes::FESpace FESpace;
  KN<K> * x=p.first->x();
  if ( !x) {  // defined 
    FESpace * Vh= p.first->newVh();     
    throwassert( Vh);
    *p.first = x = new KN<K>(Vh->NbOfDF);
    *x=K(); 
  }
  return x;}



void init_lgmesh3() {
   if(verbosity)  cout <<"lg_mesh3 ";

    //   Global.Add("buildmesh","(",new OneOperatorCode<classBuildMesh3>);
    // Global.Add("buildmesh","(",new OneOperatorCode<BuildMeshFile3>);

 atype<pmesh3>()->AddCast( new E_F1_funcT<pmesh3,pmesh3*>(UnRef<pmesh3 >)); 
 atype<pfes3 >()->AddCast(  new E_F1_funcT<pfes3,pfes3*>(UnRef<pfes3>));
 
 atype<pf3rbase>()->AddCast(  new E_F1_funcT<pf3rbase,pf3rbase>(UnRef<pf3rbase>));
 atype<pf3cbase>()->AddCast(  new E_F1_funcT<pf3cbase,pf3cbase>(UnRef<pf3cbase>));
 
 Add<pf3r>("[]",".",new OneOperator1<KN<double> *,pf3r>(pf3r2vect<R,v_fes3>));
 Add<pf3c>("[]",".",new OneOperator1<KN<Complex> *,pf3c>(pf3r2vect<Complex,v_fes3>));
 Add<pf3r>("(","",new OneQuadOperator<Op4_pf32K<R,v_fes3>,Op4_pf32K<R,v_fes3>::Op> );
 Add<pf3c>("(","",new OneQuadOperator<Op4_pf32K<Complex,v_fes3>,Op4_pf32K<Complex,v_fes3>::Op> );
 Add<double>("(","",new OneQuadOperator<Op4_K2R<R>,Op4_K2R<R>::Op> );
// Add<long>("(","",new OneTernaryOperator<Op3_K2R<long>,Op3_K2R<long>::Op> ); // FH stupide 
 Add<Complex>("(","",new OneQuadOperator<Op4_K2R<Complex>,Op4_K2R<Complex>::Op> );
 // Add<pmesh3 *>("(","",new OneTernaryOperator<Op3_Mesh2mp,Op3_Mesh2mp::Op> );
 
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
   Add<GlgElement<Mesh3> >("mesure",".",new OneOperator1_<double,GlgElement<Mesh3> >(getmes));
   Add<pmesh3*>("mesure",".",new OneOperator1<double,pmesh3*>(pmesh_mes));
   Add<pmesh3*>("bordermesure",".",new OneOperator1<double,pmesh3*>(pmesh_mesb));
   Add<pmesh3*>("nt",".",new OneOperator1<long,pmesh3*>(pmesh_nt));
   Add<pmesh3*>("nv",".",new OneOperator1<long,pmesh3*>(pmesh_nv));
   Add<pmesh3*>("nbe",".",new OneOperator1<long,pmesh3*>(pmesh_nbe));

 TheOperators->Add("<-",
       new OneOperator2_<pf3rbase*,pf3rbase*,pfes3* >(MakePtrFE3_2),
       new OneOperator3_<pf3rbasearray*,pf3rbasearray*,pfes3*,long >(MakePtrFE3_3),  


       new OneOperator2_<pf3cbase*,pf3cbase*,pfes3* >(MakePtrFE3_2),
       new OneOperator3_<pf3cbasearray*,pf3cbasearray*,pfes3*,long >(MakePtrFE3_3) //,
     //  new OneOperator2_<pmesharray*,pmesharray*,long >(MakePtr)
       
       
       );
 TheOperators->Add("<-",
		   new OneOperatorMakePtrFE3<double,v_fes3>(atype<double>()),  //  scalar case
		   new OneOperatorMakePtrFE3<double,v_fes3>(atype<E_Array>()),  //  vect case
		   new OneOperatorMakePtrFE3<Complex,v_fes3>(atype<Complex>()),  //  scalar complex  case
		   new OneOperatorMakePtrFE3<Complex,v_fes3>(atype<E_Array>())  //  vect complex case
       );
 TheOperators->Add("<-",
       new OneOperator2_<pfes3*,pfes3*,pfes3>(&set_copy_incr));

 TheOperators->Add("=",
       new OneOperator2_<pf3r,pf3r,double,E_F_StackF0F0opt2<double> >(set_fe3<double>) ,
       new OneOperator2_<pf3c,pf3c,Complex,E_F_StackF0F0opt2<Complex> >(set_fe3<Complex>)        
       ) ;     


 map_type[typeid(double).name()]->AddCast(
   new E_F1_funcT<double,pf3r>(pf3r2R<R,0,v_fes3>)
   );
   
 map_type[typeid(Complex).name()]->AddCast(
   new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,0,v_fes3>)
   );

 Global.Add("dz","(",new OneOperatorCode<CODE_Diff<Ftest,op_dz> >);
 Global.Add("dxz","(",new OneOperatorCode<CODE_Diff<Ftest,op_dxz> >);
 Global.Add("dyz","(",new OneOperatorCode<CODE_Diff<Ftest,op_dyz> >);
 Global.Add("dzx","(",new OneOperatorCode<CODE_Diff<Ftest,op_dzx> >);
 Global.Add("dzx","(",new OneOperatorCode<CODE_Diff<Ftest,op_dzy> >);
 Global.Add("dzz","(",new OneOperatorCode<CODE_Diff<Ftest,op_dzz> >);


 Global.Add("dz","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dz> >);
 Global.Add("dxz","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dxz> >);
 Global.Add("dyz","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dyz> >);
 Global.Add("dzx","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dzx> >);
 Global.Add("dzx","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dzy> >);
 Global.Add("dzz","(",new OneOperatorCode<CODE_Diff<Finconnue,op_dzz> >);


 
// bof  
 Global.Add("dx","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dx,v_fes3>));
 Global.Add("dy","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dy,v_fes3>));
 Global.Add("dz","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dz,v_fes3>));
 Global.Add("dxx","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dxx,v_fes3>));
 Global.Add("dyy","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dyy,v_fes3>));
 Global.Add("dxy","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dxy,v_fes3>));
 Global.Add("dyx","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dyx,v_fes3>));
 Global.Add("dzx","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dzx,v_fes3>));
 Global.Add("dzy","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dzy,v_fes3>));
 Global.Add("dzz","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dzz,v_fes3>));
 Global.Add("dxz","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dxz,v_fes3>));
 Global.Add("dyz","(",new E_F1_funcT<double,pf3r>(pf3r2R<R,op_dyz,v_fes3>));



 Global.Add("dx","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dx,v_fes3>));
 Global.Add("dy","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dy,v_fes3>));
 Global.Add("dz","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dz,v_fes3>));
 Global.Add("dxx","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dxx,v_fes3>));
 Global.Add("dyy","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dyy,v_fes3>));
 Global.Add("dxy","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dxy,v_fes3>));
 Global.Add("dyx","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dyx,v_fes3>));
 Global.Add("dzx","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dzx,v_fes3>));
 Global.Add("dzy","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dzy,v_fes3>));
 Global.Add("dzz","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dzz,v_fes3>));
 Global.Add("dxz","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dxz,v_fes3>));
 Global.Add("dyz","(",new E_F1_funcT<Complex,pf3c>(pf3r2R<Complex,op_dyz,v_fes3>));


 Global.Add("int3d","(",new OneOperatorCode<CDomainOfIntegration3d>);
 Global.Add("int2d","(",new OneOperatorCode<CDomainOfIntegrationBorder3d>);
 Global.Add("intallfaces","(",new OneOperatorCode<CDomainOfIntegrationAllFaces>);


 /*
 Add<pfer>("n",".",new OneOperator1<long,pfer>(pfer_nbdf<R>));
 Add<pfec>("n",".",new OneOperator1<long,pfec>(pfer_nbdf<Complex>));
 Add<pmesh*>("area",".",new OneOperator1<double,pmesh*>(pmesh_area));
 Add<pmesh*>("nt",".",new OneOperator1<long,pmesh*>(pmesh_nt));
 Add<pmesh*>("nv",".",new OneOperator1<long,pmesh*>(pmesh_nv));
 Add<pfes*>("ndof",".",new OneOperator1<long,pfes*>(pVh_ndof));
 Add<pfes*>("nt",".",new OneOperator1<long,pfes*>(pVh_nt));
 Add<pfes*>("ndofK",".",new OneOperator1<long,pfes*>(pVh_ndofK));
 Add<pfes*>("(","", new OneTernaryOperator<pVh_ndf,pVh_ndf::Op>  );
 */
}
//#include "InitFunct.hpp"
//static addingInitFunct TheaddingInitFunct(-10,init_lgmesh);
template E_set_fev3<double,v_fes3>::E_set_fev3(const E_Array * a,Expression pp) ;
template E_set_fev3<Complex,v_fes3>::E_set_fev3(const E_Array * a,Expression pp) ;
