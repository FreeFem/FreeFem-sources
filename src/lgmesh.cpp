#include  <iostream>
#include  <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"

#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp" 

//#include "fem3.hpp"
#include "MeshPoint.hpp"
#include <complex>
#include "Operator.hpp" 

#include <set>
#include <vector>
#include <fstream>

#include "lex.hpp"

#include "Mesh2.h"

#include "BamgFreeFem.hpp"
#include "lgfem.hpp"
using Fem2D::Mesh;
using Fem2D::MeshPoint;
extern bool NoWait; 

typedef Mesh * pmesh;
class MoveMesh :  public E_F0mps { public:  
 
   typedef pmesh  Result;
   Expression getmesh;
   Expression U,V;
   int nbsol;    
    vector<Expression> sol;
   
    MoveMesh(const basicAC_F0 & args) :nbsol(args.size()-2),sol(args.size()-2)
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
      
      throwassert(a);
      if (a->size() !=2) CompileError("movemesh(Th,[u,v],...) need 2 componate in array ",atype<pmesh>());
      U=to<double>( (*a)[0]);
      V=to<double>( (*a)[1]);
      
      for (int i=2;i<args.size();i++)
        sol[i-2]=to<double>(args[i]);      
    }   
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<E_Array>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new MoveMesh(args);} 
    AnyType operator()(Stack s) const ;
      operator aType () const { return atype<Result>();} 

};

class SplitMesh :  public E_F0mps { public:  
 
   typedef pmesh  Result;
   Expression getmesh;
   Expression U;
   
    SplitMesh(const basicAC_F0 & args) 
    {   
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]); 
      U = to<long>(args[1]);      
     }   
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<long>());}
    static  E_F0 * f(const basicAC_F0 & args){ return new SplitMesh(args);} 
    AnyType operator()(Stack s) const ;
    operator aType () const { return atype<Result>();} 

};

class SaveMesh :  public E_F0 { public:  
 
   typedef pmesh  Result;
   Expression getmesh;
   Expression filename; 
   Expression xx,yy,zz;  
   SaveMesh(const basicAC_F0 & args) 
    {   
      xx=0;
      yy=0;
      zz=0;
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]); 
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
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<string*>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new SaveMesh(args);} 
    AnyType operator()(Stack s) const ;
  
};


class Adaptation :   public E_F0mps { public:
    typedef pmesh  Result;

   static basicAC_F0::name_and_type name_param[] ;
   static const int n_name_param =26;
   
   int nbsol;    
    Expression nargs[n_name_param];
    Expression getmesh;
    Expression em11,em22,em12;
    int  typesol[100];
    vector<Expression> sol;
    
    double arg(int i,Stack stack,double a) const { return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): a;}
    long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
    bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
    
    Adaptation(const basicAC_F0 & args) :nbsol(args.size()-1),sol(args.size()-1)
    {   
      em11=0;
      em22=0;
      em12=0;
      
      args.SetNameParam(n_name_param,name_param,nargs);
      getmesh=to<pmesh>(args[0]); 
      int ksol=0; 
      throwassert(nbsol<100);
      for (int i=1;i<nbsol+1;i++)       
         if (args[i].left()==atype<E_Array>())
          {
            const E_Array * a = dynamic_cast<const E_Array *>(args[i].LeftValue());
            throwassert(a);
            ksol+=a->size(); 
          }
         else
           ksol++;
      sol.resize(ksol); 
      ksol=0; 
      for (int i=1;i<nbsol+1;i++)       
         if (args[i].left()==atype<E_Array>())
          {
            const E_Array * a = dynamic_cast<const E_Array *>(args[i].LeftValue());
            throwassert(a);
             int N=a->size();
            typesol[i-1]=N-1; // ok en 2D
            if (N<=4) {
             for (int j=0;j<N;j++)             
              sol[ksol++]=to<double>((*a)[j]); }
            else {
             lgerror(" Adaptation vecteur a plus de 4 compossantes inconnue");
            }
              
          }
         else {
          typesol[i-1]=0;
          sol[ksol++]=to<double>(args[i]);
          }
    const E_Array * expmetrix = dynamic_cast<const E_Array *>(nargs[24]);
  
  if(expmetrix)
   {
      if(expmetrix->nbitem()!=3)
        ExecError("\nSorry we wait an array with 3 componants in: metrix=[m11,m12,m22]");
        
      em11=(*expmetrix)[0];
      em12=(*expmetrix)[1];
      em22=(*expmetrix)[2];
      
      if(  (*expmetrix)[0].left()!= atype<KN<double> *>() )
          CompileError("Sorry the fist  array componant in metrix=[m11,m12,m22] must be vector");
      if( (*expmetrix)[1].left()!= atype<KN<double> *>() )
          CompileError("Sorry the second  array componant in metrix=[m11,m12,m22] must be vector");
      if( (*expmetrix)[2].left()!= atype<KN<double> *>() ) 
          CompileError("Sorry the third  array componant in metrix=[m11,m12,m22] must be vector");
          
          
     }
   }  
    
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new Adaptation(args);} 
    AnyType operator()(Stack s) const ;    
    operator aType () const { return atype<pmesh>();}         
};


 basicAC_F0::name_and_type Adaptation::name_param[Adaptation::n_name_param] = {
          "hmin",             &typeid(double),
          "hmax",             &typeid(double),
          "err",              &typeid(double), 
          "errg",             &typeid(double), 
          "nbvx",             &typeid(long),        // 4
          "nbsmooth",         &typeid(long),
          "nbjacoby",         &typeid(long),
          "ratio",            &typeid(double), 
          "omega",            &typeid(double),        
          "iso",              &typeid(bool),         // 9
          "abserror",         &typeid(bool), 
          "cutoff",           &typeid(double),  
          "verbosity",        &typeid(long),
          "inquire",          &typeid(bool),         
          "splitpbedge",      &typeid(bool), // 14
          "maxsubdiv",        &typeid(double),
          "anisomax",         &typeid(double),
          "rescaling",        &typeid(bool),
          "keepbackvertices", &typeid(bool),
          "IsMetric",         &typeid(bool),    // 19
          "power",            &typeid(double),    // 20 
          "thetamax",         &typeid(double), 
          "splitin2",         &typeid(bool) ,
          "nomeshgeneration", &typeid(bool) ,
          "metric"          ,  &typeid(E_Array),  // 24
          "periodic"        ,  &typeid(E_Array) // 25 
    };

struct Op_trunc_mesh : public OneOperator {

    class Op: public E_F0mps   { public:
      static basicAC_F0::name_and_type name_param[] ;
      static const int n_name_param =2;
      Expression nargs[n_name_param];
    
      Expression getmesh,bbb;
      long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
      Op(const basicAC_F0 &  args,Expression t,Expression b) : getmesh(t),bbb(b) 
        { args.SetNameParam(n_name_param,name_param,nargs); }
      AnyType operator()(Stack s)  const ;
     };
    E_F0 * code(const basicAC_F0 & args) const 
     { return new Op(args,to<pmesh>(args[0]),to<bool>(args[1])) ;}
   Op_trunc_mesh() : 
     OneOperator(atype<pmesh>(),atype<pmesh>(),atype<bool>()) {}     
};

basicAC_F0::name_and_type Op_trunc_mesh::Op::name_param[Op_trunc_mesh::Op::n_name_param] =
 {
     "split",             &typeid(long),
     "label",             &typeid(long),
 
 };


AnyType Op_trunc_mesh::Op::operator()(Stack stack)  const { 
    using namespace    Fem2D;
    Mesh & Th = *GetAny<pmesh>((*getmesh)(stack));
    long kkksplit =arg(0,stack,1L);
    long label =arg(1,stack,2L);
    KN<int> split(Th.nt);
    split=kkksplit;
    MeshPoint *mp= MeshPointStack(stack),mps=*mp;
    long kk=0;
    for (int k=0;k<Th.nt;k++)
     { 
       Triangle & K(Th[k]);
       R2 B(1./3.,1./3.);
       mp->set(Th,K(B),B,K,0);
       if (  GetAny<bool>((*bbb)(stack))  ) kk++;
       else  split[k]=0  ;    
     }
     *mp=mps;
     if (verbosity>2) 
     cout << " -- Trunc mesh: Nb of Triangle = " << kk << " label=" <<label <<endl;
  pmesh pmsh = new Mesh(Th,split,false,label);
  pmsh->renum();
   /* deja fait  dans bamg2msh
  Fem2D::R2 Pn,Px;
  m->BoundingBox(Pn,Px);
  m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);*/
  pmsh->decrement();  
  return SetAny<pmesh>(pmsh); 
 };
 

AnyType SplitMesh::operator()(Stack stack) const 
{
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
   MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
   Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   throwassert(Thh);
   int label=1;
   Mesh & Th(*Thh);
   long nbv=Thh->nv;
   long nbt=Thh->nt;
   KN<int> split(nbt);
   R2 B(1./3.,1./3.);
   for (int it=0;it<nbt;it++)
    {
      Triangle & K(Th[it]);      
      mp->set(Th,K(B),B,K,K.lab);
      split[it]=GetAny<long>((*U)(stack));
    }
    
   Mesh * pth= new Mesh(*Thh,split,false,label);
   R2 Pn,Px;
   pth->BoundingBox(Pn,Px);
   if(!pth->quadtree)
   pth->quadtree=new Fem2D::FQuadTree(pth,Pn,Px,pth->nv);

   *mp=mps;
    pth->decrement();   
    return SetAny<pmesh>(pth);
}
AnyType SaveMesh::operator()(Stack stack) const 
{
  using  Fem2D::MeshPointStack;
   Fem2D::Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   string * fn =  GetAny<string*>((*filename)(stack));
   if (!xx && !yy ) {
     ::bamg::Triangles * bTh= msh2bamg(*Thh);   
     (*bTh).Write(fn->c_str(),::bamg::Triangles::AutoMesh);
     delete bTh;
     }
   else {
    MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
     ofstream fp((*fn+".points").c_str());
     ofstream ff((*fn+".faces").c_str());
     fp.precision(12);
     if (verbosity>1) 
       cout << " -- Opening files " << (*fn+".points") << " and " << (*fn+".faces") << endl;
    const   Fem2D::Mesh & Th=*Thh;
    long nbv=Thh->nv;
    long nbt=Thh->nt;
    ff << nbt << endl;
    fp << nbv << endl;
    KN<long> num(nbv);
    num=-1;
    int k=0;
    for (int it=0;it<nbt;it++)
      {
        const Fem2D::Triangle & K(Th[it]); 
        
        num[Th(K[0])]=k++;
        num[Th(K[1])]=k++;
        num[Th(K[2])]=k++;
        
         ff << " 3 " << Th(K[0])+1 << ' ' << Th(K[1])+1 << ' ' << Th(K[2])+1 << ' ' 
         	<< " 0 0 0 " << K.lab << '\n';
      }
    if( verbosity>5)
     cout << "  - end writing faces  " << endl;
    for (int iv=0;iv<nbv;iv++)
      {
 //       cout << iv << endl;
        const Fem2D::Vertex  & v(Th(iv)); 
        assert( iv == Th(num[iv]/3,num[iv]%3));
        mp->setP(Thh,num[iv]/3,num[iv]%3);
        
        fp << GetAny<double>((*xx)(stack)) << ' ';
        fp << GetAny<double>((*yy)(stack)) << ' ';
        if (zz)
         fp << GetAny<double>((*zz)(stack)) << ' ';
        else 
         fp << " 0 ";
        fp << v.lab<< '\n';
      }
    if( verbosity>5)
     cout << "  - end writing points  " << endl;

   *mp= mps;   
     
   }
   delete fn;
   return SetAny<pmesh>(Thh);
}
AnyType MoveMesh::operator()(Stack stack) const 
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
   throwassert(Thh);
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
     for (int i=0;i<sol.size();i++)
       { //  ale 
          pair<FEbase<R>,int> * s = GetAny<pair<FEbase<R>,int>*>( (*sol[i])(stack));
          assert(s->first.Vh);
          assert( &s->first.Vh->Th == Thh); // same old mesh
          throwassert(0); // a faire ????
       }
   *mp=mps;
    pth->decrement();   
    return SetAny<pmesh>(pth);

}


AnyType Adaptation::operator()(Stack stack) const 
{
  using namespace bamg;
  using bamg::Min;
  using bamg::Max;
  using bamg::Abs;
  
  int nbcperiodic=0;
  Expression *periodic=0;

  Real8 err         = arg(2,stack,0.01);    // coef in the metric
  Real8  errg       = Min(arg(3,stack,0.01),err);
  long nbsx         = Max(100L,arg(4,stack,9000L));
  long nbsmooth     =  arg(5,stack,3L);
  long nbjacobi     =  arg(6,stack,0L) ;              // if increased will be more smooth
  const Real8 raison = arg(7,stack,1.8);
  const Real8 omega =  arg(8,stack,1.0) ; 
  bool iso          =   arg(9,stack,false);
  bool AbsError     =   arg(10,stack,true);
  Real8 CutOff      =   arg(11,stack, 1.0e-6);
  verbosity         =   arg(12,stack, (long) verbosity);   
  bool inq          =   arg(13,stack,false);
  bool SplitEdgeWith2Boundary =  arg(14,stack,true);
  double maxsubdiv             = Max(Min( arg(15,stack,10.0),10.0),0.1);
  double anisomax              =   Max((double) arg(16,stack,1.0e6),1.0);
  bool rescaling               = arg(17,stack,true) ;
  bool KeepBackVertices        =   arg(18,stack,true) ;
  int givenmetric              =  arg(19,stack,false) ;
  double powerM                = arg(20,stack,1.0) ;
  double cutoffradian          = arg(21,stack,-1.0)* bamg::Pi/180. ;
  bool split                    = arg(22,stack,false) ;
  bool nomeshgeneration         = arg(23,stack,false) ;
  const E_Array * expmetrix = dynamic_cast<const E_Array *>(nargs[24]);
  GetPeriodic(nargs[25],nbcperiodic,periodic);

  KN<double> *mm11=0, *mm12=0,* mm22=0;

  using Fem2D::MeshPoint;
  using Fem2D::Mesh;
   Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
  throwassert(Thh);
    Triangles * oTh =0;
  if (nbcperiodic) {
    KN<int> ndfv(Thh->nv);
    KN<int> ndfe(Thh->neb);
    int nbdfv=0,nbdfe=0;      
    BuildPeriodic(nbcperiodic,periodic,*Thh,stack,nbdfv,ndfv,nbdfe,ndfe);
     oTh = msh2bamg(*Thh,cutoffradian,nbdfv,ndfv,nbdfe,ndfe);
    cerr << " Sorry periodic mesh adaptation is not well implemented "<< endl;
    ExecError("adaptmesh( ... )");
  }
  else
   oTh = msh2bamg(*Thh,cutoffradian);
  Triangles &Th(*oTh);
  bool mtx=em11 && em22 && em12;
  if( mtx )
   {
      mm11= GetAny<KN<double> *>( (*em11)(stack) );
      mm22= GetAny<KN<double> *>( (*em22)(stack) );
      mm12= GetAny<KN<double> *>( (*em12)(stack) );
      if (mm11->N() != Th.nbv || mm22->N() != Th.nbv || mm12->N() != Th.nbv)
        ExecError("The size of3  metrics array must be equal to nb of vertex");
      
   }
   
  KN<double> &m11=*mm11;
  KN<double> &m12=*mm12;
  KN<double> &m22=*mm22;
 
  
  Real8 hmax = 0.3*Th.MaximalHmax(); // final largest edge 
  Real8 hmin = Th.MinimalHmin();        // final smallest edge
  
  Real8 coef =1;                // a priori don't touch
  // gestion des arguments 
  hmin              = Max(hmin, arg(0,stack,hmin));
  hmax              = Min(hmax,arg(1,stack,hmax));
  if (inq) 
   if (!withrgraphique) {initgraphique();withrgraphique=true;}
  
    if (iso)  anisomax=1;
  if (verbosity>2) 
    {
      cout << endl  << endl; 
      cout << " \t\t ## adapt : nbsol= " << nbsol << ", nbsx = " << nbsx << ", err = " << err ;
      cout << ", hmin = " << hmin << ", hmax = " << hmax <<endl;
      cout << " \t\t    ratio  = " << raison << ", nbsmooth = " << nbsmooth ;
      cout << ", omega = " << omega <<  ", coef = " << coef << ", iso = " << iso << endl;
      cout << " \t\t    AbsError =" << AbsError << ", CutOff = " << CutOff << ", nbjacobi = " << nbjacobi <<endl;
      cout << " \t\t    maxsubdiv = " << maxsubdiv << " splitpbedge = " << SplitEdgeWith2Boundary  <<endl;
      cout << " \t\t    anisomax = " << anisomax << ", rescaling = " << rescaling << ", power = " << powerM
           << ", KeepBackvertices = " << KeepBackVertices << " IsMetric = " << givenmetric
      << endl << endl ; 
    }
    
 // 
 Th.ReMakeTriangleContainingTheVertex();
  MeshPoint* mp(Fem2D::MeshPointStack(stack));
   
  Int4 i,iv; 
  int ksol =0;
  for (i=0;i<nbsol;i++)
    ksol += typesol[i]+1; //  marche en 2d 
   
  double * lessol = new double [Th.nbv*ksol];
  double *ss = lessol;
  // be careful because renum --
  // the triangle was no renum 
  for ( iv=0;iv<Th.nbv;iv++) 
    Th[iv].color=1; // color 
 for (Int4  it = 0; it < Thh->nt; it++)
    for (Int4  jt = 0; jt < 3; jt++)
      { 
        bamg::Vertex & v= Th(it)[jt];
        const Fem2D::Vertex & vf = (*Thh)[it][jt];
        if (&v && v.color)
          {
            v.color =0; // uncolor
            mp->setP(Thh ,it,jt);
 
            ss = lessol + ksol* Th.Number(v);
            for (int j =0; j < ksol; j++)
              *ss++= GetAny<double>( (*sol[j])(stack) );
           
          }
      }
  mp->unset();
  // computation of the metric --- 
  // better thing -> create keyword in the language 
  //    a faire F Hecht .
  Metric Mhmax(hmax);
  for ( iv=0;iv<Th.nbv;iv++) 
    Th[iv].m = Mhmax;
    
   if (mtx) 
    for ( iv=0;iv<Th.nbv;iv++) 
      if ( ::Max(m11[iv],m12[iv],m22[iv]) > hmax) 
       Th[iv].m = MetricAnIso(m11[iv],m12[iv],m22[iv]);
   
  if ( givenmetric)
    if (ksol == 1) 
      {
        for (Int4  iv = 0,k=0; iv < Th.nbv ; iv++)
          Th[iv].m.IntersectWith(Metric(lessol[k++]));
      }
    else if (ksol == 3) 
      {
        for (Int4  iv = 0,k=0; iv < Th.nbv ; iv++, k += 3)
          {
	    Metric MM(lessol[k],lessol[k+1],lessol[k+2]);
	    MatVVP2x2 vp(MM);
	    vp.Abs();
	    Th[iv].m.IntersectWith(vp);
          }
      }
    else
      lgerror("Adapt mesh: ksol is wrong, IsMetric  and ksol != 1 or 3");
  else
  Th.IntersectConsMetric(lessol,nbsol,typesol,hmin,hmax,sqrt(err)*coef,anisomax,AbsError?0.0:CutOff,nbjacobi,rescaling,powerM,0);
#ifdef DRAWING1
  if ( (inq!=0) ) {
    if (!withrgraphique) {initgraphique();withrgraphique=true;}
    reffecran();
    Th.InitDraw();
    Th.Draw();
    if (!NoWait) Th.inquire();} 
#endif     
  delete [] lessol;
  Th.IntersectGeomMetric(errg,iso);
  Th.SmoothMetric(raison);
  Th.MaxSubDivision(maxsubdiv);
  Th.BoundAnisotropy(anisomax);
  // end of metric's computation 
   if (mtx) 
    for ( iv=0;iv<Th.nbv;iv++) 
      {
       m11[iv] = Th[iv].m.a11  ;
       m22[iv] = Th[iv].m.a22 ; 
       m12[iv] = Th[iv].m.a21;
      }
  
  Triangles* nTh = 0;
  
  if ( ! nomeshgeneration)
  {
  nTh= new Triangles(nbsx,Th,KeepBackVertices); // Adaption is here
  
  if (split)
	    Th.SplitElement(1);
 
  if(SplitEdgeWith2Boundary)
    nTh->SplitInternalEdgeWithBorderVertices();
  
  if(verbosity>2) 
    nTh->ShowHistogram();
  if (nbsmooth)
    nTh->SmoothingVertex(nbsmooth,omega);
  if(verbosity>2 && nbsmooth) 
    nTh->ShowHistogram();
   
#ifdef DRAWING
  if ((inq!=0)) {
    if (!withrgraphique) {initgraphique();withrgraphique=true;}
    
    reffecran();
    nTh->InitDraw();
    nTh->Draw();
    if(!NoWait) nTh->inquire();}
#else
  inq=0;
#endif  
  Metric M(hmax);
  for (iv=0;iv < Th.nbv;iv++)
    Th[iv].m = M;

  Mesh * g=  bamg2msh(nTh,true);

  delete nTh;
  delete oTh;
  g->decrement();
  return SetAny<pmesh>(g);
  }
 else {
   delete oTh;
   return SetAny<pmesh>(Thh);
 }
}   

 
Mesh * MoveTheMesh(const Fem2D::Mesh &Th,const KN_<double> & U,const KN_<double> &V)
{
  using  Fem2D::Triangle;
  using  Fem2D::Vertex;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
  int nbv=Th.nv;
  int nbt=Th.nt;
  int neb=Th.neb;
  Vertex * v= new Vertex[nbv];
  Triangle *t= new Triangle[nbt];
  BoundaryEdge *b= new BoundaryEdge[neb];
  Vertex *vv=v;
  for (int i=0;i<nbv;i++)
   {
     R2 P(U[i],V[i]);     
     vv->x=P.x;
     vv->y=P.y;
     vv->lab = Th(i).lab;
     vv++;      
    }   
  Triangle *tt= t; 
  int nberr=0;
   
  for (int i=0;i<nbt;i++)
    {
      int i0=Th(i,0), i1=Th(i,1),i2=Th(i,2);
      R a= Area2(v[i0],v[i1],v[i2])/2;
      if ( a < Th[i].area/100 ) 
       { nberr++;
        if (verbosity>1) 
         {
          if (nberr==1) cerr << "Erreur: MoveMesh ";
          if (nberr < verbosity*5) {
            cerr << " " <<i;
            if ( nberr % 5 )  cerr << "\n\t";}
         }}
       else 
        (*tt++).set(v,i0,i1,i2,Th[i].lab,a);
    
      if (nberr)
       { if (verbosity) 
         cerr << "Error movemesh: " << nberr << " triangles was reverse  (=> no move)" <<  endl;  
         delete []v;
         delete []t;
         delete []b;   
         throw(ErrorExec("Error move mesh triangles was reverse",1));      
         return 0;
       }
   }  
  BoundaryEdge * bb=b;
  for (int i=0;i<neb;i++)
    {        
     int i1=Th(Th.bedges[i][0]);
     int i2=Th(Th.bedges[i][1]);
     int lab=Th.bedges[i].lab;     
     *bb++ = BoundaryEdge(v,i1,i2,lab);   
    }
 {
  Mesh * m = new Mesh(nbv,nbt,neb,v,t,b);
  R2 Pn,Px;
  m->BoundingBox(Pn,Px);
  m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
//  m->decrement();  
  return m;
      
}
}
Mesh * Carre(int nx,int ny,Expression fx,Expression fy,Stack stack)
{
  using  Fem2D::Vertex;
  using  Fem2D::Triangle;
  using  Fem2D::R2;
  using  Fem2D::BoundaryEdge;
  using  Fem2D::Mesh;
 // using  Fem2D::R;
  using  Fem2D::MeshPointStack;
  int nx1=nx+1,ny1=ny+1;
  int nbv=(nx1)*(ny1);
  int nbt=(nx)*ny*2;
  int neb=(nx+ny)*2;
  Vertex * v= new Vertex[nbv];
  Triangle *t= new Triangle[nbt];
  BoundaryEdge *b= new BoundaryEdge[neb];
  Vertex *vv=v;
  for (int j=0;j<ny1;j++)
   for (int i=0;i<nx1;i++)
   {
     R2 P((R) i/nx,(R) j/ny);
     
     vv->x=P.x;
     vv->y=P.y;
     vv->lab =  0;
     vv++;
   }
  if (fx || fy)
   {
     MeshPoint *mp=MeshPointStack(stack);
    vv=v;
    for (int j=0;j<ny1;j++)
     for (int i=0;i<nx1;i++)
       {
         mp->set(vv->x,vv->y);
         if (fx)  vv->x= GetAny<R>((*fx)(stack));
         if (fy)  vv->y= GetAny<R>((*fy)(stack));
         vv++;
       } 
   }
  Triangle *tt= t;   
  
  for (int j=0;j<ny;j++)
   for (int i=0;i<nx;i++)
     { 
       int i0 = i + j*nx1;
       int i1= i0+1;
       int i2=i1+nx1;
       int i3=i2-1;   
       (tt++)->set(v,i0,i1,i2,0,0.0);
       (tt++)->set(v,i0,i2,i3,0,0.0);
     }  
  BoundaryEdge * bb=b;
  for (int i=0;i<nx;i++)
    {  // bottom 
      int j=0;
      int i1=i,i2=i1+1;
      *bb++ = BoundaryEdge(v,i1,i2,1);
      v[i1].lab=v[i2].lab=1;
    }     
  for (int j=0;j<ny;j++)
    { // right
      int i=nx;
      int i1= i + j*nx1 ,i2=i1+nx1;
      *bb++ = BoundaryEdge(v,i1,i2,2);
      v[i1].lab=v[i2].lab=2;
    }     
    
  for (int i=0;i<nx;i++)
    {  // up
      int j=ny;
      int i1=i + j*nx1,i2=i1+1;
      *bb++ = BoundaryEdge(v,i1,i2,3);
      v[i1].lab=v[i2].lab=3;

    }     
  for (int j=0;j<ny;j++)
    { // left
      int i=0;
      int i1= i + j*nx1,i2=i1+nx1;
      *bb++ = BoundaryEdge(v,i1,i2,4);
      v[i1].lab=v[i2].lab=4;
    }     
    
 {
  Mesh * m = new Mesh(nbv,nbt,neb,v,t,b);
  R2 Pn,Px;
  m->BoundingBox(Pn,Px);
  m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
/*  initgraphique();
  m->Draw();
  rattente(1);*/
  m->decrement();
  return m;
  }

}

class MeshCarre2 :   public E_F0mps { public:
    typedef pmesh  Result;
     Expression nx,ny;
     MeshCarre2(const basicAC_F0 & args) 
    { args.SetNameParam();
      nx=to<long>(args[0]); 
      ny=to<long>(args[1]); 
     }
    
    static ArrayOfaType  typeargs() { 
        return  ArrayOfaType(atype<long>(),atype<long>(),false);}
        
    static  E_F0 * f(const basicAC_F0 & args){
        return new MeshCarre2(args);} 
        
    AnyType operator()(Stack s) const { 
      return SetAny<pmesh>(Carre( GetAny<long>( (*nx)(s)) , GetAny<long>( (*ny)(s)),0,0,s ));}
    operator aType () const { return atype<pmesh>();} 
      
};
class MeshCarre2f :   public E_F0mps { public:
    typedef pmesh  Result;
     Expression nx,ny;
     Expression fx,fy;
     MeshCarre2f(const basicAC_F0 & args) 
    { args.SetNameParam();
      nx=to<long>(args[0]); 
      ny=to<long>(args[1]); 
      const E_Array *  a= dynamic_cast<const E_Array*>(args[2].LeftValue());
      throwassert(a);fx=0;fy=0;
      if (a->size()>0) fx=to<double>( (*a)[0]);
      if (a->size()>1) fy=to<double>( (*a)[1]);
     }
    
    static ArrayOfaType  typeargs() { 
        return  ArrayOfaType(atype<long>(),atype<long>(),atype<E_Array>(),false);}
        
    static  E_F0 * f(const basicAC_F0 & args){
        return new MeshCarre2f(args);} 
        
    AnyType operator()(Stack s) const { 
      return SetAny<pmesh>(Carre( GetAny<long>( (*nx)(s)) , GetAny<long>( (*ny)(s)), fx,fy,s ));}

    operator aType () const { return atype<pmesh>();} 
      
};


/*

Grid * Etruncmesh::eval()
{
  Analvar save(*an);
  throwassert(idgrid && idgrid->typesol ==Iden::maillage );
  Grid* go = idgrid->fg;
  throwassert(go);
  int * flag= new int[go->nt];
  int * bb  = new int[go->nt];
  Real xl[]={ 1./3.,1./3.,1./3.};      
  for (int i=0;i<go->nt;i++)
    {
      int oldlocal = an->local;
      const bTriangle & T = go->t[i];
      const bVertex & v0 = *T.v[0];
      const bVertex & v1 = *T.v[1];
      const bVertex & v2 = *T.v[2];
      Real x = v0.x*xl[0] + v1.x*xl[1] + v2.x*xl[2];
      Real y = v0.y*xl[0] + v1.y*xl[1] + v2.y*xl[2];                        
      an->setAn(0,x, y, T.where, xl,-1,-1,i);
      Real ee = e->eval();
      flag[i] = (int) Max((Real)-32000.0,Min(e->eval(),(Real)32000.0));
      if (b) 
        bb[i]  = (int) Max((Real)-32000.0,Min(b->eval(),(Real)32000.0));
      else 
        bb[i] = 1;
      //   cout << ee  << " " << flag[i] <<  " " << bb[i] << endl;

      an->local  = oldlocal;
    }
  Grid* g = new Grid();
  // for (int i=0;i<go->nt;i++)
  //  cout << flag[i] << (i%10 == 9 ? '\n' : ' ');
  cout << endl;
  Triangles * Th = new Triangles(*go->Th,flag,bb);
  delete [] flag;
  delete [] bb;
  if( ! Th) erreur("trunc triangulation");
  double hmax = Th->MaximalHmax();
  //  cout << " hmax = " << hmax << " ------- " << endl;
  Metric M(hmax);
  for (int iv=0;iv < Th->nbv;iv++)
    (*Th)[iv].m = M;
  
#ifdef DRAWING1
  reffecran();
  Th->InitDraw();
  Th->inquire();
#endif 
  
  g->th2t(Th);
  g->renum();
  g->prepgrid(0);
   
  if(!toScilab)
  g->draw(*an);
  // an->activeMesh=g; // set the activegrid
  *an=save;
  return g; 
}

Grid * Emovemesh::eval()
{
   Analvar save(*an); 
  int i;
  Grid* go = idmoved->fg;
  Grid* gn = new Grid(go);
  throwassert(go);
  
  Grid& t =*gn;
  an->gridxyng = go;
  Real xl[3] = {0.,0.,0.};
  for ( i = 0; i < go->nv; i++)
    {
      int oldlocal = an->local;
      an->setAn(0, go->v[i].x, go->v[i].y, go->v[i].where,xl, i); 
      t.v[i].x = ex->eval();
      t.v[i].y = ey->eval();
      an->local  = oldlocal;
    }
  an->gridxyng =0;
  t.prepgrid(1);
  //  for(i=0;i<t.nt;i++)
  //  for(int iloc=0;iloc<3;iloc++)
  //    t.t[i].e[iloc] = go->t[i].e[iloc];
  if(!toScilab)
  t.draw(*an);
  Geometry * tGh = new Geometry(go->Th->Gh);
  Triangles* tTh = new Triangles(*go->Th,tGh);// copy the Triangles
  cout << "\t\t MoveMesh Grid * " << gn << " Gh = " <<  tGh << " Th = " << tTh  << endl;
  //tTh->Write("movemesh.Th");
  
  Triangles & Th = *tTh;
  Geometry & Gh = *tGh;
  Geometry & GhO = go->Th->Gh;
  int * renu = go->NumThinGrid;
  
  cout << "\t\t renu = " << renu << " " << (renu ? renu[0] : 0) << endl;
  // move of the geometry 
  for ( i = 0; i <GhO.nbv;i++)
    {
      int oldlocal = an->local;
      an->setAn(0,GhO.vertices[i].r.x, GhO.vertices[i].r.y, GhO.vertices[i].ref(),xl);
      Gh.vertices[i].r.x  = ex->eval();
      Gh.vertices[i].r.y  = ey->eval();
      an->local  = oldlocal;
    }
  // change the tangente 

  for (i=0;i<Gh.nbe;i++) 
    {
      R2 AB = Gh.edges[i].v[1]->r - Gh.edges[i].v[0]->r;        
      Real8 lAB = Norme2(AB); // length of current edge AB
        
      for (int jj=0;jj<2;jj++) 
        if( ! Gh.edges[i].v[jj]->Corner() &&
            (    (jj==0 && Gh.edges[i].TgA()) 
                 || (jj==1 && Gh.edges[i].TgB()) ) )
          {       
                
            // recompute the tangent
            R2 tg =  Gh.edges[i].v[1-jj]->r 
              - Gh.edges[i].Adj[jj]->v[1-Gh.edges[i].SensAdj[jj]]->r;
            Real8 ltg =  Norme2(tg);
            tg =  tg *(lAB/ltg);
            if ( (tg,AB) < 0) 
              tg = -tg;
            Gh.edges[i].tg[jj] = tg;                 
          }          
    } // for (i=0;i<nbe;i++)
if(!toScilab)
{    
#ifdef DRAWING
  Gh.InitDraw();
#endif
#ifdef DRAWING2    
  reffecran();
  Gh.InitDraw();
  Gh.Draw();
#endif 
}     
  //  move of the bamg mesh 
  Th.pmin =  R2(t.v[0].x,t.v[0].y);
  Th.pmax =  Th.pmin;

  for ( i = 0; i < Th.nbv; i++)
    {   // Be carefull we do a renumbering 
      int j = renu ? renu[i] : i;
      Th.vertices[i].r =   R2(t.v[j].x,t.v[j].y);
      Th.pmin.x = Min( Th.pmin.x, Th.vertices[i].r.x);
      Th.pmin.y = Min( Th.pmin.y, Th.vertices[i].r.y);
      Th.pmax.x = Max( Th.pmax.x, Th.vertices[i].r.x);
      Th.pmax.y = Max( Th.pmax.y, Th.vertices[i].r.y);          
    }
  {R2 P10 = (Th.pmax-Th.pmin)*0.1;
  Th.pmin = Th.pmin - P10;Th.pmax = Th.pmax + P10;}

  Gh.pmin =  Gh.vertices[0].r;
  Gh.pmax =  Gh.pmin;

  for ( i = 0; i < Gh.nbv; i++)
    {
      Gh.pmin.x = Min( Gh.pmin.x, Gh.vertices[i].r.x);
      Gh.pmin.y = Min( Gh.pmin.y, Gh.vertices[i].r.y);
      Gh.pmax.x = Max( Gh.pmax.x, Gh.vertices[i].r.x);
      Gh.pmax.y = Max( Gh.pmax.y, Gh.vertices[i].r.y);          
    }
  {R2 P10 = (Gh.pmax-Gh.pmin)*0.1;
  Gh.pmin = Gh.pmin - P10;Gh.pmax = Gh.pmax + P10;}
    
  delete [] renu;
  renu=0;
  
  Gh.coefIcoor = (MaxICoor)/(Max( Gh.pmax.x- Gh.pmin.x, Gh.pmax.y- Gh.pmin.y));
  Th.coefIcoor = (MaxICoor)/(Max( Th.pmax.x- Th.pmin.x, Th.pmax.y- Th.pmin.y));

  //  for ( i = 0; i < Gh.nbv; i++)
  //    Gh.vertices[i].i = Gh.toI2(Gh.vertices[i].r);
  
  for ( i = 0; i < Th.nbv; i++)
    Th.vertices[i].i = Th.toI2(Th.vertices[i].r);
  // remove all the adj  and save the flag
  for ( i= 0; i < Th.nbt; i++)
    { 
      Triangle & t = Th(i);
      for ( int j = 0 ; j<3; j++) 
        t.SetAdj2(j,0,t.GetAllflag(j));
    }
   

  //  Th.quadtree=new QuadTree(&Th);
  Th.nbt = Th.nbt - Th.NbOutT; // remove all the  the ouside triangles 
  Th.SetIntCoor("In movemesh"); 
  Th.FillHoleInMesh(); 
  // delete Th.quadtree; // delete old 
  if (!Th.quadtree)
    Th.quadtree=new QuadTree(&Th);
  Th.ReMakeTriangleContainingTheVertex();
  

#ifdef DRAWING2 
  Th.inquire();
#endif  
  gn->Th = &Th;
  gn->Gh = &Gh;
  Gh.NbRef++;
  gn->nbholes = go->nbholes;
  *an=save;  
  return gn;
}

*/   

   void  MeshErrorIO(ios& )
{
   ExecError("Mesh IO Error ");
}

inline pmesh *  initMesh(pmesh * const & p, string * const & s) {
  Mesh * m;
  *p= m =new Mesh(*s); 
  m->MakeQuadTree();
  delete s;  
  return p;
 }

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
      
      throwassert(a);
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
   throwassert(Thh);
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
    return SetAny<double>(minarea/2.);

}


void init_lgmesh() {
    cout <<"lg_mesh ";
    bamg::MeshIstreamErrorHandler = MeshErrorIO;
   Global.Add("adaptmesh","(",new OneOperatorCode<Adaptation>);
   Global.Add("movemesh","(",new OneOperatorCode<MoveMesh>);
   Global.Add("splitmesh","(",new OneOperatorCode<SplitMesh>);
   Global.Add("checkmovemesh","(",new OneOperatorCode<CheckMoveMesh>);
   Global.Add("square","(",new OneOperatorCode<MeshCarre2>);
   Global.Add("square","(",new OneOperatorCode<MeshCarre2f>);
   Global.Add("savemesh","(",new OneOperatorCode<SaveMesh>);
   Global.Add("trunc","(", new Op_trunc_mesh);
   Global.Add("readmesh","(",new OneOperator1_<pmesh,string*>(ReadMeshbamg));
   Global.Add("triangulate","(",new OneOperator1_<pmesh,string*>(ReadTriangulate));
   TheOperators->Add("<-",
       new OneOperator2_<pmesh*,pmesh*,string* >(&initMesh));
       
// use for :   mesh Th = readmesh ( ...);       
  TheOperators->Add("<-",
       new OneOperator2_<pmesh*,pmesh*,pmesh >(&set_copy_incr));
   
}
#include "InitFunct.hpp"
static addingInitFunct TheaddingInitFunct(-10,init_lgmesh);
