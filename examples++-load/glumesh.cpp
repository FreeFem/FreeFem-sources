// $Id$

#include  <iostream>
#include  <cfloat>
#include <cmath>
#include <complex>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;  
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"


#include "FESpacen.hpp" 
#include "FESpace.hpp" 

#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "Operator.hpp" 
#include "lex.hpp"

#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
#include "LayerMesh.hpp"

#include <set>
#include <vector>
#include <fstream>


  using namespace  Fem2D;

class listMesh { 
public:
  list<Mesh *> *lth;
  void init()  { lth=new list<Mesh *>;}
  void destroy() { delete lth;}
  listMesh(Stack s,Mesh *th) : lth(Add2StackOfPtr2Free(s,new list<Mesh*>)) { lth->push_back(th);}
  listMesh(Stack s,Mesh *tha,Mesh *thb) : lth(Add2StackOfPtr2Free(s,new list<Mesh*>)) { lth->push_back(tha);lth->push_back(thb);}
  listMesh(Stack s,listMesh &l,Mesh *th) : lth(Add2StackOfPtr2Free(s,new list<Mesh*>(*l.lth))) { lth->push_back(th);}

};

Mesh * GluMesh(listMesh const & lst)
{
  int nbt=0;
  int neb=0;
  int nebx=0;
  int nbv=0;
  int nbvx=0;
  double hmin=1e100;
  R2 Pn(1e100,1e100),Px(-1e100,-1e100);
  const list<Mesh *> lth(*lst.lth);
   Mesh * th0=0;
  for(list<Mesh *>::const_iterator i=lth.begin();i != lth.end();++i)
    {
       Mesh &Th(**i);
      th0=&Th;
      if(verbosity>1)  cout << " GluMesh + "<< Th.nv << " " << Th.nt << endl;
      nbt+= Th.nt;
      nbvx += Th.nv;
      nebx += Th.neb;
      for (int k=0;k<Th.nt;k++)
	for (int e=0;e<3;e++)
	  hmin=min(hmin,Th[k].lenEdge(e));
      
      for (int i=0;i<Th.nv;i++)
	{ 	  
	  R2 P(Th(i));
	  Pn=Minc(P,Pn);
	  Px=Maxc(P,Px);     
	}
    } 

  cout << "      - hmin =" <<  hmin << " ,  Bounding Box: " << Pn << " "<< Px << endl;
  
  Vertex * v= new Vertex[nbvx];
  Triangle *t= new Triangle[nbt];
  Triangle *tt=t;
  BoundaryEdge *b= new BoundaryEdge[nebx];
  BoundaryEdge *bb= b;
  
  ffassert(hmin>Norme2(Pn-Px)/1e9);
  double hseuil =hmin/10.; 
  
  
  FQuadTree *quadtree=new Fem2D::FQuadTree(th0,Pn,Px,0);
  
  for(list<Mesh *>::const_iterator i=lth.begin();i != lth.end();++i)
    {
      const Mesh &Th(**i);
      if(verbosity>1)  cout << " GluMesh + "<< Th.nv << " " << Th.nt << endl;
      int nbv0 = nbv;
      for (int i=0;i<Th.nv;i++)
	{
	  const Vertex &vi(Th(i));
	  Vertex * pvi=quadtree->ToClose(vi,hseuil);
	  if(!pvi) { 
	    v[nbv].x = vi.x;
	    v[nbv].y = vi.y;
	    v[nbv].lab = vi.lab;
	    quadtree->Add(v[nbv++]);
	  }
	}
      
      for (int k=0;k<Th.nt;k++)
	{
	  const Triangle  &K(Th[k]);
	  int i0=quadtree->NearestVertex(K[0])-v;
	  int i1=quadtree->NearestVertex(K[1])-v;
	  int i2=quadtree->NearestVertex(K[2])-v;	  
	  (*tt++).set(v,i0,i1,i2,K.lab);
	}
      
      for (int k=0;k<Th.neb;k++)
	{
	  const BoundaryEdge & be(Th.bedges[k]);
	  int i0=quadtree->NearestVertex(be[0])-v;
	  int i1=quadtree->NearestVertex(be[1])-v;
	  if(i1<nbv0 && i0 < nbv) continue;
	  (*bb++).set(v,i0,i1,be.lab);
	  neb++;
	}
      
    } //  
  delete quadtree;
 
  if(verbosity>1)
    {
      cout << "     Nb of glu point " << nbvx -nbv;
      cout << "     Nb of glu  Boundary edge " << nebx-neb;
    }
  
  {
    Mesh * m = new Mesh(nbv,nbt,neb,v,t,b);
    R2 Pn,Px;
    m->BoundingBox(Pn,Px);
    m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
    m->decrement();
    return m;
  }
  
}

template<class RR,class AA=RR,class BB=AA> 
struct Op2_addmesh: public binary_function<AA,BB,RR> { 
  static RR f(Stack s,const AA & a,const BB & b)  
  { return RR(s, a, b );} 
};

template<class RR,class AA=RR,class BB=AA> 
struct Op2_setmesh: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  {
    ffassert(a );
    if(*a) delete *a;
    return (*a=GluMesh(b),a); 
  } 
};

class BuildLayeMesh_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression enmax,ezmin,ezmax,xx,yy,zz;
  static const int n_name_param =5; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  
public:
  BuildLayeMesh_Op(const basicAC_F0 &  args,Expression tth,Expression nmaxx) 
    : eTh(tth),enmax(nmaxx) 
  {
    args.SetNameParam(n_name_param,name_param,nargs);
    const E_Array * a2 =0, *a1=0 ;
    if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[0]);
    if(nargs[1])  a2  = dynamic_cast<const E_Array *>(nargs[1]);
    int err =0;
    if(a1) {
      if(a1->size() !=2) 
	CompileError("LayerMesh (Th,n, zbound=[zmin,zmax],) ");
      ezmin=to<double>( (*a1)[0]);
      ezmax=to<double>( (*a2)[1]);
    }
    if(a2) {
      if(a2->size() !=3) 
	CompileError("LayerMesh (Th,n, transfo=[X,Y,Z],) ");
      xx=to<double>( (*a2)[0]);
      yy=to<double>( (*a2)[1]);
      zz=to<double>( (*a2)[2]);
    }    
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type BuildLayeMesh_Op::name_param[]= {
  {  "zbound", &typeid(E_Array)},
  {  "transfo", &typeid(E_Array)},
  {  "coef", &typeid(double)},
  {  "reftet", &typeid(KN_<long> )},
  {  "refface", &typeid(KN_<long> )}
};


class  BuildLayerMesh : public OneOperator { public:  
    BuildLayerMesh() : OneOperator(atype<pmesh3>(),atype<pmesh>(),atype<long>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new BuildLayeMesh_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1])); 
  }
};

class SetMesh_Op : public E_F0mps 
{
public:
  Expression a; 
  
  static const int n_name_param =2; //  add nbiter FH 30/01/2007 11 -> 12 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}

  
public:
  SetMesh_Op(const basicAC_F0 &  args,Expression aa) : a(aa) {
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type SetMesh_Op::name_param[]= {
  {  "refe", &typeid(KN_<long> )},
  {  "reft", &typeid(KN_<long> )}
};

int  ChangeLab(const map<int,int> & m,int lab)
{
  map<int,int>::const_iterator i=m.find(lab);
  if(i != m.end())
    lab=i->second;
  return lab;
}

AnyType BuildLayeMesh_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  int nlayer = (int) GetAny<long>((*enmax)(stack));
  ffassert(pTh && nlayer>0);
  Mesh & Th=*pTh;
  Mesh *m= pTh;
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int neb=Th.neb; // nombre d'aretes fontiere
  cout << " " << nbv<< " "<< nbv << " nbe "<< neb << endl; 
  KN<double> zmin(nbv),zmax(nbv);
  KN<double> clayer(nbv);
  //  nombre de layer est nlayer*clayer
  clayer=-1;
  zmin=0.;
  zmax=1.;
  for (int it=0;it<nbt;++it)
    for(int iv=0;iv<3;++iv)      
    {
      int i=Th(it,iv);
      if(clayer[i]<0)
	{
	  mp->setP(&Th,it,iv);
	  if(ezmin)  zmin[i]=GetAny<double>((*ezmin)(stack)); 
	  if(ezmax)  zmax[i]=GetAny<double>((*ezmax)(stack)); 
	  clayer[i]=Max( 0. , Min( 1. , arg(2,stack,1.) ) ); 
	}
			     
    }
  ffassert(clayer.min() >=0);
  KN<long> zzempty;
  KN<long> nre (arg(3,stack,zzempty));  
  KN<long> nrt (arg(4,stack,zzempty));  

  if(nre.N() <=0 && nrt.N() ) return m;
  ffassert( nre.N() %2 ==0);
  ffassert( nrt.N() %2 ==0);
  map<int,int> mape,mapt;
  int z00 = false;
  for(int i=0;i<nre.N();i+=2)
    { z00 = z00 || ( nre[i]==0 && nre[i+1]==0);
      if(nre[i] != nre[i+1])
	mape[nre[i]]=nre[i+1];
    }
  for(int i=0;i<nrt.N();i+=2)
    mapt[nrt[i]]=nrt[i+1];
  int nebn =0;
  KN<int> ni(nbv);
  for(int i=0;i<nbv;++i)
    ni[i]=Max(0,Min(nlayer,(int) lrint(nlayer*clayer[i])));
  
  
  Mesh3 *Th3= build_layer(Th, nlayer,ni,zmin, zmax);

  Th3->BuildBound();
  Th3->BuildAdj();
  Th3->Buildbnormalv();  
  Th3->BuildjElementConteningVertex();
  Th3->BuildGTree();
  Th3->decrement();  
  *mp=mps;
  return Th3;
}

AnyType SetMesh_Op::operator()(Stack stack)  const 
{
  Mesh * pTh= GetAny<Mesh *>((*a)(stack));
  Mesh & Th=*pTh;
  Mesh *m= pTh;
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int neb=Th.neb; // nombre d'aretes fontiere
  cout << " " << nbv<< " "<< nbv << " nbe "<< neb << endl;  
  KN<long> zz;
  KN<long> nre (arg(0,stack,zz));  
  KN<long> nrt (arg(1,stack,zz));  

  if(nre.N() <=0 && nrt.N() ) return m;
  ffassert( nre.N() %2 ==0);
  ffassert( nrt.N() %2 ==0);
  map<int,int> mape,mapt;
  int z00 = false;
  for(int i=0;i<nre.N();i+=2)
    { z00 = z00 || ( nre[i]==0 && nre[i+1]==0);
      if(nre[i] != nre[i+1])
	mape[nre[i]]=nre[i+1];
    }
  for(int i=0;i<nrt.N();i+=2)
    mapt[nrt[i]]=nrt[i+1];
  int nebn =0;
    for(int i=0;i<neb;++i)
      {
	int l0,l1=ChangeLab(mape,l0=m->bedges[i].lab) ;
	 nebn++;
      }
	  

  Vertex * v= new Vertex[nbv];
  Triangle *t= new Triangle[nbt];
  BoundaryEdge *b= new BoundaryEdge[nebn];
  // generation des nouveaux sommets 
  Vertex *vv=v;
  // copie des anciens sommets (remarque il n'y a pas operateur de copy des sommets)
  for (int i=0;i<nbv;i++)
    {
     Vertex & V=Th(i);
     vv->x=V.x;
     vv->y=V.y;
     vv->lab = V.lab;
     vv++;      
   }

  //  generation des triangles 
  Triangle *tt= t; 
  int nberr=0;
   
  for (int i=0;i<nbt;i++)
    {
      int i0=Th(i,0), i1=Th(i,1),i2=Th(i,2);
      // les 3 triangles par triangles origines 
      (*tt++).set(v,i0,i1,i2,ChangeLab(mapt,Th[i].lab));
    }  
  
  // les arete frontieres qui n'ont pas change
  BoundaryEdge * bb=b;
  for (int i=0;i<neb;i++)
    {        
      int i1=Th(Th.bedges[i][0]);
      int i2=Th(Th.bedges[i][1]);
      int l0,l1=ChangeLab(mape,l0=m->bedges[i].lab) ;
      *bb++ = BoundaryEdge(v,i1,i2,l1);   
    }
  assert(nebn==bb-b);
  m =  new Mesh(nbv,nbt,nebn,v,t,b);

  R2 Pn,Px;                                                                                                                               
  m->BoundingBox(Pn,Px);
  m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);   
  m->decrement();
  return m;
}


class SetMesh : public OneOperator { public:  
typedef Mesh *pmesh;
    SetMesh() : OneOperator(atype<pmesh>(),atype<pmesh>() ) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new SetMesh_Op(args,t[0]->CastTo(args[0])); 
  }
};


//  truc pour que la fonction 
// Init::Init() soit appele a moment du chargement dynamique
// du fichier 
//  
class Init { public:
  Init();
};

static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  Dcl_Type<listMesh>();
  typedef Mesh *pmesh;
  if (verbosity)
    cout << " lood: glumesh  " << endl;
  TheOperators->Add("+",new OneBinaryOperator_st< Op2_addmesh<listMesh,pmesh,pmesh>  >      );
  TheOperators->Add("=",new OneBinaryOperator< Op2_setmesh<pmesh*,pmesh*,listMesh>  >     );
  TheOperators->Add("<-",new OneBinaryOperator< Op2_setmesh<pmesh*,pmesh*,listMesh>  >     );
  Global.Add("change","(",new SetMesh);
  Global.Add("buildlayers","(",new  BuildLayerMesh);
}
