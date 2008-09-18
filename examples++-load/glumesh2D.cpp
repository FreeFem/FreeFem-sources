// $Id$

#include <iostream>
#include <cfloat>
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
  listMesh(Stack s,const listMesh &l,Mesh *th) : lth(Add2StackOfPtr2Free(s,new list<Mesh*>(*l.lth))) { lth->push_back(th);}

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
  {
    map<pair<int,int>,int> bbe;
    for(list<Mesh *>::const_iterator i=lth.begin();i != lth.end();++i)
      {
	const Mesh &Th(**i);
	if(!*i) continue;
	if(verbosity>1)  cout << " GluMesh + "<< Th.nv << " " << Th.nt << endl;
	int nbv0 = nbv;
	
	for (int ii=0;ii<Th.nv;ii++)
	{
	  const Vertex &vi(Th(ii));
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
	
	// bug glumesh done after 
      
	for (int k=0;k<Th.neb;k++)
	  {
	    const BoundaryEdge & be(Th.bedges[k]);
	    int i0=quadtree->NearestVertex(be[0])-v;
	    int i1=quadtree->NearestVertex(be[1])-v;
	    int ii0=i0,ii1=i1;
	    if(ii1<ii0) Exchange(ii0,ii1);
	    pair<int,int> i01(ii0,ii1);
	    if( bbe.find(i01) == bbe.end())
	    {	    
	      (*bb++).set(v,i0,i1,be.lab);
	      bbe[i01]=	  neb++;
	    }
	    
	  }
	
      } 
  }
  /*
  Mesh::Vertex *becog = new Vertex[nebx];
  cout << "creation quadtree" << endl;
  FQuadTree *quadtree_be=new Fem2D::FQuadTree(becog,Pn,Px,-nebx);
  cout << "fin creation quadtree" << endl;
  double hseuil_border = hseuil/2.;
  // to remove common egde wrong code .... FH
  if(0) 
    for(list<Mesh *>::const_iterator i=lth.begin();i != lth.end();++i)
      {
      const Mesh &Th(**i);
      if(!*i) continue;
      
      if(verbosity>1)  cout << " GluMesh + "<< Th.nv << " " << Th.nt << endl;
      for (int k=0;k<Th.neb;k++)
	{
	  
	  const BoundaryEdge & be(Th.bedges[k]);
	  int i0 = Th.operator()(be[0]); //quadtree->NearestVertex(be[0])-v;
	  int i1 = Th.operator()(be[1]); //quadtree->NearestVertex(be[1])-v;
	  
	  R2 r2vi= ((R2) be[0] + be[1])/2.;
	  
	  //const R2 r2vi( cdgx, cdgy);
	  const Mesh::Vertex & vi(r2vi);
	  
	  Vertex * pvi=quadtree_be->ToClose(vi,hseuil_border);
	  if(!pvi){
	    becog[neb].x = vi.x;
	    becog[neb].y = vi.y;
	    becog[neb].lab = vi.lab;
	    quadtree_be->Add( becog[neb++] );
	    
	    int iglu0=quadtree->NearestVertex(be[0])-v;
	    int iglu1=quadtree->NearestVertex(be[1])-v;
	    
	    (bb++)->set(v,iglu0,iglu1,vi.lab);
	  }
	}
    
    } //  
  */
  delete quadtree;
  //delete quadtree_be;
   
  if(verbosity>1)
    {
      cout << "     Nb points : "<< nbv << " , nb edges : " << neb << endl;
      cout << "     Nb of glu point " << nbvx -nbv;
      cout << "     Nb of glu  Boundary edge " << nebx-neb;
    }
  
  {
    Mesh * m = new Mesh(nbv,nbt,neb,v,t,b);
    R2 Pn,Px;
    m->BoundingBox(Pn,Px);
    m->quadtree=new Fem2D::FQuadTree(m,Pn,Px,m->nv);
    //    m->decrement();
    return m;
  }
  
}

template<class RR,class AA=RR,class BB=AA> 
struct Op2_addmesh: public binary_function<AA,BB,RR> { 
  static RR f(Stack s,const AA & a,const BB & b)  
  { return RR(s, a, b );} 
};

template<bool INIT,class RR,class AA=RR,class BB=AA> 
struct Op2_setmesh: public binary_function<AA,BB,RR> { 
  static RR f(Stack stack, const AA & a,const BB & b)  
  {    
    ffassert(a );
    pmesh  p=GluMesh(b);
    
    if(!INIT &&  *a) delete *a;
    //  Add2StackOfPtr2FreeRC(stack,p); //  the pointer is use to set variable so no remove. 
    return *a=p,a;
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
  //  m->decrement();
  Add2StackOfPtr2FreeRC(stack,m);
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
  
  //Dcl_Type<listMesh3>();
  //typedef Mesh3 *pmesh3;
  
  if (verbosity)
    cout << " load: glumesh2D  " << endl;
  //cout << " je suis dans Init " << endl; 
  TheOperators->Add("+",new OneBinaryOperator_st< Op2_addmesh<listMesh,pmesh,pmesh>  >      );
  TheOperators->Add("+",new OneBinaryOperator_st< Op2_addmesh<listMesh,listMesh,pmesh>  >      );
  TheOperators->Add("=",new OneBinaryOperator_st< Op2_setmesh<false,pmesh*,pmesh*,listMesh>  >     );
  TheOperators->Add("<-",new OneBinaryOperator_st< Op2_setmesh<true,pmesh*,pmesh*,listMesh>  >     );
  
  Global.Add("change","(",new SetMesh);

}
