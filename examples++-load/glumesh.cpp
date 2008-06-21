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
#include "TransfoMesh_v2.hpp"
//#include "GQuadTree.hpp"

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
      assert( *i );
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
      if(!*i) continue;
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
	  if(i1<nbv0 && i0 < nbv0) continue;
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
    //    m->decrement();
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

// glumesh3D

class listMesh3 { 
public:
  list<Mesh3 *> *lth;
  void init()  { lth=new list<Mesh3 *>;}
  void destroy() { delete lth;}
  listMesh3(Stack s,Mesh3 *th) : lth(Add2StackOfPtr2Free(s,new list<Mesh3*>)) { lth->push_back(th);}
  listMesh3(Stack s,Mesh3 *tha,Mesh3 *thb) : lth(Add2StackOfPtr2Free(s,new list<Mesh3*>)) { lth->push_back(tha);lth->push_back(thb);}
  listMesh3(Stack s,const listMesh3 &l,Mesh3 *th) : lth(Add2StackOfPtr2Free(s,new list<Mesh3*>(*l.lth))) { lth->push_back(th);}

};

Mesh3 * GluMesh3(listMesh3 const & lst)
{
  int nbt=0;
  int nbe=0;
  int nbex=0;
  int nbv=0;
  int nbvx=0;
  
  double hmin=1e100;
  R3 Pn(1e100,1e100,1e100),Px(-1e100,-1e100,-1e100);
  const list<Mesh3 *> lth(*lst.lth);
  Mesh3 * th0=0;
  
  for(list<Mesh3 *>::const_iterator i=lth.begin();i != lth.end();++i)
    {
	    Mesh3 &Th3(**i);  // definis ???
	    th0=&Th3;
      	if(verbosity>1)  cout << " GluMesh3 + "<< Th3.nv << " " << Th3.nt << " "<< Th3.nbe << endl;
      	
      	nbt+= Th3.nt;
      	nbvx += Th3.nv;
      	nbex += Th3.nbe;
      	
      	for (int k=0;k<Th3.nt;k++){
			for (int e=0;e<6;e++){
				hmin=min(hmin,Th3[k].lenEdge(e));   // calcul de .lenEdge pour un Mesh3
			}
		}
		for (int i=0;i<Th3.nv;i++){ 	  
			R3 P(Th3(i));
			Pn=Minc(P,Pn);
			Px=Maxc(P,Px);     
		}
    } 
  cout << "      - hmin =" <<  hmin << " ,  Bounding Box: " << Pn << " "<< Px << endl;

  
  Vertex3  *v= new Vertex3[nbvx];
  Tet      *t= new Tet[nbt];
  Tet      *tt=t;
  Triangle3 *b= new Triangle3[nbex];
  Triangle3 *bb= b;
     
  ffassert(hmin>Norme2(Pn-Px)/1e9);
  double hseuil =hmin/10.; 
  
  cout << " creation of : BuildGTree" << endl;
  //th0->BuildBound();
  //th0->BuildGTree();   
 
  EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,Pn,Px,0);
    
  cout << " creation of : Mesh3 th3_tmp" << endl;
  /*Mesh3 th3_tmp;
  th3_tmp.set(nbvx,nbt,nbex);*/
  
  
  for(list<Mesh3 *>::const_iterator i=lth.begin();i != lth.end();++i){
      const Mesh3 &Th3(**i);
      if(!*i) continue;
      if(verbosity>1)  cout << " loop over mesh for create new mesh "<< endl;
      if(verbosity>1)  cout << " GluMesh + "<< Th3.nv << " " << Th3.nt <<" " << Th3.nbe << endl;
      int nbv0 = nbv;
      
      for (int i=0;i<Th3.nv;i++){
	      const Vertex3 &vi(Th3.vertices[i]);
	      Vertex3 * pvi=gtree->ToClose(vi,hseuil);
	      if(!pvi){
		      v[nbv].x = vi.x;
		      v[nbv].y = vi.y;
		      v[nbv].z = vi.z;
		      v[nbv].lab = vi.lab;
		      gtree->Add( v[nbv++] );
		  }
	  }
      
      for (int k=0;k<Th3.nt;k++){
	  	const Tet  &K(Th3[k]);
	  	int iv[4];
	  	iv[0]=gtree->NearestVertex(K[0])-v;
	  	iv[1]=gtree->NearestVertex(K[1])-v;
	  	iv[2]=gtree->NearestVertex(K[2])-v;  //-v;
	  	iv[3]=gtree->NearestVertex(K[3])-v;  //-v;
	  	(tt++)->set(v,iv,K.lab);
	  }
      
      for (int k=0;k<Th3.nbe;k++)
	{
	  const Triangle3 & K(Th3.be(k));//bedges[k]);
	  int iv[3];
	  iv[0]=gtree->NearestVertex(K[0])-v; //-v;
	  iv[1]=gtree->NearestVertex(K[1])-v; //-v;
	  iv[2]=gtree->NearestVertex(K[2])-v; //-v;
	  
	  if( iv[2]<nbv0 && iv[1]<nbv0 && iv[0] < nbv0 ) continue;
	  (bb++)->set(v,iv,K.lab);
	  nbe++;
	}
      
  }   
  //delete th0;
 
  if(verbosity>1)
    {
      cout << "     Nb of glu3D  point " << nbvx-nbv;
      cout << "     Nb of glu3D  Boundary edge " << nbex-nbe << endl;
    }
	
    Mesh3 *mpq= new Mesh3(nbv,nbt,nbe,v,t,b);  
	/*       
	mpq->BuildBound();
	mpq->BuildAdj();
	mpq->Buildbnormalv();  
	mpq->BuildjElementConteningVertex();
	*/
	mpq->BuildGTree();
	mpq->decrement();  
    return mpq;
  
}

template<class RR,class AA=RR,class BB=AA> 
struct Op3_addmesh: public binary_function<AA,BB,RR> { 
  static RR f(Stack s,const AA & a,const BB & b)  
  { return RR(s, a, b );} 
};

template<class RR,class AA=RR,class BB=AA> 
struct Op3_setmesh: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  {
    ffassert(a );
    if(*a) delete *a;
    return (*a=GluMesh3(b),a); 
  } 
};


class BuildLayeMesh_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression enmax,ezmin,ezmax,xx,yy,zz;
  static const int n_name_param =7; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  
public:
  BuildLayeMesh_Op(const basicAC_F0 &  args,Expression tth,Expression nmaxx) 
    : eTh(tth),enmax(nmaxx), ezmin(0),ezmax(0),xx(0),yy(0),zz(0)
  {
    cout << "construction par BuilLayeMesh_Op" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
    const E_Array * a2 =0, *a1=0 ;
    if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[0]);
    if(nargs[1])  a2  = dynamic_cast<const E_Array *>(nargs[1]);
    int err =0;
    //cout << nargs[0] << " "<< a1 << endl;
    //cout << nargs[1] << " "<< a2 << endl;
    if(a1) {
      if(a1->size() !=2) 
	CompileError("LayerMesh (Th,n, zbound=[zmin,zmax],) ");
	  //cout << "lecture de ezmin , ezmax" << endl; 
      ezmin=to<double>( (*a1)[0]);
      ezmax=to<double>( (*a1)[1]); 
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
  {  "reffacemid", &typeid(KN_<long> )},
  {  "reffaceup", &typeid(KN_<long> )},
  {  "reffacedown", &typeid(KN_<long> )}
};


class  BuildLayerMesh : public OneOperator { public:  
    BuildLayerMesh() : OneOperator(atype<pmesh3>(),atype<pmesh>(),atype<long>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	cout << " je suis dans code(const basicAC_F0 & args) const" << endl; 
	//cout << "args: " << args << endl;   
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


//// version 3D de change label

class SetMesh3D_Op : public E_F0mps 
{
public:
  Expression a; 
  
  static const int n_name_param =2; //  add nbiter FH 30/01/2007 11 -> 12 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}

  
public:
  SetMesh3D_Op(const basicAC_F0 &  args,Expression aa) : a(aa) {
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type SetMesh3D_Op::name_param[]= {
  {  "reftet", &typeid(KN_<long> )},
  {  "refface", &typeid(KN_<long> )}
};
/*  besoin en cas de fichier 2D / fichier 3D 
int  ChangeLab(const map<int,int> & m,int lab)
{
  map<int,int>::const_iterator i=m.find(lab);
  if(i != m.end())
    lab=i->second;
  return lab;
}
*/
AnyType SetMesh3D_Op::operator()(Stack stack)  const 
{
  Mesh3 * pTh= GetAny<Mesh3 *>((*a)(stack));
  Mesh3 & Th=*pTh;
  Mesh3 *m= pTh;
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.nbe; // nombre d'aretes fontiere
  cout << "Number of Vertex="<< nbv << "Number of BorderElement=" << nbe << endl;  
  KN<long> zz;
  KN<long> nrtet (arg(0,stack,zz));  
  KN<long> nrface (arg(1,stack,zz));  

  if(nrface.N() <=0 && nrtet.N() ) return m;
  ffassert( nrtet.N() %2 ==0);
  ffassert( nrface.N() %2 ==0);
  
  map<int,int> maptet,mapface;
  
  int z00 = false;
  for(int i=0;i<nrface.N();i+=2)
    { z00 = z00 || ( nrface[i]==0 && nrface[i+1]==0);

      if( nrface[i] != nrface[i+1] ){	
      mapface[nrface[i]] = nrface[i+1];
      }      
    }
      
  for(int i=0;i<nrtet.N();i+=2)
    maptet[nrtet[i]]=nrtet[i+1];
  
  // sert a quoi ???  
  int nben =0;
  for(int i=0;i<nbe;++i)
  {
  	const Triangle3 & K = Th.be(i);	
    int l0,l1=ChangeLab(mapface, l0= K.lab) ;
	nben++;
  }
	  
	  
  Vertex3   *v = new Vertex3[nbv];
  Tet       *t = new Tet[nbt];
  Triangle3 *b = new Triangle3[nben];
  // generation des nouveaux sommets 
  Vertex3 *vv=v;
  // copie des anciens sommets (remarque il n'y a pas operateur de copy des sommets)
  for (int i=0;i<nbv;i++)
    {
     const Vertex3 & V(Th.vertices[i]);
     vv->x=V.x;
     vv->y=V.y;
     vv->z=V.z;
     vv->lab = V.lab;
     vv++;      
   }

  //  generation des triangles 
  Tet *tt= t; 
  int nberr=0;
   
  for (int i=0;i<nbt;i++)
    {
      const Tet &K( Th.elements[i] );	
      int iv[4];
      //int i0=Th(i,0), i1=Th(i,1),i2=Th(i,2);
      iv[0]= Th.operator()(K[0]);
      iv[1]= Th.operator()(K[1]);
      iv[2]= Th.operator()(K[2]);
      iv[3]= Th.operator()(K[3]);
      // les 3 triangles par triangles origines 
      int lab=K.lab;
      if(lab==41) cout << "i triangle "<< i << "lab= " << lab << "changelab " <<ChangeLab(maptet,lab) << endl;
      (*tt++).set( v, iv, ChangeLab(maptet,lab));
    }  
  
  // les arete frontieres qui n'ont pas change
  
  Triangle3 * bb=b;
  for (int i=0;i<nbe;i++)
    { 
      const Triangle3 &K( Th.be(i) );
      int iv[3];       
    
      iv[0] = Th.operator()(K[0]);
      iv[1] = Th.operator()(K[1]);
      iv[2] = Th.operator()(K[2]);
      
      int l0,l1=ChangeLab(mapface,l0=K.lab) ;
      (*bb++).set( v, iv, l1);   
    }
    
  assert(nben==bb-b);
  
  Mesh3 *mpq = new Mesh3(nbv,nbt,nbe,v,t,b);
  mpq->BuildGTree();
  mpq->decrement();  

  return mpq;
}


class SetMesh3D : public OneOperator { public:  
typedef Mesh *pmesh;
    SetMesh3D() : OneOperator(atype<pmesh3>(),atype<pmesh3>() ) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new SetMesh3D_Op(args,t[0]->CastTo(args[0])); 
  }
};


AnyType BuildLayeMesh_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  int nlayer = (int) GetAny<long>((*enmax)(stack));
  ffassert(pTh && nlayer>0);
  Mesh &Th=*pTh;
  Mesh *m= pTh;   // question a quoi sert *m ??
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int neb=Th.neb; // nombre d'aretes fontiere
  cout << " " << nbv<< " "<< nbv << " nbe "<< neb << endl; 
  KN<double> zmin(nbv),zmax(nbv);
  KN<double> clayer(nbv); //  nombre de layer est nlayer*clayer
  
  clayer=-1;
  zmin=0.;
  zmax=1.;
  for (int it=0;it<nbt;++it){
    for(int iv=0;iv<3;++iv)      
    {
      int i=Th(it,iv);
       if(clayer[i]<0)
	{
	  mp->setP(&Th,it,iv);
	  //cout << "mp: fait " << endl;
	  if(ezmin){ zmin[i]=GetAny<double>((*ezmin)(stack));}
	  if(ezmax){ zmax[i]=GetAny<double>((*ezmax)(stack));}

	  clayer[i]=Max( 0. , Min( 1. , arg(2,stack,1.) ) ); 
	
	}	     
    }
  }
  ffassert(clayer.min() >=0);

  cout << "lecture valeur des references " << endl;
  
  KN<long> zzempty;
  KN<long> nrtet  (arg(3,stack,zzempty));  
  KN<long> nrfmid (arg(4,stack,zzempty));  
  KN<long> nrfup  (arg(5,stack,zzempty));  
  KN<long> nrfdown (arg(6,stack,zzempty));  

  cout << nrtet.N() <<  nrfmid.N() << nrfup.N() << nrfdown.N() << endl;
  
  //if( nrtet.N() && nrfmid.N() && nrfup.N() && nrfdown.N() ) return m;
  ffassert( nrtet.N() %2 ==0);
  ffassert( nrfmid.N() %2 ==0);
  ffassert( nrfup.N() %2 ==0);
  ffassert( nrfdown.N() %2 ==0);
  
  // realisation de la map par default
  
  map< int, int > maptet; 
  map< int, int > maptrimil, maptrizmax, maptrizmin;
  map< int, int > mapemil, mapezmax, mapezmin;
  
  build_layer_map_tetrahedra( Th, maptet );
  build_layer_map_triangle( Th, maptrimil, maptrizmax, maptrizmin );
  build_layer_map_edge( Th, mapemil, mapezmax, mapezmin );
    
  // Map utilisateur
  map< int, int > :: iterator imap;
  for( int ii=0; ii < nrtet.N(); ii+=2){
	imap = maptet.find(nrtet[ii]);
	if( imap != maptet.end()){
		imap -> second = nrtet[ii+1];
	}
  }  

  for( int ii=0; ii < nrfmid.N(); ii+=2){
	imap = maptrimil.find(nrfmid[ii]);
	if( imap != maptrimil.end()){
		imap -> second = nrfmid[ii+1];
	}
  }  
  
  for( int ii=0; ii < nrfup.N(); ii+=2){
	imap = maptrizmax.find(nrfup[ii]);
	if( imap != maptrizmax.end()){
		imap -> second = nrfup[ii+1];
	}
  }  
  
  for( int ii=0; ii < nrfdown.N(); ii+=2){
	imap = maptrizmin.find(nrfdown[ii]);
	if( imap != maptrizmin.end()){
		imap -> second = nrfdown[ii+1];
	}
  }  
    
  int nebn =0;
  KN<int> ni(nbv);
  for(int i=0;i<nbv;i++)
    ni[i]=Max(0,Min(nlayer,(int) lrint(nlayer*clayer[i])));
    
//    map< int, int > maptet; 
// 	  map< int, int > maptrimil, maptrizmax, maptrizmin;
// 	  map< int, int > mapemil, mapezmax, mapezmin;
// 	   
// 	  build_layer_map_tetrahedra( Th, maptet );
// 	  build_layer_map_triangle( Th, maptrimil, maptrizmax, maptrizmin );
// 	  build_layer_map_edge( Th, mapemil, mapezmax, mapezmin );
// 	  

	  Mesh3 *Th3= build_layer(Th, nlayer, ni, zmin, zmax, maptet, maptrimil, maptrizmax, maptrizmin, mapemil, mapezmax, mapezmin);
    
  
  if( !(xx) && !(yy) && !(zz) ){
	 /*
	  map< int, int > maptet; 
	  map< int, int > maptrimil, maptrizmax, maptrizmin;
	  map< int, int > mapemil, mapezmax, mapezmin;
	   
	  build_layer_map_tetrahedra( Th, maptet );
	  build_layer_map_triangle( Th, maptrimil, maptrizmax, maptrizmin );
	  build_layer_map_edge( Th, mapemil, mapezmax, mapezmin );
	  
	  Mesh3 *Th3= build_layer(Th, nlayer, ni, zmin, zmax, maptet, maptrimil, maptrizmax, maptrizmin, mapemil, mapezmax, mapezmin);
	  */
	  Th3->BuildBound();
	  Th3->BuildAdj();
	  Th3->Buildbnormalv();  
	  Th3->BuildjElementConteningVertex();
	  Th3->BuildGTree();
	  Th3->decrement();  
  
  	*mp=mps;
  	return Th3;
  }
  else{
	  
	  //Mesh3 *Th3= build_layer(Th, nlayer, ni, zmin, zmax);
	  
	  KN<double> txx(Th3->nv), tyy(Th3->nv), tzz(Th3->nv);
	  KN<int> takemesh(Th3->nv);
	  MeshPoint *mp3(MeshPointStack(stack)); 
	  
	  takemesh=0;  
	  Mesh3 &rTh3 = *Th3;
	  for (int it=0;it<Th3->nt;++it){
		  for( int iv=0; iv<4; ++iv){
			   int i=(*Th3)(it,iv);  
				if(takemesh[i]==0){
					mp3->setP(Th3,it,iv);
					if(xx){ txx[i]=GetAny<double>((*xx)(stack));}
					if(yy){ tyy[i]=GetAny<double>((*yy)(stack));}
					if(zz){ tzz[i]=GetAny<double>((*zz)(stack));}

					takemesh[i] = takemesh[i]+1;
				}
	 		}
		}
		
		int border_only = 0;
		int recollement_elem=0, recollement_border=1, point_confondus_ok=0;
		Mesh3 *T_Th3=Transfo_Mesh3( rTh3, txx, tyy, tzz, border_only, recollement_elem, recollement_border, point_confondus_ok);
		  
		
 		T_Th3->BuildBound();
  		T_Th3->BuildAdj();
  		T_Th3->Buildbnormalv();  
  		T_Th3->BuildjElementConteningVertex();
  		T_Th3->BuildGTree();
  		T_Th3->decrement();  
  
	 	*mp=mps;
		return T_Th3;
	}
}


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
  
  Dcl_Type<listMesh3>();
  typedef Mesh3 *pmesh3;
  
  if (verbosity)
    cout << " load: glumesh  " << endl;
  //cout << " je suis dans Init " << endl; 
  TheOperators->Add("+",new OneBinaryOperator_st< Op2_addmesh<listMesh,pmesh,pmesh>  >      );
  TheOperators->Add("+",new OneBinaryOperator_st< Op2_addmesh<listMesh,listMesh,pmesh>  >      );
  TheOperators->Add("=",new OneBinaryOperator< Op2_setmesh<pmesh*,pmesh*,listMesh>  >     );
  TheOperators->Add("<-",new OneBinaryOperator< Op2_setmesh<pmesh*,pmesh*,listMesh>  >     );
  
  TheOperators->Add("+",new OneBinaryOperator_st< Op3_addmesh<listMesh3,pmesh3,pmesh3>  >      );
  TheOperators->Add("+",new OneBinaryOperator_st< Op3_addmesh<listMesh3,listMesh3,pmesh3>  >      );
  TheOperators->Add("=",new OneBinaryOperator< Op3_setmesh<pmesh3*,pmesh3*,listMesh3>  >     );
  TheOperators->Add("<-",new OneBinaryOperator< Op3_setmesh<pmesh3*,pmesh3*,listMesh3>  >     );
  
  Global.Add("change","(",new SetMesh);
  Global.Add("change3D","(",new SetMesh3D);
  Global.Add("buildlayers","(",new  BuildLayerMesh);
  //  Global.Add("build2D3D","(",new Build2D3D);
}
