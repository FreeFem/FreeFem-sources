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
#include "LayerMesh.hpp"
#include "TransfoMesh_v2.hpp"
//#include "GQuadTree.hpp"

#include <set>
#include <vector>
#include <fstream>


using namespace  Fem2D;


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
  
  for(list<Mesh3 *>::const_iterator i=lth.begin();i != lth.end();i++)
    {
      Mesh3 &Th3(**i);  // definis ???
      th0=&Th3;
      if(verbosity>1)  cout << " determination of hmin : GluMesh3D + "<< Th3.nv << " " << Th3.nt << " "<< Th3.nbe << endl;
      
      nbt  += Th3.nt;
      nbvx += Th3.nv;
      nbex += Th3.nbe;
      
      for (int k=0;k<Th3.nt;k++){
	for (int e=0;e<6;e++){
	  hmin=min(hmin,Th3[k].lenEdge(e));   // calcul de .lenEdge pour un Mesh3
	}
      }
      
      for (int k=0;k<Th3.nbe;k++){
	for (int e=0;e<3;e++){
	  hmin=min(hmin,Th3.be(k).lenEdge(e));   // calcul de .lenEdge pour un Mesh3
	}
      }
      
      for (int ii=0;ii<Th3.nv;ii++){ 
	R3 P( Th3.vertices[ii].x, Th3.vertices[ii].y, Th3.vertices[ii].z);
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
  
  
  //cout << " creation of : BuildGTree for border elements" << endl;
  //Vertex3  *becog= new Vertex3[nbex];  
  //EF23::GTree<Vertex3> *gtree_be = new EF23::GTree<Vertex3>(becog,Pn,Px,0);
  
  cout << " creation of : Mesh3 th3_tmp" << endl;
  /*Mesh3 th3_tmp;
    th3_tmp.set(nbvx,nbt,nbex);*/
  {  
    EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,Pn,Px,0);  
    typedef SortArray<int,3>  SortFace;
    HashTable<SortFace,int> bbe(nbex,nbvx);
    typedef HashTable<SortArray<int,3>,int>::iterator HT_iterator;
    for(list<Mesh3 *>::const_iterator i=lth.begin(); i!=lth.end();i++)
      {
	const Mesh3 &Th3(**i);
	if(verbosity>1)  cout << " loop over mesh for create new mesh "<< endl;
	if(verbosity>1)  cout << " GluMesh3D + "<< Th3.nv << " " << Th3.nt <<" " << Th3.nbe << endl;
	int nbv0 = nbv;
	
	for (int ii=0;ii<Th3.nv;ii++)
	  {
	    const Vertex3 &vi(Th3.vertices[ii]);
	    Vertex3 * pvi=gtree->ToClose(vi,hseuil);
       
	    //  if( abs(vi.x) <1e-6 || abs(vi.x-1) <1e-6 || abs(vi.x+1) <1e-6 ){
	    //   cout << ii << *i << "avant test" << vi << endl;	   }
	    
	    if(!pvi)
	      {
		if(abs(vi.x) <1e-6 || abs(vi.x-1) <1e-6 || abs(vi.x+1) <1e-6){
		  cout << "nbv=" << nbv << " vi=" << vi << endl;
		}
		v[nbv].x = vi.x;
		v[nbv].y = vi.y;
		v[nbv].z = vi.z;
		v[nbv].lab = vi.lab;
		gtree->Add( v[nbv++] );
	      }
	  }
	
	for (int k=0;k<Th3.nt;k++)
	  {
	    const Tet  &K(Th3.elements[k]);
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
	    SortFace sbe(iv);
	    if(bbe.find(sbe)) continue;
	    //if( iv[2]<nbv0 && iv[1]<nbv0 && iv[0] < nbv0 ) continue;  // Faux
	    (bb++)->set(v,iv,K.lab);
	      nbe++;
	      
	  }   
	cout << "nbe="<<nbe<<endl;   
	// Remarque on a besoin 
	
      }
      delete gtree;
  }
   /* 
  if(0) 
    {// old code  ..
      cout << " creation of : BuildGTree for border elements" << endl;
      Vertex3  *becog= new Vertex3[nbex];  
      EF23::GTree<Vertex3> *gtree_be = new EF23::GTree<Vertex3>(becog,Pn,Px,0);
      
      double hseuil_border = hseuil/3.;
      
      for(list<Mesh3 *>::const_iterator i=lth.begin();i != lth.end();i++){
	const Mesh3 &Th3(**i);
	
	for (int k=0;k<Th3.nbe;k++)
	  {
	    const Triangle3 & K(Th3.be(k));//bedges[k]);
	    int iv[3];
	    iv[0]=Th3.operator()(K[0]); //gtree->NearestVertex(K[0])-v; //-v;
	    iv[1]=Th3.operator()(K[1]); //gtree->NearestVertex(K[1])-v; //-v;
	    iv[2]=Th3.operator()(K[2]); //gtree->NearestVertex(K[2])-v; //-v;
	    
	    //assert( iv[2]<nbv && iv[1]<nbv && iv[0] < nbv );  
	    
	    R cdgx,cdgy,cdgz;
	    
	    cdgx = (Th3.vertices[iv[0]].x+ Th3.vertices[iv[1]].x+ Th3.vertices[iv[2]].x)/3.;
	    cdgy = (Th3.vertices[iv[0]].y+ Th3.vertices[iv[1]].y+ Th3.vertices[iv[2]].y)/3.;
	    cdgz = (Th3.vertices[iv[0]].z+ Th3.vertices[iv[1]].z+ Th3.vertices[iv[2]].z)/3.;
	    
	    //cout << "cdg=" << cdgx << " " << cdgy << " " << cdgz << " "<<endl; 
	    const R3 r3vi( cdgx, cdgy, cdgz ); 
	    const Vertex3 &vi( r3vi);
	    
	    Vertex3 * pvi=gtree_be->ToClose(vi,hseuil_border);
	    if(!pvi){
	      becog[nbe].x = vi.x;
	      becog[nbe].y = vi.y;
	      becog[nbe].z = vi.z;
	      becog[nbe].lab = vi.lab;
	      gtree_be->Add( becog[nbe++]);
	      
	      int igluv[3];
	      igluv[0]=gtree->NearestVertex(K[0])-v; //-v;
	      igluv[1]=gtree->NearestVertex(K[1])-v; //-v;
	      igluv[2]=gtree->NearestVertex(K[2])-v; //-v;
	      
	      (bb++)->set(v,igluv,K.lab);
	    }
	  }
	
      }
      //delete th0;
    }
  */
  
  if(verbosity>1)
    {
      cout << "     Nb of glu3D  point " << nbvx-nbv;
      cout << "     Nb of glu3D  Boundary faces " << nbex-nbe << endl;
    }
  cout << " nbv="  << nbv  << endl;
  cout << " nbvx=" << nbvx << endl;
  cout << " nbt="  << nbt  << endl;
  cout << " nbe="  << nbe  << endl;
  cout << " nbex=" << nbex << endl;
  
  cout << "     Nb of glu3D  point " << nbvx-nbv;
  cout << "     Nb of glu3D  Boundary edge " << nbex-nbe << endl;
  
  
  
  Mesh3 *mpq= new Mesh3(nbv,nbt,nbe,v,t,b);  
  cout << "fin de la definition de mpq" << endl;
  
  if(nbt !=0){     
    mpq->BuildBound();
    cout << "fin de BuildBound" << endl;
    mpq->BuildAdj();
    cout << "fin de BuildAdj" << endl;
    mpq->Buildbnormalv();  
    cout << "fin de Buildnormalv()" << endl;
    mpq->BuildjElementConteningVertex();
    cout << "fin de ConteningVertex()" << endl;
    mpq->BuildGTree();
    cout << "fin de BuildGTree()" << endl;
  }
  return mpq;
  
}


template<class RR,class AA=RR,class BB=AA> 
struct Op3_addmesh: public binary_function<AA,BB,RR> { 
  static RR f(Stack s,const AA & a,const BB & b)  
  { return RR(s, a, b );} 
};

template<bool INIT,class RR,class AA=RR,class BB=AA> 
struct Op3_setmesh: public binary_function<AA,BB,RR> { 
  static RR f(const AA & a,const BB & b)  
  {
    ffassert(a );
    pmesh3  p=GluMesh3(b);
    
    if(!INIT &&  *a) delete *a;
    //  Add2StackOfPtr2FreeRC(stack,p); //  the pointer is use to set variable so no remove. 
    return *a=p,a;
  } 
};

// Movemesh3D

class Movemesh3D_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression xx,yy,zz;
  static const int n_name_param =3; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  
public:
  Movemesh3D_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth), xx(0) , yy(0) , zz(0)
  {
    args.SetNameParam(n_name_param,name_param,nargs);
    const E_Array * a1=0 ;
    if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[0]);
    int err =0;
 
    if(a1) {
      if(a1->size() !=3) 
      CompileError("Build2D3D (Th,transfo=[X,Y,Z],) ");
      xx=to<double>( (*a1)[0]);
      yy=to<double>( (*a1)[1]);
      zz=to<double>( (*a1)[2]);
    }    
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type Movemesh3D_Op::name_param[]= {
  {  "transfo", &typeid(E_Array)},
  {  "reftet", &typeid(KN_<long> )},
  {  "refface", &typeid(KN_<long> )}
  
};


class  Movemesh3D : public OneOperator { public:  
    Movemesh3D() : OneOperator(atype<pmesh3>(),atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	cout << " je suis dans code(const basicAC_F0 & args) const" << endl;
	return  new Movemesh3D_Op(args,t[0]->CastTo(args[0])); 
  }
};

AnyType Movemesh3D_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  
  ffassert(pTh);
  Mesh3 &Th=*pTh;
  Mesh3 *m= pTh;   // question a quoi sert *m ??
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.nbe; // nombre d'aretes fontiere
  cout << " " << nbv<< " "<< nbv << " nbe "<< nbe << endl; 
 
 // lecture des references
   
  KN<long> zzempty;
  KN<long> nrtet  (arg(1,stack,zzempty));  
  KN<long> nrf (arg(2,stack,zzempty)); 

  cout << nrtet.N() << nrf.N() << endl;
  
  //if( nrtet.N() && nrfmid.N() && nrfup.N() && nrfdown.N() ) return m;
  ffassert( nrtet.N() %2 ==0);
  ffassert( nrf.N() %2 ==0);
  

  // realisation de la map par default
  
    assert((xx) && (yy) && (zz) );
  
	KN<double> txx(Th.nv), tyy(Th.nv), tzz(Th.nv);
	Mesh3 &rTh3 = Th;
	{
	KN<int> takemesh(Th.nv);
	MeshPoint *mp3(MeshPointStack(stack)); 
	  
	takemesh=0;
	for (int it=0;it<Th.nt;++it){
		  for( int iv=0; iv<4; ++iv){
			   int i=Th(it,iv);  
				if(takemesh[i]==0){
					mp3->setP(&rTh3,it,iv);
					if(xx){ txx[i]=GetAny<double>((*xx)(stack));}
					if(yy){ tyy[i]=GetAny<double>((*yy)(stack));}
					if(zz){ tzz[i]=GetAny<double>((*zz)(stack));}
					takemesh[i] = takemesh[i]+1;
				}
			}
		}
	}
	
	int border_only = 0;
	int recollement_elem=0, recollement_border=1, point_confondus_ok=0;
	Mesh3 *T_Th3=Transfo_Mesh3( rTh3, txx, tyy, tzz, border_only, 
		recollement_elem, recollement_border, point_confondus_ok);
	  
	/*  changement de label   */
	
	
	T_Th3->BuildBound();
	T_Th3->BuildAdj();
	T_Th3->Buildbnormalv();  
	T_Th3->BuildjElementConteningVertex();
	T_Th3->BuildGTree();
	//	T_Th3->decrement();  
	Add2StackOfPtr2FreeRC(stack,T_Th3);
  
	*mp=mps;
	return T_Th3;
}

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
//  besoin en cas de fichier 2D / fichier 3D 

int  ChangeLab3D(const map<int,int> & m,int lab)
{
  map<int,int>::const_iterator i=m.find(lab);
  if(i != m.end())
    lab=i->second;
  return lab;
}

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
    int l0,l1=ChangeLab3D(mapface, l0= K.lab) ;
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
      
      (*tt++).set( v, iv, ChangeLab3D(maptet,lab));
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
      
      int l0,l1=ChangeLab3D(mapface,l0=K.lab) ;
      (*bb++).set( v, iv, l1);   
    }
    
  assert(nben==bb-b);
  
  Mesh3 *mpq = new Mesh3(nbv,nbt,nbe,v,t,b);
  
  mpq->BuildBound();
  mpq->BuildAdj();
  mpq->Buildbnormalv();  
  mpq->BuildjElementConteningVertex(); 
  mpq->BuildGTree();
  //mpq->decrement();
  Add2StackOfPtr2FreeRC(stack,mpq);

  
  return mpq;
}


class SetMesh3D : public OneOperator { public:  
typedef Mesh3 *pmesh3;
    SetMesh3D() : OneOperator(atype<pmesh3>(),atype<pmesh3>() ) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new SetMesh3D_Op(args,t[0]->CastTo(args[0])); 
  }
};

// ---------------------------------
// Movemesh2d_3D_surf


class Movemesh2D_3D_surf_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression xx,yy,zz; 
  static const int n_name_param =3; //  add nbiter FH 30/01/2007 11 -> 12 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}

  
public:
  Movemesh2D_3D_surf_Op(const basicAC_F0 &  args,Expression tth) : 
  eTh(tth),xx(0),yy(0),zz(0) 
  {
  
    args.SetNameParam(n_name_param,name_param,nargs);
    
    const E_Array * a1=0 ;
    if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[0]);
    int err =0;
    
    if(a1) {
      if(a1->size() !=3) 
	CompileError("Movemesh2D_3D_surf (Th,transfo=[X,Y,Z],) ");
      xx=to<double>( (*a1)[0]);
      yy=to<double>( (*a1)[1]);
      zz=to<double>( (*a1)[2]);
    }    
    
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type Movemesh2D_3D_surf_Op::name_param[]= {
  {  "transfo", &typeid(E_Array )},
  {  "refface", &typeid(KN_<long>)},
  {  "mesuremesh", &typeid(KN_<long>)}
};

AnyType Movemesh2D_3D_surf_Op::operator()(Stack stack)  const 
{
  Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
  Mesh & Th=*pTh;
  Mesh *m= pTh;
  int nbv=Th.nv;  // nombre de sommet 
  int nbt=Th.nt;  // nombre de triangles
  int nbe=Th.neb; // nombre d'aretes fontiere
  
  cout << "Vertex Triangle Edge"<< nbv << nbt << nbe << endl;  
  
  KN<long> zzempty;
  KN<long> nrface (arg(1,stack,zzempty));
  KN<long> mesuremesh (arg(2,stack,zzempty));
  
  
  if(nrface.N()<0 ) return m;
  ffassert( nrface.N() %2 ==0);
  
  map<int,int> mapface;
  
  int z00 = false;
  for(int i=0;i<nrface.N();i+=2)
    { z00 = z00 || ( nrface[i]==0 && nrface[i+1]==0);

      if( nrface[i] != nrface[i+1] ){	
      mapface[nrface[i]] = nrface[i+1];
      }      
    }
  
  int surface_orientation=1; 
  if( mesuremesh.N() >0  && mesuremesh[0] > 0){
	  surface_orientation=-1;
  }
        

  KN<double> txx(nbv), tyy(nbv), tzz(nbv);
  MeshPoint *mp3(MeshPointStack(stack)); 
  
  {
	KN<int> takemesh(nbv);  
	takemesh=0;  
	Mesh &rTh = Th;
	for (int it=0; it<nbt; ++it){
		  for( int iv=0; iv<3; ++iv){
				int i=Th(it,iv);
				if(takemesh[i]==0){
					mp3->setP(&Th,it,iv);
					if(xx){ 
						txx[i]=GetAny<double>((*xx)(stack));
					}
					if(yy){ 
						tyy[i]=GetAny<double>((*yy)(stack));
					}
					if(zz){ 
						tzz[i]=GetAny<double>((*zz)(stack));
					}
					takemesh[i] = takemesh[i]+1;
				}
			}
		}
  }
  
  Mesh3 *Th3= new Mesh3;
  
  int vertex_out=1;
  
  if( vertex_out == 1){
	/* determinate the same vertex */ 
	int border_only = 0;
	int recollement_border=1, point_confondus_ok=0;

	// faire version de Transfo_Mesh2_tetgen pour ce cas précis.
	Th3= MoveMesh2_func( Th, txx, tyy, tzz, 
		border_only, recollement_border, point_confondus_ok);
	
	// Rajouter fonction flip a l interieure
	
	for(int ii=0; ii < Th3->nbe; ii++){
		
		const Triangle3 & K(Th3->be(ii)); // const Triangle2 & K(Th2.elements[ii]); // Version Mesh2  
		int iv[3];
		int lab;
		double mes_triangle3;
		
		iv[0] = Th3->operator()(K[0]);
		iv[1] = Th3->operator()(K[1]);
		iv[2] = Th3->operator()(K[2]);
		
		map< int, int>:: const_iterator imap;
		imap = mapface.find(K.lab); 
		
		if(imap!= mapface.end()){
			 lab=imap->second;
		} 
		else{
			lab=K.lab;
		}
		
		Th3->be(ii).set( Th3->vertices, iv, lab ) ;
		mes_triangle3 = Th3->be(ii).mesure();
		
		if( surface_orientation*mes_triangle3 < 0){
			int iv_temp=iv[1];
			iv[1]=iv[2];
			iv[2]=iv_temp;
			Th3->be(ii).set( Th3->vertices, iv, lab ) ;
		}
		
		/* autre methode a tester */
		/*
		Triangle3 Kmes;
		Kmes.set( Th3->vertices, iv, lab ) ;
		mes_triangle3 = Kmes.mesure();
		if( surface_orientation*mes_triangle3) < 0){
			int iv_temp=iv[1];
			iv[1]=iv[2];
			iv[2]=iv_temp;
		}
		Th3->be(ii).set( Th3->vertices, iv, lab ) ;
		*/
	}
  }
  
  if( vertex_out == 0){
	  
	  Tet       *t;
	  Vertex3   *v = new Vertex3[nbv];
	  Triangle3 *b = new Triangle3[nbe];
	  // generation des nouveaux sommets 
	  Vertex3 *vv=v;
	 // copie des anciens sommets (remarque il n'y a pas operateur de copy des sommets)
	  for (int i=0;i<nbv;i++)
	  {
		  const Mesh::Vertex & V( Th.vertices[i]);
		  vv->x = txx[i];
		  vv->y = tyy[i];
		  vv->z = tzz[i];
		  vv->lab = V.lab;
		  vv++;      
	  }
 
	// les arete frontieres qui n'ont pas change
  
	  Triangle3 * bb=b;
	  for (int i=0;i<nbt;i++)
	  { 
      	const Mesh::Triangle &K( Th.t(i) );
      	int iv[3];       
    
      	iv[0] = Th.operator()(K[0]);
      	iv[1] = Th.operator()(K[1]);
      	iv[2] = Th.operator()(K[2]);
      	
      // calcul de la mesure 
      //  inversement si on a besoin    
      
      	(*bb++).set( v, iv, K.lab);
      }
      
      Th3 = new Mesh3(nbv,0,nbt,v,t,b);  
      
  }
  
  return Th3;
}


class Movemesh2D_3D_surf : public OneOperator { public:  
typedef Mesh *pmesh;
typedef Mesh3 *pmesh3;
    
  Movemesh2D_3D_surf() : OneOperator(atype<pmesh3>(),atype<pmesh>() ) {}
    E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new Movemesh2D_3D_surf_Op(args,t[0]->CastTo(args[0]));  // CastTo(args[]); // plus tard
  }
};


class Init { public:
  Init();
};

static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  Dcl_Type<listMesh3>();
  typedef Mesh3 *pmesh3;
  
  if (verbosity)
    cout << " load: glumesh  " << endl;
  //cout << " je suis dans Init " << endl; 
  
  TheOperators->Add("+",new OneBinaryOperator_st< Op3_addmesh<listMesh3,pmesh3,pmesh3>  >      );
  TheOperators->Add("+",new OneBinaryOperator_st< Op3_addmesh<listMesh3,listMesh3,pmesh3>  >      );
  TheOperators->Add("=",new OneBinaryOperator< Op3_setmesh<false,pmesh3*,pmesh3*,listMesh3>  >     );
  TheOperators->Add("<-",new OneBinaryOperator< Op3_setmesh<true,pmesh3*,pmesh3*,listMesh3>  >     );
  
  Global.Add("change","(",new SetMesh3D);
  Global.Add("movemesh2D3Dsurf","(",new Movemesh2D_3D_surf);
  Global.Add("movemesh3D","(",new Movemesh3D);
  
}


