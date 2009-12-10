#define INMESH3 
// now this code is in msh3.cpp 
// F. HecHt .. 
#ifndef WITH_NO_INIT
/*
#include <fstream>
#include <iostream>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <complex>
//#include <set>
#include<stdlib.h>
//#include <vector>
#include <map>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "ufunction.hpp"
using namespace std;  
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"


#include "FESpacen.hpp" 
#include "FESpace.hpp" 

#include "MatriceCreuse_tpl.hpp"
#include "MeshPoint.hpp"
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "Operator.hpp" 
#include "lex.hpp"
#include "libmesh5.h"
#include "lgfem.hpp"
#include "lgmesh3.hpp"
#include "lgsolver.hpp"
#include "problem.hpp"
*/
#include "ff++.hpp"
#endif

//  TransfoMesh_v2.cpp
using namespace std;
// LayerMesh.cpp
// buildlayer.cpp
// rajout global

#include <set>
#include <vector>
#include <list>
#include "../src/femlib/splitsimplex.hpp"
#include "msh3.hpp"

using namespace  Fem2D;

Mesh3 * truncmesh(const Mesh3 &Th,const long &kksplit,int *split, bool kk, const int newbelabel);

struct Op_trunc_mesh3 : public OneOperator {
  typedef Mesh3 *pmesh3;
  class Op: public E_F0mps   { 
  public:
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
  { return new Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1])); }
  Op_trunc_mesh3() : 
    OneOperator(atype<pmesh3>(),atype<pmesh3>(),atype<bool>()) {};     
};

basicAC_F0::name_and_type Op_trunc_mesh3::Op::name_param[Op_trunc_mesh3::Op::n_name_param] =
 {
   {  "split",             &typeid(long)},
   {  "label",             &typeid(long)}
 
 };


Mesh3 * truncmesh(const Mesh3 &Th,const long &kksplit,int *split, bool kk, const int newbelabel){
  
  static const int FaceTriangle[4]={3,0,1,2};  //={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}


  // computation of number of border elements and vertex without split
  int nbe = 0;
  int nt  = 0;
  int nv  = 0;
  int nvtrunc =0;
  int nbedge=0;
  double hmin=1e100;
  R3 bmin,bmax;
  
  const int kksplit2 = kksplit*kksplit; 
  const int kksplit3 = kksplit2*kksplit; 

  
  for (int i=0;i<Th.nt;i++)
    if(split[i]) 
      {	
	// computation of number of tetrahedrons 
	nt=nt+kksplit3;  
	// computation of number of border elements
	for (int j=0;j<4;j++)
	  {
	    int jt=j,it=Th.ElementAdj(i,jt);
	    if(it==i || it <0) nbe += kksplit2;  //on est sur la frontiere
	    else if (!split[it]) nbe += kksplit2; //le voisin ne doit pas etre decoupe
	    //else{
	    // rien a faire
	    //}
	  }

	for (int e=0;e<6;e++){
	  hmin=min(hmin,Th[i].lenEdge(e));   // calcul de .lenEdge pour un Mesh3
	}
      }
 
  double hseuil = (hmin/kksplit)/10.; 

  /* determination de bmin, bmax et hmin */ 
 
  KN<int> takevertex(Th.nv);
  for(int i=0; i<Th.nv; i++){
    takevertex[i]=-1;
  } 
  for(int i=0; i<Th.nt; i++){
    if(split[i]){
      const Tet &K(Th.elements[i]);
	  
      for(int ii=0; ii<4; ii++){
	int iv= Th.operator()( K[ii] );
	if( takevertex[iv] == -1 ){
	  R3 P( Th.vertices[iv].x, Th.vertices[iv].y, Th.vertices[iv].z);
	  bmin=Minc(P,bmin);
	  bmax=Maxc(P,bmax); 
	  takevertex[iv]=nvtrunc;
	  nvtrunc++;
	}
      }
      
    }
  }
  
  if( kksplit > 1 ){
    KN<int> NbFois(nvtrunc,0);
    for(int i=0; i<nvtrunc; i++)
      NbFois[i] = 0;
    int majedge=0;
    for(int i=0; i<Th.nt; i++){
      if(split[i]){
	const Tet &K(Th.elements[i]);
	for(int e=0; e<6; e++){
	  int e1 = Th.operator()( K[ Th[i].nvedge[e][0] ] );
	  int e2 = Th.operator()( K[ Th[i].nvedge[e][1] ] );
	  
	  if( takevertex[e1] < takevertex[e2] ){
	    NbFois[ takevertex[e1] ]++;
	    majedge++;
	  }
	  else{
	    NbFois[ takevertex[e2] ]++;
	    majedge++;
	  }
	}

      }
    }

    KN<int> first(nvtrunc+1,0);
    KN<int> current(nvtrunc,0);
    for(int i=1; i<nvtrunc; i++)
      first[i] = first[i-1] + NbFois[i-1];
    
    first[nvtrunc] = majedge;
    for(int i=0; i<nvtrunc; i++)
      current[i] = first[i];

    cout << "majedge =" << majedge << endl;

    KN<int> tableau(majedge);
    for(int i=0; i<Th.nt; i++){
      if(split[i]){
	const Tet &K(Th.elements[i]);
	for(int e=0; e<6; e++){
	  int e1 = Th.operator()( K[ Th[i].nvedge[e][0] ] );
	  int e2 = Th.operator()( K[ Th[i].nvedge[e][1] ] );
	  //int e1 = Th[i].nvedge[e][0];
	  //int e2 = Th[i].nvedge[e][1];
	  
	  if( takevertex[e1] < takevertex[e2] ){
	    tableau[ current[ takevertex[e1] ] ] = takevertex[e2]; 
	    current[ takevertex[e1] ]++;
	  }
	  else{
	    tableau[ current[ takevertex[e2] ] ] = takevertex[e1];
	    current[ takevertex[e2] ]++;
	  }
	  
	}
      }
    }
    for(int i=0; i<nvtrunc; i++)
      assert(current[i] == first[i+1]);

    // determination du nombre d'edge
    
    for(int i=0; i<nvtrunc; i++){
      list<int> list1;
      list<int>::const_iterator ilist;
      //cout << "i = "<< i  << " nvtrunc = " <<  nvtrunc << endl;
      for(int jj=first[i]; jj < first[i+1]; jj++){
	int labOk = 0;
	//cout << "tableauu[ jj ] = " <<  tableau[ jj ] << endl;
	for( ilist=list1.begin(); ilist!=list1.end(); ilist++){
	  if( *ilist == tableau[ jj ] ){ labOk = 1;   break; }	  
	}
	if(labOk == 0){
	  list1.push_back(tableau[jj]);
	  //cout << "push back :: tableauu[ jj ] = " <<  tableau[ jj ] << endl;
	  nbedge++;
	}

      }
      //cout << " nbedge="<< nbedge  << endl;
    }
  }



  /* determination des vertex, triangles et tetrahèdre obtenue après splitting dans le Simplex */ 

  int nfacesub = kksplit2;
  int ntetsub = kksplit3;
  int nvsub = (kksplit+1)*(kksplit+2)*(kksplit+3)/6;
  int ntrisub = 4*kksplit2;
   
  R3 *vertexsub; //[nvsub];
  int *tetsub;   //[4*ntetsub];
  int *trisub;   //[4*kksplit*kksplit];

  SplitSimplex<R3>( kksplit, nvsub, vertexsub, ntetsub, tetsub);
  SplitSurfaceSimplex( kksplit, ntrisub, trisub);
  
  for( int iii=0; iii< nvsub; iii++){
    cout << "vertexsub["<< iii <<"]=" << " "<< vertexsub[iii].x << " "<< vertexsub[iii].y << " "<< vertexsub[iii].z << endl;
  }
  for( int iii=0; iii< ntetsub; iii++){
    cout << "tetsub=" << tetsub[4*iii] << " "<< tetsub[4*iii+1] << " "<< tetsub[4*iii+2] << " "<< tetsub[4*iii+3] <<endl;
  }

  for( int iii=0; iii< 4*kksplit2; iii++){
    cout << iii << " tetsub=" << trisub[3*iii] << " "<< trisub[3*iii+1] << " " << trisub[3*iii+2] << endl;
  }
  
  cout << "Th.nv= " << Th.nv << "kksplit="<< kksplit << endl;
  // determination de nv 
  /*if(kksplit == 1)
    nv=nvtrunc;
  else{
  */
    int ntnosplit  = nt/kksplit3;
    int nbenosplit = nbe/kksplit2;
    int nfacenosplit = (4*ntnosplit+nbenosplit)/2;
    nv = ntnosplit*(nvsub - 4*( (kksplit+1)*(kksplit+2)/2 - 3*(kksplit-1) -3 ) - 6*( kksplit-1 ) - 4);
    cout << " nv= " << nv << endl;
    nv = nv + nfacenosplit*( (kksplit+1)*(kksplit+2)/2 - 3*(kksplit-1) -3 );
    cout << " nv= " << nv << endl;
    nv = nv + nbedge*( kksplit-1 );
    cout << " nv= " << nv << endl;
    nv = nv + nvtrunc;
    cout << " nv= " << nv << endl;
    //}
  

  int itt=0; 
  int ie=0;

 
  Vertex3 *v=new Vertex3[nv];
  Tet *t = new Tet[nt];
  Tet *tt = t;
  
  Triangle3 *b  = new Triangle3[nbe];
  Triangle3 *bb = b;

  EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,bmin,bmax,0);
  /*
  //determination des sommets
  int np=0;
  for(int i=0; i<Th.nt; i++){
    if(split[i]){      
      const Tet &K(Th.elements[i]);
      int ivv[4]; 
      for(int ii=0; ii< 4; ii++){
	ivv[ii] =Th.operator()(K[ii]);
      }

      Vertex3 vertextetsub[nvsub];
      for( int iv=0; iv<nvsub; iv++){
	double alpha=vertexsub[iv].x;
	double beta=vertexsub[iv].y;
	double gamma=vertexsub[iv].z;

	vertextetsub[iv].x = (1-alpha-beta-gamma)*Th.vertices[ivv[0]].x + alpha*Th.vertices[ivv[1]].x + beta*Th.vertices[ivv[2]].x + gamma*Th.vertices[ivv[3]].x;  
	vertextetsub[iv].y = (1-alpha-beta-gamma)*Th.vertices[ivv[0]].y + alpha*Th.vertices[ivv[1]].y + beta*Th.vertices[ivv[2]].y + gamma*Th.vertices[ivv[3]].y;
	vertextetsub[iv].z = (1-alpha-beta-gamma)*Th.vertices[ivv[0]].z + alpha*Th.vertices[ivv[1]].z + beta*Th.vertices[ivv[2]].z + gamma*Th.vertices[ivv[3]].z;  
	vertextetsub[iv].lab = K.lab;
      }
      

      int newindex[nvsub];
      for( int iv=0; iv<nvsub; iv++){
	//onst R3 viR3(vertextetsub[iv].x,vertextetsub[iv].y,vertextetsub[iv].z);
	const Vertex3 &vi( vertextetsub[iv] );
	Vertex3 * pvi=gtree->ToClose(vi,hseuil);
	 
	if(!pvi){
	  v[np].x   = vi.x;
	  v[np].y   = vi.y;
	  v[np].z   = vi.z;
	  v[np].lab = K.lab;
	  newindex[iv] = np;
	  gtree->Add( v[np] );
	  np++;
	}
	else{
	  newindex[iv] = pvi-v;	  
	}
      }
    }
  }
  
  if(np !=nv) cout << "np=" << np << " nv=" << nv << " nvtrunc="<< nvtrunc << endl; 
  assert( np == nv );
  */

  int np=0;
  for(int i=0; i<Th.nt; i++){
    if(split[i]){      
      const Tet &K(Th.elements[i]);
      int ivv[4]; 
      for(int ii=0; ii< 4; ii++){
	ivv[ii] =Th.operator()(K[ii]);
      }

      R3 vertextetsub[nvsub];
      for( int iv=0; iv<nvsub; iv++){
	double alpha=vertexsub[iv].x;
	double beta=vertexsub[iv].y;
	double gamma=vertexsub[iv].z;

	vertextetsub[iv].x = (1-alpha-beta-gamma)*Th.vertices[ivv[0]].x + alpha*Th.vertices[ivv[1]].x + beta*Th.vertices[ivv[2]].x + gamma*Th.vertices[ivv[3]].x;  
	vertextetsub[iv].y = (1-alpha-beta-gamma)*Th.vertices[ivv[0]].y + alpha*Th.vertices[ivv[1]].y + beta*Th.vertices[ivv[2]].y + gamma*Th.vertices[ivv[3]].y;
	vertextetsub[iv].z = (1-alpha-beta-gamma)*Th.vertices[ivv[0]].z + alpha*Th.vertices[ivv[1]].z + beta*Th.vertices[ivv[2]].z + gamma*Th.vertices[ivv[3]].z;  
      }
    
      int newindex[nvsub];
      for( int iv=0; iv<nvsub; iv++){
	const R3 viR3(vertextetsub[iv].x,vertextetsub[iv].y,vertextetsub[iv].z);
	const Vertex3 &vi( viR3 );
	Vertex3 * pvi=gtree->ToClose(vi,hseuil);
	
	if(!pvi){
	  v[np].x   = vi.x;
	  v[np].y   = vi.y;
	  v[np].z   = vi.z;
	  v[np].lab = K.lab;
	  newindex[iv] = np;
	  gtree->Add( v[np] );
	  np++;
	  }
	else{
	  newindex[iv] = pvi-v;
	}
	
	  //assert(pvi);
	  //newindex[iv] = pvi-v;
	if(np>nv) cout << "np=" << np << " nv=" << nv << endl; 
	assert( np <= nv );
      }

      for( int ii=0; ii<ntetsub; ii++){
	int ivt[4];
	for( int jj=0; jj< 4; jj++){
	  ivt[jj] = newindex[tetsub[4*ii+jj]];
	  assert( tetsub[4*ii+jj] < nvsub );
	  assert( ivt[jj] < np );
	}
	(tt++)->set( v, ivt, K.lab);
	itt++;
	assert( itt <= nt );
      }
      
      for (int j=0;j<4;j++)
	{
	  int jt=j,it=Th.ElementAdj(i,jt);
	 
	  /*
	  if( (it==i || it <0) ){

	  
	    int ivb[3];
	    int label = 32;

	    for( int ii=0; ii<nfacesub; ii++){

	      int iface = 3*FaceTriangle[j]*nfacesub+3*ii;
	      if(verbosity > 1) cout << "face " << ie << "iface " << iface << " " << FaceTriangle[j] 
				     << " " << trisub[iface] 
				     << " " << trisub[iface+1]  
				     << " " << trisub[iface+2] << endl; 
	      
	      for( int jjj=0; jjj<3; jjj++){
		ivb[jjj] = newindex[ trisub[iface+jjj] ]; 
		assert( trisub[ iface+jjj ] < nvsub );
		if(verbosity > 1) cout << ivb[jjj] << " np:" << np<< endl;
		assert( ivb[jjj] < np );
	      }
	    	      
	      (bb++)->set( v, ivb, label);
	      ie++;
	     
	    }
	  }
	  else 
	  */
	  if ( !(it==i || it <0)  && !split[it]) { 
	    int ivb[3];
	    
	    for( int ii=0; ii<nfacesub; ii++){
	      int iface = 3*FaceTriangle[j]*nfacesub+3*ii;

	      for( int jjj=0; jjj<3; jjj++){
		ivb[jjj] = newindex[ trisub[iface+jjj] ];
		assert( trisub[ iface+jjj ] < nvsub );
		assert( ivb[jjj] < np );
	      }
	    
	      (bb++)->set( v, ivb, newbelabel);
	      ie++;
	    } 
	  }
	  assert( ie <= nbe);
       	}
      
    }
  }

  delete [] vertexsub; //[nvsub];
  delete [] tetsub;   //[4*ntetsub];
  delete [] trisub;   //[4*kksplit*kksplit];

  // split border elements
  int nv2Dsub   = (kksplit+1)*(kksplit+2)/4;
  int ntri2Dsub = kksplit2;
  R2 *vertex2Dsub; //[nvsub];  
  int *tri2Dsub;   //[4*kksplit*kksplit];

  SplitSimplex<R2>( kksplit, nv2Dsub, vertex2Dsub, ntri2Dsub, tri2Dsub);


  for( int ibe=0; ibe < Th.nbe; ibe++){  
    int iff;
    int it=Th.BoundaryElement(ibe,iff);
    
    if( split[it] == 0 ) continue;
    const Triangle3 &K(Th.be(ibe));
    int ivv[3];
    
    ivv[0] = Th.operator()(K[0]);
    ivv[1] = Th.operator()(K[1]);
    ivv[2] = Th.operator()(K[2]);

     R3 vertextrisub[nv2Dsub];
     for( int iv=0; iv<nv2Dsub; iv++){
       double alpha=vertex2Dsub[iv].x;
       double beta=vertex2Dsub[iv].y;
       
       vertextrisub[iv].x = (1-alpha-beta)*Th.vertices[ivv[0]].x + alpha*Th.vertices[ivv[1]].x + beta*Th.vertices[ivv[2]].x;  
       vertextrisub[iv].y = (1-alpha-beta)*Th.vertices[ivv[0]].y + alpha*Th.vertices[ivv[1]].y + beta*Th.vertices[ivv[2]].y;
       vertextrisub[iv].z = (1-alpha-beta)*Th.vertices[ivv[0]].z + alpha*Th.vertices[ivv[1]].z + beta*Th.vertices[ivv[2]].z;  
       
     }
     int newindex[nv2Dsub];
     for( int iv=0; iv<nv2Dsub; iv++){
       const Vertex3 &vi( vertextrisub[iv] );
       Vertex3 * pvi=gtree->ToClose(vi,hseuil);
       assert(pvi);
       newindex[iv] = pvi-v;
     }

     for( int ii=0; ii<nfacesub; ii++){
       int ivb[3];
       for( int jjj=0; jjj<3; jjj++){
	 ivb[jjj] = newindex[ tri2Dsub[3*ii+jjj] ]; 
	 assert( tri2Dsub[ 3*ii+jjj  ] < nvsub );
	 if(verbosity > 1) cout << ivb[jjj] << " np:" << np<< endl;
	 assert( ivb[jjj] < np );
       }
       
       (bb++)->set( v, ivb, K.lab);
       ie++;
       
     }

     
  }

  delete [] vertex2Dsub;   //[4*ntetsub];
  delete [] tri2Dsub;   //[4*kksplit*kksplit];

  

  cout << "nbofv initial" << Th.nv << endl; 
  cout << "nv=" << nv << " np=" << np << endl; 
  assert( nv == np );
  cout << "itt=" << itt << " nt=" << nt << endl;
  assert( itt == nt );
  cout << "ie=" << ie << " nbe=" << nbe << endl;
  assert( ie ==nbe);

  //delete gtree;

  Mesh3 *Tht = new Mesh3( nv, nt, nbe, v, t, b); 
  
  delete gtree;
  
 
  return Tht;
}


AnyType Op_trunc_mesh3::Op::operator()(Stack stack)  const { 
    
  Mesh3 *pTh = GetAny<Mesh3 *>((*getmesh)(stack));
  Mesh3 &Th = *pTh;
  long kkksplit =arg(0,stack,1L);
  long label =arg(1,stack,2L);
  KN<int> split(Th.nt);
  split=kkksplit;
  MeshPoint *mp= MeshPointStack(stack),mps=*mp;
  long kk=0;
    
  for (int k=0;k<Th.nt;k++)
    { 
      const Tet & K( Th.elements[k] );
      R3 B(1./3.,1./3.,1./3.);
      mp->set(Th,K(B),B,K,0);
      if (  GetAny<bool>( (*bbb)(stack) )  ) kk++;
      else  split[k]=0  ;    
    }
  //*mp=mps;
  //if (verbosity>1) 
    cout << "  -- Trunc mesh: Nb of Tetrahedrons = " << kk << " label=" <<label <<endl;
  Mesh3 * Tht = truncmesh(Th,kkksplit,split,false,label);


 //  cout << Tht->nv << " " << Tht->nt << " " << Tht->nbe << " " << endl;
//   cout << "==================================" <<  Tht << endl;
//   exit(1);
//   for( int jjj=0; jjj<Tht->nv; jjj++){
//     cout << "apres trunc mesh :: vertex" << jjj+1 <<" " <<  Tht->vertices[jjj].x << " " << Tht->vertices[jjj].y << " " << Tht->vertices[jjj].z << " " << Tht->vertices[jjj].lab << endl;
    //if( abs(Tht->vertices[jjj].lab) > 1 ) exit(1); 
  //  }
  string filename("Thtpp_res.mesh");
  Tht->Save(filename); 

  cout << "==================================" <<  Tht << endl;


  Add2StackOfPtr2FreeRC(stack,Tht);//  07/2008 FH 
  *mp=mps;
  return Tht; 
 };
 
#ifndef WITH_NO_INIT
class Init { public:
  Init();
};

static Init init;  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  typedef Mesh3 *pmesh3;
  Global.Add("trunc","(", new Op_trunc_mesh3);

}

#endif
#endif
