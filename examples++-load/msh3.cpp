// ORIG-DATE: Novembre 2008
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
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

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */

//  FH   July 2009
//   comment all
//       Th3_t->BuildBound();
//       Th3_t->BuildAdj();
//       Th3_t->Buildbnormalv();  
//       Th3_t->BuildjElementConteningVertex();
//   is now in the constructor of Mesh3 to be consistante. 
//

#ifndef WITH_NO_INIT // cf [[WITH_NO_INIT]]
#include "ff++.hpp"
#endif
#include "AFunction_ext.hpp" // [[file:../src/fflib/AFunction_ext.hpp]]

//  TransfoMesh_v2.cpp
using namespace std;
// LayerMesh.cpp
// buildlayer.cpp
// trunc3d.cpp
// rajout global
#include <climits>
#include <set>
#include <vector>
#include "msh3.hpp" // [[file:msh3.hpp]]
#include "splitsimplex.hpp" // [[file:../src/femlib/splitsimplex.hpp]]

using namespace  Fem2D;

int  ChangeLab3D(const map<int,int> & m,int lab);

void TestSameVertexMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int & nv_t, int *Numero_Som){
  
  Vertex3  *v=new Vertex3[Th3.nv];
  nv_t=0;
  
  EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,Pinf,Psup,0);
  
  // creation of octree
  for (int ii=0;ii<Th3.nv;ii++){
    const R3 r3vi( Th3.vertices[ii].x, Th3.vertices[ii].y, Th3.vertices[ii].z );
    const Vertex3 &vi(r3vi);
    
    Vertex3 * pvi=gtree->ToClose(vi,hseuil);
    
    if(!pvi){
      v[nv_t].x = vi.x;
      v[nv_t].y = vi.y;
      v[nv_t].z = vi.z;
      v[nv_t].lab = Th3.vertices[ii].lab; // lab mis a zero par default
      Numero_Som[ii] = nv_t; 
      gtree->Add( v[nv_t] );
      nv_t=nv_t+1;
    }
    else{
      Numero_Som[ii] = pvi-v;
    }
  }
  
  delete gtree;
  delete [] v;
}

void TestSameTetrahedraMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int & nt_t ){

  Vertex3 *vt=new Vertex3[Th3.nt];
  EF23::GTree<Vertex3> *gtree_t = new EF23::GTree<Vertex3>(vt,Pinf,Psup,0);
  
  nt_t=0;
  // creation of octree
  for (int ii=0;ii<Th3.nt;ii++){
    const Tet & K(Th3.elements[ii]);
    int iv[4];

    for(int jj=0; jj <4; jj++)
      iv[jj] = Th3.operator()(K[jj]);
    
    const double Cdg_x = ( Th3.vertices[iv[0]].x + Th3.vertices[iv[1]].x + Th3.vertices[iv[2]].x + Th3.vertices[iv[3]].x )/4.;
    const double Cdg_y = ( Th3.vertices[iv[0]].y + Th3.vertices[iv[1]].y + Th3.vertices[iv[2]].y + Th3.vertices[iv[3]].y )/4.;
    const double Cdg_z = ( Th3.vertices[iv[0]].z + Th3.vertices[iv[1]].z + Th3.vertices[iv[2]].z + Th3.vertices[iv[3]].z )/4.;

    const R3 r3vi( Cdg_x, Cdg_y, Cdg_z );
    const Vertex3 &vi(r3vi);
    
    Vertex3 * pvi=gtree_t->ToClose(vi,hseuil);
    
    if(!pvi){
      vt[nt_t].x = vi.x;
      vt[nt_t].y = vi.y;
      vt[nt_t].z = vi.z;
      vt[nt_t].lab = K.lab ; // lab mis a zero par default
      gtree_t->Add( vt[nt_t] );
      nt_t=nt_t+1;
    }
  }
  
  delete gtree_t;
  delete [] vt;
} 

void TestSameTetrahedraMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int *Elem_ok, int & nt_t ){

  Vertex3 *vt=new Vertex3[Th3.nt];
  EF23::GTree<Vertex3> *gtree_t = new EF23::GTree<Vertex3>(vt,Pinf,Psup,0);
  
  nt_t=0;
  // creation of octree
  for (int ii=0;ii<Th3.nt;ii++){
    if(Elem_ok[ii]!=1) continue;
    const Tet & K(Th3.elements[ii]);
    int iv[4];

    for(int jj=0; jj <4; jj++)
      iv[jj] = Th3.operator()(K[jj]);
    
    const double Cdg_x = ( Th3.vertices[iv[0]].x + Th3.vertices[iv[1]].x + Th3.vertices[iv[2]].x + Th3.vertices[iv[3]].x )/4.;
    const double Cdg_y = ( Th3.vertices[iv[0]].y + Th3.vertices[iv[1]].y + Th3.vertices[iv[2]].y + Th3.vertices[iv[3]].y )/4.;
    const double Cdg_z = ( Th3.vertices[iv[0]].z + Th3.vertices[iv[1]].z + Th3.vertices[iv[2]].z + Th3.vertices[iv[3]].z )/4.;

    const R3 r3vi( Cdg_x, Cdg_y, Cdg_z );
    const Vertex3 &vi(r3vi);
    
    Vertex3 * pvi=gtree_t->ToClose(vi,hseuil);
    
    if(!pvi){
      vt[nt_t].x = vi.x;
      vt[nt_t].y = vi.y;
      vt[nt_t].z = vi.z;
      vt[nt_t].lab = K.lab ; // lab mis a zero par default
      gtree_t->Add( vt[nt_t] );
      nt_t=nt_t+1;
    }
    else{
      Elem_ok[ii]=0;
    }
  }
  
  delete gtree_t;
  delete [] vt;
} 


void TestSameTriangleMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int & nbe_t){

  Vertex3  *vbe= new Vertex3[Th3.nbe];
  EF23::GTree<Vertex3> *gtree_be = new EF23::GTree<Vertex3>(vbe,Pinf,Psup,0);
  
  nbe_t=0;
  // creation of octree
  for (int ii=0;ii<Th3.nbe;ii++){
    const Triangle3 & K(Th3.be(ii));
    int iv[3];

    for(int jj=0; jj<3; jj++)
      iv[jj] = Th3.operator()(K[jj]);
    
    const double Cdg_x = ( Th3.vertices[iv[0]].x + Th3.vertices[iv[1]].x + Th3.vertices[iv[2]].x )/3.;
    const double Cdg_y = ( Th3.vertices[iv[0]].y + Th3.vertices[iv[1]].y + Th3.vertices[iv[2]].y )/3.;
    const double Cdg_z = ( Th3.vertices[iv[0]].z + Th3.vertices[iv[1]].z + Th3.vertices[iv[2]].z )/3.;

    const R3 r3vi( Cdg_x, Cdg_y, Cdg_z );
    const Vertex3 &vi(r3vi);
    
    Vertex3 * pvi=gtree_be->ToClose(vi,hseuil);
    
    if(!pvi){
      vbe[nbe_t].x = vi.x;
      vbe[nbe_t].y = vi.y;
      vbe[nbe_t].z = vi.z;
      vbe[nbe_t].lab = K.lab ; // lab mis a zero par default
      gtree_be->Add( vbe[nbe_t] );
      nbe_t=nbe_t+1;
    }
   
  }
  
  delete gtree_be;
  delete [] vbe;
} 

void TestSameTriangleMesh3( const Mesh3 & Th3, const double & hseuil, const R3 & Psup, const R3 &Pinf, int *Border_ok ,int & nbe_t ){

  Vertex3 *vbe=new Vertex3 [Th3.nbe];
  EF23::GTree<Vertex3> *gtree_be = new EF23::GTree<Vertex3>(vbe,Pinf,Psup,0);
  
  nbe_t=0;
  // creation of octree
  for (int ii=0;ii<Th3.nbe;ii++){
    if(Border_ok[ii]!=1) continue;
    const Triangle3 & K(Th3.be(ii));
    int iv[3];

    for(int jj=0; jj<3; jj++)
      iv[jj] = Th3.operator()(K[jj]);
    
    const double Cdg_x = ( Th3.vertices[iv[0]].x + Th3.vertices[iv[1]].x + Th3.vertices[iv[2]].x )/3.;
    const double Cdg_y = ( Th3.vertices[iv[0]].y + Th3.vertices[iv[1]].y + Th3.vertices[iv[2]].y )/3.;
    const double Cdg_z = ( Th3.vertices[iv[0]].z + Th3.vertices[iv[1]].z + Th3.vertices[iv[2]].z )/3.;

    const R3 r3vi( Cdg_x, Cdg_y, Cdg_z );
    const Vertex3 &vi(r3vi);
    
    Vertex3 * pvi=gtree_be->ToClose(vi,hseuil);
    
    if(!pvi){
      vbe[nbe_t].x = vi.x;
      vbe[nbe_t].y = vi.y;
      vbe[nbe_t].z = vi.z;
      vbe[nbe_t].lab = K.lab ; // lab mis a zero par default
      gtree_be->Add( vbe[nbe_t] );
      nbe_t=nbe_t+1;
    }
    else{      
      if(K.lab == vbe[pvi-vbe].lab )  Border_ok[ii] = 0;
    }
  }
  
  delete gtree_be;
  delete [] vbe;
} 

int TestElementMesh3( const Mesh3 & Th3 ) 
// Test si le maillage à des éléments communs : Sommet, triangle, ...
{
//  FH 31/09/2009:  Change  int* to KN<int> to remove pb of missing free in some case 
  R3 Pinf(1e100,1e100,1e100),Psup(-1e100,-1e100,-1e100);   // Extremité de la boîte englobante
  double hmin=1e10;   // longueur minimal des arrêtes
  double hseuil;
  KN<int> Numero_Som(Th3.nv);
  int nv_t,nt_t,nbe_t;
  
  // calcul de la boite englobante
  for (int ii=0;ii<Th3.nv;ii++){ 
    R3 P( Th3.vertices[ii].x, Th3.vertices[ii].y, Th3.vertices[ii].z);
    Pinf=Minc(P,Pinf);
    Psup=Maxc(P,Psup);     
  }


  // calcul de la longueur minimal des arrêtes
  for (int k=0;k<Th3.nt;k++){
    for (int e=0;e<6;e++){
      if( Th3[k].lenEdge(e) < Norme2(Psup-Pinf)/1e9 )
	{
	  const Tet & K(Th3.elements[k]);
	  int iv[4];
	  for(int jj=0; jj <4; jj++){
	    iv[jj] = Th3.operator()(K[jj]);
	  }
	  for (int eh=0;eh<6;eh++){
	    cout << "tetrahedra: " << k << " edge : " << eh << " lenght "<<  Th3[k].lenEdge(eh) << endl;
	    cout << " Tet vertices : " <<iv[0] << " " << iv[1] << " " << iv[2] << " " << iv[3] << " " << endl;
	  }
	  cout << " A tetrahedra with a very small edge was created " << endl;
	 
	  return 1;
	}
      hmin=min(hmin,Th3[k].lenEdge(e));   // calcul de .lenEdge pour un Mesh3
    }
  }
  
  for (int k=0;k<Th3.nbe;k++){
    for (int e=0;e<3;e++){
      if( Th3.be(k).lenEdge(e) < Norme2(Psup-Pinf)/1e9 )
	{
	  for (int eh=0;eh<3;e++){
	    cout << "triangles: " << k << " edge : " << eh << " lenght "<<  Th3.be(k).lenEdge(e) << endl;
	  }
	  cout << " A triangle with a very small edges was created " << endl;
	  return 1;
	}
      hmin=min(hmin,Th3.be(k).lenEdge(e));   // calcul de .lenEdge pour un Mesh3
    }
  }

  if(verbosity > 1) cout << "      - hmin =" <<  hmin << " ,  Bounding Box: " << Pinf << " "<< Psup << endl;

  ffassert(hmin>Norme2(Psup-Pinf)/1e9);

  // determination du nombre de sommets confondus
  hseuil = hmin/10.;     
  
  if(verbosity >1) cout << "TestSameVertexMesh3 " << hseuil << " size" <<Th3.nv<< endl;
  TestSameVertexMesh3( Th3, hseuil, Psup, Pinf, nv_t, Numero_Som );
  
  if(verbosity >1) cout << "hseuil=" << hseuil << endl; 
  if(verbosity >1) cout << "NbVertexRecollement " << nv_t << " / " << "NbVertex(anc)" << Th3.nv <<endl;   

  if(nv_t != Th3.nv) {
   // delete [] Numero_Som;
    cout << " A vertex was referenced twice or more " << endl;
    return 1;
  }
  /* degenerate element ??? */ 
  KN<int> Elem_ok(Th3.nt);
  int i_elem=0;
  for(int ii=0; ii< Th3.nt; ii++){
    const Tet & K(Th3.elements[ii]);
    int iv[4];
    
    Elem_ok[ii] = 1;
    
    for(int jj=0; jj <4; jj++){
      iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
    }
      
    for(int jj=0; jj<4; jj++){
      for(int kk=jj+1; kk<4; kk++){
	if( iv[jj]==iv[kk] ){
	  Elem_ok[ii] = 0;
	}
      }
    }
    i_elem = i_elem + Elem_ok[ii];
  }
   
  if( i_elem != Th3.nt ){
    cout << "There are a false tetrahedra in the mesh" << endl;
    assert( i_elem == Th3.nt);
  }

  KN<int> Border_ok(Th3.nbe);
  int i_border= 0;
  for( int ii=0; ii< Th3.nbe; ii++){
    Border_ok[ii]=1;
    
    const Triangle3 & K(Th3.be(ii));
    int iv[3];
    
    for(int jj=0; jj <3; jj++){
      iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
      assert( iv[jj] >= 0 && iv[jj] < nv_t);
    }
    
    for(int jj=0; jj<3; jj++){
      for(int kk=jj+1; kk<3; kk++){
	if( iv[jj]==iv[kk] ) Border_ok[ii]=0;
      }
    }
    i_border = i_border + Border_ok[ii];
  }
  
  if( i_border != Th3.nbe){
    cout << "There are a false tetrahedra in the mesh" << endl;
    assert( i_elem == Th3.nt);
  }

  /* determination du nombre de tetrahedre confondus */ 
  hseuil = hmin/10.;
  hseuil = hseuil/4.;
  
  if(verbosity >1) cout << "TestSameTetrahedraMesh3 " << hseuil << " size "<< Th3.nt <<endl;
  TestSameTetrahedraMesh3( Th3, hseuil, Psup, Pinf, Elem_ok, nt_t );

  if(verbosity >1) cout << "hseuil=" << hseuil << endl; 
  if(verbosity >1) cout << "NbTetrahedraRecollement " << nt_t << " / " << "NbVertex(anc)" << Th3.nt <<endl;  

  if(nt_t != Th3.nt){
    cout << " a tetrahedra was referenced twice or more " << endl;
    return 1;
  }
  /* determination du nombre de triangles confondus */ 
  hseuil = hmin/10.;
  hseuil = hseuil/3.;
  if(verbosity >1) cout << "TestSameTriangleMesh3 " << hseuil << endl;
  TestSameTriangleMesh3( Th3, hseuil, Psup, Pinf, Border_ok,  nbe_t );
 
  if(verbosity >1) cout << "hseuil=" << hseuil << endl; 
  if(verbosity >1) cout << "NbVertexRecollement " << nbe_t << " / " << "NbVertex(anc)" << Th3.nbe <<endl;  

  if(nbe_t != Th3.nbe){
    cout << " a triangle was referenced twice or more " << endl;
    return 1;
  }
    
  return 0;
}


Mesh3 *TestElementMesh3_patch( const Mesh3 & Th3 ) 
// Test si le maillage à des éléments communs : Sommet, triangle, ...
{
  R3 Pinf(1e100,1e100,1e100),Psup(-1e100,-1e100,-1e100);   // Extremité de la boîte englobante
  double hmin=1e10;   // longueur minimal des arrêtes
  double hseuil;
  int *Numero_Som = new int [Th3.nv];
  int nv_t,nt_t,nbe_t;
  
  // calcul de la boite englobante
  for (int ii=0;ii<Th3.nv;ii++){ 
    R3 P( Th3.vertices[ii].x, Th3.vertices[ii].y, Th3.vertices[ii].z);
    Pinf=Minc(P,Pinf);
    Psup=Maxc(P,Psup);     
  }

  // calcul de la longueur minimal des arrêtes
  for (int k=0;k<Th3.nt;k++){
    for (int e=0;e<6;e++){
      if( Th3[k].lenEdge(e) < Norme2(Psup-Pinf)/1e9 ) continue;
      hmin=min(hmin,Th3[k].lenEdge(e));   // calcul de .lenEdge pour un Mesh3
    }
  }
  
  for (int k=0;k<Th3.nbe;k++){
    for (int e=0;e<3;e++){
      if( Th3[k].lenEdge(e) < Norme2(Psup-Pinf)/1e9 ) continue;
      hmin=min(hmin,Th3.be(k).lenEdge(e));   // calcul de .lenEdge pour un Mesh3
    }
  }

  if(verbosity > 1) cout << "      - hmin =" <<  hmin << " ,  Bounding Box: " << Pinf << " "<< Psup << endl;

  ffassert(hmin>Norme2(Psup-Pinf)/1e9);

  // determination du nombre de sommets confondus
  hseuil = hmin/10.;     
  if(verbosity >1) cout << "TestSameVertexMesh3" << endl;
  TestSameVertexMesh3( Th3, hseuil, Psup, Pinf, nv_t, Numero_Som );
  
  if(verbosity >1) cout << "hseuil=" << hseuil << endl; 
  if(verbosity >1) cout << "NbVertexRecollement " << nv_t << " / " << "NbVertex(anc)" << Th3.nv <<endl;   

  /* degenerate element ??? */ 
  int *Elem_ok=new int[Th3.nt];
  int i_elem=0;
  for(int ii=0; ii< Th3.nt; ii++){
    const Tet & K(Th3.elements[ii]);
    int iv[4];
    
    Elem_ok[ii] = 1;
    
    for(int jj=0; jj <4; jj++){
      iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
    }
      
    for(int jj=0; jj<4; jj++){
      for(int kk=jj+1; kk<4; kk++){
	if( iv[jj]==iv[kk] ){
	  Elem_ok[ii] = 0;
	}
      }
    }
    i_elem = i_elem + Elem_ok[ii];
  }
   
  if( i_elem != Th3.nt ){
    cout << "There are a false tetrahedra in the mesh" << endl;
    //assert( i_elem == Th3.nt);
  }

  int *Border_ok=new int[Th3.nbe];
  int i_border= 0;
  for( int ii=0; ii< Th3.nbe; ii++){
    Border_ok[ii]=1;
    
    const Triangle3 & K(Th3.be(ii));
    int iv[3];
    
    for(int jj=0; jj <3; jj++){
      iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
      assert( iv[jj] >= 0 && iv[jj] < nv_t);
    }
    
    for(int jj=0; jj<3; jj++){
      for(int kk=jj+1; kk<3; kk++){
	if( iv[jj]==iv[kk] ) Border_ok[ii]=0;
      }
    }
    i_border = i_border + Border_ok[ii];
  }
  
  if( i_border != Th3.nbe){
    cout << "There are a false tetrahedra in the mesh" << endl;
    //assert( i_elem == Th3.nt);
  }

  /* determination du nombre de tetrahedre confondus */ 
  hseuil = hmin/10.;
  hseuil = hseuil/4.;
  nt_t=0;
  TestSameTetrahedraMesh3( Th3, hseuil, Psup, Pinf, Elem_ok, nt_t );

  if(verbosity >1) cout << "hseuil=" << hseuil << endl; 
  if(verbosity >1) cout << "NbVertexRecollement " << nt_t << " / " << "NbVertex(anc)" << Th3.nt <<endl;  

  /* determination du nombre de triangles confondus */ 
  hseuil = hmin/10.;
  hseuil = hseuil/3.;
  TestSameTriangleMesh3( Th3, hseuil, Psup, Pinf, Border_ok,  nbe_t );
 
  if(verbosity >1) cout << "hseuil=" << hseuil << endl; 
  if(verbosity >1) cout << "NbVertexRecollement " << nbe_t << " / " << "NbVertex(anc)" << Th3.nbe <<endl;  

  Vertex3   *v= new Vertex3[nv_t];
  Tet       *t= new Tet[nt_t];
  Triangle3 *b= new Triangle3[nbe_t]; 
  Tet       *tt=t;
  Triangle3 *bb=b; 

  EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,Pinf,Psup,0);

  // determination des nouveaux sommets
  int nbv = 0; 
  hseuil = hmin/10.;
  for (int ii=0;ii<Th3.nv;ii++){
    const Vertex3 &vi(Th3.vertices[ii]); 
    Vertex3 * pvi=gtree->ToClose(vi,hseuil);
      
    if(!pvi){
      v[nbv].x = vi.x;
      v[nbv].y = vi.y;
      v[nbv].z = vi.z;
      v[nbv].lab = vi.lab;
      gtree->Add( v[nbv++] );
    }
   
  }
  delete gtree;
  assert(nbv == nv_t);

  // determination des nouveaux tetrahedres
  int nbt = 0; 
  hseuil = hmin/10.;
  hseuil = hseuil/4.;
  for (int ii=0;ii<Th3.nt;ii++){
    if(Elem_ok[ii] == 0) continue; 
    const Tet  &K(Th3.elements[ii]);
    int iv[4];
    iv[0]=Numero_Som[Th3.operator()(K[0])];
    iv[1]=Numero_Som[Th3.operator()(K[1])];
    iv[2]=Numero_Som[Th3.operator()(K[2])];
    iv[3]=Numero_Som[Th3.operator()(K[3])];
    (tt++)->set(v,iv,K.lab);
 
  }
  assert(nbv == nv_t);

  // determination des nouveaux trianglesxs
  int nbbe = 0; 
  hseuil = hmin/10.;
  hseuil = hseuil/4.;
  for (int ii=0;ii<Th3.nbe;ii++){
    if(Border_ok[ii] == 0) continue; 
    const Triangle3 & K(Th3.be(ii));
    int iv[3];
    iv[0]=Numero_Som[Th3.operator()(K[0])]; 
    iv[1]=Numero_Som[Th3.operator()(K[1])]; 
    iv[2]=Numero_Som[Th3.operator()(K[2])]; 
  
    (bb++)->set(v,iv,K.lab);
  }
  assert(nbbe == nbe_t);

  delete [] Numero_Som;
  delete [] Border_ok;
  delete [] Elem_ok;
  
  Mesh3 *Th3_new = new Mesh3(nv_t,nt_t,nbe_t,v,t,b); 
  return Th3_new;
}





// TransfoMesh_v2.cpp


// LayerMesh.cpp

// remarque choix 2 est a encore a determiner
double  zmin_func_mesh( const int choix, const double x, const double y  )
{
  
  switch(choix){
  case 0: 
    return 0.;
    break;
  case 1:
    return 0.;   
    break;
  case 2:
    return sqrt(pow(x,2)+pow(y,2));
    break;
  default :
    cout << "zmin_func no defined" << endl;
    return 0.;
  }
}

double  zmax_func_mesh( const int choix, const double x, const double y  ){
  
  switch(choix){
  case 0: 
    return 1.;
    break;
  case 1:
    return 1.;   
    break;
  case 2:
    return 3.+sqrt(pow(x,2)+pow(y,2));
    break;
  default :
    cout << "zmaxfunc no defined" << endl;
    return 0.;
  }
}

int Ni_func_mesh( const int choix, const double x, const double y  ){
  const int multi=1;
  int res;
  switch(choix){
  case 0:
    if( x==0. && y==0.){
      res = 3;
    }
    if( x==1. && y==0.){
      res = 5;
    }
    if( x==0. && y==1.){
      res = 7;
    }
    if( x==0.5 && y==0.5){
      res = 6;
    }
    return res ;
    //return multi;
    break;
  case 1:
    return 2;   
    break;
  case 2:
    return int(multi*(3+sqrt(pow(x,2)+pow(y,2))));
    break;
  default :
    cout << "Ni_func no defined" << endl;
    return 0;
  }
}

void discretisation_max_mesh(const int choix,  const Mesh & Th2, int & Nmax){  
  int Ni;

  Nmax = 0;  

  /*for(int ii=0; ii < A2D.NbSommet2D;ii++){
    Ni   = Ni_func( choix, A2D.CoorSommet2D[ii][0], A2D.CoorSommet2D[ii][1]); 
  	Nmax = max(Ni,Nmax);
  }
  Nmax=4;*/
  for(int ii=0; ii < Th2.nv; ii++){
     const  Mesh::Vertex & P = Th2.vertices[ii];
     Ni   = Ni_func_mesh( choix, P.x, P.y); 
  	 Nmax = max(Ni,Nmax);
  }
}

void tab_zmin_zmax_Ni_mesh(const int choix, const Mesh & Th2, int & Nmax,double *tab_zmin, double *tab_zmax,int *tab_Ni){
   Nmax = 0;	
   for(int ii=0; ii < Th2.nv; ii++){
	 const  Mesh::Vertex & P = Th2.vertices[ii];
     tab_Ni[ii] = Ni_func_mesh( choix, P.x, P.y);
	 tab_zmin[ii] = zmin_func_mesh( choix, P.x, P.y );
     tab_zmax[ii] = zmax_func_mesh( choix, P.x, P.y );    
     Nmax = max(tab_Ni[ii],Nmax);
   }
}
/* Fonction permettant de transformer maillage 2D en maillage 3D*/

void Tet_mesh3_mes_neg(Mesh3 & Th3){
    int iv[4];
    int lab;
    
    for(int ii=0; ii< Th3.nt; ii++){
      const Tet & K(Th3.t(ii));
      lab   = K.lab;
 
      iv[0] = Th3.operator()(K[0]);
      iv[2] = Th3.operator()(K[1]);
      iv[1] = Th3.operator()(K[2]);
      iv[3] = Th3.operator()(K[3]);
      R3 A(Th3.vertices[iv[0]]);
      R3 B(Th3.vertices[iv[1]]);
      R3 C(Th3.vertices[iv[2]]);
      R3 D(Th3.vertices[iv[3]]);
      double mes=det(A,B,C,D)/6.;
      Th3.t(ii).set(Th3.vertices, iv, lab,mes);	
    }
}

//=======================================================================//
//   Rajout pour s'assurer un unique label pour les vertices
//=======================================================================//

void build_layer_map_tetrahedra(const Mesh &Th2, map<int, int> &maptet ){
  
  int numero_label=0;
  //cout << "in: buil_layer_map_tetrahedra" << endl;
  for(int ii=0; ii< Th2.nt; ii++){
    //cout << "ii= " << ii  << "Th2.nt=" << Th2.nt <<endl; 
    const  Mesh::Triangle & K(Th2.t(ii));
    map<int,int>::const_iterator imap=maptet.find(K.lab);
    //cout << "K.lab= " << K.lab << endl; 
    if(imap == maptet.end()){
	maptet[ K.lab ] = K.lab; //  modif FH .. numero_label;
      numero_label = numero_label+1;
    }
  }
  //cout << "number of tetraedra label=" << numero_label << endl;
}
void build_layer_map_triangle(const Mesh &Th2, map<int, int> &maptrimil, map<int, int> &maptrizmax, map<int, int> &maptrizmin ){
  
  int numero_label=0;
  //cout << "in: buil_layer_map_triangle" << endl;
  for(int ii=0; ii< Th2.nt; ii++){
    const Mesh::Triangle & K(Th2.t(ii));
    map<int,int>::const_iterator imap=maptrizmax.find(K.lab);
    
    if(imap == maptrizmax.end()){
	maptrizmax[ K.lab ] = K.lab;// modif FH jan 2010  numero_label;
      numero_label = numero_label+1;
    }
  }
  
  for(int ii=0; ii< Th2.nt; ii++){
    const Mesh::Triangle & K(Th2.t(ii));
    map<int,int>::const_iterator imap=maptrizmin.find(K.lab);
    
    if(imap == maptrizmin.end()){
	maptrizmin[ K.lab ] =  K.lab;// modif FH jan 2010 numero_label;
      numero_label = numero_label+1;
    }
  }
  
  for(int ii=0; ii< Th2.neb; ii++){
    const Mesh::BorderElement & K(Th2.be(ii));
    map<int,int>::const_iterator imap=maptrimil.find(K.lab);
    
    if(imap == maptrimil.end()){
	maptrimil[ K.lab ] = K.lab  ;//modif FH jan 2010 numero_label;
      numero_label = numero_label+1;
    }
  }
  
}
	
void build_layer_map_edge(const Mesh &Th2, map<int, int> &mapemil, map<int, int> &mapezmax, map<int, int> &mapezmin ){
  
  int numero_label=0;
		
  for(int ii=0; ii< Th2.neb; ii++){
    const Mesh::BorderElement & K(Th2.be(ii));
    map<int,int>::const_iterator imap1=mapezmax.find(K.lab);
    map<int,int>::const_iterator imap2=mapemil.find(K.lab);
    map<int,int>::const_iterator imap3=mapezmin.find(K.lab);
    
    if(imap1 == mapezmax.end()){
	mapezmax[ K.lab ] = K.lab ;//modif FH jan 2010 numero_label;
      numero_label = numero_label+1;
    }
			
    if(imap2 == mapemil.end()){
      mapemil[ K.lab ] = K.lab ;//modif FH jan 2010 numero_label;numero_label;
      numero_label = numero_label+1;
    }
    
    if(imap3 == mapezmin.end()){
      mapezmin[ K.lab ] = K.lab ;//modif FH jan 2010 numero_label;numero_label;
      numero_label = numero_label+1;
    }
    
  }
 	
}

	
Mesh3 * build_layer (const Mesh & Th2, const int Nmax, const int *tab_Ni, 
		const double *tab_zmin, const double *tab_zmax, 
		const map<int, int> &maptet, const map<int, int> &maptrimil, const map<int, int> &maptrizmax, const map<int, int> &maptrizmin, 
		const map<int, int> &mapemil, const map<int, int> &mapezmax, const map<int, int> &mapezmin ){

  int MajSom, MajElem, MajBord2D;     
  Mesh3 *Th3=new Mesh3;
  NbSom3D_NbElem3D_NbBord2D_mesh_product_mesh_tab( Nmax, tab_Ni, Th2, MajSom, MajElem, MajBord2D);   
  if(verbosity > 1) cout << "MajSom = " <<  MajSom  << "  "  << "MajElem = " <<  MajElem  << " " << "MajBord2D =" << MajBord2D << endl;
  
  if(verbosity > 1) cout << "debut :   Th3.set(MajSom, MajElem, MajBord2D);     "<< endl;
  Th3->set(MajSom,MajElem,MajBord2D);
  
  if(verbosity > 1) cout << "debut :   Som3D_mesh_product_Version_Sommet_mesh_tab( Nmax, tab_Ni, tab_zmin, tab_zmax, Th2, Th3);   "<< endl;
  Som3D_mesh_product_Version_Sommet_mesh_tab( Nmax, tab_Ni, tab_zmin, tab_zmax, Th2, 
			maptet, maptrimil, maptrizmax, maptrizmin, mapemil, mapezmax, mapezmin, *Th3);    
  
  
  //  Add FH because remove in call function.. 
  
  Th3->BuildBound();
  Th3->BuildAdj();
  Th3->Buildbnormalv();  
  Th3->BuildjElementConteningVertex();
  
    
  return Th3;
}

void NbSom3D_NbElem3D_NbBord2D_mesh_product_mesh_tab(const int Nmax, const int *tab_Ni, const Mesh &Th2,  int &MajSom, int &MajElem, int &MajBord2D){  
  int i;
  
  MajSom = 0;  
  for(int ii=0; ii < Th2.nv;ii++){  
    MajSom = MajSom + (tab_Ni[ii]+1);        
  	assert(tab_Ni[ii]<=Nmax);
  }

  MajElem = 0;  
  for(int ii=0; ii < Th2.nt; ii++){   
    const Mesh::Triangle & K(Th2.t(ii));
    for(int jj=0; jj < 3; jj++){ 
      //i  = A2D.ElemPoint2D[ii][jj]; 
      i  = Th2.operator()(K[jj]); 
      MajElem = MajElem + tab_Ni[i]; 
    }
  }
  
  // determination of NbBord2D
  MajBord2D = 2*Th2.nt;  
 
  for(int ii=0; ii < Th2.neb;ii++)
    {
      const Mesh::BorderElement  & K(Th2.be(ii));
      for(int jj=0; jj < 2; jj++)
	{  
	  // i  = A2D.ElemBord1D[ii][jj];
	  i=Th2.operator()(K[jj]); 
	  
	  MajBord2D = MajBord2D + tab_Ni[i];
	  assert( tab_Ni[i] <= Nmax);	
	}
    }
  //exit(1);
}

void Som3D_mesh_product_Version_Sommet_mesh_tab(const int Nmax, 
			const int *tab_Ni, const double *tab_zmin, const double *tab_zmax, const Mesh &Th2, 
			const map<int, int> &maptet, const map<int, int> &maptrimil, const map<int, int> &maptrizmax, const map<int, int> &maptrizmin, 
			const map<int, int> &mapemil, const map<int, int> &mapezmax, const map<int, int> &mapezmin, 
			Mesh3 & Th3){
  // intent(in)  Nmax,Mesh &A2D
  // intent(out) Mesh3 &A3D
  
  double val_zmin,val_zmax,val_dz;
  int    Ni;
  int    NumSommet;
  int    NumElement;
  KN<int>  tab_NumSommet(Th2.nv+1);

  // variable tet
  int    SommetPrisme[6];
 
  // variable creer pour le bord
  int    i_ind1,Ni_ind1;
  int    i_ind2,Ni_ind2; 
  int    i_recoll_1pp,i_recoll_2pp; 
  int    i_recoll_1, i_recoll_2; 
  //int    pas_recoll_1, pas_recoll_2; 
  int    type_dec_border; 
  
  // avec data
  int i_recoll_jMax,i_recoll_jMaxpp;
  int cas_decoupage; //, cas_data; 
  int int_decoup[3] = {1,2,4};
  int Ni_elem[3]; 
  int DiagMax1,DiagMax2;

  // determination of maximum label for vertices 
  NumSommet = 0;

  for( int ii=0; ii < Th2.nv; ii++){
    const  Mesh::Vertex & P = Th2.vertices[ii];

    val_zmin = tab_zmin[ii];
    val_zmax = tab_zmax[ii];    
    Ni       = tab_Ni[ii];

   
    //val_dz = (val_zmax - val_zmin)/Ni;
    if( Ni == 0){
      val_dz = 0.;
    }
    else{
      val_dz = (val_zmax - val_zmin)/Ni;
      //if( abs(val_dz) < 1e-9 ) Ni=0; 
    }
    

    tab_NumSommet[ii] = NumSommet; // Numero du premier sommet 3D associé au sommet 2D ii.
    //cout << "ii, tab_NumSommet[ii]= "<< ii <<" "<< tab_NumSommet[ii] << endl;
    
    for(int j=0; j <= Ni; j++){ //changer
      Th3.vertices[NumSommet].x = P.x; 
      Th3.vertices[NumSommet].y = P.y; 
      Th3.vertices[NumSommet].z = val_zmin + val_dz*j;

      Th3.vertices[NumSommet].lab = P.lab; 
      // cas maillage du bas et du haut, on un nouveau label
      if(j==0)  Th3.vertices[NumSommet].lab = P.lab ;
      if(j==Ni) Th3.vertices[NumSommet].lab = P.lab ;      
      NumSommet = NumSommet+1;
    }
  
  }
  tab_NumSommet[Th2.nv] = NumSommet;

  assert( NumSommet == Th3.nv );
  
  /*********************************/
  /*  new label for edges of cube  */
  /*********************************/ 
  /*
  cout << " new label for edges of cubes  " << endl;
  for(int ii=0; ii < Th2.neb; ii++){ 
	const Mesh::BorderElement & K(Th2.be(ii));
	int ib[2];
	ib[0] = Th2.operator()(K[0]);
    ib[1] = Th2.operator()(K[1]);  
    //map<int,int>:: const_iterator imap;
    
    for(int kk=0; kk<2;kk++){
	    // label zmin
	    map<int,int>:: const_iterator imap1;
	    imap1=mapezmin.find( K.lab );

	    assert( imap1!=mapezmin.end() );
	    Th3.vertices[ tab_NumSommet[ib[kk]] ].lab = imap1->second;  
	       
	    // label zmax
	    map<int,int>:: const_iterator imap2;
	    imap2=mapezmax.find( K.lab );
	    assert( imap2!=mapezmax.end() );
	       
	    Th3.vertices[ tab_NumSommet[ ib[kk] ] + tab_Ni[ib[kk]] ].lab = imap2->second; 
	       
	       // label côté
	     map<int,int>:: const_iterator imap3;
	     imap3=mapemil.find ( K.lab );
	     assert( imap3!=mapemil.end() );
	       
	     for(int jj=1; jj < tab_Ni[ib[kk]]; jj++){
		       Th3.vertices[ tab_NumSommet[ ib[kk] ] + jj ].lab = imap3->second;
		 }
		 }
  }
  */
  
  
  //=======================================================================
  // creation des bord du maillage 3D a partir du bord 1D et du maillage 2D
  //=======================================================================

  if(verbosity > 1) cout << "calcul element du bord " << endl;
 
  // A mettre plus haut
  int ElemBord;

  ElemBord = 0;

  // bord définies en zmax
  
  for(int ii=0; ii < Th2.nt; ii++){
    int ijj[3];
    const Mesh::Element & K(Th2.t(ii));
    int lab;
    map<int,int>::const_iterator imap=maptrizmax.find(K.lab);
    assert( imap!=maptrizmax.end() );
    lab=imap->second;
   
    for(int kk=0; kk < 3; kk++){
      ijj[kk] = Th2.operator()(K[kk]);
      ijj[kk] = tab_NumSommet[ijj[kk]+1]-1;
    }
     
    Th3.be(ElemBord).set(Th3.vertices,ijj,lab);
    
    ElemBord = ElemBord+1;
  }

  //cout << "bord en zmin" << endl;
 
  for(int ii=0; ii < Th2.nt; ii++){
    int ijj[3];//bjj[3];
    const Mesh::Element & K(Th2.t(ii));
    int lab; 
    map<int,int>::const_iterator imap=maptrizmin.find(K.lab);
    assert( imap!=maptrizmin.end() );
    lab = imap->second; 
    
  
    for(int kk=0; kk < 3; kk++){
      ijj[2-kk] = Th2.operator()(K[kk]);
      //bjj[2-kk] = ijj[2-kk] ;
      ijj[2-kk] = tab_NumSommet[ijj[2-kk]];
    }

    Th3.be(ElemBord).set(Th3.vertices,ijj,lab);

    ElemBord = ElemBord+1;
  }
  
  //cout << "bord sur le cote" << endl;

  for(int ii=0; ii < Th2.neb; ii++){ // Th2.neb ??
    int ijj[3];
	
    const Mesh::BorderElement & K(Th2.be(ii));
    int lab;
    
    map<int,int>::const_iterator imap=maptrimil.find(K.lab);
    assert( imap!=maptrimil.end() );
    lab=imap->second;

    int edgebid ;
    int ffbid   = Th2.BoundaryElement( ii, edgebid );     // ii : number of edge => sortie :: ffbid = numero triangles, edgebid = numero edges
    int j0bid,j1bid;
    Th2.VerticesNumberOfEdge( Th2.t(ffbid), edgebid, j0bid, j1bid);

    //bool ffsens = Th2.SensOfEdge( Th2.t(ffbid), edgebid ); // sens du parcours de la edge correcte ou non

    /*
      if( ffsens == true){
      i_ind1  = Th2.operator()(K[0]);
      i_ind2  = Th2.operator()(K[1]);
      }
      else{
      i_ind1  = Th2.operator()(K[1]);
      i_ind2  = Th2.operator()(K[0]);
      }
      

      printf("value of vertex edge (verticesNumberOfEdge) :: %d--%d \n", j0bid, j1bid );
      printf("value of vertex edge ( Th2.operator() ) :: %d--%d \n",  Th2.operator()(K[0]), Th2.operator()(K[1]) );
      printf("value of vertex edge ( bool sens  ) :: %d--%d \n",  i_ind1, i_ind2 );
    */
    i_ind1 = j0bid;
    i_ind2 = j1bid;


    Ni_ind1 =  tab_Ni[i_ind1]; 
    Ni_ind2 =  tab_Ni[i_ind2]; 
	  
    assert( Ni_ind1 <= Nmax);
    assert( Ni_ind2 <= Nmax);
	
    for(int jNmax=Nmax-1; jNmax >=0; jNmax--){
	    
      /*
      i_recoll_1pp = int((jNmax+1)*Ni_ind1/Nmax);
      i_recoll_2pp = int((jNmax+1)*Ni_ind2/Nmax); 
      
      i_recoll_1 = int(jNmax*Ni_ind1/Nmax);
      i_recoll_2 = int(jNmax*Ni_ind2/Nmax); 
      */
      
      i_recoll_1 = int((jNmax+1)*Ni_ind1/Nmax);
      i_recoll_2 = int((jNmax+1)*Ni_ind2/Nmax); 
      
      i_recoll_1pp = int(jNmax*Ni_ind1/Nmax);
      i_recoll_2pp = int(jNmax*Ni_ind2/Nmax);
      
//       if( (i_ind1== 11 ||  i_ind1== 0) && (i_ind2==11 || i_ind2==0) ) {
// 	printf("i_recoll1   %d,    i_recoll2 %d\n", i_recoll_1, i_recoll_2);
// 	printf("i_recoll1pp %d,  i_recoll2pp %d\n", i_recoll_1pp, i_recoll_2pp);
//       }
      /*

	1     ===   2 
	
	|           |
	|           |
	
	1pp   ===   2pp   
	
	sens 2D : 1pp => 2pp et 1 => 2

	type_dec_border = 0  tous les points sont confondus
	type_dec_border = 1  les points 1pp et 1 sont differents
	type_dec_border = 2  les points 2pp et 2 sont differents
	type_dec_border = 3  les points 1pp et 1 et les points 2pp et 2 sont differents 

	rappel :  1pp(0) 2pp(1) 2(2) 1(3) 
	data_dec_border 1 :   {3 1 0}
	data_dec_border 2 :   {2 1 0}
	data_dec_border 3 : type1 : { {2 1 0}{0 3 2} }
	                  : type2 : { {3 1 0}{1 3 2} }
      */
      type_dec_border = 0; 
      
      if( i_recoll_1pp != i_recoll_1){ 
	type_dec_border = type_dec_border + 1;
      }
      
      if( i_recoll_2pp != i_recoll_2){ 
	type_dec_border = type_dec_border + 2;
      }

//   if( (i_ind1== 11 ||  i_ind1== 0) && (i_ind2==11 || i_ind2==0) ) 
// 	cout << "type decoupage bord= " <<  type_dec_border <<endl;
           
      switch( type_dec_border ){
      case 0:
	// rien n a faire
	break;
      case 1:
	// 2pp = 2
	// avant 2,1,0 --> 0,1,2
	ijj[0] = tab_NumSommet[i_ind1]+i_recoll_1pp;   
	ijj[1] = tab_NumSommet[i_ind2]+i_recoll_2pp;
	ijj[2] = tab_NumSommet[i_ind1]+i_recoll_1; 
	
	Th3.be(ElemBord).set(Th3.vertices,ijj,lab);
	
	ElemBord = ElemBord+1;
	break;
      case 2:
	// 1pp = 1
	// avant 2,1,0 --> 0,1,2
	ijj[0] = tab_NumSommet[i_ind1]+i_recoll_1pp;
	ijj[1] = tab_NumSommet[i_ind2]+i_recoll_2pp;
	ijj[2] = tab_NumSommet[i_ind2]+i_recoll_2; 
	
	Th3.be(ElemBord).set(Th3.vertices,ijj,lab);
	
	ElemBord = ElemBord+1;
	break;
      case 3:
	int idl;
	// determination de la diagonale Max
	DiagMax1 = max( tab_NumSommet[i_ind1]+i_recoll_1pp, tab_NumSommet[i_ind2]+i_recoll_2 );
	DiagMax2 = max( tab_NumSommet[i_ind2]+i_recoll_2pp, tab_NumSommet[i_ind1]+i_recoll_1 );	

	
	if(DiagMax1 > DiagMax2){  
	  idl = 1; 
      
	  ijj[0] = tab_NumSommet[i_ind1]+i_recoll_1pp;
	  ijj[1] = tab_NumSommet[i_ind2]+i_recoll_2pp;
	  ijj[2] = tab_NumSommet[i_ind2]+i_recoll_2; 
	
	  Th3.be(ElemBord).set(Th3.vertices,ijj,lab);

	  ijj[0] = tab_NumSommet[i_ind2]+i_recoll_2;
	  ijj[1] = tab_NumSommet[i_ind1]+i_recoll_1;
	  ijj[2] = tab_NumSommet[i_ind1]+i_recoll_1pp; 
	
	  Th3.be(ElemBord+1).set(Th3.vertices,ijj,lab);

	}
	else{
	  idl = 2;
	  
	  ijj[0] = tab_NumSommet[i_ind1]+i_recoll_1pp;
	  ijj[1] = tab_NumSommet[i_ind2]+i_recoll_2pp;
	  ijj[2] = tab_NumSommet[i_ind1]+i_recoll_1; 

	  Th3.be(ElemBord).set(Th3.vertices,ijj,lab);
	  
	  ijj[0] = tab_NumSommet[i_ind2]+i_recoll_2;
	  ijj[1] = tab_NumSommet[i_ind1]+i_recoll_1;
	  ijj[2] = tab_NumSommet[i_ind2]+i_recoll_2pp; 

	  Th3.be(ElemBord+1).set(Th3.vertices,ijj,lab);
	}
	//cout << "idl=" << idl << endl; 
	ElemBord = ElemBord+2;
	break;
      default:
	break;
      }     
    }    
  }
  
  assert( ElemBord == Th3.nbe );
  //=========================================
  // Creation + determination tetraedre

  if(verbosity > 1) cout << "calcul element tetraedre " << endl;
  
  NumElement =  0;

  for(int ii=0; ii < Th2.nt; ii++){
      /*  
	  nouvelle numerotation : 
	  -----------------------
	  Valeur de cas_deoupage
	  -----------------------
	  1 : sommet 0 et 3 differents
	  2 : sommet 1 et 4 differents
	  4 : sommet 2 et 5 differents
	  ============================
	  3 : sommet 0 et 3 differents + sommet 1 et 4 differents
	  5 : sommet 0 et 3 differents + sommet 2 et 5 differents
	  6 : sommet 1 et 4 differents + sommet 2 et 5 differents
	  ============================
	  7 : aucun sommets confondus

	  data_tetraedre 
	  ==============
	  1: 0++,1++,2++,SomDiff :  {0 1 2 3}        :: data 1
	  2: 0++,1++,2++,SomDiff :  {0 1 2 4}        :: data 2
	  4: 0++,1++,2++,SomDiff :  {0 1 2 5}        :: data 3
	  ==============
	  = deux cas possible depend du sommet le plus grand : Sommet le plus grand est un ++
	  = 0++,1++,2++,SomDiffMin  ||  SomDiffMax++, j_SomDiff_py[j_SomEgal][0], j_SomDiff_py[j_SomEgal][1], Som_Egal 
	  3:a: SommetMax diag 04  {0,1,2,4} {5,4,3,0} :: data 4
	  3:b: SommetMax diag 13  {0,1,2,3} {5,4,3,1} :: data 5
	  =============================================
	  5:a: SommetMax diag 05  {0,1,2,5} {5,4,3,0} :: data 6
	  5:b: SommetMax diag 23  {0,1,2,3} {5,4,3,2} :: data 7
	  =============================================
	  6:a: SommetMax diag 15  {0,1,2,5} {5,4,3,1} :: data 8
	  6:b: SommetMax diag 24  {0,1,2,4} {5,4,3,2} :: data 9
	  =============================================
	  7: aller chercher dans la fonction          :: data 10 a data 
	  == voir hecht routine

      */
      const Mesh::Element & K(Th2.t(ii));
      int somv[4];
      int K_jj[3];
      int lab;
      
      map<int,int>::const_iterator imap=maptet.find(K.lab); 
      assert( imap != maptet.end() );
      lab=imap->second;
     
      // valeur de Nombre de points
      for(int jj=0; jj <3; jj++){
        K_jj[jj] = Th2.operator()(K[jj]);   
        Ni_elem[jj] = tab_Ni[ K_jj[jj] ];    
      }
     
      for(int jNmax=Nmax-1; jNmax >=0; jNmax--){
	// determination des sommets + cas decoupage
	cas_decoupage = 0;
	for(int jj=0; jj<3; jj++){
	  
	  i_recoll_jMax   = int( (jNmax)*Ni_elem[jj]/Nmax );
	  i_recoll_jMaxpp = int( (jNmax+1)*Ni_elem[jj]/Nmax );
	  
	  SommetPrisme[jj+3] = tab_NumSommet[ K_jj[jj] ] + i_recoll_jMaxpp;
	  SommetPrisme[jj] = tab_NumSommet[ K_jj[jj] ] + i_recoll_jMax;
	  
	  assert( SommetPrisme[jj]   <= Th3.nv);
	  assert( SommetPrisme[jj+3] <= Th3.nv);
	  if( i_recoll_jMax != i_recoll_jMaxpp) cas_decoupage = cas_decoupage + int_decoup[jj];
	}
	
	//cout << "cas du decoupage= " << cas_decoupage << endl;
	
	switch( cas_decoupage ){
	  
	case 0 : 
	// les points sont tous confondus pas d ajout element : rien a faire
	  break;
	  /*
	    CAS CREATION D UN TETRAEDRE : cas decoupage 1 2 4
	    
	  */
	case 1 :
	  // On a un tetraedre
	
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1]; 
	  somv[2] = SommetPrisme[2]; 
	  somv[3] = SommetPrisme[3];
	  
	  Th3.elements[NumElement].set(Th3.vertices, somv, lab);	
	  
	  NumElement = NumElement+1;
	  break;
	case 2 :
	  // On a un tetraedre
	  
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1]; 
	  somv[2] = SommetPrisme[2]; 
	  somv[3] = SommetPrisme[4];

	  Th3.elements[NumElement].set(Th3.vertices, somv, lab);	
	  
	  NumElement = NumElement+1;
	  break;
	case 4 :
	  // On a un tetraedre
	  
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1]; 
	  somv[2] = SommetPrisme[2]; 
	  somv[3] = SommetPrisme[5];
	  
	  Th3.elements[NumElement].set(Th3.vertices, somv, lab);	
	
	  NumElement = NumElement+1;
	  break;
	/*
	  On a une pyramide a base rectangle: decoupe deux tetraedres
	  cas decoupage 3 5 6
	*/
      case 3 :
	// determination de la diagonale dominante
	DiagMax1 = max( SommetPrisme[0], SommetPrisme[4] );
	DiagMax2 = max( SommetPrisme[1], SommetPrisme[3] );
	
	//cout << "DiagMax1=" << DiagMax1 << " "<< SommetPrisme[0]<<" " <<SommetPrisme[4] << endl;

	if( DiagMax1 > DiagMax2){
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[4]; 
	  
	  Th3.elements[NumElement].set(Th3.vertices, somv, lab);	
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[0]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}			  
	else{
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[3];
	  
	  Th3.elements[NumElement].set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[1]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}
      
	NumElement = NumElement+2;
	break;

      case 5 :
	// determination de la diagonale dominante
	DiagMax1 = max( SommetPrisme[0], SommetPrisme[5] );
	DiagMax2 = max( SommetPrisme[2], SommetPrisme[3] );
	
	//cout << "DiagMax1=" << DiagMax1 << " "<< SommetPrisme[0]<<" " <<SommetPrisme[5] << endl;

	if( DiagMax1 > DiagMax2){
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[5];
	  
	  Th3.elements[NumElement].set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[0]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}			  
	else{
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[3];
	  
	  Th3.elements[NumElement].set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[2]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}
      
	NumElement = NumElement+2;
	break;

      case 6 :
	// determination de la diagonale dominante
	DiagMax1 = max( SommetPrisme[1], SommetPrisme[5] );
	DiagMax2 = max( SommetPrisme[2], SommetPrisme[4] );

	//cout << "DiagMax1=" << DiagMax1 << " "<< SommetPrisme[1]<<" " <<SommetPrisme[5] << endl;
	
	if( DiagMax1 > DiagMax2){
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[5];
	  
	  Th3.elements[NumElement].set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[1]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}			  
	else{
	  //------------------
	  // premier tetraedre 
	  somv[0] = SommetPrisme[0];
	  somv[1] = SommetPrisme[1];	
	  somv[2] = SommetPrisme[2];
	  somv[3] = SommetPrisme[4];
	  
	  Th3.elements[NumElement].set(Th3.vertices, somv, lab);	 
	  // deuxieme tetraedre
	  somv[0] = SommetPrisme[5];
	  somv[1] = SommetPrisme[4];	
	  somv[2] = SommetPrisme[3];
	  somv[3] = SommetPrisme[2]; 
	  
	  Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	
	}
      
	NumElement = NumElement+2;
	break;

      case 7 :
	// on a un prisme 
	int nbe;
	int option=1;
	int idl[3];
	int nu[12];
	
	DiagMax1 = max( SommetPrisme[0], SommetPrisme[5] );
	DiagMax2 = max( SommetPrisme[2], SommetPrisme[3] );	
	
	// determination de idl
	// idl[0] : choix sommet 0 ou 2 (dpent1 equivalent 1 ou 3)
	
	if(DiagMax1 > DiagMax2){  
	  idl[0]=1;
	}
	else{
	  idl[0]=2;
	}
    
	DiagMax1 = max( SommetPrisme[0], SommetPrisme[4] );
	DiagMax2 = max( SommetPrisme[1], SommetPrisme[3] );	
	
	// idl[1] : choix sommet 0 ou 1 (dpent1 equivalent 1 ou 2)	
    if(DiagMax1 > DiagMax2){
	  idl[1]=1;
	}
	else{
	  idl[1]=2;
	}
	
	DiagMax1 = max( SommetPrisme[1], SommetPrisme[5] );
	DiagMax2 = max( SommetPrisme[2], SommetPrisme[4] );	
	
	// idl[2] : choix sommet 1 ou 2 (dpent1 equivalent 2 ou 3)	
	if(DiagMax1 > DiagMax2){
	  idl[2]=1;
	}
	else{
	  idl[2]=2;
	}
	//cout << "idl[0] << << idl[1] << << idl[2]" << endl;
	//cout << idl[0] << " " << idl[1] << "  "<< idl[2] << endl;
	
	nbe = 0;

	dpent1_mesh( idl, nu, nbe, option);
      
	if(nbe!=3){cout << nbe << endl; cerr << "probleme dans dpent1_mesh" << endl;  };
	
	//------------------
	// premier tetraedre 
	somv[0] = SommetPrisme[nu[0]];
	somv[1] = SommetPrisme[nu[1]];	
	somv[2] = SommetPrisme[nu[2]];	
	somv[3] = SommetPrisme[nu[3]];
	
	Th3.elements[NumElement].set(Th3.vertices, somv, lab);		
	// deuxieme tetraedre
	somv[0] = SommetPrisme[nu[4]]; 
	somv[1] = SommetPrisme[nu[5]]; 
	somv[2] = SommetPrisme[nu[6]]; 
	somv[3] = SommetPrisme[nu[7]];
	
	Th3.elements[NumElement+1].set(Th3.vertices, somv, lab);	 	
	// troisieme tetraedre
	somv[0] = SommetPrisme[nu[8]]; 
	somv[1] = SommetPrisme[nu[9]]; 
	somv[2] = SommetPrisme[nu[10]]; 
	somv[3] = SommetPrisme[nu[11]]; 
	
	Th3.elements[NumElement+2].set(Th3.vertices, somv, lab);	
	
	NumElement = NumElement+3;
	break;
      }
      
   }
  // Au final : les sommers des tetraedres et la conectivité des tetraedres finaux
  assert(NumElement <= Th3.nt);
  }

}

void dpent1_mesh(int idl[3],int nu[12],int &nbe,int &option){
  // intent(inout)  :: idl
  // intent(out)    :: nu,nbe,option
  // option ne sert à rien
  //* version simplifie pour le mailleur par couche 2D 3D
  //-----------------------------------------------------------------------
  //      subroutine dpent1 (idl,nu,nbe,option)
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //                s.p. dpent1
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  but : decoupe un pentaedre en 3 tetreadres suivant la decoupe des 3
  // ---   faces frontieres a 4 cotes
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  parametres en entre :
  //   idl : parametre de decoupe de face calculer comme ceci :
  //       si idl(i) = 0 alors la face n'est pas  decoupee
  //       idl(1)= 1 si la face 1463 est decoupe par l'arete 16 ,sinon 2
  //       idl(2)= 1 si la face 1254 est decoupe par l'arete 15 ,sinon 2
  //       idl(3)= 1 si la face 2365 est decoupe par l'arete 26 ,sinon 2
  //          id = i1 + i2 * 2 + i3 * 4
  //  parametres en sortie :
  //   nbe : nbe de tetraedre de la decoupe
  //         nbe = 0 => decoup impossible
  //         nbe = 3 => decoup possible le tableau nu est genere
  //   nu(1:4,1:nbe) : tableau de numero des sommet 3 tetraedres dans le
  //                   pentaedre
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  //  programmation : f77 ->c++ subroutine de f. hecht upmc 

  int idp[8];
  int i1,i2,i3,i,nbdp,idf,idecou;
  const int pdd[8] = {1,0,2,3,4,5,0,6};

  int mu[6][12];
  const int mu0[12] = {1,6,2,3, 1,5,2,6, 1,6,4,5};
  const int mu1[12] = {1,6,2,3, 1,4,2,6, 2,6,4,5};
  const int mu2[12] = {1,4,2,3, 2,6,3,4, 2,6,4,5};
  const int mu3[12] = {1,5,2,3, 1,5,3,6, 1,6,4,5};
  const int mu4[12] = {1,5,2,3, 1,5,3,4, 3,6,4,5};
  const int mu5[12] = {1,4,2,3, 2,5,3,4, 3,6,4,5};

  for(int jj=0; jj<12; jj++){
    mu[0][jj] = mu0[jj];
    mu[1][jj] = mu1[jj];
    mu[2][jj] = mu2[jj];
    mu[3][jj] = mu3[jj];
    mu[4][jj] = mu4[jj];
    mu[5][jj] = mu5[jj];
  }

  // calcul des descoupes possible du pentaedre
  idf  = -1;
  nbdp =  0;

  for(i3=1; i3<=2; i3++){
    for(i2=1; i2<=2; i2++){
      for(i1=1; i1<=2; i1++){
	idf=idf+1;
	if( (pdd[idf] != 0) 
	    && (  idl[0]==0  || idl[0]==i1 )
	    && (  idl[1]==0  || idl[1]==i2 )
	    && (  idl[2]==0  || idl[2]==i3 ) ){
	  //nbdp=nbdp+1;
	  idp[nbdp]=idf;
	  nbdp=nbdp+1;
	}
      }
    }
  }
  
  if(nbdp == 0){
    nbe=0;
  }
  else{
    nbe=3;
    idf=idp[0];
    idecou=pdd[idf];
   /* i=idf;  
    j=i/4;
    i=i-4*j;
    idl[2]=j+1;
    j=i/2;
    idl[1]=j+1;
    idl[0]=i-2*j+1;
    //cout << "idecou= " << idecou << endl;*/
    for(i=0; i<12;i++){
      nu[i]=mu[idecou-1][i]-1;
      //cout << "i, nu[i] "<< i <<" " << nu[i] << endl;
    }
  }

  }
//----------------------------------------------------------------------- 






// glumesh3D

class listMesh3 { 
public:
  list<const Mesh3 *> *lth;
  void init()  { lth=new list<const Mesh3 *>;}
  void destroy() { delete lth;}
  listMesh3(Stack s,const Mesh3 *th) : lth(Add2StackOfPtr2Free(s,new list<const Mesh3*>)) { lth->push_back(th);}
  listMesh3(Stack s,const Mesh3 *tha,const Mesh3 *thb) : lth(Add2StackOfPtr2Free(s,new list<const Mesh3*>)) { lth->push_back(tha);lth->push_back(thb);}
  listMesh3(Stack s,const listMesh3 &l,const Mesh3 *th) : lth(Add2StackOfPtr2Free(s,new list<const Mesh3*>(*l.lth))) { lth->push_back(th);}

};

Mesh3 * GluMesh3(listMesh3 const & lst)
{ 
  int flagsurfaceall = 0;

  int nbt=0;
  int nbe=0;
  int nbex=0;
  int nbv=0;
  int nbvx=0;
  
  double hmin=1e100;
  R3 Pn(1e100,1e100,1e100),Px(-1e100,-1e100,-1e100);
  const list<const Mesh3 *> lth(*lst.lth);
  const Mesh3 * th0=0;
  int kk=0; 
  for(list<const Mesh3 *>::const_iterator i=lth.begin();i != lth.end();i++)
    {
      if( ! *i) continue ;
      kk++;
      const Mesh3 &Th3(**i);  // definis ???
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
  if(kk==0) return 0; // no mesh ....

  if(verbosity > 1) cout << "      - hmin =" <<  hmin << " ,  Bounding Box: " << Pn << " "<< Px << endl;
  
  // probleme memoire
  Vertex3  *v= new Vertex3[nbvx];
  Tet      *t;
  if(nbt!=0) t= new Tet[nbt];
  Tet      *tt=t;
  Triangle3 *b= new Triangle3[nbex];
  Triangle3 *bb= b;
  
  ffassert(hmin>Norme2(Pn-Px)/1e9);
  double hseuil =hmin/10.; 

  //int *NumSom= new int[nbvx];

  // VERSION morice
  if(verbosity > 1) cout << " creation of : BuildGTree" << endl;   
  EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,Pn,Px,0);  
  
  nbv=0;
  //int nbv0=0;
  for(list<const Mesh3 *>::const_iterator i=lth.begin(); i!=lth.end();i++)
    {
      if( ! *i) continue ;
      const Mesh3 &Th3(**i);
      if(verbosity>1)  cout << " loop over mesh for create new mesh "<< endl;
      if(verbosity>1)  cout << " GluMesh3D + "<< Th3.nv << " " << Th3.nt <<" " << Th3.nbe << endl;
      //nbv0 =+Th3.nv;
      
      for (int ii=0;ii<Th3.nv;ii++){
	const Vertex3 &vi(Th3.vertices[ii]);
	Vertex3 * pvi=gtree->ToClose(vi,hseuil);

	   
	if(!pvi){
	  v[nbv].x = vi.x;
	  v[nbv].y = vi.y;
	  v[nbv].z = vi.z;
	  v[nbv].lab = vi.lab;
	  //NumSom[ii+nbv0] = nbv;
	  gtree->Add( v[nbv] );
	  nbv++;
	}
	/*
	  else{
	  NumSom[ii+nbv0] = pvi-v;
	  assert(pvi-v <nbv); 
	  }
	*/
      }
	  
      for (int k=0;k<Th3.nt;k++){
	const Tet  &K(Th3.elements[k]);
	int iv[4];
	iv[0]=gtree->ToClose(K[0],hseuil)-v;
	iv[1]=gtree->ToClose(K[1],hseuil)-v;
	iv[2]=gtree->ToClose(K[2],hseuil)-v;  
	iv[3]=gtree->ToClose(K[3],hseuil)-v;  
	(tt++)->set(v,iv,K.lab);
      }
      //nbv0 =+Th3.nv;
    }
  
  if(verbosity > 1) cout << " creation of : BuildGTree for border elements" << endl;
  Vertex3  *becog= new Vertex3[nbex];  
  //Vertex3  becog[nbex]; 
  EF23::GTree<Vertex3> *gtree_be = new EF23::GTree<Vertex3>(becog,Pn,Px,0);
  
  double hseuil_border = hseuil/3.;
  //nbv0=0;
  for(list<const Mesh3 *>::const_iterator i=lth.begin();i != lth.end();i++)
    {
      if( ! *i) continue ;
      const Mesh3 &Th3(**i);
      
    for (int k=0;k<Th3.nbe;k++)
      {
	const Triangle3 & K(Th3.be(k));
	int iv[3];
	iv[0]=Th3.operator()(K[0]); 
	iv[1]=Th3.operator()(K[1]); 
	iv[2]=Th3.operator()(K[2]); 

	R cdgx,cdgy,cdgz;
	  
	cdgx = (Th3.vertices[iv[0]].x+ Th3.vertices[iv[1]].x+ Th3.vertices[iv[2]].x)/3.;
	cdgy = (Th3.vertices[iv[0]].y+ Th3.vertices[iv[1]].y+ Th3.vertices[iv[2]].y)/3.;
	cdgz = (Th3.vertices[iv[0]].z+ Th3.vertices[iv[1]].z+ Th3.vertices[iv[2]].z)/3.;
	 
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
	  igluv[0]= gtree->ToClose(K[0],hseuil)-v; //NumSom[iv[0]+nbv0];  
	  igluv[1]= gtree->ToClose(K[1],hseuil)-v; //NumSom[iv[1]+nbv0]; 
	  igluv[2]= gtree->ToClose(K[2],hseuil)-v; //NumSom[iv[2]+nbv0]; 
	 
	  (bb++)->set(v,igluv,K.lab);
	}
      }
    //nbv0 =+Th3.nv;
  }
  delete gtree;
  delete gtree_be;
  delete [] becog;
  
  if(verbosity > 2) cout << " nbv="  << nbv  << endl;
  if(verbosity > 2) cout << " nbvx=" << nbvx << endl;
  if(verbosity > 2) cout << " nbt="  << nbt  << endl;
  if(verbosity > 2) cout << " nbe="  << nbe  << endl;
  if(verbosity > 2) cout << " nbex=" << nbex << endl;
  if(verbosity>1)
    {
      cout << "     Nb of glu3D  point " << nbvx-nbv;
      cout << "     Nb of glu3D  Boundary faces " << nbex-nbe << endl;
    }

  if(nbt==0){
    Mesh3 *mpq= new Mesh3(nbv,nbe,v,b);  
    if(flagsurfaceall==1) mpq->BuildBoundaryElementAdj();
    return mpq;
  }
  else{
    Mesh3 *mpq= new Mesh3(nbv,nbt,nbe,v,t,b);  
 /* 
    mpq->BuildBound();
    if(verbosity > 1) cout << "fin de BuildBound" << endl;
    mpq->BuildAdj();
    if(verbosity > 1) cout << "fin de BuildAdj" << endl;
    mpq->Buildbnormalv();  
    if(verbosity > 1) cout << "fin de Buildnormalv()" << endl;
    mpq->BuildjElementConteningVertex();
    if(verbosity > 1) cout << "fin de ConteningVertex()" << endl;
  */
    mpq->BuildGTree();
    if(verbosity > 2) cout << "fin de BuildGTree()" << endl;
    
    //Add2StackOfPtr2FreeRC(stack,mpq);
  
    return mpq;
  }
}


template<class RR,class AA=RR,class BB=AA> 
struct Op3_addmesh: public binary_function<AA,BB,RR> { 
  static RR f(Stack s,const AA & a,const BB & b)  
  { return RR(s, a, b );} 
};

template<bool INIT,class RR,class AA=RR,class BB=AA> 
struct Op3_setmesh: public binary_function<AA,BB,RR> { 
  static RR f(Stack stack,const AA & a,const BB & b)  
  {
    ffassert(a );
    const pmesh3  p=GluMesh3(b);
    
    if(!INIT && *a){
      //Add2StackOfPtr2FreeRC(stack,*a);
       (**a).destroy();
      cout << "destruction du pointeur" << endl;
    }
    //Add2StackOfPtr2FreeRC(stack,p); //  the pointer is use to set variable so no remove. 
    *a=p;
    return a;
  } 
};

// Movemesh3D

class Movemesh3D_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression xx,yy,zz;
  //Expression  lab,reg;  
  static const int n_name_param =7+2; // add FH for cleanning name  //  "+1" add to reorient tetrahedrons
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,int ii,Stack stack,KN_<long> a ) const
    {  ffassert( ! (nargs[i] && nargs[ii]) ); i= nargs[i] ? i : ii;
	return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long  arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
public:
/*  Movemesh3D_Op(const basicAC_F0 &  args,Expression tth)
    : eTh(tth), xx(0) , yy(0) , zz(0)
  {
    args.SetNameParam(n_name_param,name_param,nargs);
    const E_Array * a1=0 ;
    if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[0]);
    int err =0;
    if( nargs[1] && nargs[7] ) 
	CompileError("uncompatible movemesh3 (Th, region= , reftet=  ");
    if( nargs[2] && nargs[8] ) 
	CompileError("uncompatible movemesh3 (Th, label= , refface=  ");

    if(a1) {
      if(a1->size() !=3) 
	CompileError("movemesh3 (Th,transfo=[X,Y,Z],) ");
      xx=to<double>( (*a1)[0]); 
      yy=to<double>( (*a1)[1]);
      zz=to<double>( (*a1)[2]);
    }    
  }*/
    Movemesh3D_Op(const basicAC_F0 &  args,Expression tth,Expression xxx=0,Expression yyy=0,Expression zzz=0  )
    : eTh(tth), xx(xxx) , yy(yyy) , zz(zzz)
    {
        args.SetNameParam(n_name_param,name_param,nargs);
        const E_Array * a1=0 ;
        if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[0]);
        int err =0;
        if( nargs[1] && nargs[7] )
            CompileError("uncompatible movemesh3 (Th, region= , reftet=  ");
        if( nargs[2] && nargs[8] )
            CompileError("uncompatible movemesh3 (Th, label= , refface=  ");
        
        if(a1) {
            if(a1->size() !=3 ||  xx || yy || zz)
                CompileError("movemesh3 (Th,transfo=[X,Y,Z],) ");
            xx=to<double>( (*a1)[0]);
            yy=to<double>( (*a1)[1]);
            zz=to<double>( (*a1)[2]);
        }    
    } 
    
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type Movemesh3D_Op::name_param[]= {
  {  "transfo", &typeid(E_Array)}, //0
  {  "reftet", &typeid(KN_<long>)},// 1
  {  "refface", &typeid(KN_<long>)},
  {  "ptmerge", &typeid(double)},
  {  "facemerge",&typeid(long)},
  {  "boolsurface",&typeid(long)}, // 5
  {  "orientation",&typeid(long)},
  {  "region", &typeid(KN_<long> )}, //7 
  {  "label", &typeid(KN_<long> )} // 8
  // option a rajouter
  // facemerge 0,1 + label
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
  if(verbosity>5)  cout << "before movemesh: Vertex " << nbv<< " Tetrahedra " << nbt << " triangles "<< nbe << endl; 
 
 // lecture des references
   
  KN<long> zzempty;
  KN<long> nrtet  (arg(1,7,stack,zzempty));  
  KN<long> nrf (arg(2,8,stack,zzempty)); 
  double precis_mesh( arg(3,stack,1e-7));
  long  mergefacemesh( arg(4,stack,1L) );
  long  flagsurfaceall( arg(5,stack,0L) );
    long orientationelement( arg(6,stack,1L) );


  //if( nrtet.N() && nrfmid.N() && nrfup.N() && nrfdown.N() ) return m;
  ffassert( nrtet.N() %2 ==0);
  ffassert( nrf.N() %2 ==0);
  
  map<int,int> mapface;
  for(int i=0;i<nrf.N();i+=2)
    {
      if( nrf[i] != nrf[i+1] ){	
	mapface[nrf[i]] = nrf[i+1];
      }      
    }
  
  map<int,int> maptet;
  for(int i=0;i<nrtet.N();i+=2)
    {
      if( nrtet[i] != nrtet[i+1] ){	
	maptet[nrtet[i]] = nrtet[i+1];
      }      
    }



  // realisation de la map par default
  
  assert((xx) && (yy) && (zz) );
  
  KN<double> txx(Th.nv), tyy(Th.nv), tzz(Th.nv);
  
  Mesh3 &rTh3 = Th;
 
  KN<int> takemesh(Th.nv);
  MeshPoint *mp3(MeshPointStack(stack)); 
  
  takemesh=0;
  // loop over tetrahedral 
  for (int it=0;it<Th.nt;++it){
    for( int iv=0; iv<4; ++iv){
      int i=Th(it,iv);  
      
      if(takemesh[i]==0){
	mp3->setP(&Th,it,iv);
          if(xx){ txx[i]=GetAny<double>((*xx)(stack));} else txx[i]=mp3->P.x;
	if(yy){ tyy[i]=GetAny<double>((*yy)(stack));}else tyy[i]=mp3->P.y;
	if(zz){ tzz[i]=GetAny<double>((*zz)(stack));}else tzz[i]=mp3->P.z;
	takemesh[i] = takemesh[i]+1;
      }
    }
  }

  // loop over border elements
   // loop over tetrahedral 
  for (int it=0;it<Th.nbe;++it){
    const Triangle3 &K(Th.be(it));
    int iv[3];
    iv[0]=Th.operator()(K[0]); 
    iv[1]=Th.operator()(K[1]);
    iv[2]=Th.operator()(K[2]);
    
    R coordx,coordy,coordz;
    for(int jj=0; jj< 3; jj++){
      int i=iv[jj];
      if(takemesh[i]==0){	 
	mp3->set( Th.vertices[i].x, Th.vertices[i].y, Th.vertices[i].z );
	if(xx){ txx[i]=GetAny<double>((*xx)(stack));}
	if(yy){ tyy[i]=GetAny<double>((*yy)(stack));}
	if(zz){ tzz[i]=GetAny<double>((*zz)(stack));}
	takemesh[i] = takemesh[i]+1;
      }
    }
  }
  
  // option (Transfo_Mesh3) :: 
  
  // border_only = 0, recollement_border=1, point_confondus_ok=0;   == > 1900 triangles
  // border_only = 0, recollement_border=0, point_confondus_ok=0;   == > 1980 triangles
  // border_only = 0, recollement_border=1, point_confondus_ok=1;   == > 1820 triangles

  // border_only = 1, recollement_border=1, point_confondus_ok=0;   == > 1900 triangles
  // border_only = 1, recollement_border=0, point_confondus_ok=0;   == > 1980 triangles
  // border_only = 1, recollement_border=1, point_confondus_ok=1;   == > 1820 triangles


  int border_only=0; // ne sert a rien !!!!! A enlever
  int recollement_elem=0;
  int recollement_border, point_confondus_ok;
  
  if(mergefacemesh == 0)
    {
      recollement_border=0;
      point_confondus_ok=0;
    }
  if(mergefacemesh == 1)
    { 
      recollement_border=1;
      point_confondus_ok=0;
    }
  if(mergefacemesh == 2)
    { 
      recollement_border=1;
      point_confondus_ok=1;
    }
  
  Mesh3 *T_Th3=Transfo_Mesh3( precis_mesh,rTh3, txx, tyy, tzz, border_only, 
			      recollement_elem, recollement_border, point_confondus_ok,orientationelement);
  if((T_Th3->mes) <= 0 && (T_Th3->nt > 0))
  {
      cerr << " Erreur bad orientation in movemesh add parmetre orientation=,1 mesure="<<T_Th3->mes<< endl;
      ExecError(" movemesh(3d): mesh with neg vol");
  }
  if( nrtet.N() >0)
   {
    for (int i=0;i<nbt;i++)
      {
	const Tet &K( T_Th3->elements[i] );	
	int lab=K.lab;
	T_Th3->elements[i].lab = ChangeLab3D(maptet,lab);
      }  
  }
  // les arete frontieres qui n'ont pas change
    
  if( nrf.N()>0)
    {
    for (int i=0;i<nbe;i++)
      { 
	const Triangle3 &K( T_Th3->be(i) );
	int l0,l1=ChangeLab3D(mapface,l0=K.lab) ;
	T_Th3->be(i).lab = l1;
      }
   }
  

   if(flagsurfaceall==1) T_Th3->BuildBoundaryElementAdj();

   T_Th3->BuildGTree();
   Add2StackOfPtr2FreeRC(stack,T_Th3);
 
  *mp=mps;
  return T_Th3;
}

class  Movemesh3D : public OneOperator { public:
    int cas;
    Movemesh3D() : OneOperator(atype<pmesh3>(),atype<pmesh3>()), cas(0) {}
    Movemesh3D(int ) : OneOperator(atype<pmesh3>(),atype<pmesh3>(),atype<E_Array>()), cas(1) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	if(cas==0)
            return  new Movemesh3D_Op(args,t[0]->CastTo(args[0]));
        else if (cas == 1)
        {
            const E_Array * a = dynamic_cast<const E_Array *>(args[1].LeftValue());
            
            ffassert(a);
            if (a->size() !=3) CompileError("movemesh(Th,[X,Y,Z],...) need 3 componates in array ",atype<pmesh>());
            Expression X=to<double>( (*a)[0]);
            Expression Y=to<double>( (*a)[1]);
            Expression Z=to<double>( (*a)[2]);
            return   new Movemesh3D_Op(args,t[0]->CastTo(args[0]),X,Y,Z);
        }
      else return 0;
  }
};










//// version 3D de change label

class SetMesh3D_Op : public E_F0mps 
{
public:
  Expression a; 
  
  static const int n_name_param =2+2+2+2; //
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
    { ffassert( !(nargs[i] && nargs[i+2]));
      i = nargs[i] ? i : i+2;
      return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
    long  arg(int i,Stack stack, long  a ) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
    bool   arg(int i,Stack stack, bool   a ) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}

  
public:
  SetMesh3D_Op(const basicAC_F0 &  args,Expression aa) : a(aa) {
    args.SetNameParam(n_name_param,name_param,nargs);
      if( nargs[0] && nargs[2] ) 
	  CompileError("uncompatible change(... region= , reftet=  ");
      if( nargs[1] && nargs[3] ) 
	  CompileError("uncompatible  change(...label= , refface=  ");
      
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type SetMesh3D_Op::name_param[]= {
  {  "reftet", &typeid(KN_<long> )},
  {  "refface", &typeid(KN_<long> )},
  {  "region", &typeid(KN_<long> )},
  {  "label", &typeid(KN_<long> )},
  {  "fregion", &typeid(long)},
  {  "flabel", &typeid(long )},
  {  "rmlfaces", &typeid(long)},
  {  "rmInternalFaces", &typeid(bool)}

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
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*a)(stack));
  Mesh3 & Th=*pTh;
  if(!pTh) return pTh;
  Mesh3 *m= pTh;
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.nbe; // nombre d'aretes fontiere
  //cout << "Number of Vertex="<< nbv << "Number of BorderElement=" << nbe << endl;  
  KN<long> zz;
  KN<long> nrtet (arg(0,stack,zz));  
  KN<long> nrface (arg(1,stack,zz));  
  Expression freg = nargs[4];
  Expression flab = nargs[5];
    bool  rm_faces = nargs[6];
    long  rmlabfaces (arg(6,stack,0L));
    bool  rm_i_faces (arg(7,stack,false));

 // cout << " Chnage " << freg << " " << flab << endl;   
  if(nrface.N() <=0 && nrtet.N() <=0 && (!freg) && (!flab) && !rmlabfaces && !rm_i_faces ) return m; // modf J.M. oct 2010
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
  Tet       *t;
  if(nbt!=0) t=new Tet[nbt];
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
  int lmn= 2000000000;
  int lmx= -2000000000;
  int nberr=0;
   R3 PtHat(1./4.,1./4.,1./4.);
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
      
      tt->set( v, iv, ChangeLab3D(maptet,lab));
      if(freg)
	{//      R3 B(1./4.,1./4.,1./4.);  // 27/09/10 : J.Morice error in msh3.cpp
	      mp->set(Th,K(PtHat),PtHat,K,0);
	      tt->lab =GetAny<long>( (* freg)(stack)) ;  
	      lmn= min (lmn,tt->lab);
	      lmx= max (lmx,tt->lab);
	    }
      tt++;
	
    }  
    if(freg && verbosity> 1 ) cout << "    -- Change : new region number bound : " << lmn << " "<< lmx << endl;
    // les arete frontieres qui n'ont pas change
     lmn= 2000000000;
     lmx= -2000000000;
  
  Triangle3 * bb=b;
     R2 PtHat2(1./3.,1./3.);
    int nrmf=0; 
  for (int i=0;i<nbe;i++)
    { 
      const Triangle3 &K( Th.be(i) );
      int fk,ke = Th.BoundaryElement(i,fk); // element co
      int fkk,kke = Th.ElementAdj(ke,fkk=fk); // element co
      bool onborder = (kke==ke) || (kke <0) ;
      const Tet & KE(Th[ke]);
      R3 B= KE.PBord(fk,PtHat2);
      int iv[3];       
      bool  rmf = rm_i_faces && ! onborder; 
      iv[0] = Th.operator()(K[0]);
      iv[1] = Th.operator()(K[1]);
      iv[2] = Th.operator()(K[2]);
      
      int l0,l1=ChangeLab3D(mapface,l0=K.lab) ;
      if(flab)
	  {//      R3 B(1./4.,1./4.,1./4.);  // 27/09/10 : J.Morice error in msh3.cpp
              R3 NN= KE.N(fk);
              double mes=NN.norme();
              NN/=mes;
	      mp->set(Th,KE(B),B,KE,K.lab,NN,fk);
	      l1 =GetAny<long>( (* flab)(stack)) ;
	      lmn= min (lmn,bb->lab);
	      lmx= max (lmx,bb->lab);
	  }
        if( !rmf && rm_faces)
          rmf =  !onborder &&  ( l1 == rmlabfaces  );
        if(rmf)
            nrmf++;
        else
         (*bb++).set( v, iv, l1);

	
    }
    if(nrmf && verbosity > 2) cout << "   change  mesh3 : number of removed  internal faces " << nrmf << " == " << nbe - (bb-b) << endl;
    
  nben -= nrmf;
  nbe -= nrmf;
  assert(nben==bb-b);
  *mp=mps; 
  if(nbt != 0)
    {
      Mesh3 *mpq = new Mesh3(nbv,nbt,nbe,v,t,b);
      
	//mpq->BuildBound();
     // mpq->BuildAdj();
     // mpq->Buildbnormalv();  
     // mpq->BuildjElementConteningVertex(); 
      mpq->BuildGTree();
      //mpq->decrement();   // ?? decrement enlever ???
      Add2StackOfPtr2FreeRC(stack,mpq);  
      
      return mpq;
    }
  if(nbt == 0)
    {
      Mesh3 *mpq = new Mesh3(nbv,nbe,v,b);
      
     // mpq->BuildBound();
      Add2StackOfPtr2FreeRC(stack,mpq);  
      
      return mpq;
    }

  Mesh3 *mpq = NULL;
  return mpq;
}


class SetMesh3D : public OneOperator { public:  
typedef const Mesh3 *pmesh3;
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
  static const int n_name_param =5+1; 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const{ 
      return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  long arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
  double arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny<double>( (*nargs[i])(stack) ): a;}
public:
  Movemesh2D_3D_surf_Op(const basicAC_F0 &  args,Expression tth) : 
  eTh(tth),xx(0),yy(0),zz(0) 
  {

    args.SetNameParam(n_name_param,name_param,nargs);
    
    const E_Array * a1=0 ;
    if(nargs[0])  a1  = dynamic_cast<const E_Array *>(nargs[0]);
    int err =0;
    if( nargs[2] && nargs[5] ) 
    	CompileError("uncompatible movemesh23 (Th, label= , refface=  ");
    if(a1) {
      if(a1->size() !=3) 
	CompileError("movemesh23 (Th,transfo=[X,Y,Z],) ");
      xx=to<double>( (*a1)[0]);
      yy=to<double>( (*a1)[1]);
      zz=to<double>( (*a1)[2]);
    }    
    
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type Movemesh2D_3D_surf_Op::name_param[]= {
  {  "transfo", &typeid(E_Array )},
  {  "orientation", &typeid(long)},
  {  "refface", &typeid(KN_<long>)},
  {  "ptmerge", &typeid(double)},
  {  "boolsurface",&typeid(long)},
  {  "label", &typeid(KN_<long> )}
    
};

AnyType Movemesh2D_3D_surf_Op::operator()(Stack stack)  const 
{
  const Mesh * pTh= GetAny<const Mesh *>((*eTh)(stack));
  const Mesh & Th=*pTh;
  const Mesh *m= pTh;
  int nbv=Th.nv;  // nombre de sommet 
  int nbt=Th.nt;  // nombre de triangles
  int nbe=Th.neb; // nombre d'aretes fontiere
  if(verbosity>5)
    cout << "before movemesh: Vertex Triangle Edge"<< nbv << " "  << nbt << " " << nbe << endl;  
  
  KN<long> zzempty;
  //int intempty=0;
  int mesureM (arg(1,stack,0L));
  KN<long> nrface (arg(2,stack,arg(5,stack,zzempty)));
  double precis_mesh(arg(3,stack,-1.));
  long flagsurfaceall(arg(4,stack,-1L));
  
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
  if( mesureM <0 ){
    surface_orientation=-1;
  }
        

  KN<double> txx(nbv), tyy(nbv), tzz(nbv);
  MeshPoint *mp3(MeshPointStack(stack)); 
  
  {
    KN<int> takemesh(nbv);  
    takemesh=0;  
    const Mesh &rTh = Th;
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
  
  //Mesh3 *Th3; //= new Mesh3;
  
  int vertex_out=1;
  
  if( vertex_out == 1){
    /* determinate the same vertex */ 
    int border_only = 0;
    int recollement_border=1, point_confondus_ok=0;

    // faire version de Transfo_Mesh2_tetgen pour ce cas précis.
    Mesh3 *Th3= MoveMesh2_func( precis_mesh, Th, txx, tyy, tzz, 
			 border_only, recollement_border, point_confondus_ok);
	
    // Rajouter fonction flip a l interieure
    int nbflip=0;
    for(int ii=0; ii < Th3->nbe; ii++){
      const Triangle3 & K(Th3->be(ii)); 
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
	nbflip++;
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
    
    assert(nbflip==0 || nbflip== Th3->nbe);
    if(flagsurfaceall==1) Th3->BuildBoundaryElementAdj();
    Add2StackOfPtr2FreeRC(stack,Th3);
   
    return Th3;
  }
  else   if( vertex_out == 0){
	  
    //Tet       *t = new Tet[1];
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
      	
      	(*bb++).set( v, iv, K.lab);
		
      }
      
    //Mesh3 *Th3 = new Mesh3(nbv,0,nbt,v,t,b);
    Mesh3 *Th3 = new Mesh3(nbv,nbt,v,b);

    int nbflip=0;
    for (int i=0;i<Th3->nbe;i++)
      { 
	double mes_triangle3= Th3->be(i).mesure();
	
	if( surface_orientation*mes_triangle3 < 0){
	  const Triangle3 &K( Th3->be(i) );
	  int iv[3];       
	  
	  iv[0] = Th3->operator()(K[0]);
	  iv[1] = Th3->operator()(K[1]);
	  iv[2] = Th3->operator()(K[2]);
	  
	  int iv_temp=iv[1];
	  iv[1]=iv[2];
	  iv[2]=iv_temp;
	  Th3->be(i).set( Th3->vertices, iv, K.lab ) ;
	  nbflip++;
	}
      }
    assert(nbflip==0 || nbflip== Th3->nbe);
    if(flagsurfaceall==1) Th3->BuildBoundaryElementAdj();
    Add2StackOfPtr2FreeRC(stack,Th3);
    return Th3;
  }
  return (Mesh3*) 0;
}


class Movemesh2D_3D_surf : public OneOperator { public:  
typedef const Mesh *pmesh;
typedef const Mesh3 *pmesh3;
    
  Movemesh2D_3D_surf() : OneOperator(atype<pmesh3>(),atype<pmesh>() ) {}
    E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new Movemesh2D_3D_surf_Op(args,t[0]->CastTo(args[0]));  // CastTo(args[]); // plus tard
  }
};


/* ancien fichier de TransfoMesh */

Mesh3 * Transfo_Mesh3(const double &precis_mesh,const Mesh3 & Th3, const double *tab_XX, const double *tab_YY, const double *tab_ZZ,
		      int &border_only, int &recollement_element, int &recollement_border, int &point_confondus_ok,int orientation)
{
    // cas besoin memoire important
    
    //Mesh3 *T_Th3=new Mesh3;
    int nv_t,nt_t,nbe_t;
    
    int* Numero_Som;
    
    int* ind_nv_t;
    int* ind_nt_t;
    int* ind_nbe_t;
    
    int* label_nt_t;
    int* label_nbe_t;
    
    int i_som, i_elem, i_border;
    
    Numero_Som = new int[Th3.nv];
    
    ind_nv_t   = new int[Th3.nv];
    ind_nt_t   = new int[Th3.nt];
    ind_nbe_t  = new int[Th3.nbe];
    
    label_nt_t   = new int[Th3.nt];
    label_nbe_t  = new int[Th3.nbe];
    
    
    //cout << "Vertex, Tetrahedra, Border : "<<Th3.nv << ", "<<Th3.nt<< ", " << Th3.nbe<< endl;
    
    for(int ii=0; ii<Th3.nv; ii++){
        Numero_Som[ii]=ii;
    }
    
    if(verbosity > 1) cout <<" debut: SamePointElement " <<endl;
    
    SamePointElement( precis_mesh, tab_XX, tab_YY, tab_ZZ, Th3, recollement_element, recollement_border, point_confondus_ok,
                     Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, label_nt_t, label_nbe_t, nv_t, nt_t, nbe_t);
    
    if(verbosity > 1) cout <<" fin: SamePointElement " <<endl;
    
    
    
    // set size of Mesh T_Th3
    //T_Th3->set(nv_t,nt_t,nbe_t);
    Vertex3 *v = new Vertex3[nv_t];
    Tet     *t = new Tet[nt_t];
    Tet     *tt=t;
    Triangle3 *b= new Triangle3[nbe_t];
    Triangle3 *bb=b;
    double mes=0,mesb=0;
    cout << "Transfo TH3 : Vertex, Tetrahedra, Border : "<< "nv_t="<< nv_t << " nt_t=" << nt_t << " nbe_t=" << nbe_t << endl;
    
    // determination of vertex
    i_som = 0;
    for(int i=0; i<nv_t; i++){
        
        int & ii = ind_nv_t[i];
        assert( Numero_Som[ii] == i_som );
        
        const Vertex3 & K(Th3.vertices[ii]);
        
        v[i_som].x = tab_XX[ii];
        v[i_som].y = tab_YY[ii];
        v[i_som].z = tab_ZZ[ii];
        v[i_som].lab = K.lab;
        
        
        i_som = i_som + 1;
    }
    assert( i_som == nv_t);
    
    
    //cout << " Transfo volume elements " << endl;
    // determination of volume elements
    i_elem = 0;
    for( int i=0; i< nt_t; i++){
        int & ii=ind_nt_t[i];
        
        // creation of elements
        
        const Tet & K(Th3.elements[ii]);
        int iv[4];
        int lab;
        lab = label_nt_t[i];
        
        for(int jj=0; jj <4; jj++){
            iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
            assert( iv[jj] >= 0 && iv[jj] < nv_t);
        }
        if(orientation<0) swap(iv[1],iv[2]) ;
        (tt)->set(v, iv, lab);
        mes += tt++->mesure();
        i_elem++;
    }
    
    assert( i_elem == nt_t);
    
    //cout << " Transfo border elements " << endl;
    // determination of border elements
    i_border= 0;
    for( int i=0; i< nbe_t; i++){
        int & ii=ind_nbe_t[i];
        // creation of elements
        const Triangle3 & K(Th3.be(ii));
        int iv[3];
        int lab;
        lab = label_nbe_t[i];
        
        for(int jj=0; jj <3; jj++)
        {
            iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
            assert( iv[jj] >= 0 && iv[jj] < nv_t);
        }
        if(orientation<0) swap(iv[1],iv[2]) ;
        
        bb->set(v, iv, lab);
        mesb += bb++->mesure();
        i_border=i_border+1;
    }
    assert( i_border == nbe_t);
    if  ( mes  < 0 )
    {
        cerr << " E rror of mesh orientation , current orientation = "<< orientation << endl;
        cerr << " volume mesh = " << mes << endl;
        cerr << " surface border mesh = " << mesb << endl;
        ErrorExec(" movemesh 3d ",1);
    }
    
    delete [] Numero_Som;
    delete [] ind_nv_t;   
    delete [] ind_nt_t;  
    delete [] ind_nbe_t;
    delete [] label_nt_t;
    delete [] label_nbe_t;
    
    if( nt_t !=0){
        Mesh3 *T_Th3 = new Mesh3(nv_t,nt_t,nbe_t,v,t,b);
        
        return T_Th3;
    }
    else{
        Mesh3 *T_Th3 = new Mesh3(nv_t,nbe_t,v,b);
        
        delete [] t;
        return T_Th3;
    }
    
    
}
void SamePointElement( const double &precis_mesh, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, 
	int &recollement_element, int &recollement_border, int &point_confondus_ok,
	int *Numero_Som, int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, 
	int *label_nt_t, int *label_nbe_t, int & nv_t, int & nt_t,int & nbe_t ){
		
  int Elem_ok, Border_ok;
  double hmin,hmin_elem,hmin_border;
  R3 bmin,bmax;
  //int recollement_element=1,recollement_border=1;
  
  if(verbosity > 2) cout << "    BuilBound " <<endl;
  BuildBoundMinDist_th3( precis_mesh, tab_XX, tab_YY, tab_ZZ, Th3, bmin, bmax, hmin);
  if(verbosity > 2) cout << "   =============================== " << endl;
		
  double bmin3[3], bmax3[3];
  bmin3[0] = bmin.x;
  bmin3[1] = bmin.y;
  bmin3[2] = bmin.z;
		
  bmax3[0] = bmax.x;
  bmax3[1] = bmax.y;
  bmax3[2] = bmax.z;
 


  if(verbosity > 2) cout << "    OrderVertexTransfo_hcode gtree " << endl;
  OrderVertexTransfo_hcode_nv_gtree( Th3.nv, bmin, bmax, hmin, tab_XX, tab_YY, tab_ZZ, Numero_Som, ind_nv_t, nv_t );
  if(verbosity > 2) cout << "    fin order vertex gtree: nv_t=" << nv_t << endl;
  if(verbosity > 2) cout << "   =============================== " << endl;
		
  /* determination de nt_t et de nbe_t*/ 
  int i_elem, i_border;
		
  i_elem = 0;
  for(int ii=0; ii< Th3.nt; ii++){
    const Tet & K(Th3.elements[ii]);
    int iv[4];
			
    Elem_ok = 1;
			
    for(int jj=0; jj <4; jj++){
      iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
    }
			
    for(int jj=0; jj<4; jj++){
      for(int kk=jj+1; kk<4; kk++){
	if( iv[jj]==iv[kk] ){
	  Elem_ok = 0;
	}
      }
    }
			
    if(Elem_ok==1){
      ind_nt_t[i_elem]= ii;
      label_nt_t[i_elem] = K.lab;
      i_elem = i_elem + 1;
    }
  } 
  nt_t=i_elem;
		
  if(recollement_element ==1){
    //int point_confondus_ok_e = 0;
    if(verbosity > 1) cout << "debut recollement : nt_t= "<< nt_t << endl; 
			
    int np,dim=3;
    int *ind_np = new int [nt_t];
    int *label_t = new int [nt_t];
    double **Cdg_t=new double *[nt_t];
    for(int i=0; i<nt_t; i++) Cdg_t[i] = new double[dim];
			
    for( int i_elem=0; i_elem< nt_t; i_elem++){
      int & ii=ind_nt_t[i_elem];
      const Tet & K(Th3.elements[ii]);
      int iv[4];
      for(int jj=0; jj <4; jj++){
	iv[jj] = Th3.operator()(K[jj]) ;
      }
      Cdg_t[i_elem][0] = ( tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]] + tab_XX[iv[3]] )/4.;
      Cdg_t[i_elem][1] = ( tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]] + tab_YY[iv[3]] )/4.;
      Cdg_t[i_elem][2] = ( tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]] + tab_ZZ[iv[3]] )/4.;
      label_t[i_elem]  = K.lab;
    }
			
    hmin_elem = hmin/4;
    //PointCommun_hcode( dim, nt_t, 0, Cdg_t, bmin3, bmax3, hmin_elem, ind_np, np); //ancien
    PointCommun_hcode_gtree( dim, nt_t, 0, Cdg_t, label_t, bmin, bmax, hmin_elem, 
			     ind_np, label_nt_t, np); // nv
			
    assert( np <= nt_t );
			
    int *ind_nt_t_tmp= new int [np];
			
    for( int i_elem=0; i_elem< np; i_elem++){
      assert( ind_np[i_elem] >=0 && ind_np[i_elem] <= nt_t );
      ind_nt_t_tmp[ i_elem ] = ind_nt_t[ ind_np[i_elem] ]; 
    }
    for( int i_elem=0; i_elem< np; i_elem++){
      ind_nt_t[ i_elem ] = ind_nt_t_tmp[ i_elem ]; 
    }
			
    
    delete [] ind_np; 
    delete [] label_t;
    for(int i=0; i<nt_t; i++) delete [ ] Cdg_t[i];
    delete [] Cdg_t;

    delete [] ind_nt_t_tmp;
		
    nt_t = np;
    if(verbosity > 1) cout << "fin recollement : nt_t= "<< nt_t << endl; 
  }
		
  // determination of border elements
  i_border= 0;
  for( int ii=0; ii< Th3.nbe; ii++){
    Border_ok=1;
		
    const Triangle3 & K(Th3.be(ii));
    int iv[3];
			
    for(int jj=0; jj <3; jj++){
      iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
      assert( iv[jj] >= 0 && iv[jj] < nv_t);
    }
			
    for(int jj=0; jj<3; jj++){
      for(int kk=jj+1; kk<3; kk++){
	if( iv[jj]==iv[kk] ) Border_ok=0;
      }
    }
    if(Border_ok==1){
      ind_nbe_t[i_border]   = ii;
      label_nbe_t[i_border] = K.lab;
      i_border=i_border+1;
    }
  }
  nbe_t = i_border;
		
  if( recollement_border == 1){
    //int point_confondus_ok = 1;
    if(verbosity > 1) cout << "debut recollement : nbe_t= "<< nbe_t << endl; 
			
    int np,dim=3;
    int *ind_np = new int [nbe_t];
    double **Cdg_be=new double *[nbe_t];
    int *label_be = new int [nbe_t];
    for(int i=0; i<nbe_t; i++) Cdg_be[i] = new double[dim];
			
    for( int i_border=0; i_border< nbe_t; i_border++){
				
      int & ii=ind_nbe_t[i_border];
      const Triangle3 & K(Th3.be(ii));
      int iv[3];
				
      for(int jj=0; jj <3; jj++){
	iv[jj] = Th3.operator()(K[jj]);
      }
      Cdg_be[i_border][0] = ( tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]] )/3.; //( Th3.vertices[iv[0]].x + Th3.vertices[iv[1]].x + Th3.vertices[iv[2]].x )/3.;
      Cdg_be[i_border][1] = ( tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]] )/3.; //( Th3.vertices[iv[0]].y + Th3.vertices[iv[1]].y + Th3.vertices[iv[2]].y )/3.;
      Cdg_be[i_border][2] = ( tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]] )/3.; //( Th3.vertices[iv[0]].z + Th3.vertices[iv[1]].z + Th3.vertices[iv[2]].z )/3.;		
				
      label_be[i_border] = K.lab;
    }
    hmin_border=hmin/3.;
    if(verbosity > 1) cout << "hmin_border=" << hmin_border << endl;
			
    if(verbosity > 1) cout << "appele de PointCommun_hcode := " << point_confondus_ok<< endl;
    //PointCommun_hcode( dim, nbe_t, point_confondus_ok, Cdg_be, bmin3, bmax3, hmin_border, ind_np, np);
    PointCommun_hcode_gtree( dim, nbe_t, point_confondus_ok, Cdg_be, label_be, 
			     bmin, bmax, hmin_border, ind_np, label_nbe_t, np); 
    if(verbosity > 1) cout << "fin appele de PointCommun_hcode" << endl;
			
    assert( np <= nbe_t );
			
    int *ind_nbe_t_tmp= new int [np];
		
    for( int i_border=0; i_border<np; i_border++){
      ind_nbe_t_tmp[ i_border ] = ind_nbe_t[ ind_np[i_border] ]; 
    }
		
    for( int i_border=0; i_border< np; i_border++){
      ind_nbe_t[ i_border ] = ind_nbe_t_tmp[ i_border ]; 
    }
   
    delete [] ind_np; 
    delete [] label_be;
    for(int i=0; i<nbe_t; i++) delete [ ] Cdg_be[i];
    delete [] Cdg_be;

    delete [] ind_nbe_t_tmp;
			
    nbe_t = np;
    if(verbosity > 1) cout << "fin recollement : nbe_t= "<< nbe_t << endl; 
			
    // Affectation de la nouvelle valeur du label
			
  }
}

// 3D surface

Mesh3 * Transfo_Mesh3_surf(const double &precis_mesh, const Mesh3 & Th3, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
	int &recollement_border, int &point_confondus_ok){
	// cas besoin memoire important
	
	//Mesh3 *T_Th3=new Mesh3;
	int nv_t,nbe_t;
	int nt_t=0;
	
	int* Numero_Som;
	int* ind_nv_t;
	int* ind_nbe_t;
	int* label_nbe_t;
	
	int i_som, i_elem, i_border;
	
	assert( Th3.nt == 0);
	
	Numero_Som = new int[Th3.nv];
	ind_nv_t   = new int[Th3.nv];
	ind_nbe_t  = new int[Th3.nbe];
	label_nbe_t  = new int[Th3.nbe];
	
    
	if(verbosity > 1) cout << "Vertex, Tetrahedra, Border : "<<Th3.nv << ", "<<Th3.nt<< ", " << Th3.nbe<< endl;
    
	for(int ii=0; ii<Th3.nv; ii++){
		Numero_Som[ii]=ii;
	}
	
	if(verbosity > 1) cout <<" debut: SamePointElement " <<endl;
	
	SamePointElement_surf( precis_mesh, tab_XX, tab_YY, tab_ZZ, Th3, 
			       recollement_border, point_confondus_ok, Numero_Som,
			       ind_nv_t, ind_nbe_t, label_nbe_t, nv_t, nbe_t);
	
	if(verbosity > 1) cout <<" fin: SamePointElement " <<endl;
	
	// set size of Mesh T_Th3 


	//T_Th3->set(nv_t,nt_t,nbe_t);
	Vertex3 *v = new Vertex3[nv_t];
	//Tet     *t;
	Triangle3 *b= new Triangle3[nbe_t];
	Triangle3 *bb=b;


	if(verbosity > 1) cout << "Transfo TH3 : Vertex, Tetrahedra, Border : "<< "nv_t="<< nv_t << " nt_t=" << nt_t << " nbe_t=" << nbe_t << endl;
		
	// determination of vertex		
	i_som = 0;
	for(int i=0; i<nv_t; i++){
		
	  int & ii = ind_nv_t[i];
	  assert( Numero_Som[ii] == i_som );
	  
	  const Vertex3 & K(Th3.vertices[ii]);
	  /*
	    T_Th3->vertices[i_som].x = tab_XX[ii];
	    T_Th3->vertices[i_som].y = tab_YY[ii];
	    T_Th3->vertices[i_som].z = tab_ZZ[ii];
	    T_Th3->vertices[i_som].lab = K.lab; 
	  */
	   v[i_som].x = tab_XX[ii];
	   v[i_som].y = tab_YY[ii];
	   v[i_som].z = tab_ZZ[ii];
	   v[i_som].lab = K.lab; 

	  i_som = i_som + 1;		
	}	
	if(verbosity > 1) cout << "i_som, nv_t=" <<i_som << " "<<nv_t << endl;
	assert( i_som == nv_t);
	
	if(verbosity > 1) cout << " Transfo border elements " << endl;
	// determination of border elements
	i_border= 0;
	for( int i=0; i< nbe_t; i++){
	  int & ii=ind_nbe_t[i];
	  
	  // creation of elements
	  const Triangle3 & K(Th3.be(ii));
	  int iv[3];
	  int lab;
	  
	  //lab = K.lab; 
	  lab = label_nbe_t[i];
	  
	  for(int jj=0; jj <3; jj++){
	    iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
	    assert( iv[jj] >= 0 && iv[jj] <= nv_t);
	  }
	  //T_Th3->be(i_border).set(T_Th3->vertices, iv, lab);
	  (bb++)->set(v, iv, lab);
	  i_border=i_border+1;
	} 
	assert( i_border == nbe_t);

	delete [] Numero_Som;
	delete [] ind_nv_t;
	delete [] ind_nbe_t;  
	delete [] label_nbe_t;

	//Mesh3* T_Th3 = new Mesh3(nv_t,nt_t,nbe_t,v,t,b); 
	Mesh3* T_Th3 = new Mesh3(nv_t,nbe_t,v,b);

	return T_Th3;
}

void SamePointElement_surf( const double &precis_mesh, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3 & Th3, 
	int &recollement_border, int &point_confondus_ok, int *Numero_Som, 
	int *ind_nv_t, int *ind_nbe_t, int *label_nbe_t, int & nv_t,int & nbe_t ){
		
  int Elem_ok, Border_ok;
  double hmin,hmin_elem,hmin_border;
  R3 bmin,bmax;
  //int recollement_element=1,recollement_border=1;
  
  if(verbosity > 1) cout << "  OrderVertexTransfo_hcode gtree " <<endl;
  BuildBoundMinDist_th3( precis_mesh, tab_XX, tab_YY, tab_ZZ, Th3, bmin, bmax, hmin);
  if(verbosity > 1) cout << " =============================== " << endl;
  
  double bmin3[3], bmax3[3];
  bmin3[0] = bmin.x;
  bmin3[1] = bmin.y;
  bmin3[2] = bmin.z;
  
  bmax3[0] = bmax.x;
  bmax3[1] = bmax.y;
  bmax3[2] = bmax.z;
  
  /*
    cout << "  OrderVertexTransfo_hcode " << endl;
    OrderVertexTransfo_hcode_nv( Th3.nv, tab_XX, tab_YY, tab_ZZ, bmin3, bmax3, hmin, Numero_Som, ind_nv_t, nv_t );
    cout << "fin order vertex: nv_t=" << nv_t << endl;
  */
  if(verbosity > 1) cout << "  OrderVertexTransfo_hcode gtree " << endl;
  OrderVertexTransfo_hcode_nv_gtree( Th3.nv, bmin, bmax, hmin, tab_XX, tab_YY, tab_ZZ, Numero_Som, ind_nv_t, nv_t );
  if(verbosity > 1) cout << "fin order vertex gtree: nv_t=" << nv_t << endl;
  
  if(verbosity > 1) cout << " =============================== " << endl;
  
  /* determination de nt_t et de nbe_t*/ 
  int i_border;
  
  // determination of border elements
  i_border= 0;
  for( int ii=0; ii< Th3.nbe; ii++){
    Border_ok=1;
    
    const Triangle3 & K(Th3.be(ii));
    int iv[3];
			
    for(int jj=0; jj <3; jj++){
      iv[jj] = Numero_Som[ Th3.operator()(K[jj]) ];
    }
			
    for(int jj=0; jj<3; jj++){
      for(int kk=jj+1; kk<3; kk++){
	if( iv[jj]==iv[kk] ) Border_ok=0;
      }
    }
    if(Border_ok==1){
      ind_nbe_t[i_border]   = ii;
      label_nbe_t[i_border] = K.lab;
      i_border=i_border+1;
    }
  }
  nbe_t = i_border;
		
  if( recollement_border == 1){
    //int point_confondus_ok = 1;
    if(verbosity > 1) cout << "debut recollement : nbe_t= "<< nbe_t << endl; 
			
    int np,dim=3;
    int *ind_np = new int [nbe_t];
    int *label_be = new int [nbe_t];

    double **Cdg_be=new double *[nbe_t];
    for(int i=0; i<nbe_t; i++) Cdg_be[i] = new double[dim];
  
		
    for( int i_border=0; i_border< nbe_t; i_border++){
				
      int & ii=ind_nbe_t[i_border];
      const Triangle3 & K(Th3.be(ii));
      int iv[3];
				
      for(int jj=0; jj <3; jj++){
	iv[jj] = Th3.operator()(K[jj]);
      }
      Cdg_be[i_border][0] = ( tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]] )/3.; //( Th3.vertices[iv[0]].x + Th3.vertices[iv[1]].x + Th3.vertices[iv[2]].x )/3.;
      Cdg_be[i_border][1] = ( tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]] )/3.; //( Th3.vertices[iv[0]].y + Th3.vertices[iv[1]].y + Th3.vertices[iv[2]].y )/3.;
      Cdg_be[i_border][2] = ( tab_ZZ[iv[0]] + tab_ZZ[iv[1]] +  tab_ZZ[iv[2]] )/3.; //( Th3.vertices[iv[0]].z + Th3.vertices[iv[1]].z + Th3.vertices[iv[2]].z )/3.;		
				
      label_be[i_border] = K.lab;
    }
    hmin_border=hmin/3.;
    if(verbosity > 1) cout << "hmin_border=" << hmin_border << endl;
			
    if(verbosity > 1) cout << "appele de PointCommun_hcode := " << point_confondus_ok<< endl;
    //PointCommun_hcode( dim, nbe_t, point_confondus_ok, Cdg_be, bmin3, bmax3, hmin_border, ind_np, np);
    PointCommun_hcode_gtree( dim, nbe_t, point_confondus_ok, Cdg_be, label_be, 
			     bmin, bmax, hmin_border, ind_np, label_nbe_t, np); 
    if(verbosity > 1) cout << "fin appele de PointCommun_hcode" << endl;
			
    assert( np <= nbe_t );
			
    int *ind_nbe_t_tmp= new int [np];
		
    for( int i_border=0; i_border<np; i_border++){
      ind_nbe_t_tmp[ i_border ] = ind_nbe_t[ ind_np[i_border] ]; 
    }
		
    for( int i_border=0; i_border< np; i_border++){
      ind_nbe_t[ i_border ] = ind_nbe_t_tmp[ i_border ]; 
    }
  

    delete [] ind_np; 
    delete [] label_be;
    delete [] ind_nbe_t_tmp;
    for(int i=0; i<nbe_t; i++) delete [] Cdg_be[i];
    delete [] Cdg_be;
	
    nbe_t = np;
    if(verbosity > 1) cout << "fin recollement : nbe_t= "<< nbe_t << endl; 
			
    // Affectation de la nouvelle valeur du label
		
  }
}

void Transfo_Mesh2_map_face(const Mesh &Th2, map<int, int> &maptri ){
		
  int numero_label=0;
  for(int ii=0; ii< Th2.nt; ii++){
    const Mesh::Triangle & K(Th2.t(ii));
    map<int,int>::const_iterator imap=maptri.find(K.lab);
			
    if(imap == maptri.end()){
      maptri[ K.lab ] = numero_label;
      numero_label = numero_label+1;
    }
  }
}


Mesh3 * MoveMesh2_func( const double &precis_mesh, const Mesh & Th2, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, 
	int &border_only, int &recollement_border, int &point_confondus_ok){
	
  //Mesh3 *T_Th3= new Mesh3;
	int nv_t,nt_t,nbe_t;
	int* Numero_Som;
	
	int* ind_nv_t;
	int* ind_nt_t=0;
	int* ind_nbe_t;
	int* label_nbe_t;
	
	//int i_som;
	Numero_Som = new int[Th2.nv];
	ind_nv_t   = new int[Th2.nv];
	ind_nbe_t  = new int[Th2.nt];   
	label_nbe_t = new int[Th2.nt];
	if(verbosity>5)
	cout << "before movemesh::Vertex  triangle2  border " << Th2.nv << " "<<Th2.nt<< " " << Th2.neb<< endl;
	
	for(int ii=0; ii<Th2.nv; ii++){ 
		Numero_Som[ii]=ii;
	}
	
	if(verbosity > 1) cout <<" debut: SamePointElement " <<endl;
	
	SamePointElement_Mesh2( precis_mesh, tab_XX, tab_YY, tab_ZZ, Th2, recollement_border, point_confondus_ok, 
		Numero_Som, ind_nv_t, ind_nt_t, ind_nbe_t, label_nbe_t, nv_t, nt_t, nbe_t);
	
	if(verbosity > 1) cout <<" fin: SamePointElement " <<endl;
	
	cout << "After movemesh::Vertex  triangle2  border " << nv_t << " "<< nt_t << " " << nbe_t<< endl;
	
	Vertex3 *v=new Vertex3[nv_t];
	Tet *t;
	Triangle3 *b=new Triangle3[nbe_t];
	Triangle3 *bb=b;

	//T_Th3->set(nv_t,0,nbe_t);
	
	for(int nnv=0; nnv < nv_t; nnv++)
	{
		int ii = ind_nv_t[nnv];
		assert( Numero_Som[ii] == nnv );
		const Mesh::Vertex & K = Th2.vertices[ii];//const Vertex2 & K(Th2.vertices[ii]); //Version Mesh2   
		/*
		  T_Th3->vertices[nnv].x = tab_XX[ii];
		  T_Th3->vertices[nnv].y = tab_YY[ii];
		  T_Th3->vertices[nnv].z = tab_ZZ[ii];       
		  T_Th3->vertices[nnv].lab =  K.lab;   
		*/
		v[nnv].x = tab_XX[ii];
		v[nnv].y = tab_YY[ii];
		v[nnv].z = tab_ZZ[ii];       
		v[nnv].lab =  K.lab;
		
	}
	  
	for(int ibe=0; ibe < nbe_t; ibe++){
		int lab;
		int iv[3];
		int ii=ind_nbe_t[ibe];
		// creation of elements
		const Mesh::Triangle & K(Th2.t(ii)); // const Triangle2 & K(Th2.elements[ii]); // Version Mesh2  
		iv[0] = Numero_Som[ Th2.operator()(K[0]) ];
		iv[1] = Numero_Som[ Th2.operator()(K[1]) ];
		iv[2] = Numero_Som[ Th2.operator()(K[2]) ];
		
		/*
		map< int, int>:: const_iterator imap;
		imap = maptri.find(K.lab); // imap= maptri.find( label_nbe_t[ibe] );
		assert( imap != maptri.end());
		lab = imap->second; // K.lab; // before 
		*/
		//T_Th3->be(ibe).set(T_Th3->vertices,iv,K.lab);
		(bb++)->set(v,iv,K.lab);
	}  


	//Mesh3 *T_Th3 = new Mesh3(nv_t,0,nbe_t,v,t,b);
	Mesh3 *T_Th3 = new Mesh3(nv_t,nbe_t,v,b);

 	delete [ ] Numero_Som;
	delete [ ] ind_nv_t;
	delete [ ] ind_nbe_t;
	delete [ ] label_nbe_t;

	return T_Th3;
}

void SamePointElement_Mesh2( const double &precis_mesh, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh & Th2, 
	int &recollement_border, int &point_confondus_ok, int *Numero_Som, 
	int *ind_nv_t, int *ind_nt_t, int *ind_nbe_t, int* label_nbe_t,
	int & nv_t, int & nt_t,int & nbe_t ){

  int Border_ok;
  //int recollement_border=0;
  R3 bmin,bmax;
  double hmin,hmin_border;
		
  if(verbosity > 1) cout << "calculus of bound and minimal distance" << endl;
  BuildBoundMinDist_th2( precis_mesh, tab_XX, tab_YY, tab_ZZ, Th2, bmin, bmax, hmin);
  // assertion pour la taille de l octree
  assert(hmin>Norme2(bmin-bmax)/1e9);
		
  double bmin3[3], bmax3[3];
  bmin3[0] = bmin.x;
  bmin3[1] = bmin.y;
  bmin3[2] = bmin.z;
		
  bmax3[0] = bmax.x;
  bmax3[1] = bmax.y;
  bmax3[2] = bmax.z;
  /*
    cout << "debut: OrderVertexTransfo_hcode " <<endl;
    OrderVertexTransfo_hcode_nv( Th2.nv, tab_XX, tab_YY, tab_ZZ, bmin3, bmax3, hmin, Numero_Som, ind_nv_t, nv_t );
    cout << "fin order vertex: nv_t=" << nv_t << endl;
  */
  if(verbosity > 1) cout << "debut: OrderVertexTransfo_hcode_gtree " <<endl;
  OrderVertexTransfo_hcode_nv_gtree( Th2.nv, bmin, bmax, hmin, tab_XX, tab_YY, tab_ZZ, Numero_Som, ind_nv_t, nv_t);
  if(verbosity > 1) cout << "fin: OrderVertexTransfo_hcode_gtree " <<endl;
		
		
  /* determination de nt_t et de nbe_t*/ 
  nt_t = 0;
  int i_border;
		
  // determination of border elements
  i_border= 0;
  for( int
	 ii=0; ii< Th2.nt; ii++){
    Border_ok=1;
    const Mesh::Triangle & K(Th2.t(ii)); // const Triangle2 & K(Th2.elements[ii]); // avant Mesh2
    int iv[3];
			
    for(int jj=0; jj <3; jj++){
      iv[jj] = Numero_Som[ Th2.operator()(K[jj]) ];
    }
			
    for(int jj=0; jj<3; jj++){
      for(int kk=jj+1; kk<3; kk++){
	if( iv[jj]==iv[kk] ) Border_ok=0;
      }
    }
    if(Border_ok==1){
      ind_nbe_t[i_border] = ii;
      label_nbe_t[i_border] = K.lab;
      i_border=i_border+1;
    }
  }
  nbe_t = i_border;
		
  if( recollement_border == 1){
    //int point_confondus_ok=1;
    if(verbosity > 1) cout << "debut recollement : nbe_t= "<< nbe_t << endl; 
			
    int np,dim=3;
    int *ind_np = new int [nbe_t];
    int *label_be = new int [nbe_t ];
    double **Cdg_be=new double *[nbe_t];
		
    for(int i=0; i<nbe_t; i++) Cdg_be[i] = new double[dim];
		
    for( int i_border=0; i_border< nbe_t; i_border++){
      int & ii=ind_nbe_t[i_border];
      const Mesh::Triangle & K(Th2.t(ii)); //const Triangle2 & K(Th2.elements[ii]);  // avant Mesh2
      int iv[3];
				
      for(int jj=0; jj <3; jj++){
	iv[jj] = Th2.operator()(K[jj]);
      }
      
      Cdg_be[i_border][0] = ( tab_XX[iv[0]] + tab_XX[iv[1]] + tab_XX[iv[2]] )/3.; 
      Cdg_be[i_border][1] = ( tab_YY[iv[0]] + tab_YY[iv[1]] + tab_YY[iv[2]] )/3.; 
      Cdg_be[i_border][2] = ( tab_ZZ[iv[0]] + tab_ZZ[iv[1]] + tab_ZZ[iv[2]] )/3.; 
				
      label_be[i_border] = K.lab;
    }
			
    hmin_border=hmin/3.;
    if(verbosity > 1) cout << "points commun " << endl;
    //PointCommun_hcode( dim, nbe_t, point_confondus_ok, Cdg_be, bmin3, bmax3, hmin_border, ind_np, np); // ancien
    PointCommun_hcode_gtree( dim, nbe_t, point_confondus_ok, Cdg_be, label_be, bmin, bmax, hmin_border, 
			     ind_np, label_nbe_t,np); // new
    if(verbosity > 1) cout << "points commun finis " <<endl;
    assert( np <= nbe_t );

    	
    //int *ind_nbe_t_tmp= new int [np];
    int ind_nbe_t_tmp[np];
	
    for( int i_border=0; i_border<np; i_border++){
      ind_nbe_t_tmp[ i_border ] = ind_nbe_t[ ind_np[i_border] ]; 
    }
			
    for( int i_border=0; i_border< np; i_border++){
      ind_nbe_t[ i_border ] = ind_nbe_t_tmp[ i_border ]; 
    }	
		
   
    delete [ ] ind_np; //= new int [nbe_t];
    delete [ ] label_be;// = new int [nbe_t ];
    for(int i=0; i<nbe_t; i++) delete [ ] Cdg_be[i];
    delete [ ] Cdg_be; //=new double *[nbe_t];
    
    nbe_t = np;
    if(verbosity > 1) cout << "fin recollement : nbe_t= "<< nbe_t << endl; 

  }
}


//======================
//    Fin cas 2D
//======================
// version Mesh2
void BuildBoundMinDist_th2( const double &precis_mesh, const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh & Th2, R3 &bmin, R3 &bmax, double &hmin){
  //determination de la boite englobante
  //R3 bmin,bmax;
  double precispt;
		
  bmin.x = tab_XX[0];
  bmin.y = tab_YY[0];
  bmin.z = tab_ZZ[0];
		
  bmax.x = bmin.x;
  bmax.y = bmin.y;
  bmax.z = bmin.z;
		
  //R3 bmax = new R3(bmin);
  
  if(verbosity >1) cout << " determination of bmin and bmax" << endl;
		
  for(int ii=1; ii<Th2.nv; ii++){
    bmin.x = min(bmin.x,tab_XX[ii]);
    bmin.y = min(bmin.y,tab_YY[ii]);
    bmin.z = min(bmin.z,tab_ZZ[ii]);
		
    bmax.x = max(bmax.x,tab_XX[ii]);
    bmax.y = max(bmax.y,tab_YY[ii]);
    bmax.z = max(bmax.z,tab_ZZ[ii]);
  }
  double longmini_box=1e10;
		
  longmini_box = pow(bmax.x-bmin.x,2)+pow(bmax.y-bmin.y,2)+ pow(bmax.z-bmin.z,2);
  longmini_box = sqrt(longmini_box);
		
  // determination de hmin 
  if(precis_mesh< 0){
    precispt=longmini_box*1e-7;
  }
  else{
    precispt=precis_mesh;
  }

  hmin = 1e10;
  for( int ii=0; ii< Th2.nt; ii++){
    const Mesh :: Triangle & K(Th2.t(ii)); // const Triangle2 & K(Th2.elements[ii]);
    double longedge;
    int iv[3];
    for(int jj=0; jj<3; jj++){
      iv[jj] = Th2.operator()(K[jj]) ;
    }
			
    for( int jj=0; jj<3; jj++){
      for( int kk=jj+1; kk<3; kk++){ 
	int & i1= iv[jj];
	int & i2= iv[kk];
	longedge = pow(tab_XX[i1]-tab_XX[i2],2) 
	  + pow(tab_YY[i1]-tab_YY[i2],2) 
	  + pow(tab_ZZ[i1]-tab_ZZ[i2],2);
	longedge = sqrt(longedge);
	//cout << "longedge=" << longedge << endl; 
	if( longedge > precispt ) hmin = min( hmin, longedge);
      }
    }
  }
  if(verbosity >5) cout << "    longmin_box=" << longmini_box << endl;
  if(verbosity >5) cout << "    hmin =" << hmin << endl;
  if(verbosity >5) cout << "    Norme2(bmin-bmax)=" <<  Norme2(bmin-bmax) << endl;
  assert( hmin < longmini_box);

  // assertion pour la taille de l octree
  assert(hmin>Norme2(bmin-bmax)/1e9);

  /*  // ?????????
    hmin = 1e10;
    for( int ii=0; ii< Th2.nt; ii++){
      const Mesh :: Triangle & K(Th2.t(ii)); // const Triangle2 & K(Th2.elements[ii]);
      double longedge;
      int iv[3];
      for(int jj=0; jj<3; jj++){
	iv[jj] = Th2.operator()(K[jj]) ;
      }
			
      for( int jj=0; jj<3; jj++){
	for( int kk=jj+1; kk<3; kk++){ 
	  int & i1= iv[jj];
	  int & i2= iv[kk];
	  longedge = pow(tab_XX[i1]-tab_XX[i2],2) 
	    + pow(tab_YY[i1]-tab_YY[i2],2) 
	    + pow(tab_ZZ[i1]-tab_ZZ[i2],2);
	  longedge = sqrt(longedge);
	  //cout << "longedge=" << longedge << endl; 
	  if( longedge > longmini_box*1e-7 ) hmin = min( hmin, longedge);
	}
      }
    }
    cout << "longmin_box=" << longmini_box << endl;
    cout << "hmin =" << hmin << endl;
    cout << "Norme2(bmin-bmax)=" <<  Norme2(bmin-bmax) << endl;
    assert( hmin < longmini_box);
    // assertion pour la taille de l octree
    assert(hmin>Norme2(bmin-bmax)/1e9);
  */
}

// version Mesh3

void BuildBoundMinDist_th3(  const double &precis_mesh,  const double *tab_XX, const double *tab_YY, const double *tab_ZZ, const Mesh3& Th3, R3 &bmin, R3 &bmax, double &hmin){
  //determination de la boite englobante
  //R3 bmin,bmax;
  double precispt;
  bmin.x = tab_XX[0];
  bmin.y = tab_YY[0];
  bmin.z = tab_ZZ[0];
		
  bmax.x = bmin.x;
  bmax.y = bmin.y;
  bmax.z = bmin.z;
		
  //R3 bmax = new R3(bmin);
		
  if(verbosity >1) cout << " determination of bmin and bmax" << endl;
		
  for(int ii=1; ii<Th3.nv; ii++){
    bmin.x = min(bmin.x,tab_XX[ii]);
    bmin.y = min(bmin.y,tab_YY[ii]);
    bmin.z = min(bmin.z,tab_ZZ[ii]);
		
    bmax.x = max(bmax.x,tab_XX[ii]);
    bmax.y = max(bmax.y,tab_YY[ii]);
    bmax.z = max(bmax.z,tab_ZZ[ii]);
  }
		
  double longmini_box;
		
  //longmini_box = min(bmax.x-bmin.x, bmax.y-bmin.y);
  //longmini_box = min(longmini_box, bmax.z-bmin.z);
		
  longmini_box = pow(bmax.x-bmin.x,2)+pow(bmax.y-bmin.y,2)+ pow(bmax.z-bmin.z,2);
  longmini_box = sqrt(longmini_box);
		
		
  if(verbosity >1) cout << " bmin := " << bmin.x << " " << bmin.y << " " << bmin.z << endl;
  if(verbosity >1) cout << " bmax := " << bmax.x << " " << bmax.y << " " << bmax.z << endl;
  if(verbosity >1) cout << " box volume :=" << longmini_box << endl;
	

  if(precis_mesh< 0){
    precispt=longmini_box*1e-7;
  }
  else{
    precispt=precis_mesh;
  }
  // determination de hmin 
		
  hmin = 1e10;
  for( int ii=0; ii< Th3.nt; ii++){
    const Tet & K(Th3.elements[ii]);
    double longedge;
    int iv[4];
    for(int jj=0; jj <4; jj++){
      iv[jj] = Th3.operator()(K[jj]) ;
    }
			
    for( int jj=0; jj<4; jj++){
      for( int kk=jj+1; kk<4; kk++){ 
	int & i1= iv[jj];
	int & i2= iv[kk];
	longedge = pow(tab_XX[i1]-tab_XX[i2],2) 
	  + pow(tab_YY[i1]-tab_YY[i2],2) 
	  + pow(tab_ZZ[i1]-tab_ZZ[i2],2);
	longedge = sqrt(longedge);
	if(longedge > precispt ) hmin = min( hmin, longedge);
      }
    }
  }

  if( Th3.nt == 0){
    for( int ii=0; ii< Th3.nbe; ii++){
      if(verbosity >10) cout << "border " << ii <<" hmin =" << hmin << endl;
      const Triangle3 & K(Th3.be(ii));
      double longedge;
      int iv[3];
      for(int jj=0; jj <3; jj++){
	iv[jj] = Th3.operator()(K[jj]) ;
      }
			
      for( int jj=0; jj<3; jj++){
	for( int kk=jj+1; kk<3; kk++){ 
	  int & i1= iv[jj];
	  int & i2= iv[kk];
	  longedge = pow(tab_XX[i1]-tab_XX[i2],2) 
	    + pow(tab_YY[i1]-tab_YY[i2],2) 
	    + pow(tab_ZZ[i1]-tab_ZZ[i2],2);
	  longedge = sqrt(longedge);
	  if(longedge > precispt ) hmin = min( hmin, longedge);
	}
      }
    }
  }

  if(verbosity >5) cout << "    longmini_box" << longmini_box << endl;
  if(verbosity >5) cout << "    hmin =" << hmin << endl;
  assert( hmin < longmini_box);
  if(verbosity >5) cout << "    Norme2(bmin-bmax)=" <<  Norme2(bmin-bmax) << endl;
  // assertion pour la taille de l octree
  assert(hmin>Norme2(bmin-bmax)/1e9);
}

//======================
//
//======================
void OrderVertexTransfo_hcode_nv( const int &tab_nv, const double *tab_XX, const double *tab_YY, const double *tab_ZZ,  
				  const double *bmin, const double *bmax, const double hmin,int *Numero_Som, int * ind_nv_t, int & nv_t){
  size_t i;
  size_t j[3];
  size_t k[3];
  size_t NbCode = 100000;
  int *tcode; //= new int[NbCode];	

  int *posv = new int[tab_nv];
  double epsilon= hmin/10.;
  /*
    double epsilon=0.001;
			
    // determination de boite englobante
    double bmin[3],bmax[3];
	
    bmin[0] = tab_XX[0];
    bmin[1] = tab_YY[0];
    bmin[2] = tab_ZZ[0];
	
    bmax[0] = bmin[0];
    bmax[1] = bmin[1];
    bmax[2] = bmin[2];
	
    cout << " determination bmin et bmax" << endl;
	
    for(int ii=1; ii<tab_nv; ii++){
    bmin[0] = min(bmin[0],tab_XX[ii]);
    bmin[1] = min(bmin[1],tab_YY[ii]);
    bmin[2] = min(bmin[2],tab_ZZ[ii]);
 		
    bmax[0] = max(bmax[0],tab_XX[ii]);
    bmax[1] = max(bmax[1],tab_YY[ii]);
    bmax[2] = max(bmax[2],tab_ZZ[ii]);
    }
  */
	
  k[0] = int( (bmax[0]-bmin[0])/epsilon ); 
  k[1] = int( (bmax[1]-bmin[1])/epsilon ); 
  k[2] = int( (bmax[2]-bmin[2])/epsilon ); 
	
  int numberofpoints=0;
  int numberofpointsdiff;
	
  for(int ii=0; ii<tab_nv; ii++){
    numberofpointsdiff=0;
    for(int jj=ii+1; jj<tab_nv; jj++){
      double dist = 0.;
      dist = pow(tab_XX[jj]-tab_XX[ii],2)+pow(tab_YY[jj]-tab_YY[ii],2)+pow(tab_ZZ[jj]-tab_ZZ[ii],2); //pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
      if( sqrt(dist) < epsilon ){
	numberofpointsdiff=1;
      } 
    }
    if( numberofpointsdiff==0) numberofpoints=numberofpoints+1;
  }
	
  if(verbosity >4) cout << "   -- numberofpoints " << numberofpoints << endl;
  if(verbosity >4) cout << "   -- taille boite englobante =" << endl;
  if(verbosity >4)
    {
      for(int ii=0; ii<3; ii++){
	cout << "ii=" << ii << " " << bmin[ii] << " " << bmax[ii] <<endl;
      }
      
      for(int ii=0; ii<3; ii++){
	cout << "k[" << ii << "]= " << k[ii]<<endl;
      }
    }
  NbCode = min( 4*(k[0]+k[1]+k[2]), NbCode );
  tcode = new int[NbCode];
	
	
  /* initialisation des codes */
  for(int ii=0; ii< NbCode; ii++){
    tcode[ii] = -1;
  }
	
	
  for(int ii=0; ii < tab_nv; ii++){
    // boucle dans l autre sens pour assurer l'ordre des elements pour la suite
    //cout << "vertex ii " << ii << "  max : " << tab_nv;  
    j[0] = int( (tab_XX[ii]-bmin[0])/epsilon  );
    j[1] = int( (tab_YY[ii]-bmin[1])/epsilon  );
    j[2] = int( (tab_ZZ[ii]-bmin[2])/epsilon  );  
		
    assert( j[0] <=k[0] && j[0]>=0);
    assert( j[1] <=k[1] && j[1]>=0);
    assert( j[2] <=k[2] && j[2]>=0);
    i = (j[2]*(k[1]+1)+j[1]*(k[0]+1)+j[0]);
    //cout << i << endl;	
    i = i%NbCode;    
    assert( i < NbCode );
    posv[ii] = tcode[i];
    tcode[i] = ii;
  }
	
  if(verbosity >1) cout << " boucle numero de Sommet " << endl;
  for(int ii=0; ii<tab_nv; ii++){
    Numero_Som[ii]=-1;
  }
	
  if(verbosity >1) cout << " determinations des points confondus et numerotation " << endl;
	
  nv_t=0;
  for(int	icode =0; icode < NbCode; icode++){
    //int ii,jj;		
    double dist;		
		
    for(int ii=tcode[icode]; ii!=-1; ii=posv[ii]){
      if(Numero_Som[ii] != -1) continue; 
      Numero_Som[ii] = nv_t;				
      for(int jj=posv[ii]; jj!=-1; jj=posv[jj]){
	if(Numero_Som[jj] != -1) continue; 
	dist=pow(tab_XX[jj]-tab_XX[ii],2)+pow(tab_YY[jj]-tab_YY[ii],2)+pow(tab_ZZ[jj]-tab_ZZ[ii],2);	
			
	if( sqrt(dist) < epsilon ){
	  // point semblable
	  Numero_Som[jj] = Numero_Som[ii];
	  //cout << "point semblable" << endl;
	  //exit(-1); 
	}
				
      }
      ind_nv_t[nv_t] = ii; // Remarque on donne a nv_t le plus grand
      nv_t++; //nv_t = nvt+1;
    }
  }
  if(verbosity >1) cout << "      nv_t = " << nv_t << " / " << "nv_t(anc)" << tab_nv <<endl;
  assert( nv_t == numberofpoints);

  delete [] posv;
  delete [] tcode;
}

void PointCommun_hcode( const int &dim, const int &NbPoints, const int &point_confondus_ok, double **Coord_Point, 
			const double *bmin, const double *bmax, const double hmin, int * ind_np, int & np){
	
  size_t i;
  size_t j[dim];
  size_t k[dim];
  size_t NbCode = 100000;
  int *tcode; //= new int[NbCode];
  int *posv = new int[NbPoints];
  int *Numero_Som =new int[NbPoints];
	
  double epsilon=hmin/10.;
  /*
    double epsilon=0.0001;
    double bmin[dim],bmax[dim];
	
    for(int jj=0; jj<dim; jj++){ 
    bmin[jj] = Coord_Point[0][jj];
    bmax[jj] = bmin[jj];
    }
    for(int ii=1; ii<NbPoints; ii++){
    for(int jj=0; jj<dim; jj++){ 
    bmin[jj] = min(bmin[jj],Coord_Point[ii][jj]);
    bmax[jj] = max(bmax[jj],Coord_Point[ii][jj]);
    }
    }
  */
  assert( dim > 1);
	
  for(int jj=0; jj<dim; jj++){ 
    k[jj] = int( (bmax[jj]-bmin[jj])/epsilon );
  } 
	
  int numberofpoints=0;
  int numberofpointsdiff;
	
  for(int ii=0; ii<NbPoints; ii++){
    numberofpointsdiff=0;
    for(int jj=ii+1; jj<NbPoints; jj++){
      double dist = 0.;
      for( int kk=0; kk<3; kk++){
	dist = dist +  pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
      }
      if( sqrt(dist) < 1e-10){
	numberofpointsdiff=1;
      } 
    }
    if( numberofpointsdiff==0) numberofpoints=numberofpoints+1;
  }
	
  if(verbosity >1) cout << "numberofpoints " << numberofpoints << endl;
	
  NbCode = min( 4*(k[0]+k[1]+k[2]), NbCode);
  if(verbosity >1) cout << "NbCode=" << NbCode << endl;
  tcode = new int[NbCode];
  /* initialisation des codes */
  for(int ii=0; ii< NbCode; ii++){
    tcode[ii] = -1;
  }
	
  for(int ii=0; ii < NbPoints; ii++){
    // boucle dans l autre sens pour assurer l'ordre des elements pour la suite
		
    for( int jj=0; jj<dim; jj++){
      j[jj] = int( (Coord_Point[ii][jj]-bmin[jj])/epsilon  );
    }
	
    assert( j[0] <=k[0] && j[0]>=0);
    assert( j[1] <=k[1] && j[1]>=0);
    assert( j[2] <=k[2] && j[2]>=0);
		
    i = j[0];
    for(int jj=1; jj<dim; jj++){
      i=i+j[jj]*(k[jj-1]+1);
    }
    i = i%NbCode;    
		
    assert( i < NbCode );
    posv[ii] = tcode[i];
    tcode[i] = ii;
  }
  for(int ii=0; ii < NbPoints; ii++){
    ind_np[ii] = -1;
    Numero_Som[ii] = -1;
  }
	
  /* Resolution probleme dans le cas où le maillage se colle */
	
  /* maintenant determinations des points confondus et numerotation*/
	
  switch( point_confondus_ok ){
	
  case 0:
    np=0;
    for(int	icode =0; icode < NbCode; icode++){
      //int ii,jj;		
      double dist;		
		
      for(int ii=tcode[icode]; ii!=-1; ii=posv[ii]){
			
	if( Numero_Som[ii]!= -1 ) continue;
	//minimum_np=ii;
	Numero_Som[ii] = np;
			
	for(int jj=posv[ii]; jj!=-1; jj=posv[jj]){
				
	  if(Numero_Som[jj] != -1) continue; 
	  dist = 0.;
				
	  for( int kk=0; kk<dim; kk++){
	    dist = dist +  pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
	  }
					
	  if( sqrt(dist) < epsilon ){
	    // point semblable
	    Numero_Som[jj] = Numero_Som[ii];
	    //minimum_np = min( jj, minimum_np);
	  }
	}
	ind_np[np] = ii; //min(ii,minimum_np);	// Remarque on donne a np le plus petit element
	np++; //nv_t = nvt+1;
      }
    }
    break;
		
  case 1:
    int point_multiple;
    np=0;
    for(int icode =0; icode < NbCode; icode++){
      //int ii,jj;		
      double dist;		
		
      for(int ii=tcode[icode]; ii!=-1; ii=posv[ii]){
			
	if( Numero_Som[ii]!= -1 ) continue;										
	//minimum_np=ii;
	Numero_Som[ii] = np;
	point_multiple = 0; 
				
				
	for(int jj=posv[ii]; jj!=-1; jj=posv[jj]){
				
	  if(Numero_Som[jj] != -1) continue; 				
	  dist = 0.;
				
	  for( int kk=0; kk<dim; kk++){
	    dist = dist +  pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
	  }
							
	  if( sqrt(dist) < epsilon ){
	    // point semblable
	    Numero_Som[jj] = Numero_Som[ii];
	    point_multiple = 1; 
	    //minimum_np = min( jj, minimum_np);
	  }				
	}
	if(point_multiple ==0){
	  ind_np[np] = ii; //min(ii,minimum_np);	// Remarque on donne a np le plus petit element
	  np++; //nv_t = nvt+1;
	}
      }
    }
    break;
  default:
    cout << " point_confondus_ok dans fonction PointCommun_hcode vaut 1 ou 0." << endl;
    exit(-1);
  }

  delete [] tcode; 
  delete [] posv;
  delete [] Numero_Som;

}

void OrderVertexTransfo_hcode_nv_gtree( const int & tab_nv, const R3 &bmin, const R3 &bmax, const double &hmin,  
					const double *tab_XX, const double *tab_YY, const double *tab_ZZ, int *Numero_Som, int * ind_nv_t, int & nv_t){

  size_t i;
  size_t j[3];
  size_t k[3];
	
  // parametre interne pour debugger le code
  int verifnumberofpoints;
  verifnumberofpoints = 1;
	
  // hmin a determiner plus haut
  assert(hmin>Norme2(bmin-bmax)/1e9);
  double hseuil =hmin/10.; 

  //hseuil = hseuil/10.;
	
  Vertex3  *v= new Vertex3[tab_nv];
  //Vertex3  v[tab_nv];
  
  EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,bmin,bmax,0);
	
  if(verbosity >2){
    cout << "  -- taille de la boite " << endl;
    cout << "\t" <<  bmin.x << " " << bmin.y << " " << bmin.z <<  endl;
    cout << "\t"  << bmax.x << " " << bmax.y << " " << bmax.z <<  endl;
  }
	
  // creation of octree
  nv_t = 0;
  for (int ii=0;ii<tab_nv;ii++){
    const R3 r3vi( tab_XX[ii], tab_YY[ii],tab_ZZ[ii]);
    /*vi.x = tab_XX[ii];
      vi.y = tab_YY[ii];
      vi.z = tab_ZZ[ii];*/
    const Vertex3 &vi(r3vi);
    /*vi.x = tab_XX[ii];
      vi.y = tab_YY[ii];
      vi.z = tab_ZZ[ii];*/
	

    Vertex3 * pvi=gtree->ToClose(vi,hseuil);
    if(!pvi){
      v[nv_t].x = vi.x;
      v[nv_t].y = vi.y;
      v[nv_t].z = vi.z;
      v[nv_t].lab = vi.lab; // lab mis a zero par default
      ind_nv_t[nv_t] = ii;
      Numero_Som[ii] = nv_t;
      gtree->Add( v[nv_t] );
      nv_t=nv_t+1;
    }

    else{
      Numero_Som[ii] = pvi-v;
    }
  }

  delete gtree;
  delete [] v;
	
  if(verbosity >3) cout << "    hseuil=" << hseuil <<endl;
  if(verbosity >3) cout << "    nv_t = " << nv_t << " / " << "nv_t(anc)" << tab_nv <<endl;
  
  if(verifnumberofpoints ==1){
    int numberofpoints=0;
    int numberofpointsdiff;
	
    for(int ii=0; ii<tab_nv; ii++){
      numberofpointsdiff=0;
      for(int jj=ii+1; jj<tab_nv; jj++){
	double dist = 0.;
	dist = pow(tab_XX[jj]-tab_XX[ii],2)+pow(tab_YY[jj]-tab_YY[ii],2)+pow(tab_ZZ[jj]-tab_ZZ[ii],2);
	if( sqrt(dist) < hseuil){
	  numberofpointsdiff=1;
	} 
      }
      if( numberofpointsdiff==0) numberofpoints=numberofpoints+1;
    }
    if(verbosity >2) cout << "  -- numberofpoints " << numberofpoints << endl;
//    if(verbosity >2) cout << "  -- taille boite englobante =" << endl;
    //    assert(nv_t==numberofpoints);
  }
}

void PointCommun_hcode_gtree( const int &dim, const int &NbPoints, const int &point_confondus_ok, 
			      double **Coord_Point, const int * label_point, 
			      const R3 & bmin, const R3 & bmax, const double &hmin, int * ind_np, int * ind_label, int & np){
	
  double hseuil =hmin/10.; 
  Vertex3  *v= new Vertex3[NbPoints];
  //Vertex3 v[NbPoints];
  EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,bmin,bmax,0);
	
  if(verbosity >1) cout<< "verif hmin vertex3 GTree switch: " << point_confondus_ok << endl;
	
  int int_point_confondus_ok = point_confondus_ok;
	
  if(int_point_confondus_ok == 0){
    // accepte les points double 
    np = 0;
    for (int ii=0;ii<NbPoints;ii++){
      const R3 r3vi( Coord_Point[ii][0], Coord_Point[ii][1], Coord_Point[ii][2] ); 
      const Vertex3 &vi( r3vi);
      Vertex3 * pvi=gtree->ToClose(vi,hseuil);
			
      if(!pvi){
	v[np].x = vi.x;
	v[np].y = vi.y;
	v[np].z = vi.z;
	v[np].lab = vi.lab; // lab mis a zero par default
	ind_np[np] = ii;
	ind_label[np] = label_point[ii];
	gtree->Add( v[np++] );
      }
      else{
	ind_label[pvi-v] = min( ind_label[pvi-v], label_point[ii] );
      }
    }
    if(verbosity >1) cout << "np="<< np << endl;
		
  }
  if(int_point_confondus_ok == 1){
    // accepte les points double sont enleves
    np = 0;
    for (int ii=0;ii<NbPoints;ii++){
      const R3 r3vi( Coord_Point[ii][0], Coord_Point[ii][1], Coord_Point[ii][2] ); 
      //int label = label_point[ii];
      const Vertex3 &vi( r3vi);
      Vertex3 * pvi=gtree->ToClose(vi,hseuil);
			
      if(!pvi){
	v[np].x = vi.x;
	v[np].y = vi.y;
	v[np].z = vi.z;
	v[np].lab = vi.lab; // lab mis a zero par default
	ind_np[np] = ii;
	ind_label[np] = label_point[ii];
	gtree->Add( v[np++] );
      }
      else{
	ind_label[pvi-v] = min( ind_label[pvi-v], label_point[ii] );
      }
    }
		
    int ind_multiple[np];
		
    for(int ii=0; ii<np; ii++){
      ind_multiple[ii]=-1;
    }
    for (int ii=0;ii<NbPoints;ii++){
      const R3 r3vi(Coord_Point[ii][0], Coord_Point[ii][1], Coord_Point[ii][2] );
      //int label =  label_point[ii];
      const Vertex3 & vi(r3vi);
      Vertex3 * pvi=gtree->ToClose(vi,hseuil);
      ind_multiple[pvi-v]=ind_multiple[pvi-v]+1; 
    }
		
    int jnp;
    jnp=0;
    for(int ii=0; ii<np; ii++){
      if(ind_multiple[ii]==0){
	assert( jnp <= ii);
	ind_np[jnp] = ind_np[ii];
	ind_label[jnp] = ind_label[ii];
	jnp++;
      }
    }
    np=jnp;
  }
  if( int_point_confondus_ok != 0 && int_point_confondus_ok != 1){
    cout << " point_confondus_ok dans fonction PointCommun_hcode vaut 1 ou 0." << endl;
    exit(1);
  }
	
  delete gtree;
  delete [] v; 

  /*
    int z_verifnumberofpoints;
    z_verifnumberofpoints = 0;
    if(z_verifnumberofpoints ==1){
    int numberofpoints=0;
    int numberofpointsdiff;
    for(int ii=0; ii<NbPoints; ii++){
    numberofpointsdiff=0;
    for(int jj=ii+1; jj<NbPoints; jj++){
    double dist = 0.;
    for( int kk=0; kk<3; kk++){
    dist = dist +  pow(Coord_Point[jj][kk]-Coord_Point[ii][kk],2);
    }
    if( sqrt(dist) < hseuil/10){
    numberofpointsdiff=1;
    } 
    }
    if( numberofpointsdiff==0) numberofpoints=numberofpoints+1;
    if( point_confondus_ok==1 && numberofpointsdiff==1) numberofpoints=numberofpoints-1;
    }
    cout << "numberofpoints =" << numberofpoints<< endl;
    cout << "np =" << np<< endl;
    //assert( numberofpoints == np);
    }
  */
}

/* fin TransfoMesh_v2.cpp*/

/* debut buildlayer.cpp */
class BuildLayeMesh_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression enmax,ezmin,ezmax,xx,yy,zz;
  static const int n_name_param =9+4; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
  //int    arg(int i,Stack stack,int a ) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
public:
  BuildLayeMesh_Op(const basicAC_F0 &  args,Expression tth,Expression nmaxx) 
    : eTh(tth),enmax(nmaxx), ezmin(0),ezmax(0),xx(0),yy(0),zz(0)
  {
    if(verbosity >1) cout << "construction par BuilLayeMesh_Op" << endl;
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
    if( nargs[3] && nargs[9] ) 
	CompileError("uncompatible buildlayer (Th, region= , reftet=  ");
    if( nargs[4] && nargs[10] ) 
	CompileError("uncompatible buildlayer (Th, midlabel= , reffacemid=  ");
    if( nargs[5] && nargs[11] ) 
	CompileError("uncompatible buildlayer (Th, toplabel= , reffaceup=  ");
    if( nargs[6] && nargs[12] ) 
	CompileError("uncompatible buildlayer (Th, downlabel= , reffacelow=  ");
    
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type BuildLayeMesh_Op::name_param[]= {
  {  "zbound", &typeid(E_Array)},
  {  "transfo", &typeid(E_Array)},
  {  "coef", &typeid(double)},
  {  "reftet", &typeid(KN_<long>)}, // 3
  {  "reffacemid", &typeid(KN_<long> )},
  {  "reffaceup", &typeid(KN_<long> )},
  {  "reffacelow", &typeid(KN_<long> )},
  {  "facemerge", &typeid(long)},
  {  "ptmerge",&typeid(double)},
  {  "region", &typeid(KN_<long>)}, // 9
  {  "labelmid", &typeid(KN_<long> )},
  {  "labelup", &typeid(KN_<long> )},
  {  "labeldown", &typeid(KN_<long> )}, // 12
    
};


class  BuildLayerMesh : public OneOperator { public:  
    BuildLayerMesh() : OneOperator(atype<pmesh3>(),atype<pmesh>(),atype<long>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
    if(verbosity >1) cout << " je suis dans code(const basicAC_F0 & args) const" << endl; 
    //cout << "args: " << args << endl;   
    return  new BuildLayeMesh_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1])); 
  }
};

/* debut buildlayer.cpp */
class cubeMesh_Op : public E_F0mps
{
public:
    Expression nx,ny,nz;
    Expression xx,yy,zz;
    static const int n_name_param =3; //
    static basicAC_F0::name_and_type name_param[] ;
    Expression nargs[n_name_param];
    KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
    { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
    double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
    long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
    //int    arg(int i,Stack stack,int a ) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
public:
    cubeMesh_Op(const basicAC_F0 &  args,Expression nnx,Expression nny,Expression nnz,Expression transfo=0)
    : nx(nnx),ny(nny),nz(nnz),xx(0),yy(0),zz(0)
    {
        if(verbosity >1) cout << "construction par cubeMesh_Op" << endl;
        args.SetNameParam(n_name_param,name_param,nargs);
        const E_Array * a2 = dynamic_cast<const E_Array *>(transfo);
        int err =0;
        //cout << nargs[0] << " "<< a1 << endl;
        //cout << nargs[1] << " "<< a2 << endl;
         if(a2) {
            if(a2->size() !=3)
                CompileError("cube (n1,n2,n3, [X,Y,Z]) ");
            xx=to<double>( (*a2)[0]);
            yy=to<double>( (*a2)[1]);
            zz=to<double>( (*a2)[2]);
        }
        
    } 
    
    AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type cubeMesh_Op::name_param[]= {
    {  "region", &typeid(long )}, // 0
    {  "label", &typeid(KN_<long> )},
    {  "flags", &typeid(long )}
    
};


class  cubeMesh : public OneOperator { public:
    int xyz;
    cubeMesh() : OneOperator(atype<pmesh3>(),atype<long>(),atype<long>(),atype<long>()), xyz(0) {}
    cubeMesh(int ) : OneOperator(atype<pmesh3>(),atype<long>(),atype<long>(),atype<long>(),atype<E_Array>()) ,xyz(1) {}
    
    E_F0 * code(const basicAC_F0 & args) const
    {
        //cout << "args: " << args << endl;
        if(xyz)
            return  new cubeMesh_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]),t[2]->CastTo(args[2]),t[3]->CastTo(args[3]));
        else
            return  new cubeMesh_Op(args,t[0]->CastTo(args[0]),t[1]->CastTo(args[1]),t[2]->CastTo(args[2]));
    }
};
extern Mesh * Carre_(int nx,int ny,Expression fx,Expression fy,Stack stack,int flags,KN_<long> lab,long reg);

AnyType cubeMesh_Op::operator()(Stack stack)  const
{
   
    
    int n1 = (int) GetAny<long>((*nx)(stack));
    int n2 = (int) GetAny<long>((*ny)(stack));
    int nlayer = (int) GetAny<long>((*nz)(stack));
    MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
    long flags=arg(2,stack,0L);;
    KN<long> label2;
    long ll3[]={1,2,3,4,5,6};
    KN_<long>l3_(ll3,6);
    KN<long> l3=arg(1,stack,l3_);
    ffassert(l3.N()==6);
    long region=arg(0,stack,0L);
 
    Mesh * pTh=  Carre_( n1 , n2,0,0,stack,flags,label2,0L );// WARNING no clean of teh stack in this case (durdur)
    ffassert(pTh && nlayer>0);
    Mesh &Th=*pTh;
    int nbv=Th.nv; // nombre de sommet
    int nbt=Th.nt; // nombre de triangles
    int neb=Th.neb; // nombre d'aretes fontiere
    if(verbosity>2)
        cout << "  -- cubeMesh_Op input: nv" << nbv<< "  nt: "<< nbt << " nbe "<< neb << endl;
    KN<double> zmin(nbv),zmax(nbv);
    KN<double> clayer(nbv); //  nombre de layer est nlayer*clayer
    
    clayer=1.;
    zmin=0.;
    zmax=1.;
    double maxdz = 1.;
    
    if(verbosity >3) cout << "lecture valeur des references " << endl;
    long rup[]={0,l3[6-1]}, rdown[]={0,l3[5-1]}, r0[]={0,region};
    long rmid[]={1,l3[1-1],2,l3[2-1],3,l3[3-1],4,l3[4-1]};
  //  KN_<long> zzempty;
    KN_<long> nrtet(r0,2);
    KN_<long> nrfmid(rmid,8);
    KN_<long> nrfup(rup,2);
    KN_<long> nrfdown(rdown,2);
    
    int point_confondus_ok (0L);
    double precis_mesh (-1L);
    
    
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
    double epsz = maxdz *1e-6;
    if(verbosity>9999)    cout << ":: epsz " << epsz <<endl;
    for(int i=0;i<nbv;i++)
    {
        ni[i]=Max(0,Min(nlayer,(int) lrint(nlayer*clayer[i])));
        if(abs(zmin[i]-zmax[i]) < epsz) ni[i]=0; // Corr FH aug. 2014...
    }
    if(verbosity>9999)    cout << " cubeMesh_Op: ni = " << ni << endl;
    // triangle
    for (int it=0;it<nbt;++it){
        const Mesh::Element &K(Th.t(it));
        int i0 = Th.operator()(K[0]);
        int i1 = Th.operator()(K[1]);
        int i2 = Th.operator()(K[2]);
        
        if( ni[i0] == 0 && ni[i1] == 0 && ni[i2] == 0 ){
            cout << "A tetrahedra with null volume will be created with triangle " << it << " of 2D Mesh " << endl;
            cout << "stop procedure of buildlayer" << endl;
            exit(1);
        }
        
    }
    
    // cas maillage volumique + surfacique
    Mesh3 *Th3= build_layer(Th, nlayer, ni, zmin, zmax, maptet, maptrimil, maptrizmax, maptrizmin, mapemil, mapezmax, mapezmin);
    // cas maillage surfacique simplement // A construire Jacques + donner le numero des edges que l'on veut pas creer � l'int�rieure
    
    delete pTh;
    
   
    if( xx && yy && zz )
    {
        //Mesh3 *Th3= build_layer(Th, nlayer, ni, zmin, zmax);
        
        KN<double> txx(Th3->nv), tyy(Th3->nv), tzz(Th3->nv);
        KN<int> takemesh(Th3->nv);
        //MeshPoint *mp3(MeshPointStack(stack));
        
        takemesh=0;
        for (int it=0;it<Th3->nt;++it){
            for( int iv=0; iv<4; ++iv){
                int i=(*Th3)(it,iv);
                if(takemesh[i]==0){
                    mp->setP(Th3,it,iv);
                    { txx[i]=GetAny<double>((*xx)(stack));}
                    { tyy[i]=GetAny<double>((*yy)(stack));}
                    { tzz[i]=GetAny<double>((*zz)(stack));}
                    takemesh[i] = takemesh[i]+1;
                }
            }
        }
        
        int border_only = 0;
        int recollement_elem=0, recollement_border=1; 
        if(point_confondus_ok == 2){
            recollement_border = 0;
            point_confondus_ok = 1;
        }
        
        Mesh3 *T_Th3=Transfo_Mesh3( precis_mesh, *Th3, txx, tyy, tzz, border_only, recollement_elem, recollement_border, point_confondus_ok,1);
        delete Th3;
        Th3 =T_Th3;
        
    }
    Th3->BuildGTree(); //A decommenter
    if(verbosity >10) cout << " Cube %%% " << Th3 << endl;
    Add2StackOfPtr2FreeRC(stack,Th3);
    *mp=mps;
    return Th3;
    
}

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
  if(verbosity>2)
  cout << "  -- BuildLayeMesh_Op input: nv" << nbv<< "  nt: "<< nbt << " nbe "<< neb << endl; 
  KN<double> zmin(nbv),zmax(nbv);
  KN<double> clayer(nbv); //  nombre de layer est nlayer*clayer
  
  clayer=-1;
  zmin=0.;
  zmax=1.;
    double maxdz = 0;   
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
	  maxdz = max(maxdz, abs(zmin[i]-zmax[i]));
	  clayer[i]=Max( 0. , Min( 1. , arg(2,stack,1.) ) ); 
	
	}	     
    }
  }
  
  ffassert(clayer.min() >=0);

  if(verbosity >1) cout << "lecture valeur des references " << endl;
  
  KN<long> zzempty;
  KN<long> nrtet  (arg(3,stack,arg(3+6,stack,zzempty)));  
  KN<long> nrfmid (arg(4,stack,arg(4+6,stack,zzempty)));  
  KN<long> nrfup  (arg(5,stack,arg(5+6,stack,zzempty)));  
  KN<long> nrfdown (arg(6,stack,arg(6+6,stack,zzempty)));  
  int point_confondus_ok (arg(7,stack,0L));
  double precis_mesh (arg(8,stack,-1L));

  
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
    double epsz = maxdz *1e-6;
 if(verbosity>9999)    cout << "BuildLayeMesh_Op:: epsz " << epsz <<endl; 
  for(int i=0;i<nbv;i++)
    {
      ni[i]=Max(0,Min(nlayer,(int) lrint(nlayer*clayer[i])));
      if(abs(zmin[i]-zmax[i]) < epsz) ni[i]=0; // Corr FH aug. 2014...
    }
  if(verbosity>9999)    cout << " BuildLayeMesh_Op: ni = " << ni << endl; 
  // triangle 
  for (int it=0;it<nbt;++it){
    const Mesh::Element &K(Th.t(it));
    int i0 = Th.operator()(K[0]); 
    int i1 = Th.operator()(K[1]); 
    int i2 = Th.operator()(K[2]); 
    
    if( ni[i0] == 0 && ni[i1] == 0 && ni[i2] == 0 ){
      cout << "A tetrahedra with null volume will be created with triangle " << it << " of 2D Mesh " << endl;
      cout << "stop procedure of buildlayer" << endl;
      exit(1);
    }

  }

  // cas maillage volumique + surfacique
  Mesh3 *Th3= build_layer(Th, nlayer, ni, zmin, zmax, maptet, maptrimil, maptrizmax, maptrizmin, mapemil, mapezmax, mapezmin);
  // cas maillage surfacique simplement // A construire Jacques + donner le numero des edges que l'on veut pas creer à l'intérieure
    
  
  if( !(xx) && !(yy) && !(zz) )
    {
	 /*
	  map< int, int > maptet; 
	  map< int, int > maptrimil, maptrizmax, maptrizmin;
	  map< int, int > mapemil, mapezmax, mapezmin;
	   
	  build_layer_map_tetrahedra( Th, maptet );
	  build_layer_map_triangle( Th, maptrimil, maptrizmax, maptrizmin );
	  build_layer_map_edge( Th, mapemil, mapezmax, mapezmin );
	  
	  Mesh3 *Th3= build_layer(Th, nlayer, ni, zmin, zmax, maptet, maptrimil, maptrizmax, maptrizmin, mapemil, mapezmax, mapezmin);
	  */
      
     // Th3->BuildBound();
     // Th3->BuildAdj();
     // Th3->Buildbnormalv();  
     // Th3->BuildjElementConteningVertex();
      
      Th3->BuildGTree(); //A decommenter

      Add2StackOfPtr2FreeRC(stack,Th3);
      *mp=mps;
      return Th3;
    }
  else
    {  
      //Mesh3 *Th3= build_layer(Th, nlayer, ni, zmin, zmax);
	  
      KN<double> txx(Th3->nv), tyy(Th3->nv), tzz(Th3->nv);
      KN<int> takemesh(Th3->nv);
      //MeshPoint *mp3(MeshPointStack(stack)); 
	  
      takemesh=0;  
      Mesh3 &rTh3 = *Th3;
      for (int it=0;it<Th3->nt;++it){
	for( int iv=0; iv<4; ++iv){
	  int i=(*Th3)(it,iv);  
	  if(takemesh[i]==0){
	    mp->setP(Th3,it,iv);
	    if(xx){ txx[i]=GetAny<double>((*xx)(stack));}
	    if(yy){ tyy[i]=GetAny<double>((*yy)(stack));}
	    if(zz){ tzz[i]=GetAny<double>((*zz)(stack));}

	    takemesh[i] = takemesh[i]+1;
	  }
	}
      }
		
      int border_only = 0;
      int recollement_elem=0, recollement_border=1; 
      if(point_confondus_ok == 2){
	recollement_border = 0;
	point_confondus_ok = 1;
      }

        Mesh3 *T_Th3=Transfo_Mesh3( precis_mesh, rTh3, txx, tyy, tzz, border_only, recollement_elem, recollement_border, point_confondus_ok,1);
		  
      
      // T_Th3->BuildBound();
      //  T_Th3->BuildAdj();
      // T_Th3->Buildbnormalv();  
      // T_Th3->BuildjElementConteningVertex();
      
      
      T_Th3->BuildGTree(); //A decommenter
      
      delete Th3;
      Add2StackOfPtr2FreeRC(stack,T_Th3);
      *mp=mps;
      return T_Th3;

    }
}


// function nouveau nom de fonction

class Movemesh2D_3D_surf_cout_Op : public E_F0mps 
{
public:
  Movemesh2D_3D_surf_cout_Op(const basicAC_F0 &  args,Expression tth)  
  {
    CompileError("The keyword movemesh2D3Dsurf is remplaced now by the keyword movemesh23 (see Manual) ::: Moreover, the parameter mesuremesh are denoted now orientation ");
  } 
   AnyType operator()(Stack stack)const { return 0L; }
  
};




class Movemesh2D_3D_surf_cout : public OneOperator { public:  
typedef const Mesh *pmesh;
typedef const Mesh3 *pmesh3;
    
  Movemesh2D_3D_surf_cout() : OneOperator(atype<pmesh3>(),atype<pmesh>() ) {}
    E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new Movemesh2D_3D_surf_cout_Op(args,t[0]->CastTo(args[0]));  // CastTo(args[]); // plus tard
  }
};


/***********************************************/

class Movemesh3D_cout_Op : public E_F0mps 
{
 
public:
  Movemesh3D_cout_Op(const basicAC_F0 &  args,Expression tth)  
 
  {
    CompileError("The keyword movemesh3D is remplaced in this new version of freefem++ by the keyword movemesh3 (see manual)");
  } 
  
  AnyType operator()(Stack stack)  const {return 0L;}
  
};




class Movemesh3D_cout : public OneOperator { public:  
typedef const Mesh *pmesh;
typedef const Mesh3 *pmesh3;
    
  Movemesh3D_cout() : OneOperator(atype<pmesh3>(),atype<pmesh>() ) {}
    E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new Movemesh3D_cout_Op(args,t[0]->CastTo(args[0]));  // CastTo(args[]); // plus tard
  }
};


//

class DeplacementTab_Op : public E_F0mps 
{
public:
  Expression eTh;
  //Expression xx,yy,zz;
  static const int n_name_param =6; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long  arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
public:
  DeplacementTab_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth)  //, xx(0) , yy(0) , zz(0)
  {
    args.SetNameParam(n_name_param,name_param,nargs);
   
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type DeplacementTab_Op::name_param[]= {
  {  "deltax", &typeid(KN_<double>)},
  {  "deltay", &typeid(KN_<double>)},
  {  "deltaz", &typeid(KN_<double>)},
  {  "ptmerge", &typeid(double)},
  {  "facemerge", &typeid(long)},
  {  "boolsurface",&typeid(long)}
};

AnyType DeplacementTab_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));
  
  ffassert(pTh);
  Mesh3 &Th=*pTh;
  Mesh3 *m= pTh;   // question a quoi sert *m ??
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.nbe; // nombre d'aretes fontiere
  if(verbosity>5)
  cout << "before movemesh: Vertex " << nbv<< " Tetrahedra " << nbt << " triangles "<< nbe << endl; 
 
 // lecture des references
   
  KN<double> zzempty;
  KN<double> dx (arg(0,stack,zzempty));
  KN<double> dy (arg(1,stack,zzempty));
  KN<double> dz (arg(2,stack,zzempty));
  double precis_mesh( arg(3,stack,1e-7));

  ffassert( dx.N() == Th.nv);
  ffassert( dy.N() == Th.nv);
  ffassert( dz.N() == Th.nv);
    
  // realisation de la map par default
 
  KN<double> txx(Th.nv), tyy(Th.nv), tzz(Th.nv);
  // loop over tetrahedral 
  for (int i=0;i<Th.nv;++i){   
    txx[i]=Th.vertices[i].x+dx[i];
    tyy[i]=Th.vertices[i].y+dy[i];
    tzz[i]=Th.vertices[i].z+dz[i];
  }
  
  int border_only = 0;
  int recollement_elem=0;
  int recollement_border,point_confondus_ok;
  
  int  mergefacemesh( arg(4,stack,0L) );
  long  flagsurfaceall( arg(5,stack,1L) );

   if(mergefacemesh == 0)
    {
      recollement_border=0;
      point_confondus_ok=0;
    }
  if(mergefacemesh == 1)
    { 
      recollement_border=1;
      point_confondus_ok=0;
    }
  if(mergefacemesh == 2)
    { 
      recollement_border=1;
      point_confondus_ok=1;
    }

  Mesh3 *T_Th3=Transfo_Mesh3( precis_mesh,Th, txx, tyy, tzz, border_only, 
			      recollement_elem, recollement_border, point_confondus_ok,1);
  
  if(nbt != 0)
    {
      //T_Th3->BuildBound();
    
      //T_Th3->BuildAdj();

      if(flagsurfaceall==1) T_Th3->BuildBoundaryElementAdj();
      
      //T_Th3->Buildbnormalv();  

     // T_Th3->BuildjElementConteningVertex();
      
      T_Th3->BuildGTree();
      
      //	T_Th3->decrement(); 
    }
  else
    {
       if(flagsurfaceall==1) T_Th3->BuildBoundaryElementAdj();
    }
  Add2StackOfPtr2FreeRC(stack,T_Th3);
 
  *mp=mps;
  return T_Th3;
}

class  DeplacementTab : public OneOperator { public:  
    DeplacementTab() : OneOperator(atype<pmesh3>(),atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	return  new DeplacementTab_Op(args,t[0]->CastTo(args[0])); 
  }

};

// CheckSurfaceMesh ==> 

int  GetBEManifold( Expression bb,  Expression &label, Expression &orient);
void GetNumberBEManifold( Expression surf, int & mani_nbe);

void GetManifolds( Expression mani, int & nbcmanifold,  int * &mani_nbe, Expression * &manifold)
{
  if ( mani ) 
    {
      int i,j;
      const E_Array * a= dynamic_cast<const  E_Array *>(mani);
      ffassert(a);
      int n = a->size();
      if( verbosity>1) 
	cout << "    the number of manifold " << n << endl;
      
      nbcmanifold = n;  // nombre de manifold définis
      
      //manifold = new Expression[n]; 
      mani_nbe = new int[n];
      int size = 0;
      for ( i=0; i<n; i++){
	GetNumberBEManifold( (*a)[i], mani_nbe[i]);
	cout << "number of manifold = " << n << "manifold i=" << i << "nb BE label=" << mani_nbe[i] << endl;
	size=size+mani_nbe[i];
      }
   
      manifold = new Expression[size*2];
      int count=0;
      for ( i=0; i<n; i++){
	Expression tmp=(*a)[i];
	const E_Array * aa = dynamic_cast<const  E_Array *>( tmp );
	for( j=0; j< mani_nbe[i]; j++){
	  if(GetBEManifold( (*aa)[j], manifold[count], manifold[count+1] ) ==0)   
	    CompileError(" a manifold is defined by a pair of [label, orientation ]");
	  count=count+2;
	}
      }
      assert(count == 2*size);
    }  
}

void GetNumberBEManifold( Expression surf, int & mani_nbe)
{
  if ( surf ) 
    {
      int i,j;
      if( verbosity>1) 
	cout << "  -- Manifoldal Condition to do" << endl;
      const E_Array * a= dynamic_cast<const  E_Array *>(surf);
      ffassert(a);
      mani_nbe = a->size();
      
    }
}


// void GetManifold( Expression surf, int & mani_nbe, Expression * &manifold)
// {
//   if ( surf ) 
//     {
//       int i,j;
//       if( verbosity>1) 
// 	cout << "  -- Manifoldal Condition to do" << endl;
//       const E_Array * a= dynamic_cast<const  E_Array *>(surf);
//       ffassert(a);
//       int n = a->size()/2;
//       mani_nbe = n;
//       if( verbosity>1) 
// 	cout << "    the number of face label in a manifold " << n << endl;
//       if( n*2 != a->size() )
// 	CompileError(" a manifold is defined by a pair of [label, orientation ]");
//       manifold = new Expression[n*2]; 
//       for ( i=0,j=0;i<n;i++,j+=2)
// 	if (GetBEManifold((*a)[i],manifold[j],manifold[j+1])==0)
// 	  CompileError(" a sub array of a sub manifold must be [label, orientation ]");
//     }
// }

int GetBEManifold( Expression bb, Expression &label, Expression &orient)
{
  
  const E_Array * a= dynamic_cast<const  E_Array *>(bb);
  if( a &&  a->size() == 2 ){
    label  =  to<long>((*a)[0]);
    orient =  to<long>((*a)[1]);
    
    return 1;
  }
  else
    return 0;
}

class CheckManifoldMesh_Op : public E_F0mps 
{
public:
  Expression eTh;
  static const int n_name_param =1; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  int nbmanifold;
  int *mani_nbe;
  Expression *manifolds;

  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long  arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
public:
  CheckManifoldMesh_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth)
  {
    args.SetNameParam(n_name_param,name_param,nargs);
    if(nargs[0])   
      GetManifolds(nargs[0],nbmanifold, mani_nbe, manifolds);
    else
      CompileError("check ::: no definition of manifold");
   
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type CheckManifoldMesh_Op::name_param[]= {
  {  "manifolds", &typeid(E_Array)}
  // option a rajouter
  // facemerge 0,1 + label
};

AnyType CheckManifoldMesh_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh= GetAny<Mesh3 *>((*eTh)(stack));


  int size=0;
  KN<int> BeginManifold(nbmanifold+1);
  
  for (int i=0; i< nbmanifold; i++){
    BeginManifold[i]=size;
    size=size+mani_nbe[i];
  }
  BeginManifold[nbmanifold]=size;
  
  KN<int> TabLabelManifold(size ), OrientLabelManifold(size );
    
  int count=0;
  for (int i=0; i< nbmanifold; i++){
    for(int j=0; j< mani_nbe[i]; j++){
      TabLabelManifold   [count]  = GetAny< long > ( ( *manifolds[2*count] )(stack) ); 
      OrientLabelManifold [count] = GetAny< long > ( ( *manifolds[2*count+1] )(stack) ); 
      count++;
    }
  }
  assert(count == size);

//   int count=0;
//   for(int ii=0; ii<nbvariete; ii++)
//     {
//       beginvariete[ii]=count;
//       for(int jj=0; jj<labelvariete[ii]; jj++)
// 	{
// 	  TabLabelVariete    [count] = GetAny< long > ( (*surface[2*ii+1][2*jj])(stack) ); 
// 	  OrientLabelVariete [count] = GetAny< long > ( (*surface[2*ii+1][2*jj+1])(stack) ); 
// 	  count++;
// 	}
//     }
//   beginvariete[nbvariete]=count;
  
  long resultat=1;
  pTh->BuildBoundaryElementAdj( nbmanifold, BeginManifold,TabLabelManifold,OrientLabelManifold); // nbvariete, beginvariete, TabLabelVariete, OrientLabelVariete);
  
  cout << "utilisation V2" << endl;
  *mp=mps;
  return resultat;
}
      

class  CheckManifoldMesh : public OneOperator { public:  
    CheckManifoldMesh() : OneOperator(atype<long>(),atype<pmesh3>() ) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
	return  new CheckManifoldMesh_Op(args,t[0]->CastTo(args[0])); 
  }
};




Mesh3 * truncmesh(const Mesh3 &Th,const long &kksplit,int *split, bool kk, const int newbelabel);

struct Op_trunc_mesh3 : public OneOperator {
  typedef const Mesh3 *pmesh3;
  class Op: public E_F0mps   { 
  public:
    static basicAC_F0::name_and_type name_param[] ;
    static const int n_name_param =5;
    Expression nargs[n_name_param];
    
    Expression getmesh,bbb;
    long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
    bool arg(int i,Stack stack,bool a) const{ return nargs[i] ? GetAny<bool>( (*nargs[i])(stack) ): a;}
      KN<long> *  arg(int i,Stack stack) const{ return nargs[i] ? GetAny<KN<long> *>( (*nargs[i])(stack) ): 0;}
      
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
   {  "label",             &typeid(long)},
     { "new2old", &typeid(KN<long>*)},  //  ajout FH pour P. Jovilet jan 2014
     { "old2new", &typeid(KN<long>*)},   //  ajout FH pour P. Jovilet jan 2014
     { "renum",&typeid(bool)}
     
 };


Mesh3 * truncmesh(const Mesh3 &Th,const long &kksplit,int *split, bool kk, const int newbelabel)
{
    
    static const int FaceTriangle[4]={3,0,1,2};  //={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}
    
    
    // computation of number of border elements and vertex without split
    int nbe = 0;
    int nt  = 0;
    int nv  = 0;
    int nvtrunc =0;
    int nedge=0;
    int nface=0;
    double hmin=1e100;
    R3 bmin,bmax;
    int nbeee=0,nbfi=0;
    const int kksplit2 = kksplit*kksplit;
    const int kksplit3 = kksplit2*kksplit;
    int ntsplit =0;
    int tagb[4]={1,2,4,8} ;
    KN<int> tagTonB(Th.nt);
    tagTonB=0;
    
    for( int ibe=0; ibe < Th.nbe; ibe++)
    {
        int iff;
        int it=Th.BoundaryElement(ibe,iff);
        tagTonB[it]|= tagb[iff];
        int ifff=iff,itt=Th.ElementAdj(it,ifff);
        if(itt >=0 &&  itt != it)
            tagTonB[itt]|= tagb[ifff];
    }
    
    for (int i=0;i<Th.nt;i++)
        if(split[i])
        {
            ++ntsplit;
            // computation of number of tetrahedrons
            nt=nt+kksplit3;
            // computation of number of border elements
            for (int j=0;j<4;j++)
            {
                int jt=j,it=Th.ElementAdj(i,jt);
                if ( (it==i || it <0)  || ! split[it]) nbeee++;// boundary face ...
                else nbfi++; // internal face count 2 times ...
                if(it==i || it <0) nbe += kksplit2;  //on est sur la frontiere
                else if (!split[it]) nbe += kksplit2; //le voisin ne doit pas etre decoupe
                else if ( (tagTonB[i]&tagb[j] ) != 0 && i<it) nbe += kksplit2; // internal boundary ..
            }
            
            for (int e=0;e<6;e++){
                hmin=min(hmin,Th[i].lenEdge(e));   // calcul de .lenEdge pour un Mesh3
            }
        }
    ffassert( nbfi %2 ==0) ;
    nface = nbeee + nbfi/2;
    double hseuil = (hmin/kksplit)/1000.;
    if(verbosity>5)
        cout << "  number of  not intern boundary faces = " << nbeee << ",  all faces  =  " << nbe << ", hseuil=" << hseuil <<endl;
    
    /* determination de bmin, bmax et hmin */
    
    KN<int> takevertex(Th.nv,-1);
    
    for(int i=0; i<Th.nt; i++){
        if(split[i])
        {
            const Tet &K(Th.elements[i]);
            for(int ii=0; ii<4; ii++)
            {
                int iv= Th.operator()( K[ii] );
                if( takevertex[iv] == -1 )
                {
                    bmin=Minc(Th.vertices[iv],bmin);
                    bmax=Maxc(Th.vertices[iv],bmax);
                    takevertex[iv]=nvtrunc++;
                }
            }
            
        }
    }
    //  take same numbering
    for(int i=0,k=0; i<Th.nv; i++)
        if(takevertex[i]>=0) takevertex[i]=k++;
    
    if( kksplit > 1 )
    { // compute the number of slip edge ...
        nedge=0;
        HashTable<SortArray<int,2>,int> edges(3*nface,nface);
        for(int i=0; i<Th.nt; i++){
            if(split[i])
            {
                const Tet &K(Th.elements[i]);
                for(int e=0;e<6;++e)
                {
                    
                    int e1 = Th( K[ Th[i].nvedge[e][0] ] );
                    int e2 = Th( K[ Th[i].nvedge[e][1] ] );
                    SortArray<int,2> key(e1,e2);
                    if(!edges.find(key) )
                        edges.add(key,nedge++);
                }
            }
        }
    }
    if(verbosity>10) cout << "    -- nvertex  " << nvtrunc << ", nedges = "<< nedge
        << ", nfaces = " << nface << " ntet =" << ntsplit
        << endl
        << "    -- Euler/Poincare constante = " << nvtrunc-nedge+nface-ntsplit
        << endl;
    
    /* determination des vertex, triangles et tetrahedre obtenue apres splitting dans le Simplex */
    
    int nfacesub = kksplit2;
    int ntetsub = kksplit3;
    int nvsub = (kksplit+1)*(kksplit+2)*(kksplit+3)/6;
    int ntrisub = 4*kksplit2;
    
    R3 *vertexsub; //[nvsub];
    int *tetsub;   //[4*ntetsub];
    int *trisub;   //[4*kksplit*kksplit];
    
    SplitSimplex<R3>( kksplit, nvsub, vertexsub, ntetsub, tetsub);
    SplitSurfaceSimplex( kksplit, ntrisub, trisub);
    
    if(verbosity>3)
        cout << "  -- trunc (3d) : Th.nv= " << Th.nv << "kksplit="<< kksplit << endl;
    
    int ntnosplit  = nt/kksplit3;
    int nbenosplit = nbe/kksplit2;
    int nfacenosplit = (4*ntnosplit+nbenosplit)/2;
    nv = ntnosplit*(nvsub - 4*( (kksplit+1)*(kksplit+2)/2 - 3*(kksplit-1) -3 ) - 6*( kksplit-1 ) - 4);
    if(verbosity>100) cout << "       1) nv= " << nv << endl;
    nv = nv + nfacenosplit*( (kksplit+1)*(kksplit+2)/2 - 3*(kksplit-1) -3 );
    if(verbosity>100) cout << "       2) nv= " << nv << endl;
    nv = nv + nedge*( kksplit-1 );
    if(verbosity>100) cout << "       3) nv= " << nv << endl;
    nv = nv + nvtrunc;
    if(verbosity>100) cout << "       4) nv= " << nv << endl;
    
    
    
    int itt=0;
    int ie=0;
    
    
    Vertex3 *v=new Vertex3[nv];
    Tet *t = new Tet[nt];
    Tet *tt = t;
    
    Triangle3 *b  = new Triangle3[nbe];
    Triangle3 *bb = b;
    R3 hh = (bmax-bmax)/10.;
    EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,bmin-hh,bmax+hh,0);
    const R3 * pP[4];
    int np=0; // nb of new points ..
    // first build old point to keep the numbering order for DDM ...
    for(int i=0,k=0; i<Th.nv; i++)
        if(takevertex[i]>=0)
        {
            Vertex3 * pvi=gtree->ToClose(Th(i),hseuil);
            if(!pvi)
            {
                (R3&) v[np]= Th(i);
                v[np].lab=Th(i).lab;
                gtree->Add( v[np] );
                np++;
            }
            else ffassert(0);
        }

    {
        KN<R3>  vertextetsub(nvsub);
        KN<int> newindex (nvsub);
        
        for(int i=0; i<Th.nt; i++)
            if(split[i])
            {
                const Tet &K(Th.elements[i]);
                
                for(int ii=0; ii< 4; ii++)
                    pP[ii] = & K[ii];
                
                for( int iv=0; iv<nvsub; iv++)
                    (R3&) vertextetsub[iv]= vertexsub[iv].Bary(pP);
                
                for( int iv=0; iv<nvsub; iv++)
                {
                    Vertex3 * pvi=gtree->ToClose(vertextetsub[iv],hseuil);
                    
                    if(!pvi)
                    {
                        (R3&) v[np]   = vertextetsub[iv];
                        v[np].lab = K.lab;
                        newindex[iv] = np;
                        gtree->Add( v[np] );
                        np++;
                    }
                    else
                        newindex[iv] = pvi-v;
                    
                    ffassert( np <= nv );
                }
                
                for( int ii=0; ii<ntetsub; ii++)
                {
                    int ivt[4];
                    for( int jj=0; jj< 4; jj++)
                    {
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
                    
                    if ( ( (tagTonB[i]&tagb[j]) ==0 ) &&  !(it==i || it <0)  && !split[it])
                    {
                        // new border not on boundary
                        int ivb[3];
                        
                        for( int ii=0; ii<nfacesub; ii++)
                        {
                            int iface = 3*FaceTriangle[j]*nfacesub+3*ii;
                            
                            for( int jjj=0; jjj<3; jjj++)
                            {
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
    if(verbosity>10)
        cout  << "    ++ np=" << np << "==  nv=" << nv << endl;
    ffassert( np == nv); 
    if(verbosity>8)
        cout << "   -- Number of new  border face not on Border " << ie << endl;
    delete [] vertexsub; //[nvsub];
    delete [] tetsub;   //[4*ntetsub];
    delete [] trisub;   //[4*kksplit*kksplit];
    
    // split border elements
    int nv2Dsub   = (kksplit+1)*(kksplit+2)/4;
    int ntri2Dsub = kksplit2;
    R2 *vertex2Dsub; //[nvsub];
    int *tri2Dsub;   //[4*kksplit*kksplit];
    
    SplitSimplex<R2>( kksplit, nv2Dsub, vertex2Dsub, ntri2Dsub, tri2Dsub);
    
    
    for( int ibe=0; ibe < Th.nbe; ibe++)
    {
        int iff;
        int it=Th.BoundaryElement(ibe,iff);
        int ifff=iff,itt=Th.ElementAdj(it,ifff);
        if(itt<0) itt=it;
        if( split[it] == 0 && split[itt] == 0) continue; // boundary not on one element
        
        const Triangle3 &K(Th.be(ibe));
        int ivv[3];
        
        ivv[0] = Th.operator()(K[0]);
        ivv[1] = Th.operator()(K[1]);
        ivv[2] = Th.operator()(K[2]);
        
        R3 *vertextrisub = new R3  [nv2Dsub];
        int *newindex = new int[nv2Dsub];
        for( int iv=0; iv<nv2Dsub; iv++)
        {
            double alpha=vertex2Dsub[iv].x;
            double beta=vertex2Dsub[iv].y;
            
            vertextrisub[iv].x = (1-alpha-beta)*Th.vertices[ivv[0]].x + alpha*Th.vertices[ivv[1]].x + beta*Th.vertices[ivv[2]].x;
            vertextrisub[iv].y = (1-alpha-beta)*Th.vertices[ivv[0]].y + alpha*Th.vertices[ivv[1]].y + beta*Th.vertices[ivv[2]].y;
            vertextrisub[iv].z = (1-alpha-beta)*Th.vertices[ivv[0]].z + alpha*Th.vertices[ivv[1]].z + beta*Th.vertices[ivv[2]].z;
            
        }
        
        for( int iv=0; iv<nv2Dsub; iv++)
        {
            const Vertex3 &vi( vertextrisub[iv] );
            Vertex3 * pvi=gtree->ToClose(vi,hseuil);
            assert(pvi);
            newindex[iv] = pvi-v;
        }
        
        for( int ii=0; ii<nfacesub; ii++)
        {
            int ivb[3];
            for( int jjj=0; jjj<3; jjj++)
            {
                ivb[jjj] = newindex[ tri2Dsub[3*ii+jjj] ];
                assert( tri2Dsub[ 3*ii+jjj  ] < nvsub );
                if(verbosity > 199 ) cout << "        " << ivb[jjj] << " np:" << np<< endl;
                assert( ivb[jjj] < np );
            }
            
            (bb++)->set( v, ivb, K.lab);
            ie++;
            assert( ie <= nbe);
        }
        delete [] vertextrisub;
        delete [] newindex;
        
        
    }
    
    delete [] vertex2Dsub;   //[4*ntetsub];
    delete [] tri2Dsub;   //[4*kksplit*kksplit];
    
    
    if(verbosity>99)
    {
        cout << "nbofv initial" << Th.nv << endl;
        cout << "nv=" << nv << " np=" << np << endl;
        cout << "itt=" << itt << " nt=" << nt << endl;
        cout << "ie=" << ie << " nbe=" << nbe << endl;
    }
    ffassert( nv == np );
    ffassert( ie ==nbe);
    ffassert( itt == nt );
    
    //delete gtree;
    
    Mesh3 *Tht = new Mesh3( nv, nt, nbe, v, t, b);
    Tht->BuildGTree(); // Add JM. Oct 2010
    delete gtree;
    
    
    return Tht;
}


AnyType Op_trunc_mesh3::Op::operator()(Stack stack)  const {
    
  Mesh3 *pTh = GetAny<Mesh3 *>((*getmesh)(stack));
  Mesh3 &Th = *pTh;
  long kkksplit =arg(0,stack,1L);
  long label =arg(1,stack,2L);
   KN<long> * pn2o =  arg(2,stack);
KN<long> * po2n =  arg(3,stack);
  bool renum=arg(4,stack,true);// not use to day to by compatible with 2d version ...
    
  KN<int> split(Th.nt);
  split=kkksplit;
  MeshPoint *mp= MeshPointStack(stack),mps=*mp;
  long kk=0;
  long ks=kkksplit*kkksplit*kkksplit;
  for (int k=0;k<Th.nt;k++)
    { 
      const Tet & K( Th.elements[k] );
      R3 B(1./4.,1./4.,1./4.);  // 27/09/10 : J.Morice error in msh3.cpp
      mp->set(Th,K(B),B,K,0);
      if (  GetAny<bool>( (*bbb)(stack) )  ) kk++;
      else  split[k]=0  ;    
    }
  //*mp=mps;
  if (verbosity>1) 
    cout << "  -- Trunc mesh: Nb of Tetrahedrons = " << kk << " label=" <<label <<endl;
  Mesh3 * Tht = truncmesh(Th,kkksplit,split,false,label);
  
    if(pn2o)
    {
        pn2o->resize(kk*ks);
        KN<long> &n2o(*pn2o);
        int l=0;
        for(int k=0; k< Th.nt; ++k)
            if( split[k] )
                for(int i=0; i< ks; ++i)
                    n2o[l++] = k;
    }
    if(po2n)
    {
        po2n->resize(Th.nt);
        KN<long> &o2n(*po2n);
        int l=0;
        for(int k=0; k< Th.nt; ++k)
            if( split[k] )
            {
                o2n[k] = l;
                l+=ks;
            }
            else o2n[k]=-1;
    }
  // if(renum) Tht->renum(); do not  exist ???? FH ....
  Add2StackOfPtr2FreeRC(stack,Tht);//  07/2008 FH
  *mp=mps;
  return Tht;
 };
 
//////////////////
// new functions added by J. Morice 05/10
//  --  extractmesh2D
//  --  extractmesh
//  --  movemesh32       // projection

class ExtractMesh2D_Op : public E_F0mps 
	{
	public:
		Expression eTh;
		static const int n_name_param = 4; // 
		static basicAC_F0::name_and_type name_param[] ;
		Expression nargs[n_name_param];
		KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
		{ return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
	        KN_<long>  arg(int i,int j,Stack stack,KN_<long> a ) const
	        { 
		  if( nargs[i]) 
		    return GetAny<KN_<long> >( (*nargs[i])(stack) );
		  else 
		    return nargs[j] ? GetAny<KN_<long> >( (*nargs[j])(stack) ): a;
		}
		double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
		//long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny< long >( (*nargs[i])(stack) ): a;}
		long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
	public:
		ExtractMesh2D_Op(const basicAC_F0 &  args,Expression tth) 
		: eTh(tth)
		{
			if(verbosity >1) cout << "construction par ExtractMesh_Op" << endl;
			args.SetNameParam(n_name_param,name_param,nargs); 
			if( nargs[1] && nargs[3] ) 
			  CompileError("uncompatible extractmesh (Th, region= , reft=  ");
			if( nargs[0] && nargs[2] ) 
			  CompileError("uncompatible extractmesh (Th, label= , refe=  ");
		} 
		
		AnyType operator()(Stack stack)  const ;
	};

basicAC_F0::name_and_type ExtractMesh2D_Op::name_param[]= {
{  "refe", &typeid(KN_<long>)},
{  "reft", &typeid(KN_<long>)},
{  "label", &typeid(KN_<long>)},
{  "region", &typeid(KN_<long>)}
};


class   ExtractMesh2D : public OneOperator { public:  
    ExtractMesh2D() : OneOperator(atype<pmesh>(),atype<pmesh>()) {}
	
	E_F0 * code(const basicAC_F0 & args) const 
	{
		return  new ExtractMesh2D_Op(args,t[0]->CastTo(args[0])); 
	}
};

AnyType ExtractMesh2D_Op::operator()(Stack stack)  const 
{
	MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
	Mesh * pTh= GetAny<Mesh *>((*eTh)(stack));
	Mesh &Th=*pTh;
	
	int nbv=Th.nv; // nombre de sommet 
	int nbt=Th.nt; // nombre de triangles
	int nbe=Th.neb; // nombre d'aretes fontiere
	
	KN<long> zzempty(0);
	// recuperation des labels pour les surfaces et les elements
	KN<long> labelface (arg(0,2,stack,zzempty));  
	KN<long> labelelement (arg(1,3,stack,zzempty));  
	
	// a trier les tableaux d'entier 
	
	int nv=0, nt=0, ns=0;
	int boolsurf=0;
	int boolelement=0;
	if( nargs[0] || nargs[2]){
	  boolsurf=1;
	}
	if( nargs[1] || nargs[3]){
	  boolelement=1;
	  if( boolsurf == 0){ 
	    cerr << "we need edge label to extract the mesh "<< " (This will be changed in a new version of freefem++)" << endl;
	    exit(1);
	  }
	  // a prevoir une fonction pour recalculer les differents labels des surfaces	  
	}
	
	if( boolsurf==1 && boolelement==1 ){
		// Volumetric mesh
		
		KN<int> takevertex(Th.nv,-1);
		KN<int> takeelem(Th.nt,0);
		KN<int> takebe(Th.neb,0);
		int neLab=0;
		for(int ie=0; ie<Th.nt; ie++){
		  const Mesh::Triangle & K(Th.t(ie));			
		  for(int ii=0; ii<labelelement.N(); ii++){
		    if( K.lab == labelelement[ii] ){
		      neLab++;
		      takeelem[ie] = 1;
		      for(int jj=0; jj<3; jj++){
			if( takevertex[ Th.operator()(K[jj]) ] != -1) continue;
			takevertex[ Th.operator()(K[jj]) ] = nv;
			nv++;
		      }
		      break;		      
		    }
		  }
		}
		
		
		int nbeLab=0;
		for(int ibe=0; ibe<Th.neb; ibe++){
		  const Mesh::BorderElement & K(Th.be(ibe));			
		  for(int ii=0; ii<labelface.N(); ii++){
		    if( K.lab == labelface[ii] ){
		      nbeLab++;
		      takebe[ibe] = 1;
		      for(int jj=0; jj<2; jj++){
			if( takevertex[ Th.operator()(K[jj]) ] == -1){
			  cerr << "Volumetric extract mesh: error in the definition of surface mesh " <<endl;
			  cerr << "The edges of label "<< labelface[ii]  << " is not the extract mesh " << endl;
			}
		      }
		      break;
		    }
		  }
		}
		
		ns = nbeLab;
		nt = neLab;
		
		if( nt == 0) { cerr << " The label of 2D element is not correct :: there is no element" << endl; exit(1); } 
		Mesh::Vertex   *v  = new Mesh::Vertex[nv]; 
		Mesh::Triangle *t  = new Mesh::Triangle[nt];
		Mesh::Triangle *tt = t;
		Mesh::BorderElement *b  = new Mesh::BorderElement[ns];
		Mesh::BorderElement *bb = b;
		
		for(int ii=0; ii<Th.nv; ii++){
		  if( takevertex[ ii ] == -1) continue;
		  int iv = takevertex[ii];
		  assert(iv>=0 && iv< nv);
		  v[iv].x = Th.vertices[ii].x;
		  v[iv].y = Th.vertices[ii].y;
		  v[iv].lab = Th.vertices[ii].lab;
		}
		
		for(int ie=0; ie<Th.nt; ie++){
		  if( takeelem[ ie ] != 1) continue;
		  const Mesh::Triangle & K(Th.t(ie));	
		  int ivv[3];
		  for(int jj=0; jj<3; jj++) ivv[jj] =  takevertex[ Th.operator()(K[jj]) ];      
		  (tt++)->set(v,ivv[0],ivv[1],ivv[2],K.lab);
		}
		
		
		for(int ibe=0; ibe<Th.neb; ibe++){
		  if( takebe[ ibe ] != 1) continue;
		  const Mesh::BorderElement &K( Th.be(ibe) );
		  int ivv[2];
		  for(int jj=0; jj<2; jj++) ivv[jj] =  takevertex[ Th.operator()(K[jj]) ];      
		  (bb++)->set(v,ivv[0],ivv[1],K.lab);
		}
		
		
		Mesh *pThnew = new Mesh(nv,nt,ns,v,t,b);  // attention aux composantes connexes.
		//Lorenzo
		R2 Pn,Px;
		pThnew->BoundingBox(Pn,Px);
		if(!pThnew->quadtree)
			pThnew->quadtree=new Fem2D::FQuadTree(pTh,Pn,Px,pTh->nv);
		//Lorenzo
				
		return pThnew;
		
	}
	
	
	
	
	if( boolsurf==1 && boolelement==0 ){
		// case only surface mesh
		// demander à F. pour la pertinence
		KN<int> takevertex(Th.nv,-1);
		KN<int> takebe(Th.neb,0);
		int nbeLab=0;
		for(int ibe=0; ibe<Th.neb; ibe++){
		  const Mesh::BorderElement & K(Th.be(ibe));			
		  for(int ii=0; ii<labelface.N(); ii++){
		    if( K.lab == labelface[ii] ){
		      nbeLab++;
		      takebe[ibe] = 1;
		      for(int jj=0; jj<2; jj++){
			if( takevertex[ Th.operator()(K[jj]) ] != -1) continue;
			takevertex[ Th.operator()(K[jj]) ] = nv;
			nv++;
		      }
		      break;
		    }
		  }
		}
		
		ns = nbeLab;
		
		Mesh::Vertex        *v  = new Mesh::Vertex[nv]; 
		Mesh::BorderElement *b  = new Mesh::BorderElement[ns];
		Mesh::BorderElement *bb = b;
		
		for(int ii=0; ii<Th.nv; ii++){
			if( takevertex[ ii ] == -1) continue;
			int iv = takevertex[ii];
			assert(iv>=0 && iv< nv);
			v[iv].x = Th.vertices[ii].x;
			v[iv].y = Th.vertices[ii].y;
			v[iv].lab = Th.vertices[ii].lab;
		}
		
		for(int ibe=0; ibe<Th.neb; ibe++){
		  if( takebe[ ibe ] != 1) continue;
		  const Mesh::BorderElement &K( Th.be(ibe) );
		  int ivv[2];
		  for(int jj=0; jj<2; jj++) ivv[jj] =  takevertex[ Th.operator()(K[jj]) ];      
		  (bb++)->set(v,ivv[0],ivv[1],K.lab);
		}
		
		
		Mesh *pThnew =0; //= new Mesh(nv,ns,v,b);  // a definir
		// attention aux composantes connexes.
		ffassert(0); 
		
		return pThnew;
		
	}
	return (Mesh*) 0;
	
}

//////////////////////////////////
//
//
//      Fin du 2D
//
//
//////////////////////////////////

//////////////////////////////////
//
//
//      Debut du 3D
//
//
///////////////////////////////////////////////////////////
//  recuperer une partie du maillage à l'aide des labels

class ExtractMesh_Op : public E_F0mps 
{
public:
  Expression eTh;
  static const int n_name_param = 4; // 
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<long>  arg(int i,int j,Stack stack,KN_<long> a ) const
  { 
    if( nargs[i]) 
      return GetAny<KN_<long> >( (*nargs[i])(stack) );
    else 
      return nargs[j] ? GetAny<KN_<long> >( (*nargs[j])(stack) ): a;
  }
  double  arg(int i,Stack stack,double a ) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  long    arg(int i,Stack stack,long a ) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
public:
  ExtractMesh_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth)
  {
    if(verbosity >1) cout << "construction par ExtractMesh_Op" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);    

    if( nargs[1] && nargs[3] ) 
      CompileError("uncompatible extractmesh (Th, region= , reftet=  ");
    if( nargs[0] && nargs[2] ) 
      CompileError("uncompatible extractmesh (Th, label= , refface=  ");
  } 
  
  AnyType operator()(Stack stack)  const ;
};

basicAC_F0::name_and_type ExtractMesh_Op::name_param[]= {
  {  "refface", &typeid(KN_<long>)},
  {  "reftet", &typeid(KN_<long>)},
  {  "label", &typeid(KN_<long>)},
  {  "region", &typeid(KN_<long>)},
};


class   ExtractMesh : public OneOperator { public:  
    ExtractMesh() : OneOperator(atype<pmesh3>(),atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  {
    return  new ExtractMesh_Op(args,t[0]->CastTo(args[0])); 
  }
};

AnyType ExtractMesh_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 *pTh= GetAny<Mesh3 *>((*eTh)(stack));
  Mesh3  &Th=*pTh;
  
  int nbv=Th.nv; // nombre de sommet 
  int nbt=Th.nt; // nombre de triangles
  int nbe=Th.nbe; // nombre d'aretes fontiere
  
  KN<long> zzempty(0);
  KN<long> labelface (arg(0,2,stack,zzempty));  
  KN<long> labelelement  (arg(1,3,stack,zzempty));  
 
  // a trier les tableaux d'entier 

  int nv=0, nt=0, ns=0;
  int boolsurf=0;
  int boolelement=0;
  if( nargs[0] || nargs[2] ){
    boolsurf=1;
  }
  if( nargs[1] || nargs[3] ){
    boolelement=1;
  }
  cout << " labelface.N()  " << labelface.N() << endl;
  for(int ii=0; ii<labelface.N(); ii++)
    cout << ii << " " << labelface[ii] << endl;
 
  cout << " labelelement.N()  " << labelelement.N() << endl;
  for(int ii=0; ii<labelelement.N(); ii++)
    cout << ii << " " << labelelement[ii] << endl;


  if( boolsurf==1 && boolelement==1 ){
    // case only surface mesh
    KN<int> takevertex(Th.nv,-1);
    KN<int> takeelem(Th.nt,0);
    KN<int> takebe(Th.nbe,0);
    int neLab=0;

    for(int ie=0; ie<Th.nt; ie++){
      const Tet & K(Th.elements[ie]);		
      for(int ii=0; ii<labelelement.N(); ii++){
	if( K.lab == labelelement[ii] ){	  
	  neLab++;
	  takeelem[ie] = 1;
	  for(int jj=0; jj<4; jj++){
	    if( takevertex[ Th.operator()(K[jj]) ] != -1) continue;
	    takevertex[ Th.operator()(K[jj]) ] = nv;
	    nv++;
	  }
	  break;
	}
      }
    }
        
    int nbeLab=0;
    for(int ibe=0; ibe<Th.nbe; ibe++){
      const Triangle3 & K(Th.be(ibe));			
      for(int ii=0; ii<labelface.N(); ii++){
	if( K.lab == labelface[ii] ){
	  nbeLab++;
	  takebe[ibe] = 1;
	  for(int jj=0; jj<3; jj++){
	    if( takevertex[ Th.operator()(K[jj]) ] == -1){
	      int iivv=Th.operator()(K[jj]);
	      cerr << "Error(Extract Mesh):: the vertex ("<< Th.vertices[iivv].x  << ", " <<Th.vertices[iivv].y << ", " <<Th.vertices[iivv].z <<") is not in the volume mesh" << endl;
	    }
	    //continue;
	    //takevertex[ Th.operator()(K[jj]) ] = nv;
	    //nv++;
	  }
	  break;
	}
      }
    }
    
    ns = nbeLab;
    nt = neLab;
    cout << nt <<" " << ns << endl;
    if( nt == 0) { cerr << " The label of 3D element is not correct :: there is no element" << endl; exit(1); } 
    Vertex3 *v  = new Vertex3[nv]; 
    Tet *t  = new Tet[nt];
    Tet *tt = t;
    Triangle3 *b  = new Triangle3[ns];
    Triangle3 *bb = b;
    
    for(int ii=0; ii<Th.nv; ii++){
      if( takevertex[ ii ] == -1) continue;
      int iv = takevertex[ii];
      assert(iv>=0 && iv< nv);
      v[iv].x = Th.vertices[ii].x;
      v[iv].y = Th.vertices[ii].y;
      v[iv].z = Th.vertices[ii].z;
      v[iv].lab = Th.vertices[ii].lab;
    }
    
    for(int ie=0; ie<Th.nt; ie++){
      if( takeelem[ ie ] != 1) continue;
      const Tet & K(Th.elements[ie]);	
      int ivv[4];
      for(int jj=0; jj<4; jj++) ivv[jj] =  takevertex[ Th.operator()(K[jj]) ];      
      (tt++)->set(v,ivv,K.lab);
    }
		
    
    for(int ibe=0; ibe<Th.nbe; ibe++){
      if( takebe[ ibe ] != 1) continue;
      const Triangle3 &K( Th.be(ibe) );
      int ivv[3];
      for(int jj=0; jj<3; jj++) ivv[jj] =  takevertex[ Th.operator()(K[jj]) ];      
      (bb++)->set(v,ivv,K.lab);
    }
    
    
    Mesh3 *pThnew = new Mesh3(nv,nt,ns,v,t,b);  // peut etre a dÃ©finir ???
    // attention aux composantes connexes.
	pThnew->BuildGTree();  //Lorenzo
    
    return pThnew;
    
  }
  
  if( boolsurf==1 && boolelement==0 ){
		
		
    KN<int> takevertex(Th.nv,-1);
    for(int iv=0; iv<Th.nv; iv++) takevertex[iv]=-1;
    KN<int> takebe(Th.nbe,-1);
    for(int iv=0; iv<Th.nv; iv++) takebe[iv]=-1;
    int nbeLab=0;
    for(int ibe=0; ibe<Th.nbe; ibe++){
      const Triangle3 &K( Th.be(ibe) );
      
      for(int ii=0; ii<labelface.N(); ii++){
	
	if( K.lab == labelface[ii] ){
	  nbeLab++;
	  cout << "takeface" << endl;
	  takebe[ibe] = 1;
	  for(int jj=0; jj<3; jj++){
	    if( takevertex[ Th.operator()(K[jj]) ] != -1) continue;
	    takevertex[ Th.operator()(K[jj]) ] = nv;
	    nv++;
	  }
					
	  //break;
	}
      }
    }
    
    ns = nbeLab;
    
    Vertex3   *v  = new Vertex3[nv]; 
    Triangle3 *b  = new Triangle3[ns];
    Triangle3 *bb = b;
    
    for(int ii=0; ii<Th.nv; ii++){
      if( takevertex[ ii ] == -1) continue;
      int iv = takevertex[ii];
      assert(iv>=0 && iv< nv);
      v[iv].x = Th.vertices[ii].x;
      v[iv].y = Th.vertices[ii].y;
      v[iv].z = Th.vertices[ii].z;
		  v[iv].lab = Th.vertices[ii].lab;
    }
    
    for(int ibe=0; ibe<Th.nbe; ibe++){
      if( takebe[ ibe ] != 1) continue;
      const Triangle3 &K( Th.be(ibe) );
      int ivv[3];
      for(int jj=0; jj<3; jj++) ivv[jj] =  takevertex[ Th.operator()(K[jj]) ];      
      (bb++)->set(v,ivv,K.lab);
    }
    
    
    cout <<" nv" << nv << " ns " << endl;  
    Mesh3 *pThnew = new Mesh3(nv,ns,v,b);  
	pThnew->BuildGTree();  //Lorenzo	
    return pThnew;    
  }
   
  Mesh3 *pThnew = NULL;
  return pThnew;    
}

bool AddLayers(Mesh3 const  * const & pTh, KN<double> * const & psupp, long const & nlayer,KN<double> * const & pphi)
{
    ffassert(pTh && psupp && pphi);
    const int nve = Mesh3::Element::nv;
    const Mesh3 & Th= *pTh;
    const int nt = Th.nt;
    const int nv = Th.nv;
    
    KN<double> & supp(*psupp);
    KN<double> u(nv), s(nt);
    KN<double> & phi(*pphi);
    ffassert(supp.N()==nt);//P0
    ffassert(phi.N()==nv); // P1
    s = supp;
    phi=0.;
    // supp = 0.;
    // cout << " s  " << s << endl;
    
    for(int step=0; step < nlayer; ++ step)
    {
        
        
        u = 0.;
        for(int k=0; k<nt; ++k)
            for(int i=0; i<nve; ++i)
                u[Th(k,i)] += s[k];
        
        for(int v=0; v < nv; ++v)
            u[v] = u[v] >0.;
        // cout << " u  " << u << endl;
        
        phi += u;
        
        s = 0.;
        for(int k=0; k<nt; ++k)
            for(int i=0; i<nve; ++i)
                s[k] += u[Th(k,i)];
        
        for(int k=0; k < nt; ++k)
            s[k] = s[k] > 0.;
        supp += s;
        // cout << " s  " << s << endl;
    }
    // cout << " phi  " << phi << endl;
    phi *= (1./nlayer);
    // supp =s;
    return true;
}

Mesh3 * GluMesh3tab(KN<pmesh3> * const & tab, long const & lab_delete)
{
    int flagsurfaceall = 0;
    
    int nbt=0;
    int nbe=0;
    int nbex=0;
    int nbv=0;
    int nbvx=0;
    
    double hmin=1e100;
    R3 Pn(1e100,1e100,1e100),Px(-1e100,-1e100,-1e100);
    const Mesh3 * th0=0;
    
    for(int i = 0;i<tab->n;i++)
    {
        const Mesh3 &Th3(*tab->operator[](i));
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
    
    if(verbosity > 1) cout << "      - hmin =" <<  hmin << " ,  Bounding Box: " << Pn << " "<< Px << endl;
    
    // probleme memoire
    Vertex3  *v= new Vertex3[nbvx];
    Tet      *t;
    if(nbt!=0) t= new Tet[nbt];
    Tet      *tt=t;
    Triangle3 *b= new Triangle3[nbex];
    Triangle3 *bb= b;
    
    ffassert(hmin>Norme2(Pn-Px)/1e9);
    double hseuil =hmin/10.;
    
    //int *NumSom= new int[nbvx];
    
    // VERSION morice
    if(verbosity > 1) cout << " creation of : BuildGTree" << endl;
    EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(v,Pn,Px,0);
    
    nbv=0;
    //int nbv0=0;
    for(int i = 0;i<tab->n;i++)
    {
        const Mesh3 &Th3(*tab->operator[](i));
        
        if(verbosity>1)  cout << " loop over mesh for create new mesh "<< endl;
        if(verbosity>1)  cout << " GluMesh3D + "<< Th3.nv << " " << Th3.nt <<" " << Th3.nbe << endl;
        //nbv0 =+Th3.nv;
        
        for (int ii=0;ii<Th3.nv;ii++){
            const Vertex3 &vi(Th3.vertices[ii]);
            Vertex3 * pvi=gtree->ToClose(vi,hseuil);
            
            
            if(!pvi){
                v[nbv].x = vi.x;
                v[nbv].y = vi.y;
                v[nbv].z = vi.z;
                v[nbv].lab = vi.lab;
                //NumSom[ii+nbv0] = nbv;
                gtree->Add( v[nbv] );
                nbv++;
            }
            /*
             else{
             NumSom[ii+nbv0] = pvi-v;
             assert(pvi-v <nbv);
             }
             */
        }
        
        for (int k=0;k<Th3.nt;k++){
            const Tet  &K(Th3.elements[k]);
            int iv[4];
            iv[0]=gtree->ToClose(K[0],hseuil)-v;
            iv[1]=gtree->ToClose(K[1],hseuil)-v;
            iv[2]=gtree->ToClose(K[2],hseuil)-v;
            iv[3]=gtree->ToClose(K[3],hseuil)-v;
            (tt++)->set(v,iv,K.lab);
        }
        //nbv0 =+Th3.nv;
    }
    
    
    if(verbosity > 1) cout << " creation of : BuildGTree for border elements" << endl;
    Vertex3  *becog= new Vertex3[nbex];
    //Vertex3  becog[nbex];
    EF23::GTree<Vertex3> *gtree_be = new EF23::GTree<Vertex3>(becog,Pn,Px,0);
    
    double hseuil_border = hseuil/3;
    //nbv0=0;
    for(int i = 0;i<tab->n;i++)
    {
        const Mesh3 &Th3(*tab->operator[](i));
        
        for (int k=0;k<Th3.nbe;k++)
        {
            const Triangle3 & K(Th3.be(k));
            if ((K.lab != lab_delete))//&&(K.lab != 3))
            {
                
                int iv[3];
                iv[0]=Th3.operator()(K[0]);
                iv[1]=Th3.operator()(K[1]);
                iv[2]=Th3.operator()(K[2]);
                
                R cdgx,cdgy,cdgz;
                
                cdgx = (Th3.vertices[iv[0]].x+ Th3.vertices[iv[1]].x+ Th3.vertices[iv[2]].x)/3.;
                cdgy = (Th3.vertices[iv[0]].y+ Th3.vertices[iv[1]].y+ Th3.vertices[iv[2]].y)/3.;
                cdgz = (Th3.vertices[iv[0]].z+ Th3.vertices[iv[1]].z+ Th3.vertices[iv[2]].z)/3.;
                
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
                    igluv[0]= gtree->ToClose(K[0],hseuil)-v; //NumSom[iv[0]+nbv0];
                    igluv[1]= gtree->ToClose(K[1],hseuil)-v; //NumSom[iv[1]+nbv0];
                    igluv[2]= gtree->ToClose(K[2],hseuil)-v; //NumSom[iv[2]+nbv0];
                    
                    (bb++)->set(v,igluv,K.lab);
                }
            }
        }
        
        //nbv0 =+Th3.nv;
    }
    delete gtree;
    delete gtree_be;
    delete [] becog;
    
    if(verbosity > 2) cout << " nbv="  << nbv  << endl;
    if(verbosity > 2) cout << " nbvx=" << nbvx << endl;
    if(verbosity > 2) cout << " nbt="  << nbt  << endl;
    if(verbosity > 2) cout << " nbe="  << nbe  << endl;
    if(verbosity > 2) cout << " nbex=" << nbex << endl;
    if(verbosity>1)
    {
        cout << "     Nb of glu3D  point " << nbvx-nbv;
        cout << "     Nb of glu3D  Boundary faces " << nbex-nbe << endl;
    }
    
    if(nbt==0){
        Mesh3 *mpq= new Mesh3(nbv,nbe,v,b);  
        if(flagsurfaceall==1) mpq->BuildBoundaryElementAdj();
        return mpq;
    }
    else{
        Mesh3 *mpq= new Mesh3(nbv,nbt,nbe,v,t,b); 
        /* 
         mpq->BuildBound();
         if(verbosity > 1) cout << "fin de BuildBound" << endl;
         mpq->BuildAdj();
         if(verbosity > 1) cout << "fin de BuildAdj" << endl;
         mpq->Buildbnormalv();  
         if(verbosity > 1) cout << "fin de Buildnormalv()" << endl;
         mpq->BuildjElementConteningVertex();
         if(verbosity > 1) cout << "fin de ConteningVertex()" << endl;
         */
        mpq->BuildGTree();
        if(verbosity > 2) cout << "fin de BuildGTree()" << endl;
        
        //Add2StackOfPtr2FreeRC(stack,mpq);
        
        return mpq;
    }
}

struct Op_GluMesh3tab : public OneOperator {
    typedef const Mesh3 *pmesh3;
    class Op: public E_F0mps   {
    public:
        static basicAC_F0::name_and_type name_param[] ;
        static const int n_name_param =1;
        Expression nargs[n_name_param];
        
        Expression getmeshtab;
        long arg(int i,Stack stack,long a) const{ return nargs[i] ? GetAny<long>( (*nargs[i])(stack) ): a;}
        
        Op(const basicAC_F0 &  args,Expression t) : getmeshtab(t)
        { args.SetNameParam(n_name_param,name_param,nargs); }
        AnyType operator()(Stack s)  const ;
    };
    
    E_F0 * code(const basicAC_F0 & args) const
    { return new Op(args,t[0]->CastTo(args[0])); }
    Op_GluMesh3tab() :
    OneOperator(atype<const pmesh3>(),atype<KN<pmesh3>*>()) {};
};

basicAC_F0::name_and_type Op_GluMesh3tab::Op::name_param[Op_GluMesh3tab::Op::n_name_param] =
{
    {  "labtodel",             &typeid(long)}
};

AnyType Op_GluMesh3tab::Op::operator()(Stack stack)  const {
    
    KN<const Mesh3*> *tab = GetAny<KN<const Mesh3*> *>((*getmeshtab)(stack));
    long labtodel = arg(0,stack,LONG_MIN);
    
    Mesh3 * Tht = GluMesh3tab(tab,labtodel);
    
    Add2StackOfPtr2FreeRC(stack,Tht);
    return Tht;
}



// <<WITH_NO_INIT>> because i include this file in tetgen.cpp (very bad) [[file:tetgen.cpp::WITH_NO_INIT]]
#ifndef WITH_NO_INIT

// <<dynamic_loading>>

static void Load_Init()
{  
  
  Dcl_Type<listMesh3>();
  typedef const Mesh *pmesh;
  typedef const Mesh3 *pmesh3;
  
  if (verbosity>1  && mpirank == 0)
    cout << " load: msh3  "  <<endl;
  
  TheOperators->Add("+",new OneBinaryOperator_st< Op3_addmesh<listMesh3,pmesh3,pmesh3>  >      );
  TheOperators->Add("+",new OneBinaryOperator_st< Op3_addmesh<listMesh3,listMesh3,pmesh3>  >      );
  //TheOperators->Add("=",new OneBinaryOperator< Op3_setmesh<false,pmesh3*,pmesh3*,listMesh3>  >     );
  //TheOperators->Add("<-",new OneBinaryOperator< Op3_setmesh<true,pmesh3*,pmesh3*,listMesh3>  >     );
  
  TheOperators->Add("=",new OneBinaryOperator_st< Op3_setmesh<false,pmesh3*,pmesh3*,listMesh3>  >     );
  TheOperators->Add("<-",new OneBinaryOperator_st< Op3_setmesh<true,pmesh3*,pmesh3*,listMesh3>  >     );

  Global.Add("change","(",new SetMesh3D);
  Global.Add("movemesh23","(",new Movemesh2D_3D_surf);
  Global.Add("movemesh2D3Dsurf","(",new Movemesh2D_3D_surf_cout);
  Global.Add("movemesh3","(",new Movemesh3D);
  Global.Add("movemesh","(",new Movemesh3D(1));
  Global.Add("movemesh3D","(", new Movemesh3D_cout);
  Global.Add("deplacement","(",new DeplacementTab);
  Global.Add("checkbemesh","(",new CheckManifoldMesh);  
  Global.Add("buildlayers","(",new  BuildLayerMesh);  
  Global.Add("cube","(",new  cubeMesh);
  Global.Add("cube","(",new  cubeMesh(1));
  Global.Add("trunc","(", new Op_trunc_mesh3);
  Global.Add("gluemesh","(",new Op_GluMesh3tab);

  Global.Add("extract","(",new ExtractMesh);
  Global.Add("extract","(",new ExtractMesh2D);
    
  Global.Add("AddLayers","(",new OneOperator4_<bool, const Mesh3 * , KN<double> *,long, KN<double> * >(AddLayers));
  typedef const Mesh3 *pmesh3;
  // Global.Add("trunc","(", new Op_trunc_mesh3);
}

// <<msh3_load_init>> static loading: calling Load_Init() from a function which is accessible from
// [[file:~/ff/src/fflib/load.cpp::static_load_msh3]]

void msh3_Load_Init(){Load_Init();}

// dynamic loading: calling [[file:../src/fflib/InitFunct.hpp::LOADFUNC]] on [[Load_Init]]
LOADFUNC(Load_Init)
    
#endif // [[WITH_NO_INIT]]
