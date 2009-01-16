// ORIG-DATE:     Janvier 2009
// -*- Mode : c++ -*-
//
// SUMMARY  : interface avec le logiciel du domaine public netgen    
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
//#include "LayerMesh.hpp"
//#include "TransfoMesh_v2.hpp"
#include "msh3.hpp"
//#include "GQuadTree.hpp"

#include <set>
#include <vector>
#include <fstream>

namespace nglib {
#include "./includenetgen/nglib.h"
}
using namespace nglib;

using namespace Fem2D;


class RemplissageNetgen_Op : public E_F0mps 
{
public:
  Expression eTh;
  Expression xx,yy,zz;
  static const int n_name_param =5; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  int  arg(int i,Stack stack, int a) const{ return nargs[i] ? GetAny< int >( (*nargs[i])(stack) ): a;}
  string*  arg(int i,Stack stack, string* a) const{ return nargs[i] ? GetAny< string* >( (*nargs[i])(stack) ): a;}
public:
  RemplissageNetgen_Op(const basicAC_F0 &  args,Expression tth) 
    : eTh(tth)
  {
    if(verbosity) cout << "construction par RemplissageNetgen_Op" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
  
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type RemplissageNetgen_Op::name_param[]= {
  {  "reftet", &typeid(long)}, 
  {  "refface", &typeid(KN_<long>)},
  //  Netgen Options
  {  "maxh", &typeid(double)},
  {  "fineness", &typeid(double)},
  {  "secondorder", &typeid(long)}
};

class  RemplissageNetgen : public OneOperator { public:  
    RemplissageNetgen() : OneOperator(atype<pmesh3>(),atype<pmesh3>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
	return  new RemplissageNetgen_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

AnyType RemplissageNetgen_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  Mesh3 * pTh3= GetAny<Mesh3 *>((*eTh)(stack));
  ffassert( pTh3 );
  Mesh3 &Th3=*pTh3;
  Mesh3 *m= pTh3;   // question a quoi sert *m ??

  // lecture des arguments
  KN<long> zzempty;
  long     nrtet (arg(0, stack, 1));
  KN<long> nrf (arg(1,stack,zzempty)); 
  double   netgen_maxh (arg( 2, stack, 1.)); 
  double  netgen_fineness (arg( 3, stack, 1.)); 
  int  netgen_secondorder (arg( 4, stack, 0));

  ffassert( nrf.N() %2 ==0);

  int nbv=Th3.nv; // nombre de sommet 
  int nbt=Th3.nt; // nombre de triangles
  int nbe=Th3.nbe; // nombre d'aretes fontiere

  if(nbt!=0) {
    cerr << " The mesh must be a 3D surface mesh " << endl;
    exit(1);
  }
  if(verbosity >1) cout <<" ======================= " << endl;
  if(verbosity >1) cout <<" == RemplissageNetgen == " << endl;
 

  Ng_Mesh * netgen_mesh;
  Ng_Init();
  // creates mesh structure
  netgen_mesh = Ng_NewMesh ();

  double point[3];
  int trig[3], tet[4];

  for (int ii = 0; ii < Th3.nv; ii++)
    {
      point[0] = Th3.vertices[ii].x;
      point[1] = Th3.vertices[ii].y;
      point[2] = Th3.vertices[ii].z;
      
      Ng_AddPoint (netgen_mesh, point);    
    }
  
  for (int ii = 0; ii < Th3.nbe; ii++)
    {
      const Triangle3 & K(Th3.be(ii));
      for(int jj=0; jj <3; jj++)
	trig[jj] = Th3.operator()(K[jj])+1;
      
      Ng_AddSurfaceElement (netgen_mesh, NG_TRIG, trig);       
    }
  
  Ng_Meshing_Parameters netgen_mp;
  netgen_mp.maxh     = 1;     //netgen_maxh;
  netgen_mp.fineness = 1.;    //netgen_fineness;
  netgen_mp.secondorder = 0;  //netgen_secondorder;
 

  //cout << "start meshing" << endl;
  //Ng_GenerateVolumeMesh (netgen_mesh, &netgen_mp);
  //cout << "meshing done" << endl;
  
  netgen_mp.maxh     = 0.09;     //netgen_maxh;
  netgen_mp.fineness = 1.;    //netgen_fineness;
  netgen_mp.secondorder = 0;  //netgen_secondorder;
  cout << "start remeshing" << endl;
  Ng_GenerateVolumeMesh (netgen_mesh, &netgen_mp);
  cout << "meshing done" << endl;

  /* Transformation netgen -> freefem++ */

  map<int,int> mapface;
  for(int i=0;i<nrf.N();i+=2)
    { 
      if( nrf[i] != nrf[i+1] )	
	mapface[nrf[i]] = nrf[i+1];     
    }
      
  // read information of netgen
  int netgen_nv = Ng_GetNP(netgen_mesh);
  int netgen_nt = Ng_GetNE(netgen_mesh);
  int netgen_nbe = Ng_GetNSE(netgen_mesh);
  
  Vertex3   *v = new Vertex3[netgen_nv];
  Tet       *t = new Tet[netgen_nt];
  Triangle3 *b = new Triangle3[netgen_nbe];
  // generation des nouveaux sommets 
  Vertex3   *vv=v;
  Tet       *tt=t;
  Triangle3 *bb=b;

  cout << " donnee sortie netgen:  Vertex" << netgen_nv << " Tetrahedre " << netgen_nt << "  " << netgen_nbe << endl;

  for (int ii = 0; ii < netgen_nv; ii++)
    {
      Ng_GetPoint (netgen_mesh, ii+1, point);
      vv->x=point[0];
      vv->y=point[1];
      vv->z=point[2];
      vv->lab = 1;
      vv++;      
    }
  for (int ii = 0; ii < netgen_nt; ii++)
    {
      int iv[4];
      Ng_GetVolumeElement (netgen_mesh, ii+1, iv);
      for(int jj=0;jj<4;jj++)
	iv[jj]=iv[jj]-1;

      (*tt++).set( v, iv, nrtet);
    }
  for (int ii = 0; ii < netgen_nbe; ii++)
    {
      int label=0;
      int iv[3];
      Ng_GetSurfaceElement (netgen_mesh, ii+1, iv);
      for(int jj=0;jj<3;jj++)
	iv[jj]=iv[jj]-1;
      (*bb++).set(v, iv ,label);
    }


  Ng_DeleteMesh (netgen_mesh);
  Ng_Exit();

  Mesh3 * Th3_t = new Mesh3(netgen_nv,netgen_nt,netgen_nbe,v,t,b);

  Th3_t->BuildBound();
  Th3_t->BuildAdj();
  Th3_t->Buildbnormalv();  
  Th3_t->BuildjElementConteningVertex();
  Th3_t->BuildGTree();
  //Th3->decrement();    
  Add2StackOfPtr2FreeRC(stack,Th3_t);

		
  *mp=mps;
  
  return Th3_t;
}

class Netgen_STL_Op : public E_F0mps 
{
public:
  Expression filename;
  static const int n_name_param =4; // 
  static basicAC_F0::name_and_type name_param[] ;
  Expression nargs[n_name_param];
  KN_<long>  arg(int i,Stack stack,KN_<long> a ) const
  { return nargs[i] ? GetAny<KN_<long> >( (*nargs[i])(stack) ): a;}
  KN_<double>  arg(int i,Stack stack,KN_<double> a ) const
  { return nargs[i] ? GetAny<KN_<double> >( (*nargs[i])(stack) ): a;}
  double  arg(int i,Stack stack,double a) const{ return nargs[i] ? GetAny< double >( (*nargs[i])(stack) ): a;}
  int  arg(int i,Stack stack, int a) const{ return nargs[i] ? GetAny< int >( (*nargs[i])(stack) ): a;}
  string*  arg(int i,Stack stack, string* a) const{ return nargs[i] ? GetAny< string* >( (*nargs[i])(stack) ): a;}
public:
  Netgen_STL_Op(const basicAC_F0 &  args,Expression ffname) 
    : filename(ffname)
  {
    if(verbosity) cout << "construction par RemplissageNetgen_Op" << endl;
    args.SetNameParam(n_name_param,name_param,nargs);
  } 
  
  AnyType operator()(Stack stack)  const ;
};


basicAC_F0::name_and_type Netgen_STL_Op::name_param[]= {
  //  Netgen Options
  {  "maxh", &typeid(double)},
  {  "fineness", &typeid(double)},
  {  "secondorder", &typeid(long)},
  {  "meshsizefilename", &typeid(long)}
};

class  Netgen_STL : public OneOperator { public:  
    Netgen_STL() : OneOperator(atype<pmesh3>(),atype<string *>()) {}
  
  E_F0 * code(const basicAC_F0 & args) const 
  { 
    return  new Netgen_STL_Op( args,t[0]->CastTo(args[0]) ); 
  }
};

AnyType Netgen_STL_Op::operator()(Stack stack)  const 
{
  MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
  string * pffname= GetAny<string *>((*filename)(stack));
  size_t size_pffname = pffname->size()+1;
  char* char_pffname =new char[size_pffname];
  strncpy( char_pffname, pffname->c_str(), size_pffname); 

  // lecture des arguments
  double   netgen_maxh (arg( 2, stack, 1.)); 
  double  netgen_fineness (arg( 3, stack, 1.)); 
  int  netgen_secondorder (arg( 4, stack, 0));
  if( nargs[3] ) string *netgen_meshsize_filename =  GetAny<string *>( (*(nargs[3]))(stack) );

  if(verbosity >1) cout <<" ===================== " << endl;
  if(verbosity >1) cout <<" ==   Netgen_STL    == " << endl;
  
  int i, j, rv;

  Ng_Mesh * netgen_mesh;
  Ng_STL_Geometry * netgen_geom;

  Ng_Meshing_Parameters netgen_mp;
  netgen_mp.maxh=100000;
  netgen_mp.fineness = 0.5;
  netgen_mp.secondorder = 0;
  netgen_mp.meshsize_filename = "hinge.msz";

  Ng_Init();

  netgen_geom = Ng_STL_LoadGeometry (char_pffname);
  if (!netgen_geom)
    {
      cerr << "Ng_STL_LoadGeometry return NULL" << endl;
      exit(1);
    }

  rv = Ng_STL_InitSTLGeometry(netgen_geom);
  cout << "InitSTLGeometry: NG_result=" << rv << endl;

  netgen_mesh = Ng_NewMesh ();

  rv = Ng_STL_MakeEdges (netgen_geom, netgen_mesh, &netgen_mp);
  cout << "Make Edges: Ng_result=" << rv << endl;

  rv = Ng_STL_GenerateSurfaceMesh (netgen_geom, netgen_mesh, &netgen_mp);
  cout << "Generate Surface Mesh: Ng_result=" << rv << endl;
	
  //Ng_SaveMesh (mesh, "surface.vol");
  
  rv = Ng_GenerateVolumeMesh(netgen_mesh,&netgen_mp);
  cout << "Generate Volume Mesh: Ng_result=" << rv << endl;

  //Ng_SaveMesh (mesh, "volume.vol");
  // read information of netgen
  int netgen_nv = Ng_GetNP(netgen_mesh);
  int netgen_nt = Ng_GetNE(netgen_mesh);
  int netgen_nbe = Ng_GetNSE(netgen_mesh);
  
  Vertex3   *v = new Vertex3[netgen_nv];
  Tet       *t = new Tet[netgen_nt];
  Triangle3 *b = new Triangle3[netgen_nbe];
  // generation des nouveaux sommets 
  Vertex3   *vv=v;
  Tet       *tt=t;
  Triangle3 *bb=b;

  cout << " donnee sortie netgen:  Vertex" << netgen_nv << " Tetrahedre " << netgen_nt << "  " << netgen_nbe << endl;

  for (int ii = 0; ii < netgen_nv; ii++)
    {
      double point[3];
      Ng_GetPoint (netgen_mesh, ii+1, point);
      vv->x=point[0];
      vv->y=point[1];
      vv->z=point[2];
      vv->lab = 1;
      vv++;      
    }
  for (int ii = 0; ii < netgen_nt; ii++)
    {
      int nrtet=1;
      int iv[4];
      Ng_GetVolumeElement (netgen_mesh, ii+1, iv);
      for(int jj=0;jj<4;jj++)
	iv[jj]=iv[jj]-1;

      (*tt++).set( v, iv, nrtet);
    }
  for (int ii = 0; ii < netgen_nbe; ii++)
    {
      int label=0;
      int iv[3];
      Ng_GetSurfaceElement (netgen_mesh, ii+1, iv);
      for(int jj=0;jj<3;jj++)
	iv[jj]=iv[jj]-1;
      (*bb++).set(v, iv ,label);
    }


  Ng_DeleteMesh (netgen_mesh);
  Ng_Exit();
  Mesh3 * Th3_t = new Mesh3(netgen_nv,netgen_nt,netgen_nbe,v,t,b);

  Th3_t->BuildBound();
  Th3_t->BuildAdj();
  Th3_t->Buildbnormalv();  
  Th3_t->BuildjElementConteningVertex();
  Th3_t->BuildGTree();
  //Th3->decrement();    
  Add2StackOfPtr2FreeRC(stack,Th3_t);
		
  *mp=mps;
  
  return Th3_t;
}

class Init1 { public:
  Init1();
};

static Init1 init1;  //  une variable globale qui serat construite  au chargement dynamique 

Init1::Init1(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  
  //if (verbosity)
  if(verbosity) cout << " load: netgen  " << endl;
  Global.Add("netg","(",new RemplissageNetgen);
  Global.Add("netgstl","(",new Netgen_STL);
  
}
