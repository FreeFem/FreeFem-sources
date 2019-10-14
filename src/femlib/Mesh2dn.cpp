// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model  mesh 2d
// USAGE    : LGPL
// ORG      : LJLL Universite Pierre et Marie Curi, Paris,  FRANCE
// AUTHOR   : Frederic Hecht
// E-MAIL   : frederic.hecht@ann.jussieu.fr
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

#include <fstream>
#include <iostream>
#include "ufunction.hpp"
#include "error.hpp"
#include "RNM.hpp"
#include "libmesh5.h"


#include "Mesh2dn.hpp"
//  for plotStream (a change)
#include "Mesh3dn.hpp"
#include "MeshSn.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "PlotStream.hpp"

namespace Fem2D {

  static const int  nvfaceTet[4][3]  = { {2,1,3},{0,2,3},{1,0,3},{0,1,2} };
  static const int  nvedgeTet[6][2] = { {0,1},{0,2},{0,3},{0,1},{1,2},{2,3} };

  static const int  nvfaceTria[1][3]  = { {0,1,2} };
  static const int  nvedgeTria[3][2] = { {1,2},{2,0},{0,1}}; //  tourne de le sens trigo  donc Normal ext   vect(1,0) ^ perp

  static const int   nvfaceSeg[1][3]  = {{-1,-1,1}};
  static const int  nvedgeSeg[1][2] = { {0,1} };

  static const int  nvadjSeg[2][1] = { {0},{1} };


  template<> const int (* const GenericElement<DataSeg2>::nvface)[3] = 0 ;
  template<> const int (* const GenericElement<DataSeg2>::nvedge)[2] = nvedgeSeg; //nvedgeTria ;
  template<> const int (* const GenericElement<DataSeg2>::nvadj)[1] = nvadjSeg ;


  template<> const int (* const GenericElement<DataTriangle2>::nvface)[3] = nvfaceTria ;
  template<>  const int (* const GenericElement<DataTriangle2>::nvedge)[2] = nvedgeTria ;
  template<>  const int (* const GenericElement<DataTriangle2>::nvadj)[2] = nvedgeTria ;
  template<> const int  GenericElement<DataTriangle2>::nitemdim[4] = {3,3,1,0 }  ;


  static const int onWhatIsEdge2d[3][7] = {  {0,1,3, 2,0,0, 0}, // edge 0
					   {3,0,1, 0,2,0, 0},
					   {1,3,0, 0,0,2, 0}};

  template<>
  const int (* const GenericElement<DataTriangle2>::onWhatBorder)[7] = onWhatIsEdge2d ;

  template<> int   GenericMesh<Triangle2,BoundaryEdge2,Vertex2>::kfind=0;
  template<> int   GenericMesh<Triangle2,BoundaryEdge2,Vertex2>::kthrough=0;


Mesh2::Mesh2(const char * filename)
{ // read the mesh

  int nt,nv,nbe;
  int ok=load(filename);
  if(ok)
    {
      ifstream f(filename);
      if(!f) {cerr << "Mesh2::Mesh2 Erreur openning " << filename<<endl;exit(1);}
      if(verbosity)
      cout << " Read On file \"" <<filename<<"\""<<  endl;
      f >> nv >> nt >> nbe ;
      this->set(nv,nt,nbe);
      if(verbosity)
      cout << "  -- Nb of Vertex " << nv << " " << " Nb of Triangles " << nt
	   << " , Nb of border edges " << nbe <<  endl;
      assert(f.good() && nt && nv) ;
      for (int i=0;i<nv;i++)
	{
	  f >> this->vertices[i];
	  assert(f.good());
	}
      mes=0;
      for (int i=0;i<nt;i++)
	{
	  this->t(i).Read1(f,this->vertices,nv);
	  mes += t(i).mesure();
	}
      mesb=0.;
      for (int i=0;i<nbe;i++)
	{
	  this->be(i).Read1(f,this->vertices,nv);
	  mesb += be(i).mesure();
	}
    }
  BuildBound();
  BuildAdj();
  Buildbnormalv();
  BuildjElementConteningVertex();

  if(verbosity)
  cout << "   - mesh mesure = " << mes << " border mesure: " << mesb << endl;
}

int Mesh2::load(const string & filename)
{

  int bin;
  int ver,inm,dim;
  int lf=filename.size()+20;
  KN<char>  fileb(lf),filef(lf);
  char * pfile;
  strcpy(filef,filename.c_str());
  strcpy(fileb,filef);
  strcat(filef,".mesh");
  strcat(fileb,".meshb");
  if( (inm=GmfOpenMesh(pfile=fileb, GmfRead,&ver,&dim)) )
    bin=true;
  else if( (inm=GmfOpenMesh(pfile=filef, GmfRead,&ver,&dim)) )
    bin=false;
  else
    { cerr << " Erreur ouverture file " << (char *) fileb << " " << (char *) filef << endl;
      return   1;
    }
  int nv,nt,neb;
  nv = GmfStatKwd(inm,GmfVertices);
  nt = GmfStatKwd(inm,GmfTriangles);
  neb=GmfStatKwd(inm,GmfEdges);
  this->set(nv,nt,neb);

  if(verbosity)
  cout << pfile <<": ver " << ver << ", d "<< dim  << ", nt " << nt << ", nv " << nv << " nbe:  = " << nbe << endl;
  if(dim  != Rd::d) {
    cerr << "Err dim == " << dim << " != " << Rd::d << endl;
    return 2; }
  if( nv<=0 && nt <=0 ) {
    cerr << " missing data "<< endl;
    return 3;
  }
  int iv[4],lab;
  float cr[3]={0,0,0};
  // read vertices
  GmfGotoKwd(inm,GmfVertices);
  int mxlab=0;
  int mnlab=0;
  for(int i=0;i<nv;++i)
    {
      if(ver<2) {
	GmfGetLin(inm,GmfVertices,&cr[0],&cr[1],&lab);
	vertices[i].x=cr[0];
	vertices[i].y=cr[1];
    }
      else
	GmfGetLin(inm,GmfVertices,&vertices[i].x,&vertices[i].y,&lab);
      vertices[i].lab=lab;

      mxlab= max(mxlab,lab);
      mnlab= min(mnlab,lab);
    }


  //    /* read mesh triangles */
    {
      mes=0;
      GmfGotoKwd(inm,GmfTriangles);
      for(int i=0;i<nt;++i)
	{
	  GmfGetLin(inm,GmfTriangles,&iv[0],&iv[1],&iv[2],&lab);
	  for (int j=0;j<3;++j)
	    iv[j]--;
	  this->elements[i].set(this->vertices,iv,lab);
	  mes += this->elements[i].mesure();
	}

    }


  /* read mesh segement */
  mesb=0;
  GmfGotoKwd(inm,GmfEdges);
  for(int i=0;i<nbe;++i)
    {
      GmfGetLin(inm,GmfEdges,&iv[0],&iv[1],&lab);
      assert( iv[0]>0 && iv[0]<=nv && iv[1]>0 && iv[1]<=nv);
      for (int j=0;j<2;j++) iv[j]--;
      this->borderelements[i].set(vertices,iv,lab);
      mesb += this->borderelements[i].mesure();
    }

  GmfCloseMesh(inm);
    if(verbosity>1)  cout << "   mesure :  "<< mes << " border mesure : " << mesb<< endl;
  return(0); // OK

}

int Mesh2::Save(const string & filename)
{
  int ver = GmfDouble, outm;
  if ( !(outm = GmfOpenMesh(filename.c_str(),GmfWrite,ver,2)) ) {
    cerr <<"  -- Mesh2::Save  UNABLE TO OPEN  :"<< filename << endl;
    return(1);
  }
  double fx,fy;
  GmfSetKwd(outm,GmfVertices,this->nv);
  for (int k=0; k<nv; k++) {
    const  Vertex & P = this->vertices[k];
    GmfSetLin(outm,GmfVertices,fx=P.x,fy=P.y,P.lab);
  }

  GmfSetKwd(outm,GmfTriangles,nt);
  for (int k=0; k<nt; k++) {
    const Element & K(this->elements[k]);
    int i0=this->operator()(K[0])+1;
    int i1=this->operator()(K[1])+1;
    int i2=this->operator()(K[2])+1;
    int lab=K.lab;
    GmfSetLin(outm,GmfTriangles,i0,i1,i2,lab);
  }

  GmfSetKwd(outm,GmfEdges,nbe);
  for (int k=0; k<nbe; k++) {
    const BorderElement & K(this->borderelements[k]);
    int i0=this->operator()(K[0])+1;
    int i1=this->operator()(K[1])+1;
    int lab=K.lab;
    GmfSetLin(outm,GmfEdges,i0,i1,lab);
  }

  GmfCloseMesh(outm);
  return (0);

}

const     string Gsbegin="Mesh2::GSave v0",Gsend="end";

template<class Mesh>
void GSave2(FILE * ff,const Mesh & Th)
    {
	PlotStream f(ff);

	f <<  Gsbegin ;
	int nbe=Th.nbBrdElmts();
	f << Th.nv << Th.nt << nbe;
	for (int k=0; k<Th.nv; k++) {
	    const  typename Mesh::Vertex & P = Th(k);
	    f << P.x <<P.y  << P.lab ;
	}

	    for (int k=0; k<Th.nt; k++) {
		const typename Mesh::Element & K(Th[k]);
		int i0=Th(K[0]);
		int i1=Th(K[1]);
		int i2=Th(K[2]);

		int lab=K.lab;
		f << i0 << i1 << i2  << lab;
	    }



	for (int k=0; k<nbe; k++) {
	    const typename Mesh::BorderElement & K(Th.be(k));
	    int i0=Th(K[0]);
	    int i1=Th(K[1]);
	    int lab=K.lab;
	    f << i0 << i1  << lab;
	}
	f << Gsend;
    }

 template   void GSave2<Mesh>(FILE * ff,const Mesh & Th) ;


    Mesh2::Mesh2(FILE *f)
    {

	GRead(f);
	assert( (nt >= 0 || nbe>=0)  && nv>0) ;
	BuildBound();
	if(verbosity>1)
	    cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;

	if(nt > 0){
	    BuildAdj();
	    Buildbnormalv();
	    BuildjElementConteningVertex();
	}

	if(verbosity>1)
	    cout << "  -- Mesh2  (File *), d "<< 2  << ", n Tet " << nt << ", n Vtx "
	    << nv << " n Bord " << nbe << endl;

    }

    void Mesh2::GRead(FILE * ff)
    {
	PlotStream f(ff);
	string s;
	f >> s;
	ffassert( s== Gsbegin);
	f >> nv >> nt >> nbe;
	if(verbosity>1)
	    cout << " GRead : nv " << nv << " " << nt << " " << nbe << endl;
	this->vertices = new Vertex[nv];
	this->elements = new Element [nt];
	this->borderelements = new BorderElement[nbe];
	for (int k=0; k<nv; k++) {
	    Vertex & P = this->vertices[k];
	    f >> P.x >>P.y >> P.lab ;
	}
	mes=0.;
	mesb=0.;

	if(nt != 0)
	  {

	      for (int k=0; k<nt; k++) {
		  int i[4],lab;
		  Element & K(this->elements[k]);
		  f >> i[0] >> i[1] >> i[2] >> lab;
		  K.set(this->vertices,i,lab);
		  mes += K.mesure();

	      }
	  }


	for (int k=0; k<nbe; k++) {
	    int i[4],lab;
	    BorderElement & K(this->borderelements[k]);
	    f >> i[0] >> i[1]   >> lab;
	    K.set(this->vertices,i,lab);
	    mesb += K.mesure();

	}
	f >> s;
	ffassert( s== Gsend);
    }

Mesh2::Mesh2(int nnv, int nnt, int nnbe, Vertex2 *vv, Triangle2 *tt, BoundaryEdge2 *bb)
{

	nv = nnv;
	nt = nnt;
	nbe =nnbe;

	vertices = vv;
	elements = tt;
	borderelements = bb;

	mes=0.;
	mesb=0.;

	for (int i=0;i<nt;i++)
	  mes += this->elements[i].mesure();

	for (int i=0;i<nbe;i++)
	  mesb += this->be(i).mesure();


  if(verbosity>1)
  cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;

}

}
