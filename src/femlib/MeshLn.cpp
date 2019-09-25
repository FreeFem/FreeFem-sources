// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model  mesh 3dL
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
#include <cstring>
#include "libmesh5.h"
#include "ufunction.hpp"
#include "error.hpp"
#include "RNM.hpp"
namespace Fem2D
{
}
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "MeshSn.hpp"
#include "MeshLn.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "PlotStream.hpp"



namespace Fem2D
{
   
    
  template<> const int (* const GenericElement<DataPoint3>::nvface)[3] = 0 ;
  template<> const int (* const GenericElement<DataPoint3>::nvedge)[2] = 0 ;
    
  template<> int  GenericMesh<EdgeL,BoundaryPointL,Vertex3>::kfind=0;
  template<> int  GenericMesh<EdgeL,BoundaryPointL,Vertex3>::kthrough=0;
    
   
  /*   void Add(int *p,int n,int o) {
       for(int i=0;i<n;++i)
       p[i] += o;
       }
  */
    
  const string GsbeginL="MeshS::GSave v0",GsendL="end";
  void MeshL::GSave(FILE * ff,int offset) const
  {
    PlotStream f(ff);
        
    f <<  GsbeginL ;
    f << nv << nt << nbe;
    for (int k=0; k<nv; k++) {
      const  Vertex & P = this->vertices[k];
      f << P.x <<P.y << P.z << P.lab ;
    }
        
    for (int k=0; k<nt; ++k) {
      int iv[EdgeL::nv];
      const Element & K(this->elements[k]);
      for(int i=0;i<EdgeL::nv;++i)
	iv[i]=this->operator()(K[i])+offset;
      int lab=K.lab;
      f << iv[0] << iv[1] << lab;
    }
    for (int k=0; k<nbe; k++) {
      const BorderElement & K(this->borderelements[k]);
      int iv=this->operator()(K[0])+offset;
      int lab=K.lab;
      f << iv << lab;
    }
    f << GsendL;
  }
    
    
    
  void  MeshL::read(istream &f)
  { // read the mesh
    int i;
    //    f >> nv >> nt >> neb ;
    string str;
    int err=0;
    while(!f.eof())
      {
	f >> str;
	//cout << str << endl;
	if( str== "Vertices")
	  {
	    f >> nv;
	    assert(!this->vertices );
	    if(verbosity>2)
	      cout << "  -- Nb of Vertex " << nv << endl;
	    this->vertices = new Vertex[nv];
	    for (i=0;i<nv;i++)
	      {
		f >> this->vertices[i];
		assert(f.good());
	      }
	  }
	else if (str=="Edges")
	  {
	    f >> nt;
	    assert(this->vertices && !this->elements);
	    this->elements = new Element [nt];
	    mes=0;
	    assert(this->elements);
	    if(verbosity>2)
	      cout <<   "  -- Nb of Elements " << nt << endl;
	    for (int i=0;i<nt;i++)
	      {
		this->t(i).Read1(f,this->vertices,this->nv);
		if(this->t(i).mesure()<=0) err++; // Modif FH nov 2014
		mes += this->t(i).mesure();
	      }
	  }
	else if (str=="Points")
	  {
	    mesb=0;
	    int kmv=0,ij;
	    f >> nbe;
	    assert(vertices);
	    this->borderelements = new BorderElement[nbe];
	    if(verbosity>2)
	      cout <<   "  -- Nb of border Triangles " << nbe << endl;
	    for (i=0;i<nbe;i++)
	      {
		this->be(i).Read1(f,this->vertices,this->nv);
		//mesb += this->be(i).mesure();
		for(int j=0;j<BorderElement::nv;++j)
		  if(!vertices[ij=this->be(i,j)].lab)
		    {
		      vertices[ij].lab=1;
		      kmv++;
		    }
	      }
	  }
	else if(str[0]=='#')
	  {// on mange la ligne
	    int c;
	    while ( (c=f.get()) != '\n' &&  c != EOF)
	      //cout << c << endl;
	      ;
	  }
      }
    assert( (nt >= 0 || nbe>=0)  && nv>0) ;
    /*   done at up level ...
         BuildBound();
         
         if(nt > 0){
         BuildAdj();
         Buildbnormalv();
         BuildjElementConteningVertex();
         }
    */
    if(err!=0)
      {
	cerr << " MeshS::read: sorry bad mesh. Number of negative Triangles " << err << endl;
	this->~MeshL();
	ffassert(0);
      }
  }
    
    
  void MeshL::readmsh(ifstream & f,int offset)
  {
    int err=0;
    f >> nv >> nt >> nbe;
    if(verbosity>2)
      cout << " GRead : nv " << nv << " " << nt << " " << nbe << endl;
    this->vertices = new Vertex[nv];
    this->elements = new Element[nt];
    this->borderelements = new BorderElement[nbe];
    for (int k=0; k<nv; k++) {
      Vertex & P = this->vertices[k];
      f >> P.x >>P.y >> P.z >> P.lab ;
    }
    mes=0.;
    mesb=0.;
        
    if(nt == 0) {
      cerr << "  A meshS type must have elements  " << endl;
      ffassert(0);exit(1);
            
    }
        
    for (int k=0; k<nt; k++) {
      int i[2],lab;
      Element & K(this->elements[k]);
      f >> i[0] >> i[1] >> lab;
      Add(i,2,offset);
      K.set(this->vertices,i,lab);
      mes += K.mesure();
      err += K.mesure() <0;
                
    }
        
        
        
    for (int k=0; k<nbe; k++) {
      int i[1],lab;
      BorderElement & K(this->borderelements[0]);
      f >> i[0] >> lab;
      Add(i,1,offset);
      K.set(this->vertices,i,lab);
      mesb += K.mesure();
            
    }
    if(err!=0)
      {
	cerr << " MeshS::readmsh : sorry bad mesh. Number of negative Tri " << err << endl;
	this->~MeshL();
	ffassert(0);
      }
        
  }
    
    
  MeshL::MeshL(const string filename)
    :mapSurf2Curv(0),mapCurv2Surf(0)    {
    int ok=load(filename);
    if(verbosity) {
      cout << "read meshL ok " << ok ;
      cout << "surface Mesh, num Edges:= " << nt << ", num Vertice:= " << nv << " num boundary Points:= " << nbe << endl;
    }
        
    if (ok) {
      ifstream f(filename.c_str());
      if(!f) {
	cerr << "  --  MeshL::MeshL Erreur openning " << filename<<endl;ffassert(0);exit(1);}
      if(verbosity>2)
	cout << "  -- MeshL:  Read On file \"" <<filename<<"\""<<  endl;
      if(filename.rfind(".msh")==filename.length()-4)
	readmsh(f,-1);
      else
	read(f);
    }
        
    BuildBound();
    BuildAdj();
    Buildbnormalv();
    BuildjElementConteningVertex();
        
        
    if(verbosity>2) cout << "  -- End of read: MeshS mesure = " << mes << " border mesure " << mesb << endl;
        
    if(verbosity) cout << "  -- MeshL : "<<filename  << ", space dimension "<< 3  << ", num Edges elts " << nt << ", num Vertice "
		       << nv << " num Bondary Points " << nbe << endl;
        
        
    ffassert(mes>=0);
        
  }
    
    
    
    
  int MeshL::load(const string & filename)
  {
    int bin;
    int ver,inm,dim;
    int lf=filename.size()+20;
    KN<char>  fileb(lf),filef(lf);
    char *data = new char[filename.size()+1];
    size_t ssize = filename.size()+1;
    char *ptr;
    char *pfile=data;
    strncpy( data, filename.c_str(),ssize);
    ptr = strstr(data,".mesh");
    if( !ptr ){
      strcpy(filef,filename.c_str());
      strcpy(fileb,filef);
      strcat(filef,".mesh");
      strcat(fileb,".meshb");
      if( (inm=GmfOpenMesh(pfile=fileb, GmfRead,&ver,&dim)) )
	bin=true;
      else if( (inm=GmfOpenMesh(pfile=filef, GmfRead,&ver,&dim)) )
	bin=false;
      else
	if(verbosity>5){
	  cerr << " Erreur ouverture file " << (char *) fileb  << " " << (char *) filef  <<endl;
	  return   1;
	}
    }
    else{
      if( !(inm=GmfOpenMesh(data, GmfRead,&ver,&dim)) ){
	if(verbosity>5)
	  cerr << " Erreur ouverture file " << (char *) data  << endl;
	return   1;
      }
    }
    // data file is readed and the meshes are initilized
    int nv=-1,nEdge=-1,nPts=-1;
    nv=GmfStatKwd(inm,GmfVertices);  // vertice
    nEdge=GmfStatKwd(inm,GmfEdges); // segment elements
    ///  nPts=GmfStatKwd(inm,GmfVertices);  // points border element
    this->set(nv,nEdge,nPts);
       
    if(nEdge == 0) {
      cerr << "  A meshL type must have elements  " << endl;
      ffassert(0);exit(1);}
        
    if(verbosity>1)
      cout << "  -- MeshL(load): "<< (char *) data <<  ", MeshVersionFormatted:= " << ver << ", space dimension:= "<< dim
	   << ", num Edge elts:= " << nEdge << ", num vertice:= " << nv << " num Points boundaries:= " << nPts << endl;
        
    if(dim  != 3) {
      cerr << "Err dim == " << dim << " !=3 " <<endl;
      return 2; }
    if( nv<=0 && (nEdge <=0 || nPts <0) ) {
      cerr << " missing data "<< endl;
      return 3;
    }
    int iv[3],lab;
    float cr[3];
    int mxlab=0, mnlab=0;
    // read vertices
    GmfGotoKwd(inm,GmfVertices);
    for(int i=0;i<this->nv;++i) {
      if(ver<2) {
	GmfGetLin(inm,GmfVertices,&cr[0],&cr[1],&cr[2],&lab);
	vertices[i].x=cr[0];
	vertices[i].y=cr[1];
	vertices[i].z=cr[2];}
      else
	GmfGetLin(inm,GmfVertices,&vertices[i].x,&vertices[i].y,&vertices[i].z,&lab);
      vertices[i].lab=lab;
      mxlab= max(mxlab,lab);
      mnlab= min(mnlab,lab);
    }
    // read triangles (element meshS)
    if(mnlab==0 && mxlab==0 ) {
      int kmv=0;
      mes=0;
      GmfGotoKwd(inm,GmfEdges);
      for(int i=0;i<nEdge;++i) {
	GmfGetLin(inm,GmfEdges,&iv[0],&iv[1],&lab);
	assert( iv[0]>0 && iv[0]<=nv && iv[1]>0 && iv[1]<=nv);
	for(int j=0;j<2;++j)
	  if(!vertices[iv[j]-1].lab) {
	    vertices[iv[j]-1].lab=1;
	    kmv++;
	  }
	for (int j=0;j<2;++j) iv[j]--;
	elements[i].set(vertices,iv,lab);
	mes += elements[i].mesure();
      }
      if(kmv&& verbosity>1) cout << "    Aucun label Hack (FH)  ??? => 1 sur les triangle frontiere "<<endl;
    }
    else {
      mes=0;
      GmfGotoKwd(inm,GmfEdges);
      for(int i=0;i<nEdge;++i) {
	GmfGetLin(inm,GmfEdges,&iv[0],&iv[1],&lab);
	assert( iv[0]>0 && iv[0]<=nv && iv[1]>0 && iv[1]<=nv);
	for (int j=0;j<2;++j) iv[j]--;
	elements[i].set(this->vertices,iv,lab);
	mes += elements[i].mesure();
      }
    }
    // read edges (boundary elements meshS)
    mesb=0;
    GmfGotoKwd(inm,GmfVertices);
    for(int i=0;i<nPts;++i) {
      GmfGetLin(inm,GmfVertices,&iv[0],&lab);
      assert( iv[0]>0 && iv[0]<=nv);
      borderelements[i].set(this->vertices,iv,lab);   //element
      mesb += this->borderelements[i].mesure();
    }
        
        
    if(verbosity>1)
      cout << "  -- MeshL(load): "<< (char *) data <<  ", MeshVersionFormatted:= " << ver << ", space dimension:= "<< dim
	   << ", Edges elts:= " << nt << ", num vertice:= " << nv << ", num Points boundaries:= " << nbe << endl;
        
    GmfCloseMesh(inm);
    delete[] data;
    return 0; // OK
  }
    
    
  MeshL::MeshL(const string filename, const long change)
    :mapSurf2Curv(0),mapCurv2Surf(0)  {
        
        
    int ok=load(filename);
    if(verbosity) {
      cout << "read mesh ok " << ok  << endl;
      cout << ", nt " << nt << ", nv " << nv << " nbe:  = " << nbe << endl;
    }
    if(ok)
      {
	ifstream f(filename.c_str());
	if(!f) {
	  cerr << "  --  Mesh3::Mesh3 Erreur openning " << filename<<endl;ffassert(0);exit(1);}
	if(verbosity>2)
	  cout << "  -- Mesh3:  Read On file \"" <<filename<<"\""<<  endl;
	if(filename.rfind(".msh")==filename.length()-4)
	  readmsh(f,-1);
	else
	  read(f);
      }
        
    if(change){
      // verification multiple points
      double hseuil=hmin();
      hseuil = hseuil/10;
      cout << " hseuil = " << hseuil << endl;
      KN<int> Numero_Som(this->nv);
      Vertex *vv=new Vertex[this->nv];
      int nv_t=0;
      {
	R3 Pinf(1e100,1e100,1e100),Psup(-1e100,-1e100,-1e100);
	for (int ii=0;ii< this->nv;ii++){
	  R3 P( vertices[ii].x, vertices[ii].y, vertices[ii].z);
	  Pinf=Minc(P,Pinf);
	  Psup=Maxc(P,Psup);
	}
	EF23::GTree<Vertex3> *gtree = new EF23::GTree<Vertex3>(vv,Pinf,Psup,0);
	// creation of octree
	for (int ii=0;ii<this->nv;ii++){
	  const R3 r3vi( this->vertices[ii].x, this->vertices[ii].y, this->vertices[ii].z );
	  const Vertex3 &vi(r3vi);
                    
	  Vertex3 * pvi=gtree->ToClose(vi,hseuil);
                    
	  if(!pvi){
	    vv[nv_t].x = vi.x;
	    vv[nv_t].y = vi.y;
	    vv[nv_t].z = vi.z;
	    vv[nv_t].lab = this->vertices[ii].lab; // lab mis a zero par default
	    Numero_Som[ii] = nv_t;
	    gtree->Add( vv[nv_t] );
	    nv_t=nv_t+1;
	  }
	  else{
	    Numero_Som[ii] = pvi-vv;
	  }
	}
                
	delete gtree;
	//delete [] vv;
      }
            
      // general case
            
      KN<int> takevertex(nv_t,0);
      for (int k=0; k<nbe; k++) {
	const BorderElement & K(this->borderelements[k]);
	for(int jj=0; jj<BorderElement::nv; jj++){
	  takevertex[ Numero_Som[this->operator()(K[jj])] ] = 1;
	}
      }
      for(int k=0; k< this->nt; k++){
	const Element & K(this->elements[k]);
	for(int jj=0; jj<Element::nv; jj++){
	  takevertex[ Numero_Som[this->operator()(K[jj])] ] = 1;
	}
      }
            
      int newvertex=0;
      for(int iv=0; iv<nv_t; iv++){
	newvertex+=takevertex[iv];
      }
            
      if( newvertex != this->nv){
                
	// determination of vertex
	Vertex *vvv = new Vertex[ newvertex ];
	KN<int> newNumero_Som(nv_t);
	int iii=0;
	for(int iv=0;  iv< nv_t; iv++){
	  if( takevertex[iv ] == 1  ){
	    vvv[iii].x = vv[iv].x;
	    vvv[iii].y = vv[iv].y;
	    vvv[iii].z = vv[iv].z;
	    vvv[iii].lab = vv[iv].lab; // lab mis a zero par default
	    newNumero_Som[iv] = iii;
	    iii++;
	  }
	}
	ffassert( newvertex== iii );
                
	Element *tt=new Element[this->nt];
	BorderElement *bb = new BorderElement[this->nbe];
                
	Element *ttt=tt;
	BorderElement *bbb=bb;
                
	for (int k=0; k<this->nbe; k++) {
	  const BorderElement & K(this->borderelements[k]);
	  int iv[2];
	  for(int jj=0; jj<2; jj++){
	    iv[jj] = Numero_Som[this->operator()(K[jj])];
	    iv[jj] = newNumero_Som[iv[jj]];
	  }
	  (bbb++)->set(vvv,iv,K.lab);
	}
                
	for(int k=0; k< this->nt; k++){
	  const Element & K(this->elements[k]);
	  int iv[3];
	  for(int jj=0; jj<3; jj++){
	    iv[jj] = Numero_Som[this->operator()(K[jj])];
	    iv[jj] = newNumero_Som[iv[jj]];
	  }
	  (ttt++)->set(vvv,iv,K.lab);
	}
	cout << " delete vertices + autre " << endl;
	delete [] vertices;
	delete [] elements;
	delete [] borderelements;
                
	nv = newvertex;
                
	vertices = vvv;
	elements = tt;
	borderelements = bb;
                
	//&this = new Mesh3(newvertex,this->nt,this-nbe,vvv,tt,bb);
                
	delete [] newNumero_Som;
      }
      else{
	cout << " no need to change the mesh " << endl;
      }
      delete [] Numero_Som;
    }
        
    BuildBound();
    BuildAdj();
    Buildbnormalv();
    BuildjElementConteningVertex();
        
    if(verbosity>2)
      cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;
    if(verbosity)
      cout << "  -- MeshS : "<<filename  << ", d "<< 3  << ", n Tri " << nt << ", n Vtx "
	   << nv << " n Bord Edges " << nbe << endl;
    ffassert(mes>=0); // add F. Hecht sep 2009.
  }
    
 
  /*
    MeshS::MeshS(FILE *f,int offset)
    :mapSurf2Vol(0),mapVol2Surf(0)
    {
    GRead(f,offset);// remove 1
    assert( (nt >= 0 || nbe>=0)  && nv>0) ;
    BuildBound();
    if(verbosity>2)
    cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;
      
    BuildAdj();
    Buildbnormalv();
    BuildjElementConteningVertex();
        
    // if not edges then build the edges - need access to the old adjacensce to build eges and rebuild the new adj
    if (nbe==0) {
    BuildEdges();
    delete [] TheAdjacencesLink;
    delete [] BoundaryElementHeadLink;
    TheAdjacencesLink=0;
    BoundaryElementHeadLink=0;
    BuildBound();
    BuildAdj();
    Buildbnormalv();
    BuildjElementConteningVertex();
    }
        
    if(verbosity>2)
    cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;
        
    if(verbosity>1)
    cout << "  -- MeshS  (File *), d "<< 3  << ", n Tri " << nt << ", n Vtx "
    << nv << " n Bord " << nbe << endl;
    ffassert(mes>=0); // add F. Hecht sep 2009.
    }
  */
    
  double MeshL::hmin() const
  {
    R3 Pinf(1e100,1e100,1e100),Psup(-1e100,-1e100,-1e100);   // Extremite de la boite englobante
    double hmin=1e10;
        
    for (int ii=0;ii< this->nv;ii++) {
      R3 P( vertices[ii].x, vertices[ii].y, vertices[ii].z);
      Pinf=Minc(P,Pinf);
      Psup=Maxc(P,Psup);
    }
        
    for (int k=0;k<this->nt;k++) {
        
      if( this->elements[k].lenEdge(0) < Norme2(Psup-Pinf)/1e9 ){
	const EdgeL & K(this->elements[k]);
	int iv[3];
	for(int jj=0; jj <2; jj++)
	  iv[jj] = this->operator()(K[jj]);
	if(verbosity>2)
	  cout << "EdgeL: " << k << " lenght "<<  this->elements[k].lenEdge(0) << endl;
	if(verbosity>2) cout << " A triangleS with a very small edge was created " << endl;
	return 1;
      }
      hmin=min(hmin,this->elements[k].lenEdge(0));   // calcul de .lenEdge pour un Mesh3
          
    }
    ffassert(hmin>Norme2(Psup-Pinf)/1e9);
    return hmin;
  }
    
    
  // brute force method
  void MeshL::GRead(FILE * ff,int offset)
  {
    PlotStream f(ff);
    string s;
    f >> s;
    ffassert( s== GsbeginL);
    f >> nv >> nt >> nbe;
    if(verbosity>2)
      cout << " GRead : nv " << nv << " " << nt << " " << nbe << endl;
    this->vertices = new Vertex[nv];
    this->elements = new Element [nt];
    this->borderelements = new BorderElement[nbe];
    for (int k=0; k<nv; k++) {
      Vertex & P = this->vertices[k];
      f >> P.x >>P.y >> P.z >> P.lab ;
    }
    mes=0.;
    mesb=0.;
        
    if(nt == 0) {
      cerr << "  A meshL type must have elements  " << endl;
      ffassert(0);exit(1);}
 
            
    for (int k=0; k<nt; k++) {
      int i[2],lab;
      Element & K(this->elements[k]);
      f >> i[0] >> i[1] >> lab;
      Add(i,2,offset);
      K.set(this->vertices,i,lab);
      mes += K.mesure();
    
    }
    for (int k=0; k<nbe; k++) {
      int i[2],lab;
      BorderElement & K(this->borderelements[k]);
      f >> i[0] >> lab;
      Add(i,1,offset);
      K.set(this->vertices,i,lab);
      mesb += K.mesure();
            
    }
    f >> s;
    ffassert( s== GsendL);
  }
    
    
  const MeshL::Element * MeshL::Find( Rd P, R1 & Phat,bool & outside,const Element * tstart) const
    
  {
    ffassert(0);
  }
    
    
    
  MeshL::MeshL(int nnv, int nnt, int nnbe, Vertex3 *vv, EdgeL *tt, BoundaryPointL *bb, bool cleanmesh, bool removeduplicate, bool rebuildboundary, int orientation, double precis_mesh)
    :mapSurf2Curv(0),mapCurv2Surf(0)
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
    if (cleanmesh) {
      if(verbosity>3)
	cout << "before clean meshL, nv: " <<nv << " nt:" << nt << " nbe:" << nbe << endl;
      clean_mesh(precis_mesh, nv, nt, nbe, vertices, elements, borderelements, removeduplicate, rebuildboundary, orientation);
      if(verbosity>3)
	cout << "after clean meshL, nv: " <<nv << " nt:" << nt << " nbe:" << nbe << endl;
    }
    BuildBound();
    BuildAdj();
    Buildbnormalv();
    BuildjElementConteningVertex();
        
    if(verbosity>1)
      cout << "  -- End of read meshL: mesure = " << mes << " border mesure " << mesb << endl;
        
    assert(mes>=0.);
    assert(mesb==nbe);
  }
   
    
    
    
  int MeshL::Save(const string & filename) const
  {
    int ver = GmfDouble, outm;
    if ( !(outm = GmfOpenMesh(filename.c_str(),GmfWrite,ver,3)) ) {
      cerr <<"  -- MeshS**::Save  UNABLE TO OPEN  :"<< filename << endl;
      return(1);
    }
    float fx,fy,fz;
    // write vertice (meshL)
    GmfSetKwd(outm,GmfVertices,nv);
    for (int k=0; k<nv; k++) {
      const  Vertex & P = vertices[k];
      GmfSetLin(outm,GmfVertices,fx=P.x,fy=P.y,fz=P.z,P.lab);
    }
    // write triangles (meshS)
    GmfSetKwd(outm,GmfEdges,nt);
    for (int k=0; k<nt; k++) {
      const EdgeL & K(elements[k]);
      int i0=this->operator()(K[0])+1;
      int i1=this->operator()(K[1])+1;
      int lab=K.lab;
      GmfSetLin(outm,GmfEdges,i0,i1,lab);
    }
        
    // no boundary points ?
   
    GmfCloseMesh(outm);
    return (0);
  }
    
    
}
