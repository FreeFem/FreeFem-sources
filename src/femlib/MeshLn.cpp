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
  static const int  nvfaceSeg[1][3]  = {{-1,-1,1}};
  static const int  nvedgeSeg[1][2] = { {0,1} };
  static const int  nvadjSeg[2][1] = { {0},{1} };
    
  // geometry element for segment ( boundary elements in surface mesh, Rd=3 RdHat=1 )
  template<> const int (* const GenericElement<DataSeg3>::nvface)[3] = 0 ;
  template<> const int (* const GenericElement<DataSeg3>::nvedge)[2] = nvedgeSeg; //nvedgeTria ;
  template<> const int (* const GenericElement<DataSeg3>::nvadj)[1] = nvadjSeg ;
    
    
    
  template<> const int (* const GenericElement<DataPoint3>::nvface)[3] = 0 ;
  template<> const int (* const GenericElement<DataPoint3>::nvedge)[2] = 0 ;
  template<> const int (* const GenericElement<DataPoint3>::nvadj)[1] = 0 ;
    
  template<> int  GenericMesh<EdgeL,BoundaryPointL,Vertex3>::kfind=0;
  template<> int  GenericMesh<EdgeL,BoundaryPointL,Vertex3>::kthrough=0;
    
  static const int onWhatIsVertex[2][3] = {  {1,0,0}, // sommet  0
        {0,1,0}}; // sommet 1
    
  template<>
  const int (* const GenericElement<DataSeg3>::onWhatBorder)[3] = onWhatIsVertex ;
  template<> const int  GenericElement<DataSeg3>::nitemdim[4] = {2,1,0,0 }  ;
   
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
    // np present for the moment
	else if (str=="Points")
	  {
	    mesb=0;
        cout <<   "  -- No boundary points in ff .mesh format  " << endl;
	    /*int kmv=0,ij;
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
	      }*/
	  }
	else if(str[0]=='#') {
	    int c;
	    while ( (c=f.get()) != '\n' &&  c != EOF)
	      ;
	  }
      }
    assert( (nt >= 0 || nbe>=0)  && nv>0) ;
    if(err!=0) {
	cerr << " MeshL::read: sorry bad mesh. Number of negative Edges " << err << endl;
	this->~MeshL();
	ffassert(0);
      }
  }
    
  // no points border in this format
  void MeshL::readmsh(ifstream & f,int offset)
  {
    int err=0;
      f >> nv >> nt ; //>> nbe;
    if(verbosity>2)
        cout << " GRead : nv " << nv << " " << nt << endl; //" " << nbe << endl;
    this->vertices = new Vertex[nv];
    this->elements = new Element[nt];
    //this->borderelements = new BorderElement[nbe];
    for (int k=0; k<nv; k++) {
      Vertex & P = this->vertices[k];
      f >> P.x >>P.y >> P.z >> P.lab ;
    }
    mes=0.;
    mesb=0.;
        
    if(nt == 0) {
      cerr << "  A meshL type must have elements  " << endl;
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
    cout <<   "  -- No boundary points in .msh format  " << endl;
    /*for (int k=0; k<nbe; k++) {
      int i[1],lab;
      BorderElement & K(this->borderelements[k]);
      f >> i[0] >> lab;
      Add(i,1,offset);
      K.set(this->vertices,i,lab);
      mesb += K.mesure();
            
    }*/
    if(err!=0)
      {
	cerr << " MeshL::readmsh : sorry bad mesh. Number of negative Edges " << err << endl;
	this->~MeshL();
	ffassert(0);
      }
        
  }
    
    
  MeshL::MeshL(const string filename)
    :mapSurf2Curv(0),mapCurv2Surf(0)    {
    int ok=load(filename);
    if(verbosity) {
      cout << "read meshL ok " << ok ;
      cout << "curve Mesh, num element Edges:= " << nt << ", num Vertice:= " << nv << " num boundary Points:= " << nbe << endl;
    }
        
    if (ok) {
      ifstream f(filename.c_str());
      if(!f) {
	cerr << "  --  MeshL::MeshL Erreur opening " << filename<<endl;ffassert(0);exit(1);}
      if(verbosity>2)
	cout << "  -- MeshL:  Read On file \"" <<filename<<"\""<<  endl;
      if(filename.rfind(".msh")==filename.length()-4)
	readmsh(f,-1);
      else
	read(f);
    }
        
    BuildBound();
    BuildAdj();
    //Buildbnormalv();
    BuildjElementConteningVertex();
        
        
    if(verbosity>2) cout << "  -- End of read: MeshL mesure = " << mes << endl;
        
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
    nPts=0; // GmfStatKwd(inm,GmfVertices);  // points border element
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
    // unused loop...not border points in .mesh
    mesb=0;
    GmfGotoKwd(inm,GmfVertices);
    for(int i=0;i<nPts;++i) {
        GmfGetLin(inm,GmfVertices,&iv[0],&lab);
        assert( iv[0]>0 && iv[0]<=nv);
        borderelements[i].set(this->vertices,iv,lab);
        mesb += this->borderelements[i].mesure();
    }
        
        
    if(verbosity>1)
        cout << "  -- MeshL(load): "<< (char *) data <<  ", MeshVersionFormatted:= " << ver << ", space dimension:= "<< dim
	    << ", Edges elts:= " << nt << ", num vertice:= " << nv << ", num Points boundaries:= " << nbe << endl;
        
    GmfCloseMesh(inm);
    delete[] data;
    return 0; // OK
  }
    
    
  MeshL::MeshL(const string filename, bool cleanmesh, bool removeduplicate, bool rebuildboundary, int orientation, double precis_mesh, double ridgeangledetection)
    :mapSurf2Curv(0),mapCurv2Surf(0)  {
        
        
    int ok=load(filename);
    if(verbosity) {
      cout << "read meshL ok " << ok  << endl;
      cout << ", nt " << nt << ", nv " << nv << " nbe:  = " << nbe << endl;
    }
    if(ok) {
        ifstream f(filename.c_str());
        if(!f)
            cerr << "  --  MeshL Erreur opening " << filename<<endl;ffassert(0);exit(1);
        if(verbosity>2)
                cout << "  -- MeshL:  Read On file \"" <<filename<<"\""<<  endl;
        if(filename.rfind(".msh")==filename.length()-4)
            readmsh(f,-1);
        else
            read(f);
    }
     
    if (cleanmesh) {
        if(verbosity>3)
            cout << "before clean meshL, nv: " <<nv << " nt:" << nt << " nbe:" << nbe << endl;
        clean_mesh(precis_mesh, nv, nt, nbe, vertices, elements, borderelements, removeduplicate, rebuildboundary, orientation);
        if(verbosity>3)
            cout << "after clean meshL, nv: " <<nv << " nt:" << nt << " nbe:" << nbe << endl;
    }
        
    BuildBound();
    BuildAdj();
    //Buildbnormalv();
    BuildjElementConteningVertex();
        
    if(verbosity>2)
      cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;
    if(verbosity)
      cout << "  -- MeshL : "<<filename  << ", d "<< 3  << ", n Edges " << nt << ", n Vtx "
	   << nv << " n Border points " << nbe << endl;
    ffassert(mes>=0); // add F. Hecht sep 2009.
  }
    
 
  
  MeshL::MeshL(FILE *f,int offset)
  :mapSurf2Curv(0),mapCurv2Surf(0)     {
    GRead(f,offset);// remove 1
    assert( (nt >= 0 || nbe>=0)  && nv>0) ;
    BuildBound();
    if(verbosity>2)
    cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;
      
    BuildAdj();
    Buildbnormalv();
    BuildjElementConteningVertex();
      
    if(verbosity>2)
        cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;
        
    if(verbosity>1)
        cout << "  -- MeshL  (File *), d "<< 3  << ", n Tri " << nt << ", n Vtx " << nv << " n Bord " << nbe << endl;
    ffassert(mes>=0); // add F. Hecht sep 2009.
  }
  
    
  double MeshL::hmin() const {
    R3 Pinf(1e100,1e100,1e100),Psup(-1e100,-1e100,-1e100);   // Extremite de la boite englobante
    double hmin=1e10;
        
    for (int ii=0;ii< this->nv;ii++) {
      R3 P( vertices[ii].x, vertices[ii].y, vertices[ii].z);
      Pinf=Minc(P,Pinf);
      Psup=Maxc(P,Psup);
    }
        
    for (int k=0;k<this->nt;k++) {
        
        if( this->elements[k].mesure() < Norme2(Psup-Pinf)/1e9 ) {
            const EdgeL & K(this->elements[k]);
            int iv[2];
            for(int jj=0; jj <2; jj++)
                iv[jj] = this->operator()(K[jj]);
            if(verbosity>2)
                cout << "EdgeL: " << k << " lenght "<<  this->elements[k].mesure() << endl;
            if(verbosity>2) cout << " A triangleS with a very small edge was created " << endl;
            return 1;
        }
        hmin=min(hmin,this->elements[k].mesure());   // calcul de .lenEdge pour un Mesh3
          
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
      // rewritie FH 31 jan 2020 ..
      static int count =0;
      if( verbosity && count++< 5  )
          cerr << " MeshL::Find warning brute force to day " << endl;
      //  find the neast points ..
      double dmin2 = 1e200;
      R1 Phm;
      bool out = true;
      int n=-1;
      for (int i=0;i<nt;i++) {
          kthrough++;
          const EdgeL & K(this->elements[i]);
          R3 A(K[0]),B(K[1]), AB(A,B);
          R3 AP(A,P);
          double lab2 = AB.norme2();
          double l = min(1.,max(0.,(AB,AP)/lab2));
          R1 Ph(l);
          R3 Pt=K(Ph);
          double d2=R3(P,Pt).norme2();
          if(dmin2>d2)
          {
              dmin2 = d2;
              Phm=Ph;
              n =i;
              out = d2 < lab2*1e-2; // BofBof FH ...
          }
   
      }
      if( n<0) return  0; 
      Phat=Phm;
      return this->elements+n; // outside
  }
    
    
  MeshL::MeshL(int nnv, int nnt, int nnbe, Vertex3 *vv, EdgeL *tt, BoundaryPointL *bb, bool cleanmesh, bool removeduplicate, bool rebuildboundary, int orientation, double precis_mesh, double ridgeangledetection)
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
    //Buildbnormalv();
    BuildjElementConteningVertex();
        
    if(verbosity>1)
      cout << "  -- End of read meshL: mesure = " << mes << " border mesure " << mesb << endl;
        
    assert(mes>=0.);
    assert(mesb==0.);
  }
   
 
  int MeshL::Save(const string & filename) const
  {
    int ver = GmfDouble, outm;
    if ( !(outm = GmfOpenMesh(filename.c_str(),GmfWrite,ver,3)) ) {
        cerr <<"  -- MeshL**::Save  UNABLE TO OPEN  :"<< filename << endl;
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

  
    
    // determine the bounder points list for meshL
   void MeshL::BuildBdElem(double angle) {
      
        
       delete [] borderelements; // to remove the previous pointers
       borderelements = new BoundaryPointL[2 * nt]; // 2 * nt upper evaluated
       
       HashTable<SortArray<int, 1>, int> pointI(2 * nt, nt);
       int* AdjLink = new int[2 * nt];
       
       int nbeL=0,nbiL=0,nk=0;
       // Build border points from the edge list
       for (int i = 0; i < nt; i++)
           for (int j = 0; j < 2; j++) {
               int jt = j, it = ElementAdj(i, jt);
               EdgeL &K(elements[i]);  // current element
               // True border point -> no adjacence / on domain border
               if ((it == i || it < 0)) {
                   int iv[1];
                       iv[0] = this->operator () (K [EdgeL::nvedge[0][j]]);
                   if(verbosity>15)
                       cout << " the edge " << iv[0] << " is a boundary " << endl;
                   be(nbeL++).set(vertices,iv,K.lab);
                   
               }
               // internal point -- check angular and no manifold
               else {
                   EdgeL &K_adj(elements[it]); // adjacence element
                   int iv[1];
                   iv[0] = this->operator () (K [EdgeL::nvedge[0][j]]);
                   SortArray<int, 1> key(iv[0]);
                   typename HashTable<SortArray<int,1>,int>::iterator p= pointI.find(key);
                   if (!p) {
                       //edge element
                       R3 A(K[0]),B(K[1]);
                       R3 E(B-A);
                       E/=E.norme();
                       // adj edge element
                       R3 A_adj(K_adj[0]),B_adj(K_adj[1]);
                       R3 E_adj(B_adj-A_adj);
                          E_adj/=E_adj.norme();
                       
                       R pdt = (E,E_adj); // scalar product
                       pdt = acos(pdt); // radian angle (Normal,Normal_adj)
                       if(verbosity>15)
                           cout << "Element num: " << i << " N " << E << " Element adjacent num: " << it << " E_adj " << E_adj << " angle between N N_adj = " << pdt <<endl;
                        
                       if(pdt >= angle) {
                           if(verbosity>15)
                               cout << " the edge " <<nbeL <<": [" << iv[0] << " " << iv[1] << "] is a boundary with the angular criteria" << endl;
                           int lab = min(K.lab, K_adj.lab);
                           be(nbeL).set(vertices,iv,lab);
                           pointI.add(key, nbeL++);
                       }
                   }
               }
               nk++;  // increment the total edge jump --- nt * 2
                
           }
       assert(nt*2==nk);
       delete [] AdjLink;
       // update the number of border points
       nbe = nbeL;
       if (verbosity>5)
           cout << " Building border point from meshS nbe: "<< nbeL << " nbi: " << nbiL << endl;
        
       BuildBound();
       delete []TheAdjacencesLink;
       delete [] BoundaryElementHeadLink;
       TheAdjacencesLink=0;
       BoundaryElementHeadLink=0;
       BuildAdj();
       //Buildbnormalv();
       BuildjElementConteningVertex();
       
   }

  MeshL::MeshL(const  Serialize &serialized)
  :GenericMesh<EdgeL,BoundaryPointL,Vertex3> (serialized),mapSurf2Curv(0),mapCurv2Surf(0) {
   BuildBound();
   if(verbosity>1)
       cout << "  -- End of serialized: mesure = " << mes << " border mesure " << mesb << endl;
 
   if(nt > 0){
       BuildAdj();
       //Buildbnormalv();
       BuildjElementConteningVertex();
   }
  if(verbosity>1)
       cout << "  -- MeshL  (serialized), d "<< 3  << ", n Edges " << nt << ", n Vtx "
       << nv << " n Bord " << nbe << endl;
   ffassert(mes>=0);
  }

    
    
}
