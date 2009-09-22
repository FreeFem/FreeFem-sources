// ORIG-DATE:     Dec 2007
// -*- Mode : c++ -*-
//
// SUMMARY  :  Model  mesh 3d   
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
#include "Mesh2dn.hpp"
#include "Mesh3dn.hpp"
#include "rgraph.hpp"
#include "fem.hpp"
#include "PlotStream.hpp"

namespace Fem2D 
{
  
  
  
  
  /*
    const short int v_tet_face[4][3]=  {{3,2,1},{0,2,3},{ 3,1,0},{ 0,1,2}};
    const short int a_tet_face[4][3]=  {{ 0,1,0},{ 0,2,0},{ 0,3,1},{ 1,2,1}};
    const bool  sig_tet_face[4][3]=  {{ 0,1,0},{ 1,0,1},{ 0,1,0},{ 1,0,1}};
    const short int v_tet_edge[6][2]= {{ 1,2},{1,3},{1,4},{2,3},{2,4},{3,4}}; 
    const short int fadj_tet_edge[6][2]= {{4,3},{2,4},{3,2},{4,1},{1,3},{2,1}};
    const short int op_tet_edge[6]={ 6, 5, 4, 3, 2, 1};  
  */
  
  //  Attention  nvfaceTet  donnne les faces  les 4 faces de tet telle que la
  // tel que  le tet forme des trois sommet  + l'autre sommet soit positif.
  //  donc  le  produit vectoriel des 2 aretes  (0,1) (0,2)  donne une  normale interieur.
  //  Ok pour les gradients des $\lambda_i$  
  static const int  nvfaceTet[4][3]  ={{3,2,1}, {0,2,3},{ 3,1,0},{ 0,1,2}}  ;//{ {2,1,3},{0,2,3},{1,0,3},{0,1,2} };
  static const int  nvedgeTet[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
  
  static const int  nvfaceTria[1][3]  = { {0,1,2} };
  static const int  nvedgeTria[3][2] = { {1,2},{2,0},{0,1}}; 
  
  static const int   nvfaceSeg[1][3]  = {{-1,-1,1}};
  static const int  nvedgeSeg[1][2] = { {0,1} };
  
  
  template<>
  const int (* const GenericElement<DataTriangle3>::nvface)[3] = nvfaceTria ;
  template<>
  const int (* const GenericElement<DataTriangle3>::nvedge)[2] = nvedgeTria ;
  template<>
  const int (* const GenericElement<DataTriangle3>::nvadj)[2] = nvedgeTria ;
  template<> const int  GenericElement<DataTriangle3>::nitemdim[4] = {3,3,1,0 }  ;
  static const int onWhatIsEdge3[3][7] = {
    {0,1,3, 2,0,0, 0}, // edge 0 
    {3,0,1, 0,2,0, 0},
    {1,3,0, 0,0,2, 0} };
  
  template<>
  const int (* const GenericElement<DataTriangle3>::onWhatBorder)[7] = onWhatIsEdge3 ;
  
  
  template<>
  const int (* const GenericElement<DataTet>::nvface)[3] = nvfaceTet ;
  template<>
  const int (* const GenericElement<DataTet>::nvedge)[2] = nvedgeTet ;
  template<>
  const int (* const GenericElement<DataTet>::nvadj)[3] = nvfaceTet ;
  template<> const int  GenericElement<DataTet>::nitemdim[4] = {4,6,4,1 }  ;
  
  int onWhatIsFace[4][15] ; 
    typedef const int   (*const PtrConst15int) [15]; //  a pointeur on  const arry of 15 int. (to help xcode) 
 // static const int (* const SetonWhatIsFace(int  onWhatIsFace[4][15] ,const int  nvfaceTet[4][3],const int nvedgeTet[6][2]))[15];
    static PtrConst15int SetonWhatIsFace(int  onWhatIsFace[4][15] ,const int  nvfaceTet[4][3],const int nvedgeTet[6][2]);
    
  template<>
  const int (* const GenericElement<DataTet>::onWhatBorder)[15] = SetonWhatIsFace(onWhatIsFace,nvfaceTet,nvedgeTet) ;
  
  template<> int   GenericMesh<Tet,Triangle3,Vertex3>::kfind=0;
  template<> int   GenericMesh<Tet,Triangle3,Vertex3>::kthrough=0;
  

//  const int (* const SetonWhatIsFace(int  onWhatIsFace[4][15] ,const int  nvfaceTet[4][3],const int nvedgeTet[6][2]))[15]
  PtrConst15int  SetonWhatIsFace(int  onWhatIsFace[4][15] ,const int  nvfaceTet[4][3],const int nvedgeTet[6][2])
  {
    for(int i=0;i<15;++i)
      for(int j=0;j<4;++j)
	onWhatIsFace[j][i]=0;
    for(int j=0;j<4;++j)
      for(int i=0;i<3;++i)
	onWhatIsFace[j][nvfaceTet[j][i]]=1;
    for(int j=0;j<4;++j)
      { 
	onWhatIsFace[j][j+4+6]=3;
	int ff[]={0,0,0,0};
	int jo=1+2+3-nvfaceTet[j][0]-nvfaceTet[j][1]-nvfaceTet[j][2];
	ff[jo]=1;
	for(int i=0;i<6;++i)
	  if(ff[nvedgeTet[i][0]]+ff[nvedgeTet[i][1]]==0)
	    onWhatIsFace[j][i+4]=2;	
      }
    if(0)
      for(int j=0;j<4;++j)
	{
	  for(int i=0;i<15;++i)
	    cout << onWhatIsFace[j][i] << " ";
	  cout << endl;
	}
    
    return onWhatIsFace;
  }
  
  Mesh3::Mesh3(const string  filename)
  {
    int ok=load(filename);		
    if(ok)
      {
	ifstream f(filename.c_str());
	if(!f) {	
	  cerr << "  --  Mesh3::Mesh3 Erreur openning " << filename<<endl;ffassert(0);exit(1);}	
	if(verbosity>1)
	  cout << "  -- Mesh3:  Read On file \"" <<filename<<"\""<<  endl;
	if(filename.rfind(".msh")==filename.length()-4) 
	    readmsh(f);
        else 
	    read(f);
      }
    
    BuildBound();
    if(nt > 0){ 
      BuildAdj();
      Buildbnormalv();  
      BuildjElementConteningVertex();  
    }

    
    if(verbosity>1)
      cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;  
    if(verbosity)
      cout << "  -- Mesh3 : "<<filename  << ", d "<< 3  << ", n Tet " << nt << ", n Vtx "
	   << nv << " n Bord " << nbe << endl;
    
  }

  
  void  Mesh3::read(istream &f)
  { // read the mesh
    int i;
    //	f >> nv >> nt >> neb ;
    string str;
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
	else if (str=="Tetrahedra")
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
		mes += this->t(i).mesure();}
	  }
	else if (str=="Triangles")
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
		mesb += this->be(i).mesure();
		for(int j=0;j<3;++j)
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
  }
  
  int Mesh3::load(const string & filename)
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
      { 
	if(verbosity>5)
	  cerr << " Erreur ouverture file " << (char *) fileb << " " << (char *) filef << endl;
	return   1;
      }
    int nv,nt,neb;
    nv = GmfStatKwd(inm,GmfVertices);
    nt = GmfStatKwd(inm,GmfTetrahedra);
    neb= GmfStatKwd(inm,GmfTriangles);
    this->set(nv,nt,neb);
    if(verbosity>1)
      cout << "  -- Mesh3(load): "<<pfile <<", ver " << ver << ", d "<< dim  
	   << ", nt " << nt << ", nv " << nv << " nbe:  = " << nbe << endl;
    if(dim  != 3) { 
      cerr << "Err dim == " << dim << " !=3 " <<endl;
      return 2; }
    if( nv<=0 && (nt <0 || nbe <=0)  ) {
      cerr << " missing data "<< endl;
      return 3;
    }
    int iv[4],lab;
    float cr[3];
    // read vertices 
    GmfGotoKwd(inm,GmfVertices);
    int mxlab=0;
    int mnlab=0;
    for(int i=0;i<nv;++i)
      {  
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
    
    
    //    /* read mesh triangles */
    if(nbe > 0) {
      if(mnlab==0 && mxlab==0 )
	{
	  int kmv=0;
	  mesb=0;
	  GmfGotoKwd(inm,GmfTriangles);
	  for(int i=0;i<nbe;++i)
	    {  
	      GmfGetLin(inm,GmfTriangles,&iv[0],&iv[1],&iv[2],&lab);
	      for(int j=0;j<3;++j)
		if(!vertices[iv[j]-1].lab)
		  {
		    vertices[iv[j]-1].lab=1;
		    kmv++;
		  }
	      for (int j=0;j<3;++j)  
		iv[j]--;
	      this->be(i).set(this->vertices,iv,lab);
	      mesb += this->be(i).mesure();
	    }
	  
	  if(kmv&& verbosity>1)
	    cout << "    Aucun label Hack (FH)  ??? => 1 sur les triangle frontiere "<<endl;
	}
      else
	{
	  mesb=0;
	  GmfGotoKwd(inm,GmfTriangles);
	  for(int i=0;i<nbe;++i)
	    {  
	      GmfGetLin(inm,GmfTriangles,&iv[0],&iv[1],&iv[2],&lab);
	      for (int j=0;j<3;++j)  
		iv[j]--;
	      this->be(i).set(this->vertices,iv,lab);
	      mesb += this->be(i).mesure();
	    }
	}
    }
    
    if(nt>0)
      {
	/* read mesh tetrahedra */
	GmfGotoKwd(inm,GmfTetrahedra);
	for(int i=0;i<nt;++i)
	  {  
	    GmfGetLin(inm,GmfTetrahedra,&iv[0],&iv[1],&iv[2],&iv[3],&lab);
	    assert( iv[0]>0 && iv[0]<=nv && iv[1]>0 && iv[1]<=nv && iv[2]>0 && iv[2]<=nv && iv[3]>0 && iv[3]<=nv);
	    for (int j=0;j<4;j++) iv[j]--;
	    this->elements[i].set(vertices,iv,lab); 
	    mes += this->elements[i].mesure();	    
	  }
      }
    GmfCloseMesh(inm);    
    return(0); // OK
    
  }
  
const     string Gsbegin="Mesh3::GSave v0",Gsend="end";  
  void Mesh3::GSave(FILE * ff) const
  {  
    PlotStream f(ff);
    
    f <<  Gsbegin ;
    f << nv << nt << nbe;
    for (int k=0; k<nv; k++) {
      const  Vertex & P = this->vertices[k];		
      f << P.x <<P.y << P.z << P.lab ;
    }
    
    if(nt != 0){
      
      for (int k=0; k<nt; k++) {
	const Element & K(this->elements[k]);
	int i0=this->operator()(K[0]);
	int i1=this->operator()(K[1]);
	int i2=this->operator()(K[2]);
	int i3=this->operator()(K[3]);
	int lab=K.lab;
	f << i0 << i1 << i2 << i3 << lab;
      }
    }
    
    
    for (int k=0; k<nbe; k++) {
      const BorderElement & K(this->borderelements[k]);
      int i0=this->operator()(K[0]);
      int i1=this->operator()(K[1]);
      int i2=this->operator()(K[2]);
      int lab=K.lab;
      f << i0 << i1 << i2  << lab;
    }
    f << Gsend;
  }
	      
   Mesh3::Mesh3(const  Serialize &serialized)
   :GenericMesh<Tet,Triangle3,Vertex3> (serialized) 
    {
	BuildBound();
	if(verbosity>1)
	    cout << "  -- End of serialized: mesure = " << mes << " border mesure " << mesb << endl;  
	
	if(nt > 0){ 
	    BuildAdj();
	    Buildbnormalv();  
	    BuildjElementConteningVertex();  
	}
	
	if(verbosity)
	    cout << "  -- Mesh3  (serialized), d "<< 3  << ", n Tet " << nt << ", n Vtx "
	    << nv << " n Bord " << nbe << endl;
    }
  Mesh3::Mesh3(FILE *f)
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
    
    if(verbosity)
      cout << "  -- Mesh3  (File *), d "<< 3  << ", n Tet " << nt << ", n Vtx "
	   << nv << " n Bord " << nbe << endl;
    
  }
  
  void Mesh3::GRead(FILE * ff)
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
      f >> P.x >>P.y >> P.z >> P.lab ;
    }
    mes=0.;
    mesb=0.;
    
    if(nt != 0)
      {
	
	for (int k=0; k<nt; k++) {
	  int i[4],lab;
	  Element & K(this->elements[k]);
	  f >> i[0] >> i[1] >> i[2] >> i[3] >> lab;
	  K.set(this->vertices,i,lab);
	  mes += K.mesure();	    
	  
	}
      }
    
    
    for (int k=0; k<nbe; k++) {
      int i[4],lab;
      BorderElement & K(this->borderelements[k]);
      f >> i[0] >> i[1] >> i[2]  >> lab;
      K.set(this->vertices,i,lab);
      mesb += K.mesure();	    
      
    }
    f >> s;
    ffassert( s== Gsend);
  }
    void Mesh3::readmsh(ifstream & f)
    {  

	f >> nv >> nt >> nbe;
	if(verbosity>1)
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
	
	if(nt != 0)
	  {
	      
	      for (int k=0; k<nt; k++) {
		  int i[4],lab;
		  Element & K(this->elements[k]);
		  f >> i[0] >> i[1] >> i[2] >> i[3] >> lab;
		  K.set(this->vertices,i,lab);
		  mes += K.mesure();	    
		  
	      }
	  }
	
	
	for (int k=0; k<nbe; k++) {
	    int i[4],lab;
	    BorderElement & K(this->borderelements[k]);
	    f >> i[0] >> i[1] >> i[2]  >> lab;
	    K.set(this->vertices,i,lab);
	    mesb += K.mesure();	    
	    
	}

    }
    
  
  
  int Mesh3::Save(const string & filename) const
  {
    int ver = GmfFloat, outm;
    if ( !(outm = GmfOpenMesh(filename.c_str(),GmfWrite,ver,3)) ) {
      cerr <<"  -- Mesh3::Save  UNABLE TO OPEN  :"<< filename << endl;
      return(1);
    }
    
    float fx,fy,fz;
    GmfSetKwd(outm,GmfVertices,this->nv);
    for (int k=0; k<nv; k++) {
      const  Vertex & P = this->vertices[k];
      GmfSetLin(outm,GmfVertices,fx=P.x,fy=P.y,fz=P.z,P.lab);
    }
    
    if(nt != 0){
      GmfSetKwd(outm,GmfTetrahedra,nt);
      for (int k=0; k<nt; k++) {
	const Element & K(this->elements[k]);
	int i0=this->operator()(K[0])+1;
	int i1=this->operator()(K[1])+1;
	int i2=this->operator()(K[2])+1;
	int i3=this->operator()(K[3])+1;
	int lab=K.lab;
	GmfSetLin(outm,GmfTetrahedra,i0,i1,i2,i3,lab);
      }
    }
    
    GmfSetKwd(outm,GmfTriangles,nbe);
    for (int k=0; k<nbe; k++) {
      const BorderElement & K(this->borderelements[k]);
      int i0=this->operator()(K[0])+1;
      int i1=this->operator()(K[1])+1;
      int i2=this->operator()(K[2])+1;
      int lab=K.lab;
      GmfSetLin(outm,GmfTriangles,i0,i1,i2,lab);
    }
    
    GmfCloseMesh(outm);
    return (0);
    
  }

   int Mesh3::SaveSurface(const string & filename) const
  {
    int ver = GmfFloat, outm;
    if ( !(outm = GmfOpenMesh(filename.c_str(),GmfWrite,ver,3)) ) {
      cerr <<"  -- Mesh3::Save  UNABLE TO OPEN  :"<< filename << endl;
      return(1);
    }

    // Number of Vertex in the surface
    int *v_num_surf=new int[nv];
    int *liste_v_num_surf=new int[nv];
    for (int k=0; k<nv; k++){ 
      v_num_surf[k]=-1;
      liste_v_num_surf[k]=0;
    }
    // Search Vertex on the surface
    int nbv_surf=0;
    for (int k=0; k<nbe; k++) {
      const BorderElement & K(this->borderelements[k]);     
      for(int jj=0; jj<3; jj++){
	int i0=this->operator()(K[jj]);
	if( v_num_surf[i0] == -1 ){
	  v_num_surf[i0] = nbv_surf;
	  liste_v_num_surf[nbv_surf]= i0;
	  nbv_surf++;
	}
      }
    }

    float fx,fy,fz;
    GmfSetKwd(outm,GmfVertices,nbv_surf);
    for (int k=0; k<nbv_surf; k++) {
      int k0 = liste_v_num_surf[k];
      const  Vertex & P = this->vertices[k0];
      GmfSetLin(outm,GmfVertices,fx=P.x,fy=P.y,fz=P.z,P.lab);
    }
    
    GmfSetKwd(outm,GmfTriangles,nbe);
    for (int k=0; k<nbe; k++) {
      const BorderElement & K(this->borderelements[k]);
      int i0=v_num_surf[this->operator()(K[0])]+1;
      int i1=v_num_surf[this->operator()(K[1])]+1;
      int i2=v_num_surf[this->operator()(K[2])]+1;
      int lab=K.lab;

      assert( i0-1 < nbv_surf &&  i1-1 < nbv_surf &&  i2-1 < nbv_surf );
      assert( 0<i0 && 0<i1 && 0<i2 );

      GmfSetLin(outm,GmfTriangles,i0,i1,i2,lab);
    }
    
    GmfCloseMesh(outm);

    delete [ ] v_num_surf;
    delete [ ] liste_v_num_surf;

    return (0); 
  }


  int Mesh3::SaveSurface(const string & filename1,const string & filename2) const
  {    
    // Number of Vertex in the surface
    int *v_num_surf=new int[nv];
    int *liste_v_num_surf=new int[nv];
    for (int k=0; k<nv; k++){ 
      v_num_surf[k]=-1;
      liste_v_num_surf[k]=0;
    }
    // Search Vertex on the surface
    int nbv_surf=0;
    for (int k=0; k<nbe; k++) {
      const BorderElement & K(this->borderelements[k]);     
      for(int jj=0; jj<3; jj++){
	int i0=this->operator()(K[jj]);
	if( v_num_surf[i0] == -1){
	  v_num_surf[i0] = nbv_surf;
	  nbv_surf++;
	}
      }
    }

    // file .points
    FILE *fpoints = fopen(filename1.c_str(),"w");
    fprintf(fpoints,"%i\n",nbv_surf);
    
    for (int k=0; k<nbv_surf; k++) {
      int k0 = liste_v_num_surf[k];
      const  Vertex & P = this->vertices[k];
      fprintf(fpoints,"%f %f %f %i\n",P.x,P.y,P.z,P.lab);
    }
    fclose(fpoints);
    
    // file .faces
    FILE *ffaces = fopen(filename2.c_str(),"w");
    fprintf(ffaces,"%i\n",nbe);
    for (int k=0; k<nbe; k++) {
      const BorderElement & K(this->borderelements[k]);
      int i0=this->operator()(K[0]);
      int i1=this->operator()(K[1]);
      int i2=this->operator()(K[2]);
      int lab=K.lab;
      int label0= this->vertices[i0].lab; 
      int label1= this->vertices[i1].lab; 
      int label2= this->vertices[i2].lab;
      //GmfSetLin(outm,GmfTriangles,i0,i1,i2,lab);
      int nature=3;
      int i0v=v_num_surf[i0]+1;
      int i1v=v_num_surf[i1]+1;
      int i2v=v_num_surf[i2]+1;
      assert( i0v-1 < nbv_surf &&  i1v-1 < nbv_surf &&  i2v-1 < nbv_surf );
      assert( 0<i0v && 0<i1v && 0<i2v );

      fprintf(ffaces,"%i %i %i %i %i %i %i %i\n", nature, i0v, i1v, i2v, lab, label0, label1, label2);
    }
    fclose(ffaces);
    
    delete [ ] v_num_surf;
    delete [ ] liste_v_num_surf;
    return (0);  
  }


  Mesh3::Mesh3(int nnv, int nnt, int nnbe, Vertex3 *vv, Tet *tt, Triangle3 *bb)
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
    
//  Add FH to be consitant we all constructor ...  July 09
      BuildBound();
      if(nt > 0){ 
	  BuildAdj();
	  Buildbnormalv();  
	  BuildjElementConteningVertex();  
      }
//  end add       
          
    if(verbosity>1)
      cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;  
    
  }
  
  Mesh3::Mesh3(int nnv, int nnbe, Vertex3 *vv, Triangle3 *bb)
  {
  
    nv = nnv;
    nbe =nnbe;
    
    vertices = vv;
    borderelements = bb;
    
    mes=0.;
    mesb=0.;
    
    for (int i=0;i<nbe;i++)  
      mesb += this->be(i).mesure();  

//  Add FH to be consitant we all constructor ...  July 09
      BuildBound();
      if(nt > 0){ 
	  BuildAdj();
	  Buildbnormalv();  
	  BuildjElementConteningVertex();  
      }
//  end add       
      
    if(verbosity>1)
      cout << "  -- End of Construct  mesh3: mesure = " << mes << " border mesure " << mesb << endl;  
  }

  void Mesh3::flipSurfaceMesh3(int surface_orientation)
  {
    /* inverse the orientation of the surface if necessary*/
    /* and control that all surfaces are oriented in the same way*/
    int nbflip=0;
    for (int i=0;i<this->nbe;i++)
      { 
	double mes_triangle3= this->be(i).mesure();
	
	if( surface_orientation*mes_triangle3 < 0){
	  const Triangle3 &K( this->be(i) );
	  int iv[3];       
	  
	  iv[0] = this->operator()(K[0]);
	  iv[1] = this->operator()(K[1]);
	  iv[2] = this->operator()(K[2]);
	  
	  int iv_temp=iv[1];
	  iv[1]=iv[2];
	  iv[2]=iv_temp;
	  this->be(i).set( this->vertices, iv, K.lab ) ;
	  nbflip++;
	}
      }
    assert(nbflip==0 || nbflip== this->nbe); 
  }


  int  signe_permutation(int i0,int i1,int i2,int i3)
  {
    int p=1;
    if(i0>i1) Exchange(i0,i1), p = -p;
    if(i0>i2) Exchange(i0,i2), p = -p;
    if(i0>i3) Exchange(i0,i3), p = -p;
    if(i1>i2) Exchange(i1,i2), p = -p;
    if(i1>i3) Exchange(i1,i3), p = -p;
    if(i2>i3) Exchange(i2,i3), p = -p;
    return p;
  }


  int  WalkInTet(const Mesh3 & Th,int it, R3 & Phat,const R3 & U, R & dt)
  {
    bool ddd=verbosity>200;
      bool nomove=true;
    R lambda[4];
    Phat.toBary(lambda);
    if(ddd) cout << "\n\n\n   WT: "  << Phat << " :  "  << lambda[0] << " " <<lambda[1] <<" " <<lambda[2] << " " <<lambda[3] << " u = "<< U << " dt " << dt <<endl;
#ifndef NDEBUG      
      for(int i=0;i<4;++i)
      assert(lambda[i]<1.000001 && lambda[i]>-0.0000001);
#endif
    typedef R3 Rd;
    const Mesh3::Element & T(Th[it]);
    const int nve = 4;
    const Rd  Q[nve]={(const R3) T[0],(const R3) T[1],(const R3) T[2],(const R3) T[3]};
    
    Rd P  =T(Phat);
    
    //  cout << " " << u << " " << v ;
    Rd PF = P + U*dt;
    
    //  couleur(15);MoveTo( P); LineTo( PF);
    R l[nve];
    double Det=T.mesure()*6;
    l[0] = det(PF  ,Q[1],Q[2],Q[3]);
    l[1] = det(Q[0],PF  ,Q[2],Q[3]); 
    l[2] = det(Q[0],Q[1],PF  ,Q[3]); 
    l[3] = Det - l[0]-l[1]-l[2];
    l[0] /= Det;
    l[1] /= Det;
    l[2] /= Det;
    l[3] /= Det;
     if(ddd)  cout << "\t\t\tWT " << it << ", " << Phat << ",  PF=" << PF
		   << " :  "  << l[0] << " " <<l[1] <<" " <<l[2] << " " <<l[3] 
	           << " == " << det(Q[0],Q[1],Q[2],PF)/Det
		   << " : l (in) "  << lambda[0] << " " <<lambda[1] <<" " <<lambda[2] << " " <<lambda[3] 
		   << " PF= K(l) = " << Th[it](R3(l+1)) 
		   <<endl ;
		  
    const R eps = 1e-8;
    int neg[nve],k=0;
    int kk=-1;
    if (l[0]>-eps && l[1]>-eps && l[2]>-eps && l[3]>-eps) 
      {
	dt =0;
	Phat=R3(l+1);
	nomove=false;
	return -1;
      }
    else 
      {
	// on regarde de les reelement negatif 
        // on ne veut pas des points sur les faces.
        // car sinon il va y avoir un probleme ans on va projete sur la face
	//  et remettre le point dans le tetraedre.
	if (l[0]<=-eps ) neg[k++]=0;
	if (l[1]<=-eps ) neg[k++]=1;
	if (l[2]<=-eps ) neg[k++]=2;
	if (l[3]<=-eps ) neg[k++]=3;
	
	R eps1 = T.mesure()   * 1.e-5;
	   if(ddd)  cout << " k= " << k << endl;
    
	if (k==3) //  3 face de sortie possible 
	  {
	    // let j be the vertex beetween the 3 faces 
	    int j = 6-neg[0]-neg[1]-neg[2]; //  sommet intersection des 3 faces.
	    int i0 = Tet::nvface[j][0];
	    int i1 = Tet::nvface[j][1];
	    int i2 = Tet::nvface[j][2];
	     if(ddd)  cout << "  -------- " << j << " " << i0 << " " << i1 << " " << i2  << endl;
	    //  le tet i0,i1,i2,j est positif. 
	    assert(signe_permutation(i0,i1,i2,j)==1);
	    // 
	    R v0= det(Q[i0],Q[j],P,PF); 
	    R v1= det(Q[i1],Q[j],P,PF); 
	    R v2= det(Q[i2],Q[j],P,PF); 
	     if(ddd)   cout << "\t\t\t " << j << " v0123 =" << v0 << " "<< v1 << " " << v2 << endl;
	    if( v0 > eps && v1 < -eps ) 
	      kk= i1 ;// on sort par la face j i0, j1
	    else if( v1 > eps && v2 < -eps ) 
	      kk= i0 ;
	    else if( v2 > eps && v0 < -eps ) 
	      kk= i1 ;
	    else 
	      {  // on ne sort pas par une face 
		int nul[3], nn=0, mm=3;
		  if (Abs(v0) <=eps) nul[nn++]=i0; else nul[--mm]=i0;
		  if (Abs(v1) <=eps) nul[nn++]=i1; else nul[--mm]=i1;
		  if (Abs(v2) <=eps) nul[nn++]=i2; else nul[--mm]=i2;
		assert(nn>0); 
		if(nn == 1) // on sort par l'arete nul[0] entre le face   nul[1] et nul[2]
		  kk =  nul[1+(rand()/(RAND_MAX/2))%2];		  
		else // on sort par le sommet j.  on choisi la face alleatoirement 
		  kk = nul[(rand()/(RAND_MAX/3))%3];
		  
	      }
	  }
	else if (k==2)
	  {
	    //  numero des l'arete entre les 2 faces
	    int i0=neg[0],i1=neg[1];
	    int e = i0 + i1 - (i0==0); 
	    // on a:
	    //   e      =     0        1       2      3        4      5
	    //  (i0,i1) =   (0,1)  , (0,2), (0,3) , (1,2)  , (1,3), (2,3) 
	    // avec i0,i1 sont les sommets qui ne sont pas dans l'arete
	    int   jj0[6] = {2,3,1,0,2,0};
	    int   jj1[6] = {3,1,2,3,0,1};
	    int j0 = jj0[e];
	    int j1 = jj1[e];
	     if(ddd)   cout << " e " << e << " i0 " << i0 << " " << i1 << " j0 =" << j0 << " " << j1 << endl;
	    // le tet  j0,j1,i0,i1  doit est positif (ie. la pemutation est positive)
	    // de meme  i0,i1,j0,j1
	    assert(signe_permutation(j0,j1,i0,i1)==1);
	    R v0= det(Q[j0],Q[j1],P,PF); 
            if(ddd) cout << " v0 =" << v0 <<endl;
	    if( Abs(v0) < eps  ) 
	      {
	      // on sort par l'arete  j0,j1
	      // on choisi aleatoirement la face de sortie 
	      kk = (rand()/(RAND_MAX/2)) ? i0 : i1; 
	      if(ddd)
		  cout << " rand choose  2 :  " << kk << endl;
	      }
	    else 
	      kk= v0 >0 ? i0 : i1; // Attention dyslexie ici durdur FH....
	    
	  }
	else if (k==1) //  une face possible de sortie (cas simple)
	  kk = neg[0];
	
	if(kk>=0)
	  {
	    R d=lambda[kk]-l[kk];
	    if ( l[kk] )
	     {
	    throwassert(d);
	    R coef =  lambda[kk]/d;
	    R coef1 = 1-coef;
	    nomove= (coef<1.e-6);
	    dt        = dt*coef1;
	    lambda[0] = lambda[0]*coef1 + coef *l[0];
	    lambda[1] = lambda[1]*coef1 + coef *l[1];
	    lambda[2] = lambda[2]*coef1 + coef *l[2];
	    lambda[3] = lambda[3]*coef1 + coef *l[3];
            if(ddd) cout << "   \t\t -> kk=" << kk << " d=" << d << " , l= "<< lambda[0]  << " " 
			 <<lambda[1] << " " <<lambda[2] << " " << lambda[3] << endl;
	    lambda[kk] =0;
	     }
	      
	 
	  
	  }
	 
      }
    if(nomove )
	// on ne bouge pas on utilse Find ... 
	  {
	      R dx2= (U,U)*dt*dt;
	      R ddt=dt, dc=1;
	      // if Udt < h/2 => recherche un point final  
	      if(dx2*dx2*dx2 > Det*Det/4)
		  dt=0;       
	      else 
		{ 
		    dc=0.25;
		    ddt=dt*dc;
		    PF= P + U*ddt; // on avance que d'un 1/4
		    dt -= ddt;
		}
	      bool outside;
	      const Mesh3::Element  *K=Th.Find(PF, Phat,outside,&Th[it]);
	      if(outside) dt=0; // on a fini 
	      if(ddd) cout << "   \t ***** WT :  Lock -> Find P+U*ddt*c "<< it<< " " << " -> "<< Th(K) 
		  << " dt = " << dt << " c = " << dc << " outside: "<< outside <<" , PF " << PF << endl;
	      return 4+Th(K);
     }
      
    //  on remet le point dans le tet. 
    int jj=0;
    R lmx=lambda[0];
    if (lmx<lambda[1])  jj=1, lmx=lambda[1];
    if (lmx<lambda[2])  jj=2, lmx=lambda[2];
    if (lmx<lambda[3])  jj=3, lmx=lambda[3];
    if(lambda[0]<0) lambda[jj] += lambda[0],lambda[0]=0;
    if(lambda[1]<0) lambda[jj] += lambda[1],lambda[1]=0;
    if(lambda[2]<0) lambda[jj] += lambda[2],lambda[2]=0;
    if(lambda[3]<0) lambda[jj] += lambda[3],lambda[3]=0;
    Phat=R3(lambda+1);
    if(ddd) cout  << "\t\t\t -> "<< dt << " : "  << Phat << " K(Phat) ="<< Th[it](Phat) <<  ", " << kk << " jj= "<< jj << " "<< lmx << endl; 
    assert(kk<0 || lambda[kk]==0);
    return kk;
  }        
	      
    
} //   End  namespace Fem2D
  
	      
