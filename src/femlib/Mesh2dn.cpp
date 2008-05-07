#include <fstream>
#include <iostream>
#include "ufunction.hpp"
#include "error.hpp"
#include "RNM.hpp"
#include "libmesh5.h"


#include "Mesh2dn.hpp"

 namespace Fem2D {
long verbosity=1;

const R1 R1::KHat[2]={R1(0),R1(1)};
const R2 R2::KHat[3]={R2(0,0),R2(1,0),R2(0,1)};
const  R3 R3::KHat[4]={R3(0,0,0),R3(1,0,0),R3(0,1,0),R3(0,0,1)};

static const int  nvfaceTet[4][3]  = { {2,1,3},{0,2,3},{1,0,3},{0,1,2} };
static const int  nvedgeTet[6][2] = { {0,1},{0,2},{0,3},{0,1},{1,2},{2,3} };

static const int  nvfaceTria[1][3]  = { {0,1,2} };
static const int  nvedgeTria[3][2] = { {1,2},{2,0},{0,1}}; 

static const int   nvfaceSeg[1][3]  = {{-1,-1,1}};
static const int  nvedgeSeg[1][2] = { {0,1} };


template<>
const int (* const GenericElement<DataTriangle2>::nvface)[3] = nvfaceTria ;
template<>
const int (* const GenericElement<DataTriangle2>::nvedge)[2] = nvedgeTria ;
template<>
const int (* const GenericElement<DataTriangle2>::nvadj)[2] = nvedgeTria ;
template<> const int  GenericElement<DataTriangle2>::nitemdim[4] = {3,3,1,0 }  ;


static const int onWhatIsEdge[3][7] = {  {0,1,3, 2,0,0, 0}, // edge 0 
    {3,0,1, 0,2,0, 0},
    {1,3,0, 0,0,2, 0}};

template<>
const int (* const GenericElement<DataTriangle2>::onWhatBorder)[7] = onWhatIsEdge ;


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
	//	vertices[i].z=cr[2];
    }
      else
	GmfGetLin(inm,GmfVertices,&vertices[i].x,&vertices[i].y,&lab);	    
      vertices[i].lab=lab;
      
      mxlab= max(mxlab,lab);
      mnlab= min(mnlab,lab);
    }
  
  
  //    /* read mesh triangles */
  if(mnlab==0 &&mxlab==0 )
    {
      mes=0;
      GmfGotoKwd(inm,GmfTriangles);
      for(int i=0;i<nbe;++i)
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
  return(0); // OK
  
}

int Mesh2::Save(const string & filename)
{
  int ver = GmfFloat, outm;
  if ( !(outm = GmfOpenMesh(filename.c_str(),GmfWrite,ver,2)) ) {
    cerr <<"  -- Mesh2::Save  UNABLE TO OPEN  :"<< filename << endl;
    return(1);
  }
  float fx,fy;
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
}
