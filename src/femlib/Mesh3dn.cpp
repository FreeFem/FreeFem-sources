#include <fstream>
#include <iostream>
#include <cstring>
#include "libmesh5.h"
#include "ufunction.hpp"
#include "error.hpp"
#include "RNM.hpp"
#include "Mesh3dn.hpp"


namespace Fem2D {
  
  /*
    const short int v_tet_face[4][3]=  {{3,2,1},{0,2,3},{ 3,1,0},{ 0,1,2}};
    const short int a_tet_face[4][3]=  {{ 0,1,0},{ 0,2,0},{ 0,3,1},{ 1,2,1}};
    const bool  sig_tet_face[4][3]=  {{ 0,1,0},{ 1,0,1},{ 0,1,0},{ 1,0,1}};
    const short int v_tet_edge[6][2]= {{ 1,2},{1,3},{1,4},{2,3},{2,4},{3,4}}; 
    const short int fadj_tet_edge[6][2]= {{4,3},{2,4},{3,2},{4,1},{1,3},{2,1}};
    const short int op_tet_edge[6]={ 6, 5, 4, 3, 2, 1};  
  */
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
  static const int onWhatIsEdge[3][7] = {
    {0,1,3, 2,0,0, 0}, // edge 0 
    {3,0,1, 0,2,0, 0},
    {1,3,0, 0,0,2, 0} };
  
  template<>
  const int (* const GenericElement<DataTriangle3>::onWhatBorder)[7] = onWhatIsEdge ;
  
  
  template<>
  const int (* const GenericElement<DataTet>::nvface)[3] = nvfaceTet ;
  template<>
  const int (* const GenericElement<DataTet>::nvedge)[2] = nvedgeTet ;
  template<>
  const int (* const GenericElement<DataTet>::nvadj)[3] = nvfaceTet ;
  template<> const int  GenericElement<DataTet>::nitemdim[4] = {4,6,4,1 }  ;
  
  int onWhatIsFace[4][15] ; 
  
  static const int (* const SetonWhatIsFace(int  onWhatIsFace[4][15] ,const int  nvfaceTet[4][3],const int nvedgeTet[6][2]))[15];
  
  template<>
  const int (* const GenericElement<DataTet>::onWhatBorder)[15] = SetonWhatIsFace(onWhatIsFace,nvfaceTet,nvedgeTet) ;
  
  template<> int   GenericMesh<Tet,Triangle3,Vertex3>::kfind=0;
  template<> int   GenericMesh<Tet,Triangle3,Vertex3>::kthrough=0;

  
  const int (* const SetonWhatIsFace(int  onWhatIsFace[4][15] ,const int  nvfaceTet[4][3],const int nvedgeTet[6][2]))[15]
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
{ // read the mesh
  int i;
  int ok=load(filename);
  if(ok)
    {
      ifstream f(filename.c_str());
      if(!f) {	
	cerr << "  --  Mesh3::Mesh3 Erreur openning " << filename<<endl;ffassert(0);exit(1);}
      if(verbosity>1)
      cout << "  -- Mesh3:  Read On file \"" <<filename<<"\""<<  endl;
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
      assert( nt >0 && nv>0) ;
    }
  BuildBound();
  BuildAdj();
  Buildbnormalv();  
  BuildjElementConteningVertex();  
    
  if(verbosity>1)
  cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;  
  if(verbosity)
    cout << "  -- Mesh3 : "<<filename  << ", d "<< 3  << ", n Tet " << nt << ", n Vtx "
	 << nv << " n Bord " << nbe << endl;
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
    { if(verbosity>5)
	cerr << " Erreur ouverture file " << (char *) fileb << " " << (char *) filef << endl;
      return   1;
    }
  int nv,nt,neb;
  nv = GmfStatKwd(inm,GmfVertices);
  nt = GmfStatKwd(inm,GmfTetrahedra);
  neb=GmfStatKwd(inm,GmfTriangles);
  this->set(nv,nt,neb);
  cout << "  -- Mesh3(load): "<<pfile <<", ver " << ver << ", d "<< dim  << ", nt " << nt << ", nv " << nv << " nbe:  = " << nbe << endl;
  if(dim  != 3) { 
    cerr << "Err dim == " << dim << " !=3 " <<endl;
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
	vertices[i].z=cr[2];}
      else
	GmfGetLin(inm,GmfVertices,&vertices[i].x,&vertices[i].y,&vertices[i].z,&lab);	    
      vertices[i].lab=lab;
      
      mxlab= max(mxlab,lab);
      mnlab= min(mnlab,lab);
    }
  
  
  //    /* read mesh triangles */
  if(mnlab==0 &&mxlab==0 )
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
  
  GmfCloseMesh(inm);    
  return(0); // OK
  
}


int Mesh3::Save(const string & filename)
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
	
	
  BuildBound();
  BuildAdj();
  Buildbnormalv();  
  BuildjElementConteningVertex();  
    
  if(verbosity>1)
  cout << "  -- End of read: mesure = " << mes << " border mesure " << mesb << endl;  
	
} 
}