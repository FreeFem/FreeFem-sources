#include <fstream>
#include <iostream>
#include "ufunction.hpp"
#include "error.hpp"
#include "RNM.hpp"
#include "libmesh5.h"

#include "Mesh1dn.hpp"

   //long verbosity=1;

namespace Fem2D {

const R1 R1::KHat[2]={R1(0),R1(1)};
const R2 R2::KHat[3]={R2(0,0),R2(1,0),R2(0,1)};
const  R3 R3::KHat[4]={R3(0,0,0),R3(1,0,0),R3(0,1,0),R3(0,0,1)};

static const int  nvfaceTet[4][3]  = { {2,1,3},{0,2,3},{1,0,3},{0,1,2} };
static const int  nvedgeTet[6][2] = { {0,1},{0,2},{0,3},{0,1},{1,2},{2,3} };

static const int  nvfaceTria[1][3]  = { {0,1,2} };
static const int  nvedgeTria[3][2] = { {1,2},{2,0},{0,1}}; 

static const int   nvfaceSeg[1][3]  = {{-1,-1,-1}};
static const int  nvedgeSeg[1][2] = { {0,1} };
static const int  nvadjSeg[2][1] = { {0},{1} };


template<>
const int (* const GenericElement<DataSeg1>::nvface)[3] = 0 ;
template<>
const int (* const GenericElement<DataSeg1>::nvedge)[2] = nvedgeSeg; //nvedgeTria ;
template<>
const int (* const GenericElement<DataSeg1>::nvadj)[1] = nvadjSeg ;
template<> const int  GenericElement<DataSeg1>::nitemdim[4] = {2,1,0,0 }  ;


static const int onWhatIsVertex[2][3] = {  {1,0,0}, // sommet  0 
					   {0,1,0} // sommet 1 
				      };

template<>
const int (* const GenericElement<DataSeg1>::onWhatBorder)[3] = onWhatIsVertex ;


Mesh1::Mesh1(const char * filename)
{ // read the mesh  

  int nt,nv,nbe;
  int ok=0;//load(filename);
  if(ok)
    {
      ifstream f(filename);
      if(!f) {cerr << "Mesh1::Mesh1 Erreur openning " << filename<<endl;exit(1);}
      if(verbosity)
      cout << " Read On file \"" <<filename<<"\""<<  endl;
      f >> nv >> nt >> nbe ;
      this->set(nv,nt,nbe);
      if(verbosity)
      cout << "  -- Nb of Vertex " << nv << " " << " Nb of Seg " << nt 
	   << " , Nb of border Vertex " << nbe <<  endl;
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
}

 
