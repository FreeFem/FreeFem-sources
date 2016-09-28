// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// AUTHOR:   F. Hecht,    
// ORG    :  UPMC
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
/*
 
 This file is part of Freefem++ or Bamg
 
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
 */
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97
// ORIG-DATE:     Dec 97
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <setjmp.h>
#include <new>
#include <cassert>
#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"

using namespace std;
using namespace bamg;
//long verbosity=2;
#ifdef DRAWING 
bool withrgraphique=  initgraphique;
#else
bool withrgraphique=false; 
#endif
#ifdef HAVE_UNISTD_H 
#include <unistd.h>
#else
#ifdef __MWERKS__
#define R_OK 0
#define F_OK 0
#define W_OK 0
#include <unix.h>
inline int access( char *fileName, int notUsed )
{
	struct stat		statRec;
	return  stat( fileName, &statRec );
}
#endif

#ifdef PURE_WIN32
#define R_OK 0
#define F_OK 0
#define W_OK 0
#define access(i,j) 1
#endif
#endif
int initgraph=0;
void NewHandler();
void NewHandler()
{ cerr << " Not enought memory " << endl;
  exit(1);
}


void  MeshErrorIO(ios& )
{
  MeshError(999);
  exit(1);
}

int main(int argc, char **argv)
{ 
  Real8 cutoffradian = 30.0*Pi/180.0;
  char * fin=0,*fout=0,*fgeom=0;
  verbosity = 2;
  int ok=0,i=0;
  if (argc >= 3)
    {
      fin=argv[1];
      fout=argv[2];
      ok= !access(fin,R_OK) && (access(fout,F_OK) || !access(fout,W_OK)) ;
    }
  
  for (i=3;i<argc && ok ;i+=2)
    if (!strcmp(argv[i],"-g")) 
      {
	fgeom=argv[i+1];
	ok= ok &&  (access(fgeom,F_OK) || !access(fgeom,W_OK) );
      }
    else if (!strcmp(argv[i],"-thetamax"))
      cutoffradian = atof(argv[i+1])*Pi/180.0;
    else if (!strcmp(argv[i],"-v"))
      verbosity = atoi(argv[i+1]);
    else
      ok=0;
      
  
  
  if (!ok)
    {
      cout << "2D mesh convertisor of type   .am_fmt .amdba .am .nopo .msh .ftq + bd mesh\n";
      cout << "                    in type   .am_fmt .amdba .am .nopo .msh .ftq + bd mesh\n";
      cout << endl;
      cout << " usage: cnmsh2 InMesh OutMesh \n";
      cout << " or   : cnmsh2 Inmeshfile  OutBdmeshfile  [-g OutGeomfile ] [ -thetamax theta ] [-v level]\n";
      cout << endl;
      cout << " where  theta is the angular limit for  smooth curve in  degre \n";
      cout << "        level is level of verbosity in [0..99] 0 => no message, 99 to much \n";
      cout << endl;
      cout << " remark: ";
      cout << "   the  type of OutMesh is defined with the suffix\n";
      cout << "   si the suffixe is not in .am_fmt .amdba .am .nopo .msh .ftq, then \n";
      cout << "   the mesh is a bd mesh (cf .mesh)" << endl;
      return 1;
    }
  MeshIstreamErrorHandler = MeshErrorIO;
  set_new_handler( &NewHandler);
#ifdef DRAWING
  if(initgraphique)
    { 
      initgraphique(); 
      initgraph=1;
    }
#endif
  //    cout << "open file in " << fin<< endl;
  Triangles Th(fin,cutoffradian);
  if (!fgeom)
    ((void )  0); //    Th.Write( (const char *) fout);
  else
    {
      //      cout << " geom " << fgeom << endl;
      Th.Gh.Write(fgeom);
      Th.Write(fout,Triangles::BDMesh);
    }

   return 0;
}
