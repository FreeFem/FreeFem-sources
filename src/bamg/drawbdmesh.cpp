// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33) 1 39 63 55 14       
//
// AUTHOR:   F. Hecht,    
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <setjmp.h>
#include <new.h>
#include <assert.h>
#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"
//#include <fstream.h>
int initgraph=0;
void NewHandler();
void NewHandler()
{ cerr << " Not enought memory " << endl;
  exit(1);
}

void forDebug();

void forDebug()
{ 
#ifdef DRAWING 
if(initgraph) {
  rattente(1);
  closegraphique();}
#endif
//  DelAllocData();
  cout << "on a fini" << endl;
}
void  MeshErrorIO(ios& )
{
  MeshError(999);
  exit(1);
}

int main(int argc, char **argv)
{ 
  MeshIstreamErrorHandler = MeshErrorIO;
  //  atexit( forDebug);
  double raison=0;
  verbosity = 9;
    set_new_handler( &NewHandler);
   char * fmeshback= "1.mesh";
   char * fmetrix= 0; 
    ostream *f_metr(&cout);
    if (argc >=2) 
     fmeshback = argv[1];
     if (argc >=3) 
     fmetrix   = argv[2];
     if (argc >=4) 
      raison   =  atof(argv[3]);
   
#ifdef DRAWING 
   initgraphique(); 
   initgraph=1;

#endif
   cout << "open " << fmeshback << " raison = " << raison <<  endl;
   Triangles TTh(fmeshback);
   Real8 hmin = TTh.MinimalHmin();
   Real8 hmax = TTh.MaximalHmax();
   Real8 coef = 1;
   if (fmetrix)
         TTh.ReadMetric(fmetrix,hmin,hmax,coef);
        
#ifdef DRAWING 

     reffecran();
     TTh.InitDraw();
     TTh.inquire();
	if(raison&& fmetrix)
	  {
	    TTh.inquire();
	    TTh.SmoothMetric(raison);
	  }
     if (fmetrix)
       for(Int4 i=0;i<TTh.nbv;i++)
	   ((Metric) TTh.vertices[i]).Draw(TTh.vertices[i]);    
     TTh.Draw();
     TTh.inquire();
     closegraphique();
#endif
   return(0); 
}
