// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0
//
// AUTHOR:   F. Hecht,
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr
//
// ORIG-DATE:     Dec 97
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
 */
// E-MAIL :   Frederic.Hecht@Inria.fr
//
// ORIG-DATE:     Dec 97#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <setjmp.h>
#include <new>
#include <cassert>
#include "Meshio.h"
#include "Mesh2.h"
#include "QuadTree.h"
using namespace bamg;
extern long verbosity;

#ifdef DRAWING
bool withrgraphique = true;
#else
bool withrgraphique = false;
#endif
//#include <fstream.h>
int initgraph = 0;
void NewHandler( );
void NewHandler( ) {
  cerr << " Not enought memory " << endl;
  exit(1);
}

void forDebug( );

void forDebug( ) {
#ifdef DRAWING
  if (initgraph) {
    rattente(1);
    closegraphique( );
  }
#endif
  //  DelAllocData();
  cout << "on a fini" << endl;
}
void MeshErrorIO(ios &) {
  MeshError(999);
  exit(1);
}

int mymain(int argc, char **argv) {
  MeshIstreamErrorHandler = MeshErrorIO;
  //  atexit( forDebug);
  double raison = 0;
  verbosity = 9;
  set_new_handler(&NewHandler);
  char *fmeshback = "1.mesh";
  char *fmetrix = 0;
  // ostream *f_metr(&cout);
  if (argc >= 2) fmeshback = argv[1];
  if (argc >= 3) fmetrix = argv[2];
  if (argc >= 4) raison = atof(argv[3]);

#ifdef DRAWING
  initgraphique( );
  initgraph = 1;

#endif
  cout << "open " << fmeshback << " raison = " << raison << endl;
  Triangles TTh(fmeshback);
  Real8 hmin = TTh.MinimalHmin( );
  Real8 hmax = TTh.MaximalHmax( );
  Real8 coef = 1;
  if (fmetrix) TTh.ReadMetric(fmetrix, hmin, hmax, coef);

#ifdef DRAWING

  reffecran( );
  TTh.InitDraw( );
  TTh.inquire( );
  if (raison && fmetrix) {
    TTh.inquire( );
    TTh.SmoothMetric(raison);
  }
  if (fmetrix)
    for (Int4 i = 0; i < TTh.nbv; i++) ((Metric)TTh.vertices[i]).Draw(TTh.vertices[i]);
  TTh.Draw( );
  TTh.inquire( );
  closegraphique( );
#endif
  return (0);
}
