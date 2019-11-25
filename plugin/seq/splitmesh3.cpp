/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include <iostream>
#include <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;
#include "rgraph.hpp"
#include "RNM.hpp"
#include <fem.hpp>
#include <cmath>

using namespace Fem2D;

Mesh const *SplitMesh3(Stack stack, Fem2D::Mesh const *const &pTh) {
  assert(pTh);
  const Mesh &Th(*pTh);    // le maillage d'origne a decoupe
  using Fem2D::BoundaryEdge;
  using Fem2D::Mesh;
  using Fem2D::R2;
  using Fem2D::Triangle;
  using Fem2D::Vertex;
  int nbv = Th.nv;     // nombre de sommet
  int nbt = Th.nt;     // nombre de triangles
  int neb = Th.neb;    // nombre d'aretes fontiere
  // allocation des nouveaux items du maillage
  Vertex *v = new Vertex[nbv + nbt];
  Triangle *t = new Triangle[nbt * 3];
  BoundaryEdge *b = new BoundaryEdge[neb];
  // generation des nouveaus sommets
  Vertex *vv = v;

  // copie des anciens sommets (remarque il n'y a pas operateur de copy des sommets)
  for (int i = 0; i < nbv; i++) {
    Vertex &V = Th(i);
    vv->x = V.x;
    vv->y = V.y;
    vv->lab = V.lab;
    vv++;
  }

  // generation des points barycentre de trianngles
  for (int k = 0; k < nbt; k++) {
    Triangle &K = Th[k];
    R2 G = ((R2)K[0] + K[1] + K[2]) / 3.;
    vv->x = G.x;
    vv->y = G.y;
    vv->lab = 0;
    vv++;
  }

  // generation des triangles
  Triangle *tt = t;
  int nberr = 0;

  for (int i = 0; i < nbt; i++) {
    int i0 = Th(i, 0), i1 = Th(i, 1), i2 = Th(i, 2);
    int ii = nbv + i;    // numero du
    // les 3 triangles par triangles origines
    (*tt++).set(v, ii, i1, i2, Th[i].lab);
    (*tt++).set(v, i0, ii, i2, Th[i].lab);
    (*tt++).set(v, i0, i1, ii, Th[i].lab);
  }

  // les arete frontieres qui n'ont pas change
  BoundaryEdge *bb = b;

  for (int i = 0; i < neb; i++) {
    int i1 = Th(Th.bedges[i][0]);
    int i2 = Th(Th.bedges[i][1]);
    int lab = Th.bedges[i].lab;
    *bb++ = BoundaryEdge(v, i1, i2, lab);
  }

  // generation de la class Mesh a partir des 3 tableaux : v,t,b
  {
    Mesh *m = new Mesh(nbv + nbt, nbt * 3, neb, v, t, b);
    R2 Pn, Px;
    m->BoundingBox(Pn, Px);
    m->quadtree = new Fem2D::FQuadTree(m, Pn, Px, m->nv);
    // m->decrement();
    Add2StackOfPtr2FreeRC(stack, m);
    return m;
  }
}

// truc pour que la fonction
// static void Load_Init() soit appele a moment du chargement dynamique
// du fichier
//

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  if (verbosity) {
    cout << " lood: Split3  " << endl;
  }

  Global.Add("splitmesh3", "(", new OneOperator1s_< Mesh const *, Mesh const * >(SplitMesh3));
  // utilisation
  // mesh Th,Th3;
  // ... construction du maillage Th ici
  // Th3=splitmesh3(Th);
  /*  example complet : splitmesh3.edp
   *  load "splitmesh3"
   *  mesh Th=square(5,5);
   *  mesh Th3=splitmesh3(Th);
   *  plot(Th3,wait=1);
   */
}

LOADFUNC(Load_Init)
