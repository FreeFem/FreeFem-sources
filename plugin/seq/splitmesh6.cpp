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
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep:
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

Mesh const *SplitMesh6(Stack stack, Fem2D::Mesh const *const &pTh) {
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
  int nbe = 0;

  for (int k = 0; k < nbt; ++k) {
    for (int e = 0; e < 3; ++e) {
      int ee = e;
      int kk = Th.ElementAdj(k, ee);
      if (kk <= k) {
        nbe++;
      }
    }
  }

  Vertex *v = new Vertex[nbv + nbt + nbe];
  Triangle *t = new Triangle[nbt * 6];
  BoundaryEdge *b = new BoundaryEdge[neb * 2];
  // generation des nouveaus sommets
  Vertex *vv = v;
  KN< int > mm(3 * nbt);

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

  // generation des milieux des cote
  int nn = 0;

  for (int k = 0; k < nbt; ++k) {
    for (int e = 0; e < 3; ++e) {
      int ee = e;
      int kk = Th.ElementAdj(k, ee);
      if ((kk >= k) || (kk < 0)) {
        int v0 = Th(k, EdgesVertexTriangle[e][0]);
        int v1 = Th(k, EdgesVertexTriangle[e][1]);
        R2 M = ((R2)Th(v0) + Th(v1)) / 2.;
        int lab = 0;
        BoundaryEdge *be = Th.TheBoundaryEdge(v0, v1);
        if (be) {
          lab = be->lab;
        }

        vv->x = M.x;
        vv->y = M.y;
        vv->lab = lab;

        mm[k * 3 + e] = vv - v;    // numero du sommet
        vv++;
        nn++;
      } else {
        mm[k * 3 + e] = mm[kk * 3 + ee];
      }
    }
  }

  cout << " nb edge = " << nbe << " == " << nn << endl;
  ffassert(nbe == nn);

  // generation des triangles
  Triangle *tt = t;
  int nberr = 0;

  for (int i = 0; i < nbt; i++) {
    int i0 = Th(i, 0), i1 = Th(i, 1), i2 = Th(i, 2);
    int j0 = mm[i * 3], j1 = mm[i * 3 + 1], j2 = mm[i * 3 + 2];
    int ii = nbv + i;    // numero du
    // les 3 triangles par triangles origines
    (*tt++).set(v, ii, i1, j0, Th[i].lab);
    (*tt++).set(v, ii, j0, i2, Th[i].lab);
    (*tt++).set(v, i0, ii, j1, Th[i].lab);
    (*tt++).set(v, j1, ii, i2, Th[i].lab);
    (*tt++).set(v, i0, j2, ii, Th[i].lab);
    (*tt++).set(v, j2, i1, ii, Th[i].lab);
  }

  // les arete frontieres qui n'ont pas change
  BoundaryEdge *bb = b;

  for (int i = 0; i < neb; i++) {
    int ki;
    int k = Th.BoundaryElement(i, ki);
    int i1 = Th(Th.bedges[i][0]);
    int i2 = Th(Th.bedges[i][1]);
    int ii = mm[3 * k + ki];
    int lab = Th.bedges[i].lab;
    *bb++ = BoundaryEdge(v, i1, ii, lab);
    *bb++ = BoundaryEdge(v, ii, i2, lab);
  }

  // generation de la class Mesh a partir des 3 tableaux : v,t,b
  {
    Mesh *m = new Mesh(nbv + nbt + nbe, nbt * 6, neb * 2, v, t, b);
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
    cout << " lood: Split6  " << endl;
  }

  Global.Add("splitmesh6", "(", new OneOperator1s_< Mesh const *, Mesh const * >(SplitMesh6));
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
