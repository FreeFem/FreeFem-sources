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

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

#include <iostream>
#include <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
using namespace std;
#include "rgraph.hpp"
#include "RNM.hpp"
#include "ff++.hpp"
#include "AFunction_ext.hpp"    // [[file:../src/fflib/AFunction_ext.hpp]]
#include "msh3.hpp"
#include <fem.hpp>
#include <cmath>

using namespace Fem2D;

static void findPerm(int *ivt, int *pivt, Vertex3 *v) {
  std::copy(ivt, ivt + 4, pivt);
  R3 AB(v[ivt[0]], v[ivt[1]]);
  R3 AC(v[ivt[0]], v[ivt[2]]);
  R3 AD(v[ivt[0]], v[ivt[3]]);
  if (det(AB, AC, AD) > 0) {
    return;
  } else if (det(AB, AD, AC) > 0) {
    std::swap(pivt[2], pivt[3]);
    return;
  } else if (det(AC, AB, AD) > 0) {
    std::swap(pivt[1], pivt[2]);
    return;
  }
}

Mesh3 const *SplitMesh12(Stack stack, Fem2D::Mesh3 const *const &pTh) {
  assert(pTh);
  const Mesh3 &Th(*pTh);    // le maillage d'origne a decoupe
  using Fem2D::BoundaryEdge;
  using Fem2D::Mesh;
  using Fem2D::R2;
  using Fem2D::Triangle;
  using Fem2D::Vertex;
  int nbv = Th.nv;     // nombre de sommet
  int nbt = Th.nt;     // nombre de triangles
  int nbe = Th.nbe;    // nombre d'aretes fontiere
  // allocation des nouveaux items du maillage
  int nbnv = 0;

  for (int k = 0; k < nbt; ++k) {
    for (int e = 0; e < 4; ++e) {
      int ee = e;
      int kk = Th.ElementAdj(k, ee);
      if (kk > k) {
        nbnv++;
      }
    }
  }

  Vertex3 *v = new Vertex3[nbv + nbt + nbe + nbnv];
  Tet *t = new Tet[nbt * 12];
  Triangle3 *b = new Triangle3[nbe * 3];
  // generation des nouveaus sommets
  Vertex3 *vv = v;

  // copie des anciens sommets (remarque il n'y a pas operateur de copy des sommets)
  for (int i = 0; i < nbv; i++) {
    const Vertex3 &V = Th(i);
    vv->x = V.x;
    vv->y = V.y;
    vv->z = V.z;
    vv->lab = V.lab;
    vv++;
  }

  // generation des points barycentre de trianngles
  for (int k = 0; k < nbt; k++) {
    const Tet &K = Th[k];
    R3 G = ((R3)K[0] + K[1] + K[2] + K[3]) / 4.;
    vv->x = G.x;
    vv->y = G.y;
    vv->z = G.z;
    vv->lab = 0;
    vv++;
  }

  for (int i = 0; i < nbe; i++) {
    const Triangle3 &K(Th.be(i));
    R3 G = ((R3)K[0] + K[1] + K[2]) / 3.;
    vv->x = G.x;
    vv->y = G.y;
    vv->z = G.z;
    vv->lab = Th.be(i).lab;
    vv++;
  }

  // generation des triangles
  Tet *tt = t;

  for (int i = 0; i < nbe; i++) {
    int ki;
    int k = Th.BoundaryElement(i, ki);
    const Triangle3 &K(Th.be(i));
    int i0 = Th.operator( )(K[0]), i1 = Th.operator( )(K[1]), i2 = Th.operator( )(K[2]);
    int ii = nbv + nbt + i;    // numero du
    int jj = nbv + k;          // numero du
    int ivt[4] = {jj, i1, i2, ii};
    int pivt[4];
    findPerm(ivt, pivt, v);
    (*tt++).set(v, pivt, K.lab);
    ivt[0] = i0;
    ivt[1] = jj;
    findPerm(ivt, pivt, v);
    (*tt++).set(v, pivt, K.lab);
    ivt[1] = i1;
    ivt[2] = jj;
    findPerm(ivt, pivt, v);
    (*tt++).set(v, pivt, K.lab);
  }

  const int nvfaceTet[4][3] = {{3, 2, 1}, {0, 2, 3}, {3, 1, 0}, {0, 1, 2}};

  for (int k = 0; k < nbt; k++) {
    for (int e = 0; e < 4; ++e) {
      int ee = e;
      int kk = Th.ElementAdj(k, ee);
      if (kk > k) {
        const Tet &K = Th[k];
        const Tet &KAdj = Th[kk];
        {
          R3 uuu = K[nvfaceTet[e][1]] - K[nvfaceTet[e][0]],
             vvv = K[nvfaceTet[e][2]] - K[nvfaceTet[e][0]];
          R3 n = uuu ^ vvv;
          Vertex3 dir(v[nbv + k] - v[nbv + kk]);
          Vertex3 w0(v[nbv + kk] - K[nvfaceTet[e][0]]);
          double aa = -(n, w0);
          double bb = (n, dir);
          double rr = aa / bb;
          vv->x = v[nbv + kk].x + rr * dir.x;
          vv->y = v[nbv + kk].y + rr * dir.y;
          vv->z = v[nbv + kk].z + rr * dir.z;
          vv->lab = 0;
        }
        int i0 = Th.operator( )(K[nvfaceTet[e][0]]), i1 = Th.operator( )(K[nvfaceTet[e][1]]),
            i2 = Th.operator( )(K[nvfaceTet[e][2]]);

        for (int ij = 0; ij < 2; ++ij) {
          int ivt[4] = {ij == 0 ? nbv + k : nbv + kk, i1, i2, static_cast< int >(vv - v)};
          int lab = ij == 0 ? K.lab : KAdj.lab;
          int pivt[4];
          findPerm(ivt, pivt, v);
          (*tt++).set(v, pivt, lab);
          ivt[1] = ivt[0];
          ivt[0] = i0;
          findPerm(ivt, pivt, v);
          (*tt++).set(v, pivt, lab);
          ivt[2] = ivt[1];
          ivt[1] = i1;
          findPerm(ivt, pivt, v);
          (*tt++).set(v, pivt, lab);
        }

        vv++;
      }
    }
  }

  // les arete frontieres qui n'ont pas change
  Triangle3 *bb = b;

  for (int i = 0; i < nbe; i++) {
    const Triangle3 &K(Th.be(i));
    int orig[3];
    int ivv[3];

    orig[0] = Th.operator( )(K[0]);
    orig[1] = Th.operator( )(K[1]);
    orig[2] = Th.operator( )(K[2]);
    ivv[0] = nbv + nbt + i;
    ivv[1] = orig[1];
    ivv[2] = orig[2];
    (bb++)->set(v, ivv, K.lab);
    ivv[1] = ivv[0];
    ivv[0] = orig[0];
    (bb++)->set(v, ivv, K.lab);
    ivv[2] = ivv[1];
    ivv[1] = orig[1];
    (bb++)->set(v, ivv, K.lab);
  }

  // generation de la class Mesh a partir des 3 tableaux : v,t,b
  {
    Mesh3 *m = new Mesh3(nbv + nbt + nbe + nbnv, nbt * 12, nbe * 3, v, t, b);
    m->BuildGTree( );
    // m->decrement();
    Add2StackOfPtr2FreeRC(stack, m);
    return m;
  }
}

// truc pour que la fonction
// static void Load_Init() soit appele a moment du chargement dynamique
// du fichier
//

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh12"  a freefem++
  if (verbosity > 1) {
    cout << " load: Split12  " << endl;
  }

  Global.Add("splitmesh12", "(", new OneOperator1s_< Mesh3 const *, Mesh3 const * >(SplitMesh12));
}

LOADFUNC(Load_Init)
