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
/* add
 add Powell-Sabin division:
 – Center of triangle is taken the center of the inscribed circle
 – Points on edges are obtained by intersecting the edge with the line joining the centers of the two and on boundary just barycenter .


 
 */
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
// seg intersection in 2d
R2 Intersection(R2 B,R2 BB,R2 A,R2 AA)
{
    double b = -det(A, AA, B);
    double bb = det(A, AA, BB);
    double s = b + bb;
    assert( b/s > 0 && bb/s >=0 );
    R2 P = BB * b / s + B* bb / s;
    return P;

}
// c
R2 incircleCenter(R2 A, R2 B, R2 C)
{
    R2 AB(A,B),CA(C,A),BC(B,C);
    double lc=AB.norme(),lb=CA.norme(),la=BC.norme();
    R2 ABu =AB/lc, CAu=CA/lb, BCu = BC/la;
    R2 AA = A + (ABu-CAu);
    R2 BB = B + (BCu-ABu);
    R2 CC = C + (CAu-BCu);
    //  inter (A,AA) et ((B,BB)
    // intersection ddroite [a,aa] and [b,bb]
    R2 G= Intersection(A,AA,B,BB);
    assert(det(A,B,G)>0);
    assert(det(A,G,C)>0);
    assert(det(G,B,C)>0);
    return G;
 }

Mesh const *SplitMesh6New(Stack stack, Fem2D::Mesh const *const &pTh,int flags) {
  //  flags == 0 basis
    //  flags == 1  Powell-Sabin refinend
    
  assert(pTh);
    bool PowellSabin = flags ==1;
    if(verbosity>1)
         cout << "SplitMesh6New "<<flags <<endl;
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

  // generation des points barycentre de trianngles o c
  for (int k = 0; k < nbt; k++) {
    Triangle &K = Th[k];
      R2 G;
      if(PowellSabin)
          G = incircleCenter(K[0] ,K[1] ,K[2] );
      else
          G = ((R2)K[0] + K[1] + K[2]) / 3.;
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
          R2 M;
        if(PowellSabin && kk >0 && k !=kk)
            M = Intersection(Th(v0),Th(v1),v[nbv+k],v[nbv+kk]);
        else
          M = ((R2)Th(v0) + Th(v1)) / 2.;
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
  if(verbosity>9)
  cout << " SplitMesh6New:: nb edge = " << nbe << " == " << nn << endl;
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
    m->renum(); 
    R2 Pn, Px;
    m->BoundingBox(Pn, Px);
    m->quadtree = new Fem2D::FQuadTree(m, Pn, Px, m->nv);
    // m->decrement();
    Add2StackOfPtr2FreeRC(stack, m);
    return m;
  }
}
Mesh const *SplitMesh6(Stack stack, Fem2D::Mesh const *const &pTh) 
{
    return SplitMesh6New(stack,pTh,0);
}
Mesh const *SplitPowellSabin(Stack stack, Fem2D::Mesh const *const &pTh)
{
    return SplitMesh6New(stack,pTh,1);
}

// truc pour que la fonction
// static void Load_Init() soit appele a moment du chargement dynamique
// du fichier
//
//  PowellSabin
static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  if (verbosity) {
    cout << " load: splitmesh6  " << endl;
  }

    Global.Add("splitmesh6", "(", new OneOperator1s_< Mesh const *, Mesh const * >(SplitMesh6));
    Global.Add("splitmesh6PowellSabin", "(", new OneOperator1s_< Mesh const *, Mesh const * >(SplitPowellSabin));
   
     // utilisation
  // mesh Th,Th3;
  // ... construction du maillage Th ici
  // Th3=splitmesh3(Th);
  /*  example complet : splitmesh3.edp
   *  load "splitmesh3"
   *  mesh Th=square(5,5);
   *  mesh Th6=splitmesh6(Th); //  split with barycentrer
   *  mesh ThPS=splitmesh6PowellSabin(Th); //  split with PowellSabin point
   
 
   *  plot(Th3,wait=1);
   */
}

LOADFUNC(Load_Init)
