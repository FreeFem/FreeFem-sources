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

#include "ff++.hpp"

using namespace Fem2D;

const Mesh *Split_Edges(Stack stack, Fem2D::Mesh const *const &pTh, long *dK) {
  bool v10 = verbosity > 10;
  bool vp = verbosity > 1;

  assert(pTh);
  const Mesh &Th(*pTh);    // le maillage d'origne a decoupe
  using Fem2D::BoundaryEdge;
  using Fem2D::Mesh;
  using Fem2D::R2;
  using Fem2D::Triangle;
  using Fem2D::Vertex;
  int nbtn = Th.nt;
  int nbvn = Th.nv;
  int nebn = Th.neb;

  for (int k = 0; k < Th.nt; ++k) {
    for (int e = 0; e < 3; ++e) {
      if (dK[k] & (1 << e)) {
        nbtn++;
        int ee(e), kk;
        kk = Th.ElementAdj(k, ee);
        if (kk <= k) {
          nbvn++;
        }
      }
    }
  }

  // bug corrige ????
  for (int i = 0; i < Th.neb; i++) {
    int ek, k = Th.BoundaryElement(i, ek);
    if (dK[k] & (1 << ek)) {
      nebn++;
    }
  }

  if (vp) {
    cout << " Split_Edges:  nbv " << nbvn << " nbtn = " << nbtn << " nebn = " << nebn << endl;
  }

  int nbv = Th.nv;     // nombre de sommet
  int nbt = Th.nt;     // nombre de triangles
  int neb = Th.neb;    // nombre d'aretes fontiere
  // allocation des nouveaux items du maillage
  Vertex *v = new Vertex[nbvn];
  Triangle *t = new Triangle[nbtn];
  BoundaryEdge *b = new BoundaryEdge[nebn];
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

  KNM< int > NGP(3, nbt);

  // generation des points barycentre dearete a decoupe
  for (int k = 0; k < nbt; k++) {
    for (int e = 0; e < 3; ++e) {
      if (dK[k] & (1 << e)) {
        int ee(e), kk;
        kk = Th.ElementAdj(k, ee);
        if (kk <= k) {
          Triangle &K = Th[k];
          int i0 = (e + 1) % 3;
          int i1 = (e + 2) % 3;
          R2 A = ((R2)K[i0] + K[i1]) / 2.;
          vv->x = A.x;
          vv->y = A.y;
          vv->lab = 0;
          int j = vv - v;
          NGP(e, k) = j;
          if ((kk != k) && (kk >= 0)) {
            NGP(ee, kk) = j;
          }

          vv++;
        }
      }
    }
  }

  ffassert(vv - v == nbvn);
  // generation des triangles
  Triangle *tt = t;
  int nberr = 0;

  for (int k = 0; k < nbt; k++) {
    Triangle &K = Th[k];
    int j0 = Th(k, 0), j1 = Th(k, 1), j2 = Th(k, 2);
    int p[] = {0, 1, 2};
    R Le[] = {K.lenEdge2(0), K.lenEdge2(1), K.lenEdge2(2)};
    if (Le[p[0]] < Le[p[1]]) {
      Exchange(p[0], p[1]);
    }

    if (Le[p[1]] < Le[p[2]]) {
      Exchange(p[1], p[2]);
    }

    if (Le[p[0]] < Le[p[1]]) {
      Exchange(p[0], p[1]);
    }

    if (v10) {
      cout << k << " \t " << j0 << " " << j1 << " " << j2 << " ------ " << tt - t << endl;
    }

    Triangle *t0[] = {tt, tt, tt};
    (*tt++).set(v, j0, j1, j2, Th[k].lab);

    for (int ie = 0; ie < 3; ++ie) {
      int e = p[ie];
      int e1 = (e + 1) % 3, e2 = (e + 2) % 3;
      Triangle *td = t0[e], &Kd = *td;
      ffassert(td);

      if (dK[k] & (1 << e)) {
        Triangle *tn = tt++;
        int iee = NGP(e, k);
        int id[] = {int(&(Kd[0]) - v), int(&(Kd[1]) - v), int(&(Kd[2]) - v)};
        int in[] = {id[0], id[1], id[2]};
        id[e1] = iee;
        in[e2] = iee;

        if (v10) {
          cout << k << " \t " << in[0] << " " << in[1] << " " << in[2] << "  nn " << tn - t << " "
               << Le[e] << endl;
        }

        (*td).set(v, id[0], id[1], id[2], Th[k].lab);
        if (v10) {
          cout << k << " \t " << id[0] << " " << id[1] << " " << id[2] << "  dd " << td - t << endl;
        }

        (*tn).set(v, in[0], in[1], in[2], Th[k].lab);
        t0[e] = 0;    // done ..
        if (t0[e1] && t0[e2]) {
          t0[e1] = td;
          t0[e2] = tn;
        }
      }
    }
  }

  ffassert(tt - t == nbtn);
  // les arete frontieres qui n'ont pas change
  BoundaryEdge *bb = b;

  for (int i = 0; i < neb; i++) {
    int ek, k = Th.BoundaryElement(i, ek);
    int i1 = Th(Th.bedges[i][0]);
    int i2 = Th(Th.bedges[i][1]);
    int lab = Th.bedges[i].lab;
    if (dK[k] & (1 << ek)) {
      int iee = NGP(ek, k);
      assert(iee > 0);
      if (v10) {
        cout << " " << i1 << " " << iee << " " << i2 << " " << lab << " " << endl;
      }

      *bb++ = BoundaryEdge(v, i1, iee, lab);
      *bb++ = BoundaryEdge(v, iee, i2, lab);
    } else {
      *bb++ = BoundaryEdge(v, i1, i2, lab);
    }
  }

  ffassert(bb - b == nebn);
  // generation de la class Mesh a partir des 3 tableaux : v,t,b
  {
    Mesh *m = new Mesh(nbvn, nbtn, nebn, v, t, b);
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
class SplitEdges : public E_F0mps {
 public:
  typedef pmesh Result;
  Expression expTh;
  Expression spt;

  SplitEdges(const basicAC_F0 &args) {
    args.SetNameParam( );
    expTh = to< pmesh >(args[0]);    // a the expression to get the mesh
    spt = to< double >(args[1]);     // a the expression to get the mesh
  }

  static ArrayOfaType typeargs( ) { return ArrayOfaType(atype< pmesh >( ), atype< double >( )); }

  static E_F0 *f(const basicAC_F0 &args) { return new SplitEdges(args); }

  AnyType operator( )(Stack s) const;
};

AnyType SplitEdges::operator( )(Stack stack) const {
  const Mesh *pTh = GetAny< pmesh >((*expTh)(stack));
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  bool v10 = verbosity > 10;

  ffassert(pTh);
  const Mesh &Th(*pTh);
  KN< long > dK(Th.nt);
  dK = 0L;
  int ne = 0;

  for (int k = 0; k < Th.nt; k++) {
    for (int e = 0; e < 3; e++) {
      Triangle &K = Th[k];
      int e1 = (e + 1) % 3;
      int e2 = (e + 2) % 3;
      R2 P1 = K[e1], P2 = K[e2];
      R2 P = (P1 + P2) / 2.;
      MeshPointStack(stack)->set(P.x, P.y);
      double de = fabs(GetAny< double >((*spt)(stack)));
      bool be = fabs(de) > 1e-30;
      if (be) {
        dK[k] += (1 << e);
        ne++;
      }

      if (v10) {
        cout << k << " " << e << "   f " << P << " = " << de << " " << be << " " << dK[k] << " "
             << (1 << e) << endl;
      }

      int ee(e), kk;
      kk = Th.ElementAdj(k, ee);
      if ((kk < k) && (kk >= 0)) {
        bool bee = dK[kk] & (1 << ee);
        if (bee != be) {
          cout << " Bizarre edge right != compatible left " << k << " " << e << " P = " << P
               << " kk " << kk << " " << ee << " " << dK[kk] << endl;
          dK[k] = dK[k] | (1 << e);
          dK[kk] = dK[kk] | (1 << ee);
        }
      }
    }
  }

  if (verbosity > 0) {
    cout << "  SplitEdges: nb split edge = " << ne << endl;
  }

  *mp = mps;
  return SetAny< pmesh >(Split_Edges(stack, pTh, (long *)dK));
}

static void Load_Init( ) {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
  if (verbosity) {
    cout << " lood: Split3  " << endl;
  }

  Global.Add("SplitEdgeMesh", "(", new OneOperatorCode< SplitEdges >( ));
  // utilisation
  // mesh Th,Th3;
  // ... construction du maillage Th ici
  // Th3=splitmesh3(Th);
  /*  example complet : splitmesh3.edp
   *  load "splitedge"
   *  mesh Th=square(5,5);
   *  mesh Th3=SplitEdgeMesh(Th,x<0.51 && y < 0.49 );
   *  plot(Th3,wait=1);
   */
}

LOADFUNC(Load_Init)
