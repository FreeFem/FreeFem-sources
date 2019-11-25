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
// AUTHORS : Jacques Morice
// Frederic Hecht
// E-MAIL  : jacques.morice@ann.jussieu.fr
// frederic.hecht@sorbonne-universite.fr

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

/*
 * Thank to the ARN ()  FF2A3 grant
 * ref:ANR-07-CIS7-002-01
 */

/*
 * calcul demander par F. Hecht
 */

#ifndef WITH_NO_INIT
#include "ff++.hpp"
#include "AFunction_ext.hpp"
#endif

using namespace std;

#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <queue>
// #include "msh3.hpp"

using namespace Fem2D;

// fonction determinant les points d'intersection
static int debug = 0;
class FINDLOCALMIN_P1_Op : public E_F0mps {
 public:
  Expression eTh, eu, er;
  static const int n_name_param = 2;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];
  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  bool arg(int i, Stack stack, bool a) const {
    return nargs[i] ? GetAny< bool >((*nargs[i])(stack)) : a;
  }

  KN< long > *arg(int i, Stack stack, KN< long > *a) const {
    return nargs[i] ? GetAny< KN< long > * >((*nargs[i])(stack)) : a;
  }

  string *arg(int i, Stack stack, string *a) const {
    return nargs[i] ? GetAny< string * >((*nargs[i])(stack)) : a;
  }

 public:
  FINDLOCALMIN_P1_Op(const basicAC_F0 &args, Expression tth, Expression fu, Expression fr)
    : eTh(tth), eu(fu), er(fr) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type FINDLOCALMIN_P1_Op::name_param[] = {{"eps", &typeid(double)},
                                                              {"convex", &typeid(long)}};

class ISOLINE_P1_Op : public E_F0mps {
 public:
  Expression eTh, eff, emat, exx, eyy, exy, iso;
  static const int n_name_param = 7;    //
  static basicAC_F0::name_and_type name_param[];
  Expression nargs[n_name_param];

  double arg(int i, Stack stack, double a) const {
    return nargs[i] ? GetAny< double >((*nargs[i])(stack)) : a;
  }

  long arg(int i, Stack stack, long a) const {
    return nargs[i] ? GetAny< long >((*nargs[i])(stack)) : a;
  }

  KN< long > *arg(int i, Stack stack, KN< long > *a) const {
    return nargs[i] ? GetAny< KN< long > * >((*nargs[i])(stack)) : a;
  }

  string *arg(int i, Stack stack, string *a) const {
    return nargs[i] ? GetAny< string * >((*nargs[i])(stack)) : a;
  }

 public:
  ISOLINE_P1_Op(const basicAC_F0 &args, Expression tth, Expression fff, Expression xxx,
                Expression yyy)
    : eTh(tth), eff(fff), emat(0), exx(xxx), eyy(yyy), exy(0) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  ISOLINE_P1_Op(const basicAC_F0 &args, Expression tth, Expression fff, Expression xxyy)
    : eTh(tth), eff(fff), emat(0), exx(0), eyy(0), exy(xxyy) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type ISOLINE_P1_Op::name_param[] = {
  {"iso", &typeid(double)},   {"close", &typeid(long)}, {"smoothing", &typeid(double)},
  {"ratio", &typeid(double)}, {"eps", &typeid(double)}, {"beginend", &typeid(KN< long > *)},
  {"file", &typeid(string *)}};

int IsoLineK(R2 *P, double *f, R2 *Q, int *i0, int *i1, double eps) {
  int kv = 0, ke = 0, e = 3;
  int tv[3], te[3], vk[3];

  for (int i = 0; i < 3; ++i) {
    if (abs(f[i]) <= eps) {
      e -= tv[kv++] = i;
      vk[i] = 1;
    } else {
      vk[i] = 0;
    }
  }

  if (debug) {
    cout << " ** " << kv << endl;
  }

  if (kv > 1) {    // on 2  vertex on the isoline ....
    if (kv == 2) {
      if (f[e] > 0.) {
        int j0 = (e + 1) % 3;
        int j1 = (e + 2) % 3;
        te[ke] = e + 3, i0[ke] = j0, i1[ke] = j0, ++ke;
        te[ke] = e, i0[ke] = j1, i1[ke] = j1, ++ke;
        // pb d'unicity, need to see the adj triangle ...
        // return 10+e ; // edge number + 10
      } else {
        return 0;    // skip edge ...
      }
    } else {
      return 0;    // const funct...
    }
  } else {    // see internal edge ..
    for (int e = 0; e < 3; ++e) {
      int j0 = (e + 1) % 3;
      int j1 = (e + 2) % 3;
      if (vk[j0]) {    // the intial  point on iso line
        if (0. < f[j1]) {
          te[ke] = e, i0[ke] = j0, i1[ke] = j0, ++ke;
        } else {
          te[ke] = e + 3, i0[ke] = j0, i1[ke] = j0, ++ke;
        }
      } else if (vk[j1]) {
        ;                                       // skip the final point on iso line
      } else if (f[j0] < 0. && 0. < f[j1]) {    // good  sens
        te[ke] = e, i0[ke] = j0, i1[ke] = j1, ++ke;
      } else if (f[j0] > 0. && 0. > f[j1]) {    // inverse  sens
        te[ke] = e + 3, i0[ke] = j1, i1[ke] = j0, ++ke;
      }
    }
  }

  if (ke == 2) {
    // the  K[i1[0]] , Q[0], Q[1] must be direct ...
    // the  K[i0[1]] , Q[0], Q[1] must be direct ...
    // Warning   no trivail case ..  make a plot to see
    // with is good
    // the first edge must be

    if (te[0] < 3) {    // oriente the line
      assert(te[1] >= 3);
      std::swap(te[0], te[1]);
      std::swap(i0[0], i0[1]);
      std::swap(i1[0], i1[1]);
      if (debug) {
        cout << " swap " << endl;
      }
    }

    for (int i = 0; i < 2; ++i) {
      int j0 = i0[i], j1 = i1[i];
      if (j0 == j1) {
        Q[i] = P[j0];
      } else {
        Q[i] = (P[j0] * (f[j1]) - P[j1] * (f[j0])) / (f[j1] - f[j0]);
      }

      if (debug) {
        cout << i << " " << j0 << " " << j1 << " : " << Q[i] << "***" << endl;
      }
    }

    if (debug) {
      cout << "i0 " << i0[0] << " " << i0[1] << " " << det(P[i1[0]], Q[0], Q[1]) << endl;
      cout << "i1 " << i1[0] << " " << i1[1] << " " << det(P[i0[1]], Q[1], Q[0]) << endl;
      cout << "f " << f[0] << " " << f[1] << " " << f[2] << endl;
      cout << "P " << P[0] << ", " << P[1] << ", " << P[2] << endl;
      cout << "Q " << Q[0] << ", " << Q[1] << endl;
    }

    if (!vk[i1[0]]) {
      assert(det(P[i1[0]], Q[0], Q[1]) > 0);
    }

    if (!vk[i0[1]]) {
      assert(det(P[i0[1]], Q[1], Q[0]) > 0);
    }

    return 2;
  }

  // remark, the left of the line is upper .
  return 0;
}

int LineBorder(R2 *P, double *f, long close, R2 *Q, int *i1, int *i2, double eps) {
  int np = 0;

  if (close) {
    if (f[0] > -eps) {
      Q[np] = P[0];
      i1[np] = i2[np] = 0, np++;
    }

    if (f[0] * f[1] <= -eps * eps) {
      Q[np] = (P[0] * (f[1]) - P[1] * (f[0])) / (f[1] - f[0]);
      i1[np] = 0, i2[np] = 1, np++;
    }

    if (f[1] > -eps) {
      Q[np] = P[1];
      i1[np] = i2[np] = 1, np++;
    }
  }

  return np;
}

struct R2_I2 {
  R2 P;
  int nx;
  R2_I2(R2 A, int nxx = -1) : P(A), nx(nxx) {}

  bool add(int k0, int k1, multimap< int, int > &L) {
    if (nx == -1) {
      nx = k1;
    } else {
      if (nx > 0) {    // more than 2 seg ... put data in  the multi map ..
        L.insert(make_pair(k0, nx));
        L.insert(make_pair(k0, k1));
        nx = -2;
      } else {
        L.insert(make_pair(k0, k1));
      }

      return false;
    }

    return true;
  }

  int next(int k0, multimap< int, int > &L, int rm = 0) {
    int nxx = -1;

    if (nx >= 0) {
      nxx = nx;
      if (rm) {
        nx = -2;
      }
    } else {
      typedef multimap< int, int >::iterator IT;
      IT f = L.find(k0);
      if (f == L.end( )) {
        nxx = -1;    //
      } else {
        nxx = f->second;
        if (rm) {
          L.erase(f);
        }
      }
    }

    return nxx;
  }

  int count(int k0, const multimap< int, int > &L) const {
    if (nx >= 0) {
      return 1;
    } else {
      return L.count(k0);
    }
  }
};
// Absact mesh for data on grid .... FH ..
int Th_Grid(const KNM_< double > *g, int k, int ii) {
  int N = g->N( ) - 1;
  int kq = k / 2;    // number of the quad
  int k0 = k % 2;    // up or down
  int I = kq % N + (k0 ? (ii % 2) : (ii != 0));
  int J = kq / N + (k0 ? (ii != 0) : (ii == 2));

  return J * (N + 1) + I;
}

R2 V_Grid(const KNM_< double > *g, int k) {
  int i = k % g->N( ), j = k / g->N( );

  return R2(i, j);
}

int EA_Grid(const KNM_< double > *g, int k, int &e) {
  int kq = k / 2;    // number of the quad
  int k0 = k % 2;    // up or down
  bool intern = k0 ? (e == 0) : (e == 2);

  if (intern) {
    e = 2 - e;
    return 2 * kq + 1 - k0;
  }

  ffassert(0);
  return 0;
}

struct SMesh {
  const Mesh *pTh;
  const KNM_< double > *g;
  int nv, nt, neb;
  int operator( )(int k, int i) const { return pTh ? (*pTh)(k, i) : Th_Grid(g, k, i); }

  R2 operator( )(int i) const { return pTh ? (*pTh)(i) : V_Grid(g, i); }

  int ElementAdj(int k, int &e) { return pTh ? pTh->ElementAdj(k, e) : EA_Grid(g, k, e); }

  SMesh(const Mesh *PTh) : pTh(PTh), g(0), nv(pTh->nv), nt(pTh->nt), neb(pTh->neb) {}

  SMesh(KNM_< double > *gg)
    : pTh(0), g(gg), nv(gg->N( ) * gg->M( )), nt((gg->N( ) - 1) * (gg->M( ) - 1) * 2),
      neb((gg->N( ) + gg->M( ) - 2) * 2) {}
};
::AnyType FINDLOCALMIN_P1_Op::operator( )(Stack stack) const {
  typedef std::pair< double, int > KEY;
  typedef std::priority_queue< KEY, std::vector< KEY >, std::greater< KEY > > myPQ;
  typedef std::priority_queue< KEY > myPQL;

  const Mesh *pTh = GetAny< const Mesh * >((*eTh)(stack));
  ffassert(pTh);
  const Mesh &Th = *pTh;
  int ddd1 = verbosity > 9;
  int ddd0 = verbosity > 2;
  int nbv = Th.nv;    // nombre de sommet
  int nbt = Th.nt;    // nombre de triangles
  int convex = arg(1, stack, 0L);
  double eps = arg(0, stack, 0.);
  if (verbosity > 2) {
    cout << "    -- findlocalmin: convex = " << convex << " eps= " << eps << endl;
  }

  KN< long > rs(nbv);
  KN< double > *pu = GetAny< KN< double > * >((*eu)(stack));
  KN< double > *pr = GetAny< KN< double > * >((*er)(stack));
  KN< double > &U(*pu), &R(*pr);
  if (U.N( ) != nbv) {
    ffassert(0);
  }

  if (R.N( ) != nbt) {
    ffassert(0);
  }

  R = -1.;    //
  // 1 recher de min local
  rs = 1;
  // voisanage de sommet
  KN< long > head(nbv), next(3 * nbt);
  head = -1;

  for (int p = 0; p < next.N( ); ++p) {
    int j = Th(p / 3, p % 3);
    next[p] = head[j];
    head[j] = p;
  }

  KN< int > nbrv(nbv);    // nb traingle around  vertices
  nbrv = 0;

  for (int k = 0; k < nbt; ++k) {
    const Mesh::Element &K = Th[k];

    for (int i = 0; i < 3; ++i) {
      ++nbrv[Th(K[i])];
      int i0 = Th(K[(i + 1) % 3]);
      int i1 = Th(K[(i + 2) % 3]);
      if (U[i0] > U[i1]) {
        rs[i0] = 0;
      } else if (U[i1] > U[i0]) {
        rs[i1] = 0;
      }
    }
  }

  if (ddd0) {
    cout << "rq nb min =  " << rs.sum( ) << endl;
  }

  myPQ lvm;

  for (int i = 0; i < nbv; ++i) {
    if (rs[i] == 1) {
      lvm.push(make_pair(U[i], i));
    }
  }

  // analyse des minimals locaux
  KN< long > cv(nbv);
  cv = 0L;
  long col = 0;
  int nmin = 0;
  KN< long > sm(nbv);
  long rkm = 0;    // Numero de region

  while (!lvm.empty( )) {
    KEY tlvm = lvm.top( );
    lvm.pop( );
    double Ui = tlvm.first;
    int i = tlvm.second;

    if (rs[i] == 1 && nbrv[i] > 0) {
      ++col;    // change the color
      myPQ pqv;
      pqv.push(make_pair(Ui, i));

      if (ddd1) {
        cout << " ** " << i << " "
             << "ui =" << Ui << endl;
      }

      double Uvp = Ui;

      while (!pqv.empty( )) {
        KEY tp = pqv.top( );
        pqv.pop( );
        double Uv = tp.first;
        int iv = tp.second;
        if (ddd1) {
          cout << "\t\t" << iv << " " << Uv << endl;
        }

        assert(Uv >= Uvp);    // verif piority queue
        Uvp = Uv;

        for (int p = head[iv]; p >= 0; p = next[p]) {    // lpour les triangle autour de iv
          int k = p / 3, ik = p % 3, i1 = (ik + 1) % 3, i2 = (ik + 2) % 3;
          if (R[k] < 0) {    // nouveau triangle
            int j1 = Th(k, i1);
            int j2 = Th(k, i2);
            if (Uv <= U[j1] && Uv <= U[j2]) {
              R[k] = rkm;
              --nbrv[iv];
              --nbrv[j1];
              --nbrv[j2];

              if (ddd1) {
                cout << "\t\t\t Add " << k << " " << j1 << " " << U[j1] << " / " << j2 << " "
                     << U[j2] << endl;
              }

              if (nbrv[j1] > 0) {
                pqv.push(make_pair(U[j1], j1));
              }

              if (nbrv[j2] > 0) {
                pqv.push(make_pair(U[j2], j2));
              }
            }
          }
        }

        if (convex > 1) {
          break;
        }
      }

      sm[rkm++] = i;
    }
  }

  if (ddd1) {
    cout << " R = " << R << endl;
  }

  if (verbosity > 2) {
    cout << "    -- FindlocalMin nb min =" << nmin << endl;
  }

  sm.resize(rkm);
  KN< long > *ppr = new KN< long >(sm);
  return Add2StackOfPtr2Free(stack, ppr);
}

AnyType ISOLINE_P1_Op::operator( )(Stack stack) const {
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;

  KNM< double > *pxy = 0;
  KN< double > *pxx = 0;
  KN< double > *pyy = 0;
  if (exy) {
    pxy = GetAny< KNM< double > * >((*exy)(stack));
  }

  if (exx) {
    pxx = GetAny< KN< double > * >((*exx)(stack));
  }

  if (eyy) {
    pyy = GetAny< KN< double > * >((*eyy)(stack));
  }

  ffassert((pxx || pyy) == !pxy);
  const Mesh *pTh = GetAny< const Mesh * >((*eTh)(stack));
  ffassert(pTh);
  SMesh Th(pTh);
  int nbv = Th.nv;    // nombre de sommet
  int nbt = Th.nt;    // nombre de triangles
  // int nbe=Th.neb; // nombre d'aretes fontiere
  long nbc;
  // value of isoline
  double isovalue = arg(0, stack, 0.);
  long close = arg(1, stack, 1L);
  double smoothing = arg(2, stack, 0.);
  double ratio = arg(3, stack, 1.);
  double epsr = arg(4, stack, 1e-10);
  KN< long > *pbeginend = arg(5, stack, (KN< long > *)0);
  string *file = arg(6, stack, (string *)0);
  vector< R2_I2 > P;
  multimap< int, int > L;
  if (verbosity >= 1000) {
    debug = verbosity / 1000;
  } else {
    debug = 0;
  }

  map< pair< int, int >, int > FP;
  const double unset = -1e-100;
  KN< double > tff(nbv, unset);

  // loop over triangle
  if (pTh) {
    for (int it = 0; it < Th.nt; ++it) {
      for (int iv = 0; iv < 3; ++iv) {
        int i = Th(it, iv);
        if (tff[i] == unset) {
          mp->setP(pTh, it, iv);
          tff[i] = GetAny< double >((*eff)(stack)) - isovalue;
        }
      }
    }
  } else {
    ffassert(0);
  }

  if (close < 0) {
    tff = -tff;
    close = -close;
  }

  *mp = mps;    // restore the stat of local variable ...

  double tffmax = tff.max( ), tffmin = tff.min( );
  if (verbosity) {
    cout << " -- isoline close=" << close << " iso= " << isovalue << " " << epsr << endl
         << "    bound  isovalue :" << tffmin << " " << tffmax << endl;
  }

  double eps = (tffmax - tffmin) * epsr;
  if ((tffmax < 0.) || (tffmin > 0.)) {
    return 0L;
  }

  if (epsr < 0) {
    eps = -epsr;
  }

  ostream *fff = 0;
  if (debug > 9) {
    fff = new ofstream("g-iso");
  }

  for (int k = 0; k < Th.nt; ++k) {
    int iK[3] = {Th(k, 0), Th(k, 1), Th(k, 2)};
    R2 Pk[3] = {Th(iK[0]), Th(iK[1]), Th(iK[2])};
    R fk[3] = {tff[iK[0]], tff[iK[1]], tff[iK[2]]};
    R2 Qk[6];
    int i1[6], i2[6];
    int np = IsoLineK(Pk, fk, Qk, i1, i2, eps);
    if (np == 2) {
      for (int i = 0; i < np; ++i) {    // sort i1,i2 ..
        i1[i] = iK[i1[i]];
        i2[i] = iK[i2[i]];
        if (i2[i] < i1[i]) {
          std::swap(i1[i], i2[i]);
        }
      }

      int p[2];    // point number

      for (int i = 0; i < 2; ++i) {
        pair< map< pair< int, int >, int >::iterator, bool > ii;
        pair< int, int > e(i1[i], i2[i]);
        ii = FP.insert(make_pair(e, P.size( )));
        if (ii.second) {
          P.push_back(R2_I2(Qk[i]));
        }

        if (debug) {
          cout << i1[i] << " ---  " << i2[i] << " ;     " << ii.second << endl;
        }

        p[i] = ii.first->second;
      }

      // add line k[0], k[1]
      P[p[0]].add(p[0], p[1], L);
      if (debug) {
        cout << " +++ " << Qk[0] << " ->  " << Qk[1] << " :: " << p[0] << " -> " << p[1] << endl;
      }

      if (fff) {
        *fff << Qk[0] << "\n"
             << Qk[1] << "\n"
             << ((Qk[0] * 0.4 + Qk[1] * .6) + R2(Qk[0], Qk[1]).perp( ) * .4) << "\n\n";
      }
    }
  }

  if (fff) {
    delete fff;
  }

  if (close) {
    if (debug) {
      cout << " Close path " << endl;
    }

    for (int k = 0; k < Th.nt; ++k) {
      // Triangle &K=Th[k];
      for (int e = 0; e < 3; ++e) {
        int ee, kk = Th.ElementAdj(k, ee = e);
        if (kk == k || kk < 0) {    // true border element edge
          int iK[2] = {Th(k, (e + 1) % 3), Th(k, (e + 2) % 3)};
          R2 Pk[2] = {Th(iK[0]), Th(iK[1])};
          R fk[2] = {tff[iK[0]], tff[iK[1]]};
          R2 Qk[2];
          int i1[2], i2[2];
          if (debug) {
            cout << " LB : " << Pk[0] << ", " << fk[0] << " ->  " << Pk[1] << ", " << fk[1] << " : "
                 << iK[0] << " " << iK[1] << endl;
          }

          int np = LineBorder(Pk, fk, close, Qk, i1, i2, eps);
          if (np >= 10) {    // full edge
            int ke = 0;
          } else if (np == 2) {
            for (int i = 0; i < 2; ++i) {
              i1[i] = iK[i1[i]];
              i2[i] = iK[i2[i]];
              if (i2[i] < i1[i]) {
                std::swap(i1[i], i2[i]);
              }
            }

            if (debug) {
              cout << " add  : " << Qk[0] << ", " << i1[0] << ',' << i2[0] << " ->   " << Qk[1]
                   << ", " << i1[1] << ',' << i2[1] << endl;
            }

            int p[2];    // point number

            for (int i = 0; i < 2; ++i) {
              pair< map< pair< int, int >, int >::iterator, bool > ii;
              pair< int, int > ee(i1[i], i2[i]);
              ii = FP.insert(make_pair(ee, P.size( )));
              if (ii.second) {
                P.push_back(R2_I2(Qk[i]));
              }

              if (debug) {
                cout << i1[i] << " ---  " << i2[i] << " ;     " << ii.second << endl;
              }

              p[i] = ii.first->second;
            }

            // add line k[0], k[1]
            P[p[0]].add(p[0], p[1], L);
            if (debug) {
              cout << " +++ " << Qk[0] << " ->  " << Qk[1] << " :: " << p[0] << " -> " << p[1]
                   << endl;
            }
          }
        }
      }
    }

    if (debug) {
      cout << " End Close path " << endl;
    }
  }

  if (verbosity > 99) {    // dump the data base
    cout << " IsolineP1 " << endl;

    for (int i = 0; i < P.size( ); ++i) {
      cout << "\t" << i << " :  " << P[i].P << " ->  " << P[i].nx << endl;
    }

    cout << " multmap for execption " << endl;

    for (multimap< int, int >::const_iterator i = L.begin( ); i != L.end( ); ++i) {
      cout << "\t" << i->first << " -> " << i->second << endl;
    }

    cout << " End multmap for execption " << endl;
  }

  vector< int > iQ, QQ;    // QQ the list of curve
  int np = P.size( );
  KN< int > start(np);
  start = -1;
  int kkk = 0;

  for (int i = 0; i < np; ++i) {
    int nx = P[i].next(i, L, 0);
    if (nx >= 0) {
      start[nx] = i;
    }
  }

  vector< int > starting;
  int iss = 0;

  for (multimap< int, int >::const_iterator i = L.begin( ); i != L.end( ); ++i) {
    starting.push_back(i->first);
    start[i->second] = i->first;
  }

  for (int i = 0; i < np; ++i) {
    if (start[i] == -1) {
      starting.push_back(i);
    }
  }

  while (1) {
    ffassert(kkk++ < 100000);
    int kk = 0, k = 0;

    for (int i = 0; i < np; ++i) {    // correction FH 18/9/2012 ....
      kk += P[i].count(i, L);
      if (start[i] == -1) {
        k++;
      }
    }

    if (kk == 0) {
      break;
    }

    int i = -1;
    if (iss < starting.size( )) {
      start[i = starting[iss++]] = -1;
    } else {
      for (int ii = 0; ii < np; ++ii) {
        if ((start[ii] >= 0)) {
          i = ii;
        }
      }
    }

    if (i < 0) {
      break;
    }

    {
      if (verbosity > 9) {
        cout << "  isolineP1: start curve  = " << i << " -> " << P[i].next(i, L, 0);
      }

      if (verbosity > 99) {
        cout << " (" << i << endl;
      }

      iQ.push_back(QQ.size( ));
      QQ.push_back(i);
      start[i] = -2;
      int i0 = i, i1 = 0, ie = i;

      while (1) {
        i1 = P[i0].next(i0, L, 1);
        if (i1 < 0) {
          break;
        }

        if (start[i1] < 0) {
          if (verbosity > 99) {
            cout << " -- " << i1;
          }

          QQ.push_back(i1);
          break;
        }

        QQ.push_back(i1);
        if (verbosity > 99) {
          cout << " " << i1;
        }

        start[i1] = -2;
        i0 = i1;
      }

      if (verbosity > 99) {
        cout << ") " << endl;
      } else if (verbosity > 9) {
        cout << endl;
      }

      iQ.push_back(QQ.size( ));
    }
  }

  // sort iQ
  if (iQ.size( ) > 2) {
    vector< pair< int, pair< int, int > > > sQ(iQ.size( ) / 2);

    for (int i = 0, j = 0; i < iQ.size( ); i += 2, ++j) {
      int i0 = iQ[i];
      int i1 = iQ[i + 1];
      sQ[j] = make_pair(i0 - i1, make_pair(i0, i1));
    }

    std::sort(sQ.begin( ), sQ.end( ));

    for (int i = 0, j = 0; i < iQ.size( ); i += 2, ++j) {
      iQ[i] = sQ[j].second.first;
      iQ[i + 1] = sQ[j].second.second;
    }
  }

  if (smoothing > 0) {
    KN< R2 > P1(QQ.size( )), P2(QQ.size( ));

    for (int i = 0; i < QQ.size( ); ++i) {
      P1[i] = P[QQ[i]].P;
    }

    // Smoothing the curve
    double c1 = 1, c0 = 4, ct = 2 * c1 + c0;
    c1 /= ct;
    c0 /= ct;

    for (int i = 0; i < iQ.size( );) {
      int i0 = iQ[i++];
      int i1 = iQ[i++] - 1;
      int nbsmoothing = pow((i1 - i0), ratio) * smoothing;
      if (verbosity > 2) {
        cout << "     curve " << i << " size = " << i1 - i0 << " nbsmoothing = " << nbsmoothing
             << " " << i0 << " " << i1 << endl;
      }

      P2 = P1;

      for (int step = 0; step < nbsmoothing; ++step) {
        for (int j = i0 + 1; j < i1; ++j) {
          P2[j] = c0 * P1[j] + c1 * P1[j - 1] + c1 * P1[j + 1];
        }

        if (QQ[i0] == QQ[i1]) {    // close curve
          int j0 = i0 + 1;         // prec
          int j1 = i1 - 1;         // next
          P2[i0] = P2[i1] = c0 * P1[i0] + c1 * P1[j0] + c1 * P1[j1];
        }

        P1 = P2;
      }
    }

    for (int i = 0; i < QQ.size( ); ++i) {
      P[QQ[i]].P = P1[i];
    }
  }

  if (pbeginend) {
    pbeginend->resize(iQ.size( ));

    for (int i = 0; i < iQ.size( ); ++i) {
      (*pbeginend)[i] = iQ[i];
    }
  }

  if (pxx && pyy) {
    pxx->resize(QQ.size( ));
    pyy->resize(QQ.size( ));

    for (int i = 0; i < QQ.size( ); ++i) {
      int j = QQ[i];
      (*pxx)[i] = P[j].P.x;
      (*pyy)[i] = P[j].P.y;
    }
  } else if (pxy) {
    pxy->resize(3, QQ.size( ));

    for (int k = 0; k < iQ.size( ); k += 2) {
      int i0 = iQ[k], i1 = iQ[k + 1];
      double lg = 0;
      R2 Po = P[QQ[i0]].P;

      for (int i = i0; i < i1; ++i) {
        int j = QQ[i];
        (*pxy)(0, i) = P[j].P.x;
        (*pxy)(1, i) = P[j].P.y;
        lg += R2(P[j].P, Po).norme( );
        (*pxy)(2, i) = lg;
        Po = P[j].P;
      }
    }
  } else {
    ffassert(0);
  }

  nbc = iQ.size( ) / 2;
  if (file) {
    ofstream fqq(file->c_str( ));
    int i = 0, i0, i1, n2 = iQ.size( ), k = 0;

    while (i < n2) {
      k++;
      i0 = iQ[i++];
      i1 = iQ[i++];

      for (int l = i0; l < i1; ++l) {
        int j = QQ[l];
        fqq << P[j].P.x << " " << P[j].P.y << " " << k << " " << j << endl;
      }

      fqq << endl;
    }
  }

  /*
   * int err=0;
   * int pe[10];
   * for(int i=0;i< P.size();++i)
   * {
   * int pb = P[i].count(i,L);
   * if(pb)
   * if(err<10)
   * pe[err++]=i;
   * else
   * err++;
   * }
   *
   * if(err>0)
   * {
   * for(int i=0;i<10;i++)
   * cout << " PB point = " << pe[i] << " " << P[pe[i]].P << " odd count = " <<
   * P[pe[i]].count(pe[i],L)  << endl; ffassert(0);
   * }
   */
  // construction des courble

  return nbc;
}

class ISOLINE_P1 : public OneOperator {
 public:
  typedef const Mesh *pmesh;
  int cas;

  ISOLINE_P1( )
    : OneOperator(atype< long >( ), atype< pmesh >( ), atype< double >( ),
                  atype< KN< double > * >( ), atype< KN< double > * >( )),
      cas(4) {}

  ISOLINE_P1(int)
    : OneOperator(atype< long >( ), atype< pmesh >( ), atype< double >( ),
                  atype< KNM< double > * >( )),
      cas(3) {}

  E_F0 *code(const basicAC_F0 &args) const {
    if (cas == 4) {
      return new ISOLINE_P1_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]),
                               t[2]->CastTo(args[2]), t[3]->CastTo(args[3]));
    } else if (cas == 3) {
      return new ISOLINE_P1_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]),
                               t[2]->CastTo(args[2]));
    } else {
      ffassert(0);    // bug
    }
  }
};

class FINDLOCALMIN_P1 : public OneOperator {
 public:
  typedef const Mesh *pmesh;
  int cas;

  FINDLOCALMIN_P1( )
    : OneOperator(atype< KN< long > * >( ), atype< pmesh >( ), atype< KN< double > * >( ),
                  atype< KN< double > * >( )),
      cas(1) {}

  E_F0 *code(const basicAC_F0 &args) const {
    if (cas == 1) {
      return new FINDLOCALMIN_P1_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]),
                                    t[2]->CastTo(args[2]));
    } else {
      ffassert(0);    // bug
    }
  }
};

R3 *Curve(Stack stack, const KNM_< double > &b, const long &li0, const long &li1, const double &ss,
          long *const &pi) {
  assert(b.N( ) >= 3);
  int i0 = li0, i1 = li1, im;
  if (i0 < 0) {
    i0 = 0;
  }

  if (i1 < 0) {
    i1 = b.M( ) - 1;
  }

  double lg = b(2, i1);
  R3 Q;
  ffassert(lg > 0 && b(2, 0) == 0.);
  double s = ss * lg;
  int k = 0, k1 = i1;

  while (i0 < i1 - 1) {
    ffassert(k++ < k1);
    im = (i0 + i1) / 2;
    if (s < b(2, im)) {
      i1 = im;
    } else if (s > b(2, im)) {
      i0 = im;
    } else {
      Q = R3(b(0, im), b(1, im), 0);
      i0 = i1 = im;
      break;
    }
  }

  if (i0 < i1) {
    ffassert(b(2, i0) <= s);
    ffassert(b(2, i1) >= s);
    R2 A(b(0, i0), b(1, i0));
    R2 B(b(0, i1), b(1, i1));
    double l1 = (b(2, i1) - s);
    double l0 = s - b(2, i0);
    Q = (l1 * A + l0 * B) / (l1 + l0);
  }

  if (pi) {
    *pi = i0;
  }

  R3 *pQ = Add2StackOfPtr2Free(stack, new R3(Q));
  return pQ;
}

R3 *Curve(Stack stack, const KNM_< double > &b, const long &li0, const long &li1,
          const double &ss) {
  return Curve(stack, b, li0, li1, ss, 0);
}

double mesure(Stack stack, const KNM_< double > &b, const KN_< long > &be) {
  double mes = 0;
  int nbc2 = be.N( );

  for (int k = 0; k < nbc2;) {
    int i0 = be[k++];
    int i1 = be[k++];
    R2 A(b(0, i0), b(1, i0));
    double mk = 0;

    for (int i = i0 + 1; i < i1; ++i) {
      R2 B(b(0, i - 1), b(1, i - 1));
      R2 C(b(0, i), b(1, i));
      mk += det(A, B, C);
    }

    if (verbosity > 9) {
      cout << " mesure: composante " << k / 2 << "  mesure  " << mk / 2. << endl;
    }

    mes += mk;
  }

  return mes / 2.;
}

R3 *Curve(Stack stack, const KNM_< double > &b, const double &ss) {
  return Curve(stack, b, -1, -1, ss);
}

template< class R, class A0, class A1, class A2, class A3, class A4,
          class E = E_F0 >    // extend (4th arg.)
class E_F_F0F0F0F0F0s_ : public E {
 public:    // extend
  typedef R (*func)(Stack, const A0 &, const A1 &, const A2 &, const A3 &,
                    const A4 &);    // extend (4th arg.)
  func f;
  Expression a0, a1, a2, a3, a4;    // extend
  E_F_F0F0F0F0F0s_(func ff, Expression aa0, Expression aa1, Expression aa2, Expression aa3,
                   Expression aa4)                             // extend
    : f(ff), a0(aa0), a1(aa1), a2(aa2), a3(aa3), a4(aa4) {}    // extend (4th arg.)

  AnyType operator( )(Stack s) const {
    return SetAny< R >(f(s, GetAny< A0 >((*a0)(s)), GetAny< A1 >((*a1)(s)), GetAny< A2 >((*a2)(s)),
                         GetAny< A3 >((*a3)(s)), GetAny< A4 >((*a4)(s))));
  }    // extend (4th arg.)

  virtual size_t nbitem( ) const { return a4->nbitem( ); }    // modif

  bool MeshIndependent( ) const {
    return a0->MeshIndependent( ) && a1->MeshIndependent( ) && a2->MeshIndependent( ) &&
           a3->MeshIndependent( ) && a4->MeshIndependent( );
  }    // extend (4th arg.)
};

template< class R, class A = R, class B = A, class C = B, class D = C, class E = D,
          class CODE = E_F_F0F0F0F0F0s_< R, A, B, C, D, E, E_F0 > >    // extend (4th arg.)
class OneOperator5s_ : public OneOperator {
  aType r;    // return type
  typedef typename CODE::func func;
  func f;

 public:
  E_F0 *code(const basicAC_F0 &args) const {
    if (args.named_parameter && !args.named_parameter->empty( )) {
      CompileError(" They are used Named parameter ");
    }

    return new CODE(f, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]), t[2]->CastTo(args[2]),
                    t[3]->CastTo(args[3]), t[4]->CastTo(args[4]));
  }    // extend

  OneOperator5s_(func ff)
    :    // 3->4
      OneOperator(map_type[typeid(R).name( )], map_type[typeid(A).name( )],
                  map_type[typeid(B).name( )], map_type[typeid(C).name( )],
                  map_type[typeid(D).name( )],
                  map_type[typeid(E).name( )]),    // extens
      f(ff) {}
};

static void finit( ) {
  typedef const Mesh *pmesh;

  Global.Add("isoline", "(", new ISOLINE_P1);
  Global.Add("isoline", "(", new ISOLINE_P1(1));

  Global.Add("Curve", "(", new OneOperator2s_< R3 *, KNM_< double >, double >(Curve));
  Global.Add("Curve", "(", new OneOperator4s_< R3 *, KNM_< double >, long, long, double >(Curve));
  Global.Add("Curve", "(",
             new OneOperator5s_< R3 *, KNM_< double >, long, long, double, long * >(Curve));

  Global.Add("Area", "(", new OneOperator2s_< double, KNM_< double >, KN_< long > >(mesure));
  Global.Add("findalllocalmin", "(", new FINDLOCALMIN_P1);
}

LOADFUNC(finit);    // une variable globale qui serat construite  au chargement dynamique
