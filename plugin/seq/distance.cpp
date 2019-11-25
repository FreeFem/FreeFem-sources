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
// SUMMARY : Signed distance function
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

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
#include <limits>
static int debug = 0, debug1 = 0;
using namespace Fem2D;

template< class Rd >
Rd zero(double a, const Rd &A, double b, const Rd &B) {    // zero affine
  double la = b / (b - a), lb = a / (a - b);

  return la * A + lb * B;
}

template< class Rd >
double distmin(const Rd &A, const Rd &B, const Rd &Q) {
  double d = 0;
  Rd AB(A, B), AQ(A, Q);
  double ab2 = (AB, AB);
  double lc = (AQ, AB) / ab2;
  Rd CQ = AQ - lc * AB;    // ( CQ , AB) = 0
  Rd C = A + lc * AB;      // or Q - CQ

  if (lc < 0) {
    d = Norme2(AQ);
  } else if (lc > 1.) {
    d = Norme2(Rd(B, Q));
  } else {
    d = Norme2(CQ);
  }

  if (verbosity > 9999) {
    cout << " distmin: d =" << d << " /" << lc << " :: " << A << " " << B << " " << Q << " C" << C
         << endl;
  }

  return d;
}

template< class Rd >
double distmin(const Rd &A, const Rd &B, const Rd &Q, double aq, double bq) {
  double d = 0;
  Rd AB(A, B), AQ(A, Q);
  double ab2 = (AB, AB);
  double lp = (AQ, AB) / ab2;
  Rd PQ = AQ - lp * AB;    // ( PQ , AB) = 0

  if (lp < 0) {
    d = aq;
  } else if (lp > 1.) {
    d = bq;
  } else {
    d = Norme2(PQ);
  }

  if (verbosity > 9999) {
    cout << " distmin:AB/Q: d =" << d << " /" << lp << " :: A " << A << " B " << B << " Q " << Q
         << "  P " << (A + lp * AB) << endl;
  }

  return d;
}

template< class Rd >
double distmin(const Rd &A, double a, const Rd &B, double b, const Rd &Q, double aq, double bq) {
  // compule the min_P in [A,B] of  p + || PQ ||
  // where p = affine interpolation of a,b on [A,B]
  // P = (1-l) A + l B ; p =  (1-l) a + l b  for l in [0,1]
  double dmin = min(a + aq, b + bq);
  double ab = b - a;
  Rd AB(A, B), AQ(A, Q);
  double ab2 = (AB, AB);
  Rd H = ab * AB / ab2;    // H = d p/ dP
  double h2 = (H, H);
  int cas = 0;

  if (h2 < 1) {
    cas = 1;
    double lc = (AQ, AB) / ab2;
    Rd CQ = AQ - lc * AB;    // ( CQ , AB) = 0
    Rd C = A + lc * AB;      // or Q - CQ
    assert(abs((CQ, AB)) < 1e-6);
    double r2 = (CQ, CQ) / ab2;
    double lpm = lc + copysign(sqrt(r2 * h2 / (1 - h2)), -ab);
    // lpm in [0,1]
    if (verbosity > 999) {
      Rd M = A + lpm * AB;
      cout << " lgm " << lpm << " r= " << sqrt(r2) << " M= " << M << " Q =" << Q
           << " ::" << a + lpm * ab << " " << ab << endl;
    }

    lpm = max(0., min(1., lpm));

    if (lpm > 0. && lpm < 1) {
      cas = 2;
      Rd M = A + lpm * AB, MQ(M, Q);
      dmin = a + lpm * ab + sqrt((MQ, MQ));
      /*  // check ???? min ....
       * lpm += 0.001;
       * MQ=Rd(A +lpm*AB,Q);
       * double dmin1 = a + lpm*ab + sqrt((MQ,MQ));
       * lpm -= 0.002;
       * MQ=Rd(A +lpm*AB,Q);
       * double dmin2 = a + lpm*ab + sqrt((MQ,MQ));
       * if(verbosity>99) cout << " ### "<< dmin1 << " " << dmin << " " << dmin2 << " " << endl;
       * ffassert(dmin1 >= dmin1&& dmin2 >= dmin );
       */
    }
  }

  if (verbosity > 99) {
    cout << " distmin/ AaBaQ " << A << " " << a << " / " << B << " " << b << " / " << Q
         << " / dmin= " << dmin << " cas =" << cas << endl;
  }

  return dmin;
}

double distmin(const R3 &A, double a, const R3 &B, double b, const R3 &C, double c, const R3 &Q,
               double aq, double bq, double cq) {
  // let fA the affine function on ABC
  const int in = 1;
  int cas = 0, flat = 0;
  double dmin = min(min(a + aq, b + bq), c + cq);
  R3 AB(A, B), AC(A, C), AQ(A, Q);
  double ab2 = (AB, AB), acab = (AC, AB), ac2 = (AC, AC);
  double aqab = (AQ, AB), aqac = (AQ, AC);
  double pdet = ab2 * ac2 - acab * acab;
  double pb = (aqab * ac2 - aqac * acab) / pdet;
  double pc = (ab2 * aqac - aqab * acab) / pdet;
  double pa = 1 - pb - pc;
  R3 P = pa * A + pb * B + pc * C;    // proj of Q on plan ABC.
  R3 PQ(P, Q);
  // PG = grad(fA)
  double fab = b - a, fac = c - a;

  if (abs(fab) + abs(fac) < 1e-16) {    // const function fA (flat)
    flat = 1;
    if (a >= 0. && b >= 0. && c >= 0.) {
      dmin = a + Norme2(PQ);
      cas = in;
    }
  } else {
    // Calcul de d FA / dP =
    R3 AZ = fab * AC - fac * AB;
    // fA in constant on line AZ
    R3 AG = AZ ^ PQ;    // otho of  AZ in ABC :
    // AG = gb AB + gc AC ;
    double agab = (AG, AB), agac = (AG, AC);
    double gab = (agab * ac2 - agac * acab) / pdet;
    double gac = (ab2 * agac - agab * acab) / pdet;
    R3 AGG = gab * AB + gac * AC;
    ffassert(Norme2(AGG - AG) < 1e-6);    // verif ..
    double fag = fab * gab + fac * gac;
    R3 H = AG / fag;
    double r2 = (PQ, PQ);
    double h2 = (H, H);
    double lpm = -sqrt(r2 * h2 / (1 - h2));
    double gb = gab / fag;
    double gc = gac / fag;
    double ga = -gb - gc;
    double am = pa + lpm * ga, bm = pb + lpm * gb, cm = pc + lpm * gc;
    if (am >= 0 && bm >= 0 && cm > 0.) {
      R3 M = am * A + bm * B + cm * C;
      R3 QM(Q, M);
      dmin = am * a + bm * b + cm * c + Norme2(QM);
      cas = in;
      if (debug1) {
        // verif .
        {
          double df1 = ga * a + gb * b + gc * c;
          cout << "df1 " << df1 << " " << fag << endl;

          ffassert(abs(df1 - 1) < 1e-6);
          ffassert(abs((H, PQ)) < 1e-6);
          ffassert(abs((AZ, PQ)) < 1e-6);
          ffassert(abs((AG, PQ)) < 1e-6);
        }
        {
          am += 0.001, bm -= 0.001;
          R3 M = am * A + bm * B + cm * C, QM(Q, M);
          double dmin1 = am * a + bm * b + cm * c + Norme2(QM);
          ffassert(dmin <= dmin1);
        }
        {
          bm += 0.001, cm -= 0.001;
          R3 M = am * A + bm * B + cm * C, QM(Q, M);
          double dmin2 = am * a + bm * b + cm * c + Norme2(QM);
          ffassert(dmin <= dmin2);
        }
      }
    }
  }

  if (cas != in) {
    if (flat) {    // externe => test 3 bord ...
      double dminab = a + distmin(A, B, Q, aq, bq);
      double dminac = a + distmin(A, C, Q, aq, cq);
      double dminbc = a + distmin(B, C, Q, bq, cq);
      dmin = min(min(dminab, dminac), min(dminbc, dmin));
    } else {    // externe => test 3 bord ...
      double dminab = distmin(A, a, B, b, Q, aq, bq);
      double dminac = distmin(A, a, C, c, Q, aq, cq);
      double dminbc = distmin(B, b, C, c, Q, bq, cq);
      dmin = min(min(dminab, dminac), min(dminbc, dmin));
    }
  }

  if (debug) {
    cout << "       AaBbCc/q  " << dmin << " " << cas << flat << endl;
  }

  return dmin;
}

double distmin(const R3 &A, double a, const R3 &B, double b, const R3 &C, double c, const R3 &Q) {
  R3 AQ(A, Q), BQ(B, Q), CQ(C, Q);
  double aq = Norme2(AQ), bq = Norme2(BQ), cq = Norme2(CQ);

  return distmin(A, a, B, b, C, c, Q, aq, bq, cq);
}

template< class Rd >
double distmin(const Rd &A, double a, const Rd &B, double b, const Rd &Q) {
  double aq = Norme2(Rd(A, Q));
  double bq = Norme2(Rd(B, Q));

  return distmin(A, a, B, b, Q, aq, bq);
}

double distmin(const R3 &A, const R3 &B, const R3 &C, const R3 &Q) {
  double d = 0;
  R3 AB(A, B), AC(A, C), AQ(A, Q);
  double ab2 = (AB, AB), acab = (AC, AB), ac2 = (AC, AC);
  double aqab = (AQ, AB), aqac = (AQ, AC);
  double det = ab2 * ac2 - acab * acab;
  double b = (aqab * ac2 - aqac * acab) / det;
  double c = (ab2 * aqac - aqab * acab) / det;
  double a = 1 - b - c;
  R3 P = a * A + b * B + c * C;    // proj of Q on plan ABC.

  if (debug) {
    cout << " distmin ABC/q " << a << " " << b << " " << c << endl;
  }

  R3 PQ(P, Q);
  assert(abs((PQ, AB)) < 1e-7);
  assert(abs((PQ, AC)) < 1e-7);
  if (a >= 0. && b >= 0. && c >= 0.) {
    d = Norme2(PQ);
  } else {
    double d1 = distmin(A, B, Q);
    double d2 = distmin(B, C, Q);
    double d3 = distmin(C, A, Q);
    d = min(min(d1, d2), d3);
  }

  return d;
}

typedef Mesh3::Element Tet;

int DistanceIso0(const Tet &K, double *f, double *fK) {
  int ret = 1;
  // fk value f
  const double eps = 1e-16;
  R3 P[10];
  int np = 0;

  if (abs(f[0]) < eps) {
    f[0] = 0.;
    P[np++] = K[0];
  }

  if (abs(f[1]) < eps) {
    f[1] = 0.;
    P[np++] = K[1];
  }

  if (abs(f[2]) < eps) {
    f[2] = 0.;
    P[np++] = K[2];
  }

  if (abs(f[3]) < eps) {
    f[3] = 0.;
    P[np++] = K[3];
  }

  for (int e = 0; e < 6; ++e) {
    int i1 = Tet::nvedge[e][0];
    int i2 = Tet::nvedge[e][1];

    if ((f[i1] < 0. && f[i2] > 0.) or (f[i1] > 0. && f[i2] < 0.)) {
      P[np++] = zero< R3 >(f[i1], K[i1], f[i2], K[i2]);
    }
  }

  if (np && debug) {
    cout << " np " << np << " " << P[0] << " " << P[1] << " :: " << f[0] << " " << f[1] << " "
         << f[2] << " " << f[3] << endl;
  }

  if (np == 0) {
    ret = 0;
  } else if (np == 1) {
    for (int j = 0; j < 4; ++j) {
      fK[j] = Norme2(R3(P[0], K[j]));
    }
  } else if (np == 2) {
    for (int j = 0; j < 4; ++j) {
      fK[j] = distmin(P[0], P[1], (R3)K[j]);
    }
  } else if (np == 3 || np == 4) {
    for (int j = 0; j < 4; ++j) {
      fK[j] = distmin(P[0], P[1], P[2], (R3)K[j]);
    }
  } else {
    fK[0] = fK[1] = fK[2] = fK[3] = 0.;    //
  }

  if (debug) {
    cout << ret << " 3d DistanceIso0  " << np << " " << fK[0] << " " << fK[1] << fK[2] << " "
         << fK[3] << endl;
  }

  return ret;
}

int DistanceIso0(const Triangle &K, double *f, double *fK) {
  // fk value f

  const double eps = 1e-16;
  R2 P[6];
  int np = 0;

  if (abs(f[0]) < eps) {
    f[0] = 0.;
  }

  if (abs(f[1]) < eps) {
    f[1] = 0.;
  }

  if (abs(f[2]) < eps) {
    f[2] = 0.;
  }

  int ke[6];

  for (int e = 0; e < 3; ++e) {
    int i1 = (e + 1) % 3;
    int i2 = (e + 2) % 3;
    if (f[i1] == 0) {
      ke[np] = i1, P[np++] = K[i1];
    } else if ((f[i1] < 0. && f[i2] > 0.) or (f[i1] > 0. && f[i2] < 0.)) {
      ke[np] = e;
      P[np++] = zero< R2 >(f[i1], K[i1], f[i2], K[i2]);
    }
  }

  if (np && debug) {
    cout << " np " << np << " " << P[0] << " " << P[1] << " :: " << f[0] << " " << f[1] << " "
         << f[2] << endl;
  }

  if (np == 0) {
    return 0;
  } else if (np == 1) {
    for (int j = 0; j < 3; ++j) {
      fK[j] = Norme2(R2(P[0], K[j]));
    }
  } else if (np == 2) {
    for (int j = 0; j < 3; ++j) {
      fK[j] = distmin(P[0], P[1], (R2)K[j]);
    }
  } else {
    fK[0] = fK[1] = fK[2] = 0.;    //
  }

  if (debug) {
    cout << np << " DistanceIso0  np="
         << " " << fK[0] << " " << fK[1] << " " << fK[2] << endl;
  }

  return np;
}

double CheckDist(double a, double b) {
  for (int i = 0; i < 30; ++i) {
    double a = 1, b = 1.1, c = 1.5;
    R3 A(-0.5, 0.001, 0.002), B(0.5, -0.001, 0.0001), C(0.0001, 1., -0.0003), Q(i * 0.1, 0.001, 1.);
    double d = distmin(A, a, B, b, C, c, Q);
    cout << " d = " << i << " == " << d << endl;
  }

  return 0;
}

// fonction determinant les points d'intersection

class Distance2d_Op : public E_F0mps {
 public:
  Expression eTh, eff, exx;
  static const int n_name_param = 1;    //
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
  Distance2d_Op(const basicAC_F0 &args, Expression tth, Expression fff, Expression xxx)
    : eTh(tth), eff(fff), exx(xxx) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type Distance2d_Op::name_param[] = {{"distmax", &typeid(double)}};

pair< double, long > Add(const Mesh &Th, int kk, int ee, double *fv) {
  const Triangle &K = Th[kk];
  int a = (ee + 1) % 3;
  int b = (ee + 2) % 3;
  int q = ee;
  int ia = Th(kk, a), ib = Th(kk, b), iq = Th(kk, q);
  double fq = distmin< R2 >(K[a], fv[ia], K[b], fv[ib], K[q]);

  if (debug) {
    cout << iq << " ** add " << kk << " " << ee << " ; " << fq << " :: " << fv[ia] << " " << fv[ib]
         << " || " << fv[iq] << endl;
  }

  return pair< double, long >(fq, kk * 3 + ee);
}

pair< double, long > Add(const Mesh3 &Th, int kk, int ee, double *fv) {
  typedef Mesh3::Element Tet;
  const Tet &K = Th[kk];
  int a = Tet::nvface[ee][0];
  int b = Tet::nvface[ee][1];
  int c = Tet::nvface[ee][2];
  int q = ee;
  int ia = Th(kk, a), ib = Th(kk, b), ic = Th(kk, c), iq = Th(kk, q);
  double fq = distmin(K[a], fv[ia], K[b], fv[ib], K[c], fv[ic], K[q]);
  if (debug) {
    cout << " ** add " << kk << " " << ee << " ; " << fq << " :: " << fv[ia] << " " << fv[ib] << " "
         << fv[ic] << " || " << fv[iq] << endl;
  }

  return pair< double, long >(fq, kk * 4 + ee);
}

int DistanceIso0(const Mesh &Th, int k, double *f, double *fv) {
  typedef Mesh::Element Elem;
  const int nbve = 3;
  const Elem &K = Th[k];
  int iK[nbve] = {Th(k, 0), Th(k, 1), Th(k, 2)};
  R fk[nbve] = {f[iK[0]], f[iK[1]], f[iK[2]]};
  // cut here ..
  double FK[nbve] = {fv[iK[0]], fv[iK[1]], fv[iK[2]]};
  int cas = DistanceIso0(K, fk, FK);
  if (cas > 1) {    // OK iso cut triangle
    //

    fv[iK[0]] = min(fv[iK[0]], FK[0]);
    fv[iK[1]] = min(fv[iK[1]], FK[1]);
    fv[iK[2]] = min(fv[iK[2]], FK[2]);
    if (debug) {
      cout << " DistanceIso0 set K" << cas << " " << iK[0] << " " << iK[1] << " " << iK[2] << " "
           << fv[iK[0]] << " " << fv[iK[1]] << " " << fv[iK[2]] << endl;
    }
  }

  return cas > 1;
}

int DistanceIso0(const Mesh3 &Th, int k, double *f, double *fv) {
  typedef Mesh3::Element Elem;
  const int nbve = 4;
  const Elem &K = Th[k];
  int iK[nbve] = {Th(k, 0), Th(k, 1), Th(k, 2), Th(k, 3)};
  R fk[nbve] = {f[iK[0]], f[iK[1]], f[iK[2]], f[iK[3]]};
  // cut here ..
  double FK[nbve] = {fv[iK[0]], fv[iK[1]], fv[iK[2]], fv[iK[3]]};
  int cas = DistanceIso0(K, fk, FK);
  if (cas > 0) {    // OK iso cut triangle
    fv[iK[0]] = min(fv[iK[0]], FK[0]);
    fv[iK[1]] = min(fv[iK[1]], FK[1]);
    fv[iK[2]] = min(fv[iK[2]], FK[2]);
    fv[iK[3]] = min(fv[iK[3]], FK[3]);
  }

  return cas;
}

template< class Mesh >
AnyType Distance(Stack stack, const Mesh *pTh, Expression eff, KN< double > *pxx, double dmax) {
  typedef typename Mesh::Element Elem;
  const int nbve = Mesh::Rd::d + 1;
  debug = 0;
  if (verbosity > 99) {
    debug = 1;
  }

  double unset = -std::numeric_limits< double >::max( );
  double distinf = std::numeric_limits< double >::max( );
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  double isovalue = 0.;
  ffassert(pTh);
  const Mesh &Th(*pTh);
  long nbv = Th.nv;    // nombre de sommet
  long nbt = Th.nt;    // nombre de triangles
  // long  nbe=Th.neb; // nombre d'aretes fontiere
  typedef KN< double > Rn;
  typedef KN< long > Zn;
  typedef pair< double, long > KEY;    // Distance max , vertex i/triangle k= 3*k+i
  if (verbosity > 2) {
    cout << "      distance max = " << dmax << " nt =" << nbt << endl;
  }

  Rn tK(nbt), tV(nbv), f(nbv);
  pxx->resize(nbv);
  Rn &fv = *pxx;
  Zn mK(nbt);
  mK = 0L;
  f = unset;
  fv = distinf;
  typedef std::priority_queue< KEY, std::vector< KEY >, std::greater< KEY > > myPQ;
  myPQ pqs;
  vector< long > markT(Th.nt, 1);

  for (int it = 0; it < Th.nt; ++it) {
    for (int iv = 0; iv < nbve; ++iv) {
      int i = Th(it, iv);
      if (f[i] == unset) {
        mp->setP(pTh, it, iv);
        f[i] = GetAny< double >((*eff)(stack)) - isovalue;
      }
    }
  }

  long err = 0, nt0 = 0;

  for (long k = 0; k < Th.nt; ++k) {
    int cas = DistanceIso0(Th, k, f, fv);
    if (cas) {
      ++nt0;
      markT[k] = 0;
    }    // init
  }

  for (long k = 0; k < Th.nt; ++k) {
    if (markT[k] == 0) {
      for (int e = 0; e < nbve; ++e) {
        int ee = e, kk = Th.ElementAdj(k, ee);
        if (kk >= 0 && (markT[kk] != 0)) {
          pqs.push(Add(Th, kk, ee, fv));
        }
      }
    }
  }

  //
  if (verbosity > 3) {
    cout << "    Distance: nb elemets in queue after init" << pqs.size( ) << " / nb set " << nt0
         << endl;
  }

  double fm;

  while (!pqs.empty( )) {
    KEY t = pqs.top( );
    pqs.pop( );
    fm = t.first;
    if (fm > dmax) {
      break;
    }

    int e = t.second % nbve;
    int k = t.second / nbve;
    int iq = Th(k, e);

    if (debug) {
      cout << iq << " " << fm << " -- k=" << k << " e=" << e << " fv:" << fv[iq] << " / "
           << markT[k] << endl;
    }

    if (markT[k] != 0) {    // already done, skeep
      markT[k] = 0;
      fv[iq] = min(fm, fv[iq]);

      for (int e = 0; e < nbve; ++e) {
        int ee = e, kk = Th.ElementAdj(k, ee);
        if (kk >= 0 && (markT[kk] != 0)) {
          pqs.push(Add(Th, kk, ee, fv));
        }
      }
    }
  }

  // clean
  err = 0;

  for (int i = 0; i < nbv; ++i) {
    if (fv[i] == distinf) {
      fv[i] = dmax;
      err++;
    }

    fv[i] = copysign(fv[i], f[i]);
  }

  if (err && fm < dmax) {
    if (verbosity) {
      cout << " Distance 2d : we miss some point " << err << endl;
    }
  }

  return err;
}

AnyType Distance2d_Op::operator( )(Stack stack) const {
  double distinf = std::numeric_limits< double >::max( );
  double dmax = arg(0, stack, distinf);

  KN< double > *pxx = 0;
  pxx = GetAny< KN< double > * >((*exx)(stack));
  const Mesh *pTh = GetAny< const Mesh * >((*eTh)(stack));

  return Distance< Mesh >(stack, pTh, eff, pxx, dmax);
}

class Distance3d_Op : public E_F0mps {
 public:
  Expression eTh, eff, exx;
  static const int n_name_param = 1;
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
  Distance3d_Op(const basicAC_F0 &args, Expression tth, Expression fff, Expression xxx)
    : eTh(tth), eff(fff), exx(xxx) {
    args.SetNameParam(n_name_param, name_param, nargs);
  }

  AnyType operator( )(Stack stack) const;
};

basicAC_F0::name_and_type Distance3d_Op::name_param[] = {{"distmax", &typeid(double)}};
inline double max(double a, double b, double c) { return max(max(a, b), c); }

inline double min(double a, double b, double c) { return min(min(a, b), c); }

AnyType Distance3d_Op::operator( )(Stack stack) const {
  double distinf = std::numeric_limits< double >::max( );
  double dmax = arg(0, stack, distinf);

  KN< double > *pxx = 0;
  pxx = GetAny< KN< double > * >((*exx)(stack));
  const Mesh3 *pTh = GetAny< const Mesh3 * >((*eTh)(stack));

  return Distance< Mesh3 >(stack, pTh, eff, pxx, dmax);
}

class Distance2d_P1 : public OneOperator {
 public:
  typedef const Mesh *pmesh;
  Distance2d_P1( )
    : OneOperator(atype< long >( ), atype< pmesh >( ), atype< double >( ),
                  atype< KN< double > * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new Distance2d_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]),
                             t[2]->CastTo(args[2]));
  }
};

class Distance3d_P1 : public OneOperator {
 public:
  typedef const Mesh3 *pmesh3;
  Distance3d_P1( )
    : OneOperator(atype< long >( ), atype< pmesh3 >( ), atype< double >( ),
                  atype< KN< double > * >( )) {}

  E_F0 *code(const basicAC_F0 &args) const {
    return new Distance3d_Op(args, t[0]->CastTo(args[0]), t[1]->CastTo(args[1]),
                             t[2]->CastTo(args[2]));
  }
};

static void finit( ) {
  typedef const Mesh *pmesh;
  typedef const Mesh3 *pmesh3;

  Global.Add("distance", "(", new Distance2d_P1);
  Global.Add("distance", "(", new Distance3d_P1);
  Global.Add("checkdist", "(", new OneOperator2< double, double, double >(CheckDist));
}

LOADFUNC(finit);    // une variable globale qui serat construite  au chargement dynamique
