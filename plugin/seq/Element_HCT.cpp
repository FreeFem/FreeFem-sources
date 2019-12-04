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
// SUMMARY : HCT finite element
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// Hanen Narje
// E-MAIL  : frederic.hecht@sorbonne-universite.fr
// ferchichihanen@gmail.com

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

// Reference
// Sur l'implementation des elements finis de Hsieh-Clough-Tocher complet et reduit
// https://hal.inria.fr/file/index/docid/76557/filename/RR-0004.pdf

// Related files:
// to check  and validate: testFEHCT.edp
// to get a real example: bilapHCT.edp

#include "ff++.hpp"
#include "AddNewFE.h"

namespace Fem2D {
  class TypeOfFE_HCT : public TypeOfFE {
   public:
    static int Data[];

    TypeOfFE_HCT( )
      : TypeOfFE(12,
                 3,    // hack   u, u_x, u_y for interpolation
                 Data, 2, 1,
                 9 + 6,    // nb coef to build interpolation
                 6,        // np point to build interpolation
                 0) {
      const R2 Pt[] = {R2(0, 0), R2(1, 0), R2(0, 1), R2(0.5, 0.5), R2(0, 0.5), R2(0.5, 0)};
      // for the 3 vertices 3 coef => 9 coef ..
      int kk = 0;

      for (int p = 0; p < 3; p++) {
        P_Pi_h[p] = Pt[p];
        pij_alpha[kk] = IPJ(kk, p, 0);
        kk++;    // VALUE
        pij_alpha[kk] = IPJ(kk, p, 1);
        kk++;    // DX
        pij_alpha[kk] = IPJ(kk, p, 2);
        kk++;    // DY
      }

      int p = 3;

      for (int e = 0; e < 3; ++e) {    // point d'integration sur l'arete e
        P_Pi_h[p] = Pt[p];
        pij_alpha[kk++] = IPJ(9 + e, p, 1);    // coef =  ne_x * sge
        pij_alpha[kk++] = IPJ(9 + e, p, 2);    // coef =  ne_y * sge
        p++;
      }

      assert(P_Pi_h.N( ) == p);
      assert(pij_alpha.N( ) == kk);
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const R2 &P, RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const;
  };
  // on what     nu df on node node of df
  int TypeOfFE_HCT::Data[] = {0, 0, 0, 1,  1,  1, 2, 2, 2, 3, 4,  5,     // support on what
                              0, 1, 2, 0,  1,  2, 0, 1, 2, 0, 0,  0,     // df on node
                              0, 0, 0, 1,  1,  1, 2, 2, 2, 3, 4,  5,     // th node
                              0, 0, 0, 0,  0,  0, 0, 0, 0, 0, 0,  0,     // df of prevoius FE
                              0, 1, 2, 3,  4,  5, 6, 7, 8, 9, 10, 11,    // ++ which df on prevoiu
                              0, 0, 0,                                   // ???
                              0, 0, 0, 12, 12, 12};
  void TypeOfFE_HCT::Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
    const Triangle &T(K.T);
    int k = 0;

    // coef pour les 3 sommets  fois le 3 composantes
    for (int i = 0; i < 9; i++) {
      v[k++] = 1;
    }

    // integration sur les aretes
    for (int i = 0; i < 3; i++) {
      R2 N(T.Edge(i).perp( ));
      N *= T.EdgeOrientation(i) / N.norme( );
      v[k++] = N.x;
      v[k++] = N.y;
    }

    ffassert(v.N( ) == k);
  }

  void set2zero(double *p, int k) {
    for (int i = 0; i < k; ++i) {
      p[i] = 0.;
    }
  }

#define P3(a, b, c) a *b *c
#define P3abcx(x, a, b, c) a##x *b *c + a *b##x *c + a *b *c##x
#define P3abcxy(x, y, a, b, c) \
  a##x *b##y *c + a##y *b##x *c + a##y *b *c##x + a##x *b *c##y + a *b##x *c##y + a *b##y *c##x
#define P3abcxyz(x, y, z, a, b, c)                                                               \
  a##x *b##y *c##z + a##y *b##x *c##z + a##y *b##z *c##x + a##x *b##z *c##y + a##z *b##x *c##y + \
    a##z *b##y *c##x

#define P3X(a, b, c) P3abcx(x, a, b, c)
#define P3Y(a, b, c) P3abcx(y, a, b, c)
#define P3XY(a, b, c) P3abcxy(x, y, a, b, c)
#define P3XX(a, b, c) P3abcxy(x, x, a, b, c)
#define P3YY(a, b, c) P3abcxy(y, y, a, b, c)
#define P3XXX(a, b, c) P3abcxyz(x, x, x, a, b, c)
#define P3YYY(a, b, c) P3abcxyz(y, y, y, a, b, c)
#define P3XXY(a, b, c) P3abcxyz(x, x, y, a, b, c)
#define P3XYY(a, b, c) P3abcxyz(x, y, y, a, b, c)

#define LL10(P3)                                                                                 \
  {                                                                                              \
    P3(li, li, li), P3(li1, li1, li1), P3(li2, li2, li2), P3(li, li, li2), P3(li, li, li1),      \
      P3(li1, li1, li), P3(li1, li1, li2), P3(li2, li2, li1), P3(li2, li2, li), P3(li, li1, li2) \
  }

  void TypeOfFE_HCT::FB(const bool *whatd, const Mesh &, const Triangle &K, const R2 &P,
                        RNMK_ &val) const {
    typedef double R;
    double area = K.area;
    int Nop = val.K( );
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l[3] = {1 - P.x - P.y, P.x, P.y};
    R2 Dl[3] = {K.H(0), K.H(1), K.H(2)};
    R2 E[3] = {K.Edge(0), K.Edge(1), K.Edge(2)};
    R lg2[3] = {E[0].norme2( ), E[1].norme2( ), E[2].norme2( )};
    R lg[3] = {sqrt(lg2[0]), sqrt(lg2[1]), sqrt(lg2[2])};
    R eta[3] = {(lg2[2] - lg2[1]) / lg2[0], (lg2[0] - lg2[2]) / lg2[1], (lg2[1] - lg2[0]) / lg2[2]};
    double sgE[3] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
    val = 0;

    throwassert(val.N( ) >= 6);
    throwassert(val.M( ) == 3);

    int i0 = 0;
    if (l[1] < l[i0]) {
      i0 = 1;
    }

    if (l[2] < l[i0]) {
      i0 = 2;
    }

    int i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;
    double etai = eta[i0], etai1 = eta[i1], etai2 = eta[i2];
    double li = l[i0], li1 = l[i1], li2 = l[i2];
    double lix = Dl[i0].x, li1x = Dl[i1].x, li2x = Dl[i2].x;
    double liy = Dl[i0].y, li1y = Dl[i1].y, li2y = Dl[i2].y;
    // i0,i1,i2,
    int p12[12] = {
      3 * i0,     3 * i1,     3 * i2, 3 * i0 + 2, 3 * i0 + 1, 3 * i1 + 2, 3 * i1 + 1,
      3 * i2 + 2, 3 * i2 + 1, 9 + i0, 9 + i1,     9 + i2};    // renumerotation DL .. ff-> paper

    /*
     * double ll[10]={ P3(li,li,li),   P3(li1,li1,li1),   P3(li2,li2,li2),
     * P3(li,li,li2),  P3(li,li,li1),   P3(li1,li1,li),
     * P3(li1,li1,li2), P3(li2,li2,li1),  P3(li2,li2,li),
     * P3(li,li1,li2) };
     */
    // paper DOF fig 1.5.1 corresponding array Ai:
    // ( P(q_ij), Dp(q_ij-1 - q_ij), Dp(q_ij+1 - q_ij), (j=0,1,2) (9 DOF)
    // Dp(b_ij)(h_ij)  ou bij = middle of edge ij (j=0,1,2) (3 DOF)
    // Warning
    // double ccc[] = { 1./(Dl[0],Ne[0]), 1./(Dl[1],Ne[1]), 1./(Dl[2],Ne[2]) };
    double c12 = 1. / 12.;
    double Ai[12][10] = {
      {(-0.5) * (etai1 - etai2), 0, 0, (1.5) * (3 + etai1), (1.5) * (3 - etai2), 0, 0, 0, 0, 0},
      {(0.5) * (1 - 2 * etai - etai2), 1, 0, (-1.5) * (1 - etai), (1.5) * (etai + etai2), 3, 3, 0,
       0, 3 * (1 - etai)},
      {(0.5) * (1 + 2 * etai + etai1), 0, 1, (-1.5) * (etai + etai1), (-1.5) * (1 + etai), 0, 0, 3,
       3, 3 * (1 + etai)},
      {(-c12) * (1 + etai1), 0, 0, (0.25) * (7 + etai1), (-0.5), 0, 0, 0, 0, 0},
      {(-c12) * (1 - etai2), 0, 0, (-0.5), (0.25) * (7 - etai2), 0, 0, 0, 0, 0},
      {(-c12) * (7 + etai2), 0, 0, (0.5), (0.25) * (5 + etai2), 1, 0, 0, 0, -1},
      {(1. / 6.) * (4 - etai), 0, 0, (-0.25) * (3 - etai), (-0.25) * (5 - etai), 0, 1, 0, 0,
       (0.5) * (3 - etai)},
      {(1. / 6.) * (4 + etai), 0, 0, (-0.25) * (5 + etai), (-0.25) * (3 + etai), 0, 0, 1, 0,
       (0.5) * (3 + etai)},
      {(-c12) * (7 - etai1), 0, 0, (0.25) * (5 - etai1), (0.5), 0, 0, 0, 1, -1},
      {(4. / 3.), 0, 0, -2, -2, 0, 0, 0, 0, 4},
      {(-2. / 3.), 0, 0, 2, 0, 0, 0, 0, 0, 0},
      {(-2. / 3.), 0, 0, 0, 2, 0, 0, 0, 0, 0}};
    const int nnzdd = 6 * 3;    // nb coef..
    double add[] = {1., 1., 0., 0., 1., 1., 1., 0., 0., 1., 1., 1., 0., 0., 1., 1., 1., 1.};
    int idd[] = {0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 7, 7, 8, 8, 9, 10, 11};
    int jdd[] = {0, 1, 2, 1, 2, 3, 4, 5, 4, 5, 6, 7, 8, 7, 8, 9, 10, 11};

    for (int i = 0, kb = 0, kc = 15; i < 3; ++i, kb += 5, ++kc) {
      int ii = i * 5 + 1;
      int ip = (i + 1) % 3;
      int is = (i + 2) % 3;
      R2 Es = E[is];
      R2 Ep = -E[ip];
      double dd[4] = {Es.x, Es.y, Ep.x, Ep.y};

      add[ii] = dd[0];
      add[ii + 1] = dd[2];
      add[ii + 2] = dd[1];
      add[ii + 3] = dd[3];

      add[kc] = (sgE[i] * 2 * area / lg[i]);    // hauteur
    }

    double AAA[12][10];
    set2zero(&AAA[0][0], 120);
    double AA[12][10];
    set2zero(&AA[0][0], 120);

    for (int jj = 0; jj < 10; ++jj) {
      for (int i = 0; i < 12; ++i) {
        AAA[p12[i]][jj] += Ai[i][jj];
      }
    }

    for (int k = 0; k < nnzdd; ++k) {
      int i = idd[k];
      int j = jdd[k];
      double dij = add[k];

      for (int jj = 0; jj < 10; ++jj) {
        AA[i][jj] += dij * AAA[j][jj];
      }
    }

    if (whatd[op_id] || whatd[op_dx] || whatd[op_dy]) {

      double ll[10] = LL10(P3);
      double llx[10] = LL10(P3X);
      double lly[10] = LL10(P3Y);
      RN_ f(val('.', 0, op_id));
      RN_ fx(val('.', 1, op_id));
      RN_ fy(val('.', 2, op_id));

      for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 10; ++j) {
          f[i] += AA[i][j] * ll[j];
          fx[i] += AA[i][j] * llx[j];
          fy[i] += AA[i][j] * lly[j];
        }
      }

      if (whatd[op_dx]) {
        val('.', 0, op_dx) = fx;
      }

      if (whatd[op_dy]) {
        val('.', 0, op_dy) = fy;
      }
    }

    if (whatd[op_dx] || whatd[op_dxx] || whatd[op_dxy]) {
      // cout << "dx dxx dxy"<< endl;
      double ll[10] = LL10(P3X);
      double llx[10] = LL10(P3XX);
      double lly[10] = LL10(P3XY);
      RN_ f(val('.', 0, op_dx));
      RN_ fx(val('.', 1, op_dx));
      RN_ fy(val('.', 2, op_dx));
      f = 0.;

      for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 10; ++j) {
          f[i] += AA[i][j] * ll[j];
          fx[i] += AA[i][j] * llx[j];
          fy[i] += AA[i][j] * lly[j];
        }
      }

      if (whatd[op_dxx]) {
        val('.', 0, op_dxx) = fx;
      }

      if (whatd[op_dxy]) {
        val('.', 0, op_dxy) = fy;
      }
    }

    if (whatd[op_dy] || whatd[op_dyy]) {
      double ll[10] = LL10(P3Y);
      double llx[10] = LL10(P3XY);
      double lly[10] = LL10(P3YY);
      RN_ f(val('.', 0, op_dy));
      RN_ fx(val('.', 1, op_dy));
      RN_ fy(val('.', 2, op_dy));
      f = 0.;

      for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 10; ++j) {
          f[i] += AA[i][j] * ll[j];
          fx[i] += AA[i][j] * llx[j];
          fy[i] += AA[i][j] * lly[j];
        }
      }

      if (whatd[op_dyy]) {
        val('.', 0, op_dyy) = fy;
      }
    }

    if (Nop > op_dxx) {
      val('.', 1, op_dxx) = NAN;
      val('.', 2, op_dxx) = NAN;
    }

    if (Nop > op_dyy) {
      val('.', 1, op_dyy) = NAN;
      val('.', 2, op_dyy) = NAN;
    }

    if (Nop > op_dxy) {
      val('.', 1, op_dxy) = NAN;
      val('.', 2, op_dxy) = NAN;
    }
  }

  // a static variable to add the finite element to freefem++
  static TypeOfFE_HCT Lagrange_HCT;
  static AddNewFE FE_HCT("HCT", &Lagrange_HCT);
}    // namespace Fem2D
// --- fin --
