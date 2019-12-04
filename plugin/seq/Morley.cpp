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
// SUMMARY : Morley finite element
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

// the Polynomial space is P2, an the degree of freedom are
// the 3 values a the 3 vertices and the three
// normal derivative at the middle at the tree the edges
// remark;
// to compute the interpolante, we need
// the value , plus the value of the normal derivative
// so I use the following hack, I say the is a tree dim vectorial
// finite element with give  the value, x derivative ,and the y derivative
// Ref: chapter VII section 50  fig 50.2  of
// Ciarlet,   HandBook of Numerical Analysis, Volume II Finite elemet methodes (parts 1),
// NORTH-HOLLAND

// Related files:
// to check  and validate: testFEMorlay.edp
// to get a real example: bilapMorley.edp

#include "ff++.hpp"
#include "AddNewFE.h"

namespace Fem2D {
  // ------ P2 Morley
  class TypeOfFE_P2Morley : public TypeOfFE {
   public:
    static int Data[];

    TypeOfFE_P2Morley( )
      : TypeOfFE(3 + 3 + 0,
                 3,    // hack   u, u_x, u_y for interpolation
                 Data, 2, 1,
                 3 + 6,    // nb coef to build interpolation
                 6,        // np point to build interpolation
                 0) {
      const R2 Pt[] = {R2(0, 0), R2(1, 0), R2(0, 1), R2(0.5, 0.5), R2(0, 0.5), R2(0.5, 0)};
      // for the 3 vertices 6 coef
      int kk = 0;

      for (int p = 0; p < 3; p++) {
        P_Pi_h[p] = Pt[p];
        pij_alpha[kk] = IPJ(kk, p, 0);
        kk++;
      }

      // for
      int p = 3;

      for (int e = 0; e < 3; ++e) {    // point d'integration sur l'arete e
        P_Pi_h[p] = Pt[p];
        pij_alpha[kk++] = IPJ(3 + e, p, 1);    // coef = 0.5* l_e *ne_x * sge
        pij_alpha[kk++] = IPJ(3 + e, p, 2);    // coef = 0.5* l_e *ne_y * sge
        p++;
      }

      assert(P_Pi_h.N( ) == p);
      assert(pij_alpha.N( ) == kk);
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const R2 &P, RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const;
  };

  // on what     nu df on node node of df
  int TypeOfFE_P2Morley::Data[] = {0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0,
                                   0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 6, 6, 6};
  void TypeOfFE_P2Morley::Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
    const Triangle &T(K.T);
    int k = 0;

    // coef pour les 3 sommets  fois le 2 composantes
    for (int i = 0; i < 3; i++) {
      v[k++] = 1;
    }

    // integration sur les aretes
    for (int i = 0; i < 3; i++) {
      R2 N(T.Edge(i).perp( ));
      N *= T.EdgeOrientation(i);
      v[k++] = N.x;
      v[k++] = N.y;
    }
  }

  void TypeOfFE_P2Morley::FB(const bool *whatd, const Mesh &, const Triangle &K, const R2 &P,
                             RNMK_ &val) const {
    typedef double R;
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - P.x - P.y, l1 = P.x, l2 = P.y;
    R2 Dl[3] = {K.H(0), K.H(1), K.H(2)};
    R2 E[3] = {K.Edge(0), K.Edge(1), K.Edge(2)};
    double sgE[3] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
    R2 Ne[3] = {E[0].perp( ) * sgE[0], E[1].perp( ) * sgE[1], E[2].perp( ) * sgE[2]};
    double ccc[] = {1. / (Dl[0], Ne[0]), 1. / (Dl[1], Ne[1]), 1. / (Dl[2], Ne[2])};
    R l3 = (l0 - l0 * l0) * ccc[0];
    R l4 = (l1 - l1 * l1) * ccc[1];
    R l5 = (l2 - l2 * l2) * ccc[2];
    R dl3 = (1 - 2 * l0) * ccc[0];
    R dl4 = (1 - 2 * l1) * ccc[1];
    R dl5 = (1 - 2 * l2) * ccc[2];

    KN< bool > wd(KN_< const bool >(whatd, last_operatortype));
    KN< int > wd_op(last_operatortype), wd_j(last_operatortype);
    wd_j = 0;
    wd_op = -1;
    if (whatd[op_id]) {
      wd[op_dx] = wd[op_dy] = true;
      wd_op[op_dx] = wd_op[op_dy] = op_id;
      wd_j[op_dx] = 1;
      wd_j[op_dy] = 2;
    }

    if (whatd[op_dx]) {
      wd[op_dxx] = wd[op_dxy] = true;
      wd_op[op_dxx] = wd_op[op_dxy] = op_dx;
      wd_j[op_dxx] = 1;
      wd_j[op_dxy] = 2;
    }

    if (whatd[op_dy]) {
      wd[op_dyy] = wd[op_dxy] = true;
      wd_op[op_dxy] = wd_op[op_dyy] = op_dy;
      wd_j[op_dxy] = 1;
      wd_j[op_dyy] = 2;
    }

    // on previligie les originaux
    for (int i = 0; i < last_operatortype; ++i) {
      if (whatd[i]) {
        wd_op[i] = i;
        wd_j[i] = 0;
      }
    }

    val = 0;

    throwassert(val.N( ) >= 6);
    throwassert(val.M( ) == 3);

    val = 0;

    if (wd[op_id]) {
      RN_ f(val('.', 0, op_id));

      f[0] = l0;
      f[1] = l1;
      f[2] = l2;
      // remark  \int_O \Delta u = \int_G \dn(u)
      // \Delta l1^2 = Div ( 2 l1 \nalba l1) = 2 (\nalba l1,\nalba l1)
      f[3] = l3;    // (1-2 lO) \nabla l0
      f[4] = l4;
      f[5] = l5;
    }

    if (wd[op_dx]) {
      RN_ fx(val('.', wd_j[op_dx], wd_op[op_dx]));

      fx[0] = Dl[0].x;
      fx[1] = Dl[1].x;
      fx[2] = Dl[2].x;

      fx[3] = dl3 * Dl[0].x;
      fx[4] = dl4 * Dl[1].x;
      fx[5] = dl5 * Dl[2].x;
    }

    if (wd[op_dy]) {
      // RN_ fy(val('.',0,op_dy));
      RN_ fy(val('.', wd_j[op_dy], wd_op[op_dy]));
      fy[0] = Dl[0].y;
      fy[1] = Dl[1].y;
      fy[2] = Dl[2].y;

      fy[3] = dl3 * Dl[0].y;
      fy[4] = dl4 * Dl[1].y;
      fy[5] = dl5 * Dl[2].y;
    }

    if (wd[op_dxx]) {
      // RN_ fxx(val('.',0,op_dxx));
      RN_ fxx(val('.', wd_j[op_dxx], wd_op[op_dxx]));
      fxx[3] = Dl[0].x * Dl[0].x * ccc[0] * -2.;
      fxx[4] = Dl[1].x * Dl[1].x * ccc[1] * -2.;
      fxx[5] = Dl[2].x * Dl[2].x * ccc[2] * -2.;
    }

    if (wd[op_dyy]) {
      // RN_ fyy(val('.',0,op_dyy));
      RN_ fyy(val('.', wd_j[op_dyy], wd_op[op_dyy]));

      fyy[3] = Dl[0].y * Dl[0].y * ccc[0] * -2.;
      fyy[4] = Dl[1].y * Dl[1].y * ccc[1] * -2.;
      fyy[5] = Dl[2].y * Dl[2].y * ccc[2] * -2.;
    }

    if (wd[op_dxy]) {
      assert(val.K( ) > wd_op[op_dxy]);

      // RN_ fxy(val('.',0,op_dxy));
      RN_ fxy(val('.', wd_j[op_dxy], wd_op[op_dxy]));
      fxy[3] = Dl[0].x * Dl[0].y * ccc[0] * -2.;
      fxy[4] = Dl[1].x * Dl[1].y * ccc[1] * -2.;
      fxy[5] = Dl[2].x * Dl[2].y * ccc[2] * -2.;
    }

    {
      int vop[last_operatortype] = {}, nop = 0;

      for (int j = 0; j < last_operatortype; j++) {
        if (wd[j]) {
          vop[nop++] = j;
        }
      }

      // mise a zero des flux of int_e_j dn(w_i) pour i,j = 0,..,2

      for (int i = 0; i < 3; ++i) {
        for (int e = 0; e < 3; ++e) {
          // 0= int_e  w_i = dn(li) - a = a = (Dn[i],N)
          double a = (Dl[i], Ne[e]);

          // but   = 0
          for (int opp = 0; opp < nop; ++opp) {
            int op = vop[opp];
            int k = wd_op[op];
            int j = wd_j[op];
            val(i, j, k) -= a * val(e + 3, j, k);
          }
        }
      }

      ;
      // copie of the derivative for the hack.
      if (whatd[op_id]) {
        if (wd_op[op_dx] != op_id) {
          val('.', 1, op_id) = val('.', wd_j[op_dx], wd_op[op_dx]);
        }

        if (wd_op[op_dy] != op_id) {
          val('.', 2, op_id) = val('.', wd_j[op_dy], wd_op[op_dy]);
        }
      }

      if (whatd[op_dx]) {
        if (wd_op[op_dxx] != op_dx) {
          val('.', 1, op_dx) = val('.', wd_j[op_dxx], wd_op[op_dxx]);
        }

        if (wd_op[op_dxy] != op_dx) {
          val('.', 2, op_dx) = val('.', wd_j[op_dxy], wd_op[op_dxy]);
        }
      }

      if (whatd[op_dy]) {
        if (wd_op[op_dxy] != op_dy) {
          val('.', 1, op_dy) = val('.', wd_j[op_dxy], wd_op[op_dxy]);
        }

        if (wd_op[op_dyy] != op_dy) {
          val('.', 2, op_dy) = val('.', wd_j[op_dyy], wd_op[op_dyy]);
        }
      }
    }
  }

  // a static variable to add the finite element to freefem++
  static TypeOfFE_P2Morley P2LagrangeP2Morley;
  static AddNewFE P2Morley("P2Morley", &P2LagrangeP2Morley);
}    // namespace Fem2D
     // --- fin --
