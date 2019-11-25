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

#include "ff++.hpp"
#include "AddNewFE.h"

// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend
// de l'orientation des aretes
//
/// ---------------------------------------------------------------
namespace Fem2D {
  static const QuadratureFormular1d &QFE = QF_GaussLegendre3;
  static const GQuadratureFormular< R2 > &QFK = QuadratureFormular_T_5;

  class TypeOfFE_P2pnc : public TypeOfFE {
   public:
    static int Data[];

    TypeOfFE_P2pnc( ) : TypeOfFE(7, 1, Data, 3, 1, 2 * 3 * QFE.n + QFK.n, 3 * QFE.n + QFK.n, 0) {
      int p = 0, k = 0;
      const R2 Ph[] = {R2(0, 0), R2(1, 0), R2(0, 1)};

      for (int i = 0; i < 3; i++) {
        R2 A = Ph[(i + 1) % 3], B = Ph[(i + 2) % 3];

        for (int j = 0; j < QFE.n; j++) {
          R l0 = QFE[j].x, l1 = 1 - l0;
          pij_alpha[k++] = IPJ(2 * i, p, 0);
          pij_alpha[k++] = IPJ(2 * i + 1, p, 0);
          P_Pi_h[p++] = A * l0 + B * l1;
        }
      }

      for (int i = 0; i < QFK.n; i++) {
        pij_alpha[k++] = IPJ(6, p, 0);
        P_Pi_h[p++] = QFK[i];
      }

      ffassert(k == this->pij_alpha.N( ));
      ffassert(p == this->P_Pi_h.N( ));
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const;
  };

  // const QuadratureFormular1d & TypeOfFE_P2pnc::QFE=QF_GaussLegendre3;
  // const GQuadratureFormular<R2> &  TypeOfFE_P2pnc::QFK=QuadratureFormular_T_5;

  int TypeOfFE_P2pnc::Data[] = {3, 3, 4, 4, 5, 5, 6,    // the support number  of the node of the df
                                0, 1, 0, 1, 0, 1, 0,    // the number of the df on  the node
                                0, 0, 1, 1, 2, 2, 3,    // the node of the df
                                0, 0, 0, 0, 0, 0, 0,    // the df come from which FE (generaly 0)
                                0, 1, 2, 3, 4, 5, 6,    // which are de df on sub FE
                                0,    // for each compontant $j=0,N-1$ it give the sub FE associated
                                0, 7};
  void TypeOfFE_P2pnc::Pi_h_alpha(
    const baseFElement &K, KN_< double > &v) const {    // compute the coef of interpolation ...
    const Triangle &T(K.T);
    int k = 0;
    R oe[3] = {T.EdgeOrientation(0), T.EdgeOrientation(1), T.EdgeOrientation(2)};
    static int ddd = 10;

    ++ddd;

    for (int i = 0; i < 3; i++) {
      for (int p = 0; p < QFE.n; ++p) {
        R l1 = QFE[p].x, l0 = 1 - QFE[p].x;
        if (oe[i] < 0) {
          swap(l0, l1);
        }

        if (ddd < 3) {
          cout << p << " " << oe[i] << " " << l0 << " " << l1 << endl;
        }

        R p0 = l0, p1 = l1;
        R cc1 = p0 * QFE[p].a;
        R cc0 = p1 * QFE[p].a;
        v[k++] = cc0;
        v[k++] = cc1;
      }
    }

    for (int p = 0; p < QFK.n; ++p) {
      double w = QFK[p].a;
      R l1 = QFK[p].x, l2 = QFK[p].y, l0 = 1 - l1 - l2;
      R b = 1;
      v[k++] = b * w;
    }

    ffassert(k == this->pij_alpha.N( ));
  }

  void TypeOfFE_P2pnc::FB(const bool *whatd, const Mesh &TH, const Triangle &K, const RdHat &PHat,
                          RNMK_ &val) const {
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R oe[3] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
    R l7[] = {
      l0, l1, l2, l1 * l2, l0 * l2, l0 * l1, (l0 - l1) * (l1 - l2) * (l2 - l0)};    // 7 monome
    double C1[][7] = {{1, 1, -1, 3, 3, -1, -2} /* 0 */,    {3, -1, 1, 1, -1, 3, -2} /* 1 */,
                      {-1, 3, 3, -1, 1, 1, -2} /* 2 */,    {0, 0, -12, 0, 0, -12, 12} /* 3 */,
                      {0, -12, 0, -0, -12, 0, 12} /* 4 */, {-12, 0, -0, -12, -0, 0, 12} /* 5 */,
                      {-10, 10, -10, 10, -10, 10, 0} /* 6 */};

    if (val.N( ) < 7) {
      throwassert(val.N( ) >= 7);
    }

    throwassert(val.M( ) == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));
    int p[] = {0, 1, 2, 3, 4, 5, 6};
    if (oe[0] < 0) {
      swap(p[0], p[1]);
    }

    if (oe[1] < 0) {
      swap(p[2], p[3]);
    }

    if (oe[2] < 0) {
      swap(p[4], p[5]);
    }

    if (whatd[op_id]) {
      for (int i = 0; i < 7; ++i) {
        for (int j = 0; j < 7; ++j) {
          f0[p[i]] += C1[j][i] * l7[j];
        }
      }
    }

    if (whatd[op_dx] || whatd[op_dy]) {
      R2 D0 = K.H(0), D1 = K.H(1), D2 = K.H(2);
      R2 D7[] = {D0,
                 D1,
                 D2,
                 D1 * l2 + D2 * l1,
                 D0 * l2 + D2 * l0,
                 D0 * l1 + D1 * l0,
                 ((D0 - D1) * (l1 - l2) * (l2 - l0) + (l0 - l1) * (D1 - D2) * (l2 - l0) +
                  (l0 - l1) * (l1 - l2) * (D2 - D0))};

      if (whatd[op_dx]) {
        RN_ f0x(val('.', 0, op_dx));

        for (int i = 0; i < 7; ++i) {
          for (int j = 0; j < 7; ++j) {
            f0x[p[i]] += C1[j][i] * D7[j].x;
          }
        }
      }

      if (whatd[op_dy]) {
        RN_ f0y(val('.', 0, op_dy));

        for (int i = 0; i < 7; ++i) {
          for (int j = 0; j < 7; ++j) {
            f0y[p[i]] += C1[j][i] * D7[j].y;
          }
        }
      }
    }
  }
}    // namespace Fem2D
static TypeOfFE_P2pnc VTypeOfFE_P2pnc;
static AddNewFE P1dcLagrange("P2pnc", &VTypeOfFE_P2pnc);

static void finit( ) {
  // link with FreeFem++
}

LOADFUNC(finit)

// --- fin --
