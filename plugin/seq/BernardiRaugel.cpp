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
// SUMMARY : Bernardi Raugel finite element
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

// RELATED FILES:
// VALIDATION: test_P2BR.edp
// EXAMPLE: NS_P2BR_P0.edp

// *INDENT-OFF* //
// ff-c++-LIBRARY-dep:
// ff-c++-cpp-dep:
// *INDENT-ON* //

#include <ff++.hpp>
#include "AddNewFE.h"

namespace Fem2D {
  /*!
   * The P2BR finite element: the Bernardi Raugel Finite Element
   * See Bernardi, C., Raugel, G.: Analysis of some finite elements for the Stokes problem. Math.
   * Comp. 44, 71-79 (1985) It is a 2d coupled finite element The Polynomial space is P1^2 + 3
   * normals bubbles edges function (P_2) The degree of freedom is 6 values of the 2 components at
   * the 3 vertices and the 3 flux on the 3 edges So 9 degrees of freedom and N = 2
   */
  class TypeOfFE_P2BRLagrange : public TypeOfFE {
   public:
    static int Data[];

    TypeOfFE_P2BRLagrange( )
      : TypeOfFE(6 + 3 + 0,          // Number of degrees of freedom
                 2,                  // Dimension N
                 Data,               // Data array
                 4,                  // Number of subdivisions for plotting
                 1,                  // Number of sub finite element
                 6 + 3 * (2 + 2),    // Number kPi of coefficients to build the interpolation
                 9,                  // number nPi of integration points to build the interpolation
                 0    // Array to store the coefficient alpha_k to build the interpolator
        ) {
      const double gauss1 = (1. - sqrt(1. / 3.)) / 2;
      const double gauss2 = 1. - gauss1;
      const R2 Pt[] = {R2(0, 0), R2(1, 0), R2(0, 1)};

      // For the 3 vertices: 6 coefficents
      int kk = 0;

      for (int p = 0; p < 3; p++) {
        P_Pi_h[p] = Pt[p];
        pij_alpha[kk] = IPJ(kk, p, 0);
        ++kk;
        pij_alpha[kk] = IPJ(kk, p, 1);
        ++kk;
      }

      // Integration point on edge e
      int p = 3;

      for (int e = 0; e < 3; ++e) {
        R2 A = Pt[VerticesOfTriangularEdge[e][0]];
        R2 B = Pt[VerticesOfTriangularEdge[e][1]];
        P_Pi_h[p] = A * gauss1 + B * gauss2;
        pij_alpha[kk++] = IPJ(6 + e, p, 0);    // coef = 0.5 * l_e * ne_x * sge
        pij_alpha[kk++] = IPJ(6 + e, p, 1);    // coef = 0.5 * l_e * ne_y * sge
        p++;
        P_Pi_h[p] = A * gauss2 + B * gauss1;
        pij_alpha[kk++] = IPJ(6 + e, p, 0);    // coef = 0.5 * l_e * ne_x * sge
        pij_alpha[kk++] = IPJ(6 + e, p, 1);    // coef = 0.5 * l_e * ne_y * sge
        p++;
      }

      assert(P_Pi_h.N( ) == p);
      assert(pij_alpha.N( ) == kk);
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const;
  };

  // Data array
  int TypeOfFE_P2BRLagrange::Data[] = {0, 0, 1, 1, 2, 2, 3, 4, 5, 0, 1, 0, 1, 0, 1, 0, 0,
                                       0, 0, 0, 1, 1, 2, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0, 9, 9};

  /*!
   * \brief Define alpha k
   * \param K const baseFElement &
   * \param v KN<double> &
   */
  // Define alpha_k
  void TypeOfFE_P2BRLagrange::Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
    const Triangle &T(K.T);
    int k = 0;

    // Coefficents for the 3 vertices time the 2 components
    for (int i = 0; i < 6; i++) {
      v[k++] = 1;
    }

    // Integration on edges
    for (int i = 0; i < 3; i++) {
      R2 N(T.Edge(i).perp( ));
      N *= T.EdgeOrientation(i) * 0.5;
      v[k++] = N.x;
      v[k++] = N.y;
      v[k++] = N.x;
      v[k++] = N.y;
    }
  }

  /*!
   * \brief Shape function
   * \param whatd const bool*
   * \param const Mesh &
   * \param K const Triangle &
   * \param PHat const RdHat &
   * \param val RNMK_
   */
  // Shape function
  void TypeOfFE_P2BRLagrange::FB(const bool *whatd, const Mesh &, const Triangle &K,
                                 const RdHat &PHat, RNMK_ &val) const {
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    // int_e_1 l0*l0 = |e_1|/3 and int_e_1 l0*l1 = |e_1|/6
    // to get the flux = 1
    R2 E[3] = {K.Edge(0), K.Edge(1), K.Edge(2)};
    double l2E[3] = {(E[0], E[0]), (E[1], E[1]), (E[2], E[2])};
    // double lE[3] = {sqrt(l2E[0]), sqrt(l2E[1]), sqrt(l2E[2])};
    double sgE[3] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
    R2 cN[3] = {E[0].perp( ) * (6. * sgE[0] / l2E[0]), E[1].perp( ) * (6. * sgE[1] / l2E[1]),
                E[2].perp( ) * (6. * sgE[2] / l2E[2])};

    val = 0;

    throwassert(val.N( ) >= 9);
    throwassert(val.M( ) == 2);

    val = 0;

    if (whatd[op_id]) {
      RN_ f0(val('.', 0, op_id));
      RN_ f1(val('.', 1, op_id));

      f1[1] = f0[0] = l0;
      f1[3] = f0[2] = l1;
      f1[5] = f0[4] = l2;

      f0[6] = cN[0].x * l1 * l2;    // opposite to the vertex 0
      f0[7] = cN[1].x * l0 * l2;    // opposite to the vertex 1
      f0[8] = cN[2].x * l1 * l0;    // opposite to the vertex 2

      f1[6] = cN[0].y * l1 * l2;    // opposite to the vertex 0
      f1[7] = cN[1].y * l0 * l2;    // opposite to the vertex 1
      f1[8] = cN[2].y * l1 * l0;    // opposite to the vertex 2
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));

      if (whatd[op_dx]) {
        RN_ f0x(val('.', 0, op_dx));
        RN_ f1x(val('.', 1, op_dx));

        f1x[1] = f0x[0] = Dl0.x;
        f1x[3] = f0x[2] = Dl1.x;
        f1x[5] = f0x[4] = Dl2.x;

        f0x[6] = cN[0].x * (Dl1.x * l2 + Dl2.x * l1);
        f0x[7] = cN[1].x * (Dl2.x * l0 + Dl0.x * l2);
        f0x[8] = cN[2].x * (Dl0.x * l1 + Dl1.x * l0);

        f1x[6] = cN[0].y * (Dl1.x * l2 + Dl2.x * l1);
        f1x[7] = cN[1].y * (Dl2.x * l0 + Dl0.x * l2);
        f1x[8] = cN[2].y * (Dl0.x * l1 + Dl1.x * l0);
      }

      if (whatd[op_dy]) {
        RN_ f0y(val('.', 0, op_dy));
        RN_ f1y(val('.', 1, op_dy));

        f1y[1] = f0y[0] = Dl0.y;
        f1y[3] = f0y[2] = Dl1.y;
        f1y[5] = f0y[4] = Dl2.y;

        f0y[6] = cN[0].x * (Dl1.y * l2 + Dl2.y * l1);
        f0y[7] = cN[1].x * (Dl2.y * l0 + Dl0.y * l2);
        f0y[8] = cN[2].x * (Dl0.y * l1 + Dl1.y * l0);

        f1y[6] = cN[0].y * (Dl1.y * l2 + Dl2.y * l1);
        f1y[7] = cN[1].y * (Dl2.y * l0 + Dl0.y * l2);
        f1y[8] = cN[2].y * (Dl0.y * l1 + Dl1.y * l0);
      }

      if (whatd[op_dxx]) {
        RN_ f0xx(val('.', 0, op_dxx));
        RN_ f1xx(val('.', 1, op_dxx));

        f0xx[6] = 2 * cN[0].x * Dl1.x * Dl2.x;
        f0xx[7] = 2 * cN[1].x * Dl0.x * Dl2.x;
        f0xx[8] = 2 * cN[2].x * Dl0.x * Dl1.x;
        f1xx[6] = 2 * cN[0].y * Dl1.x * Dl2.x;
        f1xx[7] = 2 * cN[1].y * Dl0.x * Dl2.x;
        f1xx[8] = 2 * cN[2].y * Dl0.x * Dl1.x;
      }

      if (whatd[op_dyy]) {
        RN_ f0yy(val('.', 0, op_dyy));
        RN_ f1yy(val('.', 1, op_dyy));

        f0yy[6] = 2 * cN[0].x * Dl1.y * Dl2.y;
        f0yy[7] = 2 * cN[1].x * Dl0.y * Dl2.y;
        f0yy[8] = 2 * cN[2].x * Dl0.y * Dl1.y;
        f1yy[6] = 2 * cN[0].y * Dl1.y * Dl2.y;
        f1yy[7] = 2 * cN[1].y * Dl0.y * Dl2.y;
        f1yy[8] = 2 * cN[2].y * Dl0.y * Dl1.y;
      }

      if (whatd[op_dxy]) {
        assert(val.K( ) > op_dxy);
        RN_ f0xy(val('.', 0, op_dxy));
        RN_ f1xy(val('.', 1, op_dxy));

        f0xy[6] = cN[0].x * (Dl1.x * Dl2.y + Dl1.y * Dl2.x);
        f0xy[7] = cN[1].x * (Dl0.x * Dl2.y + Dl0.y * Dl2.x);
        f0xy[8] = cN[2].x * (Dl0.x * Dl1.y + Dl0.y * Dl1.x);
        f1xy[6] = cN[0].y * (Dl1.x * Dl2.y + Dl1.y * Dl2.x);
        f1xy[7] = cN[1].y * (Dl0.x * Dl2.y + Dl0.y * Dl2.x);
        f1xy[8] = cN[2].y * (Dl0.x * Dl1.y + Dl0.y * Dl1.x);
      }
    }

    // Now, remove the flux part on 6 first DOF
    // w_i = w_i - a_i w_{k_i} - b_i w_{l_i}
    {
      int k[6] = {6 + 1, 6 + 1, 6 + 2, 6 + 2, 6 + 0, 6 + 0};
      int l[6] = {6 + 2, 6 + 2, 6 + 0, 6 + 0, 6 + 1, 6 + 1};
      R2 eN[3] = {E[0].perp( ) * (0.5 * sgE[0]), E[1].perp( ) * (0.5 * sgE[1]),
                  E[2].perp( ) * (0.5 * sgE[2])};
      double a[6] = {eN[1].x, eN[1].y, eN[2].x, eN[2].y, eN[0].x, eN[0].y};
      double b[6] = {eN[2].x, eN[2].y, eN[0].x, eN[0].y, eN[1].x, eN[1].y};
      int nop = 0;
      int vop[last_operatortype] = {};

      for (int j = 0; j < last_operatortype; j++) {
        if (whatd[j]) {
          vop[nop++] = j;
        }
      }

      for (int i = 0; i < 6; ++i) {
        for (int jj = 0; jj < nop; ++jj) {
          int j = vop[jj];
          val(i, 0, j) -= a[i] * val(k[i], 0, j) + b[i] * val(l[i], 0, j);
          val(i, 1, j) -= a[i] * val(k[i], 1, j) + b[i] * val(l[i], 1, j);
        }
      }
    }
  }

  // Add the finite element to the FreeFem++ table
  // a static variable to define the finite element
  static TypeOfFE_P2BRLagrange P2LagrangeP2BR;
  // now, adding FE in FreeFem++ table
  static AddNewFE P2BR("P2BR", &P2LagrangeP2BR);
}    // namespace Fem2D
