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
#include "AddNewFE.h"

// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend
// de l'orientation des aretes
//
/// ---------------------------------------------------------------
namespace Fem2D {
  struct InitTypeOfFE_PkEdge {
    int k;             // order poly on edge
    int npe;           // nb point on edge
    int ndf;           // nb dof
    KN< R > X;         // point on edge
    KN< int > Data;    // data of TypeOfFE
    InitTypeOfFE_PkEdge(int KK) : k(KK), npe(k + 1), ndf(3 * npe), X(npe), Data(5 * ndf + 3) {
      const QuadratureFormular1d QF(-1 + 2 * npe, npe, GaussLegendre(npe), true);

      for (int i = 0; i < npe; ++i) {
        X[i] = QF[i].x;
      }

      HeapSort((R *)X, npe);
      // int j = 0;
      int o[6];
      o[0] = 0;

      for (int i = 1; i < 6; ++i) {
        o[i] = o[i - 1] + ndf;
      }

      for (int df = 0; df < ndf; ++df) {
        int e = df / npe;
        int n = df % npe;
        Data[o[0] + df] = 3 + e;
        Data[o[1] + df] = n;
        Data[o[2] + df] = e;
        Data[o[3] + df] = 0;
        Data[o[4] + df] = df;
      }

      Data[o[5]] = 0;
      Data[o[5] + 1] = 0;
      Data[o[5] + 2] = ndf;    // end_dfcomp
    }
  };

  class TypeOfFE_PkEdge : public InitTypeOfFE_PkEdge, public TypeOfFE {
   public:
    static double Pi_h_coef[];

    TypeOfFE_PkEdge(int KK)
      : InitTypeOfFE_PkEdge(KK), TypeOfFE(ndf, 1, Data, -k, 1, ndf * 2, ndf, 0) {
      int kkk = 0;

      for (int i = 0; i < NbDoF; i++) {
        int e = i / npe;
        int j = i % npe;
        int ii = e * npe + npe - j - 1;
        R2 A(TriangleHat[VerticesOfTriangularEdge[e][0]]);
        R2 B(TriangleHat[VerticesOfTriangularEdge[e][1]]);
        pij_alpha[kkk++] = IPJ(i, i, 0);
        pij_alpha[kkk++] = IPJ(i, ii, 0);
        P_Pi_h[i] = A * (1. - X[j]) + B * (X[j]);    // X=0 => A  X=1 => B;
      }
    }

    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
      int kkk = 0;

      for (int e = 0; e < 3; ++e) {
        int i0 = 0;
        if (K.EdgeOrientation(e) < 0.) {
          i0 = 1 - i0;
        }

        int i1 = 1 - i0;

        for (int p = 0; p < npe; ++p) {
          v[kkk + i0] = 0;
          v[kkk + i1] = 1;
          kkk += 2;
        }
      }
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
  };

  // ENDOFCLASS TypeOfFE_PkEdge

  void TypeOfFE_PkEdge::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat,
                           RNMK_ &val) const {
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R L[3] = {l0, l1, l2};

    assert(val.N( ) >= ndf);
    assert(val.M( ) == 1);
    int ee = 0;
    if (L[0] <= min(L[1], L[2])) {
      ee = 0;    // arete
    } else if (L[1] <= min(L[0], L[2])) {
      ee = 1;
    } else {
      ee = 2;
    }

    int e3 = ee * npe;
    double s = 1. - L[ee];
    R xe = L[VerticesOfTriangularEdge[ee][0]] / s;    // go from 0 to 1 on edge
    R dxe = -1;
    if (K.EdgeOrientation(ee) < 0.) {
      xe = 1 - xe, dxe = -1;
    }

    assert(s);
    val = 0;
    if (whatd[op_id]) {
      RN_ f0(val('.', 0, op_id));

      for (int l = 0; l < npe; l++) {
        int dof = e3 + l;
        R f = 1.;

        for (int i = 0; i < npe; ++i) {
          if (i != l) {
            f *= (xe - X[i]) / (X[l] - X[i]);
          }
        }

        f0[dof] = f;
      }
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      R2 E = K.Edge(ee);
      R lE2 = (E, E);
      dxe /= lE2;

      for (int l = 0; l < npe; l++) {
        int dof = e3 + l;
        R f = 1., df = 0., ddf = 0.;

        for (int i = 0; i < npe; ++i) {
          if (i != l) {
            R xx = (xe - X[i]) / (X[l] - X[i]);
            R dxx = dxe / (X[l] - X[i]);
            ddf = ddf * xx + 2 * df * dxx;
            df = df * xx + f * dxx;
            f *= xx;
          }
        }

        if (whatd[op_id]) {
          val(dof, 0, op_id) = f;
        }

        if (whatd[op_dx]) {
          val(dof, 0, op_dx) = df * E.x;
        }

        if (whatd[op_dy]) {
          val(dof, 0, op_dy) = df * E.y;
        }

        if (whatd[op_dxx]) {
          val(dof, 0, op_dxx) = ddf * E.x * E.x;
        }

        if (whatd[op_dyy]) {
          val(dof, 0, op_dyy) = ddf * E.y * E.y;
        }

        if (whatd[op_dxy]) {
          val(dof, 0, op_dxy) = ddf * E.x * E.y;
        }
      }
    }
  }

  struct InitTypeOfFE_PkEdgedc {
    int k;             // order poly on edge
    int npe;           // nb point on edge
    int ndf;           // nb dof
    KN< R > X;         // point on edge
    KN< int > Data;    // data of TypeOfFE
    InitTypeOfFE_PkEdgedc(int KK) : k(KK), npe(k + 1), ndf(3 * npe), X(npe), Data(5 * ndf + 3) {
      const QuadratureFormular1d QF(-1 + 2 * npe, npe, GaussLegendre(npe), true);

      for (int i = 0; i < npe; ++i) {
        X[i] = QF[i].x;
      }

      HeapSort((R *)X, npe);
      int o[6];
      o[0] = 0;

      for (int i = 1; i < 6; ++i) {
        o[i] = o[i - 1] + ndf;
      }

      for (int df = 0; df < ndf; ++df) {
        Data[o[0] + df] = 6;    // 3+e;
        Data[o[1] + df] = df;
        Data[o[2] + df] = 0;
        Data[o[3] + df] = 0;
        Data[o[4] + df] = df;
      }

      Data[o[5]] = 0;
      Data[o[5] + 1] = 0;
      Data[o[5] + 2] = ndf;    // end_dfcomp
    }
  };

  class TypeOfFE_PkEdgedc : public InitTypeOfFE_PkEdgedc, public TypeOfFE {
   public:
    static double Pi_h_coef[];

    TypeOfFE_PkEdgedc(int KK)
      : InitTypeOfFE_PkEdgedc(KK), TypeOfFE(ndf, 1, Data, -2, 1, ndf * 2, ndf, 0) {
      int kkk = 0;

      for (int i = 0; i < NbDoF; i++) {
        int e = i / npe;
        int j = i % npe;
        int ii = e * npe + npe - j - 1;
        R2 A(TriangleHat[VerticesOfTriangularEdge[e][0]]);
        R2 B(TriangleHat[VerticesOfTriangularEdge[e][1]]);
        pij_alpha[kkk++] = IPJ(i, i, 0);
        pij_alpha[kkk++] = IPJ(i, ii, 0);
        P_Pi_h[i] = A * (1. - X[j]) + B * (X[j]);    // X=0 => A  X=1 => B;
      }
    }

    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
      int kkk = 0;

      for (int e = 0; e < 3; ++e) {
        int i0 = 0;
        if (K.EdgeOrientation(e) < 0.) {
          i0 = 1 - i0;
        }

        int i1 = 1 - i0;

        for (int p = 0; p < npe; ++p) {
          v[kkk + i0] = 0;
          v[kkk + i1] = 1;
          kkk += 2;
        }
      }
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
  };

  // ENDOFCLASS TypeOfFE_PkEdge

  void TypeOfFE_PkEdgedc::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat,
                             RNMK_ &val) const {
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R L[3] = {l0, l1, l2};

    assert(val.N( ) >= ndf);
    assert(val.M( ) == 1);
    int ee = 0;
    if (L[0] <= min(L[1], L[2])) {
      ee = 0;    // arete
    } else if (L[1] <= min(L[0], L[2])) {
      ee = 1;
    } else {
      ee = 2;
    }

    int e3 = ee * npe;
    double s = 1. - L[ee];
    R xe = L[VerticesOfTriangularEdge[ee][0]] / s;    // go from 0 to 1 on edge
    if (K.EdgeOrientation(ee) < 0.) {
      xe = 1 - xe;
    }

    assert(s);
    val = 0;
    if (whatd[op_id]) {
      RN_ f0(val('.', 0, op_id));

      for (int l = 0; l < npe; l++) {
        int df = e3 + l;
        R f = 1.;

        for (int i = 0; i < npe; ++i) {
          if (i != l) {
            f *= (xe - X[i]) / (X[l] - X[i]);
          }
        }

        f0[df] = f;
      }
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      cerr << " TO DO ???  FH " << endl;
      ffassert(0);
    }
  }

  // link with FreeFem++
  static TypeOfFE_PkEdgedc PkEdgeP0dc(0);
  static TypeOfFE_PkEdgedc PkEdgeP1dc(1);
  static TypeOfFE_PkEdgedc PkEdgeP2dc(2);
  static TypeOfFE_PkEdgedc PkEdgeP3dc(3);
  static TypeOfFE_PkEdgedc PkEdgeP4dc(4);
  static TypeOfFE_PkEdgedc PkEdgeP5dc(
    5);    // a static variable to add the finite element to freefem++
  static AddNewFE P0Edgedc("P0edgedc", &PkEdgeP0dc);
  static AddNewFE P1Edgedc("P1edgedc", &PkEdgeP1dc);
  static AddNewFE P2Edgedc("P2edgedc", &PkEdgeP2dc);
  static AddNewFE P3Edgedc("P3edgedc", &PkEdgeP3dc);
  static AddNewFE P4Edgedc("P4edgedc", &PkEdgeP4dc);
  static AddNewFE P5Edgedc("P5edgedc", &PkEdgeP5dc);

  // link with FreeFem++
  static TypeOfFE_PkEdge PkEdgeP1(1);
  static TypeOfFE_PkEdge PkEdgeP2(2);
  static TypeOfFE_PkEdge PkEdgeP3(3);
  static TypeOfFE_PkEdge PkEdgeP4(4);
  static TypeOfFE_PkEdge PkEdgeP5(5);
  // a static variable to add the finite element to freefem++
  static AddNewFE P1Edge("P1edge", &PkEdgeP1);
  static AddNewFE P2Edge("P2edge", &PkEdgeP2);
  static AddNewFE P3Edge("P3edge", &PkEdgeP3);
  static AddNewFE P4Edge("P4edge", &PkEdgeP4);
  static AddNewFE P5Edge("P5edge", &PkEdgeP5);
}    // namespace Fem2D

// --- fin --
