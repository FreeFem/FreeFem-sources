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
// SUMMARY : P2 bubble finite element
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

/*
 * P2 bulle 3 bulle per face + on e full bull
 *
 */

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "ff++.hpp"
#include "AddNewFE.h"
#include <iostream>

namespace Fem2D {
  // Author: F. Hecht , P-H Tournier, J-H Tang
  // Jan 2017
  // in tets
  class TypeOfFE_P2_bulle3_3d : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    typedef GFElement< Mesh3 > FElement;

    static int dfon[];
    static const int d = Mesh::Rd::d;
    static const GQuadratureFormular< R1 > QFe;    // quadrature formula on an edge
    static const GQuadratureFormular< R2 > QFf;    // quadrature formula on a face
    TypeOfFE_P2_bulle3_3d( );                      // constructor
    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
            RNMK_ &val) const;
    void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
             int *nump) const;
  };

  int TypeOfFE_P2_bulle3_3d::dfon[] = {1, 1, 3, 1};    // 2 dofs on each edge, 2 dofs on each face

  TypeOfFE_P2_bulle3_3d::TypeOfFE_P2_bulle3_3d( )
    : GTypeOfFE< Mesh >(TypeOfFE_P2_bulle3_3d::dfon, 1, 3, false, false) {
    typedef Element E;
    const R alpha = 0.1885804846964451;    // ; (7.-sqrt(13.))/18.;
    const double a1 = alpha, a0 = (1 - 2 * a1);
    int n = this->NbDoF;
    const int d = E::Rd::d;
    if (verbosity > 9) {
      cout << "\n +++ P2 3bulle : ndof : " << n << endl;
    }

    R3 *Pt = this->PtInterpolation;
    // construction of interpolation ppoint

    {
      int k = 0;

      for (int i = 0; i <= d; ++i) {
        Pt[k++] = Rd( );
      }

      for (int i = 0; i < d; ++i) {
        Pt[i + 1][i] = 1.;
      }

      for (int i = 0; i < E::ne; ++i) {
        Pt[k++] = (Pt[E::nvedge[i][0]] + Pt[E::nvedge[i][1]]) * 0.5;
      }

      for (int i = 0; i < E::nf; ++i) {
        for (int j = 0; j < 3; ++j) {
          int j1 = (j + 1) % 3, j2 = (j + 2) % 3;
          Pt[k++] =
            Pt[E::nvface[i][j]] * a0 + Pt[E::nvface[i][j1]] * a1 + Pt[E::nvface[i][j2]] * a1;
        }
      }

      Pt[k++] = Rd(0.25, 0.25, 0.25);
      ffassert(k == n);
    }
    if (verbosity > 9) {
      cout << this->PtInterpolation << endl;
    }

    for (int i = 0; i < n; i++) {
      this->pInterpolation[i] = i;
      this->cInterpolation[i] = 0;
      this->dofInterpolation[i] = i;
      this->coefInterpolation[i] = 1.;
    }
  }

  void TypeOfFE_P2_bulle3_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
                                  int ocoef, int odf, int *nump) const {
    int n = this->NbDoF;
    int *p = M.p;

    for (int i = 0; i < n; ++i) {
      M.p[i] = i;
    }

    int k = 10;
    if (verbosity > 9) {
      cout << " P2 3 bulle set:";
    }

    for (int ff = 0; ff < Element::nf; ff++, k += 3) {
      // oriantation de la face  a endroit
      int fp = K.facePermutation(ff);
      if (fp & 1) {
        Exchange(p[k], p[k + 1]);
      }

      if (fp & 2) {
        Exchange(p[k + 1], p[k + 2]);
      }

      if (fp & 4) {
        Exchange(p[k], p[k + 1]);
      }
    }
  }

  void TypeOfFE_P2_bulle3_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                                 const RdHat &PHat, RNMK_ &val) const {
    assert(val.N( ) >= 10 + 3 * 4 + 1);    // 23 degrees of freedom
    assert(val.M( ) == 1);                 // 3 components
    // -------------
    // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL
    // number (i.e. perm[0] is the local number of the vertex with the smallest global number,
    // perm[3] is the local number of the vertex with the biggest global number.)
    const Element::Vertex *tV[4] = {&K.at(0), &K.at(1), &K.at(2), &K.at(3)};
    static const int nvf[4][3] = {{3, 2, 1}, {0, 2, 3}, {3, 1, 0}, {0, 1, 2}};
    static const int nve[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    // -------------
    // -------------
    // the 4 barycentric coordinates for the reference tetrahedron evaluated at the point P
    // (they have the same value at the real tetrahedron's point corresponding to the reference
    // tetrahedron's point P)

    // ===================================================================================================================================================
    //
    // permutation of vertices
    //
    // ===================================================================================================================================================
    int k = 10;
    int p[23] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22};

    for (int ff = 0; ff < Element::nf; ff++, k += 3) {
      // orientation de la face a envert
      int fp = K.facePermutation(ff);
      // fp=0; // No perm
      if (fp & 4) {
        Exchange(p[k], p[k + 1]);
      }

      if (fp & 2) {
        Exchange(p[k + 1], p[k + 2]);
      }

      if (fp & 1) {
        Exchange(p[k], p[k + 1]);
      }
    }

    assert(val.N( ) >= E::nv + E::ne);
    assert(val.M( ) == 1);

    // ===================================================================================================================================================
    //
    // some constants
    //
    // ===================================================================================================================================================
    R alpha = 0.1885804846964451;    // alpha = (7-sqrt(13))/18;
    R c = 1. / (alpha * alpha * (2.0 * alpha - 1.0) * (3.0 * alpha - 1.0));
    R beta1 = -0.1530178854880987;
    R beta2 = 0.1174552862797528;
    R beta3 = 0.1004882060140154;
    R beta4 = -0.1422503968333846;
    R beta5 = -0.4698211451190111;
    R beta6 = -0.0341147839626256;
    R beta7 = -0.0997720100233590;
    R w = 1. - PHat.sum( );
    // ===================================================================================================================================================
    //
    // basis function
    //
    // ===================================================================================================================================================
    RdHat P(PHat);
    R l[] = {w * (2 * w - 1),
             P.x * (2 * P.x - 1),
             P.y * (2 * P.y - 1),
             P.z * (2 * P.z - 1),
             4 * w * P.x,    //
             4 * w * P.y,
             4 * w * P.z,
             4 * P.x * P.y,
             4 * P.x * P.z,
             4 * P.y * P.z,

             (P.x * P.y * P.z * (P.z - alpha)) * c,    // 10 , v3
             (P.x * P.y * P.z * (P.y - alpha)) * c,    // v2
             (P.x * P.y * P.z * (P.x - alpha)) * c,    // v1

             (w * P.y * P.z * (w - alpha)) * c,      // 13,  v0
             (w * P.y * P.z * (P.y - alpha)) * c,    // 2
             (w * P.y * P.z * (P.z - alpha)) * c,    // 3

             (w * P.x * P.z * (P.z - alpha)) * c,    // 16,  v3
             (w * P.x * P.z * (P.x - alpha)) * c,    // v1
             (w * P.x * P.z * (w - alpha)) * c,      // v0

             (w * P.x * P.y * (w - alpha)) * c,      // 19,   0
             (w * P.x * P.y * (P.x - alpha)) * c,    // 1
             (w * P.x * P.y * (P.y - alpha)) * c,    // 2
             256 * w * P.x * P.y * P.z};
    // 1 PB numerotation des arete

    val = 0;
    RN_ f0(val('.', 0, op_id));
    if (whatd & Fop_D0) {
      // corner
      f0[p[0]] = l[0] + beta1 * (l[13] + l[18] + l[19]) +
                 beta2 * (l[14] + l[15] + l[16] + l[17] + l[20] + l[21]) + beta3 * l[22];
      f0[p[1]] = l[1] + beta1 * (l[12] + l[17] + l[20]) +
                 beta2 * (l[10] + l[11] + l[16] + l[18] + l[19] + l[21]) + beta3 * l[22];
      f0[p[2]] = l[2] + beta1 * (l[11] + l[14] + l[21]) +
                 beta2 * (l[10] + l[12] + l[13] + l[15] + l[19] + l[20]) + beta3 * l[22];
      f0[p[3]] = l[3] + beta1 * (l[10] + l[15] + l[16]) +
                 beta2 * (l[11] + l[12] + l[13] + l[14] + l[17] + l[18]) + beta3 * l[22];
      // edge
      f0[p[4]] = l[4] + beta4 * (l[16] + l[21]) + beta5 * (l[17] + l[18] + l[19] + l[20]) +
                 beta6 * l[22];    // 01 faces 2, 3
      f0[p[5]] = l[5] + beta4 * (l[15] + l[20]) + beta5 * (l[13] + l[14] + l[19] + l[21]) +
                 beta6 * l[22];    // 02 faces 1, 3
      f0[p[6]] = l[6] + beta4 * (l[14] + l[17]) + beta5 * (l[13] + l[15] + l[16] + l[18]) +
                 beta6 * l[22];    // 03 faces 1, 2
      f0[p[7]] = l[7] + beta4 * (l[10] + l[19]) + beta5 * (l[11] + l[12] + l[20] + l[21]) +
                 beta6 * l[22];    // 12 faces 0, 3
      f0[p[8]] = l[8] + beta4 * (l[11] + l[18]) + beta5 * (l[10] + l[12] + l[16] + l[17]) +
                 beta6 * l[22];    // 13 faces 0, 2
      f0[p[9]] = l[9] + beta4 * (l[12] + l[13]) + beta5 * (l[10] + l[11] + l[14] + l[15]) +
                 beta6 * l[22];    // 23 faces 0, 1

      // face
      for (int i = 10; i < 22; i++) {
        f0[p[i]] = l[i] + beta7 * l[22];
      }

      // center
      f0[p[22]] = l[22];
    }

    // ===================================================================================================================================================
    //
    // first derivative
    //
    // ===================================================================================================================================================
    if (whatd & (Fop_D1 | Fop_D2)) {
      R3 Dl[4];
      K.Gradlambda(Dl);
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // * * * * * * * * * * * * * * * * * * * * * * * * * * derivative of lambda
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // * * * * * * * * * * * * * * * * * * * * * * * * * *
      R lx[] = {
        4 * P.x + 4 * P.y + 4 * P.z - 3,
        4 * P.x - 1,
        0,
        0,
        4 - 4 * P.y - 4 * P.z - 8 * P.x,
        -4 * P.y,
        -4 * P.z,
        4 * P.y,
        4 * P.z,
        0,
        -c * P.y * P.z * (alpha - P.z),
        -c * P.y * P.z * (alpha - P.y),
        c * P.x * P.y * P.z - c * P.y * P.z * (alpha - P.x),
        c * P.y * P.z * (P.x + P.y + P.z - 1) + c * P.y * P.z * (alpha + P.x + P.y + P.z - 1),
        c * P.y * P.z * (alpha - P.y),
        c * P.y * P.z * (alpha - P.z),
        c * P.x * P.z * (alpha - P.z) + c * P.z * (alpha - P.z) * (P.x + P.y + P.z - 1),
        c * P.x * P.z * (alpha - P.x) - c * P.x * P.z * (P.x + P.y + P.z - 1) +
          c * P.z * (alpha - P.x) * (P.x + P.y + P.z - 1),
        c * P.z * (P.x + P.y + P.z - 1) * (alpha + P.x + P.y + P.z - 1) +
          c * P.x * P.z * (P.x + P.y + P.z - 1) + c * P.x * P.z * (alpha + P.x + P.y + P.z - 1),
        c * P.y * (P.x + P.y + P.z - 1) * (alpha + P.x + P.y + P.z - 1) +
          c * P.x * P.y * (P.x + P.y + P.z - 1) + c * P.x * P.y * (alpha + P.x + P.y + P.z - 1),
        c * P.x * P.y * (alpha - P.x) - c * P.x * P.y * (P.x + P.y + P.z - 1) +
          c * P.y * (alpha - P.x) * (P.x + P.y + P.z - 1),
        c * P.x * P.y * (alpha - P.y) + c * P.y * (alpha - P.y) * (P.x + P.y + P.z - 1),
        -P.y * P.z * (256 * P.x + 256 * P.y + 256 * P.z - 256) - 256 * P.x * P.y * P.z,
      };
      R ly[] = {
        4 * P.x + 4 * P.y + 4 * P.z - 3,
        0,
        4 * P.y - 1,
        0,
        -4 * P.x,
        4 - 8 * P.y - 4 * P.z - 4 * P.x,
        -4 * P.z,
        4 * P.x,
        0,
        4 * P.z,
        -c * P.x * P.z * (alpha - P.z),
        c * P.x * P.y * P.z - c * P.x * P.z * (alpha - P.y),
        -c * P.x * P.z * (alpha - P.x),
        c * P.z * (P.x + P.y + P.z - 1) * (alpha + P.x + P.y + P.z - 1) +
          c * P.y * P.z * (P.x + P.y + P.z - 1) + c * P.y * P.z * (alpha + P.x + P.y + P.z - 1),
        c * P.y * P.z * (alpha - P.y) - c * P.y * P.z * (P.x + P.y + P.z - 1) +
          c * P.z * (alpha - P.y) * (P.x + P.y + P.z - 1),
        c * P.y * P.z * (alpha - P.z) + c * P.z * (alpha - P.z) * (P.x + P.y + P.z - 1),
        c * P.x * P.z * (alpha - P.z),
        c * P.x * P.z * (alpha - P.x),
        c * P.x * P.z * (P.x + P.y + P.z - 1) + c * P.x * P.z * (alpha + P.x + P.y + P.z - 1),
        c * P.x * (P.x + P.y + P.z - 1) * (alpha + P.x + P.y + P.z - 1) +
          c * P.x * P.y * (P.x + P.y + P.z - 1) + c * P.x * P.y * (alpha + P.x + P.y + P.z - 1),
        c * P.x * P.y * (alpha - P.x) + c * P.x * (alpha - P.x) * (P.x + P.y + P.z - 1),
        c * P.x * P.y * (alpha - P.y) - c * P.x * P.y * (P.x + P.y + P.z - 1) +
          c * P.x * (alpha - P.y) * (P.x + P.y + P.z - 1),
        -P.x * P.z * (256 * P.x + 256 * P.y + 256 * P.z - 256) - 256 * P.x * P.y * P.z,
      };
      R lz[] = {
        4 * P.x + 4 * P.y + 4 * P.z - 3,
        0,
        0,
        4 * P.z - 1,
        -4 * P.x,
        -4 * P.y,
        4 - 4 * P.y - 8 * P.z - 4 * P.x,
        0,
        4 * P.x,
        4 * P.y,
        c * P.x * P.y * P.z - c * P.x * P.y * (alpha - P.z),
        -c * P.x * P.y * (alpha - P.y),
        -c * P.x * P.y * (alpha - P.x),
        c * P.y * (P.x + P.y + P.z - 1) * (alpha + P.x + P.y + P.z - 1) +
          c * P.y * P.z * (P.x + P.y + P.z - 1) + c * P.y * P.z * (alpha + P.x + P.y + P.z - 1),
        c * P.y * P.z * (alpha - P.y) + c * P.y * (alpha - P.y) * (P.x + P.y + P.z - 1),
        c * P.y * P.z * (alpha - P.z) - c * P.y * P.z * (P.x + P.y + P.z - 1) +
          c * P.y * (alpha - P.z) * (P.x + P.y + P.z - 1),
        c * P.x * P.z * (alpha - P.z) - c * P.x * P.z * (P.x + P.y + P.z - 1) +
          c * P.x * (alpha - P.z) * (P.x + P.y + P.z - 1),
        c * P.x * P.z * (alpha - P.x) + c * P.x * (alpha - P.x) * (P.x + P.y + P.z - 1),
        c * P.x * (P.x + P.y + P.z - 1) * (alpha + P.x + P.y + P.z - 1) +
          c * P.x * P.z * (P.x + P.y + P.z - 1) + c * P.x * P.z * (alpha + P.x + P.y + P.z - 1),
        c * P.x * P.y * (P.x + P.y + P.z - 1) + c * P.x * P.y * (alpha + P.x + P.y + P.z - 1),
        c * P.x * P.y * (alpha - P.x),
        c * P.x * P.y * (alpha - P.y),
        -P.x * P.y * (256 * P.x + 256 * P.y + 256 * P.z - 256) - 256 * P.x * P.y * P.z,
      };

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // * * * * * * * * * * * * * * * * * * * * * * * * * * derivative of reference basis function
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // * * * * * * * * * * * * * * * * * * * * * * * * * *
      R ff0x[23], ff0y[23], ff0z[23];
      // corner
      ff0x[p[0]] = lx[0] + beta1 * (lx[13] + lx[18] + lx[19]) +
                   beta2 * (lx[14] + lx[15] + lx[16] + lx[17] + lx[20] + lx[21]) + beta3 * lx[22];
      ff0y[p[0]] = ly[0] + beta1 * (ly[13] + ly[18] + ly[19]) +
                   beta2 * (ly[14] + ly[15] + ly[16] + ly[17] + ly[20] + ly[21]) + beta3 * ly[22];
      ff0z[p[0]] = lz[0] + beta1 * (lz[13] + lz[18] + lz[19]) +
                   beta2 * (lz[14] + lz[15] + lz[16] + lz[17] + lz[20] + lz[21]) + beta3 * lz[22];
      ff0x[p[1]] = lx[1] + beta1 * (lx[12] + lx[17] + lx[20]) +
                   beta2 * (lx[10] + lx[11] + lx[16] + lx[18] + lx[19] + lx[21]) + beta3 * lx[22];
      ff0y[p[1]] = ly[1] + beta1 * (ly[12] + ly[17] + ly[20]) +
                   beta2 * (ly[10] + ly[11] + ly[16] + ly[18] + ly[19] + ly[21]) + beta3 * ly[22];
      ff0z[p[1]] = lz[1] + beta1 * (lz[12] + lz[17] + lz[20]) +
                   beta2 * (lz[10] + lz[11] + lz[16] + lz[18] + lz[19] + lz[21]) + beta3 * lz[22];
      ff0x[p[2]] = lx[2] + beta1 * (lx[11] + lx[14] + lx[21]) +
                   beta2 * (lx[10] + lx[12] + lx[13] + lx[15] + lx[19] + lx[20]) + beta3 * lx[22];
      ff0y[p[2]] = ly[2] + beta1 * (ly[11] + ly[14] + ly[21]) +
                   beta2 * (ly[10] + ly[12] + ly[13] + ly[15] + ly[19] + ly[20]) + beta3 * ly[22];
      ff0z[p[2]] = lz[2] + beta1 * (lz[11] + lz[14] + lz[21]) +
                   beta2 * (lz[10] + lz[12] + lz[13] + lz[15] + lz[19] + lz[20]) + beta3 * lz[22];
      ff0x[p[3]] = lx[3] + beta1 * (lx[10] + lx[15] + lx[16]) +
                   beta2 * (lx[11] + lx[12] + lx[13] + lx[14] + lx[17] + lx[18]) + beta3 * lx[22];
      ff0y[p[3]] = ly[3] + beta1 * (ly[10] + ly[15] + ly[16]) +
                   beta2 * (ly[11] + ly[12] + ly[13] + ly[14] + ly[17] + ly[18]) + beta3 * ly[22];
      ff0z[p[3]] = lz[3] + beta1 * (lz[10] + lz[15] + lz[16]) +
                   beta2 * (lz[11] + lz[12] + lz[13] + lz[14] + lz[17] + lz[18]) + beta3 * lz[22];
      // edge
      ff0x[p[4]] = lx[4] + beta4 * (lx[16] + lx[21]) + beta5 * (lx[17] + lx[18] + lx[19] + lx[20]) +
                   beta6 * lx[22];
      ff0y[p[4]] = ly[4] + beta4 * (ly[16] + ly[21]) + beta5 * (ly[17] + ly[18] + ly[19] + ly[20]) +
                   beta6 * ly[22];
      ff0z[p[4]] = lz[4] + beta4 * (lz[16] + lz[21]) + beta5 * (lz[17] + lz[18] + lz[19] + lz[20]) +
                   beta6 * lz[22];
      ff0x[p[5]] = lx[5] + beta4 * (lx[15] + lx[20]) + beta5 * (lx[13] + lx[14] + lx[19] + lx[21]) +
                   beta6 * lx[22];
      ff0y[p[5]] = ly[5] + beta4 * (ly[15] + ly[20]) + beta5 * (ly[13] + ly[14] + ly[19] + ly[21]) +
                   beta6 * ly[22];
      ff0z[p[5]] = lz[5] + beta4 * (lz[15] + lz[20]) + beta5 * (lz[13] + lz[14] + lz[19] + lz[21]) +
                   beta6 * lz[22];
      ff0x[p[6]] = lx[6] + beta4 * (lx[14] + lx[17]) + beta5 * (lx[13] + lx[15] + lx[16] + lx[18]) +
                   beta6 * lx[22];
      ff0y[p[6]] = ly[6] + beta4 * (ly[14] + ly[17]) + beta5 * (ly[13] + ly[15] + ly[16] + ly[18]) +
                   beta6 * ly[22];
      ff0z[p[6]] = lz[6] + beta4 * (lz[14] + lz[17]) + beta5 * (lz[13] + lz[15] + lz[16] + lz[18]) +
                   beta6 * lz[22];
      ff0x[p[7]] = lx[7] + beta4 * (lx[10] + lx[19]) + beta5 * (lx[11] + lx[12] + lx[20] + lx[21]) +
                   beta6 * lx[22];
      ff0y[p[7]] = ly[7] + beta4 * (ly[10] + ly[19]) + beta5 * (ly[11] + ly[12] + ly[20] + ly[21]) +
                   beta6 * ly[22];
      ff0z[p[7]] = lz[7] + beta4 * (lz[10] + lz[19]) + beta5 * (lz[11] + lz[12] + lz[20] + lz[21]) +
                   beta6 * lz[22];
      ff0x[p[8]] = lx[8] + beta4 * (lx[11] + lx[18]) + beta5 * (lx[10] + lx[12] + lx[16] + lx[17]) +
                   beta6 * lx[22];
      ff0y[p[8]] = ly[8] + beta4 * (ly[11] + ly[18]) + beta5 * (ly[10] + ly[12] + ly[16] + ly[17]) +
                   beta6 * ly[22];
      ff0z[p[8]] = lz[8] + beta4 * (lz[11] + lz[18]) + beta5 * (lz[10] + lz[12] + lz[16] + lz[17]) +
                   beta6 * lz[22];
      ff0x[p[9]] = lx[9] + beta4 * (lx[12] + lx[13]) + beta5 * (lx[10] + lx[11] + lx[14] + lx[15]) +
                   beta6 * lx[22];
      ff0y[p[9]] = ly[9] + beta4 * (ly[12] + ly[13]) + beta5 * (ly[10] + ly[11] + ly[14] + ly[15]) +
                   beta6 * ly[22];
      ff0z[p[9]] = lz[9] + beta4 * (lz[12] + lz[13]) + beta5 * (lz[10] + lz[11] + lz[14] + lz[15]) +
                   beta6 * lz[22];

      // face
      for (int i = 10; i < 22; i++) {
        ff0x[p[i]] = lx[i] + beta7 * lx[22];
        ff0y[p[i]] = ly[i] + beta7 * ly[22];
        ff0z[p[i]] = lz[i] + beta7 * lz[22];
      }

      // center
      ff0x[p[22]] = lx[22];
      ff0y[p[22]] = ly[22];
      ff0z[p[22]] = lz[22];

      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // * * * * * * * * * * * * * * * * * * * * * * * * * * derivative of "exact" basis function
      // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      // * * * * * * * * * * * * * * * * * * * * * * * * * *
      RN_ f0x(val('.', 0, op_dx));
      RN_ f0y(val('.', 0, op_dy));
      RN_ f0z(val('.', 0, op_dz));

      for (int i = 0; i < 23; i++) {
        f0x[i] = Dl[1].x * ff0x[i] + Dl[2].x * ff0y[i] + Dl[3].x * ff0z[i];
        f0y[i] = Dl[1].y * ff0x[i] + Dl[2].y * ff0y[i] + Dl[3].y * ff0z[i];
        f0z[i] = Dl[1].z * ff0x[i] + Dl[2].z * ff0y[i] + Dl[3].z * ff0z[i];
      }
    }
  }

  static TypeOfFE_P2_bulle3_3d P2Bulle3_3d;
  GTypeOfFE< Mesh3 > &Elm_P2Bulle3_3d(P2Bulle3_3d);

  static AddNewFE3 TFE_P2Bulle3_3d("P2b3d", &Elm_P2Bulle3_3d);
}    // namespace Fem2D

// --- fin --
