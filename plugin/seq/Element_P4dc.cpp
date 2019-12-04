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
  // ------ P4dc  Hierarchical (just remove P1 node of the P2 finite element)  --------
  class TypeOfFE_P4dcLagrange : public TypeOfFE {
   public:
    static const int k = 4;
    static const int ndf = (k + 2) * (k + 1) / 2;
    static int Data[];
    static double Pi_h_coef[];
    static const int nn[15][4];
    static const int aa[15][4];
    static const int ff[15];
    static const int il[15];
    static const int jl[15];
    static const int kl[15];
    static const R2 G;
    static const R cshrink;
    static const R cshrink1;
    // (1 -1/3)*
    static R2 Shrink(const R2 &P) { return (P - G) * cshrink + G; }

    static R2 Shrink1(const R2 &P) { return (P - G) * cshrink1 + G; }

    TypeOfFE_P4dcLagrange( ) : TypeOfFE(3 + 3 * 3 + 3, 1, Data, 4, 1, 15, 15, Pi_h_coef) {
      static const R2 Pt[15] = {R2(0 / 4., 0 / 4.), R2(4 / 4., 0 / 4.), R2(0 / 4., 4 / 4.),
                                R2(3 / 4., 1 / 4.), R2(2 / 4., 2 / 4.), R2(1 / 4., 3 / 4.),
                                R2(0 / 4., 3 / 4.), R2(0 / 4., 2 / 4.), R2(0 / 4., 1 / 4.),
                                R2(1 / 4., 0 / 4.), R2(2 / 4., 0 / 4.), R2(3 / 4., 0 / 4.),
                                R2(1 / 4., 2 / 4.), R2(2 / 4., 1 / 4.), R2(1 / 4., 1 / 4.)};

      for (int i = 0; i < NbDoF; i++) {
        pij_alpha[i] = IPJ(i, i, 0);
        P_Pi_h[i] = Shrink(Pt[i]);
      }

      // 3,4,5, 6,7,8, 9,10,11,
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
  };

  const R2 TypeOfFE_P4dcLagrange::G(1. / 3., 1. / 3.);
  const R TypeOfFE_P4dcLagrange::cshrink = 1 - 1e-2;
  const R TypeOfFE_P4dcLagrange::cshrink1 = 1. / TypeOfFE_P4dcLagrange::cshrink;

  // on what     nu df on node node of df
  int TypeOfFE_P4dcLagrange::Data[] = {
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,  6,  6,  6,  6,    // the support number  of the node of the df
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,    // the number of the df on  the node
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,     // the node of the df
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,     // the df come from which FE (generaly 0)
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,    // which are de df on sub FE
    0,    // for each compontant $j=0,N-1$ it give the sub FE associated

    0, 15};
  double TypeOfFE_P4dcLagrange::Pi_h_coef[] = {1., 1., 1., 1., 1., 1., 1., 1.,
                                               1., 1., 1., 1., 1., 1., 1.};
  void TypeOfFE_P4dcLagrange::FB(const bool *whatd, const Mesh &, const Triangle &K,
                                 const RdHat &PHat, RNMK_ &val) const {
    R2 P = Shrink1(PHat);
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1. - P.x - P.y, l1 = P.x, l2 = P.y;
    R L[3] = {l0 * k, l1 * k, l2 * k};

    throwassert(val.N( ) >= 14);
    throwassert(val.M( ) == 1);
    // Attention il faut renumeroter les fonction de bases
    // car dans freefem++, il y a un node par sommet, arete or element
    // et la numerotation naturelle  mais 2 noud pas arete
    // donc p est la perumation
    // echange de numerotation si les arete sont dans le mauvais sens
    int p[15] = {};

    for (int i = 0; i < 15; ++i) {
      p[i] = i;
    }

    // if(K.EdgeOrientation(0) <0) Exchange(p[3],p[5]);// 3,4
    // if(K.EdgeOrientation(1) <0) Exchange(p[6],p[8]);// 5,6
    // if(K.EdgeOrientation(2) <0) Exchange(p[9],p[11]);// 7,8
    // cout << KN_<int>(p,10) <<endl;
    val = 0;
    /*
     * //  les fonction de base du Pk Lagrange sont
     * //
     * //
     */
    // --

    if (whatd[op_id]) {
      RN_ f0(val('.', 0, op_id));

      for (int df = 0; df < ndf; df++) {
        int pdf = p[df];
        R f = 1. / ff[df];

        for (int i = 0; i < k; ++i) {
          f *= L[nn[df][i]] - aa[df][i];
        }

        f0[pdf] = f;
      }
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      R ks = k * cshrink1;
      R2 D[] = {K.H(0) * ks, K.H(1) * ks, K.H(2) * ks};
      if (whatd[op_dx] || whatd[op_dy]) {
        for (int df = 0; df < ndf; df++) {
          int pdf = p[df];
          R fx = 0., fy = 0., f = 1. / ff[df];

          for (int i = 0; i < k; ++i) {
            int n = nn[df][i];
            R Ln = L[n] - aa[df][i];
            fx = fx * Ln + f * D[n].x;
            fy = fy * Ln + f * D[n].y;
            f = f * Ln;
          }

          if (whatd[op_dx]) {
            val(pdf, 0, op_dx) = fx;
          }

          if (whatd[op_dy]) {
            val(pdf, 0, op_dy) = fy;
          }
        }
      }

      if (whatd[op_dyy] || whatd[op_dxy] || whatd[op_dxx]) {
        for (int df = 0; df < ndf; df++) {
          int pdf = p[df];
          R fx = 0., fy = 0., f = 1. / ff[df];
          R fxx = 0., fyy = 0., fxy = 0.;

          for (int i = 0; i < k; ++i) {
            int n = nn[df][i];
            R Ln = L[n] - aa[df][i];
            fxx = fxx * Ln + 2. * fx * D[n].x;
            fyy = fyy * Ln + 2. * fy * D[n].y;
            fxy = fxy * Ln + fx * D[n].y + fy * D[n].x;
            fx = fx * Ln + f * D[n].x;
            fy = fy * Ln + f * D[n].y;
            f = f * Ln;
          }

          if (whatd[op_dxx]) {
            val(pdf, 0, op_dxx) = fxx;
          }

          if (whatd[op_dyy]) {
            val(pdf, 0, op_dyy) = fyy;
          }

          if (whatd[op_dxy]) {
            val(pdf, 0, op_dxy) = fxy;
          }
        }
      }
    }
  }

#include "Element_P4dc.hpp"

  // link with FreeFem++
  static TypeOfFE_P4dcLagrange P4dcLagrangeP4dc;
  // a static variable to add the finite element to freefem++
  static AddNewFE P4dcLagrange("P4dc", &P4dcLagrangeP4dc);
}    // namespace Fem2D

// --- fin --
