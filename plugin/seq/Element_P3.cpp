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
  // ------ P3  Hierarchical (just remove P1 node of the P2 finite element)  --------
  class TypeOfFE_P3Lagrange : public TypeOfFE {
   public:
    static const int k = 3;
    static const int ndf = (k + 2) * (k + 1) / 2;
    static int Data[];
    static double Pi_h_coef[];
    static const int nn[10][3];
    static const int aa[10][3];
    static const int ff[10];
    static const int il[10];
    static const int jl[10];
    static const int kl[10];

    TypeOfFE_P3Lagrange( ) : TypeOfFE(3 + 2 * 3 + 1, 1, Data, 4, 1, 16, 10, 0) {
      static const R2 Pt[10] = {R2(0 / 3., 0 / 3.), R2(3 / 3., 0 / 3.), R2(0 / 3., 3 / 3.),
                                R2(2 / 3., 1 / 3.), R2(1 / 3., 2 / 3.), R2(0 / 3., 2 / 3.),
                                R2(0 / 3., 1 / 3.), R2(1 / 3., 0 / 3.), R2(2 / 3., 0 / 3.),
                                R2(1 / 3., 1 / 3.)};
      // 3,4,5,6,7,8
      int other[10] = {-1, -1, -1, 4, 3, 6, 5, 8, 7, -1};
      int kk = 0;

      for (int i = 0; i < NbDoF; i++) {
        pij_alpha[kk++] = IPJ(i, i, 0);
        if (other[i] >= 0) {
          pij_alpha[kk++] = IPJ(i, other[i], 0);
        }

        P_Pi_h[i] = Pt[i];
      }

      assert(P_Pi_h.N( ) == NbDoF);
      assert(pij_alpha.N( ) == kk);
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
      for (int i = 0; i < 16; ++i) {
        v[i] = 1;
      }

      int e0 = K.EdgeOrientation(0);
      int e1 = K.EdgeOrientation(1);
      int e2 = K.EdgeOrientation(2);
      int ooo[6] = {e0, e0, e1, e1, e2, e2};
      int iii[6] = {};
      int jjj[6] = {};

      for (int i = 0; i < 6; ++i) {
        iii[i] = 3 + 2 * i;    // si  orient = 1
        jjj[i] = 4 + 2 * i;    // si orient = -1
      }

      for (int i = 0; i < 6; ++i) {
        if (ooo[i] == 1) {
          v[jjj[i]] = 0;
        } else {
          v[iii[i]] = 0;
        }
      }
    }
  };

  // on what     nu df on node node of df
  int TypeOfFE_P3Lagrange::Data[] = {
    0, 1, 2, 3, 3, 4, 4, 5, 5, 6,    // the support number  of the node of the df
    0, 0, 0, 0, 1, 0, 1, 0, 1, 0,    // the number of the df on  the node
    0, 1, 2, 3, 3, 4, 4, 5, 5, 6,    // the node of the df
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,    // the df come from which FE (generaly 0)
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9,    // which are de df on sub FE
    0,                               // for each compontant $j=0,N-1$ it give the sub FE associated
    0, 10};
  double TypeOfFE_P3Lagrange::Pi_h_coef[] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  void TypeOfFE_P3Lagrange::FB(const bool *whatd, const Mesh &, const Triangle &K,
                               const RdHat &PHat, RNMK_ &val) const {
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R L[3] = {l0 * k, l1 * k, l2 * k};

    throwassert(val.N( ) >= 10);
    throwassert(val.M( ) == 1);
    // Attention il faut renumeroter les fonction de bases
    // car dans freefem++, il y a un node par sommet, arete or element
    // et la numerotation naturelle  mais 2 noud pas arete
    // donc p est la perumation
    // echange de numerotation si les arete sont dans le mauvais sens
    int p[10] = {};

    for (int i = 0; i < 10; ++i) {
      p[i] = i;
    }

    if (K.EdgeOrientation(0) < 0) {
      Exchange(p[3], p[4]);    // 3,4
    }

    if (K.EdgeOrientation(1) < 0) {
      Exchange(p[5], p[6]);    // 5,6
    }

    if (K.EdgeOrientation(2) < 0) {
      Exchange(p[7], p[8]);    // 7,8
    }

    val = 0;

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
      R2 D[] = {K.H(0) * k, K.H(1) * k, K.H(2) * k};
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

#include "Element_P3.hpp"

  // Author: F. Hecht , P-H Tournier, Jet Hoe Tang jethoe.tang@googlemail.com
  // Jan 2017
  // in tets
  class TypeOfFE_P3_3d : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    typedef GFElement< Mesh3 > FElement;
    static const int kp = 3;    // P3
    static const int ndof = (kp + 3) * (kp + 2) * (kp + 1) / 6;
    static int dfon[];
    static int nl[20][3];
    static int cl[20][3];
    static int cp[20];
    static int pp[20][4];
    static const int d = Mesh::Rd::d;

    TypeOfFE_P3_3d( );    // constructor
    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
            RNMK_ &val) const;
    void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
             int *nump) const;
  };

  int TypeOfFE_P3_3d::nl[20][3] = {
    {0, 0, 0} /* 0 */,  {1, 1, 1} /* 1 */,  {2, 2, 2} /* 2 */,  {3, 3, 3} /* 3 */,
    {0, 0, 1} /* 4 */,  {0, 1, 1} /* 5 */,  {0, 0, 2} /* 6 */,  {0, 2, 2} /* 7 */,
    {0, 0, 3} /* 8 */,  {0, 3, 3} /* 9 */,  {1, 1, 2} /* 10 */, {1, 2, 2} /* 11 */,
    {1, 1, 3} /* 12 */, {1, 3, 3} /* 13 */, {2, 2, 3} /* 14 */, {2, 3, 3} /* 15 */,
    {1, 2, 3} /* 16 */, {0, 2, 3} /* 17 */, {0, 1, 3} /* 18 */, {0, 1, 2} /* 19 */};
  int TypeOfFE_P3_3d::cl[20][3] = {
    {0, 1, 2} /* 0 */,  {0, 1, 2} /* 1 */,  {0, 1, 2} /* 2 */,  {0, 1, 2} /* 3 */,
    {0, 1, 0} /* 4 */,  {0, 0, 1} /* 5 */,  {0, 1, 0} /* 6 */,  {0, 0, 1} /* 7 */,
    {0, 1, 0} /* 8 */,  {0, 0, 1} /* 9 */,  {0, 1, 0} /* 10 */, {0, 0, 1} /* 11 */,
    {0, 1, 0} /* 12 */, {0, 0, 1} /* 13 */, {0, 1, 0} /* 14 */, {0, 0, 1} /* 15 */,
    {0, 0, 0} /* 16 */, {0, 0, 0} /* 17 */, {0, 0, 0} /* 18 */, {0, 0, 0} /* 19 */};
  int TypeOfFE_P3_3d::cp[20] = {6, 6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1};
  int TypeOfFE_P3_3d::pp[20][4] = {
    {3, 0, 0, 0} /* 0 */,  {0, 3, 0, 0} /* 1 */,  {0, 0, 3, 0} /* 2 */,  {0, 0, 0, 3} /* 3 */,
    {2, 1, 0, 0} /* 4 */,  {1, 2, 0, 0} /* 5 */,  {2, 0, 1, 0} /* 6 */,  {1, 0, 2, 0} /* 7 */,
    {2, 0, 0, 1} /* 8 */,  {1, 0, 0, 2} /* 9 */,  {0, 2, 1, 0} /* 10 */, {0, 1, 2, 0} /* 11 */,
    {0, 2, 0, 1} /* 12 */, {0, 1, 0, 2} /* 13 */, {0, 0, 2, 1} /* 14 */, {0, 0, 1, 2} /* 15 */,
    {0, 1, 1, 1} /* 16 */, {1, 0, 1, 1} /* 17 */, {1, 1, 0, 1} /* 18 */, {1, 1, 1, 0} /* 19 */};
  int TypeOfFE_P3_3d::dfon[] = {1, 2, 1, 0};    // 2 dofs on each edge, 2 dofs on each face

  TypeOfFE_P3_3d::TypeOfFE_P3_3d( ) : GTypeOfFE< Mesh >(TypeOfFE_P3_3d::dfon, 1, 3, false, false) {
    typedef Element E;
    int n = this->NbDoF;
    bool dd = verbosity > 5;
    if (dd) {
      cout << "\n +++ P3  : ndof : " << n << " " << this->PtInterpolation.N( ) << endl;
    }

    R3 *Pt = this->PtInterpolation;
    // construction of interpolation ppoint

    {
      double cc = 1. / 3.;

      for (int i = 0; i < ndof; ++i) {
        Pt[i] = R3::KHat[0] * cc * pp[i][0] + R3::KHat[1] * cc * pp[i][1] +
                R3::KHat[2] * cc * pp[i][2] + R3::KHat[3] * cc * pp[i][3];
      }

      if (dd) {
        cout << this->PtInterpolation << endl;
      }
    }

    for (int i = 0; i < n; i++) {
      this->pInterpolation[i] = i;
      this->cInterpolation[i] = 0;
      this->dofInterpolation[i] = i;
      this->coefInterpolation[i] = 1.;
    }
  }

  void TypeOfFE_P3_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
                           int ocoef, int odf, int *nump) const {
    int n = this->NbDoF;
    int *p = M.p;

    for (int i = 0; i < n; ++i) {
      M.p[i] = i;
    }

    if (verbosity > 9) {
      cout << " P3  set:";
    }

    int dof = 4;

    for (int e = 0; e < 6; ++e) {
      int oe = K.EdgeOrientation(e);
      if (oe < 0) {
        swap(p[dof], p[dof + 1]);
      }

      dof += 2;
    }
  }

  void TypeOfFE_P3_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                          const RdHat &PHat, RNMK_ &val) const {
    assert(val.N( ) >= 20);    // 23 degrees of freedom
    assert(val.M( ) == 1);     // 3 components
    // int n = this->NbDoF;
    // -------------
    // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL
    // number (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
    // perm[3] is the local number of the vertex with the biggest global number.)
    // -------------
    R ld[4];
    PHat.toBary(ld);
    ld[0] *= 3.;
    ld[1] *= 3.;
    ld[2] *= 3.;
    ld[3] *= 3.;

    int p[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};

    {
      int dof = 4;

      for (int e = 0; e < 6; ++e) {
        int oe = K.EdgeOrientation(e);
        if (oe < 0) {
          swap(p[dof], p[dof + 1]);
        }

        dof += 2;
      }
    }
    static int ddd = 100;
    ddd++;
    val = 0.;
    RN_ f0(val('.', 0, op_id));
    if (ddd < 20) {
      cout << ld[0] << " " << ld[1] << " " << ld[2] << " " << ld[3] << " ::";
    }

    if (whatd & Fop_D0) {
      for (int i = 0; i < 20; ++i) {
        R fi = 1. / cp[i];

        for (int l = 0; l < 3; ++l) {
          fi *= ld[nl[i][l]] - cl[i][l];
        }

        if (ddd < 20) {
          cout << " " << fi;
        }

        f0[p[i]] = fi;
      }

      if (ddd < 20) {
        cout << endl;
      }
    }

    if (whatd & (Fop_D1 | Fop_D2)) {
      R3 Dld[4], Df[20];
      K.Gradlambda(Dld);
      Dld[0] *= 3.;
      Dld[1] *= 3.;
      Dld[2] *= 3.;
      Dld[3] *= 3.;

      for (int i = 0; i < 20; ++i) {
        R fi = 1. / cp[i];
        R3 &dfi = Df[p[i]];

        for (int l = 0; l < 3; ++l) {
          double ci = ld[nl[i][l]] - cl[i][l];
          dfi *= ci;
          dfi += fi * Dld[nl[i][l]];
          fi *= ci;
        }

        RN_ f0x(val('.', 0, op_dx));
        RN_ f0y(val('.', 0, op_dy));
        RN_ f0z(val('.', 0, op_dz));
        if (whatd & Fop_dx) {
          for (int i = 0; i < 20; ++i) {
            f0x[i] = Df[i].x;
          }
        }

        if (whatd & Fop_dy) {
          for (int i = 0; i < 20; ++i) {
            f0y[i] = Df[i].y;
          }
        }

        if (whatd & Fop_dz) {
          for (int i = 0; i < 20; ++i) {
            f0z[i] = Df[i].z;
          }
        }

        ffassert(!(whatd & Fop_D2));    // no D2 to do !!!
      }
    }
  }

  // link with FreeFem++
  static TypeOfFE_P3Lagrange P3LagrangeP3;
  // a static variable to add the finite element to freefem++
  static AddNewFE P3Lagrange("P3", &P3LagrangeP3);
  static TypeOfFE_P3_3d P3_3d;
  GTypeOfFE< Mesh3 > &Elm_P3_3d(P3_3d);

  static AddNewFE3 TFE_P3_3d("P33d", &Elm_P3_3d);
  static void init( ) {
    TEF2dto3d[&P3LagrangeP3] = &Elm_P3_3d;    // P3 -> P33d
  }
}    // namespace Fem2D
LOADFUNC(Fem2D::init);

// --- fin --
