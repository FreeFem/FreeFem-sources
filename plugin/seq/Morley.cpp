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
// normal derivative at the middle at the three the edges
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

//  Add version 3D on tet
// Dof  mean value on edge
// value on normal derivative on barycenter of face

#include "ff++.hpp"
#include "AddNewFE.h"

namespace Fem2D {
  // ------ P2 Morley 2d
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

//Morley 3d

class TypeOfFE_Morley_3d : public GTypeOfFE< Mesh3 > {
 public:
  typedef Mesh3 Mesh;
  typedef Mesh3::Element Element;
  typedef GFElement< Mesh3 > FElement;

  static int dfon[];
  static const int d = Mesh::Rd::d;
  static const GQuadratureFormular< R1 > QFe;    // quadrature formula on an edge
  static const GQuadratureFormular< R2 > QFf;    // quadrature formula on a face
  TypeOfFE_Morley_3d( );                          // constructor
  void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
          RNMK_ &val) const;
  void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
           int *nump) const;
};

int TypeOfFE_Morley_3d::dfon[] = {0, 1, 1, 0};    // 1 dofs on each edge, 1 dofs on each face

// Quadrature formula on an edge, exact for degree 2 (ok: int_e (deg2*t *lambda))
const GQuadratureFormular< R1 > TypeOfFE_Morley_3d::QFe(-1 + 2 * 2, 2, GaussLegendre(2), true);
// arguments: exact, num of integration pts, integration pts, clean (~GQuadratureFormular()
// {if(clean) delete [] p;}) GaussLegendre defined in QuadratureFormular.cpp

// Quadrature formula on a face, exact for degree 2 (ok: int_f (deg2*t)), internal quadrature
// points
static GQuadraturePoint< R2 > P_QuadratureFormular_T_2_intp[3] = {
  GQuadraturePoint< R2 >(1. / 3., R2(1. / 6., 4. / 6.)),
  GQuadraturePoint< R2 >(1. / 3., R2(4. / 6., 1. / 6.)),
  GQuadraturePoint< R2 >(1. / 3., R2(1. / 6., 1. / 6.))};
const GQuadratureFormular< R2 > TypeOfFE_Morley_3d::QFf(2, 3, P_QuadratureFormular_T_2_intp);

TypeOfFE_Morley_3d::TypeOfFE_Morley_3d( )
  : GTypeOfFE< Mesh3 >(TypeOfFE_Morley_3d::dfon, d+1, 1,
                       Element::ne *  QFe.n + Element::nf * 3 * QFf.n,
                       Element::ne * QFe.n + Element::nf * QFf.n, false, true) {
  assert(QFe.n);
  assert(QFf.n);
  R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.),
             R3(0., 0., 1.)};    // 4 ref tetrahedron vertices

  {
    // We build the interpolation pts on the edges of the reference tetrahedron:
    int i;
    i = 0;

    for (int e = 0; e < Element::ne; ++e) {
      for (int q = 0; q < QFe.n; ++q, ++i) {
        double x = QFe[q].x;
        this->PtInterpolation[i] =
          Pt[Element::nvedge[e][0]] * (1. - x) + Pt[Element::nvedge[e][1]] * (x);
      }
    }

    // We build the interpolation pts on the faces of the reference tetrahedron:
    // (the index i mustn't be reinitialised!)
    for (int f = 0; f < Element::nf; ++f) {
      for (int q = 0; q < QFf.n; ++q, ++i) {
        double x = QFf[q].x;
        double y = QFf[q].y;
        this->PtInterpolation[i] = Pt[Element::nvface[f][0]] * (1. - x - y) +
                                   Pt[Element::nvface[f][1]] * x + Pt[Element::nvface[f][2]] * y;
      }
    }
  }
  {
    // We build the indices in (13.1) : edge dofs
    int i = 0, p = 0;    // i is the k in (13.1) (chapter 13 ff++doc)
    int e;               // we will need e also below, in the part referred to faces

    for (e = 0; e < (Element::ne); e++) {    // loop on the 6 dofs

      for (int q = 0; q < QFe.n; ++q, ++p) {    // loop on the 2 edge quadrature pts
        {      // loop on the 3 components
          this->pInterpolation[i] = p;          // pk in (13.1)
          this->cInterpolation[i] = 0;          // jk in (13.1)
          this->dofInterpolation[i] = e;        // ik in (13.1)
          this->coefInterpolation[i] = 0.;      // alfak: we will fill them with 'set' (below)
                                                // because they depend on the tetrahedron
        }
      }
    }

    // We build the indices in (13.1) : face dofs
    // (the indices i and p mustn't be reinitialised)
    for (int f = 0; f < (Element::nf); f++) {    // loop on the 4 dofs

      for (int q = 0; q < QFf.n; ++q, ++p) {    // loop on the 3 face quadrature pts
        for (int c = 0; c < 3; c++, i++) {      // loop on the 3 components
          this->pInterpolation[i] = p;          // pk in (13.1)
          this->cInterpolation[i] = 1+c;          // jk in (13.1)
          this->dofInterpolation[i] = e + f;    // ik in (13.1)
          this->coefInterpolation[i] = 0.;      // alphak: we will fill them with 'set' (below)
                                                // because they depend on the tetrahedron
        }
      }
    }
  }
}
void TypeOfFE_Morley_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
                            int ocoef, int odf, int *nump) const {
    ffassert(0); // to do  F. Hecht
}
void TypeOfFE_Morley_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                           const RdHat &PHat, RNMK_ &val) const {
    ffassert(0); // to do  F. Hecht
}

  // a static variable to add the finite element to freefem++
  static TypeOfFE_P2Morley P2LagrangeP2Morley;
  static AddNewFE P2Morley("P2Morley", &P2LagrangeP2Morley);

static TypeOfFE_Morley_3d Morley_3d;    // TypeOfFE_Edge1_3d is the name of the class we defined
GTypeOfFE< Mesh3 > &GMorley3d(Morley_3d);    // GTypeOfFE<Mesh3> is the mother class
static AddNewFE3 AddGMorley3d("Morley3d", &GMorley3d);    // Edge13d will be the name used by the user

}    // namespace Fem2D
     // --- fin --
