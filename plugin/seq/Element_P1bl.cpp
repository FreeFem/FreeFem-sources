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

#include "PkLagrange.hpp"
#include "ff++.hpp"
#include "AddNewFE.h"

// using namespace std;
namespace Fem2D {
  class TypeOfFE_P1Bubble2 : public TypeOfFE {
   public:
    static int Data[];
    static double Pi_h_coef[];
    TypeOfFE_P1Bubble2( ) : TypeOfFE(1, 0, 1, 1, Data, 1, 1, 4, 4, Pi_h_coef) {
      const R2 Pt[] = {R2(0, 0), R2(1, 0), R2(0, 1), R2(1. / 3., 1. / 3.)};

      for (int i = 0; i < NbDoF; i++) {
        pij_alpha[i] = IPJ(i, i, 0);
        P_Pi_h[i] = Pt[i];
      }
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
  };

  int TypeOfFE_P1Bubble2::Data[] = {0, 1, 2, 6, 0, 0, 0, 0, 0, 1, 2, 3,
                                    0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 4};
  double TypeOfFE_P1Bubble2::Pi_h_coef[] = {1., 1., 1., 1.};
  void TypeOfFE_P1Bubble2::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat,
                              RNMK_ &val) const {
    // const Triangle & K(FE.T);
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l[] = {1 - PHat.x - PHat.y, PHat.x, PHat.y};    // lb=l0*l1*l2*9.;
    int i0 = 0;

    if (l[1] < l[i0]) {
      i0 = 1;
    }

    if (l[2] < l[i0]) {
      i0 = 2;
    }

    int i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;
    R2 G[3] = {K.H(0), K.H(1), K.H(2)};
    double l0[3] = {l[i0] * 3, l[i1] - l[i0], l[i2] - l[i0]};
    R2 G0[3] = {G[i0] * 3, G[i1] - G[i0], G[i2] - G[i0]};
    if (val.N( ) < 4) {
      throwassert(val.N( ) >= 4);
    }

    throwassert(val.M( ) == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd[op_id]) {
      f0[3] = l0[0];
      f0[i1] = l0[1];
      f0[i2] = l0[2];
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      if (whatd[op_dx]) {
        RN_ f0x(val('.', 0, op_dx));
        f0x[3] = G0[0].x;
        f0x[i1] = G0[1].x;
        f0x[i2] = G0[2].x;
      }

      if (whatd[op_dy]) {
        RN_ f0y(val('.', 0, op_dy));
        f0y[3] = G0[0].y;
        f0y[i1] = G0[1].y;
        f0y[i2] = G0[2].y;
      }
    }
  }

  // Version 3D

  class TypeOfFE_P1blLagrange3d : public TypeOfFE_Lagrange< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef GFElement< Mesh3 > FElement;
    TypeOfFE_P1blLagrange3d( ) : TypeOfFE_Lagrange< Mesh3 >(-1) {}

    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
            RNMK_ &val) const;
  };

  void TypeOfFE_P1blLagrange3d::FB(const What_d whatd, const Mesh &, const Element &K,
                                   const RdHat &PHat, RNMK_ &val) const {
    R l[] = {1. - PHat.sum( ), PHat.x, PHat.y, PHat.z};
    const R d1 = dHat + 1.;
    const R d13 = d1 * d1 * d1;
    int i0 = 0;

    if (l[1] < l[i0]) {
      i0 = 1;
    }

    if (l[2] < l[i0]) {
      i0 = 2;
    }

    int i1 = (i0 + 1) % 4, i2 = (i0 + 2) % 4, i3 = (i0 + 3) % 4;
    double l0[4] = {l[i0] * 4, l[i1] - l[i0], l[i2] - l[i0], l[i3] - l[i0]};

    assert(val.N( ) >= Element::nv);
    assert(val.M( ) == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {
      f0[4] = l0[0];
      f0[i1] = l0[1];
      f0[i2] = l0[2];
      f0[i3] = l0[3];
    }

    if (whatd & Fop_D1) {
      R3 G[4];
      K.Gradlambda(G);
      R3 G0[4] = {G[i0] * 4, G[i1] - G[i0], G[i2] - G[i0], G[i3] - G[i0]};

      if (whatd & Fop_dx) {
        RN_ f0x(val('.', 0, op_dx));
        f0x[4] = G0[0].x;
        f0x[i1] = G0[1].x;
        f0x[i2] = G0[2].x;
        f0x[i3] = G0[3].x;
      }

      if (whatd & Fop_dy) {
        RN_ f0y(val('.', 0, op_dy));
        f0y[4] = G0[0].y;
        f0y[i1] = G0[1].y;
        f0y[i2] = G0[2].y;
        f0y[i3] = G0[3].y;
      }

      if (whatd & Fop_dz) {
        RN_ f0z(val('.', 0, op_dz));
        f0z[4] = G0[0].z;
        f0z[i1] = G0[1].z;
        f0z[i2] = G0[2].z;
        f0z[i3] = G0[3].z;
      }
    }
  }

  static TypeOfFE_P1blLagrange3d P1blLagrange3d;                     //
  GTypeOfFE< Mesh3 > &GP1blLagrange3d(P1blLagrange3d);               //
  static AddNewFE3 TypeOfFE_Edge1_3d("P1bl3d", &GP1blLagrange3d);    //

  // link with FreeFem++
  static TypeOfFE_P1Bubble2 P1Bulle2;
  // a static variable to add the finite element to freefem++
  static AddNewFE stP1Bulle2("P1bl", &P1Bulle2);
}    // namespace Fem2D
