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
  // -------------------
  // ttdcnc1_ finite element fully discontinue.
  // -------------------
  class TypeOfFE_P1ttdcnc1_ : public TypeOfFE {
   public:
    static int Data[];
    static double Pi_h_coef[];

    TypeOfFE_P1ttdcnc1_( ) : TypeOfFE(0, 0, 3, 1, Data, 2, 1, 3, 3, Pi_h_coef) {
      const R2 Pt[] = {R2(0.5, 0.5), R2(0.0, 0.5), R2(0.5, 0.0)};

      for (int i = 0; i < NbDoF; i++) {
        pij_alpha[i] = IPJ(i, i, 0);
        P_Pi_h[i] = Pt[i];
      }
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;

    virtual R operator( )(const FElement &K, const R2 &PHat, const KN_< R > &u, int componante,
                          int op) const;
  };

  // on what   nu df on node  node of df
  int TypeOfFE_P1ttdcnc1_::Data[] = {6, 6, 6, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 3};
  ;
  double TypeOfFE_P1ttdcnc1_::Pi_h_coef[] = {1., 1., 1.};
  R TypeOfFE_P1ttdcnc1_::operator( )(const FElement &K, const R2 &PHat, const KN_< R > &u,
                                     int componante, int op) const {
    R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
    R r = 0;

    if (op == 0) {
      R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
      R ll0 = 1 - l0 * 2, ll1 = 1 - l1 * 2, ll2 = 1 - l1 * 2;
      r = u0 * ll0 + u1 * ll1 + ll2 * u2;
    } else {
      const Triangle &T = K.T;
      R2 D0 = T.H(0), D1 = T.H(1), D2 = T.H(2);
      if (op == 1) {
        r = -(D0.x * u0 + D1.x * u1 + D2.x * u2) * 2;
      } else {
        r = -(D0.y * u0 + D1.y * u1 + D2.y * u2) * 2;
      }
    }

    return r;
  }

  void TypeOfFE_P1ttdcnc1_::FB(const bool *whatd, const Mesh &, const Triangle &K,
                               const RdHat &PHat, RNMK_ &val) const {
    // const Triangle & K(FE.T);
    R2 A(K[0]), B(K[1]), C(K[2]);
    // l1(  cshrink1*(cshrink*((1,0)-G)+G)-G)+G  = 1
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;

    if (val.N( ) < 3) {
      throwassert(val.N( ) >= 3);
    }

    throwassert(val.M( ) == 1);

    val = 0;
    if (whatd[op_id]) {
      RN_ f0(val('.', 0, 0));
      f0[0] = 1 - l0 * 2;
      f0[1] = 1 - l1 * 2;
      f0[2] = 1 - l2 * 2;
    }

    if (whatd[op_dx] || whatd[op_dy]) {
      R2 Dl0(K.H(0)), Dl1(K.H(1)), Dl2(K.H(2));
      if (whatd[op_dx]) {
        RN_ f0x(val('.', 0, op_dx));
        f0x[0] = -Dl0.x * 2;
        f0x[1] = -Dl1.x * 2;
        f0x[2] = -Dl2.x * 2;
      }

      if (whatd[op_dy]) {
        RN_ f0y(val('.', 0, op_dy));
        f0y[0] = -Dl0.y * 2;
        f0y[1] = -Dl1.y * 2;
        f0y[2] = -Dl2.y * 2;
      }
    }
  }

  // link with FreeFem++
  static TypeOfFE_P1ttdcnc1_ P1dc1LagrangeP1dc1;
  // static TypeOfFE_LagrangeDC3d TypeOfFE_LagrangeDC3dtt(1);

  // a static variable to add the finite element to freefem++
  static AddNewFE P1dcLagrange("P1dcnc", &P1dc1LagrangeP1dc1);
  // static AddNewFE3  P1dttLagrange3d("P1dcnc3d",&TypeOfFE_LagrangeDC3dtt,"P1dcnc");
}    // namespace Fem2D

// --- fin --
