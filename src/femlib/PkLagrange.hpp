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
// SUMMARY : Generic Pk Lagrange finite element class (Jan. 2008)
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

/*
Thank to the ARN FF2A3 grant
ref: ANR-07-CIS7-002-01
*/

#ifndef _PK_LAGRANGE_HPP_
#define _PK_LAGRANGE_HPP_

#include "FESpacen.hpp"

namespace Fem2D {
  template<class Rd,class E>
  static void SetPtPk(Rd *Pt, const int *dfon, int nn) {
    // P0, P1, P1b & P2
    const int dHat = E::RdHat::d;
    int k = 0;

    if (dfon[0]) {
      for (int i = 0; i <= dHat; ++i)
        Pt[k++] = Rd();

      for (int i = 0; i < dHat; ++i)
        Pt[i+1][i] = 1.;
    }

    if (dfon[1] && dHat !=1)
      for(int i = 0; i < E::ne; ++i)
        Pt[k++] = (Pt[E::nvedge[i][0]] + Pt[E::nvedge[i][1]])*0.5;

    if (dfon[dHat] == 1)
      Pt[k++] = Rd::diag(1./(dHat+1));

    if (nn != k) {
      cout << " " << nn << " != " << k << " - dHat = " << dHat << " " << dfon[0] << dfon[1] << dfon[2] << dfon[3] << " " << E::ne << endl;
      assert(nn == k);
    }
    if (verbosity > 9)
      cout << " Pk = " << KN_<Rd>(Pt, nn) << "\n";
  }

  //  a class of Lagrange Pk finite element
  template<class MMesh>
  class TypeOfFE_Lagrange: public GTypeOfFE<MMesh> {
    public:
      typedef MMesh Mesh;
      typedef typename Mesh::Element Element;
      typedef typename Element::RdHat RdHat;
      static const int dHat=RdHat::d;
      struct A4 {
        int dfon[4];

        A4(int k) {
          if (k == 0) {// P0
            dfon[0] = dfon[1] = dfon[2] = dfon[3] = 0;
            dfon[dHat] = 1;
          }
          else if (k == -1) { // P1b
            dfon[0] = 1;
            dfon[1] = dfon[2] = dfon[3] = 0;
            dfon[dHat] = 1;
          }
          else if (k == -2 && dHat==2) { // P2b
              dfon[0] = 1;
              dfon[1] = 1;
              dfon[2] = dfon[3] = 0;
              dfon[dHat] = 1;
            }
          else {
            dfon[0] = 1;
            dfon[1] = max(k - 1, 0);
            dfon[2] = dHat > 1 ? max(k - 2, 0) : 0;
            dfon[3] = dHat > 2 ? max(k - 3, 0) : 0;
          }

          if (verbosity > 9)
            cout << " A4 " << k << " " << dfon[0] << dfon[1] << dfon[2] << dfon[3] << endl;
        }
        operator const int * () const {return dfon;}
      };

      RdHat *Pt;
      TypeOfFE_Lagrange(int k):
      GTypeOfFE<Mesh>(A4(k), 1, k == -1 ? -1 : Max(k, 1), k <= 2, k == 0) {
        int n = this->NbDoF;
        if (verbosity > 9)
          cout << "\n +++ P" << k << ": ndof = " << n << endl;
        SetPtPk<RdHat, Element>(this->PtInterpolation, this->ndfOn(), this->NbDoF);
        if (verbosity > 9)
          cout << this->PtInterpolation << endl;
        for (int i = 0; i < n; i++) {
          this->pInterpolation[i] = i;
          this->cInterpolation[i] = 0.;
          this->dofInterpolation[i] = i;
          this->coefInterpolation[i] = 1.;
        }
      }
      ~TypeOfFE_Lagrange(){}

    private:
      TypeOfFE_Lagrange(const TypeOfFE_Lagrange &);
      void operator = (const TypeOfFE_Lagrange &);
  };
}

#endif //_PK_LAGRANGE_HPP_
