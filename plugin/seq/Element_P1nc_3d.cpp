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
#include "PkLagrange.hpp"
#include "ff++.hpp"

#include "AddNewFE.h"

#include <iostream>
#include  <algorithm>


// Attention probleme de numerotation des inconnues
// -------------------------------------------------
// dans freefem, il y a un noeud par objets  sommet, arete, element.
// et donc la numerotation des dl dans l'element depend
// de l'orientation des aretes
//
/// ---------------------------------------------------------------
namespace Fem2D {

  class TypeOfFE_P1nc_3d : public TypeOfFE_Lagrange< Mesh3 > {
   public:
   typedef Mesh3 Mesh;
   typedef Mesh3::Element Element;
   typedef GFElement< Mesh3 > FElement;
      static int dfon[];
      static const int d = Mesh::Rd::d;
      static const GQuadratureFormular< R2 > &QFf;    // quadrature formula on a face
      static const GQuadratureFormular< R3 > &QFt;    // quadrature formula on a face
      TypeOfFE_P1nc_3d(): TypeOfFE_Lagrange<Mesh3>(-2,2) {  }

       void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
              RNMK_ &val) const;
     };






void TypeOfFE_P1nc_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                           const RdHat &PHat, RNMK_ &val) const {
    
    assert(val.N( ) >= 4);    // 4 degrees of freedom
    assert(val.M( ) == 1);     // 1 components
/*
   Dof  Numbering  1 dof / face   the dof a associated to a vertex face
  and the vertex a numbering in increase way.
 
 we have to numbering    î on ref element with
 
   So we need of a permutation p to insure the compatibilit between adjacent element
   
 */
    //  generated  with file Element_P2pnc_3d.edp
    // Warning  p(^i) =i
   // CC(i,j) = dof(j)(mo(i)); //
   // phi_k =sum_i C1(k,i) mo_i
 

 
    double l[4],mo[13];
    PHat.toBary(l);
    val=0; ;
    RN_ f0(val('.', 0, op_id));
    if (whatd & Fop_D0) {
        RN_ f0(val('.',0,0));
        f0[0] = 1-l[0]*3.;
        f0[1] = 1-l[1]*3.;
        f0[2] = 1-l[2]*3.;
        f0[3] = 1-l[3]*3.;

       
    }
    if (whatd & (Fop_D1 )) {
        const unsigned int iDx =whatd & Fop_dx , iDy =whatd & Fop_dy , iDz =whatd & Fop_dz; 
        R3 Dl[4];
        K.Gradlambda(Dl);
         int k=0;
        RN_ f0x(val('.', 0, op_dx));
        RN_ f0y(val('.', 0, op_dy));
        RN_ f0z(val('.', 0, op_dz));
        for(int i=0; i<4;++i)
        {
            if(iDx) f0x[i] = -Dl[i].x*3.;
            if(iDy) f0y[i] = -Dl[i].y*3.;
            if(iDz) f0z[i] = -Dl[i].z*3.;
        }
        
    }

}
 static TypeOfFE_P1nc_3d P1nc_3d;
 GTypeOfFE< Mesh3 > &Elm_P1nc_3d(P1nc_3d);

 static AddNewFE3 TFE_P2pnc_3d("P1nc3d", &Elm_P1nc_3d);
  
}    // namespace Fem2D

// --- fin --
