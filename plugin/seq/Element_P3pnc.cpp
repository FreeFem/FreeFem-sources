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
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE,  CEA
// AUTHORS :   BALAZI ATCHY NILLAMA Loïc et HECHT Frédéric 
// E-MAIL  : loic.balaziatchynillama@cea.fr et rederic.hecht@sorbonne-universite.fr

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */
/*
   the DoF are
   3 dof per edge e
    int_(e) \lambda^e_0 v, int_(e) \lambda^e_1 v, int_(e) \lambda^e_0 \lambda^e_1 v
    where \lambda^e_0 and \lambda^e_1 is the barycentric coordinate on oriented edge
   3 Dof on triangle K
    int_(K) \lambda_0 v, int_(e) \lambda_1 v, int_(e) \lambda_2
  where \lambda_0 , \lambda_1, \lambda_2 is the barycentric coordinate on traingle K
 
 and the Pk polynomial space is : P3 + span( bk \lambda_0 ,  bk \lambda_1)
 bk = (\lambda_0 - \lambda_1) * (\lambda_1 - \lambda_2) * (\lambda_2 - \lambda_0);
 */
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
  static const QuadratureFormular1d &QFE = QF_GaussLegendre4;
  static const GQuadratureFormular< R2 > &QFK = QuadratureFormular_T_5;

  class TypeOfFE_P3pnc : public TypeOfFE {
   public:
    static int Data[];

    TypeOfFE_P3pnc( ) : TypeOfFE(12, 1, Data, 4, 1, 3 * 3 * QFE.n + 3 * QFK.n, 3 * QFE.n + QFK.n, 0) {
      int p = 0, k = 0;
      const R2 Ph[] = {R2(0, 0), R2(1, 0), R2(0, 1)};

      for (int i = 0; i < 3; i++) {
        R2 A = Ph[(i + 1) % 3], B = Ph[(i + 2) % 3];
        for (int j = 0; j < QFE.n; j++) {
          R l0 = QFE[j].x, l1 = 1 - l0;
          cout << "k = " << k << " AND pij = " <<  pij_alpha.N( ) << endl;
          pij_alpha[k++] = IPJ(3 * i, p, 0);
          pij_alpha[k++] = IPJ(3 * i + 1, p, 0);
          pij_alpha[k++] = IPJ(3 * i + 2, p, 0);
          P_Pi_h[p++] = A * l0 + B * l1;
        }
      }

      for (int i = 0; i < QFK.n; i++) {
        cout << "k = " << k << " AND pij = " <<  pij_alpha.N( ) << endl;
        pij_alpha[k++] = IPJ(9, p, 0);
        pij_alpha[k++] = IPJ(10, p, 0);
        pij_alpha[k++] = IPJ(11, p, 0);
        P_Pi_h[p++] = QFK[i];
      }
      cout << " QFE.n " << QFE.n << endl;
      cout << " QFK.n " << QFK.n << endl;
      cout << "k = " << k << " AND pij = " <<  pij_alpha.N( ) << endl;
      cout << "p = " << p << " AND P_Pi_h = " <<  P_Pi_h.N( )<< endl;
      ffassert(k == this->pij_alpha.N( ));
      ffassert(p == this->P_Pi_h.N( ));
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const;
  };

  // const QuadratureFormular1d & TypeOfFE_P2pnc::QFE=QF_GaussLegendre3;
  // const GQuadratureFormular<R2> &  TypeOfFE_P2pnc::QFK=QuadratureFormular_T_5;

  int TypeOfFE_P3pnc::Data[] = {3, 3, 3,  4, 4, 4,  5, 5, 5,  6, 6,  6,    // the support number  of the node of the df
                                0, 1, 2,  0, 1, 2,  0, 1, 2,  0, 1, 2,    // the number of the df on  the node
                                0, 0, 0,  1, 1, 1,  2, 2, 2,  3, 3, 3,    // the node of the df
                                0, 0, 0,  0, 0, 0,  0, 0, 0,  0, 0, 0,    // the df come from which FE (generaly 0)
                                0, 1, 2,  3, 4, 5,  6, 7, 8,  9,10,11,    // which are de df on sub FE
                                0,    // for each compontant $j=0,N-1$ it give the sub FE associated
                                0, 12};
  void TypeOfFE_P3pnc::Pi_h_alpha(
    const baseFElement &K, KN_< double > &v) const {    // compute the coef of interpolation ...
    const Triangle &T(K.T);
    int k = 0;
    R oe[3] = {T.EdgeOrientation(0), T.EdgeOrientation(1), T.EdgeOrientation(2)};
    static int ddd = 10;

    ++ddd;

    for (int i = 0; i < 3; i++) {
      for (int p = 0; p < QFE.n; ++p) {
          double x = QFE[p].x;
        R L0 = QFE[p].x , L1 = 1-QFE[p].x , L2 = L0*L1 ;
        if (oe[i] < 0) {
          swap(L0, L1);  //depend on edge orientation //No need to swap L2 since symmetric
        }

        if (ddd < 3) {
          cout << p << " " << oe[i] << " " << L0 << " " << L1 << endl;
        }

        R cc0 = L0 * QFE[p].a;
        R cc1 = L1 * QFE[p].a;
        R cc2 = L2 * QFE[p].a;
        v[k++] = cc0;
        v[k++] = cc1;
        v[k++] = cc2;
      }
    }

    for (int p = 0; p < QFK.n; ++p) {
      double w = QFK[p].a;
      R l1 = QFK[p].x, l2 = QFK[p].y, l0 = 1 - l1 - l2;
      v[k++] = l0 * w;
      v[k++] = l1 * w;
      v[k++] = l2 * w;
    }

    ffassert(k == this->pij_alpha.N( ));
  }

  void TypeOfFE_P3pnc::FB(const bool *whatd, const Mesh &TH, const Triangle &K, const RdHat &PHat,
                          RNMK_ &val) const {
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R bk = (l0 - l1) * (l1 - l2) * (l2 - l0);
    R oe[3] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
    R l12[] = {
      l0 * l0 * l0, l1 * l1 * l1, l2 * l2 * l2,                                            //3
      l0 * l0 * l1, l0 * l0 * l2, l1 * l1 * l0, l1 * l1 * l2, l2 * l2 * l0, l2 * l2 * l1,  //6
      l0 * l1 * l2,                                                                        //2
      bk * l0, bk * l1                                                                     // P4 element
     };    // 12 monome

    /*
    double C1[12][12] ={
    {-1.2,   6,        1.2,    -8.4,    10.8,     10.8,    -46.8,   -10.8,    20.4,    12,     0,     -84  },
    {-1.2,   1.2,      6,      -73.2,   75.6,     73.2,    -63.6,   -73.2,    37.2,    12,     -84,   -84  },
    {1.2,    -15.6,    -15.6,  298.8,   -205.2,   -306,    349.2,    198,     -154.8,  -180,   252,   504  },
    {1.2,    -1.2,     6,      73.2,    -63.6,    -73.2,   75.6,     37.2,    -73.2,   12,     84,    84   },
    {6,      -1.2,     1.2,    10.8,    -46.8,    -8.4,    10.8,     20.4,    -10.8,   12 ,    84,    0    },
    {-15.6,  1.2,      -15.6,  -306,    349.2,    298.8,   -205.2,   -154.8,   198,    -180,  -504,   -252 },
    {6,      1.2,      -1.2,   -46.8,   10.8,     20.4,    -10.8,    -8.4,     10.8,   12,    -84,    0    },
    {1.2,    6,        -1.2,   20.4,    -10.8,    -46.8,   10.8,     10.8,     -8.4,   12,    0,      84   },
    {-15.6,  -15.6,    1.2,    97.2,    -54,      97.2,    -54,      46.8,     46.8,   -180,  252,    -252 },
    {-9.6,   4.8,      4.8,    189.6,   21.6,     -192,    86.4,     -24,     -81.6,   60,    84,     168  },
    {4.8,    -9.6,     4.8,    -192,    86.4,     189.6,   21.6,     -81.6,   -24,     60,   -168,    -84  },
    {4.8,    4.8,      -9.6,   2.4,     -108,     2.4,    -108,      105.6,   105.6,   60,   84,      -84  },
  }; */


      double C1[12][12] ={
      {  -1.2,    -1.2,      1.2,      1.2,         6.,        -15.6,       6.,        1.2,       -15.6,      -9.6,     4.8,      4.8     }, //dof 1
      {   6.,     1.2,      -15.6,     -1.2,        -1.2,       1.2,        1.2,       6.,        -15.6,      4.8,      -9.6,     4.8     }, //dof 2
      {   1.2,    6.,       -15.6,      6.,         1.2,       -15.6,      -1.2,      -1.2,       1.2,        4.8 ,     4.8,      -9.6    }, //dof 3
      {  -8.4,    -73.2,    298.8,      73.2,       10.8,      -306.,      -46.8,      20.4,      97.2,      189.6,    -192.,      2.4    },  //dof 4
      {  10.8,    75.6,     -205.2,     -63.6,      -46.8,     349.2,      10.8,      -10.8,      -54.,      21.6,      86.4,     -108.   }, //dof 5
      {  10.8,    73.2,     -306.,      -73.2,      -8.4,      298.8,      20.4,      -46.8,      97.2,      -192.,     189.6,    2.4     }, //dof 6
      { -46.8,    -63.6,    349.2,      75.6,       10.8,      -205.2,     -10.8,     10.8,      -54.,       86.4 ,     21.6,    -108.    }, //dof 7
      { -10.8,    -73.2,    198.,       37.2,       20.4,      -154.8,     -8.4,      10.8,      46.8,       -24. ,     -81.6,    105.6   }, //dof 8
      {  20.4,    37.2,     -154.8,    -73.2,       -10.8,     198.,       10.8,      -8.4,      46.8,       -81.6 ,    -24.,     105.6   }, //dof 9
      {  12.,     12.,      -180.,      12.,        12.,       -180.,      12.,       12.,       -180.,      60. ,      60.,      60.     }, //dof 10
      {   0.,     -84.,     252.,       84.,        84.,       -504.,     -84.,       0.,        252.,       84.,      -168.,     84.     }, //dof 11
      { -84.,     -84.,     504.,       84.,        0.,        -252.,     -0.,        84.,      -252.,       168.,     -84.,     -84.     }, //dof 12
    };


    //  C1 = C^1
    //  C(i,j) = dof(i,mom_j)
    //  dof(k,sum_j C1(j,i)*mom_j) = delta_ik
    //  sum_j C1(j,i)*dof(k,mom_j) = Id
   //  sum_j C_kj C1(j,i)  = Id
   //  sum_j dof(k,mom_j) C1(j,i) = Id
   //  C_kj= dof(k,mom_j)

   // formule magic :
   //  sur un simplex en dim d , l_0, , .. l_d coord barycentric
    //  int_K (\Prod_i=0^d l_i^n_i) = |K| ( d! \Prod_i n_i! / (d+ sum_i n_i)!)
    // n_i == 0 ,
      // n_i = delta_ik => int_K l_k = |K|/(d+1)
      //  int_K l_k l_j = |K|(1+delta_kj) /((d+1)((d+2))
      //

    if (val.N( ) < 12) {
      throwassert(val.N( ) >= 12);
    }

    throwassert(val.M( ) == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));
    int p[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}; //no need to swap 2, 5, 8 since symmetric


    if (oe[0] < 0) {
      swap(p[0], p[1]);
    }

    if (oe[1] < 0) {
      swap(p[3], p[4]);
    }

    if (oe[2] < 0) {
      swap(p[6], p[7]);
    }

 //No need to swap elemet 2, 5 and 8 since they are symmetric = do not depends on edge orientation

    if (whatd[op_id]) {
      for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
          f0[p[i]] += C1[j][i] * l12[j];
        }
      }
    }

    if (whatd[op_dx] || whatd[op_dy]) {
      R2 D0 = K.H(0), D1 = K.H(1), D2 = K.H(2);
      R2 Dbk =  (l1 - l2) * (l2 - l0) * (D0 - D1)
              + (l0 - l1) * (l2 - l0) * (D1 - D2)
              + (l0 - l1) * (l1 - l2) * (D2 - D0);
      R2 D12[] = {
                 D0 * (l0 * l0*3.),
                 D1 * (l1 * l1*3.),
                 D2 * (l2 * l2*3.),
                 D0 * (l0 * l1 *2.) + D1 * (l0 * l0),
                 D0 * (l0 * l2 *2.) + D2 * (l0 * l0),
                 D1 * (l1 * l0 *2.) + D0 * (l1 * l1),
                 D1 * (l1 * l2 *2.) + D2 * (l1 * l1),
                 D2 * (l2 * l0 *2.) + D0 * (l2 * l2),
                 D2 * (l2 * l1 *2.) + D1 * (l2 * l2),
                 D0 * (l1 * l2)  + D1* (l0 * l2) + D2* (l0 * l1) ,
                 Dbk*l0 + D0 *bk, Dbk*l1 + D1 *bk
                };

      if (whatd[op_dx]) {
        RN_ f0x(val('.', 0, op_dx));

        for (int i = 0; i < 12; ++i) {
          for (int j = 0; j < 12; ++j) {
            f0x[p[i]] += C1[j][i] * D12[j].x;
          }
        }
      }

      if (whatd[op_dy]) {
        RN_ f0y(val('.', 0, op_dy));

        for (int i = 0; i < 12; ++i) {
          for (int j = 0; j < 12; ++j) {
            f0y[p[i]] += C1[j][i] * D12[j].y;
          }
        }
      }
    }
  }
}    // namespace Fem2D
static TypeOfFE_P3pnc VTypeOfFE_P3pnc;
static AddNewFE P1dcLagrange("P3pnc", &VTypeOfFE_P3pnc);

static void finit( ) {
  // link with FreeFem++
}

LOADFUNC(finit)

// --- fin --
