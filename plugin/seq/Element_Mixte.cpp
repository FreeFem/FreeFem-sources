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
// SUMMARY : Add finite elements
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep: lapack blas
//ff-c++-cpp-dep:
// *INDENT-ON* //

// RT1 and BDM1 conforme in H(div) space, ie u.n  continuous (n = is the normal)
// RT1Ortho , BDM1 are conforme in H(curl) space ie u.t  continuous  (t = is the tangent)
// RT1 , RT1Ortho  ( Neledec type I order 1)
// BDM1 (Brezzi, Douglas, Marini) ,  , BDM1Ortho  ( Neledec type I� order 1)
// the
// TDNNS1   finite element for symetrix 2x2 matrix conforme in H(div div)
// : Snn are continious
// See Thesis of Astrid Sabine Sinwel, A New Family of Mixed Finite Elements for Elasticity, 2009
// http://www.numa.uni-linz.ac.at/Teaching/PhD/Finished/sinwel-diss.pdf
// JOHANNES KEPLER UNIVERSITA,T LINZ
// Thesis of Sabine Zaglmayr, High Order Finite Element Methods for Electromagnetic Field
// Computation, 2006 JOHANNES KEPLER UNIVERSITA,T LINZ
// http://www.numerik.math.tugraz.at/~zaglmayr/pub/szthesis.pdf

// Realted files:
// LaplaceRT1.edp
// lame-TD-NSS.edp
// test-ElementMixte.edp

//ff-c++-LIBRARY-dep:   lapack blas
// lapack: dgetrf and dgetri

#include "ff++.hpp"
#include "AddNewFE.h"

#ifdef __LP64__
typedef int intblas;
typedef int integer;
#else
typedef long intblas;
typedef long integer;
#endif

typedef integer logical;
typedef float LAPACK_real;
typedef double doublereal;
typedef logical (*L_fp)( );
typedef integer ftnlen;
typedef complex< float > LAPACK_complex;
typedef complex< double > doublecomplex;
typedef void VOID;
#define complex LAPACK_complex
#define real LAPACK_real

#include "clapack.h"
#undef real
#undef complex
namespace Fem2D {
  // ------ P2 TD_NNS
  class TypeOfFE_TD_NNS0 : public TypeOfFE {
   public:
    static int Data[];

    TypeOfFE_TD_NNS0( )
      : TypeOfFE(3, 3, Data, 1, 1,
                 9,    // nb coef to build interpolation
                 3,    // np point to build interpolation
                 0) {
      const R2 Pt[] = {R2(0.5, 0.5), R2(0, 0.5), R2(0.5, 0)};
      // for the 3 vertices 6 coef
      int kk = 0;

      for (int e = 0; e < 3; ++e) {    // point d'integration sur l'arete e
        P_Pi_h[e] = Pt[e];
        pij_alpha[kk++] = IPJ(e, e, 0);
        pij_alpha[kk++] = IPJ(e, e, 1);
        pij_alpha[kk++] = IPJ(e, e, 2);
      }

      assert(P_Pi_h.N( ) == 3);
      assert(pij_alpha.N( ) == kk);
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const;
  };

  // on what     nu df on node node of df
  int TypeOfFE_TD_NNS0::Data[] = {3, 4, 5,    // support on what
                                  0, 0, 0,    // df on node
                                  0, 1, 2,    // th node of df
                                  0, 0, 0,    // df previou FE
                                  0, 1, 2,    // which df on prev
                                  0, 0, 0, 0, 0, 0, 3, 3, 3};
  void TypeOfFE_TD_NNS0::Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
    const Triangle &T(K.T);
    int k = 0;

    // integration sur les aretes
    for (int i = 0; i < 3; i++) {
      R2 N(T.Edge(i).perp( ));
      v[k++] = N.x * N.x;
      v[k++] = 2 * N.y * N.x;
      v[k++] = N.y * N.y;
    }

    assert(k == 9);
  }

  void TypeOfFE_TD_NNS0::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat,
                            RNMK_ &val) const {
    typedef double R;
    R cK = 2 * K.area;
    R2 Rl[3] = {K.Edge(0) / cK, K.Edge(1) / cK, K.Edge(2) / cK};    //
    /* bubble:
     * $ B_i = ((Rot L_i+1 ) (Rot L_(i+2)' ))^s  L_i$
     * s => symetrise ..
     */
    R S[3][3], S1[3][3];

    for (int i = 0; i < 3; ++i) {
      int i1 = (i + 1) % 3;
      int i2 = (i + 2) % 3;
      S[0][i] = -Rl[i1].x * Rl[i2].x;
      S[1][i] = -(Rl[i1].x * Rl[i2].y + Rl[i1].y * Rl[i2].x) * 0.5;
      S[2][i] = -Rl[i1].y * Rl[i2].y;
    }

    val = 0;
    KN< bool > wd(KN_< const bool >(whatd, last_operatortype));

    if (wd[op_id]) {
      for (int c = 0; c < 3; ++c) {
        for (int i = 0; i < 3; ++i) {
          val(i, c, op_id) = S[c][i];    // (c3-ll[i])/c3
        }
      }
    }
  }

  // ------ P2 TD_NNS
  class TypeOfFE_TD_NNS1 : public TypeOfFE {
   public:
    static int Data[];
    // double Pi_h_coef[];
    const QuadratureFormular1d &QFE;
    const GQuadratureFormular< R2 > &QFK;

    TypeOfFE_TD_NNS1( )
      : TypeOfFE(
          3 * 2 + 3, 3, Data, 2, 1,
          3 + 6 * 3 * QF_GaussLegendre2.n,                       // nb coef to build interpolation
          QuadratureFormular_T_1.n + 3 * QF_GaussLegendre2.n,    // np point to build interpolation
          0),
        QFE(QF_GaussLegendre2), QFK(QuadratureFormular_T_1) {
      int kk = 0, kp = 0;

      for (int p = 0; p < QFK.n; ++p) {
        P_Pi_h[kp++] = QFK[p];

        for (int c = 0; c < 3; c++) {
          pij_alpha[kk++] = IPJ(3 * 2 + c, p, c);
        }
      }

      for (int e = 0; e < 3; ++e) {
        for (int p = 0; p < QFE.n; ++p) {
          R2 A(TriangleHat[VerticesOfTriangularEdge[e][0]]);
          R2 B(TriangleHat[VerticesOfTriangularEdge[e][1]]);
          P_Pi_h[kp++] = B * (QFE[p].x) + A * (1. - QFE[p].x);    // X=0 => A  X=1 => B;
        }
      }

      for (int e = 0; e < 3; ++e) {
        for (int p = 0; p < QFE.n; ++p) {
          int pp = QFK.n + e * QFE.n + p;
          pij_alpha[kk++] = IPJ(2 * e + 0, pp, 0);
          pij_alpha[kk++] = IPJ(2 * e + 1, pp, 0);
          pij_alpha[kk++] = IPJ(2 * e + 0, pp, 1);
          pij_alpha[kk++] = IPJ(2 * e + 1, pp, 1);
          pij_alpha[kk++] = IPJ(2 * e + 0, pp, 2);
          pij_alpha[kk++] = IPJ(2 * e + 1, pp, 2);
        }
      }

      ffassert(P_Pi_h.N( ) == kp);
      ffassert(pij_alpha.N( ) == kk);
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
    void Pi_h_alpha(const baseFElement &K, KN_< double > &v) const;
  };

  // on what     nu df on node node of df
  int TypeOfFE_TD_NNS1::Data[] = {3, 3, 4, 4, 5, 5, 6, 6, 6,    // support on what
                                  0, 1, 0, 1, 0, 1, 0, 1, 2,    // df on node
                                  0, 0, 1, 1, 2, 2, 3, 3, 3,    // th node of df
                                  0, 0, 0, 0, 0, 0, 0, 0, 0,    // df previou FE
                                  0, 1, 2, 3, 4, 5, 6, 7, 8,    // which df on prev
                                  0, 0, 0, 0, 0, 0, 9, 9, 9};
  void TypeOfFE_TD_NNS1::Pi_h_alpha(const baseFElement &K, KN_< double > &v) const {
    const Triangle &T(K.T);
    int k = 0;

    // coef pour les 3 sommets  fois le 2 composantes
    // for (int i=0;i<3;i++)
    for (int p = 0; p < QFK.n; ++p) {    // wrong ...
      // 3 -1 -1
      // -1 3 -1    / 4
      // -1 -1 3
      R cc = QFK[p].a * T.area;
      v[k++] = cc;    // * T.area;
      v[k++] = cc;    // * T.area;
      v[k++] = cc;    // * T.area;
    }

    // integration sur les aretes
    for (int i = 0; i < 3; i++) {
      R s = T.EdgeOrientation(i);

      for (int p = 0; p < QFE.n; ++p) {
        R l0 = QFE[p].x, l1 = 1 - QFE[p].x;
        R p0 = (2 * l0 - l1) * 2;    // poly othogonaux to \lambda_1
        R p1 = (2 * l1 - l0) * 2;    // poly othogonaux to \lambda_0
        R cc1 = p0 * QFE[p].a;       //
        R cc0 = p1 * QFE[p].a;       //
        if (s < 0) {
          Exchange(cc1, cc0);    // exch lambda0,lambda1
        }

        R2 N(T.Edge(i).perp( ));
        v[k++] = cc0 * N.x * N.x;
        v[k++] = cc1 * N.x * N.x;
        v[k++] = cc0 * 2 * N.y * N.x;
        v[k++] = cc1 * 2 * N.y * N.x;
        v[k++] = cc0 * N.y * N.y;
        v[k++] = cc1 * N.y * N.y;
      }
    }

    ffassert(pij_alpha.N( ) == k);
  }

  void TypeOfFE_TD_NNS1::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat,
                            RNMK_ &val) const {
    typedef double R;
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R ll[3] = {l0, l1, l2};
    R2 Dl[3] = {K.H(0), K.H(1), K.H(2)};
    R cK = 2 * K.area;
    R2 Rl[3] = {K.Edge(0) / cK, K.Edge(1) / cK, K.Edge(2) / cK};    //
    /* bubble:
     * $ B_i = ((Rot L_i+1 ) (Rot L_(i+2)' ))^s  L_i$
     * s => symetrise ..
     */
    R S[3][3], S1[3][3];
    int ei0[3] = {1, 2, 0};
    int ei1[3] = {2, 0, 1};

    for (int i = 0; i < 3; ++i) {
      if (K.EdgeOrientation(i) < 0) {
        Exchange(ei0[i], ei1[i]);
      }

      int i1 = (i + 1) % 3;
      int i2 = (i + 2) % 3;
      S[0][i] = -Rl[i1].x * Rl[i2].x;
      S[1][i] = -(Rl[i1].x * Rl[i2].y + Rl[i1].y * Rl[i2].x) * 0.5;
      S[2][i] = -Rl[i1].y * Rl[i2].y;
    }

    // s[.] [i] =
    {    ////  compute the inv of S with lapack
      for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
          S1[i][j] = S[i][j];
        }
      }

      intblas N = 3, LWORK = 9;
      double WORK[9];
      intblas INFO, IPIV[4];

      dgetrf_(&N, &N, &(S1[0][0]), &N, IPIV, &INFO);
      ffassert(INFO == 0);
      dgetri_(&N, &(S1[0][0]), &N, IPIV, WORK, &LWORK, &INFO);
      ffassert(INFO == 0);
    }
    R B[3][3], BB[3][3];
    R cc = 3. / K.area;

    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        B[i][j] = S[i][j] * ll[j];
      }
    }

    for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
        BB[i][k] = 0.;

        for (int j = 0; j < 3; ++j) {
          BB[i][k] += cc * S[i][j] * ll[j] * S1[j][k];
        }
      }
    }

    if (verbosity > 1000) {
      cout << endl;
      cout << Rl[0] << " " << Rl[1] << ",  " << Rl[2] << endl;

      for (int i = 0; i < 3; ++i) {
        cout << " *****    " << BB[i][0] << " " << BB[i][1] << " " << BB[i][2] << "\t l " << ll[i]
             << endl;
      }
    }

    // the basic function are
    // the space S_i * ( a_i+b_1lambda_i)
    // the base :
    // tree egde function:
    // coefe*  S_i*( lambda_i - 1/3)   :  zero a barycenter
    // coefk*BB_i ,withh     B_i =   (S_i * lambda_i),   BB = B * S1 , ok because lambda_i = 1/3 a
    // bary so BB_ij = coefk/3  delta_ij  at G the barycenter.
    //
    KN< bool > wd(KN_< const bool >(whatd, last_operatortype));
    val = 0;

    throwassert(val.N( ) >= 6);
    throwassert(val.M( ) == 3);

    val = 0;

    if (wd[op_id]) {
      for (int c = 0; c < 3; ++c) {
        for (int i = 0; i < 3; ++i) {
          val(2 * i, c, op_id) = S[c][i] * (ll[ei0[i]] - ll[i]);
          val(2 * i + 1, c, op_id) = S[c][i] * (ll[ei1[i]] - ll[i]);
          val(i + 6, c, op_id) = BB[c][i];
        }
      }
    }

    if (wd[op_dx]) {
      for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
          BB[i][k] = 0.;

          for (int j = 0; j < 3; ++j) {
            BB[i][k] += cc * S[i][j] * Dl[j].x * S1[j][k];
          }
        }
      }

      for (int c = 0; c < 3; ++c) {
        for (int i = 0; i < 3; ++i) {
          val(2 * i, c, op_dx) = S[c][i] * (Dl[ei0[i]].x - Dl[i].x);
          val(2 * i + 1, c, op_dx) = S[c][i] * (Dl[ei1[i]].x - Dl[i].x);
          val(i + 6, c, op_dx) = BB[c][i];
        }
      }
    }

    if (wd[op_dy]) {
      for (int i = 0; i < 3; ++i) {
        for (int k = 0; k < 3; ++k) {
          BB[i][k] = 0.;

          for (int j = 0; j < 3; ++j) {
            BB[i][k] += cc * S[i][j] * Dl[j].y * S1[j][k];
          }
        }
      }

      for (int c = 0; c < 3; ++c) {
        for (int i = 0; i < 3; ++i) {
          val(2 * i, c, op_dy) = S[c][i] * (Dl[ei0[i]].y - Dl[i].y);
          val(2 * i + 1, c, op_dy) = S[c][i] * (Dl[ei1[i]].y - Dl[i].y);
          val(i + 6, c, op_dy) = BB[c][i];
        }
      }
    }
  }

  struct InitTypeOfRTk_2d {
    int k;       // order poly on edge
    int ndfi;    // nb of internal dof
    int npe;     // nb point on edge
    int ndf;     // nb dof

    KN< R > X;    // point on edge
    // KN<R> Pi_h_coef; // 1
    KN< int > Data;    // data of TypeOfFE
    const QuadratureFormular1d QFE;
    const GQuadratureFormular< R2 > &QFK;
    InitTypeOfRTk_2d(int KK)
      : k(KK), ndfi((k + 1) * (k)), npe(k + 1), ndf(3 * npe + ndfi), Data(5 * ndf + 6),
        QFE(-1 + 2 * npe, npe, GaussLegendre(npe), true), QFK(QuadratureFormular_T_5) {
      // int j = 0;
      int ndfe = ndf - ndfi;    //
      int o[6];

      o[0] = 0;

      for (int i = 1; i < 6; ++i) {
        o[i] = o[i - 1] + ndf;
      }

      for (int df = 0; df < ndf; ++df) {
        if (df < ndfe) {
          int e = df / npe;
          int n = df % npe;
          Data[o[0] + df] = 3 + e;
          Data[o[1] + df] = n;
          Data[o[2] + df] = e;
          Data[o[3] + df] = 0;
          Data[o[4] + df] = df;
        } else {
          int n = df - ndfe;
          Data[o[0] + df] = 6;
          Data[o[1] + df] = n;
          Data[o[2] + df] = 3;
          Data[o[3] + df] = 0;
          Data[o[4] + df] = df;
        }
      }

      Data[o[5] + 0] = 0;
      Data[o[5] + 1] = 0;
      Data[o[5] + 2] = 0;
      Data[o[5] + 3] = 0;
      Data[o[5] + 4] = ndf;
      Data[o[5] + 5] = ndf;    // end_dfcomp
    }
  };

  class TypeOfFE_RT1_2d : public InitTypeOfRTk_2d, public TypeOfFE {
   public:
    static double Pi_h_coef[];
    bool Ortho;

    TypeOfFE_RT1_2d(bool ortho)
      : InitTypeOfRTk_2d(1), TypeOfFE(ndf, 2, Data, 2, 1,
                                      2 * 2 * 3 * QFE.n + QFK.n * 4,    // nb coef mat interpole
                                      3 * QFE.n + QFK.n,                // nb P interpolation
                                      0),
        Ortho(ortho) {
      int kkk = 0, i = 0;

      for (int e = 0; e < 3; ++e) {
        for (int p = 0; p < QFE.n; ++p) {
          R2 A(TriangleHat[VerticesOfTriangularEdge[e][0]]);
          R2 B(TriangleHat[VerticesOfTriangularEdge[e][1]]);

          pij_alpha[kkk++] = IPJ(2 * e, i, 0);
          pij_alpha[kkk++] = IPJ(2 * e, i, 1);
          pij_alpha[kkk++] = IPJ(2 * e + 1, i, 0);
          pij_alpha[kkk++] = IPJ(2 * e + 1, i, 1);

          P_Pi_h[i++] = B * (QFE[p].x) + A * (1. - QFE[p].x);    // X=0 => A  X=1 => B;
        }
      }

      int i6 = 6, i7 = 7;

      for (int p = 0; p < QFK.n; ++p) {
        pij_alpha[kkk++] = IPJ(i6, i, 0);
        pij_alpha[kkk++] = IPJ(i6, i, 1);
        pij_alpha[kkk++] = IPJ(i7, i, 0);
        pij_alpha[kkk++] = IPJ(i7, i, 1);
        P_Pi_h[i++] = QFK[p];
      }

      ffassert(kkk == this->pij_alpha.N( ));
      ffassert(i == this->P_Pi_h.N( ));
    }

    void Pi_h_alpha(const baseFElement &K,
                    KN_< double > &v) const {    // compute the coef of interpolation ...
      const Triangle &T(K.T);
      int k = 0;
      R oe[3] = {T.EdgeOrientation(0), T.EdgeOrientation(1), T.EdgeOrientation(2)};

      for (int i = 0; i < 3; i++) {
        R2 E(Ortho ? T.Edge(i) : -T.Edge(i).perp( ));
        R s = oe[i];

        for (int p = 0; p < QFE.n; ++p) {
          R l0 = QFE[p].x, l1 = 1 - QFE[p].x;
          R p0 = (2 * l0 - l1) * 2;     // poly othogonaux to \lambda_1
          R p1 = (2 * l1 - l0) * 2;     // poly othogonaux to \lambda_0
          R cc1 = s * p0 * QFE[p].a;    //
          R cc0 = s * p1 * QFE[p].a;    //
          if (s < 0) {
            Exchange(cc1, cc0);    // exch lambda0,lambda1
          }

          v[k++] = cc0 * E.x;
          v[k++] = cc0 * E.y;
          v[k++] = cc1 * E.x;
          v[k++] = cc1 * E.y;
        }
      }

      R2 B[2] = {T.Edge(1), T.Edge(2)};
      if (Ortho) {
        B[0] = -B[0];
        B[1] = -B[1];
      } else {
        B[0] = B[0].perp( );
        B[1] = B[1].perp( );
      }

      double CK = 0.5;    // dof U= [u1,u2] > |K| int_K ( B_i.U )

      for (int p = 0; p < QFK.n; ++p) {
        double w = QFK[p].a * CK;
        v[k++] = w * B[0].x;
        v[k++] = w * B[0].y;
        v[k++] = w * B[1].x;
        v[k++] = w * B[1].y;
      }

      assert(k == this->pij_alpha.N( ));
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
  };

  // ENDOFCLASS TypeOfFE_PkEdge

  void TypeOfFE_RT1_2d::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &Phat,
                           RNMK_ &val) const {
    R2 X = K(Phat);
    R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])};
    R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x, l2 = Phat.y;
    R L[3] = {l0, l1, l2};
    R eo[] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};

    /*
     *
     * THE 2 DOF k=0,1  are: on edge e   f -> \int_e f \lambda_{e+k} . n_e
     * THE 2 internal dof are : f -> \int_K f e_i  where e_i is the canonical basis of R^2
     *
     *
     * so the basis function are
     *
     * let call \phi_i the basic fonction of RT0 (without orientation) so the normal is exterior.
     * \phi_i (X) = ( X- Q_i ) / (2 |K|) =  \lambda_{i+1} Curl( \lambda_{i+2}) - \lambda_{i+2} Curl(
     * \lambda_{i+1})
     *
     * edge function j=0,1
     * i1= i+j+1, i2= i+2-j  remark : {i,i1,i2} <=> {i,i+1,i+2}
     * \phi_i ( \lambda_{i1} - 4/3 \lambda_i) + 1/3 \phi_{i1}\lambda_{i1}
     *
     *
     *
     * we have 2 bubbles functions
     * fb_6 =      8    \phi_0   \lambda_{0} +   16 \phi_1   \lambda_{1}
     * fb_7 =     - 8   \phi_0   \lambda_{0} + 8 \phi_1   \lambda_{1}
     *
     *
     * such that  for i,j =0,1
     * int_K ( Bb_j  f_{6+i))  =   \delta_{ij}
     *
     *
     * so all basic d function are the sum of 3 function
     *
     * sum_{k=0}^2  c_k  phi_{p_k} lambda_{l_k}
     *
     */

    assert(val.N( ) >= ndf);
    assert(val.M( ) == 2);

    val = 0;

    R2 phi[3] = {X - Q[0], X - Q[1], X - Q[2]};    // phi * area *2
    if (Ortho) {
      phi[0] = phi[0].perp( );
      phi[1] = phi[1].perp( );
      phi[2] = phi[2].perp( );
    }

    int pI[8][3];    // store p_k
    int lI[8][3];    // store l_k
    R cI[8][3];      // store c_k
    int df = 0;
    R CKK = 2 * K.area;

    for (int e = 0; e < 3; ++e) {
      int i = e;
      int ii[2] = {(e + 1) % 3, (e + 2) % 3};
      R s = eo[e] / CKK;
      if (s < 0) {
        Exchange(ii[0], ii[1]);
      }

      for (int k = 0; k < 2; ++k, df++) {
        pI[df][0] = i;
        lI[df][0] = ii[k];
        cI[df][0] = s;

        pI[df][1] = i;
        lI[df][1] = i;
        cI[df][1] = -s * 4. / 3.;

        pI[df][2] = ii[k];
        lI[df][2] = ii[k];
        cI[df][2] = s / 3.;
      }
    }

    // FB  (x-Q_i) l_i l_j  =
    R s8 = 8 / CKK, s01 = s8;
    R cbb[] = {s8, 2 * s01, -s01, s8};    // { [ 8, 16], [ -8, 8] }

    // the 2 bubbles
    for (int k = 0; k < 2; ++k, df++) {    // k: ligne
      pI[df][0] = 0;                       // i
      lI[df][0] = 0;
      cI[df][0] = cbb[k];    //

      pI[df][1] = 1;    // i
      lI[df][1] = 1;
      cI[df][1] = cbb[k + 2];

      pI[df][2] = 2;
      lI[df][2] = 2;
      cI[df][2] = 0;
    }

    ffassert(df == 8);

    if (whatd[op_id]) {
      for (int df = 0; df < 8; ++df) {
        R2 fd(0., 0.);

        for (int k = 0; k < 3; ++k) {
          fd += (cI[df][k] * L[lI[df][k]]) * phi[pI[df][k]];
        }

        val(df, 0, op_id) = fd.x;
        val(df, 1, op_id) = fd.y;
      }
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      R2 DL[3] = {K.H(0), K.H(1), K.H(2)};
      R2 Dphix(1, 0);
      R2 Dphiy(0, 1);
      if (Ortho) {
        Dphix = R2(0, 1);
        Dphiy = R2(-1, 0);
      }    // x,y -> (-y,x)  Correct Bug Jan 2019 FH

      if (whatd[op_dx]) {
        for (int df = 0; df < 8; ++df) {
          R2 fd(0., 0.);

          for (int k = 0; k < 3; ++k) {
            fd += cI[df][k] * (DL[lI[df][k]].x * phi[pI[df][k]] + L[lI[df][k]] * Dphix);
          }

          val(df, 0, op_dx) = fd.x;
          val(df, 1, op_dx) = fd.y;
        }
      }

      if (whatd[op_dy]) {
        for (int df = 0; df < 8; ++df) {
          R2 fd(0., 0.);

          for (int k = 0; k < 3; ++k) {
            fd += cI[df][k] * (DL[lI[df][k]].y * phi[pI[df][k]] + L[lI[df][k]] * Dphiy);
          }

          val(df, 0, op_dy) = fd.x;
          val(df, 1, op_dy) = fd.y;
        }
      }

      if (whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
        cout << " to do FH RT1 dxx, dyy dxy " << endl;
        ffassert(0);
      }
    }
  }

  class TypeOfFE_RT2_2d : public InitTypeOfRTk_2d, public TypeOfFE {
   public:
    static double Pi_h_coef[];
    bool Ortho;

    TypeOfFE_RT2_2d(bool ortho)
      : InitTypeOfRTk_2d(2), TypeOfFE(ndf, 2, Data, 3, 1,
                                      2 * 3 * 3 * QFE.n + QFK.n * 4 * 3,    // nb coef mat interpole
                                      3 * QFE.n + QFK.n,                    // nb P interpolation
                                      0),
        Ortho(ortho) {
      int dofE = this->k + 1;                 // == 3
      int dofKs = (dofE - 1) * (dofE) / 2;    //= = 3 ..(

      ffassert(dofKs == 3);
      ffassert(dofE == 3);

      int kkk = 0, i = 0;

      for (int e = 0; e < 3; ++e) {
        for (int p = 0; p < QFE.n; ++p) {
          R2 A(TriangleHat[VerticesOfTriangularEdge[e][0]]);
          R2 B(TriangleHat[VerticesOfTriangularEdge[e][1]]);

          for (int l = 0; l < dofE; ++l) {
            pij_alpha[kkk++] = IPJ(dofE * e + l, i, 0);
            pij_alpha[kkk++] = IPJ(dofE * e + l, i, 1);
          }

          P_Pi_h[i++] = B * (QFE[p].x) + A * (1. - QFE[p].x);    // X=0 => A  X=1 => B;
        }
      }

      for (int p = 0; p < QFK.n; ++p) {
        int i6 = 3 * 3, i7 = i6 + 1;

        for (int l = 0; l < dofKs; ++l) {
          pij_alpha[kkk++] = IPJ(i6, i, 0);
          pij_alpha[kkk++] = IPJ(i6, i, 1);
          pij_alpha[kkk++] = IPJ(i7, i, 0);
          pij_alpha[kkk++] = IPJ(i7, i, 1);
          i6 += 2;
          i7 += 2;
        }

        P_Pi_h[i++] = QFK[p];
      }

      ffassert(kkk == this->pij_alpha.N( ));
      ffassert(i == this->P_Pi_h.N( ));
    }

    void Pi_h_alpha(const baseFElement &K,
                    KN_< double > &v) const {    // compute the coef of interpolation ...
      const Triangle &T(K.T);
      int k = 0;
      // magic fom:
      // inv of [[4!,3!,2!2!],[3!,2!2!,2!],[2!2!,3!,4!]]/5! =
      double c1[][3] = {{9, -18, 3} /* 0 */, {-18, 84, -18} /* 1 */, {3, -18, 9} /* 2 */};

      for (int i = 0; i < 3; i++) {
        R2 E(Ortho ? T.Edge(i) : -T.Edge(i).perp( ));
        R s = T.EdgeOrientation(i);

        for (int p = 0; p < QFE.n; ++p) {
          R l1 = QFE[p].x, l2 = 1 - QFE[p].x;
          R l11 = l1 * l1;
          R l22 = l2 * l2;
          R l21 = l2 * l1;
          R p0 = c1[0][0] * l11 + c1[0][1] * l21 + c1[0][2] * l22;    //
          R p1 = c1[1][0] * l11 + c1[1][1] * l21 + c1[1][2] * l22;    //
          R p2 = c1[2][0] * l11 + c1[2][1] * l21 + c1[2][2] * l22;    //
          R sa = s * QFE[p].a;
          R cc2 = sa * p0;    //
          R cc1 = sa * p1;    //
          R cc0 = sa * p2;    //
          if (s < 0) {
            swap(cc0, cc2);
          }

          v[k++] = cc0 * E.x;
          v[k++] = cc0 * E.y;
          v[k++] = cc1 * E.x;
          v[k++] = cc1 * E.y;
          v[k++] = cc2 * E.x;
          v[k++] = cc2 * E.y;
        }
      }

      R2 B[2] = {T.Edge(1), T.Edge(2)};
      if (Ortho) {
        B[0] = -B[0];
        B[1] = -B[1];
      } else {
        B[0] = B[0].perp( );
        B[1] = B[1].perp( );
      }

      double CK = 0.5;    // dof U= [u1,u2] > |K| int_K ( B_i.U )
      R ll[3];            //, lo[3];

      for (int p = 0; p < QFK.n; ++p) {
        double w = -QFK[p].a * CK;
        QFK[p].toBary(ll);
        ll[0] *= w;
        ll[1] *= w;
        ll[2] *= w;

        for (int l = 0; l < 3; ++l) {
          v[k++] = ll[l] * B[0].x;
          v[k++] = ll[l] * B[0].y;
          v[k++] = ll[l] * B[1].x;
          v[k++] = ll[l] * B[1].y;
        }
      }

      assert(k == this->pij_alpha.N( ));
    }

    void FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &Phat,
            RNMK_ &val) const;
  };

  void TypeOfFE_RT2_2d::FB(const bool *whatd, const Mesh &, const Triangle &K, const R2 &Phat,
                           RNMK_ &val) const {
    R2 X = K(Phat);
    R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])};
    R l0 = 1 - Phat.x - Phat.y, l1 = Phat.x, l2 = Phat.y;
    R L[3] = {l0, l1, l2};
    R eo[] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};

    assert(val.N( ) >= ndf);
    assert(val.M( ) == 2);

    val = 0;
    int p[15] = {0, 1, 2,  5,  4,  3,  6, 7,
                 8, 9, 10, 11, 12, 13, 14};    // Permutation for orinatation
    R2 Pm[18];                                 // all the momome function ..
    double cf[][6] = {{0, -5.5, 0, -2.5, -0.5, -1.5} /* 0 */,
                      {-1.25, -1.25, 0.25, -1, 0.25, -1} /* 1 */,
                      {-5.5, 0, -0.5, -1.5, 0, -2.5} /* 2 */,
                      {0, -2.5, 0, -5.5, -1.5, -0.5} /* 3 */,
                      {0.25, -1, -1.25, -1.25, -1, 0.25} /* 4 */,
                      {-0.5, -1.5, -5.5, 0, -2.5, 0} /* 5 */,
                      {-2.5, 0, -1.5, -0.5, 0, -5.5} /* 6 */,
                      {-1, 0.25, -1, 0.25, -1.25, -1.25} /* 7 */,
                      {-1.5, -0.5, -2.5, 0, -5.5, 0} /* 8 */,
                      {30, 90, -30, 180, 30, 60} /* 9 */,
                      {90, 30, 30, 60, -30, 180} /* 10 */,
                      {30, -180, -30, -90, -60, -30} /* 11 */,
                      {60, -120, 60, -60, 120, -60} /* 12 */,
                      {-120, 60, 120, -60, 60, -60} /* 13 */,
                      {-180, 30, -60, -30, -30, -90} /* 14 */};
    int Bii[][3] = {{0, 0, 0} /* 0 */,  {0, 1, 1} /* 1 */,  {0, 2, 2} /* 2 */,  {0, 1, 2} /* 3 */,
                    {0, 2, 0} /* 4 */,  {0, 0, 1} /* 5 */,  {1, 0, 0} /* 6 */,  {1, 1, 1} /* 7 */,
                    {1, 2, 2} /* 8 */,  {1, 1, 2} /* 9 */,  {1, 2, 0} /* 10 */, {1, 0, 1} /* 11 */,
                    {2, 0, 0} /* 12 */, {2, 1, 1} /* 13 */, {2, 2, 2} /* 14 */, {2, 1, 2} /* 15 */,
                    {2, 2, 0} /* 16 */, {2, 0, 1} /* 17 */};
    int fe[] = {1, 3, 2, 6, 10, 8, 12, 17, 13};
    int k6[] = {4, 5, 9, 11, 15, 16};
    R CKK = K.area * 2;
    R2 phi[3] = {X - Q[0], X - Q[1], X - Q[2]};    // phi * area *2
    if (Ortho) {
      phi[0] = phi[0].perp( );
      phi[1] = phi[1].perp( );
      phi[2] = phi[2].perp( );
    }

    for (int l = 0; l < 18; ++l) {
      int i = Bii[l][0];
      int j = Bii[l][1];
      int k = Bii[l][2];
      Pm[l] = phi[i] * (L[j] * L[k] / CKK);
    }

    // static int ddd=0;

    if (eo[0] < 0) {
      Exchange(p[0], p[2]);
    }

    if (eo[1] < 0) {
      Exchange(p[3], p[5]);
    }

    if (eo[2] < 0) {
      Exchange(p[6], p[8]);
    }

    double sg[15] = {eo[0], eo[0], eo[0], eo[1], eo[1], eo[1], eo[2], eo[2],
                     eo[2], 1.,    1.,    1.,    1.,    1.,    1.};

    if (whatd[op_id]) {
      for (int pdf = 0; pdf < 15; ++pdf) {
        int df = p[pdf];
        R2 fd(0., 0.);
        if (df < 9) {
          fd = Pm[fe[df]];    // edge function ..
        }

        for (int k = 0; k < 6; ++k) {
          fd += cf[df][k] * Pm[k6[k]];
        }

        fd *= sg[df];
        val(pdf, 0, op_id) = fd.x;
        val(pdf, 1, op_id) = fd.y;
      }
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      R2 DL[3] = {K.H(0), K.H(1), K.H(2)};
      R2 Dphi1(1, 0);    // D(phi.x)
      R2 Dphi2(0, 1);    // D(phi.y)
      R2 DxPm[18];
      R2 DyPm[18];

      if (Ortho) {
        Dphi1 = R2(0, -1);
        Dphi2 = R2(1, 0);
      }    // Correction Jan 2019 FH. 	// x,y -> (-y,x)

      for (int l = 0; l < 18; ++l) {
        int i = Bii[l][0];
        int j = Bii[l][1];
        int k = Bii[l][2];
        R Ljk = L[j] * L[k];
        R2 DLjk = L[j] * DL[k] + DL[j] * L[k];
        R2 DF1 = (Dphi1 * Ljk + phi[i].x * DLjk) / CKK;    // BUG ici Ortho ????????
        R2 DF2 = (Dphi2 * Ljk + phi[i].y * DLjk) / CKK;
        DxPm[l] = R2(DF1.x, DF2.x);
        DyPm[l] = R2(DF1.y, DF2.y);
      }

      if (whatd[op_dx]) {
        for (int pdf = 0; pdf < 15; ++pdf) {
          int df = p[pdf];
          R2 fd(0., 0.);
          if (df < 9) {
            fd = DxPm[fe[df]];    // edge function ..
          }

          for (int k = 0; k < 6; ++k) {
            fd += cf[df][k] * DxPm[k6[k]];
          }

          fd *= sg[df];
          val(pdf, 0, op_dx) = fd.x;
          val(pdf, 1, op_dx) = fd.y;
        }
      }

      if (whatd[op_dy]) {
        for (int pdf = 0; pdf < 15; ++pdf) {
          int df = p[pdf];
          R2 fd(0., 0.);
          if (df < 9) {
            fd = DyPm[fe[df]];    // edge function ..
          }

          for (int k = 0; k < 6; ++k) {
            fd += cf[df][k] * DyPm[k6[k]];
          }

          fd *= sg[df];
          val(pdf, 0, op_dy) = fd.x;
          val(pdf, 1, op_dy) = fd.y;
        }
      }

      if (whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
        cout << " to do FH RT2 dxx, dyy dxy " << endl;
        ffassert(0);
      }
    }
  }

  class TypeOfFE_BDM1_2d : public TypeOfFE {
   public:
    static int Data[];
    static double Pi_h_coef[];
    bool Ortho;
    const QuadratureFormular1d &QFE;
    TypeOfFE_BDM1_2d(bool ortho)
      : TypeOfFE(6, 2, Data, 1, 1,
                 2 * 2 * 3 * 2,    // nb coef mat interpole
                 3 * 2,            // nb P interpolation
                 0),
        QFE(QF_GaussLegendre2), Ortho(ortho) {
      int kkk = 0, i = 0;

      for (int e = 0; e < 3; ++e) {
        for (int p = 0; p < QFE.n; ++p) {
          R2 A(TriangleHat[VerticesOfTriangularEdge[e][0]]);
          R2 B(TriangleHat[VerticesOfTriangularEdge[e][1]]);

          pij_alpha[kkk++] = IPJ(2 * e, i, 0);
          pij_alpha[kkk++] = IPJ(2 * e, i, 1);
          pij_alpha[kkk++] = IPJ(2 * e + 1, i, 0);
          pij_alpha[kkk++] = IPJ(2 * e + 1, i, 1);

          P_Pi_h[i++] = B * (QFE[p].x) + A * (1. - QFE[p].x);    // X=0 => A  X=1 => B;
        }
      }

      ffassert(kkk == this->pij_alpha.N( ));
      ffassert(i == this->P_Pi_h.N( ));
    }

    void Pi_h_alpha(const baseFElement &K,
                    KN_< double > &v) const {    // compute the coef of interpolation ...
      const Triangle &T(K.T);
      int k = 0;

      for (int i = 0; i < 3; i++) {
        R2 E(Ortho ? T.Edge(i) : -T.Edge(i).perp( ));
        R s = T.EdgeOrientation(i);

        for (int p = 0; p < QFE.n; ++p) {
          R l0 = QFE[p].x, l1 = 1 - QFE[p].x;
          R p0 = s;                 // poly othogonaux to \lambda_1
          R p1 = -3 * (l0 - l1);    // poly othogonaux to \lambda_0
          R cc0 = p0 * QFE[p].a;    //
          R cc1 = p1 * QFE[p].a;    //
          v[k++] = cc0 * E.x;
          v[k++] = cc0 * E.y;
          v[k++] = cc1 * E.x;
          v[k++] = cc1 * E.y;
        }
      }

      assert(k == this->pij_alpha.N( ));
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
  };

  // ENDOFCLASS TypeOfFE_PkEdge
  int TypeOfFE_BDM1_2d::Data[] = {3, 3, 4, 4, 5, 5,    // support on what
                                  0, 1, 0, 1, 0, 1,    // df on node
                                  0, 0, 1, 1, 2, 2,    // th node of df
                                  0, 0, 0, 0, 0, 0,    // df previou FE
                                  0, 1, 2, 3, 4, 5,    // which df on prev
                                  0, 0, 0, 0, 6, 6};
  void TypeOfFE_BDM1_2d::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat,
                            RNMK_ &val) const {
    R2 X = K(PHat);
    R2 Q[] = {R2(K[0]), R2(K[1]), R2(K[2])};
    R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
    R L[3] = {l0, l1, l2};
    R2 Dl[3] = {K.H(0), K.H(1), K.H(2)};

    assert(val.N( ) >= 6);
    assert(val.M( ) == 2);

    val = 0;
    R cK = 2 * K.area;
    int ortho0 = 0, ortho1 = 1;
    R s1ortho = 1;
    if (Ortho) {
      ortho0 = 1;
      ortho1 = 0;
      s1ortho = -1;
    }

    if (whatd[op_id]) {
      for (int df = 0, e = 0; e < 3; ++e) {
        int e1 = (e + 1) % 3, e2 = (e + 2) % 3;
        R s = K.EdgeOrientation(e);
        R2 f1 = (X - Q[e]) * s / cK;
        R2 f2 = -(Dl[e1] * L[e2] + Dl[e2] * L[e1]).perp( );

        val(df, ortho0, op_id) = f1.x;
        val(df++, ortho1, op_id) = s1ortho * f1.y;

        val(df, ortho0, op_id) = f2.x;
        val(df++, ortho1, op_id) = s1ortho * f2.y;
      }
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      R2 Dphix(1, 0);
      R2 Dphiy(0, 1);

      if (whatd[op_dx]) {
        for (int df = 0, e = 0; e < 3; ++e) {
          int e1 = (e + 1) % 3, e2 = (e + 2) % 3;
          R s = K.EdgeOrientation(e);
          R2 f1 = R2(s / cK, 0.);
          R2 f2 = -(Dl[e1] * Dl[e2].x + Dl[e2] * Dl[e1].x).perp( );

          val(df, ortho0, op_dx) = f1.x;
          val(df++, ortho1, op_dx) = s1ortho * f1.y;

          val(df, ortho0, op_dx) = f2.x;
          val(df++, ortho1, op_dx) = s1ortho * f2.y;
        }
      }

      if (whatd[op_dy]) {
        for (int df = 0, e = 0; e < 3; ++e) {
          int e1 = (e + 1) % 3, e2 = (e + 2) % 3;
          R s = K.EdgeOrientation(e);
          R2 f1 = R2(0., s / cK);
          R2 f2 = -(Dl[e1] * Dl[e2].y + Dl[e2] * Dl[e1].y).perp( );

          val(df, ortho0, op_dy) = f1.x;
          val(df++, ortho1, op_dy) = s1ortho * f1.y;

          val(df, ortho0, op_dy) = f2.x;
          val(df++, ortho1, op_dy) = s1ortho * f2.y;
        }
      }
    }
  }

  // a static variable to add the finite element to freefem++
  static TypeOfFE_RT1_2d Elm_TypeOfFE_RT1_2d(false);          // RT1
  static TypeOfFE_RT1_2d Elm_TypeOfFE_RT1_2dOrtho(true);      // RT1ortho
  static TypeOfFE_RT2_2d Elm_TypeOfFE_RT2_2d(false);          // RT1
  static TypeOfFE_RT2_2d Elm_TypeOfFE_RT2_2dOrtho(true);      // RT1ortho
  static TypeOfFE_BDM1_2d Elm_TypeOfFE_BDM1_2d(false);        // BDM1
  static TypeOfFE_BDM1_2d Elm_TypeOfFE_BDM1_2dOrtho(true);    // BDM1ortho
  static TypeOfFE_TD_NNS0 Elm_TD_NNS;
  static TypeOfFE_TD_NNS1 Elm_TD_NNS1;
  static AddNewFE FE__TD_NNS("TDNNS0", &Elm_TD_NNS);
  static AddNewFE FE__TD_NNS1("TDNNS1", &Elm_TD_NNS1);
  static AddNewFE Elm__TypeOfFE_RT1_2d("RT1", &Elm_TypeOfFE_RT1_2d);
  static AddNewFE Elm__TypeOfFE_RT1_2dOrtho("RT1Ortho", &Elm_TypeOfFE_RT1_2dOrtho);
  static AddNewFE Elm__TypeOfFE_RT2_2d("RT2", &Elm_TypeOfFE_RT2_2d);
  static AddNewFE Elm__TypeOfFE_RT2_2dOrtho("RT2Ortho", &Elm_TypeOfFE_RT2_2dOrtho);
  static AddNewFE Elm__TypeOfFE_BDM1_2d("BDM1", &Elm_TypeOfFE_BDM1_2d);
  static AddNewFE Elm__TypeOfFE_BDM1_2dOrtho("BDM1Ortho", &Elm_TypeOfFE_BDM1_2dOrtho);
}    // namespace Fem2D
// --- fin --
