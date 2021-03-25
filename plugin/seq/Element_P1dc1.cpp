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
  // ttdc1_ finite element fully discontinue.
  // -------------------
  class TypeOfFE_P1ttdc1_ : public TypeOfFE {
   public:
    static int Data[];
    static double Pi_h_coef[];
    static const R2 G;
    static const R cshrink;
    static const R cshrink1;
    // (1 -1/3)*
    static R2 Shrink(const R2 &P) { return (P - G) * cshrink + G; }

    static R2 Shrink1(const R2 &P) { return (P - G) * cshrink1 + G; }

    TypeOfFE_P1ttdc1_( ) : TypeOfFE(0, 0, 3, 1, Data, 1, 1, 3, 3, Pi_h_coef) {
      const R2 Pt[] = {Shrink(R2(0, 0)), Shrink(R2(1, 0)), Shrink(R2(0, 1))};

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

  const R2 TypeOfFE_P1ttdc1_::G(1. / 3., 1. / 3.);
  const R TypeOfFE_P1ttdc1_::cshrink = 1;
  const R TypeOfFE_P1ttdc1_::cshrink1 = 1. / TypeOfFE_P1ttdc1_::cshrink;

  class TypeOfFE_P2ttdc1_ : public TypeOfFE {
   public:
    static int Data[];
    static double Pi_h_coef[];
    static const R2 G;
    static const R cshrink;
    static const R cshrink1;
    static R2 Shrink(const R2 &P) { return (P - G) * cshrink + G; }

    static R2 Shrink1(const R2 &P) { return (P - G) * cshrink1 + G; }

    TypeOfFE_P2ttdc1_( ) : TypeOfFE(0, 0, 6, 1, Data, 3, 1, 6, 6, Pi_h_coef) {
      const R2 Pt[] = {Shrink(R2(0, 0)),     Shrink(R2(1, 0)),   Shrink(R2(0, 1)),
                       Shrink(R2(0.5, 0.5)), Shrink(R2(0, 0.5)), Shrink(R2(0.5, 0))};

      for (int i = 0; i < NbDoF; i++) {
        pij_alpha[i] = IPJ(i, i, 0);
        P_Pi_h[i] = Pt[i];
      }
    }

    void FB(const bool *whatd, const Mesh &Th, const Triangle &K, const RdHat &PHat,
            RNMK_ &val) const;
  };

  // on what   nu df on node  node of df
  int TypeOfFE_P1ttdc1_::Data[] = {6, 6, 6, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 3};
  int TypeOfFE_P2ttdc1_::Data[] = {6, 6, 6, 6, 6, 6, 0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0,
                                   0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0, 6};
  double TypeOfFE_P1ttdc1_::Pi_h_coef[] = {1., 1., 1.};
  double TypeOfFE_P2ttdc1_::Pi_h_coef[] = {1., 1., 1., 1., 1., 1.};
  const R2 TypeOfFE_P2ttdc1_::G(1. / 3., 1. / 3.);
  const R TypeOfFE_P2ttdc1_::cshrink = 1;
  const R TypeOfFE_P2ttdc1_::cshrink1 = 1. / TypeOfFE_P2ttdc1_::cshrink;
  R TypeOfFE_P1ttdc1_::operator( )(const FElement &K, const R2 &P1Hat, const KN_< R > &u,
                                   int componante, int op) const {
    R2 PHat = Shrink1(P1Hat);
    R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));
    R r = 0;

    if (op == 0) {
      R l0 = 1 - PHat.x - PHat.y, l1 = PHat.x, l2 = PHat.y;
      r = u0 * l0 + u1 * l1 + l2 * u2;
    } else {
      const Triangle &T = K.T;
      R2 D0 = T.H(0) * cshrink1, D1 = T.H(1) * cshrink1, D2 = T.H(2) * cshrink1;
      if (op == 1) {
        r = D0.x * u0 + D1.x * u1 + D2.x * u2;
      } else {
        r = D0.y * u0 + D1.y * u1 + D2.y * u2;
      }
    }

    return r;
  }

  void TypeOfFE_P1ttdc1_::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat,
                             RNMK_ &val) const {
    R2 P = Shrink1(PHat);

    // const Triangle & K(FE.T);
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - P.x - P.y, l1 = P.x, l2 = P.y;

    if (val.N( ) < 3) {
      throwassert(val.N( ) >= 3);
    }

    throwassert(val.M( ) == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd[op_id]) {
      f0[0] = l0;
      f0[1] = l1;
      f0[2] = l2;
    }

    if (whatd[op_dx] || whatd[op_dy]) {
      R2 Dl0(K.H(0) * cshrink1), Dl1(K.H(1) * cshrink1), Dl2(K.H(2) * cshrink1);

      if (whatd[op_dx]) {
        RN_ f0x(val('.', 0, op_dx));
        f0x[0] = Dl0.x;
        f0x[1] = Dl1.x;
        f0x[2] = Dl2.x;
      }

      if (whatd[op_dy]) {
        RN_ f0y(val('.', 0, op_dy));
        f0y[0] = Dl0.y;
        f0y[1] = Dl1.y;
        f0y[2] = Dl2.y;
      }
    }
  }

  void TypeOfFE_P2ttdc1_::FB(const bool *whatd, const Mesh &, const Triangle &K, const RdHat &PHat,
                             RNMK_ &val) const {
    R2 P = Shrink1(PHat);

    // const Triangle & K(FE.T);
    R2 A(K[0]), B(K[1]), C(K[2]);
    R l0 = 1 - P.x - P.y, l1 = P.x, l2 = P.y;
    R l4_0 = (4 * l0 - 1), l4_1 = (4 * l1 - 1), l4_2 = (4 * l2 - 1);

    throwassert(val.N( ) >= 6);
    throwassert(val.M( ) == 1);

    val = 0;
    if (whatd[op_id]) {
      RN_ f0(val('.', 0, op_id));
      f0[0] = l0 * (2 * l0 - 1);
      f0[1] = l1 * (2 * l1 - 1);
      f0[2] = l2 * (2 * l2 - 1);
      f0[3] = 4 * l1 * l2;    // oppose au sommet 0
      f0[4] = 4 * l0 * l2;    // oppose au sommet 1
      f0[5] = 4 * l1 * l0;    // oppose au sommet 3
    }

    if (whatd[op_dx] || whatd[op_dy] || whatd[op_dxx] || whatd[op_dyy] || whatd[op_dxy]) {
      R2 Dl0(K.H(0) * cshrink1), Dl1(K.H(1) * cshrink1), Dl2(K.H(2) * cshrink1);
      if (whatd[op_dx]) {
        RN_ f0x(val('.', 0, op_dx));
        f0x[0] = Dl0.x * l4_0;
        f0x[1] = Dl1.x * l4_1;
        f0x[2] = Dl2.x * l4_2;
        f0x[3] = 4 * (Dl1.x * l2 + Dl2.x * l1);
        f0x[4] = 4 * (Dl2.x * l0 + Dl0.x * l2);
        f0x[5] = 4 * (Dl0.x * l1 + Dl1.x * l0);
      }

      if (whatd[op_dy]) {
        RN_ f0y(val('.', 0, op_dy));
        f0y[0] = Dl0.y * l4_0;
        f0y[1] = Dl1.y * l4_1;
        f0y[2] = Dl2.y * l4_2;
        f0y[3] = 4 * (Dl1.y * l2 + Dl2.y * l1);
        f0y[4] = 4 * (Dl2.y * l0 + Dl0.y * l2);
        f0y[5] = 4 * (Dl0.y * l1 + Dl1.y * l0);
      }

      if (whatd[op_dxx]) {
        RN_ fxx(val('.', 0, op_dxx));

        fxx[0] = 4 * Dl0.x * Dl0.x;
        fxx[1] = 4 * Dl1.x * Dl1.x;
        fxx[2] = 4 * Dl2.x * Dl2.x;
        fxx[3] = 8 * Dl1.x * Dl2.x;
        fxx[4] = 8 * Dl0.x * Dl2.x;
        fxx[5] = 8 * Dl0.x * Dl1.x;
      }

      if (whatd[op_dyy]) {
        RN_ fyy(val('.', 0, op_dyy));
        fyy[0] = 4 * Dl0.y * Dl0.y;
        fyy[1] = 4 * Dl1.y * Dl1.y;
        fyy[2] = 4 * Dl2.y * Dl2.y;
        fyy[3] = 8 * Dl1.y * Dl2.y;
        fyy[4] = 8 * Dl0.y * Dl2.y;
        fyy[5] = 8 * Dl0.y * Dl1.y;
      }

      if (whatd[op_dxy]) {
        assert(val.K( ) > op_dxy);
        RN_ fxy(val('.', 0, op_dxy));

        fxy[0] = 4 * Dl0.x * Dl0.y;
        fxy[1] = 4 * Dl1.x * Dl1.y;
        fxy[2] = 4 * Dl2.x * Dl2.y;
        fxy[3] = 4 * (Dl1.x * Dl2.y + Dl1.y * Dl2.x);
        fxy[4] = 4 * (Dl0.x * Dl2.y + Dl0.y * Dl2.x);
        fxy[5] = 4 * (Dl0.x * Dl1.y + Dl0.y * Dl1.x);
      }
    }
  }

  //
  // end ttdc1_



  // 3d
  template<class MMesh> 
   void SetPtPkDC(typename MMesh::Element::RdHat *Pt, int kk, int nn, R cc = 1) ;
  //volume
  template<>
   void SetPtPkDC<Mesh3>(R3 *Pt, int kk, int nn, R cc) {    // P0 P1 et P2 , P1b
    const int d = 3;
    int n = 0;
    double dK = kk;
    double cc1 = 1 - cc;                    //
    const R3 G = R3::diag(1. / (d + 1));    // barycenter

    for (int i = 0; i <= kk; ++i) {
      for (int j = 0; j <= kk - i; ++j) {
        for (int k = 0; k <= kk - i - j; ++k) {
          int l = kk - i - j - k;
          ffassert(l >= 0 && l <= kk);
          Pt[n++] = R3(k / dK, j / dK, i / dK) * cc + G * cc1;
        }
      }
    }

    ffassert(n == nn);
    if (verbosity > 9) {
      cout << " Pkdc = " << KN_< R3 >(Pt, nn) << "\n";
    }
  }
  
  //surface 
  template<>
   void SetPtPkDC<MeshS>(R2 *Pt, int kk, int nn, R cc) {
    const int dHat = 2;
    int n = 0;
    double dK = kk;
    double cc1 = 1 - cc;                 
    const R2 G = R2::diag(1. / (dHat + 1));    // barycenter

    for (int i = 0; i <= kk; ++i) 
      for (int j = 0; j <= kk - i; ++j)
          Pt[n++] = R2(j / dK, i / dK) * cc + G * cc1;
    ffassert(n == nn);
    if (verbosity > 9) {
      cout << " Pkdc = " << KN_< R2 >(Pt, nn) << "\n";
    }
  }
  
 //curve
  template<>
   void SetPtPkDC<MeshL>(R1 *Pt, int kk, int nn, R cc) {
    const int dHat = 2;
    int n = 0;
    double dK = kk;
    double cc1 = 1 - cc;                 
    const R G = R(1./2.);    // barycenter

    for (int i = 0; i <= kk; ++i) 
          Pt[n++] = R1(i / dK) * cc + G * cc1;
    ffassert(n == nn);
    if (verbosity > 9) {
      cout << " Pkdc = " << KN_< R1 >(Pt, nn) << "\n";
    }
  }  
  

  
  // 3d
  template<class MMesh>
  class TypeOfFE_LagrangeDC3d : public GTypeOfFE< MMesh > {
    // typedef typename  MMesh Mesh;

   public: 
	typedef GFElement<MMesh> FElement;
    typedef typename MMesh::Element Element;
    typedef typename Element::Rd Rd;
    typedef typename Element::RdHat RdHat;
    static const int d = Rd::d;
	static const int dHat = RdHat::d;
    const R cshrink;
    const R cshrink1;
    static const RdHat G;

    RdHat Shrink(const RdHat &P) const { return (P - G) * cshrink + G; }

    RdHat Shrink1(const RdHat &P) const { return (P - G) * cshrink1 + G; }

    const int k;
    struct A4 {
      int dfon[4];

      A4(int k) {
        int ndf = (dHat == 3) ? ((k + 3) * (k + 2) * (k + 1) / 6)
                           : ((dHat == 2) ? ((k + 2) * (k + 1) / 2) : k + 1);

        dfon[0] = dfon[1] = dfon[2] = dfon[3] = 0;
        dfon[dHat] = ndf;

        if (verbosity > 9) {
          cout << "A4 " << k << " " << dfon[0] << dfon[1] << dfon[2] << dfon[3] << endl;
        }
      }

      operator const int *( ) const { return dfon; }
    };
    RdHat *Pt;
    TypeOfFE_LagrangeDC3d(int kk, R cc)
      :    // dfon ,N,nsub(graphique) ,  const mat interpolation , discontinuous
        GTypeOfFE< MMesh >(A4(kk), 1, Max(kk, 1), true, true), cshrink(cc), cshrink1(1. / cc),
        k(kk) {
      int n = this->NbDoF;

      if (verbosity > 9) {
        cout << "\n +++ Pdc" << k << " : ndof : " << n << endl;
      }

      SetPtPkDC<MMesh>(this->PtInterpolation, k, this->NbDoF, cc);
      if (verbosity > 9) 
         cout << " ppppp " << this->PtInterpolation << endl;
      

      {
        for (int i = 0; i < n; i++) {
          this->pInterpolation[i] = i;
          this->cInterpolation[i] = 0;
          this->dofInterpolation[i] = i;
          this->coefInterpolation[i] = 1.;
        }
      }
    }

    ~TypeOfFE_LagrangeDC3d( ) {}


    void FB(const What_d whatd, const MMesh &Th, const Element &K, const RdHat &PHat,RNMK_ &val) const;
    R operator( )(const FElement &K, const RdHat &PHat, const KN_< R > &u, int componante,int op) const;
	
   private:
    TypeOfFE_LagrangeDC3d(const TypeOfFE_LagrangeDC3d &);
    void operator=(const TypeOfFE_LagrangeDC3d &);
  };
  
  
 
  template<>
  void TypeOfFE_LagrangeDC3d<Mesh3>::FB(const What_d whatd, const Mesh3 &Th, const Mesh3::Element &K,const RdHat &PHat, RNMK_ &val) const {
    
	R3 P = this->Shrink1(PHat);
    R l[] = {1. - P.sum( ), P.x, P.y, P.z};

    assert(val.N( ) >= Element::nv);
    assert(val.M( ) == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {
      f0[0] = l[0];
      f0[1] = l[1];
      f0[2] = l[2];
      f0[3] = l[3];
    }

    if (whatd & Fop_D1) {
      R3 Dl[4];
      K.Gradlambda(Dl);

      for (int i = 0; i < 4; ++i) {
        Dl[i] *= cshrink1;
      }

      if (whatd & Fop_dx) {
        RN_ f0x(val('.', 0, op_dx));
        f0x[0] = Dl[0].x;
        f0x[1] = Dl[1].x;
        f0x[2] = Dl[2].x;
        f0x[3] = Dl[3].x;
      }

      if (whatd & Fop_dy) {
        RN_ f0y(val('.', 0, op_dy));
        f0y[0] = Dl[0].y;
        f0y[1] = Dl[1].y;
        f0y[2] = Dl[2].y;
        f0y[3] = Dl[3].y;
      }

      if (whatd & Fop_dz) {
        RN_ f0z(val('.', 0, op_dz));
        f0z[0] = Dl[0].z;
        f0z[1] = Dl[1].z;
        f0z[2] = Dl[2].z;
        f0z[3] = Dl[3].z;
      }
    }
  }


	  
  template<>
  R TypeOfFE_LagrangeDC3d<Mesh3>::operator( )(const FElement &K, const R3 &PHat1, const KN_< R > &u,int componante, int op) const {
    R3 PHat = Shrink1(PHat1);
    R r = 0;

    if (k == 1) {
      R u0(u(K(0))), u1(u(K(1))), u2(u(K(2))), u3(u(K(3)));
      if (op == 0) {
        R l[4];
        PHat.toBary(l);
        r = u0 * l[0] + u1 * l[1] + l[2] * u2 + l[3] * u3;
      } else if (op == op_dx || op == op_dy || op == op_dz) {
        const Element &T = K.T;
        R3 D[4];
        T.Gradlambda(D);

        for (int i = 0; i < 4; ++i) {
          D[i] *= cshrink1;
        }

        if (op == op_dx) {
          r = D[0].x * u0 + D[1].x * u1 + D[2].x * u2 + D[3].x * u3;
        } else if (op == op_dy) {
          r = D[0].y * u0 + D[1].y * u1 + D[2].y * u2 + D[3].y * u3;
        } else {
          r = D[0].z * u0 + D[1].z * u1 + D[2].z * u2 + D[3].z * u3;
        }
      }
    } else {
      ffassert(0);    // to do ..
    }

    return r;
  }

  template<>
  const R3 TypeOfFE_LagrangeDC3d<Mesh3>::G(1. / 4., 1. / 4., 1. / 4.);    
  
  

  template<>
  void TypeOfFE_LagrangeDC3d<MeshS>::FB(const What_d whatd, const MeshS &Th, const MeshS::Element &K,const RdHat &PHat, RNMK_ &val) const {
    R2 P = this->Shrink1(PHat);
    R l[] = {1. - P.sum( ), P.x, P.y};

    assert(val.N( ) >= Element::nv);
    assert(val.M( ) == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {
      f0[0] = l[0];
      f0[1] = l[1];
      f0[2] = l[2];
    }

    if (whatd & Fop_D1) {
      R3 Dl[3];
      K.Gradlambda(Dl);

      for (int i = 0; i < 3; ++i) {
        Dl[i] *= cshrink1;
      }

      if (whatd & Fop_dx) {
        RN_ f0x(val('.', 0, op_dx));
        f0x[0] = Dl[0].x;
        f0x[1] = Dl[1].x;
        f0x[2] = Dl[2].x;
      }

      if (whatd & Fop_dy) {
        RN_ f0y(val('.', 0, op_dy));
        f0y[0] = Dl[0].y;
        f0y[1] = Dl[1].y;
        f0y[2] = Dl[2].y;
      }

      if (whatd & Fop_dz) {
        RN_ f0z(val('.', 0, op_dz));
        f0z[0] = Dl[0].z;
        f0z[1] = Dl[1].z;
        f0z[2] = Dl[2].z;
      }
    }
  }

  template<>
  R TypeOfFE_LagrangeDC3d<MeshS>::operator( )(const FElement &K, const R2 &PHat1, const KN_< R > &u,int componante, int op) const {
    R2 PHat = Shrink1(PHat1);
    R r = 0;

    if (k == 1) {
      R u0(u(K(0))), u1(u(K(1))), u2(u(K(2)));

      if (op == 0) {
        R l[3];
        PHat.toBary(l);
        r = u0 * l[0] + u1 * l[1] + l[2] * u2 ;
      } else if (op == op_dx || op == op_dy || op == op_dz) {
        const Element &T = K.T;
        R3 D[3];
        T.Gradlambda(D);

        for (int i = 0; i < 3; ++i) {
          D[i] *= cshrink1;
        }

        if (op == op_dx) {
          r = D[0].x * u0 + D[1].x * u1 + D[2].x * u2 ;
        } else if (op == op_dy) {
          r = D[0].y * u0 + D[1].y * u1 + D[2].y * u2 ;
        } else {
          r = D[0].z * u0 + D[1].z * u1 + D[2].z * u2 ;
        }
      }
    } else {
      ffassert(0);    // to do ..
    }

    return r;
  }
  template<>
  const R2 TypeOfFE_LagrangeDC3d<MeshS>::G(1. / 3., 1. / 3.);
  
  
  template<>
  void TypeOfFE_LagrangeDC3d<MeshL>::FB(const What_d whatd, const MeshL &Th, const MeshL::Element &K,const RdHat &PHat, RNMK_ &val) const {
    R1 P = this->Shrink1(PHat);
    R l[] = {1. - P.x, P.x};

    assert(val.N( ) >= Element::nv);
    assert(val.M( ) == 1);

    val = 0;
    RN_ f0(val('.', 0, op_id));

    if (whatd & Fop_D0) {
      f0[0] = l[0];
      f0[1] = l[1];
    }

    if (whatd & Fop_D1) {
      R3 Dl[2];
      K.Gradlambda(Dl);

      for (int i = 0; i < 2; ++i) {
        Dl[i] *= cshrink1;
      }

      if (whatd & Fop_dx) {
        RN_ f0x(val('.', 0, op_dx));
        f0x[0] = Dl[0].x;
        f0x[1] = Dl[1].x;
      }

      if (whatd & Fop_dy) {
        RN_ f0y(val('.', 0, op_dy));
        f0y[0] = Dl[0].y;
        f0y[1] = Dl[1].y;
      }

      if (whatd & Fop_dz) {
        RN_ f0z(val('.', 0, op_dz));
        f0z[0] = Dl[0].z;
        f0z[1] = Dl[1].z;
      }
    }
  }

  template<>
  R TypeOfFE_LagrangeDC3d<MeshL>::operator( )(const FElement &K, const R1 &PHat1, const KN_< R > &u,int componante, int op) const {
    R1 PHat = Shrink1(PHat1);
    R r = 0;

    if (k == 1) {
      R u0(u(K(0))), u1(u(K(1)));

      if (op == 0) {
        R l[2];
        PHat.toBary(l);
        r = u0 * l[0] + u1 * l[1] ;
      } else if (op == op_dx || op == op_dy || op == op_dz) {
        const Element &T = K.T;
        R3 D[2];
        T.Gradlambda(D);

        for (int i = 0; i < 2; ++i) {
          D[i] *= cshrink1;
        }

        if (op == op_dx) {
          r = D[0].x * u0 + D[1].x * u1 ;
        } else if (op == op_dy) {
          r = D[0].y * u0 + D[1].y * u1 ;
        } else {
          r = D[0].z * u0 + D[1].z * u1 ;
        }
      }
    } else {
      ffassert(0);    // to do ..
    }

    return r;
  }
  
  template<>
  const R1 TypeOfFE_LagrangeDC3d<MeshL>::G(1. / 2.);
  
  
}    // namespace Fem2D

static void finit( ) {
  // link with FreeFem++
  static TypeOfFE_P1ttdc1_ P1dc1LagrangeP1dc1;
  static TypeOfFE_P2ttdc1_ P2dc1LagrangeP2dc1;
  
  static TypeOfFE_LagrangeDC3d<Mesh3> TypeOfFE_LagrangeDC3dtt(1, 0.999);
  static TypeOfFE_LagrangeDC3d<Mesh3> TypeOfFE_LagrangeDC3dtt1(1, 1.);
  
  static TypeOfFE_LagrangeDC3d<MeshS> TypeOfFE_LagrangeDC3dStt(1, 0.999);
  static TypeOfFE_LagrangeDC3d<MeshS> TypeOfFE_LagrangeDC3dStt1(1, 1.);
 
  static TypeOfFE_LagrangeDC3d<MeshL> TypeOfFE_LagrangeDC3dLtt(1, 0.999);
  static TypeOfFE_LagrangeDC3d<MeshL> TypeOfFE_LagrangeDC3dLtt1(1, 1.);
  
  // a static variable to add the finite element to freefem++
  static AddNewFE P1dcLagrange("P1dc1", &P1dc1LagrangeP1dc1);
  static AddNewFE P2dcLagrange("P2dc1", &P2dc1LagrangeP2dc1);
  
  static AddNewFE3 P1dttLagrange3d("P1dc3d2", &TypeOfFE_LagrangeDC3dtt);// deja code 
  static AddNewFE3 P1dttLagrange3d1("P1dc3d3", &TypeOfFE_LagrangeDC3dtt1);
  
  static AddNewFES P1dttLagrange3dS("P1dcS2", &TypeOfFE_LagrangeDC3dStt);
  static AddNewFES P1dttLagrange3dS1("P1dcS3", &TypeOfFE_LagrangeDC3dStt1);
  
  static AddNewFEL P1dttLagrange3dL("P1dcL2", &TypeOfFE_LagrangeDC3dLtt);
  static AddNewFEL P1dttLagrange3dL1("P1dcL3", &TypeOfFE_LagrangeDC3dLtt1);
   
}

LOADFUNC(finit)

// --- fin --
