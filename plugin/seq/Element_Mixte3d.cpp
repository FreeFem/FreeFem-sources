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
// SUMMARY : Add 3D finite elements
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

// RT1
// Pk = P1^3  + P1h (x,y,z)  ( dim Pk = 4*3 + 3) = 15

// Related files:
// LaplaceRT1.edp
// lame-TD-NSS.edp
// test-ElementMixte.edp

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

#include "ff++.hpp"
#include "AddNewFE.h"

#ifdef __INTEL_COMPILER
#pragma optimize("", off)
#endif

namespace Fem2D {
  // Author: Marcella Bonazzoli
  // Edge elements of degree 2, 3D
  // (they are called Edge13d because the Nedelec elements of degree 1 are called Edge03d)
  // TypeOfFE_Edge1_3d derived from GTypeOfFE<> (which is defined in FESpacen.hpp)
  class TypeOfFE_Edge1_3d : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    typedef GFElement< Mesh3 > FElement;

    static int dfon[];
    static const int d = Mesh::Rd::d;
    static const GQuadratureFormular< R1 > QFe;    // quadrature formula on an edge
    static const GQuadratureFormular< R2 > QFf;    // quadrature formula on a face
    TypeOfFE_Edge1_3d( );                          // constructor
    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
            RNMK_ &val) const;
    void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
             int *nump) const;
  };

  int TypeOfFE_Edge1_3d::dfon[] = {0, 2, 2, 0};    // 2 dofs on each edge, 2 dofs on each face

  // Quadrature formula on an edge, exact for degree 3 (ok: int_e (deg2*t *lambda))
  const GQuadratureFormular< R1 > TypeOfFE_Edge1_3d::QFe(-1 + 2 * 2, 2, GaussLegendre(2), true);
  // arguments: exact, num of integration pts, integration pts, clean (~GQuadratureFormular()
  // {if(clean) delete [] p;}) GaussLegendre defined in QuadratureFormular.cpp

  // Quadrature formula on a face, exact for degree 2 (ok: int_f (deg2*t)), internal quadrature
  // points
  static GQuadraturePoint< R2 > P_QuadratureFormular_T_2_intp[3] = {
    GQuadraturePoint< R2 >(1. / 3., R2(1. / 6., 4. / 6.)),
    GQuadraturePoint< R2 >(1. / 3., R2(4. / 6., 1. / 6.)),
    GQuadraturePoint< R2 >(1. / 3., R2(1. / 6., 1. / 6.))};
  const GQuadratureFormular< R2 > TypeOfFE_Edge1_3d::QFf(2, 3, P_QuadratureFormular_T_2_intp);

  // In Mesh3dn.cpp:
  // static const int  nvedgeTet[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
  // static const int  nvfaceTet[4][3] = { {3,2,1},{0,2,3},{3,1,0},{0,1,2} };
  // In GenericMesh.hpp:
  // Vertex& at(int i) {return *vertices[i];}
  // In GenericMesh.hpp:
  // Rd Edge(int i) const {ASSERTION(i>=0 && i <ne);
  // return Rd(at(nvedge[i][0]),at(nvedge[i][1]));}
  // In GenericMesh.hpp:
  // int   EdgeOrientation(int i) const  +/- 1
  // { return &at(nvedge[i][0]) < &at(nvedge[i][1]);}

  // Constructor
  // dfon, d, nsub,
  // kPi = NbcoefforInterpolation is the number of alphas in (13.1) (chapter 13 ff++doc)
  // = 3(=numComp)*QFe.n(num quad pts edge)*2*ne(2 fncts per edge) +
  // 3(=numComp)*QFf.n(num quad pts face)*2*nf(2 fncts per face)
  // npPi = NbPtforInterpolation = ne*QFe.n+nf*QFf.n
  // invariantinterpolationMatrix = false (i.e. it depends on the tetrahedron), discon=true
  TypeOfFE_Edge1_3d::TypeOfFE_Edge1_3d( )
    : GTypeOfFE< Mesh3 >(TypeOfFE_Edge1_3d::dfon, d, 1,
                         Element::ne * 2 * 3 * QFe.n + Element::nf * 2 * 3 * QFf.n,
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

      for (e = 0; e < (Element::ne)*2; e++) {    // loop on the 12 dofs
        if (e % 2 == 1) {
          p = p - QFe.n;
        }    // if I consider an 'even' dof, the quad pts are the ones of the previous dof (they
             // correspond to the same edge)

        for (int q = 0; q < QFe.n; ++q, ++p) {    // loop on the 2 edge quadrature pts
          for (int c = 0; c < 3; c++, i++) {      // loop on the 3 components
            this->pInterpolation[i] = p;          // pk in (13.1)
            this->cInterpolation[i] = c;          // jk in (13.1)
            this->dofInterpolation[i] = e;        // ik in (13.1)
            this->coefInterpolation[i] = 0.;      // alfak: we will fill them with 'set' (below)
                                                  // because they depend on the tetrahedron
          }
        }
      }

      // We build the indices in (13.1) : face dofs
      // (the indices i and p mustn't be reinitialised)
      for (int f = 0; f < (Element::nf)*2; f++) {    // loop on the 8 dofs
        if (f % 2 == 1) {
          p = p - QFf.n;
        }    // if I consider an 'even' dof, the quad pts are the ones of the previous dof (they
             // correspond to the same face)

        for (int q = 0; q < QFf.n; ++q, ++p) {    // loop on the 3 face quadrature pts
          for (int c = 0; c < 3; c++, i++) {      // loop on the 3 components
            this->pInterpolation[i] = p;          // pk in (13.1)
            this->cInterpolation[i] = c;          // jk in (13.1)
            this->dofInterpolation[i] = e + f;    // ik in (13.1)
            this->coefInterpolation[i] = 0.;      // alphak: we will fill them with 'set' (below)
                                                  // because they depend on the tetrahedron
          }
        }
      }
    }
  }

  // For the coefficients of interpolation alphak in (13.1)
  void TypeOfFE_Edge1_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
                              int ocoef, int odf, int *nump) const {
    int i = ocoef, p = 0;

    // 12 edge dofs
    int e;

    for (e = 0; e < (Element::ne)*2; e++) {
      int ee = e / 2;
      int eo = K.EdgeOrientation(ee) > 0;    //  change FH; jan 2019
      R3 E = K.Edge(ee);    // the edge local number is given by the integer division betweeand
      if (!eo) {
        E = -E;
      }

      if (e % 2 == 1) {
        p = p - QFe.n;
      }    // if I consider an 'even' dof, the quad pts are the ones of the previous dof (they
           // correspond to the same edge)

      for (int q = 0; q < QFe.n; ++q, ++p) {
        double ll = QFe[q].x;    // value of lambda_0 or lambda_1
        if ((e % 2 + eo) == 1) {
          ll = 1 - ll;    // if exactly one between e%2 and eo is equal to 1 (so the sum is equal to
                          // 1), take the other lambda
        }

        // i.e. if I'm considering the 2nd dof of the edge (or) if the edge is badly oriented, take
        // the other lambda (but not if both)
        for (int c = 0; c < 3; c++, i++) {
          M.coef[i] = E[c] * QFe[q].a * ll;
          // QFe[q].a is the weight of the integration point q
          // QFe[q].x is the x over [0,1] of the integration point q
        }
      }
    }

    // (the indices i and p mustn't be reinitialised)
    // 8 face dofs
    for (int f = 0; f < (Element::nf)*2; f++) {    // loop on the 8 face dofs
      int ff = f / 2;                              // the face number
      int iff = f % 2;                             // 1st or 2nd dof of the face
      const Element::Vertex *fV[3] = {&K.at(Element::nvface[ff][0]), &K.at(Element::nvface[ff][1]),
                                      &K.at(Element::nvface[ff][2])};
      // We 'order' the 3 vertices of the face according to their global numbers:
      // i0 will be the local number in the FACE of the vertex with the smallest global number
      // i1 will be the local number in the face of the vertex with the second smallest global
      // number i2 will be the local number in the face of the vertex with the largest global number
      int i0 = 0, i1 = 1, i2 = 2;
      if (fV[i0] > fV[i1]) {
        Exchange(i0, i1);
      }

      if (fV[i1] > fV[i2]) {
        Exchange(i1, i2);
        if (fV[i0] > fV[i1]) {
          Exchange(i0, i1);
        }
      }

      // now local numbers in the tetrahedron:
      i0 = Element::nvface[ff][i0], i1 = Element::nvface[ff][i1], i2 = Element::nvface[ff][i2];
      int ie0 = i0,
          ie1 = iff == 0 ? i1 : i2;    // edge for the face dof (its endpoints local numbers)
      R3 E(K[ie0], K[ie1]);
      if (iff) {
        p = p - QFf.n;
      }    // if I consider an 'even' dof, the quad pts are the ones of the previous dof (they
           // correspond to the same face)

      for (int q = 0; q < QFf.n; ++q, ++p) {    // loop on the 3 face quadrature pts
        for (int c = 0; c < 3; c++, i++) {      // loop on the 3 components
          M.coef[i] = E[c] * QFf[q].a;
        }
      }
    }
  }

  // In Mesh3dn.hpp:
  // void Gradlambda(R3 * GradL) const
  // {
  // R3 V1(at(0),at(1));
  // R3 V2(at(0),at(2));
  // R3 V3(at(0),at(3));
  // R det1=1./(6.*mesure());
  // GradL[1]= (V2^V3)*det1;
  // GradL[2]= (V3^V1)*det1;
  // GradL[3]= (V1^V2)*det1;
  // GradL[0]=-GradL[1]-GradL[2]-GradL[3];
  // }

  /*
   * FB: We are on a tetrahedron K
   * (the local numbering of its 4 vertices is given by how they are listed where K is described in
   * the mesh file).
   *
   * AT FIRST we BUILD the 'pre-basis' functions OMEGA in the order given the GLOBAL numbering of
   * the tetrahedron vertices, that is the 1st examined edge is from the vertex with the 1st
   * smallest GLOBAL number to the one the 2nd smallest GLOBAL number, and so on (see nvedege), and
   * the 1st examined face is the one opposite the vertex with the smallest GLOBAL number (and 2
   * edges of the face are chosen and oriented looking at its 3 global numbers), and so on. In this
   * way:
   * - a dof which is common to 2 adjacent tetrahedra is the same in the two tetrahedra
   * (since the orientation of each edge is chosen using the global numbers),
   * - we can use the coefficients giving a dual basis that were calculated for the reference
   * tetrahedron since the structure of orientation of the edges is the same (up to a rotation) as
   * the one used in the reference tetrahedron.
   *
   * BUT THEN, when we build the DUAL basis functions PHI, we use the permutation p20 to go back to
   * the FreeFem numbering of the dofs (which follows the LOCAL numbering). For instance the 1st
   * examined edge can be the 3rd edge looking at the local numbering. Here with 'dual' I mean basis
   * functions and dofs in duality.
   */

  // val contains the values of the basis functions and of their derivatives at the point of K
  // corresponding to the point P of the reference tetrahedron, by components
  void TypeOfFE_Edge1_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                             const RdHat &PHat, RNMK_ &val) const {
    assert(val.N( ) >= 20);    // 20 degrees of freedom
    assert(val.M( ) == 3);     // 3 components
    // -------------
    // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL
    // number (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
    // perm[3] is the local number of the vertex with the biggest global number.)
    const Element::Vertex *tV[4] = {&K.at(0), &K.at(1), &K.at(2), &K.at(3)};
    int k0 = 0, k1 = 1, k2 = 2, k3 = 3;
    if (tV[k0] > tV[k1]) {
      Exchange(k0, k1);
    }

    if (tV[k1] > tV[k2]) {
      Exchange(k1, k2);
    }

    if (tV[k2] > tV[k3]) {
      Exchange(k2, k3);
    }

    if (tV[k0] > tV[k1]) {
      Exchange(k0, k1);
    }

    if (tV[k1] > tV[k2]) {
      Exchange(k1, k2);
    }

    if (tV[k0] > tV[k1]) {
      Exchange(k0, k1);
    }

    int perm[4] = {k0, k1, k2, k3};
    // -------------
    // We build mynvface to be used here instead of the FreeFem nvface,
    // (in order to exploit the results of perm and write a better code),
    // in mynvface in all the triplets the numbers are increasing
    static const int mynvface[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
    // -------------
    // If [a,b] is the i-th edge (where a,b are its vertices local numbers), edgesMap[(a+1)*(b+1)] =
    // i
    int edgesMap[13] = {-1, -1, 0,  1,  2,  -1, 3,
                        -1, 4,  -1, -1, -1, 5};    // a map<int,int> would be more slow
    // -------------
    // the 4 barycentric coordinates for the reference tetrahedron evaluated at the point P
    // (they have the same value at the real tetrahedron's point corresponding to the reference
    // tetrahedron's point P)
    R l[] = {1. - PHat.sum( ), PHat.x, PHat.y, PHat.z};
    R3 D[4];
    K.Gradlambda(D);    // (riempie un array di 4 R3)
    val = 0;
    // -----
    int p20[20];    // the permutation from the dofs numbering of my tetrahedron (numbered using the
                    // GLOBAL vertex numbers) to the dofs numbering of FreeFem !!!!!

    for (int i = 0; i < 6; ++i) {    // edges
      // see below
      int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
      int i0 = perm[ii0];
      int i1 = perm[ii1];
      int iEdge = edgesMap[(i0 + 1) * (i1 + 1)];    // i of the edge [i0,i1]
      p20[i * 2] = iEdge * 2;
      p20[i * 2 + 1] = iEdge * 2 + 1;
    }

    for (int j = 0; j < 4; ++j) {    // faces
      // using perm the 1st examined face is the one opposite the vertex with the smallest GLOBAL
      // number, and so on (see below)
      int jFace = perm[j];
      p20[12 + j * 2] = 12 + jFace * 2;
      p20[12 + j * 2 + 1] = 12 + jFace * 2 + 1;
    }

    // -----

    if (whatd & Fop_D0) {    // Fop_D0 defined in FESpacen.hpp
      R3 X = K(PHat);
      // First, the functions omega (they don't constitute a dual basis! only a basis)
      R3 omega[20];

      // 12 edge functions:
      for (int i = 0; i < 6; ++i) {
        int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        // ! :
        // using perm, [i0,i1] is already from the smallest global number to the greatest global
        // number, since nvedge always gives indices which read perm from left to right
        omega[i * 2] = l[i0] * (l[i0] * D[i1] - l[i1] * D[i0]);
        omega[i * 2 + 1] = l[i1] * (l[i0] * D[i1] - l[i1] * D[i0]);
      }

      // 8 face functions:
      for (int j = 0; j < 4; ++j) {
        // In (my)nvface the face opposite the vertex k is (my)nvface[k][],
        // using perm the 1st examined face is the one opposite the vertex with the smallest GLOBAL
        // number, and so on, and, since in mynvface the numbers always increase, [i0,i1,i2] are
        // already ordered with increasing global number
        int ii0 = mynvface[j][0];
        int ii1 = mynvface[j][1];
        int ii2 = mynvface[j][2];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        int i2 = perm[ii2];
        omega[12 + j * 2] = l[i2] * (l[i0] * D[i1] - l[i1] * D[i0]);
        omega[12 + j * 2 + 1] = l[i1] * (l[i0] * D[i2] - l[i2] * D[i0]);
      }

      // Now, the functions phi that really constitute a dual basis
      R3 phi[20];
      phi[p20[0]] = +4 * omega[0] - 2 * omega[1] - 4 * omega[16] + 2 * omega[17] - 4 * omega[18] +
                    2 * omega[19];
      phi[p20[1]] = -2 * omega[0] + 4 * omega[1] - 2 * omega[16] - 2 * omega[17] - 2 * omega[18] -
                    2 * omega[19];
      phi[p20[2]] = +4 * omega[2] - 2 * omega[3] - 4 * omega[14] + 2 * omega[15] + 2 * omega[18] -
                    4 * omega[19];
      phi[p20[3]] = -2 * omega[2] + 4 * omega[3] - 2 * omega[14] - 2 * omega[15] - 2 * omega[18] -
                    2 * omega[19];
      phi[p20[4]] = +4 * omega[4] - 2 * omega[5] + 2 * omega[14] - 4 * omega[15] + 2 * omega[16] -
                    4 * omega[17];
      phi[p20[5]] = -2 * omega[4] + 4 * omega[5] - 2 * omega[14] - 2 * omega[15] - 2 * omega[16] -
                    2 * omega[17];
      phi[p20[6]] = +4 * omega[6] - 2 * omega[7] - 4 * omega[12] + 2 * omega[13] + 2 * omega[18] -
                    4 * omega[19];
      phi[p20[7]] = -2 * omega[6] + 4 * omega[7] - 2 * omega[12] - 2 * omega[13] + 4 * omega[18] -
                    2 * omega[19];
      phi[p20[8]] = +4 * omega[8] - 2 * omega[9] + 2 * omega[12] - 4 * omega[13] + 2 * omega[16] -
                    4 * omega[17];
      phi[p20[9]] = -2 * omega[8] + 4 * omega[9] - 2 * omega[12] - 2 * omega[13] + 4 * omega[16] -
                    2 * omega[17];
      phi[p20[10]] = +4 * omega[10] - 2 * omega[11] + 2 * omega[12] - 4 * omega[13] +
                     2 * omega[14] - 4 * omega[15];
      phi[p20[11]] = -2 * omega[10] + 4 * omega[11] + 4 * omega[12] - 2 * omega[13] +
                     4 * omega[14] - 2 * omega[15];
      phi[p20[12]] = +8 * omega[12] - 4 * omega[13];
      phi[p20[13]] = -4 * omega[12] + 8 * omega[13];
      phi[p20[14]] = +8 * omega[14] - 4 * omega[15];
      phi[p20[15]] = -4 * omega[14] + 8 * omega[15];
      phi[p20[16]] = +8 * omega[16] - 4 * omega[17];
      phi[p20[17]] = -4 * omega[16] + 8 * omega[17];
      phi[p20[18]] = +8 * omega[18] - 4 * omega[19];
      phi[p20[19]] = -4 * omega[18] + 8 * omega[19];

      for (int k = 0; k < 20; ++k) {
        val(k, 0, op_id) = phi[k].x;
        val(k, 1, op_id) = phi[k].y;
        val(k, 2, op_id) = phi[k].z;
      }
    }

    if (whatd & Fop_D1) {    // First Derivatives wrt x,y,z
      R3 omegadx[20];
      R3 omegady[20];
      R3 omegadz[20];

      // 12 edge functions:
      for (int i = 0; i < 6; ++i) {
        int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        // using perm, [i0,i1] is already from the smallest global number to the greatest global
        // number, since nvedge always gives indices which read perm from left to right
        if (whatd & Fop_dx) {
          omegadx[i * 2] =
            D[i0].x * (l[i0] * D[i1] - l[i1] * D[i0]) + l[i0] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omegadx[i * 2 + 1] =
            D[i1].x * (l[i0] * D[i1] - l[i1] * D[i0]) + l[i1] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
        }

        if (whatd & Fop_dy) {
          omegady[i * 2] =
            D[i0].y * (l[i0] * D[i1] - l[i1] * D[i0]) + l[i0] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omegady[i * 2 + 1] =
            D[i1].y * (l[i0] * D[i1] - l[i1] * D[i0]) + l[i1] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
        }

        if (whatd & Fop_dz) {
          omegadz[i * 2] =
            D[i0].z * (l[i0] * D[i1] - l[i1] * D[i0]) + l[i0] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omegadz[i * 2 + 1] =
            D[i1].z * (l[i0] * D[i1] - l[i1] * D[i0]) + l[i1] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
        }
      }

      // 8 face functions:
      for (int j = 0; j < 4; ++j) {
        // In (my)nvface the face opposite the vertex k is (my)nvface[k][],
        // using perm the 1st examined face is the one opposite the vertex with the smallest GLOBAL
        // number, and so on, and, since in mynvface the numbers always increase, [i0,i1,i2] are
        // already ordered with increasing global number
        int ii0 = mynvface[j][0];
        int ii1 = mynvface[j][1];
        int ii2 = mynvface[j][2];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        int i2 = perm[ii2];
        if (whatd & Fop_dx) {
          omegadx[12 + j * 2] =
            D[i2].x * (l[i0] * D[i1] - l[i1] * D[i0]) + l[i2] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omegadx[12 + j * 2 + 1] =
            D[i1].x * (l[i0] * D[i2] - l[i2] * D[i0]) + l[i1] * (D[i0].x * D[i2] - D[i2].x * D[i0]);
        }

        if (whatd & Fop_dy) {
          omegady[12 + j * 2] =
            D[i2].y * (l[i0] * D[i1] - l[i1] * D[i0]) + l[i2] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omegady[12 + j * 2 + 1] =
            D[i1].y * (l[i0] * D[i2] - l[i2] * D[i0]) + l[i1] * (D[i0].y * D[i2] - D[i2].y * D[i0]);
        }

        if (whatd & Fop_dz) {
          omegadz[12 + j * 2] =
            D[i2].z * (l[i0] * D[i1] - l[i1] * D[i0]) + l[i2] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omegadz[12 + j * 2 + 1] =
            D[i1].z * (l[i0] * D[i2] - l[i2] * D[i0]) + l[i1] * (D[i0].z * D[i2] - D[i2].z * D[i0]);
        }
      }

      R3 phidx[20];
      if (whatd & Fop_dx) {
        phidx[p20[0]] = +4 * omegadx[0] - 2 * omegadx[1] - 4 * omegadx[16] + 2 * omegadx[17] -
                        4 * omegadx[18] + 2 * omegadx[19];
        phidx[p20[1]] = -2 * omegadx[0] + 4 * omegadx[1] - 2 * omegadx[16] - 2 * omegadx[17] -
                        2 * omegadx[18] - 2 * omegadx[19];
        phidx[p20[2]] = +4 * omegadx[2] - 2 * omegadx[3] - 4 * omegadx[14] + 2 * omegadx[15] +
                        2 * omegadx[18] - 4 * omegadx[19];
        phidx[p20[3]] = -2 * omegadx[2] + 4 * omegadx[3] - 2 * omegadx[14] - 2 * omegadx[15] -
                        2 * omegadx[18] - 2 * omegadx[19];
        phidx[p20[4]] = +4 * omegadx[4] - 2 * omegadx[5] + 2 * omegadx[14] - 4 * omegadx[15] +
                        2 * omegadx[16] - 4 * omegadx[17];
        phidx[p20[5]] = -2 * omegadx[4] + 4 * omegadx[5] - 2 * omegadx[14] - 2 * omegadx[15] -
                        2 * omegadx[16] - 2 * omegadx[17];
        phidx[p20[6]] = +4 * omegadx[6] - 2 * omegadx[7] - 4 * omegadx[12] + 2 * omegadx[13] +
                        2 * omegadx[18] - 4 * omegadx[19];
        phidx[p20[7]] = -2 * omegadx[6] + 4 * omegadx[7] - 2 * omegadx[12] - 2 * omegadx[13] +
                        4 * omegadx[18] - 2 * omegadx[19];
        phidx[p20[8]] = +4 * omegadx[8] - 2 * omegadx[9] + 2 * omegadx[12] - 4 * omegadx[13] +
                        2 * omegadx[16] - 4 * omegadx[17];
        phidx[p20[9]] = -2 * omegadx[8] + 4 * omegadx[9] - 2 * omegadx[12] - 2 * omegadx[13] +
                        4 * omegadx[16] - 2 * omegadx[17];
        phidx[p20[10]] = +4 * omegadx[10] - 2 * omegadx[11] + 2 * omegadx[12] - 4 * omegadx[13] +
                         2 * omegadx[14] - 4 * omegadx[15];
        phidx[p20[11]] = -2 * omegadx[10] + 4 * omegadx[11] + 4 * omegadx[12] - 2 * omegadx[13] +
                         4 * omegadx[14] - 2 * omegadx[15];
        phidx[p20[12]] = +8 * omegadx[12] - 4 * omegadx[13];
        phidx[p20[13]] = -4 * omegadx[12] + 8 * omegadx[13];
        phidx[p20[14]] = +8 * omegadx[14] - 4 * omegadx[15];
        phidx[p20[15]] = -4 * omegadx[14] + 8 * omegadx[15];
        phidx[p20[16]] = +8 * omegadx[16] - 4 * omegadx[17];
        phidx[p20[17]] = -4 * omegadx[16] + 8 * omegadx[17];
        phidx[p20[18]] = +8 * omegadx[18] - 4 * omegadx[19];
        phidx[p20[19]] = -4 * omegadx[18] + 8 * omegadx[19];

        for (int k = 0; k < 20; ++k) {
          val(k, 0, op_dx) = phidx[k].x;
          val(k, 1, op_dx) = phidx[k].y;
          val(k, 2, op_dx) = phidx[k].z;
        }
      }

      R3 phidy[20];
      if (whatd & Fop_dy) {
        phidy[p20[0]] = +4 * omegady[0] - 2 * omegady[1] - 4 * omegady[16] + 2 * omegady[17] -
                        4 * omegady[18] + 2 * omegady[19];
        phidy[p20[1]] = -2 * omegady[0] + 4 * omegady[1] - 2 * omegady[16] - 2 * omegady[17] -
                        2 * omegady[18] - 2 * omegady[19];
        phidy[p20[2]] = +4 * omegady[2] - 2 * omegady[3] - 4 * omegady[14] + 2 * omegady[15] +
                        2 * omegady[18] - 4 * omegady[19];
        phidy[p20[3]] = -2 * omegady[2] + 4 * omegady[3] - 2 * omegady[14] - 2 * omegady[15] -
                        2 * omegady[18] - 2 * omegady[19];
        phidy[p20[4]] = +4 * omegady[4] - 2 * omegady[5] + 2 * omegady[14] - 4 * omegady[15] +
                        2 * omegady[16] - 4 * omegady[17];
        phidy[p20[5]] = -2 * omegady[4] + 4 * omegady[5] - 2 * omegady[14] - 2 * omegady[15] -
                        2 * omegady[16] - 2 * omegady[17];
        phidy[p20[6]] = +4 * omegady[6] - 2 * omegady[7] - 4 * omegady[12] + 2 * omegady[13] +
                        2 * omegady[18] - 4 * omegady[19];
        phidy[p20[7]] = -2 * omegady[6] + 4 * omegady[7] - 2 * omegady[12] - 2 * omegady[13] +
                        4 * omegady[18] - 2 * omegady[19];
        phidy[p20[8]] = +4 * omegady[8] - 2 * omegady[9] + 2 * omegady[12] - 4 * omegady[13] +
                        2 * omegady[16] - 4 * omegady[17];
        phidy[p20[9]] = -2 * omegady[8] + 4 * omegady[9] - 2 * omegady[12] - 2 * omegady[13] +
                        4 * omegady[16] - 2 * omegady[17];
        phidy[p20[10]] = +4 * omegady[10] - 2 * omegady[11] + 2 * omegady[12] - 4 * omegady[13] +
                         2 * omegady[14] - 4 * omegady[15];
        phidy[p20[11]] = -2 * omegady[10] + 4 * omegady[11] + 4 * omegady[12] - 2 * omegady[13] +
                         4 * omegady[14] - 2 * omegady[15];
        phidy[p20[12]] = +8 * omegady[12] - 4 * omegady[13];
        phidy[p20[13]] = -4 * omegady[12] + 8 * omegady[13];
        phidy[p20[14]] = +8 * omegady[14] - 4 * omegady[15];
        phidy[p20[15]] = -4 * omegady[14] + 8 * omegady[15];
        phidy[p20[16]] = +8 * omegady[16] - 4 * omegady[17];
        phidy[p20[17]] = -4 * omegady[16] + 8 * omegady[17];
        phidy[p20[18]] = +8 * omegady[18] - 4 * omegady[19];
        phidy[p20[19]] = -4 * omegady[18] + 8 * omegady[19];

        for (int k = 0; k < 20; ++k) {
          val(k, 0, op_dy) = phidy[k].x;
          val(k, 1, op_dy) = phidy[k].y;
          val(k, 2, op_dy) = phidy[k].z;
        }
      }

      R3 phidz[20];
      if (whatd & Fop_dz) {
        phidz[p20[0]] = +4 * omegadz[0] - 2 * omegadz[1] - 4 * omegadz[16] + 2 * omegadz[17] -
                        4 * omegadz[18] + 2 * omegadz[19];
        phidz[p20[1]] = -2 * omegadz[0] + 4 * omegadz[1] - 2 * omegadz[16] - 2 * omegadz[17] -
                        2 * omegadz[18] - 2 * omegadz[19];
        phidz[p20[2]] = +4 * omegadz[2] - 2 * omegadz[3] - 4 * omegadz[14] + 2 * omegadz[15] +
                        2 * omegadz[18] - 4 * omegadz[19];
        phidz[p20[3]] = -2 * omegadz[2] + 4 * omegadz[3] - 2 * omegadz[14] - 2 * omegadz[15] -
                        2 * omegadz[18] - 2 * omegadz[19];
        phidz[p20[4]] = +4 * omegadz[4] - 2 * omegadz[5] + 2 * omegadz[14] - 4 * omegadz[15] +
                        2 * omegadz[16] - 4 * omegadz[17];
        phidz[p20[5]] = -2 * omegadz[4] + 4 * omegadz[5] - 2 * omegadz[14] - 2 * omegadz[15] -
                        2 * omegadz[16] - 2 * omegadz[17];
        phidz[p20[6]] = +4 * omegadz[6] - 2 * omegadz[7] - 4 * omegadz[12] + 2 * omegadz[13] +
                        2 * omegadz[18] - 4 * omegadz[19];
        phidz[p20[7]] = -2 * omegadz[6] + 4 * omegadz[7] - 2 * omegadz[12] - 2 * omegadz[13] +
                        4 * omegadz[18] - 2 * omegadz[19];
        phidz[p20[8]] = +4 * omegadz[8] - 2 * omegadz[9] + 2 * omegadz[12] - 4 * omegadz[13] +
                        2 * omegadz[16] - 4 * omegadz[17];
        phidz[p20[9]] = -2 * omegadz[8] + 4 * omegadz[9] - 2 * omegadz[12] - 2 * omegadz[13] +
                        4 * omegadz[16] - 2 * omegadz[17];
        phidz[p20[10]] = +4 * omegadz[10] - 2 * omegadz[11] + 2 * omegadz[12] - 4 * omegadz[13] +
                         2 * omegadz[14] - 4 * omegadz[15];
        phidz[p20[11]] = -2 * omegadz[10] + 4 * omegadz[11] + 4 * omegadz[12] - 2 * omegadz[13] +
                         4 * omegadz[14] - 2 * omegadz[15];
        phidz[p20[12]] = +8 * omegadz[12] - 4 * omegadz[13];
        phidz[p20[13]] = -4 * omegadz[12] + 8 * omegadz[13];
        phidz[p20[14]] = +8 * omegadz[14] - 4 * omegadz[15];
        phidz[p20[15]] = -4 * omegadz[14] + 8 * omegadz[15];
        phidz[p20[16]] = +8 * omegadz[16] - 4 * omegadz[17];
        phidz[p20[17]] = -4 * omegadz[16] + 8 * omegadz[17];
        phidz[p20[18]] = +8 * omegadz[18] - 4 * omegadz[19];
        phidz[p20[19]] = -4 * omegadz[18] + 8 * omegadz[19];

        for (int k = 0; k < 20; ++k) {
          val(k, 0, op_dz) = phidz[k].x;
          val(k, 1, op_dz) = phidz[k].y;
          val(k, 2, op_dz) = phidz[k].z;
        }
      }
    }

    if (whatd & Fop_D2) {    // Second Derivatives
      R3 omegadxx[20];
      R3 omegadyy[20];
      R3 omegadzz[20];
      R3 omegadxy[20];
      R3 omegadxz[20];
      R3 omegadyz[20];

      // 12 edge functions:
      for (int i = 0; i < 6; ++i) {
        int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        // using perm, [i0,i1] is already from the smallest global number to the greatest global
        // number, since nvedge always gives indices which read perm from left to right
        if (whatd & Fop_dxx) {
          omegadxx[i * 2] = D[i0].x * (D[i0].x * D[i1] - D[i1].x * D[i0]) +
                            D[i0].x * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omegadxx[i * 2 + 1] = D[i1].x * (D[i0].x * D[i1] - D[i1].x * D[i0]) +
                                D[i1].x * (D[i0].x * D[i1] - D[i1].x * D[i0]);
        }

        if (whatd & Fop_dyy) {
          omegadyy[i * 2] = D[i0].y * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
                            D[i0].y * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omegadyy[i * 2 + 1] = D[i1].y * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
                                D[i1].y * (D[i0].y * D[i1] - D[i1].y * D[i0]);
        }

        if (whatd & Fop_dzz) {
          omegadzz[i * 2] = D[i0].z * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                            D[i0].z * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omegadzz[i * 2 + 1] = D[i1].z * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                                D[i1].z * (D[i0].z * D[i1] - D[i1].z * D[i0]);
        }

        if (whatd & Fop_dxy) {
          omegadxy[i * 2] = D[i0].x * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
                            D[i0].y * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omegadxy[i * 2 + 1] = D[i1].x * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
                                D[i1].y * (D[i0].x * D[i1] - D[i1].x * D[i0]);
        }

        if (whatd & Fop_dxz) {
          omegadxz[i * 2] = D[i0].x * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                            D[i0].z * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omegadxz[i * 2 + 1] = D[i1].x * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                                D[i1].z * (D[i0].x * D[i1] - D[i1].x * D[i0]);
        }

        if (whatd & Fop_dyz) {
          omegadyz[i * 2] = D[i0].y * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                            D[i0].z * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omegadyz[i * 2 + 1] = D[i1].y * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                                D[i1].z * (D[i0].y * D[i1] - D[i1].y * D[i0]);
        }
      }

      // 8 face functions:
      for (int j = 0; j < 4; ++j) {
        // In (my)nvface the face opposite the vertex k is (my)nvface[k][],
        // using perm the 1st examined face is the one opposite the vertex with the smallest GLOBAL
        // number, and so on, and, since in mynvface the numbers always increase, [i0,i1,i2] are
        // already ordered with increasing global number
        int ii0 = mynvface[j][0];
        int ii1 = mynvface[j][1];
        int ii2 = mynvface[j][2];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        int i2 = perm[ii2];
        if (whatd & Fop_dxx) {
          omegadxx[12 + j * 2] = D[i2].x * (D[i0].x * D[i1] - D[i1].x * D[i0]) +
                                 D[i2].x * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omegadxx[12 + j * 2 + 1] = D[i1].x * (D[i0].x * D[i2] - D[i2].x * D[i0]) +
                                     D[i1].x * (D[i0].x * D[i2] - D[i2].x * D[i0]);
        }

        if (whatd & Fop_dyy) {
          omegadyy[12 + j * 2] = D[i2].y * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
                                 D[i2].y * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omegadyy[12 + j * 2 + 1] = D[i1].y * (D[i0].y * D[i2] - D[i2].y * D[i0]) +
                                     D[i1].y * (D[i0].y * D[i2] - D[i2].y * D[i0]);
        }

        if (whatd & Fop_dzz) {
          omegadzz[12 + j * 2] = D[i2].z * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                                 D[i2].z * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omegadzz[12 + j * 2 + 1] = D[i1].z * (D[i0].z * D[i2] - D[i2].z * D[i0]) +
                                     D[i1].z * (D[i0].z * D[i2] - D[i2].z * D[i0]);
        }

        if (whatd & Fop_dxy) {
          omegadxy[12 + j * 2] = D[i2].x * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
                                 D[i2].y * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omegadxy[12 + j * 2 + 1] = D[i1].x * (D[i0].y * D[i2] - D[i2].y * D[i0]) +
                                     D[i1].y * (D[i0].x * D[i2] - D[i2].x * D[i0]);
        }

        if (whatd & Fop_dxz) {
          omegadxz[12 + j * 2] = D[i2].x * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                                 D[i2].z * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omegadxz[12 + j * 2 + 1] = D[i1].x * (D[i0].z * D[i2] - D[i2].z * D[i0]) +
                                     D[i1].z * (D[i0].x * D[i2] - D[i2].x * D[i0]);
        }

        if (whatd & Fop_dyz) {
          omegadyz[12 + j * 2] = D[i2].y * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                                 D[i2].z * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omegadyz[12 + j * 2 + 1] = D[i1].y * (D[i0].z * D[i2] - D[i2].z * D[i0]) +
                                     D[i1].z * (D[i0].y * D[i2] - D[i2].y * D[i0]);
        }
      }

      R3 phidxx[20];
      if (whatd & Fop_dxx) {
        phidxx[p20[0]] = +4 * omegadxx[0] - 2 * omegadxx[1] - 4 * omegadxx[16] + 2 * omegadxx[17] -
                         4 * omegadxx[18] + 2 * omegadxx[19];
        phidxx[p20[1]] = -2 * omegadxx[0] + 4 * omegadxx[1] - 2 * omegadxx[16] - 2 * omegadxx[17] -
                         2 * omegadxx[18] - 2 * omegadxx[19];
        phidxx[p20[2]] = +4 * omegadxx[2] - 2 * omegadxx[3] - 4 * omegadxx[14] + 2 * omegadxx[15] +
                         2 * omegadxx[18] - 4 * omegadxx[19];
        phidxx[p20[3]] = -2 * omegadxx[2] + 4 * omegadxx[3] - 2 * omegadxx[14] - 2 * omegadxx[15] -
                         2 * omegadxx[18] - 2 * omegadxx[19];
        phidxx[p20[4]] = +4 * omegadxx[4] - 2 * omegadxx[5] + 2 * omegadxx[14] - 4 * omegadxx[15] +
                         2 * omegadxx[16] - 4 * omegadxx[17];
        phidxx[p20[5]] = -2 * omegadxx[4] + 4 * omegadxx[5] - 2 * omegadxx[14] - 2 * omegadxx[15] -
                         2 * omegadxx[16] - 2 * omegadxx[17];
        phidxx[p20[6]] = +4 * omegadxx[6] - 2 * omegadxx[7] - 4 * omegadxx[12] + 2 * omegadxx[13] +
                         2 * omegadxx[18] - 4 * omegadxx[19];
        phidxx[p20[7]] = -2 * omegadxx[6] + 4 * omegadxx[7] - 2 * omegadxx[12] - 2 * omegadxx[13] +
                         4 * omegadxx[18] - 2 * omegadxx[19];
        phidxx[p20[8]] = +4 * omegadxx[8] - 2 * omegadxx[9] + 2 * omegadxx[12] - 4 * omegadxx[13] +
                         2 * omegadxx[16] - 4 * omegadxx[17];
        phidxx[p20[9]] = -2 * omegadxx[8] + 4 * omegadxx[9] - 2 * omegadxx[12] - 2 * omegadxx[13] +
                         4 * omegadxx[16] - 2 * omegadxx[17];
        phidxx[p20[10]] = +4 * omegadxx[10] - 2 * omegadxx[11] + 2 * omegadxx[12] -
                          4 * omegadxx[13] + 2 * omegadxx[14] - 4 * omegadxx[15];
        phidxx[p20[11]] = -2 * omegadxx[10] + 4 * omegadxx[11] + 4 * omegadxx[12] -
                          2 * omegadxx[13] + 4 * omegadxx[14] - 2 * omegadxx[15];
        phidxx[p20[12]] = +8 * omegadxx[12] - 4 * omegadxx[13];
        phidxx[p20[13]] = -4 * omegadxx[12] + 8 * omegadxx[13];
        phidxx[p20[14]] = +8 * omegadxx[14] - 4 * omegadxx[15];
        phidxx[p20[15]] = -4 * omegadxx[14] + 8 * omegadxx[15];
        phidxx[p20[16]] = +8 * omegadxx[16] - 4 * omegadxx[17];
        phidxx[p20[17]] = -4 * omegadxx[16] + 8 * omegadxx[17];
        phidxx[p20[18]] = +8 * omegadxx[18] - 4 * omegadxx[19];
        phidxx[p20[19]] = -4 * omegadxx[18] + 8 * omegadxx[19];

        for (int k = 0; k < 20; ++k) {
          val(k, 0, op_dxx) = phidxx[k].x;
          val(k, 1, op_dxx) = phidxx[k].y;
          val(k, 2, op_dxx) = phidxx[k].z;
        }
      }

      R3 phidyy[20];
      if (whatd & Fop_dyy) {
        phidyy[p20[0]] = +4 * omegadyy[0] - 2 * omegadyy[1] - 4 * omegadyy[16] + 2 * omegadyy[17] -
                         4 * omegadyy[18] + 2 * omegadyy[19];
        phidyy[p20[1]] = -2 * omegadyy[0] + 4 * omegadyy[1] - 2 * omegadyy[16] - 2 * omegadyy[17] -
                         2 * omegadyy[18] - 2 * omegadyy[19];
        phidyy[p20[2]] = +4 * omegadyy[2] - 2 * omegadyy[3] - 4 * omegadyy[14] + 2 * omegadyy[15] +
                         2 * omegadyy[18] - 4 * omegadyy[19];
        phidyy[p20[3]] = -2 * omegadyy[2] + 4 * omegadyy[3] - 2 * omegadyy[14] - 2 * omegadyy[15] -
                         2 * omegadyy[18] - 2 * omegadyy[19];
        phidyy[p20[4]] = +4 * omegadyy[4] - 2 * omegadyy[5] + 2 * omegadyy[14] - 4 * omegadyy[15] +
                         2 * omegadyy[16] - 4 * omegadyy[17];
        phidyy[p20[5]] = -2 * omegadyy[4] + 4 * omegadyy[5] - 2 * omegadyy[14] - 2 * omegadyy[15] -
                         2 * omegadyy[16] - 2 * omegadyy[17];
        phidyy[p20[6]] = +4 * omegadyy[6] - 2 * omegadyy[7] - 4 * omegadyy[12] + 2 * omegadyy[13] +
                         2 * omegadyy[18] - 4 * omegadyy[19];
        phidyy[p20[7]] = -2 * omegadyy[6] + 4 * omegadyy[7] - 2 * omegadyy[12] - 2 * omegadyy[13] +
                         4 * omegadyy[18] - 2 * omegadyy[19];
        phidyy[p20[8]] = +4 * omegadyy[8] - 2 * omegadyy[9] + 2 * omegadyy[12] - 4 * omegadyy[13] +
                         2 * omegadyy[16] - 4 * omegadyy[17];
        phidyy[p20[9]] = -2 * omegadyy[8] + 4 * omegadyy[9] - 2 * omegadyy[12] - 2 * omegadyy[13] +
                         4 * omegadyy[16] - 2 * omegadyy[17];
        phidyy[p20[10]] = +4 * omegadyy[10] - 2 * omegadyy[11] + 2 * omegadyy[12] -
                          4 * omegadyy[13] + 2 * omegadyy[14] - 4 * omegadyy[15];
        phidyy[p20[11]] = -2 * omegadyy[10] + 4 * omegadyy[11] + 4 * omegadyy[12] -
                          2 * omegadyy[13] + 4 * omegadyy[14] - 2 * omegadyy[15];
        phidyy[p20[12]] = +8 * omegadyy[12] - 4 * omegadyy[13];
        phidyy[p20[13]] = -4 * omegadyy[12] + 8 * omegadyy[13];
        phidyy[p20[14]] = +8 * omegadyy[14] - 4 * omegadyy[15];
        phidyy[p20[15]] = -4 * omegadyy[14] + 8 * omegadyy[15];
        phidyy[p20[16]] = +8 * omegadyy[16] - 4 * omegadyy[17];
        phidyy[p20[17]] = -4 * omegadyy[16] + 8 * omegadyy[17];
        phidyy[p20[18]] = +8 * omegadyy[18] - 4 * omegadyy[19];
        phidyy[p20[19]] = -4 * omegadyy[18] + 8 * omegadyy[19];

        for (int k = 0; k < 20; ++k) {
          val(k, 0, op_dyy) = phidyy[k].x;
          val(k, 1, op_dyy) = phidyy[k].y;
          val(k, 2, op_dyy) = phidyy[k].z;
        }
      }

      R3 phidzz[20];
      if (whatd & Fop_dzz) {
        phidzz[p20[0]] = +4 * omegadzz[0] - 2 * omegadzz[1] - 4 * omegadzz[16] + 2 * omegadzz[17] -
                         4 * omegadzz[18] + 2 * omegadzz[19];
        phidzz[p20[1]] = -2 * omegadzz[0] + 4 * omegadzz[1] - 2 * omegadzz[16] - 2 * omegadzz[17] -
                         2 * omegadzz[18] - 2 * omegadzz[19];
        phidzz[p20[2]] = +4 * omegadzz[2] - 2 * omegadzz[3] - 4 * omegadzz[14] + 2 * omegadzz[15] +
                         2 * omegadzz[18] - 4 * omegadzz[19];
        phidzz[p20[3]] = -2 * omegadzz[2] + 4 * omegadzz[3] - 2 * omegadzz[14] - 2 * omegadzz[15] -
                         2 * omegadzz[18] - 2 * omegadzz[19];
        phidzz[p20[4]] = +4 * omegadzz[4] - 2 * omegadzz[5] + 2 * omegadzz[14] - 4 * omegadzz[15] +
                         2 * omegadzz[16] - 4 * omegadzz[17];
        phidzz[p20[5]] = -2 * omegadzz[4] + 4 * omegadzz[5] - 2 * omegadzz[14] - 2 * omegadzz[15] -
                         2 * omegadzz[16] - 2 * omegadzz[17];
        phidzz[p20[6]] = +4 * omegadzz[6] - 2 * omegadzz[7] - 4 * omegadzz[12] + 2 * omegadzz[13] +
                         2 * omegadzz[18] - 4 * omegadzz[19];
        phidzz[p20[7]] = -2 * omegadzz[6] + 4 * omegadzz[7] - 2 * omegadzz[12] - 2 * omegadzz[13] +
                         4 * omegadzz[18] - 2 * omegadzz[19];
        phidzz[p20[8]] = +4 * omegadzz[8] - 2 * omegadzz[9] + 2 * omegadzz[12] - 4 * omegadzz[13] +
                         2 * omegadzz[16] - 4 * omegadzz[17];
        phidzz[p20[9]] = -2 * omegadzz[8] + 4 * omegadzz[9] - 2 * omegadzz[12] - 2 * omegadzz[13] +
                         4 * omegadzz[16] - 2 * omegadzz[17];
        phidzz[p20[10]] = +4 * omegadzz[10] - 2 * omegadzz[11] + 2 * omegadzz[12] -
                          4 * omegadzz[13] + 2 * omegadzz[14] - 4 * omegadzz[15];
        phidzz[p20[11]] = -2 * omegadzz[10] + 4 * omegadzz[11] + 4 * omegadzz[12] -
                          2 * omegadzz[13] + 4 * omegadzz[14] - 2 * omegadzz[15];
        phidzz[p20[12]] = +8 * omegadzz[12] - 4 * omegadzz[13];
        phidzz[p20[13]] = -4 * omegadzz[12] + 8 * omegadzz[13];
        phidzz[p20[14]] = +8 * omegadzz[14] - 4 * omegadzz[15];
        phidzz[p20[15]] = -4 * omegadzz[14] + 8 * omegadzz[15];
        phidzz[p20[16]] = +8 * omegadzz[16] - 4 * omegadzz[17];
        phidzz[p20[17]] = -4 * omegadzz[16] + 8 * omegadzz[17];
        phidzz[p20[18]] = +8 * omegadzz[18] - 4 * omegadzz[19];
        phidzz[p20[19]] = -4 * omegadzz[18] + 8 * omegadzz[19];

        for (int k = 0; k < 20; ++k) {
          val(k, 0, op_dzz) = phidzz[k].x;
          val(k, 1, op_dzz) = phidzz[k].y;
          val(k, 2, op_dzz) = phidzz[k].z;
        }
      }

      R3 phidxy[20];
      if (whatd & Fop_dxy) {
        phidxy[p20[0]] = +4 * omegadxy[0] - 2 * omegadxy[1] - 4 * omegadxy[16] + 2 * omegadxy[17] -
                         4 * omegadxy[18] + 2 * omegadxy[19];
        phidxy[p20[1]] = -2 * omegadxy[0] + 4 * omegadxy[1] - 2 * omegadxy[16] - 2 * omegadxy[17] -
                         2 * omegadxy[18] - 2 * omegadxy[19];
        phidxy[p20[2]] = +4 * omegadxy[2] - 2 * omegadxy[3] - 4 * omegadxy[14] + 2 * omegadxy[15] +
                         2 * omegadxy[18] - 4 * omegadxy[19];
        phidxy[p20[3]] = -2 * omegadxy[2] + 4 * omegadxy[3] - 2 * omegadxy[14] - 2 * omegadxy[15] -
                         2 * omegadxy[18] - 2 * omegadxy[19];
        phidxy[p20[4]] = +4 * omegadxy[4] - 2 * omegadxy[5] + 2 * omegadxy[14] - 4 * omegadxy[15] +
                         2 * omegadxy[16] - 4 * omegadxy[17];
        phidxy[p20[5]] = -2 * omegadxy[4] + 4 * omegadxy[5] - 2 * omegadxy[14] - 2 * omegadxy[15] -
                         2 * omegadxy[16] - 2 * omegadxy[17];
        phidxy[p20[6]] = +4 * omegadxy[6] - 2 * omegadxy[7] - 4 * omegadxy[12] + 2 * omegadxy[13] +
                         2 * omegadxy[18] - 4 * omegadxy[19];
        phidxy[p20[7]] = -2 * omegadxy[6] + 4 * omegadxy[7] - 2 * omegadxy[12] - 2 * omegadxy[13] +
                         4 * omegadxy[18] - 2 * omegadxy[19];
        phidxy[p20[8]] = +4 * omegadxy[8] - 2 * omegadxy[9] + 2 * omegadxy[12] - 4 * omegadxy[13] +
                         2 * omegadxy[16] - 4 * omegadxy[17];
        phidxy[p20[9]] = -2 * omegadxy[8] + 4 * omegadxy[9] - 2 * omegadxy[12] - 2 * omegadxy[13] +
                         4 * omegadxy[16] - 2 * omegadxy[17];
        phidxy[p20[10]] = +4 * omegadxy[10] - 2 * omegadxy[11] + 2 * omegadxy[12] -
                          4 * omegadxy[13] + 2 * omegadxy[14] - 4 * omegadxy[15];
        phidxy[p20[11]] = -2 * omegadxy[10] + 4 * omegadxy[11] + 4 * omegadxy[12] -
                          2 * omegadxy[13] + 4 * omegadxy[14] - 2 * omegadxy[15];
        phidxy[p20[12]] = +8 * omegadxy[12] - 4 * omegadxy[13];
        phidxy[p20[13]] = -4 * omegadxy[12] + 8 * omegadxy[13];
        phidxy[p20[14]] = +8 * omegadxy[14] - 4 * omegadxy[15];
        phidxy[p20[15]] = -4 * omegadxy[14] + 8 * omegadxy[15];
        phidxy[p20[16]] = +8 * omegadxy[16] - 4 * omegadxy[17];
        phidxy[p20[17]] = -4 * omegadxy[16] + 8 * omegadxy[17];
        phidxy[p20[18]] = +8 * omegadxy[18] - 4 * omegadxy[19];
        phidxy[p20[19]] = -4 * omegadxy[18] + 8 * omegadxy[19];

        for (int k = 0; k < 20; ++k) {
          val(k, 0, op_dxy) = phidxy[k].x;
          val(k, 1, op_dxy) = phidxy[k].y;
          val(k, 2, op_dxy) = phidxy[k].z;
        }
      }

      R3 phidxz[20];
      if (whatd & Fop_dxz) {
        phidxz[p20[0]] = +4 * omegadxz[0] - 2 * omegadxz[1] - 4 * omegadxz[16] + 2 * omegadxz[17] -
                         4 * omegadxz[18] + 2 * omegadxz[19];
        phidxz[p20[1]] = -2 * omegadxz[0] + 4 * omegadxz[1] - 2 * omegadxz[16] - 2 * omegadxz[17] -
                         2 * omegadxz[18] - 2 * omegadxz[19];
        phidxz[p20[2]] = +4 * omegadxz[2] - 2 * omegadxz[3] - 4 * omegadxz[14] + 2 * omegadxz[15] +
                         2 * omegadxz[18] - 4 * omegadxz[19];
        phidxz[p20[3]] = -2 * omegadxz[2] + 4 * omegadxz[3] - 2 * omegadxz[14] - 2 * omegadxz[15] -
                         2 * omegadxz[18] - 2 * omegadxz[19];
        phidxz[p20[4]] = +4 * omegadxz[4] - 2 * omegadxz[5] + 2 * omegadxz[14] - 4 * omegadxz[15] +
                         2 * omegadxz[16] - 4 * omegadxz[17];
        phidxz[p20[5]] = -2 * omegadxz[4] + 4 * omegadxz[5] - 2 * omegadxz[14] - 2 * omegadxz[15] -
                         2 * omegadxz[16] - 2 * omegadxz[17];
        phidxz[p20[6]] = +4 * omegadxz[6] - 2 * omegadxz[7] - 4 * omegadxz[12] + 2 * omegadxz[13] +
                         2 * omegadxz[18] - 4 * omegadxz[19];
        phidxz[p20[7]] = -2 * omegadxz[6] + 4 * omegadxz[7] - 2 * omegadxz[12] - 2 * omegadxz[13] +
                         4 * omegadxz[18] - 2 * omegadxz[19];
        phidxz[p20[8]] = +4 * omegadxz[8] - 2 * omegadxz[9] + 2 * omegadxz[12] - 4 * omegadxz[13] +
                         2 * omegadxz[16] - 4 * omegadxz[17];
        phidxz[p20[9]] = -2 * omegadxz[8] + 4 * omegadxz[9] - 2 * omegadxz[12] - 2 * omegadxz[13] +
                         4 * omegadxz[16] - 2 * omegadxz[17];
        phidxz[p20[10]] = +4 * omegadxz[10] - 2 * omegadxz[11] + 2 * omegadxz[12] -
                          4 * omegadxz[13] + 2 * omegadxz[14] - 4 * omegadxz[15];
        phidxz[p20[11]] = -2 * omegadxz[10] + 4 * omegadxz[11] + 4 * omegadxz[12] -
                          2 * omegadxz[13] + 4 * omegadxz[14] - 2 * omegadxz[15];
        phidxz[p20[12]] = +8 * omegadxz[12] - 4 * omegadxz[13];
        phidxz[p20[13]] = -4 * omegadxz[12] + 8 * omegadxz[13];
        phidxz[p20[14]] = +8 * omegadxz[14] - 4 * omegadxz[15];
        phidxz[p20[15]] = -4 * omegadxz[14] + 8 * omegadxz[15];
        phidxz[p20[16]] = +8 * omegadxz[16] - 4 * omegadxz[17];
        phidxz[p20[17]] = -4 * omegadxz[16] + 8 * omegadxz[17];
        phidxz[p20[18]] = +8 * omegadxz[18] - 4 * omegadxz[19];
        phidxz[p20[19]] = -4 * omegadxz[18] + 8 * omegadxz[19];

        for (int k = 0; k < 20; ++k) {
          val(k, 0, op_dxz) = phidxz[k].x;
          val(k, 1, op_dxz) = phidxz[k].y;
          val(k, 2, op_dxz) = phidxz[k].z;
        }
      }

      R3 phidyz[20];
      if (whatd & Fop_dyz) {
        phidyz[p20[0]] = +4 * omegadyz[0] - 2 * omegadyz[1] - 4 * omegadyz[16] + 2 * omegadyz[17] -
                         4 * omegadyz[18] + 2 * omegadyz[19];
        phidyz[p20[1]] = -2 * omegadyz[0] + 4 * omegadyz[1] - 2 * omegadyz[16] - 2 * omegadyz[17] -
                         2 * omegadyz[18] - 2 * omegadyz[19];
        phidyz[p20[2]] = +4 * omegadyz[2] - 2 * omegadyz[3] - 4 * omegadyz[14] + 2 * omegadyz[15] +
                         2 * omegadyz[18] - 4 * omegadyz[19];
        phidyz[p20[3]] = -2 * omegadyz[2] + 4 * omegadyz[3] - 2 * omegadyz[14] - 2 * omegadyz[15] -
                         2 * omegadyz[18] - 2 * omegadyz[19];
        phidyz[p20[4]] = +4 * omegadyz[4] - 2 * omegadyz[5] + 2 * omegadyz[14] - 4 * omegadyz[15] +
                         2 * omegadyz[16] - 4 * omegadyz[17];
        phidyz[p20[5]] = -2 * omegadyz[4] + 4 * omegadyz[5] - 2 * omegadyz[14] - 2 * omegadyz[15] -
                         2 * omegadyz[16] - 2 * omegadyz[17];
        phidyz[p20[6]] = +4 * omegadyz[6] - 2 * omegadyz[7] - 4 * omegadyz[12] + 2 * omegadyz[13] +
                         2 * omegadyz[18] - 4 * omegadyz[19];
        phidyz[p20[7]] = -2 * omegadyz[6] + 4 * omegadyz[7] - 2 * omegadyz[12] - 2 * omegadyz[13] +
                         4 * omegadyz[18] - 2 * omegadyz[19];
        phidyz[p20[8]] = +4 * omegadyz[8] - 2 * omegadyz[9] + 2 * omegadyz[12] - 4 * omegadyz[13] +
                         2 * omegadyz[16] - 4 * omegadyz[17];
        phidyz[p20[9]] = -2 * omegadyz[8] + 4 * omegadyz[9] - 2 * omegadyz[12] - 2 * omegadyz[13] +
                         4 * omegadyz[16] - 2 * omegadyz[17];
        phidyz[p20[10]] = +4 * omegadyz[10] - 2 * omegadyz[11] + 2 * omegadyz[12] -
                          4 * omegadyz[13] + 2 * omegadyz[14] - 4 * omegadyz[15];
        phidyz[p20[11]] = -2 * omegadyz[10] + 4 * omegadyz[11] + 4 * omegadyz[12] -
                          2 * omegadyz[13] + 4 * omegadyz[14] - 2 * omegadyz[15];
        phidyz[p20[12]] = +8 * omegadyz[12] - 4 * omegadyz[13];
        phidyz[p20[13]] = -4 * omegadyz[12] + 8 * omegadyz[13];
        phidyz[p20[14]] = +8 * omegadyz[14] - 4 * omegadyz[15];
        phidyz[p20[15]] = -4 * omegadyz[14] + 8 * omegadyz[15];
        phidyz[p20[16]] = +8 * omegadyz[16] - 4 * omegadyz[17];
        phidyz[p20[17]] = -4 * omegadyz[16] + 8 * omegadyz[17];
        phidyz[p20[18]] = +8 * omegadyz[18] - 4 * omegadyz[19];
        phidyz[p20[19]] = -4 * omegadyz[18] + 8 * omegadyz[19];

        for (int k = 0; k < 20; ++k) {
          val(k, 0, op_dyz) = phidyz[k].x;
          val(k, 1, op_dyz) = phidyz[k].y;
          val(k, 2, op_dyz) = phidyz[k].z;
        }
      }
    }
  }

  // Author: Marcella Bonazzoli
  // Edge elements of degree 3, 3D
  // (they are called Edge23d because the Nedelec elements of degree 1 are called Edge03d)
  // TypeOfFE_Edge2_3d derived from GTypeOfFE<> (which is defined in FESpacen.hpp)
  class TypeOfFE_Edge2_3d : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    typedef GFElement< Mesh3 > FElement;

    static int dfon[];
    static const int d = Mesh::Rd::d;
    static const GQuadratureFormular< R1 > QFe;    // quadrature formula on an edge
    static const GQuadratureFormular< R2 > QFf;    // quadrature formula on a face
    static const GQuadratureFormular< R3 > QFv;    // quadrature formula on a tetrahedron
    TypeOfFE_Edge2_3d( );                          // constructor
    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
            RNMK_ &val) const;
    void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
             int *nump) const;
  };

  int TypeOfFE_Edge2_3d::dfon[] = {
    0, 3, 6, 3};    // 3 dofs on each edge, 6 dofs on each face, 3 dofs on the tetrahedron

  // Quadrature formula on an edge, exact for degree 5 (ok: int_e (deg3*t *lambda*lambda))
  const GQuadratureFormular< R1 > TypeOfFE_Edge2_3d::QFe(2 * 3 - 1, 3, GaussLegendre(3), true);
  // arguments: exact, num of integration pts, integration pts, clean (~GQuadratureFormular()
  // {if(clean) delete [] p;}) GaussLegendre defined in QuadratureFormular.cpp

  // Quadrature formula on a face, exact for degree 4 (ok: int_f (deg3*t*lambda))
  // (From Mark A. Taylor, Beth A. Wingate, Len P. Bos,
  // Several new quadrature formulas for polynomial integration in the triangle)
  static GQuadraturePoint< R2 > P_QuadratureFormular_T_4_TWB[6] = {
    GQuadraturePoint< R2 >(0.2199034873106 / 2, R2(0.0915762135098, 0.0915762135098)),
    GQuadraturePoint< R2 >(0.2199034873106 / 2, R2(0.8168475729805, 0.0915762135098)),
    GQuadraturePoint< R2 >(0.2199034873106 / 2, R2(0.0915762135098, 0.8168475729805)),
    GQuadraturePoint< R2 >(0.4467631793560 / 2, R2(0.1081030181681, 0.4459484909160)),
    GQuadraturePoint< R2 >(0.4467631793560 / 2, R2(0.4459484909160, 0.1081030181681)),
    GQuadraturePoint< R2 >(0.4467631793560 / 2, R2(0.4459484909160, 0.4459484909160))};
  const GQuadratureFormular< R2 > TypeOfFE_Edge2_3d::QFf(4, 6, P_QuadratureFormular_T_4_TWB);

  // Quadrature formula on a tetrahedron, exact for degree 3 (ok: int_v (deg3*t))
  // (From Table 4.57 pag 247 Solin HOFEM, coords trasformed like on pag 68,
  // weights multiplied by 6/8 so that their sum is equal to 1)
  static GQuadraturePoint< R3 > P_QuadratureFormular_Tet_3_solin[5] = {
    GQuadraturePoint< R3 >(R3(0.25, 0.25, 0.25), -0.8),
    GQuadraturePoint< R3 >(R3(1. / 6., 1. / 6., 1. / 6.), 0.45),
    GQuadraturePoint< R3 >(R3(1. / 6., 1. / 6., 0.5), 0.45),
    GQuadraturePoint< R3 >(R3(1. / 6., 0.5, 1. / 6.), 0.45),
    GQuadraturePoint< R3 >(R3(0.5, 1. / 6., 1. / 6.), 0.45)};
  const GQuadratureFormular< R3 > TypeOfFE_Edge2_3d::QFv(3, 5, P_QuadratureFormular_Tet_3_solin);

  // In Mesh3dn.cpp:
  // static const int  nvedgeTet[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
  // static const int  nvfaceTet[4][3] = { {3,2,1},{0,2,3},{3,1,0},{0,1,2} };
  // In GenericMesh.hpp:
  // Vertex& at(int i) {return *vertices[i];}
  // In GenericMesh.hpp:
  // Rd Edge(int i) const {ASSERTION(i>=0 && i <ne);
  // return Rd(at(nvedge[i][0]),at(nvedge[i][1]));}
  // In GenericMesh.hpp:
  // bool   EdgeOrientation(int i) const
  // { return &at(nvedge[i][0]) < &at(nvedge[i][1]);}

  // Constructor
  // dfon, d, nsub,
  // kPi = NbcoefforInterpolation is the number of alphas in (13.1) (chapter 13 ff++doc)
  // = 3(=numComp)*QFe.n(num quad pts edge)*3*ne(3 fncts per edge) +
  // 3(=numComp)*QFf.n(num quad pts face)*6*nf(6 fncts per face) +
  // 3(=numComp)*QFv.n(num quad pts volume)*3(3 volume fncts)
  // npPi = NbPtforInterpolation = ne*QFe.n+nf*QFf.n+QFv.n
  // invariantinterpolationMatrix = false (i.e. it depends on the tetrahedron), discon=true
  TypeOfFE_Edge2_3d::TypeOfFE_Edge2_3d( )
    : GTypeOfFE< Mesh3 >(TypeOfFE_Edge2_3d::dfon, d, 1,
                         3 * QFe.n * 3 * Element::ne + 3 * QFf.n * 6 * Element::nf + 3 * QFv.n * 3,
                         Element::ne * QFe.n + Element::nf * QFf.n + QFv.n, false, true) {
    assert(QFe.n);
    assert(QFf.n);
    assert(QFv.n);
    R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.),
               R3(0., 0., 1.)};    // 4 ref tetrahedron vertices

    {
      // We build the interpolation pts on the edges of the reference tetrahedron:
      int p;
      p = 0;

      for (int e = 0; e < Element::ne; ++e) {
        for (int q = 0; q < QFe.n; ++q, ++p) {
          double x = QFe[q].x;
          this->PtInterpolation[p] =
            Pt[Element::nvedge[e][0]] * (1. - x) + Pt[Element::nvedge[e][1]] * (x);
        }
      }

      // We build the interpolation pts on the faces of the reference tetrahedron:
      // (the index i mustn't be reinitialised!)
      for (int f = 0; f < Element::nf; ++f) {
        for (int q = 0; q < QFf.n; ++q, ++p) {
          double x = QFf[q].x;
          double y = QFf[q].y;
          this->PtInterpolation[p] = Pt[Element::nvface[f][0]] * (1. - x - y) +
                                     Pt[Element::nvface[f][1]] * x + Pt[Element::nvface[f][2]] * y;
        }
      }

      // We build the interpolation pts on the tetrahedron:
      // (the index i mustn't be reinitialised!)
      for (int q = 0; q < QFv.n; ++q, ++p) {
        double x = QFv[q].x;
        double y = QFv[q].y;
        double z = QFv[q].z;
        this->PtInterpolation[p] = Pt[0] * (1. - x - y - z) + Pt[1] * x + Pt[2] * y + Pt[3] * z;
      }
    }
    {
      // We build the indices in (13.1) : edge dofs
      int i = 0, p = 0;    // i is the k in (13.1) (chapter 13 ff++doc)
      int e;               // we will need e also below, in the parts referred to faces and volumes

      for (e = 0; e < (Element::ne)*3; e++) {    // loop on the 18 edge dofs
        if (e % 3 != 0) {
          p = p - QFe.n;
        }    // for 3 successive dofs the quad pts are the same (they correspond to the same edge)

        for (int q = 0; q < QFe.n; ++q, ++p) {    // loop on the edge quadrature pts
          for (int c = 0; c < 3; c++, i++) {      // loop on the 3 components
            this->pInterpolation[i] = p;          // pk in (13.1)
            this->cInterpolation[i] = c;          // jk in (13.1)
            this->dofInterpolation[i] = e;        // ik in (13.1)
            this->coefInterpolation[i] = 0.;      // alfak: we will fill them with 'set' (below)
                                                  // because they depend on the tetrahedron
          }
        }
      }

      // We build the indices in (13.1) : face dofs
      // (the indices i and p mustn't be reinitialised)
      int f;    // we will need f also below, in the part referred to volumes

      for (f = 0; f < (Element::nf)*6; f++) {    // loop on the 24 face dofs
        if (f % 6 != 0) {
          p = p - QFf.n;
        }    // for 6 successive dofs the quad pts are the same (they correspond to the same face)

        for (int q = 0; q < QFf.n; ++q, ++p) {    // loop on the face quadrature pts
          for (int c = 0; c < 3; c++, i++) {      // loop on the 3 components
            this->pInterpolation[i] = p;          // pk in (13.1)
            this->cInterpolation[i] = c;          // jk in (13.1)
            this->dofInterpolation[i] = e + f;    // ik in (13.1)
            this->coefInterpolation[i] = 0.;      // alphak: we will fill them with 'set' (below)
                                                  // because they depend on the tetrahedron
          }
        }
      }

      // We build the indices in (13.1) : volume dofs
      // (the indices i and p mustn't be reinitialised)
      for (int v = 0; v < 3; v++) {    // loop on the 3 volume dofs
        if (v > 0) {
          p = p - QFv.n;
        }

        for (int q = 0; q < QFv.n; ++q, ++p) {        // loop on the volume quadrature pts
          for (int c = 0; c < 3; c++, i++) {          // loop on the 3 components
            this->pInterpolation[i] = p;              // pk in (13.1)
            this->cInterpolation[i] = c;              // jk in (13.1)
            this->dofInterpolation[i] = e + f + v;    // ik in (13.1)
            this->coefInterpolation[i] = 0.;    // alphak: we will fill them with 'set' (below)
                                                // because they depend on the tetrahedron
          }
        }
      }
    }
  }

  // For the coefficients of interpolation alphak in (13.1)
  void TypeOfFE_Edge2_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
                              int ocoef, int odf, int *nump) const {
    int i = ocoef, p = 0;
    // (p is used only to check)

    // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL
    // number (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
    // perm[3] is the local number of the vertex with the biggest global number.)
    const Element::Vertex *tV[4] = {&K.at(0), &K.at(1), &K.at(2), &K.at(3)};
    int k0 = 0, k1 = 1, k2 = 2, k3 = 3;

    if (tV[k0] > tV[k1]) {
      Exchange(k0, k1);
    }

    if (tV[k1] > tV[k2]) {
      Exchange(k1, k2);
    }

    if (tV[k2] > tV[k3]) {
      Exchange(k2, k3);
    }

    if (tV[k0] > tV[k1]) {
      Exchange(k0, k1);
    }

    if (tV[k1] > tV[k2]) {
      Exchange(k1, k2);
    }

    if (tV[k0] > tV[k1]) {
      Exchange(k0, k1);
    }

    int perm[4] = {k0, k1, k2, k3};
    // here perm is used only for volume dofs (I wrote a version in which perm is exploited also for
    // edge and face dofs, but I don't know if it costs less)

    // 18 edge dofs
    for (int ee = 0; ee < Element::ne; ee++) {    // loop on the edges
      R3 E = K.EdgeOrientation(ee) * K.Edge(ee);
      int eo = K.EdgeOrientation(ee);
      if (eo < 0) {
        E = -E;
      }

      for (int edof = 0; edof < 3; ++edof) {    // 3 dofs for each edge
        if (edof != 0) {
          p = p - QFe.n;
        }

        for (int q = 0; q < QFe.n; ++q, ++p) {    // loop on the edge quadrature pts
          double ll = 1 - QFe[q].x;               // lambda of the first vertex of the edge
          if (eo < 0) {
            ll = 1 - ll;
          }

          double prodll;    // the product of two lambdas:
          if (edof == 0) {
            prodll = ll * ll;
          } else if (edof == 1) {
            prodll = ll * (1 - ll);
          } else {
            prodll = (1 - ll) * (1 - ll);    // if (edof==2)
          }

          for (int c = 0; c < 3; c++, i++) {    // loop on the 3 components
            M.coef[i] = E[c] * QFe[q].a * prodll;
            // QFe[q].a is the weight of the integration point q
            // QFe[q].x is the x over [0,1] of the integration point q
          }
        }
      }
    }

    // (the indices i and p mustn't be reinitialised)
    for (int ff = 0; ff < Element::nf; ff++) {    // loop on the 24 face dofs
      const Element::Vertex *fV[3] = {&K.at(Element::nvface[ff][0]), &K.at(Element::nvface[ff][1]),
                                      &K.at(Element::nvface[ff][2])};
      // We 'order' the 3 vertices of the face according to their global numbers:
      // i0f will be the local number in the FACE of the vertex with the smallest global number
      // i1f will be the local number in the face of the vertex with the second smallest global
      // number i2f will be the local number in the face of the vertex with the largest global
      // number
      int i0f = 0, i1f = 1, i2f = 2;
      if (fV[i0f] > fV[i1f]) {
        Exchange(i0f, i1f);
      }

      if (fV[i1f] > fV[i2f]) {
        Exchange(i1f, i2f);
        if (fV[i0f] > fV[i1f]) {
          Exchange(i0f, i1f);
        }
      }

      // now local numbers in the tetrahedron:
      int i0 = Element::nvface[ff][i0f], i1 = Element::nvface[ff][i1f],
          i2 = Element::nvface[ff][i2f];

      for (int fdof = 0; fdof < 6; ++fdof) {    // 6 dofs for each face
        int ie0 = i0,
            ie1 = fdof < 3 ? i1 : i2;    // edge for the face dof (its endpoints local numbers)
        // [i0,i1] for the first 3 dofs, [i0,i2] for the last 3 dofs
        R3 E(K[ie0], K[ie1]);
        if (fdof != 0) {
          p = p - QFf.n;
        }

        for (int q = 0; q < QFf.n; ++q, ++p) {    // loop on the face quadrature pts
          double ll3[3] = {1. - QFf[q].x - QFf[q].y, QFf[q].x,
                           QFf[q].y};    // the 3 lambdas of the face
          double ll;
          if (fdof % 3 == 0) {
            ll = ll3[i0f];
          } else if (fdof % 3 == 1) {
            ll = ll3[i1f];
          } else {
            ll = ll3[i2f];    // if (fdof%3==2)
          }

          for (int c = 0; c < 3; c++, i++) {    // loop on the 3 components
            M.coef[i] = E[c] * QFf[q].a * ll;
          }
        }
      }
    }

    // (the indices i and p mustn't be reinitialised)
    for (int v = 0; v < 3; v++) {          // loop on the 3 volume dofs
      R3 E(K[perm[0]], K[perm[v + 1]]);    // edge for the volume dof
      if (v > 0) {
        p = p - QFv.n;
      }

      for (int q = 0; q < QFv.n; ++q, ++p) {    // loop on the volume quadrature pts
        for (int c = 0; c < 3; c++, i++) {      // loop on the 3 components
          M.coef[i] = E[c] * QFv[q].a;
        }
      }
    }

    ffassert(i == M.ncoef && M.np == p);
  }

  // In Mesh3dn.hpp:
  // void Gradlambda(R3 * GradL) const
  // {
  // R3 V1(at(0),at(1));
  // R3 V2(at(0),at(2));
  // R3 V3(at(0),at(3));
  // R det1=1./(6.*mesure());
  // GradL[1]= (V2^V3)*det1;
  // GradL[2]= (V3^V1)*det1;
  // GradL[3]= (V1^V2)*det1;
  // GradL[0]=-GradL[1]-GradL[2]-GradL[3];
  // }

  /*
   * FB: We are on a tetrahedron K
   * (the local numbering of its 4 vertices is given by how they are listed where K is described in
   * the mesh file).
   *
   * AT FIRST we BUILD the 'pre-basis' functions OMEGA in the order given by the GLOBAL numbering of
   * the tetrahedron vertices, that is the 1st examined edge is from the vertex with the 1st
   * smallest GLOBAL number to the one with the 2nd smallest GLOBAL number, and so on (see nvedege),
   * and the 1st examined face is the one opposite the vertex with the smallest GLOBAL number
   * (and 2 edges of the face are chosen and oriented looking at its 3 GLOBAL numbers), and so on.
   * In this way:
   * - a dof which is common to 2 adjacent tetrahedra is the same in the two tetrahedra
   * (since the orientation of each edge is chosen using the global numbers),
   * - we can use the coefficients giving a dual basis that were calculated for the reference
   * tetrahedron since the structure of orientation of the edges is the same (up to a rotation) as
   * the one used in the reference tetrahedron.
   *
   * BUT THEN, when we build the DUAL basis functions PHI, we use the permutation p45 to go back to
   * the FreeFem numbering of the dofs, in which edges and faces are ordered following the LOCAL
   * numbering (for instance the 1st examined edge can be the 3rd edge looking at the local
   * numbering); however inside each edge (resp face) the 3 (resp 6) related dofs are still ordered
   * according to the GLOBAL numbers. Here with 'dual' I mean basis functions and dofs in duality.
   */

  // val contains the values of the basis functions and of their derivatives at the point of K
  // corresponding to the point P of the reference tetrahedron, by components
  void TypeOfFE_Edge2_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                             const RdHat &PHat, RNMK_ &val) const {
    assert(val.N( ) >= 45);    // 45 degrees of freedom
    assert(val.M( ) == 3);     // 3 components
    // -------------
    // perm: the permutation for which the 4 tetrahedron vertices are listed with increasing GLOBAL
    // number (i.e. perm[0] is the local number of the vertex with the smallest global number, ...
    // perm[3] is the local number of the vertex with the biggest global number.)
    const Element::Vertex *tV[4] = {&K.at(0), &K.at(1), &K.at(2), &K.at(3)};
    int k0 = 0, k1 = 1, k2 = 2, k3 = 3;
    if (tV[k0] > tV[k1]) {
      Exchange(k0, k1);
    }

    if (tV[k1] > tV[k2]) {
      Exchange(k1, k2);
    }

    if (tV[k2] > tV[k3]) {
      Exchange(k2, k3);
    }

    if (tV[k0] > tV[k1]) {
      Exchange(k0, k1);
    }

    if (tV[k1] > tV[k2]) {
      Exchange(k1, k2);
    }

    if (tV[k0] > tV[k1]) {
      Exchange(k0, k1);
    }

    int perm[4] = {k0, k1, k2, k3};
    // -------------
    // We build mynvface to be used here instead of the FreeFem nvface,
    // (in order to exploit the results of perm and write a better code),
    // in mynvface in all the triplets the numbers are increasing
    static const int mynvface[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
    // -------------
    // If [a,b] is the i-th edge (where a,b are its vertices local numbers), edgesMap[(a+1)*(b+1)] =
    // i
    int edgesMap[13] = {-1, -1, 0,  1,  2,  -1, 3,
                        -1, 4,  -1, -1, -1, 5};    // a map<int,int> would be more slow
    // -------------
    // the 4 barycentric coordinates for the reference tetrahedron evaluated at the point P
    // (they have the same value at the real tetrahedron's point corresponding to the reference
    // tetrahedron's point P)
    R l[] = {1. - PHat.sum( ), PHat.x, PHat.y, PHat.z};
    R3 D[4];
    K.Gradlambda(D);
    val = 0;
    // -------------
    int p45[45];    // the permutation from the dofs numbering of my tetrahedron (numbered using the
                    // global vertex numbers)

    // to the dofs numbering of FreeFem !!!!! (see the comment above)
    for (int i = 0; i < Element::ne; ++i) {    // edges (3 dofs for each edge)
      int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
      int i0 = perm[ii0];
      int i1 = perm[ii1];
      int iEdge = edgesMap[(i0 + 1) * (i1 + 1)];    // i of the edge [i0,i1]
      p45[i * 3] = iEdge * 3;
      p45[i * 3 + 1] = iEdge * 3 + 1;
      p45[i * 3 + 2] = iEdge * 3 + 2;
    }

    for (int j = 0; j < Element::nf; ++j) {    // faces (6 dofs for each face)
      // using perm the 1st examined face is the one opposite the vertex with the smallest GLOBAL
      // number, and so on
      int jFace = perm[j];
      int dof0j = 18 + j * 6;
      int dof0jF = 18 + jFace * 6;

      for (int s = 0; s < 6; ++s) {
        p45[dof0j + s] = dof0jF + s;
      }    // p45[18+j*6+0] = 18+jFace*6+0; ... p45[18+j*6+5] = 18+jFace*6+5;
    }

    // 3 dofs for the tetrahedron
    p45[42] = 42;
    p45[43] = 43;
    p45[44] = 44;    // (there is only 1 'volume', the FF dofs numbering is the same as mine)
    // -------------

    if (whatd & Fop_D0) {    // Fop_D0 defined in FESpacen.hpp
      // First, the functions omega (om) (they don't constitute a dual basis! only a basis)
      R3 om[45];

      // 18 edge functions (3 for each edge):
      for (int i = 0; i < Element::ne; ++i) {
        int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        // ! :
        // using perm, [i0,i1] is already from the smallest global number to the greatest global
        // number, since nvedge always gives indices which read perm from left to right
        om[i * 3] = l[i0] * l[i0] * (l[i0] * D[i1] - l[i1] * D[i0]);
        om[i * 3 + 1] = l[i0] * l[i1] * (l[i0] * D[i1] - l[i1] * D[i0]);
        om[i * 3 + 2] = l[i1] * l[i1] * (l[i0] * D[i1] - l[i1] * D[i0]);
      }

      // 24 face functions (6 for each face):
      for (int j = 0; j < Element::nf; ++j) {
        // In (my)nvface the face opposite the vertex k is (my)nvface[k][],
        // using perm the 1st examined face is the one opposite the vertex with the smallest GLOBAL
        // number, and so on, and, since in mynvface the numbers always increase, [i0,i1,i2] are
        // already ordered with increasing global number
        int ii0 = mynvface[j][0];
        int ii1 = mynvface[j][1];
        int ii2 = mynvface[j][2];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        int i2 = perm[ii2];
        // (edge i0i1, which requires i2, and in front put a lambda of the face ordered vertices)
        om[18 + j * 6] = l[i0] * l[i2] * (l[i0] * D[i1] - l[i1] * D[i0]);
        om[18 + j * 6 + 1] = l[i1] * l[i2] * (l[i0] * D[i1] - l[i1] * D[i0]);
        om[18 + j * 6 + 2] = l[i2] * l[i2] * (l[i0] * D[i1] - l[i1] * D[i0]);
        // (edge i0i2, which requires i1, and in front put a lambda of the face ordered vertices)
        om[18 + j * 6 + 3] = l[i0] * l[i1] * (l[i0] * D[i2] - l[i2] * D[i0]);
        om[18 + j * 6 + 4] = l[i1] * l[i1] * (l[i0] * D[i2] - l[i2] * D[i0]);
        om[18 + j * 6 + 5] = l[i2] * l[i1] * (l[i0] * D[i2] - l[i2] * D[i0]);
      }

      // 3 volume functions
      // (Although the volume dofs are not shared between adjacent tetrahedra,
      // I have to use perm to be allowed to use the calculated coefficients giving a dual basis)
      om[42] = l[perm[2]] * l[perm[3]] * (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]);
      om[43] = l[perm[1]] * l[perm[3]] * (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]);
      om[44] = l[perm[1]] * l[perm[2]] * (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]);

      // Now, the functions phi that really constitute a dual basis
      R3 phi[45];
      phi[p45[0]] = +9 * om[0] - 18 * om[1] + 3 * om[2] - 27 * om[30] + 12 * om[31] + 9 * om[32] +
                    9 * om[33] - 6 * om[34] - 6 * om[35] - 27 * om[36] + 12 * om[37] + 9 * om[38] +
                    9 * om[39] - 6 * om[40] - 6 * om[41] + 18 * om[42] - 6 * om[43] - 6 * om[44];
      phi[p45[1]] = -18 * om[0] + 84 * om[1] - 18 * om[2] - 6 * om[30] - 36 * om[31] + 12 * om[32] -
                    30 * om[33] + 30 * om[34] - 6 * om[36] - 36 * om[37] + 12 * om[38] -
                    30 * om[39] + 30 * om[40] + 24 * om[42];
      phi[p45[2]] = +3 * om[0] - 18 * om[1] + 9 * om[2] + 6 * om[30] - 18 * om[31] + 3 * om[32] +
                    6 * om[33] - 9 * om[34] + 6 * om[35] + 6 * om[36] - 18 * om[37] + 3 * om[38] +
                    6 * om[39] - 9 * om[40] + 6 * om[41] + 6 * om[42] + 6 * om[43] + 6 * om[44];
      phi[p45[3]] = +9 * om[3] - 18 * om[4] + 3 * om[5] - 27 * om[24] + 12 * om[25] + 9 * om[26] +
                    9 * om[27] - 6 * om[28] - 6 * om[29] + 9 * om[36] - 6 * om[37] - 6 * om[38] -
                    27 * om[39] + 9 * om[40] + 12 * om[41] - 6 * om[42] + 18 * om[43] - 6 * om[44];
      phi[p45[4]] = -18 * om[3] + 84 * om[4] - 18 * om[5] - 6 * om[24] - 36 * om[25] + 12 * om[26] -
                    30 * om[27] + 30 * om[28] - 30 * om[36] + 30 * om[38] - 6 * om[39] +
                    12 * om[40] - 36 * om[41] + 24 * om[43];
      phi[p45[5]] = +3 * om[3] - 18 * om[4] + 9 * om[5] + 6 * om[24] - 18 * om[25] + 3 * om[26] +
                    6 * om[27] - 9 * om[28] + 6 * om[29] + 6 * om[36] + 6 * om[37] - 9 * om[38] +
                    6 * om[39] + 3 * om[40] - 18 * om[41] + 6 * om[42] + 6 * om[43] + 6 * om[44];
      phi[p45[6]] = +9 * om[6] - 18 * om[7] + 3 * om[8] + 9 * om[24] - 6 * om[25] - 6 * om[26] -
                    27 * om[27] + 9 * om[28] + 12 * om[29] + 9 * om[30] - 6 * om[31] - 6 * om[32] -
                    27 * om[33] + 9 * om[34] + 12 * om[35] - 6 * om[42] - 6 * om[43] + 18 * om[44];
      phi[p45[7]] = -18 * om[6] + 84 * om[7] - 18 * om[8] - 30 * om[24] + 30 * om[26] - 6 * om[27] +
                    12 * om[28] - 36 * om[29] - 30 * om[30] + 30 * om[32] - 6 * om[33] +
                    12 * om[34] - 36 * om[35] + 24 * om[44];
      phi[p45[8]] = +3 * om[6] - 18 * om[7] + 9 * om[8] + 6 * om[24] + 6 * om[25] - 9 * om[26] +
                    6 * om[27] + 3 * om[28] - 18 * om[29] + 6 * om[30] + 6 * om[31] - 9 * om[32] +
                    6 * om[33] + 3 * om[34] - 18 * om[35] + 6 * om[42] + 6 * om[43] + 6 * om[44];
      phi[p45[9]] = +9 * om[9] - 18 * om[10] + 3 * om[11] - 27 * om[18] + 12 * om[19] + 9 * om[20] +
                    9 * om[21] - 6 * om[22] - 6 * om[23] - 3 * om[36] + 18 * om[37] - 6 * om[38] +
                    9 * om[39] - 27 * om[40] + 12 * om[41] - 6 * om[42] + 18 * om[43] - 6 * om[44];
      phi[p45[10]] = -18 * om[9] + 84 * om[10] - 18 * om[11] - 6 * om[18] - 36 * om[19] +
                     12 * om[20] - 30 * om[21] + 30 * om[22] - 12 * om[36] + 36 * om[37] +
                     6 * om[38] + 12 * om[39] - 6 * om[40] - 36 * om[41] - 24 * om[42] +
                     24 * om[43];
      phi[p45[11]] = +3 * om[9] - 18 * om[10] + 9 * om[11] + 6 * om[18] - 18 * om[19] + 3 * om[20] +
                     6 * om[21] - 9 * om[22] + 6 * om[23] - 9 * om[36] - 12 * om[37] + 27 * om[38] +
                     3 * om[39] + 6 * om[40] - 18 * om[41] - 18 * om[42] + 6 * om[43] + 6 * om[44];
      phi[p45[12]] = +9 * om[12] - 18 * om[13] + 3 * om[14] + 9 * om[18] - 6 * om[19] - 6 * om[20] -
                     27 * om[21] + 9 * om[22] + 12 * om[23] - 3 * om[30] + 18 * om[31] -
                     6 * om[32] + 9 * om[33] - 27 * om[34] + 12 * om[35] - 6 * om[42] - 6 * om[43] +
                     18 * om[44];
      phi[p45[13]] = -18 * om[12] + 84 * om[13] - 18 * om[14] - 30 * om[18] + 30 * om[20] -
                     6 * om[21] + 12 * om[22] - 36 * om[23] - 12 * om[30] + 36 * om[31] +
                     6 * om[32] + 12 * om[33] - 6 * om[34] - 36 * om[35] - 24 * om[42] +
                     24 * om[44];
      phi[p45[14]] = +3 * om[12] - 18 * om[13] + 9 * om[14] + 6 * om[18] + 6 * om[19] - 9 * om[20] +
                     6 * om[21] + 3 * om[22] - 18 * om[23] - 9 * om[30] - 12 * om[31] +
                     27 * om[32] + 3 * om[33] + 6 * om[34] - 18 * om[35] - 18 * om[42] +
                     6 * om[43] + 6 * om[44];
      phi[p45[15]] = +9 * om[15] - 18 * om[16] + 3 * om[17] - 3 * om[18] + 18 * om[19] -
                     6 * om[20] + 9 * om[21] - 27 * om[22] + 12 * om[23] - 3 * om[24] +
                     18 * om[25] - 6 * om[26] + 9 * om[27] - 27 * om[28] + 12 * om[29] -
                     6 * om[42] - 6 * om[43] + 18 * om[44];
      phi[p45[16]] = -18 * om[15] + 84 * om[16] - 18 * om[17] - 12 * om[18] + 36 * om[19] +
                     6 * om[20] + 12 * om[21] - 6 * om[22] - 36 * om[23] - 12 * om[24] +
                     36 * om[25] + 6 * om[26] + 12 * om[27] - 6 * om[28] - 36 * om[29] -
                     24 * om[43] + 24 * om[44];
      phi[p45[17]] = +3 * om[15] - 18 * om[16] + 9 * om[17] - 9 * om[18] - 12 * om[19] +
                     27 * om[20] + 3 * om[21] + 6 * om[22] - 18 * om[23] - 9 * om[24] -
                     12 * om[25] + 27 * om[26] + 3 * om[27] + 6 * om[28] - 18 * om[29] +
                     6 * om[42] - 18 * om[43] + 6 * om[44];
      phi[p45[18]] = +90 * om[18] - 30 * om[19] - 45 * om[20] - 30 * om[21] + 15 * om[22] +
                     30 * om[23] + 15 * om[42] - 45 * om[43] + 15 * om[44];
      phi[p45[19]] = -30 * om[18] + 120 * om[19] - 30 * om[20] + 30 * om[21] - 60 * om[22] +
                     30 * om[42] - 30 * om[43] + 30 * om[44];
      phi[p45[20]] = -45 * om[18] - 30 * om[19] + 90 * om[20] + 15 * om[21] + 15 * om[22] -
                     60 * om[23] + 15 * om[42] - 45 * om[43] + 15 * om[44];
      phi[p45[21]] = -30 * om[18] + 30 * om[19] + 15 * om[20] + 90 * om[21] - 45 * om[22] -
                     30 * om[23] + 15 * om[42] + 15 * om[43] - 45 * om[44];
      phi[p45[22]] = +15 * om[18] - 60 * om[19] + 15 * om[20] - 45 * om[21] + 90 * om[22] -
                     30 * om[23] + 15 * om[42] + 15 * om[43] - 45 * om[44];
      phi[p45[23]] = +30 * om[18] - 60 * om[20] - 30 * om[21] - 30 * om[22] + 120 * om[23] +
                     30 * om[42] + 30 * om[43] - 30 * om[44];
      phi[p45[24]] = +90 * om[24] - 30 * om[25] - 45 * om[26] - 30 * om[27] + 15 * om[28] +
                     30 * om[29] + 15 * om[42] - 45 * om[43] + 15 * om[44];
      phi[p45[25]] = -30 * om[24] + 120 * om[25] - 30 * om[26] + 30 * om[27] - 60 * om[28] -
                     30 * om[42] - 30 * om[43] + 30 * om[44];
      phi[p45[26]] = -45 * om[24] - 30 * om[25] + 90 * om[26] + 15 * om[27] + 15 * om[28] -
                     60 * om[29] + 15 * om[42] - 45 * om[43] + 15 * om[44];
      phi[p45[27]] = -30 * om[24] + 30 * om[25] + 15 * om[26] + 90 * om[27] - 45 * om[28] -
                     30 * om[29] + 15 * om[42] + 15 * om[43] - 45 * om[44];
      phi[p45[28]] = +15 * om[24] - 60 * om[25] + 15 * om[26] - 45 * om[27] + 90 * om[28] -
                     30 * om[29] + 15 * om[42] + 15 * om[43] - 45 * om[44];
      phi[p45[29]] = +30 * om[24] - 60 * om[26] - 30 * om[27] - 30 * om[28] + 120 * om[29] -
                     30 * om[42] + 30 * om[43] - 30 * om[44];
      phi[p45[30]] = +90 * om[30] - 30 * om[31] - 45 * om[32] - 30 * om[33] + 15 * om[34] +
                     30 * om[35] - 45 * om[42] + 15 * om[43] + 15 * om[44];
      phi[p45[31]] = -30 * om[30] + 120 * om[31] - 30 * om[32] + 30 * om[33] - 60 * om[34] -
                     30 * om[42] - 30 * om[43] + 30 * om[44];
      phi[p45[32]] = -45 * om[30] - 30 * om[31] + 90 * om[32] + 15 * om[33] + 15 * om[34] -
                     60 * om[35] - 45 * om[42] + 15 * om[43] + 15 * om[44];
      phi[p45[33]] = -30 * om[30] + 30 * om[31] + 15 * om[32] + 90 * om[33] - 45 * om[34] -
                     30 * om[35] + 15 * om[42] + 15 * om[43] - 45 * om[44];
      phi[p45[34]] = +15 * om[30] - 60 * om[31] + 15 * om[32] - 45 * om[33] + 90 * om[34] -
                     30 * om[35] + 15 * om[42] + 15 * om[43] - 45 * om[44];
      phi[p45[35]] = +30 * om[30] - 60 * om[32] - 30 * om[33] - 30 * om[34] + 120 * om[35] +
                     30 * om[42] - 30 * om[43] - 30 * om[44];
      phi[p45[36]] = +90 * om[36] - 30 * om[37] - 45 * om[38] - 30 * om[39] + 15 * om[40] +
                     30 * om[41] - 45 * om[42] + 15 * om[43] + 15 * om[44];
      phi[p45[37]] = -30 * om[36] + 120 * om[37] - 30 * om[38] + 30 * om[39] - 60 * om[40] -
                     30 * om[42] + 30 * om[43] - 30 * om[44];
      phi[p45[38]] = -45 * om[36] - 30 * om[37] + 90 * om[38] + 15 * om[39] + 15 * om[40] -
                     60 * om[41] - 45 * om[42] + 15 * om[43] + 15 * om[44];
      phi[p45[39]] = -30 * om[36] + 30 * om[37] + 15 * om[38] + 90 * om[39] - 45 * om[40] -
                     30 * om[41] + 15 * om[42] - 45 * om[43] + 15 * om[44];
      phi[p45[40]] = +15 * om[36] - 60 * om[37] + 15 * om[38] - 45 * om[39] + 90 * om[40] -
                     30 * om[41] + 15 * om[42] - 45 * om[43] + 15 * om[44];
      phi[p45[41]] = +30 * om[36] - 60 * om[38] - 30 * om[39] - 30 * om[40] + 120 * om[41] +
                     30 * om[42] - 30 * om[43] - 30 * om[44];
      phi[p45[42]] = +90 * om[42] - 30 * om[43] - 30 * om[44];
      phi[p45[43]] = -30 * om[42] + 90 * om[43] - 30 * om[44];
      phi[p45[44]] = -30 * om[42] - 30 * om[43] + 90 * om[44];

      for (int k = 0; k < 45; ++k) {
        val(k, 0, op_id) = phi[k].x;
        val(k, 1, op_id) = phi[k].y;
        val(k, 2, op_id) = phi[k].z;
      }
    }

    if (whatd & Fop_D1) {    // FIRST DERIVATIVES wrt x,y,z
      R3 omdx[45];
      R3 omdy[45];
      R3 omdz[45];

      // 18 edge functions (3 for each edge):
      for (int i = 0; i < Element::ne; ++i) {
        int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        if (whatd & Fop_dx) {
          omdx[i * 3] = 2 * D[i0].x * l[i0] * (l[i0] * D[i1] - l[i1] * D[i0]) +
                        l[i0] * l[i0] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdx[i * 3 + 1] = (D[i0].x * l[i1] + l[i0] * D[i1].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
                            l[i0] * l[i1] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdx[i * 3 + 2] = 2 * D[i1].x * l[i1] * (l[i0] * D[i1] - l[i1] * D[i0]) +
                            l[i1] * l[i1] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
        }

        if (whatd & Fop_dy) {
          omdy[i * 3] = 2 * D[i0].y * l[i0] * (l[i0] * D[i1] - l[i1] * D[i0]) +
                        l[i0] * l[i0] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdy[i * 3 + 1] = (D[i0].y * l[i1] + l[i0] * D[i1].y) * (l[i0] * D[i1] - l[i1] * D[i0]) +
                            l[i0] * l[i1] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdy[i * 3 + 2] = 2 * D[i1].y * l[i1] * (l[i0] * D[i1] - l[i1] * D[i0]) +
                            l[i1] * l[i1] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
        }

        if (whatd & Fop_dz) {
          omdz[i * 3] = 2 * D[i0].z * l[i0] * (l[i0] * D[i1] - l[i1] * D[i0]) +
                        l[i0] * l[i0] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdz[i * 3 + 1] = (D[i0].z * l[i1] + l[i0] * D[i1].z) * (l[i0] * D[i1] - l[i1] * D[i0]) +
                            l[i0] * l[i1] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdz[i * 3 + 2] = 2 * D[i1].z * l[i1] * (l[i0] * D[i1] - l[i1] * D[i0]) +
                            l[i1] * l[i1] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
        }
      }

      // 24 face functions (6 for each face):
      for (int j = 0; j < Element::nf; ++j) {
        int ii0 = mynvface[j][0];
        int ii1 = mynvface[j][1];
        int ii2 = mynvface[j][2];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        int i2 = perm[ii2];
        if (whatd & Fop_dx) {
          omdx[18 + j * 6] = (D[i0].x * l[i2] + l[i0] * D[i2].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
                             l[i0] * l[i2] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdx[18 + j * 6 + 1] =
            (D[i1].x * l[i2] + l[i1] * D[i2].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            l[i1] * l[i2] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdx[18 + j * 6 + 2] =
            (D[i2].x * l[i2] + l[i2] * D[i2].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            l[i2] * l[i2] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdx[18 + j * 6 + 3] =
            (D[i0].x * l[i1] + l[i0] * D[i1].x) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            l[i0] * l[i1] * (D[i0].x * D[i2] - D[i2].x * D[i0]);
          omdx[18 + j * 6 + 4] =
            (D[i1].x * l[i1] + l[i1] * D[i1].x) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            l[i1] * l[i1] * (D[i0].x * D[i2] - D[i2].x * D[i0]);
          omdx[18 + j * 6 + 5] =
            (D[i2].x * l[i1] + l[i2] * D[i1].x) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            l[i2] * l[i1] * (D[i0].x * D[i2] - D[i2].x * D[i0]);
        }

        if (whatd & Fop_dy) {
          omdy[18 + j * 6] = (D[i0].y * l[i2] + l[i0] * D[i2].y) * (l[i0] * D[i1] - l[i1] * D[i0]) +
                             l[i0] * l[i2] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdy[18 + j * 6 + 1] =
            (D[i1].y * l[i2] + l[i1] * D[i2].y) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            l[i1] * l[i2] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdy[18 + j * 6 + 2] =
            (D[i2].y * l[i2] + l[i2] * D[i2].y) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            l[i2] * l[i2] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdy[18 + j * 6 + 3] =
            (D[i0].y * l[i1] + l[i0] * D[i1].y) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            l[i0] * l[i1] * (D[i0].y * D[i2] - D[i2].y * D[i0]);
          omdy[18 + j * 6 + 4] =
            (D[i1].y * l[i1] + l[i1] * D[i1].y) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            l[i1] * l[i1] * (D[i0].y * D[i2] - D[i2].y * D[i0]);
          omdy[18 + j * 6 + 5] =
            (D[i2].y * l[i1] + l[i2] * D[i1].y) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            l[i2] * l[i1] * (D[i0].y * D[i2] - D[i2].y * D[i0]);
        }

        if (whatd & Fop_dz) {
          omdz[18 + j * 6] = (D[i0].z * l[i2] + l[i0] * D[i2].z) * (l[i0] * D[i1] - l[i1] * D[i0]) +
                             l[i0] * l[i2] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdz[18 + j * 6 + 1] =
            (D[i1].z * l[i2] + l[i1] * D[i2].z) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            l[i1] * l[i2] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdz[18 + j * 6 + 2] =
            (D[i2].z * l[i2] + l[i2] * D[i2].z) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            l[i2] * l[i2] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdz[18 + j * 6 + 3] =
            (D[i0].z * l[i1] + l[i0] * D[i1].z) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            l[i0] * l[i1] * (D[i0].z * D[i2] - D[i2].z * D[i0]);
          omdz[18 + j * 6 + 4] =
            (D[i1].z * l[i1] + l[i1] * D[i1].z) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            l[i1] * l[i1] * (D[i0].z * D[i2] - D[i2].z * D[i0]);
          omdz[18 + j * 6 + 5] =
            (D[i2].z * l[i1] + l[i2] * D[i1].z) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            l[i2] * l[i1] * (D[i0].z * D[i2] - D[i2].z * D[i0]);
        }
      }

      // 3 volume functions
      if (whatd & Fop_dx) {
        omdx[42] =
          (D[perm[2]].x * l[perm[3]] + l[perm[2]] * D[perm[3]].x) *
            (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]) +
          l[perm[2]] * l[perm[3]] * (D[perm[0]].x * D[perm[1]] - D[perm[1]].x * D[perm[0]]);
        omdx[43] =
          (D[perm[1]].x * l[perm[3]] + l[perm[1]] * D[perm[3]].x) *
            (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]) +
          l[perm[1]] * l[perm[3]] * (D[perm[0]].x * D[perm[2]] - D[perm[2]].x * D[perm[0]]);
        omdx[44] =
          (D[perm[1]].x * l[perm[2]] + l[perm[1]] * D[perm[2]].x) *
            (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]) +
          l[perm[1]] * l[perm[2]] * (D[perm[0]].x * D[perm[3]] - D[perm[3]].x * D[perm[0]]);
      }

      if (whatd & Fop_dy) {
        omdy[42] =
          (D[perm[2]].y * l[perm[3]] + l[perm[2]] * D[perm[3]].y) *
            (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]) +
          l[perm[2]] * l[perm[3]] * (D[perm[0]].y * D[perm[1]] - D[perm[1]].y * D[perm[0]]);
        omdy[43] =
          (D[perm[1]].y * l[perm[3]] + l[perm[1]] * D[perm[3]].y) *
            (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]) +
          l[perm[1]] * l[perm[3]] * (D[perm[0]].y * D[perm[2]] - D[perm[2]].y * D[perm[0]]);
        omdy[44] =
          (D[perm[1]].y * l[perm[2]] + l[perm[1]] * D[perm[2]].y) *
            (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]) +
          l[perm[1]] * l[perm[2]] * (D[perm[0]].y * D[perm[3]] - D[perm[3]].y * D[perm[0]]);
      }

      if (whatd & Fop_dz) {
        omdz[42] =
          (D[perm[2]].z * l[perm[3]] + l[perm[2]] * D[perm[3]].z) *
            (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]) +
          l[perm[2]] * l[perm[3]] * (D[perm[0]].z * D[perm[1]] - D[perm[1]].z * D[perm[0]]);
        omdz[43] =
          (D[perm[1]].z * l[perm[3]] + l[perm[1]] * D[perm[3]].z) *
            (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]) +
          l[perm[1]] * l[perm[3]] * (D[perm[0]].z * D[perm[2]] - D[perm[2]].z * D[perm[0]]);
        omdz[44] =
          (D[perm[1]].z * l[perm[2]] + l[perm[1]] * D[perm[2]].z) *
            (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]) +
          l[perm[1]] * l[perm[2]] * (D[perm[0]].z * D[perm[3]] - D[perm[3]].z * D[perm[0]]);
      }

      R3 phidx[45];
      if (whatd & Fop_dx) {
        phidx[p45[0]] = +9 * omdx[0] - 18 * omdx[1] + 3 * omdx[2] - 27 * omdx[30] + 12 * omdx[31] +
                        9 * omdx[32] + 9 * omdx[33] - 6 * omdx[34] - 6 * omdx[35] - 27 * omdx[36] +
                        12 * omdx[37] + 9 * omdx[38] + 9 * omdx[39] - 6 * omdx[40] - 6 * omdx[41] +
                        18 * omdx[42] - 6 * omdx[43] - 6 * omdx[44];
        phidx[p45[1]] = -18 * omdx[0] + 84 * omdx[1] - 18 * omdx[2] - 6 * omdx[30] - 36 * omdx[31] +
                        12 * omdx[32] - 30 * omdx[33] + 30 * omdx[34] - 6 * omdx[36] -
                        36 * omdx[37] + 12 * omdx[38] - 30 * omdx[39] + 30 * omdx[40] +
                        24 * omdx[42];
        phidx[p45[2]] = +3 * omdx[0] - 18 * omdx[1] + 9 * omdx[2] + 6 * omdx[30] - 18 * omdx[31] +
                        3 * omdx[32] + 6 * omdx[33] - 9 * omdx[34] + 6 * omdx[35] + 6 * omdx[36] -
                        18 * omdx[37] + 3 * omdx[38] + 6 * omdx[39] - 9 * omdx[40] + 6 * omdx[41] +
                        6 * omdx[42] + 6 * omdx[43] + 6 * omdx[44];
        phidx[p45[3]] = +9 * omdx[3] - 18 * omdx[4] + 3 * omdx[5] - 27 * omdx[24] + 12 * omdx[25] +
                        9 * omdx[26] + 9 * omdx[27] - 6 * omdx[28] - 6 * omdx[29] + 9 * omdx[36] -
                        6 * omdx[37] - 6 * omdx[38] - 27 * omdx[39] + 9 * omdx[40] + 12 * omdx[41] -
                        6 * omdx[42] + 18 * omdx[43] - 6 * omdx[44];
        phidx[p45[4]] = -18 * omdx[3] + 84 * omdx[4] - 18 * omdx[5] - 6 * omdx[24] - 36 * omdx[25] +
                        12 * omdx[26] - 30 * omdx[27] + 30 * omdx[28] - 30 * omdx[36] +
                        30 * omdx[38] - 6 * omdx[39] + 12 * omdx[40] - 36 * omdx[41] +
                        24 * omdx[43];
        phidx[p45[5]] = +3 * omdx[3] - 18 * omdx[4] + 9 * omdx[5] + 6 * omdx[24] - 18 * omdx[25] +
                        3 * omdx[26] + 6 * omdx[27] - 9 * omdx[28] + 6 * omdx[29] + 6 * omdx[36] +
                        6 * omdx[37] - 9 * omdx[38] + 6 * omdx[39] + 3 * omdx[40] - 18 * omdx[41] +
                        6 * omdx[42] + 6 * omdx[43] + 6 * omdx[44];
        phidx[p45[6]] = +9 * omdx[6] - 18 * omdx[7] + 3 * omdx[8] + 9 * omdx[24] - 6 * omdx[25] -
                        6 * omdx[26] - 27 * omdx[27] + 9 * omdx[28] + 12 * omdx[29] + 9 * omdx[30] -
                        6 * omdx[31] - 6 * omdx[32] - 27 * omdx[33] + 9 * omdx[34] + 12 * omdx[35] -
                        6 * omdx[42] - 6 * omdx[43] + 18 * omdx[44];
        phidx[p45[7]] = -18 * omdx[6] + 84 * omdx[7] - 18 * omdx[8] - 30 * omdx[24] +
                        30 * omdx[26] - 6 * omdx[27] + 12 * omdx[28] - 36 * omdx[29] -
                        30 * omdx[30] + 30 * omdx[32] - 6 * omdx[33] + 12 * omdx[34] -
                        36 * omdx[35] + 24 * omdx[44];
        phidx[p45[8]] = +3 * omdx[6] - 18 * omdx[7] + 9 * omdx[8] + 6 * omdx[24] + 6 * omdx[25] -
                        9 * omdx[26] + 6 * omdx[27] + 3 * omdx[28] - 18 * omdx[29] + 6 * omdx[30] +
                        6 * omdx[31] - 9 * omdx[32] + 6 * omdx[33] + 3 * omdx[34] - 18 * omdx[35] +
                        6 * omdx[42] + 6 * omdx[43] + 6 * omdx[44];
        phidx[p45[9]] = +9 * omdx[9] - 18 * omdx[10] + 3 * omdx[11] - 27 * omdx[18] +
                        12 * omdx[19] + 9 * omdx[20] + 9 * omdx[21] - 6 * omdx[22] - 6 * omdx[23] -
                        3 * omdx[36] + 18 * omdx[37] - 6 * omdx[38] + 9 * omdx[39] - 27 * omdx[40] +
                        12 * omdx[41] - 6 * omdx[42] + 18 * omdx[43] - 6 * omdx[44];
        phidx[p45[10]] = -18 * omdx[9] + 84 * omdx[10] - 18 * omdx[11] - 6 * omdx[18] -
                         36 * omdx[19] + 12 * omdx[20] - 30 * omdx[21] + 30 * omdx[22] -
                         12 * omdx[36] + 36 * omdx[37] + 6 * omdx[38] + 12 * omdx[39] -
                         6 * omdx[40] - 36 * omdx[41] - 24 * omdx[42] + 24 * omdx[43];
        phidx[p45[11]] = +3 * omdx[9] - 18 * omdx[10] + 9 * omdx[11] + 6 * omdx[18] -
                         18 * omdx[19] + 3 * omdx[20] + 6 * omdx[21] - 9 * omdx[22] + 6 * omdx[23] -
                         9 * omdx[36] - 12 * omdx[37] + 27 * omdx[38] + 3 * omdx[39] +
                         6 * omdx[40] - 18 * omdx[41] - 18 * omdx[42] + 6 * omdx[43] + 6 * omdx[44];
        phidx[p45[12]] = +9 * omdx[12] - 18 * omdx[13] + 3 * omdx[14] + 9 * omdx[18] -
                         6 * omdx[19] - 6 * omdx[20] - 27 * omdx[21] + 9 * omdx[22] +
                         12 * omdx[23] - 3 * omdx[30] + 18 * omdx[31] - 6 * omdx[32] +
                         9 * omdx[33] - 27 * omdx[34] + 12 * omdx[35] - 6 * omdx[42] -
                         6 * omdx[43] + 18 * omdx[44];
        phidx[p45[13]] = -18 * omdx[12] + 84 * omdx[13] - 18 * omdx[14] - 30 * omdx[18] +
                         30 * omdx[20] - 6 * omdx[21] + 12 * omdx[22] - 36 * omdx[23] -
                         12 * omdx[30] + 36 * omdx[31] + 6 * omdx[32] + 12 * omdx[33] -
                         6 * omdx[34] - 36 * omdx[35] - 24 * omdx[42] + 24 * omdx[44];
        phidx[p45[14]] = +3 * omdx[12] - 18 * omdx[13] + 9 * omdx[14] + 6 * omdx[18] +
                         6 * omdx[19] - 9 * omdx[20] + 6 * omdx[21] + 3 * omdx[22] - 18 * omdx[23] -
                         9 * omdx[30] - 12 * omdx[31] + 27 * omdx[32] + 3 * omdx[33] +
                         6 * omdx[34] - 18 * omdx[35] - 18 * omdx[42] + 6 * omdx[43] + 6 * omdx[44];
        phidx[p45[15]] = +9 * omdx[15] - 18 * omdx[16] + 3 * omdx[17] - 3 * omdx[18] +
                         18 * omdx[19] - 6 * omdx[20] + 9 * omdx[21] - 27 * omdx[22] +
                         12 * omdx[23] - 3 * omdx[24] + 18 * omdx[25] - 6 * omdx[26] +
                         9 * omdx[27] - 27 * omdx[28] + 12 * omdx[29] - 6 * omdx[42] -
                         6 * omdx[43] + 18 * omdx[44];
        phidx[p45[16]] = -18 * omdx[15] + 84 * omdx[16] - 18 * omdx[17] - 12 * omdx[18] +
                         36 * omdx[19] + 6 * omdx[20] + 12 * omdx[21] - 6 * omdx[22] -
                         36 * omdx[23] - 12 * omdx[24] + 36 * omdx[25] + 6 * omdx[26] +
                         12 * omdx[27] - 6 * omdx[28] - 36 * omdx[29] - 24 * omdx[43] +
                         24 * omdx[44];
        phidx[p45[17]] = +3 * omdx[15] - 18 * omdx[16] + 9 * omdx[17] - 9 * omdx[18] -
                         12 * omdx[19] + 27 * omdx[20] + 3 * omdx[21] + 6 * omdx[22] -
                         18 * omdx[23] - 9 * omdx[24] - 12 * omdx[25] + 27 * omdx[26] +
                         3 * omdx[27] + 6 * omdx[28] - 18 * omdx[29] + 6 * omdx[42] -
                         18 * omdx[43] + 6 * omdx[44];
        phidx[p45[18]] = +90 * omdx[18] - 30 * omdx[19] - 45 * omdx[20] - 30 * omdx[21] +
                         15 * omdx[22] + 30 * omdx[23] + 15 * omdx[42] - 45 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[19]] = -30 * omdx[18] + 120 * omdx[19] - 30 * omdx[20] + 30 * omdx[21] -
                         60 * omdx[22] + 30 * omdx[42] - 30 * omdx[43] + 30 * omdx[44];
        phidx[p45[20]] = -45 * omdx[18] - 30 * omdx[19] + 90 * omdx[20] + 15 * omdx[21] +
                         15 * omdx[22] - 60 * omdx[23] + 15 * omdx[42] - 45 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[21]] = -30 * omdx[18] + 30 * omdx[19] + 15 * omdx[20] + 90 * omdx[21] -
                         45 * omdx[22] - 30 * omdx[23] + 15 * omdx[42] + 15 * omdx[43] -
                         45 * omdx[44];
        phidx[p45[22]] = +15 * omdx[18] - 60 * omdx[19] + 15 * omdx[20] - 45 * omdx[21] +
                         90 * omdx[22] - 30 * omdx[23] + 15 * omdx[42] + 15 * omdx[43] -
                         45 * omdx[44];
        phidx[p45[23]] = +30 * omdx[18] - 60 * omdx[20] - 30 * omdx[21] - 30 * omdx[22] +
                         120 * omdx[23] + 30 * omdx[42] + 30 * omdx[43] - 30 * omdx[44];
        phidx[p45[24]] = +90 * omdx[24] - 30 * omdx[25] - 45 * omdx[26] - 30 * omdx[27] +
                         15 * omdx[28] + 30 * omdx[29] + 15 * omdx[42] - 45 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[25]] = -30 * omdx[24] + 120 * omdx[25] - 30 * omdx[26] + 30 * omdx[27] -
                         60 * omdx[28] - 30 * omdx[42] - 30 * omdx[43] + 30 * omdx[44];
        phidx[p45[26]] = -45 * omdx[24] - 30 * omdx[25] + 90 * omdx[26] + 15 * omdx[27] +
                         15 * omdx[28] - 60 * omdx[29] + 15 * omdx[42] - 45 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[27]] = -30 * omdx[24] + 30 * omdx[25] + 15 * omdx[26] + 90 * omdx[27] -
                         45 * omdx[28] - 30 * omdx[29] + 15 * omdx[42] + 15 * omdx[43] -
                         45 * omdx[44];
        phidx[p45[28]] = +15 * omdx[24] - 60 * omdx[25] + 15 * omdx[26] - 45 * omdx[27] +
                         90 * omdx[28] - 30 * omdx[29] + 15 * omdx[42] + 15 * omdx[43] -
                         45 * omdx[44];
        phidx[p45[29]] = +30 * omdx[24] - 60 * omdx[26] - 30 * omdx[27] - 30 * omdx[28] +
                         120 * omdx[29] - 30 * omdx[42] + 30 * omdx[43] - 30 * omdx[44];
        phidx[p45[30]] = +90 * omdx[30] - 30 * omdx[31] - 45 * omdx[32] - 30 * omdx[33] +
                         15 * omdx[34] + 30 * omdx[35] - 45 * omdx[42] + 15 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[31]] = -30 * omdx[30] + 120 * omdx[31] - 30 * omdx[32] + 30 * omdx[33] -
                         60 * omdx[34] - 30 * omdx[42] - 30 * omdx[43] + 30 * omdx[44];
        phidx[p45[32]] = -45 * omdx[30] - 30 * omdx[31] + 90 * omdx[32] + 15 * omdx[33] +
                         15 * omdx[34] - 60 * omdx[35] - 45 * omdx[42] + 15 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[33]] = -30 * omdx[30] + 30 * omdx[31] + 15 * omdx[32] + 90 * omdx[33] -
                         45 * omdx[34] - 30 * omdx[35] + 15 * omdx[42] + 15 * omdx[43] -
                         45 * omdx[44];
        phidx[p45[34]] = +15 * omdx[30] - 60 * omdx[31] + 15 * omdx[32] - 45 * omdx[33] +
                         90 * omdx[34] - 30 * omdx[35] + 15 * omdx[42] + 15 * omdx[43] -
                         45 * omdx[44];
        phidx[p45[35]] = +30 * omdx[30] - 60 * omdx[32] - 30 * omdx[33] - 30 * omdx[34] +
                         120 * omdx[35] + 30 * omdx[42] - 30 * omdx[43] - 30 * omdx[44];
        phidx[p45[36]] = +90 * omdx[36] - 30 * omdx[37] - 45 * omdx[38] - 30 * omdx[39] +
                         15 * omdx[40] + 30 * omdx[41] - 45 * omdx[42] + 15 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[37]] = -30 * omdx[36] + 120 * omdx[37] - 30 * omdx[38] + 30 * omdx[39] -
                         60 * omdx[40] - 30 * omdx[42] + 30 * omdx[43] - 30 * omdx[44];
        phidx[p45[38]] = -45 * omdx[36] - 30 * omdx[37] + 90 * omdx[38] + 15 * omdx[39] +
                         15 * omdx[40] - 60 * omdx[41] - 45 * omdx[42] + 15 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[39]] = -30 * omdx[36] + 30 * omdx[37] + 15 * omdx[38] + 90 * omdx[39] -
                         45 * omdx[40] - 30 * omdx[41] + 15 * omdx[42] - 45 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[40]] = +15 * omdx[36] - 60 * omdx[37] + 15 * omdx[38] - 45 * omdx[39] +
                         90 * omdx[40] - 30 * omdx[41] + 15 * omdx[42] - 45 * omdx[43] +
                         15 * omdx[44];
        phidx[p45[41]] = +30 * omdx[36] - 60 * omdx[38] - 30 * omdx[39] - 30 * omdx[40] +
                         120 * omdx[41] + 30 * omdx[42] - 30 * omdx[43] - 30 * omdx[44];
        phidx[p45[42]] = +90 * omdx[42] - 30 * omdx[43] - 30 * omdx[44];
        phidx[p45[43]] = -30 * omdx[42] + 90 * omdx[43] - 30 * omdx[44];
        phidx[p45[44]] = -30 * omdx[42] - 30 * omdx[43] + 90 * omdx[44];

        for (int k = 0; k < 45; ++k) {
          val(k, 0, op_dx) = phidx[k].x;
          val(k, 1, op_dx) = phidx[k].y;
          val(k, 2, op_dx) = phidx[k].z;
        }
      }

      R3 phidy[45];
      if (whatd & Fop_dy) {
        phidy[p45[0]] = +9 * omdy[0] - 18 * omdy[1] + 3 * omdy[2] - 27 * omdy[30] + 12 * omdy[31] +
                        9 * omdy[32] + 9 * omdy[33] - 6 * omdy[34] - 6 * omdy[35] - 27 * omdy[36] +
                        12 * omdy[37] + 9 * omdy[38] + 9 * omdy[39] - 6 * omdy[40] - 6 * omdy[41] +
                        18 * omdy[42] - 6 * omdy[43] - 6 * omdy[44];
        phidy[p45[1]] = -18 * omdy[0] + 84 * omdy[1] - 18 * omdy[2] - 6 * omdy[30] - 36 * omdy[31] +
                        12 * omdy[32] - 30 * omdy[33] + 30 * omdy[34] - 6 * omdy[36] -
                        36 * omdy[37] + 12 * omdy[38] - 30 * omdy[39] + 30 * omdy[40] +
                        24 * omdy[42];
        phidy[p45[2]] = +3 * omdy[0] - 18 * omdy[1] + 9 * omdy[2] + 6 * omdy[30] - 18 * omdy[31] +
                        3 * omdy[32] + 6 * omdy[33] - 9 * omdy[34] + 6 * omdy[35] + 6 * omdy[36] -
                        18 * omdy[37] + 3 * omdy[38] + 6 * omdy[39] - 9 * omdy[40] + 6 * omdy[41] +
                        6 * omdy[42] + 6 * omdy[43] + 6 * omdy[44];
        phidy[p45[3]] = +9 * omdy[3] - 18 * omdy[4] + 3 * omdy[5] - 27 * omdy[24] + 12 * omdy[25] +
                        9 * omdy[26] + 9 * omdy[27] - 6 * omdy[28] - 6 * omdy[29] + 9 * omdy[36] -
                        6 * omdy[37] - 6 * omdy[38] - 27 * omdy[39] + 9 * omdy[40] + 12 * omdy[41] -
                        6 * omdy[42] + 18 * omdy[43] - 6 * omdy[44];
        phidy[p45[4]] = -18 * omdy[3] + 84 * omdy[4] - 18 * omdy[5] - 6 * omdy[24] - 36 * omdy[25] +
                        12 * omdy[26] - 30 * omdy[27] + 30 * omdy[28] - 30 * omdy[36] +
                        30 * omdy[38] - 6 * omdy[39] + 12 * omdy[40] - 36 * omdy[41] +
                        24 * omdy[43];
        phidy[p45[5]] = +3 * omdy[3] - 18 * omdy[4] + 9 * omdy[5] + 6 * omdy[24] - 18 * omdy[25] +
                        3 * omdy[26] + 6 * omdy[27] - 9 * omdy[28] + 6 * omdy[29] + 6 * omdy[36] +
                        6 * omdy[37] - 9 * omdy[38] + 6 * omdy[39] + 3 * omdy[40] - 18 * omdy[41] +
                        6 * omdy[42] + 6 * omdy[43] + 6 * omdy[44];
        phidy[p45[6]] = +9 * omdy[6] - 18 * omdy[7] + 3 * omdy[8] + 9 * omdy[24] - 6 * omdy[25] -
                        6 * omdy[26] - 27 * omdy[27] + 9 * omdy[28] + 12 * omdy[29] + 9 * omdy[30] -
                        6 * omdy[31] - 6 * omdy[32] - 27 * omdy[33] + 9 * omdy[34] + 12 * omdy[35] -
                        6 * omdy[42] - 6 * omdy[43] + 18 * omdy[44];
        phidy[p45[7]] = -18 * omdy[6] + 84 * omdy[7] - 18 * omdy[8] - 30 * omdy[24] +
                        30 * omdy[26] - 6 * omdy[27] + 12 * omdy[28] - 36 * omdy[29] -
                        30 * omdy[30] + 30 * omdy[32] - 6 * omdy[33] + 12 * omdy[34] -
                        36 * omdy[35] + 24 * omdy[44];
        phidy[p45[8]] = +3 * omdy[6] - 18 * omdy[7] + 9 * omdy[8] + 6 * omdy[24] + 6 * omdy[25] -
                        9 * omdy[26] + 6 * omdy[27] + 3 * omdy[28] - 18 * omdy[29] + 6 * omdy[30] +
                        6 * omdy[31] - 9 * omdy[32] + 6 * omdy[33] + 3 * omdy[34] - 18 * omdy[35] +
                        6 * omdy[42] + 6 * omdy[43] + 6 * omdy[44];
        phidy[p45[9]] = +9 * omdy[9] - 18 * omdy[10] + 3 * omdy[11] - 27 * omdy[18] +
                        12 * omdy[19] + 9 * omdy[20] + 9 * omdy[21] - 6 * omdy[22] - 6 * omdy[23] -
                        3 * omdy[36] + 18 * omdy[37] - 6 * omdy[38] + 9 * omdy[39] - 27 * omdy[40] +
                        12 * omdy[41] - 6 * omdy[42] + 18 * omdy[43] - 6 * omdy[44];
        phidy[p45[10]] = -18 * omdy[9] + 84 * omdy[10] - 18 * omdy[11] - 6 * omdy[18] -
                         36 * omdy[19] + 12 * omdy[20] - 30 * omdy[21] + 30 * omdy[22] -
                         12 * omdy[36] + 36 * omdy[37] + 6 * omdy[38] + 12 * omdy[39] -
                         6 * omdy[40] - 36 * omdy[41] - 24 * omdy[42] + 24 * omdy[43];
        phidy[p45[11]] = +3 * omdy[9] - 18 * omdy[10] + 9 * omdy[11] + 6 * omdy[18] -
                         18 * omdy[19] + 3 * omdy[20] + 6 * omdy[21] - 9 * omdy[22] + 6 * omdy[23] -
                         9 * omdy[36] - 12 * omdy[37] + 27 * omdy[38] + 3 * omdy[39] +
                         6 * omdy[40] - 18 * omdy[41] - 18 * omdy[42] + 6 * omdy[43] + 6 * omdy[44];
        phidy[p45[12]] = +9 * omdy[12] - 18 * omdy[13] + 3 * omdy[14] + 9 * omdy[18] -
                         6 * omdy[19] - 6 * omdy[20] - 27 * omdy[21] + 9 * omdy[22] +
                         12 * omdy[23] - 3 * omdy[30] + 18 * omdy[31] - 6 * omdy[32] +
                         9 * omdy[33] - 27 * omdy[34] + 12 * omdy[35] - 6 * omdy[42] -
                         6 * omdy[43] + 18 * omdy[44];
        phidy[p45[13]] = -18 * omdy[12] + 84 * omdy[13] - 18 * omdy[14] - 30 * omdy[18] +
                         30 * omdy[20] - 6 * omdy[21] + 12 * omdy[22] - 36 * omdy[23] -
                         12 * omdy[30] + 36 * omdy[31] + 6 * omdy[32] + 12 * omdy[33] -
                         6 * omdy[34] - 36 * omdy[35] - 24 * omdy[42] + 24 * omdy[44];
        phidy[p45[14]] = +3 * omdy[12] - 18 * omdy[13] + 9 * omdy[14] + 6 * omdy[18] +
                         6 * omdy[19] - 9 * omdy[20] + 6 * omdy[21] + 3 * omdy[22] - 18 * omdy[23] -
                         9 * omdy[30] - 12 * omdy[31] + 27 * omdy[32] + 3 * omdy[33] +
                         6 * omdy[34] - 18 * omdy[35] - 18 * omdy[42] + 6 * omdy[43] + 6 * omdy[44];
        phidy[p45[15]] = +9 * omdy[15] - 18 * omdy[16] + 3 * omdy[17] - 3 * omdy[18] +
                         18 * omdy[19] - 6 * omdy[20] + 9 * omdy[21] - 27 * omdy[22] +
                         12 * omdy[23] - 3 * omdy[24] + 18 * omdy[25] - 6 * omdy[26] +
                         9 * omdy[27] - 27 * omdy[28] + 12 * omdy[29] - 6 * omdy[42] -
                         6 * omdy[43] + 18 * omdy[44];
        phidy[p45[16]] = -18 * omdy[15] + 84 * omdy[16] - 18 * omdy[17] - 12 * omdy[18] +
                         36 * omdy[19] + 6 * omdy[20] + 12 * omdy[21] - 6 * omdy[22] -
                         36 * omdy[23] - 12 * omdy[24] + 36 * omdy[25] + 6 * omdy[26] +
                         12 * omdy[27] - 6 * omdy[28] - 36 * omdy[29] - 24 * omdy[43] +
                         24 * omdy[44];
        phidy[p45[17]] = +3 * omdy[15] - 18 * omdy[16] + 9 * omdy[17] - 9 * omdy[18] -
                         12 * omdy[19] + 27 * omdy[20] + 3 * omdy[21] + 6 * omdy[22] -
                         18 * omdy[23] - 9 * omdy[24] - 12 * omdy[25] + 27 * omdy[26] +
                         3 * omdy[27] + 6 * omdy[28] - 18 * omdy[29] + 6 * omdy[42] -
                         18 * omdy[43] + 6 * omdy[44];
        phidy[p45[18]] = +90 * omdy[18] - 30 * omdy[19] - 45 * omdy[20] - 30 * omdy[21] +
                         15 * omdy[22] + 30 * omdy[23] + 15 * omdy[42] - 45 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[19]] = -30 * omdy[18] + 120 * omdy[19] - 30 * omdy[20] + 30 * omdy[21] -
                         60 * omdy[22] + 30 * omdy[42] - 30 * omdy[43] + 30 * omdy[44];
        phidy[p45[20]] = -45 * omdy[18] - 30 * omdy[19] + 90 * omdy[20] + 15 * omdy[21] +
                         15 * omdy[22] - 60 * omdy[23] + 15 * omdy[42] - 45 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[21]] = -30 * omdy[18] + 30 * omdy[19] + 15 * omdy[20] + 90 * omdy[21] -
                         45 * omdy[22] - 30 * omdy[23] + 15 * omdy[42] + 15 * omdy[43] -
                         45 * omdy[44];
        phidy[p45[22]] = +15 * omdy[18] - 60 * omdy[19] + 15 * omdy[20] - 45 * omdy[21] +
                         90 * omdy[22] - 30 * omdy[23] + 15 * omdy[42] + 15 * omdy[43] -
                         45 * omdy[44];
        phidy[p45[23]] = +30 * omdy[18] - 60 * omdy[20] - 30 * omdy[21] - 30 * omdy[22] +
                         120 * omdy[23] + 30 * omdy[42] + 30 * omdy[43] - 30 * omdy[44];
        phidy[p45[24]] = +90 * omdy[24] - 30 * omdy[25] - 45 * omdy[26] - 30 * omdy[27] +
                         15 * omdy[28] + 30 * omdy[29] + 15 * omdy[42] - 45 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[25]] = -30 * omdy[24] + 120 * omdy[25] - 30 * omdy[26] + 30 * omdy[27] -
                         60 * omdy[28] - 30 * omdy[42] - 30 * omdy[43] + 30 * omdy[44];
        phidy[p45[26]] = -45 * omdy[24] - 30 * omdy[25] + 90 * omdy[26] + 15 * omdy[27] +
                         15 * omdy[28] - 60 * omdy[29] + 15 * omdy[42] - 45 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[27]] = -30 * omdy[24] + 30 * omdy[25] + 15 * omdy[26] + 90 * omdy[27] -
                         45 * omdy[28] - 30 * omdy[29] + 15 * omdy[42] + 15 * omdy[43] -
                         45 * omdy[44];
        phidy[p45[28]] = +15 * omdy[24] - 60 * omdy[25] + 15 * omdy[26] - 45 * omdy[27] +
                         90 * omdy[28] - 30 * omdy[29] + 15 * omdy[42] + 15 * omdy[43] -
                         45 * omdy[44];
        phidy[p45[29]] = +30 * omdy[24] - 60 * omdy[26] - 30 * omdy[27] - 30 * omdy[28] +
                         120 * omdy[29] - 30 * omdy[42] + 30 * omdy[43] - 30 * omdy[44];
        phidy[p45[30]] = +90 * omdy[30] - 30 * omdy[31] - 45 * omdy[32] - 30 * omdy[33] +
                         15 * omdy[34] + 30 * omdy[35] - 45 * omdy[42] + 15 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[31]] = -30 * omdy[30] + 120 * omdy[31] - 30 * omdy[32] + 30 * omdy[33] -
                         60 * omdy[34] - 30 * omdy[42] - 30 * omdy[43] + 30 * omdy[44];
        phidy[p45[32]] = -45 * omdy[30] - 30 * omdy[31] + 90 * omdy[32] + 15 * omdy[33] +
                         15 * omdy[34] - 60 * omdy[35] - 45 * omdy[42] + 15 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[33]] = -30 * omdy[30] + 30 * omdy[31] + 15 * omdy[32] + 90 * omdy[33] -
                         45 * omdy[34] - 30 * omdy[35] + 15 * omdy[42] + 15 * omdy[43] -
                         45 * omdy[44];
        phidy[p45[34]] = +15 * omdy[30] - 60 * omdy[31] + 15 * omdy[32] - 45 * omdy[33] +
                         90 * omdy[34] - 30 * omdy[35] + 15 * omdy[42] + 15 * omdy[43] -
                         45 * omdy[44];
        phidy[p45[35]] = +30 * omdy[30] - 60 * omdy[32] - 30 * omdy[33] - 30 * omdy[34] +
                         120 * omdy[35] + 30 * omdy[42] - 30 * omdy[43] - 30 * omdy[44];
        phidy[p45[36]] = +90 * omdy[36] - 30 * omdy[37] - 45 * omdy[38] - 30 * omdy[39] +
                         15 * omdy[40] + 30 * omdy[41] - 45 * omdy[42] + 15 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[37]] = -30 * omdy[36] + 120 * omdy[37] - 30 * omdy[38] + 30 * omdy[39] -
                         60 * omdy[40] - 30 * omdy[42] + 30 * omdy[43] - 30 * omdy[44];
        phidy[p45[38]] = -45 * omdy[36] - 30 * omdy[37] + 90 * omdy[38] + 15 * omdy[39] +
                         15 * omdy[40] - 60 * omdy[41] - 45 * omdy[42] + 15 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[39]] = -30 * omdy[36] + 30 * omdy[37] + 15 * omdy[38] + 90 * omdy[39] -
                         45 * omdy[40] - 30 * omdy[41] + 15 * omdy[42] - 45 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[40]] = +15 * omdy[36] - 60 * omdy[37] + 15 * omdy[38] - 45 * omdy[39] +
                         90 * omdy[40] - 30 * omdy[41] + 15 * omdy[42] - 45 * omdy[43] +
                         15 * omdy[44];
        phidy[p45[41]] = +30 * omdy[36] - 60 * omdy[38] - 30 * omdy[39] - 30 * omdy[40] +
                         120 * omdy[41] + 30 * omdy[42] - 30 * omdy[43] - 30 * omdy[44];
        phidy[p45[42]] = +90 * omdy[42] - 30 * omdy[43] - 30 * omdy[44];
        phidy[p45[43]] = -30 * omdy[42] + 90 * omdy[43] - 30 * omdy[44];
        phidy[p45[44]] = -30 * omdy[42] - 30 * omdy[43] + 90 * omdy[44];

        for (int k = 0; k < 45; ++k) {
          val(k, 0, op_dy) = phidy[k].x;
          val(k, 1, op_dy) = phidy[k].y;
          val(k, 2, op_dy) = phidy[k].z;
        }
      }

      R3 phidz[45];
      if (whatd & Fop_dz) {
        phidz[p45[0]] = +9 * omdz[0] - 18 * omdz[1] + 3 * omdz[2] - 27 * omdz[30] + 12 * omdz[31] +
                        9 * omdz[32] + 9 * omdz[33] - 6 * omdz[34] - 6 * omdz[35] - 27 * omdz[36] +
                        12 * omdz[37] + 9 * omdz[38] + 9 * omdz[39] - 6 * omdz[40] - 6 * omdz[41] +
                        18 * omdz[42] - 6 * omdz[43] - 6 * omdz[44];
        phidz[p45[1]] = -18 * omdz[0] + 84 * omdz[1] - 18 * omdz[2] - 6 * omdz[30] - 36 * omdz[31] +
                        12 * omdz[32] - 30 * omdz[33] + 30 * omdz[34] - 6 * omdz[36] -
                        36 * omdz[37] + 12 * omdz[38] - 30 * omdz[39] + 30 * omdz[40] +
                        24 * omdz[42];
        phidz[p45[2]] = +3 * omdz[0] - 18 * omdz[1] + 9 * omdz[2] + 6 * omdz[30] - 18 * omdz[31] +
                        3 * omdz[32] + 6 * omdz[33] - 9 * omdz[34] + 6 * omdz[35] + 6 * omdz[36] -
                        18 * omdz[37] + 3 * omdz[38] + 6 * omdz[39] - 9 * omdz[40] + 6 * omdz[41] +
                        6 * omdz[42] + 6 * omdz[43] + 6 * omdz[44];
        phidz[p45[3]] = +9 * omdz[3] - 18 * omdz[4] + 3 * omdz[5] - 27 * omdz[24] + 12 * omdz[25] +
                        9 * omdz[26] + 9 * omdz[27] - 6 * omdz[28] - 6 * omdz[29] + 9 * omdz[36] -
                        6 * omdz[37] - 6 * omdz[38] - 27 * omdz[39] + 9 * omdz[40] + 12 * omdz[41] -
                        6 * omdz[42] + 18 * omdz[43] - 6 * omdz[44];
        phidz[p45[4]] = -18 * omdz[3] + 84 * omdz[4] - 18 * omdz[5] - 6 * omdz[24] - 36 * omdz[25] +
                        12 * omdz[26] - 30 * omdz[27] + 30 * omdz[28] - 30 * omdz[36] +
                        30 * omdz[38] - 6 * omdz[39] + 12 * omdz[40] - 36 * omdz[41] +
                        24 * omdz[43];
        phidz[p45[5]] = +3 * omdz[3] - 18 * omdz[4] + 9 * omdz[5] + 6 * omdz[24] - 18 * omdz[25] +
                        3 * omdz[26] + 6 * omdz[27] - 9 * omdz[28] + 6 * omdz[29] + 6 * omdz[36] +
                        6 * omdz[37] - 9 * omdz[38] + 6 * omdz[39] + 3 * omdz[40] - 18 * omdz[41] +
                        6 * omdz[42] + 6 * omdz[43] + 6 * omdz[44];
        phidz[p45[6]] = +9 * omdz[6] - 18 * omdz[7] + 3 * omdz[8] + 9 * omdz[24] - 6 * omdz[25] -
                        6 * omdz[26] - 27 * omdz[27] + 9 * omdz[28] + 12 * omdz[29] + 9 * omdz[30] -
                        6 * omdz[31] - 6 * omdz[32] - 27 * omdz[33] + 9 * omdz[34] + 12 * omdz[35] -
                        6 * omdz[42] - 6 * omdz[43] + 18 * omdz[44];
        phidz[p45[7]] = -18 * omdz[6] + 84 * omdz[7] - 18 * omdz[8] - 30 * omdz[24] +
                        30 * omdz[26] - 6 * omdz[27] + 12 * omdz[28] - 36 * omdz[29] -
                        30 * omdz[30] + 30 * omdz[32] - 6 * omdz[33] + 12 * omdz[34] -
                        36 * omdz[35] + 24 * omdz[44];
        phidz[p45[8]] = +3 * omdz[6] - 18 * omdz[7] + 9 * omdz[8] + 6 * omdz[24] + 6 * omdz[25] -
                        9 * omdz[26] + 6 * omdz[27] + 3 * omdz[28] - 18 * omdz[29] + 6 * omdz[30] +
                        6 * omdz[31] - 9 * omdz[32] + 6 * omdz[33] + 3 * omdz[34] - 18 * omdz[35] +
                        6 * omdz[42] + 6 * omdz[43] + 6 * omdz[44];
        phidz[p45[9]] = +9 * omdz[9] - 18 * omdz[10] + 3 * omdz[11] - 27 * omdz[18] +
                        12 * omdz[19] + 9 * omdz[20] + 9 * omdz[21] - 6 * omdz[22] - 6 * omdz[23] -
                        3 * omdz[36] + 18 * omdz[37] - 6 * omdz[38] + 9 * omdz[39] - 27 * omdz[40] +
                        12 * omdz[41] - 6 * omdz[42] + 18 * omdz[43] - 6 * omdz[44];
        phidz[p45[10]] = -18 * omdz[9] + 84 * omdz[10] - 18 * omdz[11] - 6 * omdz[18] -
                         36 * omdz[19] + 12 * omdz[20] - 30 * omdz[21] + 30 * omdz[22] -
                         12 * omdz[36] + 36 * omdz[37] + 6 * omdz[38] + 12 * omdz[39] -
                         6 * omdz[40] - 36 * omdz[41] - 24 * omdz[42] + 24 * omdz[43];
        phidz[p45[11]] = +3 * omdz[9] - 18 * omdz[10] + 9 * omdz[11] + 6 * omdz[18] -
                         18 * omdz[19] + 3 * omdz[20] + 6 * omdz[21] - 9 * omdz[22] + 6 * omdz[23] -
                         9 * omdz[36] - 12 * omdz[37] + 27 * omdz[38] + 3 * omdz[39] +
                         6 * omdz[40] - 18 * omdz[41] - 18 * omdz[42] + 6 * omdz[43] + 6 * omdz[44];
        phidz[p45[12]] = +9 * omdz[12] - 18 * omdz[13] + 3 * omdz[14] + 9 * omdz[18] -
                         6 * omdz[19] - 6 * omdz[20] - 27 * omdz[21] + 9 * omdz[22] +
                         12 * omdz[23] - 3 * omdz[30] + 18 * omdz[31] - 6 * omdz[32] +
                         9 * omdz[33] - 27 * omdz[34] + 12 * omdz[35] - 6 * omdz[42] -
                         6 * omdz[43] + 18 * omdz[44];
        phidz[p45[13]] = -18 * omdz[12] + 84 * omdz[13] - 18 * omdz[14] - 30 * omdz[18] +
                         30 * omdz[20] - 6 * omdz[21] + 12 * omdz[22] - 36 * omdz[23] -
                         12 * omdz[30] + 36 * omdz[31] + 6 * omdz[32] + 12 * omdz[33] -
                         6 * omdz[34] - 36 * omdz[35] - 24 * omdz[42] + 24 * omdz[44];
        phidz[p45[14]] = +3 * omdz[12] - 18 * omdz[13] + 9 * omdz[14] + 6 * omdz[18] +
                         6 * omdz[19] - 9 * omdz[20] + 6 * omdz[21] + 3 * omdz[22] - 18 * omdz[23] -
                         9 * omdz[30] - 12 * omdz[31] + 27 * omdz[32] + 3 * omdz[33] +
                         6 * omdz[34] - 18 * omdz[35] - 18 * omdz[42] + 6 * omdz[43] + 6 * omdz[44];
        phidz[p45[15]] = +9 * omdz[15] - 18 * omdz[16] + 3 * omdz[17] - 3 * omdz[18] +
                         18 * omdz[19] - 6 * omdz[20] + 9 * omdz[21] - 27 * omdz[22] +
                         12 * omdz[23] - 3 * omdz[24] + 18 * omdz[25] - 6 * omdz[26] +
                         9 * omdz[27] - 27 * omdz[28] + 12 * omdz[29] - 6 * omdz[42] -
                         6 * omdz[43] + 18 * omdz[44];
        phidz[p45[16]] = -18 * omdz[15] + 84 * omdz[16] - 18 * omdz[17] - 12 * omdz[18] +
                         36 * omdz[19] + 6 * omdz[20] + 12 * omdz[21] - 6 * omdz[22] -
                         36 * omdz[23] - 12 * omdz[24] + 36 * omdz[25] + 6 * omdz[26] +
                         12 * omdz[27] - 6 * omdz[28] - 36 * omdz[29] - 24 * omdz[43] +
                         24 * omdz[44];
        phidz[p45[17]] = +3 * omdz[15] - 18 * omdz[16] + 9 * omdz[17] - 9 * omdz[18] -
                         12 * omdz[19] + 27 * omdz[20] + 3 * omdz[21] + 6 * omdz[22] -
                         18 * omdz[23] - 9 * omdz[24] - 12 * omdz[25] + 27 * omdz[26] +
                         3 * omdz[27] + 6 * omdz[28] - 18 * omdz[29] + 6 * omdz[42] -
                         18 * omdz[43] + 6 * omdz[44];
        phidz[p45[18]] = +90 * omdz[18] - 30 * omdz[19] - 45 * omdz[20] - 30 * omdz[21] +
                         15 * omdz[22] + 30 * omdz[23] + 15 * omdz[42] - 45 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[19]] = -30 * omdz[18] + 120 * omdz[19] - 30 * omdz[20] + 30 * omdz[21] -
                         60 * omdz[22] + 30 * omdz[42] - 30 * omdz[43] + 30 * omdz[44];
        phidz[p45[20]] = -45 * omdz[18] - 30 * omdz[19] + 90 * omdz[20] + 15 * omdz[21] +
                         15 * omdz[22] - 60 * omdz[23] + 15 * omdz[42] - 45 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[21]] = -30 * omdz[18] + 30 * omdz[19] + 15 * omdz[20] + 90 * omdz[21] -
                         45 * omdz[22] - 30 * omdz[23] + 15 * omdz[42] + 15 * omdz[43] -
                         45 * omdz[44];
        phidz[p45[22]] = +15 * omdz[18] - 60 * omdz[19] + 15 * omdz[20] - 45 * omdz[21] +
                         90 * omdz[22] - 30 * omdz[23] + 15 * omdz[42] + 15 * omdz[43] -
                         45 * omdz[44];
        phidz[p45[23]] = +30 * omdz[18] - 60 * omdz[20] - 30 * omdz[21] - 30 * omdz[22] +
                         120 * omdz[23] + 30 * omdz[42] + 30 * omdz[43] - 30 * omdz[44];
        phidz[p45[24]] = +90 * omdz[24] - 30 * omdz[25] - 45 * omdz[26] - 30 * omdz[27] +
                         15 * omdz[28] + 30 * omdz[29] + 15 * omdz[42] - 45 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[25]] = -30 * omdz[24] + 120 * omdz[25] - 30 * omdz[26] + 30 * omdz[27] -
                         60 * omdz[28] - 30 * omdz[42] - 30 * omdz[43] + 30 * omdz[44];
        phidz[p45[26]] = -45 * omdz[24] - 30 * omdz[25] + 90 * omdz[26] + 15 * omdz[27] +
                         15 * omdz[28] - 60 * omdz[29] + 15 * omdz[42] - 45 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[27]] = -30 * omdz[24] + 30 * omdz[25] + 15 * omdz[26] + 90 * omdz[27] -
                         45 * omdz[28] - 30 * omdz[29] + 15 * omdz[42] + 15 * omdz[43] -
                         45 * omdz[44];
        phidz[p45[28]] = +15 * omdz[24] - 60 * omdz[25] + 15 * omdz[26] - 45 * omdz[27] +
                         90 * omdz[28] - 30 * omdz[29] + 15 * omdz[42] + 15 * omdz[43] -
                         45 * omdz[44];
        phidz[p45[29]] = +30 * omdz[24] - 60 * omdz[26] - 30 * omdz[27] - 30 * omdz[28] +
                         120 * omdz[29] - 30 * omdz[42] + 30 * omdz[43] - 30 * omdz[44];
        phidz[p45[30]] = +90 * omdz[30] - 30 * omdz[31] - 45 * omdz[32] - 30 * omdz[33] +
                         15 * omdz[34] + 30 * omdz[35] - 45 * omdz[42] + 15 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[31]] = -30 * omdz[30] + 120 * omdz[31] - 30 * omdz[32] + 30 * omdz[33] -
                         60 * omdz[34] - 30 * omdz[42] - 30 * omdz[43] + 30 * omdz[44];
        phidz[p45[32]] = -45 * omdz[30] - 30 * omdz[31] + 90 * omdz[32] + 15 * omdz[33] +
                         15 * omdz[34] - 60 * omdz[35] - 45 * omdz[42] + 15 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[33]] = -30 * omdz[30] + 30 * omdz[31] + 15 * omdz[32] + 90 * omdz[33] -
                         45 * omdz[34] - 30 * omdz[35] + 15 * omdz[42] + 15 * omdz[43] -
                         45 * omdz[44];
        phidz[p45[34]] = +15 * omdz[30] - 60 * omdz[31] + 15 * omdz[32] - 45 * omdz[33] +
                         90 * omdz[34] - 30 * omdz[35] + 15 * omdz[42] + 15 * omdz[43] -
                         45 * omdz[44];
        phidz[p45[35]] = +30 * omdz[30] - 60 * omdz[32] - 30 * omdz[33] - 30 * omdz[34] +
                         120 * omdz[35] + 30 * omdz[42] - 30 * omdz[43] - 30 * omdz[44];
        phidz[p45[36]] = +90 * omdz[36] - 30 * omdz[37] - 45 * omdz[38] - 30 * omdz[39] +
                         15 * omdz[40] + 30 * omdz[41] - 45 * omdz[42] + 15 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[37]] = -30 * omdz[36] + 120 * omdz[37] - 30 * omdz[38] + 30 * omdz[39] -
                         60 * omdz[40] - 30 * omdz[42] + 30 * omdz[43] - 30 * omdz[44];
        phidz[p45[38]] = -45 * omdz[36] - 30 * omdz[37] + 90 * omdz[38] + 15 * omdz[39] +
                         15 * omdz[40] - 60 * omdz[41] - 45 * omdz[42] + 15 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[39]] = -30 * omdz[36] + 30 * omdz[37] + 15 * omdz[38] + 90 * omdz[39] -
                         45 * omdz[40] - 30 * omdz[41] + 15 * omdz[42] - 45 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[40]] = +15 * omdz[36] - 60 * omdz[37] + 15 * omdz[38] - 45 * omdz[39] +
                         90 * omdz[40] - 30 * omdz[41] + 15 * omdz[42] - 45 * omdz[43] +
                         15 * omdz[44];
        phidz[p45[41]] = +30 * omdz[36] - 60 * omdz[38] - 30 * omdz[39] - 30 * omdz[40] +
                         120 * omdz[41] + 30 * omdz[42] - 30 * omdz[43] - 30 * omdz[44];
        phidz[p45[42]] = +90 * omdz[42] - 30 * omdz[43] - 30 * omdz[44];
        phidz[p45[43]] = -30 * omdz[42] + 90 * omdz[43] - 30 * omdz[44];
        phidz[p45[44]] = -30 * omdz[42] - 30 * omdz[43] + 90 * omdz[44];

        for (int k = 0; k < 45; ++k) {
          val(k, 0, op_dz) = phidz[k].x;
          val(k, 1, op_dz) = phidz[k].y;
          val(k, 2, op_dz) = phidz[k].z;
        }
      }
    }

    if (whatd & Fop_D2) {    // SECOND DERIVATIVES
      R3 omdxx[45];
      R3 omdyy[45];
      R3 omdzz[45];
      R3 omdxy[45];
      R3 omdxz[45];
      R3 omdyz[45];

      // 18 edge functions (3 for each edge):
      for (int i = 0; i < Element::ne; ++i) {
        int ii0 = Element::nvedge[i][0], ii1 = Element::nvedge[i][1];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        if (whatd & Fop_dxx) {
          omdxx[i * 3] = 2 * D[i0].x * D[i0].x * (l[i0] * D[i1] - l[i1] * D[i0]) +
                         4 * D[i0].x * l[i0] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxx[i * 3 + 1] =
            2 * D[i0].x * D[i1].x * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i0].x * l[i1] + l[i0] * D[i1].x) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxx[i * 3 + 2] = 2 * D[i1].x * D[i1].x * (l[i0] * D[i1] - l[i1] * D[i0]) +
                             4 * D[i1].x * l[i1] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
        }

        if (whatd & Fop_dyy) {
          omdyy[i * 3] = 2 * D[i0].y * D[i0].y * (l[i0] * D[i1] - l[i1] * D[i0]) +
                         4 * D[i0].y * l[i0] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyy[i * 3 + 1] =
            2 * D[i0].y * D[i1].y * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i0].y * l[i1] + l[i0] * D[i1].y) * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyy[i * 3 + 2] = 2 * D[i1].y * D[i1].y * (l[i0] * D[i1] - l[i1] * D[i0]) +
                             4 * D[i1].y * l[i1] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
        }

        if (whatd & Fop_dzz) {
          omdzz[i * 3] = 2 * D[i0].z * D[i0].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
                         4 * D[i0].z * l[i0] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdzz[i * 3 + 1] =
            2 * D[i0].z * D[i1].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i0].z * l[i1] + l[i0] * D[i1].z) * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdzz[i * 3 + 2] = 2 * D[i1].z * D[i1].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
                             4 * D[i1].z * l[i1] * (D[i0].z * D[i1] - D[i1].z * D[i0]);
        }

        if (whatd & Fop_dxy) {
          omdxy[i * 3] = 2 * D[i0].x * D[i0].y * (l[i0] * D[i1] - l[i1] * D[i0]) +
                         2 * D[i0].x * l[i0] * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
                         2 * D[i0].y * l[i0] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxy[i * 3 + 1] =
            (D[i0].x * D[i1].y + D[i0].y * D[i1].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i0].x * l[i1] + l[i0] * D[i1].x) * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
            (D[i0].y * l[i1] + l[i0] * D[i1].y) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxy[i * 3 + 2] = 2 * D[i1].x * D[i1].y * (l[i0] * D[i1] - l[i1] * D[i0]) +
                             2 * D[i1].x * l[i1] * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
                             2 * D[i1].y * l[i1] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
        }

        if (whatd & Fop_dxz) {
          omdxz[i * 3] = 2 * D[i0].x * D[i0].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
                         2 * D[i0].x * l[i0] * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                         2 * D[i0].z * l[i0] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxz[i * 3 + 1] =
            (D[i0].x * D[i1].z + D[i0].z * D[i1].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i0].x * l[i1] + l[i0] * D[i1].x) * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
            (D[i0].z * l[i1] + l[i0] * D[i1].z) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxz[i * 3 + 2] = 2 * D[i1].x * D[i1].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
                             2 * D[i1].x * l[i1] * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                             2 * D[i1].z * l[i1] * (D[i0].x * D[i1] - D[i1].x * D[i0]);
        }

        if (whatd & Fop_dyz) {
          omdyz[i * 3] = 2 * D[i0].y * D[i0].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
                         2 * D[i0].y * l[i0] * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                         2 * D[i0].z * l[i0] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyz[i * 3 + 1] =
            (D[i0].y * D[i1].z + D[i0].z * D[i1].y) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i0].y * l[i1] + l[i0] * D[i1].y) * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
            (D[i0].z * l[i1] + l[i0] * D[i1].z) * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyz[i * 3 + 2] = 2 * D[i1].y * D[i1].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
                             2 * D[i1].y * l[i1] * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
                             2 * D[i1].z * l[i1] * (D[i0].y * D[i1] - D[i1].y * D[i0]);
        }
      }

      // 24 face functions (6 for each face):
      for (int j = 0; j < Element::nf; ++j) {
        int ii0 = mynvface[j][0];
        int ii1 = mynvface[j][1];
        int ii2 = mynvface[j][2];
        int i0 = perm[ii0];
        int i1 = perm[ii1];
        int i2 = perm[ii2];
        if (whatd & Fop_dxx) {
          omdxx[18 + j * 6] =
            2 * D[i0].x * D[i2].x * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i0].x * l[i2] + l[i0] * D[i2].x) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxx[18 + j * 6 + 1] =
            2 * D[i1].x * D[i2].x * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i1].x * l[i2] + l[i1] * D[i2].x) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxx[18 + j * 6 + 2] =
            2 * D[i2].x * D[i2].x * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i2].x * l[i2] + l[i2] * D[i2].x) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxx[18 + j * 6 + 3] =
            2 * D[i0].x * D[i1].x * (l[i0] * D[i2] - l[i2] * D[i0]) +
            2 * (D[i0].x * l[i1] + l[i0] * D[i1].x) * (D[i0].x * D[i2] - D[i2].x * D[i0]);
          omdxx[18 + j * 6 + 4] =
            2 * D[i1].x * D[i1].x * (l[i0] * D[i2] - l[i2] * D[i0]) +
            2 * (D[i1].x * l[i1] + l[i1] * D[i1].x) * (D[i0].x * D[i2] - D[i2].x * D[i0]);
          omdxx[18 + j * 6 + 5] =
            2 * D[i2].x * D[i1].x * (l[i0] * D[i2] - l[i2] * D[i0]) +
            2 * (D[i2].x * l[i1] + l[i2] * D[i1].x) * (D[i0].x * D[i2] - D[i2].x * D[i0]);
        }

        if (whatd & Fop_dyy) {
          omdyy[18 + j * 6] =
            2 * D[i0].y * D[i2].y * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i0].y * l[i2] + l[i0] * D[i2].y) * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyy[18 + j * 6 + 1] =
            2 * D[i1].y * D[i2].y * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i1].y * l[i2] + l[i1] * D[i2].y) * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyy[18 + j * 6 + 2] =
            2 * D[i2].y * D[i2].y * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i2].y * l[i2] + l[i2] * D[i2].y) * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyy[18 + j * 6 + 3] =
            2 * D[i0].y * D[i1].y * (l[i0] * D[i2] - l[i2] * D[i0]) +
            2 * (D[i0].y * l[i1] + l[i0] * D[i1].y) * (D[i0].y * D[i2] - D[i2].y * D[i0]);
          omdyy[18 + j * 6 + 4] =
            2 * D[i1].y * D[i1].y * (l[i0] * D[i2] - l[i2] * D[i0]) +
            2 * (D[i1].y * l[i1] + l[i1] * D[i1].y) * (D[i0].y * D[i2] - D[i2].y * D[i0]);
          omdyy[18 + j * 6 + 5] =
            2 * D[i2].y * D[i1].y * (l[i0] * D[i2] - l[i2] * D[i0]) +
            2 * (D[i2].y * l[i1] + l[i2] * D[i1].y) * (D[i0].y * D[i2] - D[i2].y * D[i0]);
        }

        if (whatd & Fop_dzz) {
          omdzz[18 + j * 6] =
            2 * D[i0].z * D[i2].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i0].z * l[i2] + l[i0] * D[i2].z) * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdzz[18 + j * 6 + 1] =
            2 * D[i1].z * D[i2].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i1].z * l[i2] + l[i1] * D[i2].z) * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdzz[18 + j * 6 + 2] =
            2 * D[i2].z * D[i2].z * (l[i0] * D[i1] - l[i1] * D[i0]) +
            2 * (D[i2].z * l[i2] + l[i2] * D[i2].z) * (D[i0].z * D[i1] - D[i1].z * D[i0]);
          omdzz[18 + j * 6 + 3] =
            2 * D[i0].z * D[i1].z * (l[i0] * D[i2] - l[i2] * D[i0]) +
            2 * (D[i0].z * l[i1] + l[i0] * D[i1].z) * (D[i0].z * D[i2] - D[i2].z * D[i0]);
          omdzz[18 + j * 6 + 4] =
            2 * D[i1].z * D[i1].z * (l[i0] * D[i2] - l[i2] * D[i0]) +
            2 * (D[i1].z * l[i1] + l[i1] * D[i1].z) * (D[i0].z * D[i2] - D[i2].z * D[i0]);
          omdzz[18 + j * 6 + 5] =
            2 * D[i2].z * D[i1].z * (l[i0] * D[i2] - l[i2] * D[i0]) +
            2 * (D[i2].z * l[i1] + l[i2] * D[i1].z) * (D[i0].z * D[i2] - D[i2].z * D[i0]);
        }

        if (whatd & Fop_dxy) {
          omdxy[18 + j * 6] =
            (D[i0].x * D[i2].y + D[i0].y * D[i2].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i0].x * l[i2] + l[i0] * D[i2].x) * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
            (D[i0].y * l[i2] + l[i0] * D[i2].y) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxy[18 + j * 6 + 1] =
            (D[i1].x * D[i2].y + D[i1].y * D[i2].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i1].x * l[i2] + l[i1] * D[i2].x) * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
            (D[i1].y * l[i2] + l[i1] * D[i2].y) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxy[18 + j * 6 + 2] =
            (D[i2].x * D[i2].y + D[i2].y * D[i2].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i2].x * l[i2] + l[i2] * D[i2].x) * (D[i0].y * D[i1] - D[i1].y * D[i0]) +
            (D[i2].y * l[i2] + l[i2] * D[i2].y) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxy[18 + j * 6 + 3] =
            (D[i0].x * D[i1].y + D[i0].y * D[i1].x) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            (D[i0].x * l[i1] + l[i0] * D[i1].x) * (D[i0].y * D[i2] - D[i2].y * D[i0]) +
            (D[i0].y * l[i1] + l[i0] * D[i1].y) * (D[i0].x * D[i2] - D[i2].x * D[i0]);
          omdxy[18 + j * 6 + 4] =
            (D[i1].x * D[i1].y + D[i1].y * D[i1].x) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            (D[i1].x * l[i1] + l[i1] * D[i1].x) * (D[i0].y * D[i2] - D[i2].y * D[i0]) +
            (D[i1].y * l[i1] + l[i1] * D[i1].y) * (D[i0].x * D[i2] - D[i2].x * D[i0]);
          omdxy[18 + j * 6 + 5] =
            (D[i2].x * D[i1].y + D[i2].y * D[i1].x) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            (D[i2].x * l[i1] + l[i2] * D[i1].x) * (D[i0].y * D[i2] - D[i2].y * D[i0]) +
            (D[i2].y * l[i1] + l[i2] * D[i1].y) * (D[i0].x * D[i2] - D[i2].x * D[i0]);
        }

        if (whatd & Fop_dxz) {
          omdxz[18 + j * 6] =
            (D[i0].x * D[i2].z + D[i0].z * D[i2].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i0].x * l[i2] + l[i0] * D[i2].x) * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
            (D[i0].z * l[i2] + l[i0] * D[i2].z) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxz[18 + j * 6 + 1] =
            (D[i1].x * D[i2].z + D[i1].z * D[i2].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i1].x * l[i2] + l[i1] * D[i2].x) * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
            (D[i1].z * l[i2] + l[i1] * D[i2].z) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxz[18 + j * 6 + 2] =
            (D[i2].x * D[i2].z + D[i2].z * D[i2].x) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i2].x * l[i2] + l[i2] * D[i2].x) * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
            (D[i2].z * l[i2] + l[i2] * D[i2].z) * (D[i0].x * D[i1] - D[i1].x * D[i0]);
          omdxz[18 + j * 6 + 3] =
            (D[i0].x * D[i1].z + D[i0].z * D[i1].x) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            (D[i0].x * l[i1] + l[i0] * D[i1].x) * (D[i0].z * D[i2] - D[i2].z * D[i0]) +
            (D[i0].z * l[i1] + l[i0] * D[i1].z) * (D[i0].x * D[i2] - D[i2].x * D[i0]);
          omdxz[18 + j * 6 + 4] =
            (D[i1].x * D[i1].z + D[i1].z * D[i1].x) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            (D[i1].x * l[i1] + l[i1] * D[i1].x) * (D[i0].z * D[i2] - D[i2].z * D[i0]) +
            (D[i1].z * l[i1] + l[i1] * D[i1].z) * (D[i0].x * D[i2] - D[i2].x * D[i0]);
          omdxz[18 + j * 6 + 5] =
            (D[i2].x * D[i1].z + D[i2].z * D[i1].x) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            (D[i2].x * l[i1] + l[i2] * D[i1].x) * (D[i0].z * D[i2] - D[i2].z * D[i0]) +
            (D[i2].z * l[i1] + l[i2] * D[i1].z) * (D[i0].x * D[i2] - D[i2].x * D[i0]);
        }

        if (whatd & Fop_dyz) {
          omdyz[18 + j * 6] =
            (D[i0].y * D[i2].z + D[i0].z * D[i2].y) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i0].y * l[i2] + l[i0] * D[i2].y) * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
            (D[i0].z * l[i2] + l[i0] * D[i2].z) * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyz[18 + j * 6 + 1] =
            (D[i1].y * D[i2].z + D[i1].z * D[i2].y) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i1].y * l[i2] + l[i1] * D[i2].y) * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
            (D[i1].z * l[i2] + l[i1] * D[i2].z) * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyz[18 + j * 6 + 2] =
            (D[i2].y * D[i2].z + D[i2].z * D[i2].y) * (l[i0] * D[i1] - l[i1] * D[i0]) +
            (D[i2].y * l[i2] + l[i2] * D[i2].y) * (D[i0].z * D[i1] - D[i1].z * D[i0]) +
            (D[i2].z * l[i2] + l[i2] * D[i2].z) * (D[i0].y * D[i1] - D[i1].y * D[i0]);
          omdyz[18 + j * 6 + 3] =
            (D[i0].y * D[i1].z + D[i0].z * D[i1].y) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            (D[i0].y * l[i1] + l[i0] * D[i1].y) * (D[i0].z * D[i2] - D[i2].z * D[i0]) +
            (D[i0].z * l[i1] + l[i0] * D[i1].z) * (D[i0].y * D[i2] - D[i2].y * D[i0]);
          omdyz[18 + j * 6 + 4] =
            (D[i1].y * D[i1].z + D[i1].z * D[i1].y) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            (D[i1].y * l[i1] + l[i1] * D[i1].y) * (D[i0].z * D[i2] - D[i2].z * D[i0]) +
            (D[i1].z * l[i1] + l[i1] * D[i1].z) * (D[i0].y * D[i2] - D[i2].y * D[i0]);
          omdyz[18 + j * 6 + 5] =
            (D[i2].y * D[i1].z + D[i2].z * D[i1].y) * (l[i0] * D[i2] - l[i2] * D[i0]) +
            (D[i2].y * l[i1] + l[i2] * D[i1].y) * (D[i0].z * D[i2] - D[i2].z * D[i0]) +
            (D[i2].z * l[i1] + l[i2] * D[i1].z) * (D[i0].y * D[i2] - D[i2].y * D[i0]);
        }
      }

      // 3 volume functions
      if (whatd & Fop_dxx) {
        omdxx[42] =
          2 * D[perm[2]].x * D[perm[3]].x * (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]) +
          2 * (D[perm[2]].x * l[perm[3]] + l[perm[2]] * D[perm[3]].x) *
            (D[perm[0]].x * D[perm[1]] - D[perm[1]].x * D[perm[0]]);
        omdxx[43] =
          2 * D[perm[1]].x * D[perm[3]].x * (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]) +
          2 * (D[perm[1]].x * l[perm[3]] + l[perm[1]] * D[perm[3]].x) *
            (D[perm[0]].x * D[perm[2]] - D[perm[2]].x * D[perm[0]]);
        omdxx[44] =
          2 * D[perm[1]].x * D[perm[2]].x * (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]) +
          2 * (D[perm[1]].x * l[perm[2]] + l[perm[1]] * D[perm[2]].x) *
            (D[perm[0]].x * D[perm[3]] - D[perm[3]].x * D[perm[0]]);
      }

      if (whatd & Fop_dyy) {
        omdyy[42] =
          2 * D[perm[2]].y * D[perm[3]].y * (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]) +
          2 * (D[perm[2]].y * l[perm[3]] + l[perm[2]] * D[perm[3]].y) *
            (D[perm[0]].y * D[perm[1]] - D[perm[1]].y * D[perm[0]]);
        omdyy[43] =
          2 * D[perm[1]].y * D[perm[3]].y * (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]) +
          2 * (D[perm[1]].y * l[perm[3]] + l[perm[1]] * D[perm[3]].y) *
            (D[perm[0]].y * D[perm[2]] - D[perm[2]].y * D[perm[0]]);
        omdyy[44] =
          2 * D[perm[1]].y * D[perm[2]].y * (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]) +
          2 * (D[perm[1]].y * l[perm[2]] + l[perm[1]] * D[perm[2]].y) *
            (D[perm[0]].y * D[perm[3]] - D[perm[3]].y * D[perm[0]]);
      }

      if (whatd & Fop_dzz) {
        omdzz[42] =
          2 * D[perm[2]].z * D[perm[3]].z * (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]) +
          2 * (D[perm[2]].z * l[perm[3]] + l[perm[2]] * D[perm[3]].z) *
            (D[perm[0]].z * D[perm[1]] - D[perm[1]].z * D[perm[0]]);
        omdzz[43] =
          2 * D[perm[1]].z * D[perm[3]].z * (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]) +
          2 * (D[perm[1]].z * l[perm[3]] + l[perm[1]] * D[perm[3]].z) *
            (D[perm[0]].z * D[perm[2]] - D[perm[2]].z * D[perm[0]]);
        omdzz[44] =
          2 * D[perm[1]].z * D[perm[2]].z * (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]) +
          2 * (D[perm[1]].z * l[perm[2]] + l[perm[1]] * D[perm[2]].z) *
            (D[perm[0]].z * D[perm[3]] - D[perm[3]].z * D[perm[0]]);
      }

      if (whatd & Fop_dxy) {
        omdxy[42] = (D[perm[2]].x * D[perm[3]].y + D[perm[2]].y * D[perm[3]].x) *
                      (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]) +
                    (D[perm[2]].x * l[perm[3]] + l[perm[2]] * D[perm[3]].x) *
                      (D[perm[0]].y * D[perm[1]] - D[perm[1]].y * D[perm[0]]) +
                    (D[perm[2]].y * l[perm[3]] + l[perm[2]] * D[perm[3]].y) *
                      (D[perm[0]].x * D[perm[1]] - D[perm[1]].x * D[perm[0]]);

        omdxy[43] = (D[perm[1]].x * D[perm[3]].y + D[perm[1]].y * D[perm[3]].x) *
                      (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]) +
                    (D[perm[1]].x * l[perm[3]] + l[perm[1]] * D[perm[3]].x) *
                      (D[perm[0]].y * D[perm[2]] - D[perm[2]].y * D[perm[0]]) +
                    (D[perm[1]].y * l[perm[3]] + l[perm[1]] * D[perm[3]].y) *
                      (D[perm[0]].x * D[perm[2]] - D[perm[2]].x * D[perm[0]]);

        omdxy[44] = (D[perm[1]].x * D[perm[2]].y + D[perm[1]].y * D[perm[2]].x) *
                      (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]) +
                    (D[perm[1]].x * l[perm[2]] + l[perm[1]] * D[perm[2]].x) *
                      (D[perm[0]].y * D[perm[3]] - D[perm[3]].y * D[perm[0]]) +
                    (D[perm[1]].y * l[perm[2]] + l[perm[1]] * D[perm[2]].y) *
                      (D[perm[0]].x * D[perm[3]] - D[perm[3]].x * D[perm[0]]);
      }

      if (whatd & Fop_dxz) {
        omdxz[42] = (D[perm[2]].x * D[perm[3]].z + D[perm[2]].z * D[perm[3]].x) *
                      (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]) +
                    (D[perm[2]].x * l[perm[3]] + l[perm[2]] * D[perm[3]].x) *
                      (D[perm[0]].z * D[perm[1]] - D[perm[1]].z * D[perm[0]]) +
                    (D[perm[2]].z * l[perm[3]] + l[perm[2]] * D[perm[3]].z) *
                      (D[perm[0]].x * D[perm[1]] - D[perm[1]].x * D[perm[0]]);

        omdxz[43] = (D[perm[1]].x * D[perm[3]].z + D[perm[1]].z * D[perm[3]].x) *
                      (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]) +
                    (D[perm[1]].x * l[perm[3]] + l[perm[1]] * D[perm[3]].x) *
                      (D[perm[0]].z * D[perm[2]] - D[perm[2]].z * D[perm[0]]) +
                    (D[perm[1]].z * l[perm[3]] + l[perm[1]] * D[perm[3]].z) *
                      (D[perm[0]].x * D[perm[2]] - D[perm[2]].x * D[perm[0]]);

        omdxz[44] = (D[perm[1]].x * D[perm[2]].z + D[perm[1]].z * D[perm[2]].x) *
                      (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]) +
                    (D[perm[1]].x * l[perm[2]] + l[perm[1]] * D[perm[2]].x) *
                      (D[perm[0]].z * D[perm[3]] - D[perm[3]].z * D[perm[0]]) +
                    (D[perm[1]].z * l[perm[2]] + l[perm[1]] * D[perm[2]].z) *
                      (D[perm[0]].x * D[perm[3]] - D[perm[3]].x * D[perm[0]]);
      }

      if (whatd & Fop_dyz) {
        omdyz[42] = (D[perm[2]].y * D[perm[3]].z + D[perm[2]].z * D[perm[3]].y) *
                      (l[perm[0]] * D[perm[1]] - l[perm[1]] * D[perm[0]]) +
                    (D[perm[2]].y * l[perm[3]] + l[perm[2]] * D[perm[3]].y) *
                      (D[perm[0]].z * D[perm[1]] - D[perm[1]].z * D[perm[0]]) +
                    (D[perm[2]].z * l[perm[3]] + l[perm[2]] * D[perm[3]].z) *
                      (D[perm[0]].y * D[perm[1]] - D[perm[1]].y * D[perm[0]]);

        omdyz[43] = (D[perm[1]].y * D[perm[3]].z + D[perm[1]].z * D[perm[3]].y) *
                      (l[perm[0]] * D[perm[2]] - l[perm[2]] * D[perm[0]]) +
                    (D[perm[1]].y * l[perm[3]] + l[perm[1]] * D[perm[3]].y) *
                      (D[perm[0]].z * D[perm[2]] - D[perm[2]].z * D[perm[0]]) +
                    (D[perm[1]].z * l[perm[3]] + l[perm[1]] * D[perm[3]].z) *
                      (D[perm[0]].y * D[perm[2]] - D[perm[2]].y * D[perm[0]]);

        omdyz[44] = (D[perm[1]].y * D[perm[2]].z + D[perm[1]].z * D[perm[2]].y) *
                      (l[perm[0]] * D[perm[3]] - l[perm[3]] * D[perm[0]]) +
                    (D[perm[1]].y * l[perm[2]] + l[perm[1]] * D[perm[2]].y) *
                      (D[perm[0]].z * D[perm[3]] - D[perm[3]].z * D[perm[0]]) +
                    (D[perm[1]].z * l[perm[2]] + l[perm[1]] * D[perm[2]].z) *
                      (D[perm[0]].y * D[perm[3]] - D[perm[3]].y * D[perm[0]]);
      }

      R3 phidxx[45];
      if (whatd & Fop_dxx) {
        phidxx[p45[0]] = +9 * omdxx[0] - 18 * omdxx[1] + 3 * omdxx[2] - 27 * omdxx[30] +
                         12 * omdxx[31] + 9 * omdxx[32] + 9 * omdxx[33] - 6 * omdxx[34] -
                         6 * omdxx[35] - 27 * omdxx[36] + 12 * omdxx[37] + 9 * omdxx[38] +
                         9 * omdxx[39] - 6 * omdxx[40] - 6 * omdxx[41] + 18 * omdxx[42] -
                         6 * omdxx[43] - 6 * omdxx[44];
        phidxx[p45[1]] = -18 * omdxx[0] + 84 * omdxx[1] - 18 * omdxx[2] - 6 * omdxx[30] -
                         36 * omdxx[31] + 12 * omdxx[32] - 30 * omdxx[33] + 30 * omdxx[34] -
                         6 * omdxx[36] - 36 * omdxx[37] + 12 * omdxx[38] - 30 * omdxx[39] +
                         30 * omdxx[40] + 24 * omdxx[42];
        phidxx[p45[2]] = +3 * omdxx[0] - 18 * omdxx[1] + 9 * omdxx[2] + 6 * omdxx[30] -
                         18 * omdxx[31] + 3 * omdxx[32] + 6 * omdxx[33] - 9 * omdxx[34] +
                         6 * omdxx[35] + 6 * omdxx[36] - 18 * omdxx[37] + 3 * omdxx[38] +
                         6 * omdxx[39] - 9 * omdxx[40] + 6 * omdxx[41] + 6 * omdxx[42] +
                         6 * omdxx[43] + 6 * omdxx[44];
        phidxx[p45[3]] = +9 * omdxx[3] - 18 * omdxx[4] + 3 * omdxx[5] - 27 * omdxx[24] +
                         12 * omdxx[25] + 9 * omdxx[26] + 9 * omdxx[27] - 6 * omdxx[28] -
                         6 * omdxx[29] + 9 * omdxx[36] - 6 * omdxx[37] - 6 * omdxx[38] -
                         27 * omdxx[39] + 9 * omdxx[40] + 12 * omdxx[41] - 6 * omdxx[42] +
                         18 * omdxx[43] - 6 * omdxx[44];
        phidxx[p45[4]] = -18 * omdxx[3] + 84 * omdxx[4] - 18 * omdxx[5] - 6 * omdxx[24] -
                         36 * omdxx[25] + 12 * omdxx[26] - 30 * omdxx[27] + 30 * omdxx[28] -
                         30 * omdxx[36] + 30 * omdxx[38] - 6 * omdxx[39] + 12 * omdxx[40] -
                         36 * omdxx[41] + 24 * omdxx[43];
        phidxx[p45[5]] = +3 * omdxx[3] - 18 * omdxx[4] + 9 * omdxx[5] + 6 * omdxx[24] -
                         18 * omdxx[25] + 3 * omdxx[26] + 6 * omdxx[27] - 9 * omdxx[28] +
                         6 * omdxx[29] + 6 * omdxx[36] + 6 * omdxx[37] - 9 * omdxx[38] +
                         6 * omdxx[39] + 3 * omdxx[40] - 18 * omdxx[41] + 6 * omdxx[42] +
                         6 * omdxx[43] + 6 * omdxx[44];
        phidxx[p45[6]] = +9 * omdxx[6] - 18 * omdxx[7] + 3 * omdxx[8] + 9 * omdxx[24] -
                         6 * omdxx[25] - 6 * omdxx[26] - 27 * omdxx[27] + 9 * omdxx[28] +
                         12 * omdxx[29] + 9 * omdxx[30] - 6 * omdxx[31] - 6 * omdxx[32] -
                         27 * omdxx[33] + 9 * omdxx[34] + 12 * omdxx[35] - 6 * omdxx[42] -
                         6 * omdxx[43] + 18 * omdxx[44];
        phidxx[p45[7]] = -18 * omdxx[6] + 84 * omdxx[7] - 18 * omdxx[8] - 30 * omdxx[24] +
                         30 * omdxx[26] - 6 * omdxx[27] + 12 * omdxx[28] - 36 * omdxx[29] -
                         30 * omdxx[30] + 30 * omdxx[32] - 6 * omdxx[33] + 12 * omdxx[34] -
                         36 * omdxx[35] + 24 * omdxx[44];
        phidxx[p45[8]] = +3 * omdxx[6] - 18 * omdxx[7] + 9 * omdxx[8] + 6 * omdxx[24] +
                         6 * omdxx[25] - 9 * omdxx[26] + 6 * omdxx[27] + 3 * omdxx[28] -
                         18 * omdxx[29] + 6 * omdxx[30] + 6 * omdxx[31] - 9 * omdxx[32] +
                         6 * omdxx[33] + 3 * omdxx[34] - 18 * omdxx[35] + 6 * omdxx[42] +
                         6 * omdxx[43] + 6 * omdxx[44];
        phidxx[p45[9]] = +9 * omdxx[9] - 18 * omdxx[10] + 3 * omdxx[11] - 27 * omdxx[18] +
                         12 * omdxx[19] + 9 * omdxx[20] + 9 * omdxx[21] - 6 * omdxx[22] -
                         6 * omdxx[23] - 3 * omdxx[36] + 18 * omdxx[37] - 6 * omdxx[38] +
                         9 * omdxx[39] - 27 * omdxx[40] + 12 * omdxx[41] - 6 * omdxx[42] +
                         18 * omdxx[43] - 6 * omdxx[44];
        phidxx[p45[10]] = -18 * omdxx[9] + 84 * omdxx[10] - 18 * omdxx[11] - 6 * omdxx[18] -
                          36 * omdxx[19] + 12 * omdxx[20] - 30 * omdxx[21] + 30 * omdxx[22] -
                          12 * omdxx[36] + 36 * omdxx[37] + 6 * omdxx[38] + 12 * omdxx[39] -
                          6 * omdxx[40] - 36 * omdxx[41] - 24 * omdxx[42] + 24 * omdxx[43];
        phidxx[p45[11]] = +3 * omdxx[9] - 18 * omdxx[10] + 9 * omdxx[11] + 6 * omdxx[18] -
                          18 * omdxx[19] + 3 * omdxx[20] + 6 * omdxx[21] - 9 * omdxx[22] +
                          6 * omdxx[23] - 9 * omdxx[36] - 12 * omdxx[37] + 27 * omdxx[38] +
                          3 * omdxx[39] + 6 * omdxx[40] - 18 * omdxx[41] - 18 * omdxx[42] +
                          6 * omdxx[43] + 6 * omdxx[44];
        phidxx[p45[12]] = +9 * omdxx[12] - 18 * omdxx[13] + 3 * omdxx[14] + 9 * omdxx[18] -
                          6 * omdxx[19] - 6 * omdxx[20] - 27 * omdxx[21] + 9 * omdxx[22] +
                          12 * omdxx[23] - 3 * omdxx[30] + 18 * omdxx[31] - 6 * omdxx[32] +
                          9 * omdxx[33] - 27 * omdxx[34] + 12 * omdxx[35] - 6 * omdxx[42] -
                          6 * omdxx[43] + 18 * omdxx[44];
        phidxx[p45[13]] = -18 * omdxx[12] + 84 * omdxx[13] - 18 * omdxx[14] - 30 * omdxx[18] +
                          30 * omdxx[20] - 6 * omdxx[21] + 12 * omdxx[22] - 36 * omdxx[23] -
                          12 * omdxx[30] + 36 * omdxx[31] + 6 * omdxx[32] + 12 * omdxx[33] -
                          6 * omdxx[34] - 36 * omdxx[35] - 24 * omdxx[42] + 24 * omdxx[44];
        phidxx[p45[14]] = +3 * omdxx[12] - 18 * omdxx[13] + 9 * omdxx[14] + 6 * omdxx[18] +
                          6 * omdxx[19] - 9 * omdxx[20] + 6 * omdxx[21] + 3 * omdxx[22] -
                          18 * omdxx[23] - 9 * omdxx[30] - 12 * omdxx[31] + 27 * omdxx[32] +
                          3 * omdxx[33] + 6 * omdxx[34] - 18 * omdxx[35] - 18 * omdxx[42] +
                          6 * omdxx[43] + 6 * omdxx[44];
        phidxx[p45[15]] = +9 * omdxx[15] - 18 * omdxx[16] + 3 * omdxx[17] - 3 * omdxx[18] +
                          18 * omdxx[19] - 6 * omdxx[20] + 9 * omdxx[21] - 27 * omdxx[22] +
                          12 * omdxx[23] - 3 * omdxx[24] + 18 * omdxx[25] - 6 * omdxx[26] +
                          9 * omdxx[27] - 27 * omdxx[28] + 12 * omdxx[29] - 6 * omdxx[42] -
                          6 * omdxx[43] + 18 * omdxx[44];
        phidxx[p45[16]] = -18 * omdxx[15] + 84 * omdxx[16] - 18 * omdxx[17] - 12 * omdxx[18] +
                          36 * omdxx[19] + 6 * omdxx[20] + 12 * omdxx[21] - 6 * omdxx[22] -
                          36 * omdxx[23] - 12 * omdxx[24] + 36 * omdxx[25] + 6 * omdxx[26] +
                          12 * omdxx[27] - 6 * omdxx[28] - 36 * omdxx[29] - 24 * omdxx[43] +
                          24 * omdxx[44];
        phidxx[p45[17]] = +3 * omdxx[15] - 18 * omdxx[16] + 9 * omdxx[17] - 9 * omdxx[18] -
                          12 * omdxx[19] + 27 * omdxx[20] + 3 * omdxx[21] + 6 * omdxx[22] -
                          18 * omdxx[23] - 9 * omdxx[24] - 12 * omdxx[25] + 27 * omdxx[26] +
                          3 * omdxx[27] + 6 * omdxx[28] - 18 * omdxx[29] + 6 * omdxx[42] -
                          18 * omdxx[43] + 6 * omdxx[44];
        phidxx[p45[18]] = +90 * omdxx[18] - 30 * omdxx[19] - 45 * omdxx[20] - 30 * omdxx[21] +
                          15 * omdxx[22] + 30 * omdxx[23] + 15 * omdxx[42] - 45 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[19]] = -30 * omdxx[18] + 120 * omdxx[19] - 30 * omdxx[20] + 30 * omdxx[21] -
                          60 * omdxx[22] + 30 * omdxx[42] - 30 * omdxx[43] + 30 * omdxx[44];
        phidxx[p45[20]] = -45 * omdxx[18] - 30 * omdxx[19] + 90 * omdxx[20] + 15 * omdxx[21] +
                          15 * omdxx[22] - 60 * omdxx[23] + 15 * omdxx[42] - 45 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[21]] = -30 * omdxx[18] + 30 * omdxx[19] + 15 * omdxx[20] + 90 * omdxx[21] -
                          45 * omdxx[22] - 30 * omdxx[23] + 15 * omdxx[42] + 15 * omdxx[43] -
                          45 * omdxx[44];
        phidxx[p45[22]] = +15 * omdxx[18] - 60 * omdxx[19] + 15 * omdxx[20] - 45 * omdxx[21] +
                          90 * omdxx[22] - 30 * omdxx[23] + 15 * omdxx[42] + 15 * omdxx[43] -
                          45 * omdxx[44];
        phidxx[p45[23]] = +30 * omdxx[18] - 60 * omdxx[20] - 30 * omdxx[21] - 30 * omdxx[22] +
                          120 * omdxx[23] + 30 * omdxx[42] + 30 * omdxx[43] - 30 * omdxx[44];
        phidxx[p45[24]] = +90 * omdxx[24] - 30 * omdxx[25] - 45 * omdxx[26] - 30 * omdxx[27] +
                          15 * omdxx[28] + 30 * omdxx[29] + 15 * omdxx[42] - 45 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[25]] = -30 * omdxx[24] + 120 * omdxx[25] - 30 * omdxx[26] + 30 * omdxx[27] -
                          60 * omdxx[28] - 30 * omdxx[42] - 30 * omdxx[43] + 30 * omdxx[44];
        phidxx[p45[26]] = -45 * omdxx[24] - 30 * omdxx[25] + 90 * omdxx[26] + 15 * omdxx[27] +
                          15 * omdxx[28] - 60 * omdxx[29] + 15 * omdxx[42] - 45 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[27]] = -30 * omdxx[24] + 30 * omdxx[25] + 15 * omdxx[26] + 90 * omdxx[27] -
                          45 * omdxx[28] - 30 * omdxx[29] + 15 * omdxx[42] + 15 * omdxx[43] -
                          45 * omdxx[44];
        phidxx[p45[28]] = +15 * omdxx[24] - 60 * omdxx[25] + 15 * omdxx[26] - 45 * omdxx[27] +
                          90 * omdxx[28] - 30 * omdxx[29] + 15 * omdxx[42] + 15 * omdxx[43] -
                          45 * omdxx[44];
        phidxx[p45[29]] = +30 * omdxx[24] - 60 * omdxx[26] - 30 * omdxx[27] - 30 * omdxx[28] +
                          120 * omdxx[29] - 30 * omdxx[42] + 30 * omdxx[43] - 30 * omdxx[44];
        phidxx[p45[30]] = +90 * omdxx[30] - 30 * omdxx[31] - 45 * omdxx[32] - 30 * omdxx[33] +
                          15 * omdxx[34] + 30 * omdxx[35] - 45 * omdxx[42] + 15 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[31]] = -30 * omdxx[30] + 120 * omdxx[31] - 30 * omdxx[32] + 30 * omdxx[33] -
                          60 * omdxx[34] - 30 * omdxx[42] - 30 * omdxx[43] + 30 * omdxx[44];
        phidxx[p45[32]] = -45 * omdxx[30] - 30 * omdxx[31] + 90 * omdxx[32] + 15 * omdxx[33] +
                          15 * omdxx[34] - 60 * omdxx[35] - 45 * omdxx[42] + 15 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[33]] = -30 * omdxx[30] + 30 * omdxx[31] + 15 * omdxx[32] + 90 * omdxx[33] -
                          45 * omdxx[34] - 30 * omdxx[35] + 15 * omdxx[42] + 15 * omdxx[43] -
                          45 * omdxx[44];
        phidxx[p45[34]] = +15 * omdxx[30] - 60 * omdxx[31] + 15 * omdxx[32] - 45 * omdxx[33] +
                          90 * omdxx[34] - 30 * omdxx[35] + 15 * omdxx[42] + 15 * omdxx[43] -
                          45 * omdxx[44];
        phidxx[p45[35]] = +30 * omdxx[30] - 60 * omdxx[32] - 30 * omdxx[33] - 30 * omdxx[34] +
                          120 * omdxx[35] + 30 * omdxx[42] - 30 * omdxx[43] - 30 * omdxx[44];
        phidxx[p45[36]] = +90 * omdxx[36] - 30 * omdxx[37] - 45 * omdxx[38] - 30 * omdxx[39] +
                          15 * omdxx[40] + 30 * omdxx[41] - 45 * omdxx[42] + 15 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[37]] = -30 * omdxx[36] + 120 * omdxx[37] - 30 * omdxx[38] + 30 * omdxx[39] -
                          60 * omdxx[40] - 30 * omdxx[42] + 30 * omdxx[43] - 30 * omdxx[44];
        phidxx[p45[38]] = -45 * omdxx[36] - 30 * omdxx[37] + 90 * omdxx[38] + 15 * omdxx[39] +
                          15 * omdxx[40] - 60 * omdxx[41] - 45 * omdxx[42] + 15 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[39]] = -30 * omdxx[36] + 30 * omdxx[37] + 15 * omdxx[38] + 90 * omdxx[39] -
                          45 * omdxx[40] - 30 * omdxx[41] + 15 * omdxx[42] - 45 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[40]] = +15 * omdxx[36] - 60 * omdxx[37] + 15 * omdxx[38] - 45 * omdxx[39] +
                          90 * omdxx[40] - 30 * omdxx[41] + 15 * omdxx[42] - 45 * omdxx[43] +
                          15 * omdxx[44];
        phidxx[p45[41]] = +30 * omdxx[36] - 60 * omdxx[38] - 30 * omdxx[39] - 30 * omdxx[40] +
                          120 * omdxx[41] + 30 * omdxx[42] - 30 * omdxx[43] - 30 * omdxx[44];
        phidxx[p45[42]] = +90 * omdxx[42] - 30 * omdxx[43] - 30 * omdxx[44];
        phidxx[p45[43]] = -30 * omdxx[42] + 90 * omdxx[43] - 30 * omdxx[44];
        phidxx[p45[44]] = -30 * omdxx[42] - 30 * omdxx[43] + 90 * omdxx[44];

        for (int k = 0; k < 45; ++k) {
          val(k, 0, op_dxx) = phidxx[k].x;
          val(k, 1, op_dxx) = phidxx[k].y;
          val(k, 2, op_dxx) = phidxx[k].z;
        }
      }

      R3 phidyy[45];
      if (whatd & Fop_dyy) {
        phidyy[p45[0]] = +9 * omdyy[0] - 18 * omdyy[1] + 3 * omdyy[2] - 27 * omdyy[30] +
                         12 * omdyy[31] + 9 * omdyy[32] + 9 * omdyy[33] - 6 * omdyy[34] -
                         6 * omdyy[35] - 27 * omdyy[36] + 12 * omdyy[37] + 9 * omdyy[38] +
                         9 * omdyy[39] - 6 * omdyy[40] - 6 * omdyy[41] + 18 * omdyy[42] -
                         6 * omdyy[43] - 6 * omdyy[44];
        phidyy[p45[1]] = -18 * omdyy[0] + 84 * omdyy[1] - 18 * omdyy[2] - 6 * omdyy[30] -
                         36 * omdyy[31] + 12 * omdyy[32] - 30 * omdyy[33] + 30 * omdyy[34] -
                         6 * omdyy[36] - 36 * omdyy[37] + 12 * omdyy[38] - 30 * omdyy[39] +
                         30 * omdyy[40] + 24 * omdyy[42];
        phidyy[p45[2]] = +3 * omdyy[0] - 18 * omdyy[1] + 9 * omdyy[2] + 6 * omdyy[30] -
                         18 * omdyy[31] + 3 * omdyy[32] + 6 * omdyy[33] - 9 * omdyy[34] +
                         6 * omdyy[35] + 6 * omdyy[36] - 18 * omdyy[37] + 3 * omdyy[38] +
                         6 * omdyy[39] - 9 * omdyy[40] + 6 * omdyy[41] + 6 * omdyy[42] +
                         6 * omdyy[43] + 6 * omdyy[44];
        phidyy[p45[3]] = +9 * omdyy[3] - 18 * omdyy[4] + 3 * omdyy[5] - 27 * omdyy[24] +
                         12 * omdyy[25] + 9 * omdyy[26] + 9 * omdyy[27] - 6 * omdyy[28] -
                         6 * omdyy[29] + 9 * omdyy[36] - 6 * omdyy[37] - 6 * omdyy[38] -
                         27 * omdyy[39] + 9 * omdyy[40] + 12 * omdyy[41] - 6 * omdyy[42] +
                         18 * omdyy[43] - 6 * omdyy[44];
        phidyy[p45[4]] = -18 * omdyy[3] + 84 * omdyy[4] - 18 * omdyy[5] - 6 * omdyy[24] -
                         36 * omdyy[25] + 12 * omdyy[26] - 30 * omdyy[27] + 30 * omdyy[28] -
                         30 * omdyy[36] + 30 * omdyy[38] - 6 * omdyy[39] + 12 * omdyy[40] -
                         36 * omdyy[41] + 24 * omdyy[43];
        phidyy[p45[5]] = +3 * omdyy[3] - 18 * omdyy[4] + 9 * omdyy[5] + 6 * omdyy[24] -
                         18 * omdyy[25] + 3 * omdyy[26] + 6 * omdyy[27] - 9 * omdyy[28] +
                         6 * omdyy[29] + 6 * omdyy[36] + 6 * omdyy[37] - 9 * omdyy[38] +
                         6 * omdyy[39] + 3 * omdyy[40] - 18 * omdyy[41] + 6 * omdyy[42] +
                         6 * omdyy[43] + 6 * omdyy[44];
        phidyy[p45[6]] = +9 * omdyy[6] - 18 * omdyy[7] + 3 * omdyy[8] + 9 * omdyy[24] -
                         6 * omdyy[25] - 6 * omdyy[26] - 27 * omdyy[27] + 9 * omdyy[28] +
                         12 * omdyy[29] + 9 * omdyy[30] - 6 * omdyy[31] - 6 * omdyy[32] -
                         27 * omdyy[33] + 9 * omdyy[34] + 12 * omdyy[35] - 6 * omdyy[42] -
                         6 * omdyy[43] + 18 * omdyy[44];
        phidyy[p45[7]] = -18 * omdyy[6] + 84 * omdyy[7] - 18 * omdyy[8] - 30 * omdyy[24] +
                         30 * omdyy[26] - 6 * omdyy[27] + 12 * omdyy[28] - 36 * omdyy[29] -
                         30 * omdyy[30] + 30 * omdyy[32] - 6 * omdyy[33] + 12 * omdyy[34] -
                         36 * omdyy[35] + 24 * omdyy[44];
        phidyy[p45[8]] = +3 * omdyy[6] - 18 * omdyy[7] + 9 * omdyy[8] + 6 * omdyy[24] +
                         6 * omdyy[25] - 9 * omdyy[26] + 6 * omdyy[27] + 3 * omdyy[28] -
                         18 * omdyy[29] + 6 * omdyy[30] + 6 * omdyy[31] - 9 * omdyy[32] +
                         6 * omdyy[33] + 3 * omdyy[34] - 18 * omdyy[35] + 6 * omdyy[42] +
                         6 * omdyy[43] + 6 * omdyy[44];
        phidyy[p45[9]] = +9 * omdyy[9] - 18 * omdyy[10] + 3 * omdyy[11] - 27 * omdyy[18] +
                         12 * omdyy[19] + 9 * omdyy[20] + 9 * omdyy[21] - 6 * omdyy[22] -
                         6 * omdyy[23] - 3 * omdyy[36] + 18 * omdyy[37] - 6 * omdyy[38] +
                         9 * omdyy[39] - 27 * omdyy[40] + 12 * omdyy[41] - 6 * omdyy[42] +
                         18 * omdyy[43] - 6 * omdyy[44];
        phidyy[p45[10]] = -18 * omdyy[9] + 84 * omdyy[10] - 18 * omdyy[11] - 6 * omdyy[18] -
                          36 * omdyy[19] + 12 * omdyy[20] - 30 * omdyy[21] + 30 * omdyy[22] -
                          12 * omdyy[36] + 36 * omdyy[37] + 6 * omdyy[38] + 12 * omdyy[39] -
                          6 * omdyy[40] - 36 * omdyy[41] - 24 * omdyy[42] + 24 * omdyy[43];
        phidyy[p45[11]] = +3 * omdyy[9] - 18 * omdyy[10] + 9 * omdyy[11] + 6 * omdyy[18] -
                          18 * omdyy[19] + 3 * omdyy[20] + 6 * omdyy[21] - 9 * omdyy[22] +
                          6 * omdyy[23] - 9 * omdyy[36] - 12 * omdyy[37] + 27 * omdyy[38] +
                          3 * omdyy[39] + 6 * omdyy[40] - 18 * omdyy[41] - 18 * omdyy[42] +
                          6 * omdyy[43] + 6 * omdyy[44];
        phidyy[p45[12]] = +9 * omdyy[12] - 18 * omdyy[13] + 3 * omdyy[14] + 9 * omdyy[18] -
                          6 * omdyy[19] - 6 * omdyy[20] - 27 * omdyy[21] + 9 * omdyy[22] +
                          12 * omdyy[23] - 3 * omdyy[30] + 18 * omdyy[31] - 6 * omdyy[32] +
                          9 * omdyy[33] - 27 * omdyy[34] + 12 * omdyy[35] - 6 * omdyy[42] -
                          6 * omdyy[43] + 18 * omdyy[44];
        phidyy[p45[13]] = -18 * omdyy[12] + 84 * omdyy[13] - 18 * omdyy[14] - 30 * omdyy[18] +
                          30 * omdyy[20] - 6 * omdyy[21] + 12 * omdyy[22] - 36 * omdyy[23] -
                          12 * omdyy[30] + 36 * omdyy[31] + 6 * omdyy[32] + 12 * omdyy[33] -
                          6 * omdyy[34] - 36 * omdyy[35] - 24 * omdyy[42] + 24 * omdyy[44];
        phidyy[p45[14]] = +3 * omdyy[12] - 18 * omdyy[13] + 9 * omdyy[14] + 6 * omdyy[18] +
                          6 * omdyy[19] - 9 * omdyy[20] + 6 * omdyy[21] + 3 * omdyy[22] -
                          18 * omdyy[23] - 9 * omdyy[30] - 12 * omdyy[31] + 27 * omdyy[32] +
                          3 * omdyy[33] + 6 * omdyy[34] - 18 * omdyy[35] - 18 * omdyy[42] +
                          6 * omdyy[43] + 6 * omdyy[44];
        phidyy[p45[15]] = +9 * omdyy[15] - 18 * omdyy[16] + 3 * omdyy[17] - 3 * omdyy[18] +
                          18 * omdyy[19] - 6 * omdyy[20] + 9 * omdyy[21] - 27 * omdyy[22] +
                          12 * omdyy[23] - 3 * omdyy[24] + 18 * omdyy[25] - 6 * omdyy[26] +
                          9 * omdyy[27] - 27 * omdyy[28] + 12 * omdyy[29] - 6 * omdyy[42] -
                          6 * omdyy[43] + 18 * omdyy[44];
        phidyy[p45[16]] = -18 * omdyy[15] + 84 * omdyy[16] - 18 * omdyy[17] - 12 * omdyy[18] +
                          36 * omdyy[19] + 6 * omdyy[20] + 12 * omdyy[21] - 6 * omdyy[22] -
                          36 * omdyy[23] - 12 * omdyy[24] + 36 * omdyy[25] + 6 * omdyy[26] +
                          12 * omdyy[27] - 6 * omdyy[28] - 36 * omdyy[29] - 24 * omdyy[43] +
                          24 * omdyy[44];
        phidyy[p45[17]] = +3 * omdyy[15] - 18 * omdyy[16] + 9 * omdyy[17] - 9 * omdyy[18] -
                          12 * omdyy[19] + 27 * omdyy[20] + 3 * omdyy[21] + 6 * omdyy[22] -
                          18 * omdyy[23] - 9 * omdyy[24] - 12 * omdyy[25] + 27 * omdyy[26] +
                          3 * omdyy[27] + 6 * omdyy[28] - 18 * omdyy[29] + 6 * omdyy[42] -
                          18 * omdyy[43] + 6 * omdyy[44];
        phidyy[p45[18]] = +90 * omdyy[18] - 30 * omdyy[19] - 45 * omdyy[20] - 30 * omdyy[21] +
                          15 * omdyy[22] + 30 * omdyy[23] + 15 * omdyy[42] - 45 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[19]] = -30 * omdyy[18] + 120 * omdyy[19] - 30 * omdyy[20] + 30 * omdyy[21] -
                          60 * omdyy[22] + 30 * omdyy[42] - 30 * omdyy[43] + 30 * omdyy[44];
        phidyy[p45[20]] = -45 * omdyy[18] - 30 * omdyy[19] + 90 * omdyy[20] + 15 * omdyy[21] +
                          15 * omdyy[22] - 60 * omdyy[23] + 15 * omdyy[42] - 45 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[21]] = -30 * omdyy[18] + 30 * omdyy[19] + 15 * omdyy[20] + 90 * omdyy[21] -
                          45 * omdyy[22] - 30 * omdyy[23] + 15 * omdyy[42] + 15 * omdyy[43] -
                          45 * omdyy[44];
        phidyy[p45[22]] = +15 * omdyy[18] - 60 * omdyy[19] + 15 * omdyy[20] - 45 * omdyy[21] +
                          90 * omdyy[22] - 30 * omdyy[23] + 15 * omdyy[42] + 15 * omdyy[43] -
                          45 * omdyy[44];
        phidyy[p45[23]] = +30 * omdyy[18] - 60 * omdyy[20] - 30 * omdyy[21] - 30 * omdyy[22] +
                          120 * omdyy[23] + 30 * omdyy[42] + 30 * omdyy[43] - 30 * omdyy[44];
        phidyy[p45[24]] = +90 * omdyy[24] - 30 * omdyy[25] - 45 * omdyy[26] - 30 * omdyy[27] +
                          15 * omdyy[28] + 30 * omdyy[29] + 15 * omdyy[42] - 45 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[25]] = -30 * omdyy[24] + 120 * omdyy[25] - 30 * omdyy[26] + 30 * omdyy[27] -
                          60 * omdyy[28] - 30 * omdyy[42] - 30 * omdyy[43] + 30 * omdyy[44];
        phidyy[p45[26]] = -45 * omdyy[24] - 30 * omdyy[25] + 90 * omdyy[26] + 15 * omdyy[27] +
                          15 * omdyy[28] - 60 * omdyy[29] + 15 * omdyy[42] - 45 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[27]] = -30 * omdyy[24] + 30 * omdyy[25] + 15 * omdyy[26] + 90 * omdyy[27] -
                          45 * omdyy[28] - 30 * omdyy[29] + 15 * omdyy[42] + 15 * omdyy[43] -
                          45 * omdyy[44];
        phidyy[p45[28]] = +15 * omdyy[24] - 60 * omdyy[25] + 15 * omdyy[26] - 45 * omdyy[27] +
                          90 * omdyy[28] - 30 * omdyy[29] + 15 * omdyy[42] + 15 * omdyy[43] -
                          45 * omdyy[44];
        phidyy[p45[29]] = +30 * omdyy[24] - 60 * omdyy[26] - 30 * omdyy[27] - 30 * omdyy[28] +
                          120 * omdyy[29] - 30 * omdyy[42] + 30 * omdyy[43] - 30 * omdyy[44];
        phidyy[p45[30]] = +90 * omdyy[30] - 30 * omdyy[31] - 45 * omdyy[32] - 30 * omdyy[33] +
                          15 * omdyy[34] + 30 * omdyy[35] - 45 * omdyy[42] + 15 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[31]] = -30 * omdyy[30] + 120 * omdyy[31] - 30 * omdyy[32] + 30 * omdyy[33] -
                          60 * omdyy[34] - 30 * omdyy[42] - 30 * omdyy[43] + 30 * omdyy[44];
        phidyy[p45[32]] = -45 * omdyy[30] - 30 * omdyy[31] + 90 * omdyy[32] + 15 * omdyy[33] +
                          15 * omdyy[34] - 60 * omdyy[35] - 45 * omdyy[42] + 15 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[33]] = -30 * omdyy[30] + 30 * omdyy[31] + 15 * omdyy[32] + 90 * omdyy[33] -
                          45 * omdyy[34] - 30 * omdyy[35] + 15 * omdyy[42] + 15 * omdyy[43] -
                          45 * omdyy[44];
        phidyy[p45[34]] = +15 * omdyy[30] - 60 * omdyy[31] + 15 * omdyy[32] - 45 * omdyy[33] +
                          90 * omdyy[34] - 30 * omdyy[35] + 15 * omdyy[42] + 15 * omdyy[43] -
                          45 * omdyy[44];
        phidyy[p45[35]] = +30 * omdyy[30] - 60 * omdyy[32] - 30 * omdyy[33] - 30 * omdyy[34] +
                          120 * omdyy[35] + 30 * omdyy[42] - 30 * omdyy[43] - 30 * omdyy[44];
        phidyy[p45[36]] = +90 * omdyy[36] - 30 * omdyy[37] - 45 * omdyy[38] - 30 * omdyy[39] +
                          15 * omdyy[40] + 30 * omdyy[41] - 45 * omdyy[42] + 15 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[37]] = -30 * omdyy[36] + 120 * omdyy[37] - 30 * omdyy[38] + 30 * omdyy[39] -
                          60 * omdyy[40] - 30 * omdyy[42] + 30 * omdyy[43] - 30 * omdyy[44];
        phidyy[p45[38]] = -45 * omdyy[36] - 30 * omdyy[37] + 90 * omdyy[38] + 15 * omdyy[39] +
                          15 * omdyy[40] - 60 * omdyy[41] - 45 * omdyy[42] + 15 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[39]] = -30 * omdyy[36] + 30 * omdyy[37] + 15 * omdyy[38] + 90 * omdyy[39] -
                          45 * omdyy[40] - 30 * omdyy[41] + 15 * omdyy[42] - 45 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[40]] = +15 * omdyy[36] - 60 * omdyy[37] + 15 * omdyy[38] - 45 * omdyy[39] +
                          90 * omdyy[40] - 30 * omdyy[41] + 15 * omdyy[42] - 45 * omdyy[43] +
                          15 * omdyy[44];
        phidyy[p45[41]] = +30 * omdyy[36] - 60 * omdyy[38] - 30 * omdyy[39] - 30 * omdyy[40] +
                          120 * omdyy[41] + 30 * omdyy[42] - 30 * omdyy[43] - 30 * omdyy[44];
        phidyy[p45[42]] = +90 * omdyy[42] - 30 * omdyy[43] - 30 * omdyy[44];
        phidyy[p45[43]] = -30 * omdyy[42] + 90 * omdyy[43] - 30 * omdyy[44];
        phidyy[p45[44]] = -30 * omdyy[42] - 30 * omdyy[43] + 90 * omdyy[44];

        for (int k = 0; k < 45; ++k) {
          val(k, 0, op_dyy) = phidyy[k].x;
          val(k, 1, op_dyy) = phidyy[k].y;
          val(k, 2, op_dyy) = phidyy[k].z;
        }
      }

      R3 phidzz[45];
      if (whatd & Fop_dzz) {
        phidzz[p45[0]] = +9 * omdzz[0] - 18 * omdzz[1] + 3 * omdzz[2] - 27 * omdzz[30] +
                         12 * omdzz[31] + 9 * omdzz[32] + 9 * omdzz[33] - 6 * omdzz[34] -
                         6 * omdzz[35] - 27 * omdzz[36] + 12 * omdzz[37] + 9 * omdzz[38] +
                         9 * omdzz[39] - 6 * omdzz[40] - 6 * omdzz[41] + 18 * omdzz[42] -
                         6 * omdzz[43] - 6 * omdzz[44];
        phidzz[p45[1]] = -18 * omdzz[0] + 84 * omdzz[1] - 18 * omdzz[2] - 6 * omdzz[30] -
                         36 * omdzz[31] + 12 * omdzz[32] - 30 * omdzz[33] + 30 * omdzz[34] -
                         6 * omdzz[36] - 36 * omdzz[37] + 12 * omdzz[38] - 30 * omdzz[39] +
                         30 * omdzz[40] + 24 * omdzz[42];
        phidzz[p45[2]] = +3 * omdzz[0] - 18 * omdzz[1] + 9 * omdzz[2] + 6 * omdzz[30] -
                         18 * omdzz[31] + 3 * omdzz[32] + 6 * omdzz[33] - 9 * omdzz[34] +
                         6 * omdzz[35] + 6 * omdzz[36] - 18 * omdzz[37] + 3 * omdzz[38] +
                         6 * omdzz[39] - 9 * omdzz[40] + 6 * omdzz[41] + 6 * omdzz[42] +
                         6 * omdzz[43] + 6 * omdzz[44];
        phidzz[p45[3]] = +9 * omdzz[3] - 18 * omdzz[4] + 3 * omdzz[5] - 27 * omdzz[24] +
                         12 * omdzz[25] + 9 * omdzz[26] + 9 * omdzz[27] - 6 * omdzz[28] -
                         6 * omdzz[29] + 9 * omdzz[36] - 6 * omdzz[37] - 6 * omdzz[38] -
                         27 * omdzz[39] + 9 * omdzz[40] + 12 * omdzz[41] - 6 * omdzz[42] +
                         18 * omdzz[43] - 6 * omdzz[44];
        phidzz[p45[4]] = -18 * omdzz[3] + 84 * omdzz[4] - 18 * omdzz[5] - 6 * omdzz[24] -
                         36 * omdzz[25] + 12 * omdzz[26] - 30 * omdzz[27] + 30 * omdzz[28] -
                         30 * omdzz[36] + 30 * omdzz[38] - 6 * omdzz[39] + 12 * omdzz[40] -
                         36 * omdzz[41] + 24 * omdzz[43];
        phidzz[p45[5]] = +3 * omdzz[3] - 18 * omdzz[4] + 9 * omdzz[5] + 6 * omdzz[24] -
                         18 * omdzz[25] + 3 * omdzz[26] + 6 * omdzz[27] - 9 * omdzz[28] +
                         6 * omdzz[29] + 6 * omdzz[36] + 6 * omdzz[37] - 9 * omdzz[38] +
                         6 * omdzz[39] + 3 * omdzz[40] - 18 * omdzz[41] + 6 * omdzz[42] +
                         6 * omdzz[43] + 6 * omdzz[44];
        phidzz[p45[6]] = +9 * omdzz[6] - 18 * omdzz[7] + 3 * omdzz[8] + 9 * omdzz[24] -
                         6 * omdzz[25] - 6 * omdzz[26] - 27 * omdzz[27] + 9 * omdzz[28] +
                         12 * omdzz[29] + 9 * omdzz[30] - 6 * omdzz[31] - 6 * omdzz[32] -
                         27 * omdzz[33] + 9 * omdzz[34] + 12 * omdzz[35] - 6 * omdzz[42] -
                         6 * omdzz[43] + 18 * omdzz[44];
        phidzz[p45[7]] = -18 * omdzz[6] + 84 * omdzz[7] - 18 * omdzz[8] - 30 * omdzz[24] +
                         30 * omdzz[26] - 6 * omdzz[27] + 12 * omdzz[28] - 36 * omdzz[29] -
                         30 * omdzz[30] + 30 * omdzz[32] - 6 * omdzz[33] + 12 * omdzz[34] -
                         36 * omdzz[35] + 24 * omdzz[44];
        phidzz[p45[8]] = +3 * omdzz[6] - 18 * omdzz[7] + 9 * omdzz[8] + 6 * omdzz[24] +
                         6 * omdzz[25] - 9 * omdzz[26] + 6 * omdzz[27] + 3 * omdzz[28] -
                         18 * omdzz[29] + 6 * omdzz[30] + 6 * omdzz[31] - 9 * omdzz[32] +
                         6 * omdzz[33] + 3 * omdzz[34] - 18 * omdzz[35] + 6 * omdzz[42] +
                         6 * omdzz[43] + 6 * omdzz[44];
        phidzz[p45[9]] = +9 * omdzz[9] - 18 * omdzz[10] + 3 * omdzz[11] - 27 * omdzz[18] +
                         12 * omdzz[19] + 9 * omdzz[20] + 9 * omdzz[21] - 6 * omdzz[22] -
                         6 * omdzz[23] - 3 * omdzz[36] + 18 * omdzz[37] - 6 * omdzz[38] +
                         9 * omdzz[39] - 27 * omdzz[40] + 12 * omdzz[41] - 6 * omdzz[42] +
                         18 * omdzz[43] - 6 * omdzz[44];
        phidzz[p45[10]] = -18 * omdzz[9] + 84 * omdzz[10] - 18 * omdzz[11] - 6 * omdzz[18] -
                          36 * omdzz[19] + 12 * omdzz[20] - 30 * omdzz[21] + 30 * omdzz[22] -
                          12 * omdzz[36] + 36 * omdzz[37] + 6 * omdzz[38] + 12 * omdzz[39] -
                          6 * omdzz[40] - 36 * omdzz[41] - 24 * omdzz[42] + 24 * omdzz[43];
        phidzz[p45[11]] = +3 * omdzz[9] - 18 * omdzz[10] + 9 * omdzz[11] + 6 * omdzz[18] -
                          18 * omdzz[19] + 3 * omdzz[20] + 6 * omdzz[21] - 9 * omdzz[22] +
                          6 * omdzz[23] - 9 * omdzz[36] - 12 * omdzz[37] + 27 * omdzz[38] +
                          3 * omdzz[39] + 6 * omdzz[40] - 18 * omdzz[41] - 18 * omdzz[42] +
                          6 * omdzz[43] + 6 * omdzz[44];
        phidzz[p45[12]] = +9 * omdzz[12] - 18 * omdzz[13] + 3 * omdzz[14] + 9 * omdzz[18] -
                          6 * omdzz[19] - 6 * omdzz[20] - 27 * omdzz[21] + 9 * omdzz[22] +
                          12 * omdzz[23] - 3 * omdzz[30] + 18 * omdzz[31] - 6 * omdzz[32] +
                          9 * omdzz[33] - 27 * omdzz[34] + 12 * omdzz[35] - 6 * omdzz[42] -
                          6 * omdzz[43] + 18 * omdzz[44];
        phidzz[p45[13]] = -18 * omdzz[12] + 84 * omdzz[13] - 18 * omdzz[14] - 30 * omdzz[18] +
                          30 * omdzz[20] - 6 * omdzz[21] + 12 * omdzz[22] - 36 * omdzz[23] -
                          12 * omdzz[30] + 36 * omdzz[31] + 6 * omdzz[32] + 12 * omdzz[33] -
                          6 * omdzz[34] - 36 * omdzz[35] - 24 * omdzz[42] + 24 * omdzz[44];
        phidzz[p45[14]] = +3 * omdzz[12] - 18 * omdzz[13] + 9 * omdzz[14] + 6 * omdzz[18] +
                          6 * omdzz[19] - 9 * omdzz[20] + 6 * omdzz[21] + 3 * omdzz[22] -
                          18 * omdzz[23] - 9 * omdzz[30] - 12 * omdzz[31] + 27 * omdzz[32] +
                          3 * omdzz[33] + 6 * omdzz[34] - 18 * omdzz[35] - 18 * omdzz[42] +
                          6 * omdzz[43] + 6 * omdzz[44];
        phidzz[p45[15]] = +9 * omdzz[15] - 18 * omdzz[16] + 3 * omdzz[17] - 3 * omdzz[18] +
                          18 * omdzz[19] - 6 * omdzz[20] + 9 * omdzz[21] - 27 * omdzz[22] +
                          12 * omdzz[23] - 3 * omdzz[24] + 18 * omdzz[25] - 6 * omdzz[26] +
                          9 * omdzz[27] - 27 * omdzz[28] + 12 * omdzz[29] - 6 * omdzz[42] -
                          6 * omdzz[43] + 18 * omdzz[44];
        phidzz[p45[16]] = -18 * omdzz[15] + 84 * omdzz[16] - 18 * omdzz[17] - 12 * omdzz[18] +
                          36 * omdzz[19] + 6 * omdzz[20] + 12 * omdzz[21] - 6 * omdzz[22] -
                          36 * omdzz[23] - 12 * omdzz[24] + 36 * omdzz[25] + 6 * omdzz[26] +
                          12 * omdzz[27] - 6 * omdzz[28] - 36 * omdzz[29] - 24 * omdzz[43] +
                          24 * omdzz[44];
        phidzz[p45[17]] = +3 * omdzz[15] - 18 * omdzz[16] + 9 * omdzz[17] - 9 * omdzz[18] -
                          12 * omdzz[19] + 27 * omdzz[20] + 3 * omdzz[21] + 6 * omdzz[22] -
                          18 * omdzz[23] - 9 * omdzz[24] - 12 * omdzz[25] + 27 * omdzz[26] +
                          3 * omdzz[27] + 6 * omdzz[28] - 18 * omdzz[29] + 6 * omdzz[42] -
                          18 * omdzz[43] + 6 * omdzz[44];
        phidzz[p45[18]] = +90 * omdzz[18] - 30 * omdzz[19] - 45 * omdzz[20] - 30 * omdzz[21] +
                          15 * omdzz[22] + 30 * omdzz[23] + 15 * omdzz[42] - 45 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[19]] = -30 * omdzz[18] + 120 * omdzz[19] - 30 * omdzz[20] + 30 * omdzz[21] -
                          60 * omdzz[22] + 30 * omdzz[42] - 30 * omdzz[43] + 30 * omdzz[44];
        phidzz[p45[20]] = -45 * omdzz[18] - 30 * omdzz[19] + 90 * omdzz[20] + 15 * omdzz[21] +
                          15 * omdzz[22] - 60 * omdzz[23] + 15 * omdzz[42] - 45 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[21]] = -30 * omdzz[18] + 30 * omdzz[19] + 15 * omdzz[20] + 90 * omdzz[21] -
                          45 * omdzz[22] - 30 * omdzz[23] + 15 * omdzz[42] + 15 * omdzz[43] -
                          45 * omdzz[44];
        phidzz[p45[22]] = +15 * omdzz[18] - 60 * omdzz[19] + 15 * omdzz[20] - 45 * omdzz[21] +
                          90 * omdzz[22] - 30 * omdzz[23] + 15 * omdzz[42] + 15 * omdzz[43] -
                          45 * omdzz[44];
        phidzz[p45[23]] = +30 * omdzz[18] - 60 * omdzz[20] - 30 * omdzz[21] - 30 * omdzz[22] +
                          120 * omdzz[23] + 30 * omdzz[42] + 30 * omdzz[43] - 30 * omdzz[44];
        phidzz[p45[24]] = +90 * omdzz[24] - 30 * omdzz[25] - 45 * omdzz[26] - 30 * omdzz[27] +
                          15 * omdzz[28] + 30 * omdzz[29] + 15 * omdzz[42] - 45 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[25]] = -30 * omdzz[24] + 120 * omdzz[25] - 30 * omdzz[26] + 30 * omdzz[27] -
                          60 * omdzz[28] - 30 * omdzz[42] - 30 * omdzz[43] + 30 * omdzz[44];
        phidzz[p45[26]] = -45 * omdzz[24] - 30 * omdzz[25] + 90 * omdzz[26] + 15 * omdzz[27] +
                          15 * omdzz[28] - 60 * omdzz[29] + 15 * omdzz[42] - 45 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[27]] = -30 * omdzz[24] + 30 * omdzz[25] + 15 * omdzz[26] + 90 * omdzz[27] -
                          45 * omdzz[28] - 30 * omdzz[29] + 15 * omdzz[42] + 15 * omdzz[43] -
                          45 * omdzz[44];
        phidzz[p45[28]] = +15 * omdzz[24] - 60 * omdzz[25] + 15 * omdzz[26] - 45 * omdzz[27] +
                          90 * omdzz[28] - 30 * omdzz[29] + 15 * omdzz[42] + 15 * omdzz[43] -
                          45 * omdzz[44];
        phidzz[p45[29]] = +30 * omdzz[24] - 60 * omdzz[26] - 30 * omdzz[27] - 30 * omdzz[28] +
                          120 * omdzz[29] - 30 * omdzz[42] + 30 * omdzz[43] - 30 * omdzz[44];
        phidzz[p45[30]] = +90 * omdzz[30] - 30 * omdzz[31] - 45 * omdzz[32] - 30 * omdzz[33] +
                          15 * omdzz[34] + 30 * omdzz[35] - 45 * omdzz[42] + 15 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[31]] = -30 * omdzz[30] + 120 * omdzz[31] - 30 * omdzz[32] + 30 * omdzz[33] -
                          60 * omdzz[34] - 30 * omdzz[42] - 30 * omdzz[43] + 30 * omdzz[44];
        phidzz[p45[32]] = -45 * omdzz[30] - 30 * omdzz[31] + 90 * omdzz[32] + 15 * omdzz[33] +
                          15 * omdzz[34] - 60 * omdzz[35] - 45 * omdzz[42] + 15 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[33]] = -30 * omdzz[30] + 30 * omdzz[31] + 15 * omdzz[32] + 90 * omdzz[33] -
                          45 * omdzz[34] - 30 * omdzz[35] + 15 * omdzz[42] + 15 * omdzz[43] -
                          45 * omdzz[44];
        phidzz[p45[34]] = +15 * omdzz[30] - 60 * omdzz[31] + 15 * omdzz[32] - 45 * omdzz[33] +
                          90 * omdzz[34] - 30 * omdzz[35] + 15 * omdzz[42] + 15 * omdzz[43] -
                          45 * omdzz[44];
        phidzz[p45[35]] = +30 * omdzz[30] - 60 * omdzz[32] - 30 * omdzz[33] - 30 * omdzz[34] +
                          120 * omdzz[35] + 30 * omdzz[42] - 30 * omdzz[43] - 30 * omdzz[44];
        phidzz[p45[36]] = +90 * omdzz[36] - 30 * omdzz[37] - 45 * omdzz[38] - 30 * omdzz[39] +
                          15 * omdzz[40] + 30 * omdzz[41] - 45 * omdzz[42] + 15 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[37]] = -30 * omdzz[36] + 120 * omdzz[37] - 30 * omdzz[38] + 30 * omdzz[39] -
                          60 * omdzz[40] - 30 * omdzz[42] + 30 * omdzz[43] - 30 * omdzz[44];
        phidzz[p45[38]] = -45 * omdzz[36] - 30 * omdzz[37] + 90 * omdzz[38] + 15 * omdzz[39] +
                          15 * omdzz[40] - 60 * omdzz[41] - 45 * omdzz[42] + 15 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[39]] = -30 * omdzz[36] + 30 * omdzz[37] + 15 * omdzz[38] + 90 * omdzz[39] -
                          45 * omdzz[40] - 30 * omdzz[41] + 15 * omdzz[42] - 45 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[40]] = +15 * omdzz[36] - 60 * omdzz[37] + 15 * omdzz[38] - 45 * omdzz[39] +
                          90 * omdzz[40] - 30 * omdzz[41] + 15 * omdzz[42] - 45 * omdzz[43] +
                          15 * omdzz[44];
        phidzz[p45[41]] = +30 * omdzz[36] - 60 * omdzz[38] - 30 * omdzz[39] - 30 * omdzz[40] +
                          120 * omdzz[41] + 30 * omdzz[42] - 30 * omdzz[43] - 30 * omdzz[44];
        phidzz[p45[42]] = +90 * omdzz[42] - 30 * omdzz[43] - 30 * omdzz[44];
        phidzz[p45[43]] = -30 * omdzz[42] + 90 * omdzz[43] - 30 * omdzz[44];
        phidzz[p45[44]] = -30 * omdzz[42] - 30 * omdzz[43] + 90 * omdzz[44];

        for (int k = 0; k < 45; ++k) {
          val(k, 0, op_dzz) = phidzz[k].x;
          val(k, 1, op_dzz) = phidzz[k].y;
          val(k, 2, op_dzz) = phidzz[k].z;
        }
      }

      R3 phidxy[45];
      if (whatd & Fop_dxy) {
        phidxy[p45[0]] = +9 * omdxy[0] - 18 * omdxy[1] + 3 * omdxy[2] - 27 * omdxy[30] +
                         12 * omdxy[31] + 9 * omdxy[32] + 9 * omdxy[33] - 6 * omdxy[34] -
                         6 * omdxy[35] - 27 * omdxy[36] + 12 * omdxy[37] + 9 * omdxy[38] +
                         9 * omdxy[39] - 6 * omdxy[40] - 6 * omdxy[41] + 18 * omdxy[42] -
                         6 * omdxy[43] - 6 * omdxy[44];
        phidxy[p45[1]] = -18 * omdxy[0] + 84 * omdxy[1] - 18 * omdxy[2] - 6 * omdxy[30] -
                         36 * omdxy[31] + 12 * omdxy[32] - 30 * omdxy[33] + 30 * omdxy[34] -
                         6 * omdxy[36] - 36 * omdxy[37] + 12 * omdxy[38] - 30 * omdxy[39] +
                         30 * omdxy[40] + 24 * omdxy[42];
        phidxy[p45[2]] = +3 * omdxy[0] - 18 * omdxy[1] + 9 * omdxy[2] + 6 * omdxy[30] -
                         18 * omdxy[31] + 3 * omdxy[32] + 6 * omdxy[33] - 9 * omdxy[34] +
                         6 * omdxy[35] + 6 * omdxy[36] - 18 * omdxy[37] + 3 * omdxy[38] +
                         6 * omdxy[39] - 9 * omdxy[40] + 6 * omdxy[41] + 6 * omdxy[42] +
                         6 * omdxy[43] + 6 * omdxy[44];
        phidxy[p45[3]] = +9 * omdxy[3] - 18 * omdxy[4] + 3 * omdxy[5] - 27 * omdxy[24] +
                         12 * omdxy[25] + 9 * omdxy[26] + 9 * omdxy[27] - 6 * omdxy[28] -
                         6 * omdxy[29] + 9 * omdxy[36] - 6 * omdxy[37] - 6 * omdxy[38] -
                         27 * omdxy[39] + 9 * omdxy[40] + 12 * omdxy[41] - 6 * omdxy[42] +
                         18 * omdxy[43] - 6 * omdxy[44];
        phidxy[p45[4]] = -18 * omdxy[3] + 84 * omdxy[4] - 18 * omdxy[5] - 6 * omdxy[24] -
                         36 * omdxy[25] + 12 * omdxy[26] - 30 * omdxy[27] + 30 * omdxy[28] -
                         30 * omdxy[36] + 30 * omdxy[38] - 6 * omdxy[39] + 12 * omdxy[40] -
                         36 * omdxy[41] + 24 * omdxy[43];
        phidxy[p45[5]] = +3 * omdxy[3] - 18 * omdxy[4] + 9 * omdxy[5] + 6 * omdxy[24] -
                         18 * omdxy[25] + 3 * omdxy[26] + 6 * omdxy[27] - 9 * omdxy[28] +
                         6 * omdxy[29] + 6 * omdxy[36] + 6 * omdxy[37] - 9 * omdxy[38] +
                         6 * omdxy[39] + 3 * omdxy[40] - 18 * omdxy[41] + 6 * omdxy[42] +
                         6 * omdxy[43] + 6 * omdxy[44];
        phidxy[p45[6]] = +9 * omdxy[6] - 18 * omdxy[7] + 3 * omdxy[8] + 9 * omdxy[24] -
                         6 * omdxy[25] - 6 * omdxy[26] - 27 * omdxy[27] + 9 * omdxy[28] +
                         12 * omdxy[29] + 9 * omdxy[30] - 6 * omdxy[31] - 6 * omdxy[32] -
                         27 * omdxy[33] + 9 * omdxy[34] + 12 * omdxy[35] - 6 * omdxy[42] -
                         6 * omdxy[43] + 18 * omdxy[44];
        phidxy[p45[7]] = -18 * omdxy[6] + 84 * omdxy[7] - 18 * omdxy[8] - 30 * omdxy[24] +
                         30 * omdxy[26] - 6 * omdxy[27] + 12 * omdxy[28] - 36 * omdxy[29] -
                         30 * omdxy[30] + 30 * omdxy[32] - 6 * omdxy[33] + 12 * omdxy[34] -
                         36 * omdxy[35] + 24 * omdxy[44];
        phidxy[p45[8]] = +3 * omdxy[6] - 18 * omdxy[7] + 9 * omdxy[8] + 6 * omdxy[24] +
                         6 * omdxy[25] - 9 * omdxy[26] + 6 * omdxy[27] + 3 * omdxy[28] -
                         18 * omdxy[29] + 6 * omdxy[30] + 6 * omdxy[31] - 9 * omdxy[32] +
                         6 * omdxy[33] + 3 * omdxy[34] - 18 * omdxy[35] + 6 * omdxy[42] +
                         6 * omdxy[43] + 6 * omdxy[44];
        phidxy[p45[9]] = +9 * omdxy[9] - 18 * omdxy[10] + 3 * omdxy[11] - 27 * omdxy[18] +
                         12 * omdxy[19] + 9 * omdxy[20] + 9 * omdxy[21] - 6 * omdxy[22] -
                         6 * omdxy[23] - 3 * omdxy[36] + 18 * omdxy[37] - 6 * omdxy[38] +
                         9 * omdxy[39] - 27 * omdxy[40] + 12 * omdxy[41] - 6 * omdxy[42] +
                         18 * omdxy[43] - 6 * omdxy[44];
        phidxy[p45[10]] = -18 * omdxy[9] + 84 * omdxy[10] - 18 * omdxy[11] - 6 * omdxy[18] -
                          36 * omdxy[19] + 12 * omdxy[20] - 30 * omdxy[21] + 30 * omdxy[22] -
                          12 * omdxy[36] + 36 * omdxy[37] + 6 * omdxy[38] + 12 * omdxy[39] -
                          6 * omdxy[40] - 36 * omdxy[41] - 24 * omdxy[42] + 24 * omdxy[43];
        phidxy[p45[11]] = +3 * omdxy[9] - 18 * omdxy[10] + 9 * omdxy[11] + 6 * omdxy[18] -
                          18 * omdxy[19] + 3 * omdxy[20] + 6 * omdxy[21] - 9 * omdxy[22] +
                          6 * omdxy[23] - 9 * omdxy[36] - 12 * omdxy[37] + 27 * omdxy[38] +
                          3 * omdxy[39] + 6 * omdxy[40] - 18 * omdxy[41] - 18 * omdxy[42] +
                          6 * omdxy[43] + 6 * omdxy[44];
        phidxy[p45[12]] = +9 * omdxy[12] - 18 * omdxy[13] + 3 * omdxy[14] + 9 * omdxy[18] -
                          6 * omdxy[19] - 6 * omdxy[20] - 27 * omdxy[21] + 9 * omdxy[22] +
                          12 * omdxy[23] - 3 * omdxy[30] + 18 * omdxy[31] - 6 * omdxy[32] +
                          9 * omdxy[33] - 27 * omdxy[34] + 12 * omdxy[35] - 6 * omdxy[42] -
                          6 * omdxy[43] + 18 * omdxy[44];
        phidxy[p45[13]] = -18 * omdxy[12] + 84 * omdxy[13] - 18 * omdxy[14] - 30 * omdxy[18] +
                          30 * omdxy[20] - 6 * omdxy[21] + 12 * omdxy[22] - 36 * omdxy[23] -
                          12 * omdxy[30] + 36 * omdxy[31] + 6 * omdxy[32] + 12 * omdxy[33] -
                          6 * omdxy[34] - 36 * omdxy[35] - 24 * omdxy[42] + 24 * omdxy[44];
        phidxy[p45[14]] = +3 * omdxy[12] - 18 * omdxy[13] + 9 * omdxy[14] + 6 * omdxy[18] +
                          6 * omdxy[19] - 9 * omdxy[20] + 6 * omdxy[21] + 3 * omdxy[22] -
                          18 * omdxy[23] - 9 * omdxy[30] - 12 * omdxy[31] + 27 * omdxy[32] +
                          3 * omdxy[33] + 6 * omdxy[34] - 18 * omdxy[35] - 18 * omdxy[42] +
                          6 * omdxy[43] + 6 * omdxy[44];
        phidxy[p45[15]] = +9 * omdxy[15] - 18 * omdxy[16] + 3 * omdxy[17] - 3 * omdxy[18] +
                          18 * omdxy[19] - 6 * omdxy[20] + 9 * omdxy[21] - 27 * omdxy[22] +
                          12 * omdxy[23] - 3 * omdxy[24] + 18 * omdxy[25] - 6 * omdxy[26] +
                          9 * omdxy[27] - 27 * omdxy[28] + 12 * omdxy[29] - 6 * omdxy[42] -
                          6 * omdxy[43] + 18 * omdxy[44];
        phidxy[p45[16]] = -18 * omdxy[15] + 84 * omdxy[16] - 18 * omdxy[17] - 12 * omdxy[18] +
                          36 * omdxy[19] + 6 * omdxy[20] + 12 * omdxy[21] - 6 * omdxy[22] -
                          36 * omdxy[23] - 12 * omdxy[24] + 36 * omdxy[25] + 6 * omdxy[26] +
                          12 * omdxy[27] - 6 * omdxy[28] - 36 * omdxy[29] - 24 * omdxy[43] +
                          24 * omdxy[44];
        phidxy[p45[17]] = +3 * omdxy[15] - 18 * omdxy[16] + 9 * omdxy[17] - 9 * omdxy[18] -
                          12 * omdxy[19] + 27 * omdxy[20] + 3 * omdxy[21] + 6 * omdxy[22] -
                          18 * omdxy[23] - 9 * omdxy[24] - 12 * omdxy[25] + 27 * omdxy[26] +
                          3 * omdxy[27] + 6 * omdxy[28] - 18 * omdxy[29] + 6 * omdxy[42] -
                          18 * omdxy[43] + 6 * omdxy[44];
        phidxy[p45[18]] = +90 * omdxy[18] - 30 * omdxy[19] - 45 * omdxy[20] - 30 * omdxy[21] +
                          15 * omdxy[22] + 30 * omdxy[23] + 15 * omdxy[42] - 45 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[19]] = -30 * omdxy[18] + 120 * omdxy[19] - 30 * omdxy[20] + 30 * omdxy[21] -
                          60 * omdxy[22] + 30 * omdxy[42] - 30 * omdxy[43] + 30 * omdxy[44];
        phidxy[p45[20]] = -45 * omdxy[18] - 30 * omdxy[19] + 90 * omdxy[20] + 15 * omdxy[21] +
                          15 * omdxy[22] - 60 * omdxy[23] + 15 * omdxy[42] - 45 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[21]] = -30 * omdxy[18] + 30 * omdxy[19] + 15 * omdxy[20] + 90 * omdxy[21] -
                          45 * omdxy[22] - 30 * omdxy[23] + 15 * omdxy[42] + 15 * omdxy[43] -
                          45 * omdxy[44];
        phidxy[p45[22]] = +15 * omdxy[18] - 60 * omdxy[19] + 15 * omdxy[20] - 45 * omdxy[21] +
                          90 * omdxy[22] - 30 * omdxy[23] + 15 * omdxy[42] + 15 * omdxy[43] -
                          45 * omdxy[44];
        phidxy[p45[23]] = +30 * omdxy[18] - 60 * omdxy[20] - 30 * omdxy[21] - 30 * omdxy[22] +
                          120 * omdxy[23] + 30 * omdxy[42] + 30 * omdxy[43] - 30 * omdxy[44];
        phidxy[p45[24]] = +90 * omdxy[24] - 30 * omdxy[25] - 45 * omdxy[26] - 30 * omdxy[27] +
                          15 * omdxy[28] + 30 * omdxy[29] + 15 * omdxy[42] - 45 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[25]] = -30 * omdxy[24] + 120 * omdxy[25] - 30 * omdxy[26] + 30 * omdxy[27] -
                          60 * omdxy[28] - 30 * omdxy[42] - 30 * omdxy[43] + 30 * omdxy[44];
        phidxy[p45[26]] = -45 * omdxy[24] - 30 * omdxy[25] + 90 * omdxy[26] + 15 * omdxy[27] +
                          15 * omdxy[28] - 60 * omdxy[29] + 15 * omdxy[42] - 45 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[27]] = -30 * omdxy[24] + 30 * omdxy[25] + 15 * omdxy[26] + 90 * omdxy[27] -
                          45 * omdxy[28] - 30 * omdxy[29] + 15 * omdxy[42] + 15 * omdxy[43] -
                          45 * omdxy[44];
        phidxy[p45[28]] = +15 * omdxy[24] - 60 * omdxy[25] + 15 * omdxy[26] - 45 * omdxy[27] +
                          90 * omdxy[28] - 30 * omdxy[29] + 15 * omdxy[42] + 15 * omdxy[43] -
                          45 * omdxy[44];
        phidxy[p45[29]] = +30 * omdxy[24] - 60 * omdxy[26] - 30 * omdxy[27] - 30 * omdxy[28] +
                          120 * omdxy[29] - 30 * omdxy[42] + 30 * omdxy[43] - 30 * omdxy[44];
        phidxy[p45[30]] = +90 * omdxy[30] - 30 * omdxy[31] - 45 * omdxy[32] - 30 * omdxy[33] +
                          15 * omdxy[34] + 30 * omdxy[35] - 45 * omdxy[42] + 15 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[31]] = -30 * omdxy[30] + 120 * omdxy[31] - 30 * omdxy[32] + 30 * omdxy[33] -
                          60 * omdxy[34] - 30 * omdxy[42] - 30 * omdxy[43] + 30 * omdxy[44];
        phidxy[p45[32]] = -45 * omdxy[30] - 30 * omdxy[31] + 90 * omdxy[32] + 15 * omdxy[33] +
                          15 * omdxy[34] - 60 * omdxy[35] - 45 * omdxy[42] + 15 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[33]] = -30 * omdxy[30] + 30 * omdxy[31] + 15 * omdxy[32] + 90 * omdxy[33] -
                          45 * omdxy[34] - 30 * omdxy[35] + 15 * omdxy[42] + 15 * omdxy[43] -
                          45 * omdxy[44];
        phidxy[p45[34]] = +15 * omdxy[30] - 60 * omdxy[31] + 15 * omdxy[32] - 45 * omdxy[33] +
                          90 * omdxy[34] - 30 * omdxy[35] + 15 * omdxy[42] + 15 * omdxy[43] -
                          45 * omdxy[44];
        phidxy[p45[35]] = +30 * omdxy[30] - 60 * omdxy[32] - 30 * omdxy[33] - 30 * omdxy[34] +
                          120 * omdxy[35] + 30 * omdxy[42] - 30 * omdxy[43] - 30 * omdxy[44];
        phidxy[p45[36]] = +90 * omdxy[36] - 30 * omdxy[37] - 45 * omdxy[38] - 30 * omdxy[39] +
                          15 * omdxy[40] + 30 * omdxy[41] - 45 * omdxy[42] + 15 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[37]] = -30 * omdxy[36] + 120 * omdxy[37] - 30 * omdxy[38] + 30 * omdxy[39] -
                          60 * omdxy[40] - 30 * omdxy[42] + 30 * omdxy[43] - 30 * omdxy[44];
        phidxy[p45[38]] = -45 * omdxy[36] - 30 * omdxy[37] + 90 * omdxy[38] + 15 * omdxy[39] +
                          15 * omdxy[40] - 60 * omdxy[41] - 45 * omdxy[42] + 15 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[39]] = -30 * omdxy[36] + 30 * omdxy[37] + 15 * omdxy[38] + 90 * omdxy[39] -
                          45 * omdxy[40] - 30 * omdxy[41] + 15 * omdxy[42] - 45 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[40]] = +15 * omdxy[36] - 60 * omdxy[37] + 15 * omdxy[38] - 45 * omdxy[39] +
                          90 * omdxy[40] - 30 * omdxy[41] + 15 * omdxy[42] - 45 * omdxy[43] +
                          15 * omdxy[44];
        phidxy[p45[41]] = +30 * omdxy[36] - 60 * omdxy[38] - 30 * omdxy[39] - 30 * omdxy[40] +
                          120 * omdxy[41] + 30 * omdxy[42] - 30 * omdxy[43] - 30 * omdxy[44];
        phidxy[p45[42]] = +90 * omdxy[42] - 30 * omdxy[43] - 30 * omdxy[44];
        phidxy[p45[43]] = -30 * omdxy[42] + 90 * omdxy[43] - 30 * omdxy[44];
        phidxy[p45[44]] = -30 * omdxy[42] - 30 * omdxy[43] + 90 * omdxy[44];

        for (int k = 0; k < 45; ++k) {
          val(k, 0, op_dxy) = phidxy[k].x;
          val(k, 1, op_dxy) = phidxy[k].y;
          val(k, 2, op_dxy) = phidxy[k].z;
        }
      }

      R3 phidxz[45];
      if (whatd & Fop_dxz) {
        phidxz[p45[0]] = +9 * omdxz[0] - 18 * omdxz[1] + 3 * omdxz[2] - 27 * omdxz[30] +
                         12 * omdxz[31] + 9 * omdxz[32] + 9 * omdxz[33] - 6 * omdxz[34] -
                         6 * omdxz[35] - 27 * omdxz[36] + 12 * omdxz[37] + 9 * omdxz[38] +
                         9 * omdxz[39] - 6 * omdxz[40] - 6 * omdxz[41] + 18 * omdxz[42] -
                         6 * omdxz[43] - 6 * omdxz[44];
        phidxz[p45[1]] = -18 * omdxz[0] + 84 * omdxz[1] - 18 * omdxz[2] - 6 * omdxz[30] -
                         36 * omdxz[31] + 12 * omdxz[32] - 30 * omdxz[33] + 30 * omdxz[34] -
                         6 * omdxz[36] - 36 * omdxz[37] + 12 * omdxz[38] - 30 * omdxz[39] +
                         30 * omdxz[40] + 24 * omdxz[42];
        phidxz[p45[2]] = +3 * omdxz[0] - 18 * omdxz[1] + 9 * omdxz[2] + 6 * omdxz[30] -
                         18 * omdxz[31] + 3 * omdxz[32] + 6 * omdxz[33] - 9 * omdxz[34] +
                         6 * omdxz[35] + 6 * omdxz[36] - 18 * omdxz[37] + 3 * omdxz[38] +
                         6 * omdxz[39] - 9 * omdxz[40] + 6 * omdxz[41] + 6 * omdxz[42] +
                         6 * omdxz[43] + 6 * omdxz[44];
        phidxz[p45[3]] = +9 * omdxz[3] - 18 * omdxz[4] + 3 * omdxz[5] - 27 * omdxz[24] +
                         12 * omdxz[25] + 9 * omdxz[26] + 9 * omdxz[27] - 6 * omdxz[28] -
                         6 * omdxz[29] + 9 * omdxz[36] - 6 * omdxz[37] - 6 * omdxz[38] -
                         27 * omdxz[39] + 9 * omdxz[40] + 12 * omdxz[41] - 6 * omdxz[42] +
                         18 * omdxz[43] - 6 * omdxz[44];
        phidxz[p45[4]] = -18 * omdxz[3] + 84 * omdxz[4] - 18 * omdxz[5] - 6 * omdxz[24] -
                         36 * omdxz[25] + 12 * omdxz[26] - 30 * omdxz[27] + 30 * omdxz[28] -
                         30 * omdxz[36] + 30 * omdxz[38] - 6 * omdxz[39] + 12 * omdxz[40] -
                         36 * omdxz[41] + 24 * omdxz[43];
        phidxz[p45[5]] = +3 * omdxz[3] - 18 * omdxz[4] + 9 * omdxz[5] + 6 * omdxz[24] -
                         18 * omdxz[25] + 3 * omdxz[26] + 6 * omdxz[27] - 9 * omdxz[28] +
                         6 * omdxz[29] + 6 * omdxz[36] + 6 * omdxz[37] - 9 * omdxz[38] +
                         6 * omdxz[39] + 3 * omdxz[40] - 18 * omdxz[41] + 6 * omdxz[42] +
                         6 * omdxz[43] + 6 * omdxz[44];
        phidxz[p45[6]] = +9 * omdxz[6] - 18 * omdxz[7] + 3 * omdxz[8] + 9 * omdxz[24] -
                         6 * omdxz[25] - 6 * omdxz[26] - 27 * omdxz[27] + 9 * omdxz[28] +
                         12 * omdxz[29] + 9 * omdxz[30] - 6 * omdxz[31] - 6 * omdxz[32] -
                         27 * omdxz[33] + 9 * omdxz[34] + 12 * omdxz[35] - 6 * omdxz[42] -
                         6 * omdxz[43] + 18 * omdxz[44];
        phidxz[p45[7]] = -18 * omdxz[6] + 84 * omdxz[7] - 18 * omdxz[8] - 30 * omdxz[24] +
                         30 * omdxz[26] - 6 * omdxz[27] + 12 * omdxz[28] - 36 * omdxz[29] -
                         30 * omdxz[30] + 30 * omdxz[32] - 6 * omdxz[33] + 12 * omdxz[34] -
                         36 * omdxz[35] + 24 * omdxz[44];
        phidxz[p45[8]] = +3 * omdxz[6] - 18 * omdxz[7] + 9 * omdxz[8] + 6 * omdxz[24] +
                         6 * omdxz[25] - 9 * omdxz[26] + 6 * omdxz[27] + 3 * omdxz[28] -
                         18 * omdxz[29] + 6 * omdxz[30] + 6 * omdxz[31] - 9 * omdxz[32] +
                         6 * omdxz[33] + 3 * omdxz[34] - 18 * omdxz[35] + 6 * omdxz[42] +
                         6 * omdxz[43] + 6 * omdxz[44];
        phidxz[p45[9]] = +9 * omdxz[9] - 18 * omdxz[10] + 3 * omdxz[11] - 27 * omdxz[18] +
                         12 * omdxz[19] + 9 * omdxz[20] + 9 * omdxz[21] - 6 * omdxz[22] -
                         6 * omdxz[23] - 3 * omdxz[36] + 18 * omdxz[37] - 6 * omdxz[38] +
                         9 * omdxz[39] - 27 * omdxz[40] + 12 * omdxz[41] - 6 * omdxz[42] +
                         18 * omdxz[43] - 6 * omdxz[44];
        phidxz[p45[10]] = -18 * omdxz[9] + 84 * omdxz[10] - 18 * omdxz[11] - 6 * omdxz[18] -
                          36 * omdxz[19] + 12 * omdxz[20] - 30 * omdxz[21] + 30 * omdxz[22] -
                          12 * omdxz[36] + 36 * omdxz[37] + 6 * omdxz[38] + 12 * omdxz[39] -
                          6 * omdxz[40] - 36 * omdxz[41] - 24 * omdxz[42] + 24 * omdxz[43];
        phidxz[p45[11]] = +3 * omdxz[9] - 18 * omdxz[10] + 9 * omdxz[11] + 6 * omdxz[18] -
                          18 * omdxz[19] + 3 * omdxz[20] + 6 * omdxz[21] - 9 * omdxz[22] +
                          6 * omdxz[23] - 9 * omdxz[36] - 12 * omdxz[37] + 27 * omdxz[38] +
                          3 * omdxz[39] + 6 * omdxz[40] - 18 * omdxz[41] - 18 * omdxz[42] +
                          6 * omdxz[43] + 6 * omdxz[44];
        phidxz[p45[12]] = +9 * omdxz[12] - 18 * omdxz[13] + 3 * omdxz[14] + 9 * omdxz[18] -
                          6 * omdxz[19] - 6 * omdxz[20] - 27 * omdxz[21] + 9 * omdxz[22] +
                          12 * omdxz[23] - 3 * omdxz[30] + 18 * omdxz[31] - 6 * omdxz[32] +
                          9 * omdxz[33] - 27 * omdxz[34] + 12 * omdxz[35] - 6 * omdxz[42] -
                          6 * omdxz[43] + 18 * omdxz[44];
        phidxz[p45[13]] = -18 * omdxz[12] + 84 * omdxz[13] - 18 * omdxz[14] - 30 * omdxz[18] +
                          30 * omdxz[20] - 6 * omdxz[21] + 12 * omdxz[22] - 36 * omdxz[23] -
                          12 * omdxz[30] + 36 * omdxz[31] + 6 * omdxz[32] + 12 * omdxz[33] -
                          6 * omdxz[34] - 36 * omdxz[35] - 24 * omdxz[42] + 24 * omdxz[44];
        phidxz[p45[14]] = +3 * omdxz[12] - 18 * omdxz[13] + 9 * omdxz[14] + 6 * omdxz[18] +
                          6 * omdxz[19] - 9 * omdxz[20] + 6 * omdxz[21] + 3 * omdxz[22] -
                          18 * omdxz[23] - 9 * omdxz[30] - 12 * omdxz[31] + 27 * omdxz[32] +
                          3 * omdxz[33] + 6 * omdxz[34] - 18 * omdxz[35] - 18 * omdxz[42] +
                          6 * omdxz[43] + 6 * omdxz[44];
        phidxz[p45[15]] = +9 * omdxz[15] - 18 * omdxz[16] + 3 * omdxz[17] - 3 * omdxz[18] +
                          18 * omdxz[19] - 6 * omdxz[20] + 9 * omdxz[21] - 27 * omdxz[22] +
                          12 * omdxz[23] - 3 * omdxz[24] + 18 * omdxz[25] - 6 * omdxz[26] +
                          9 * omdxz[27] - 27 * omdxz[28] + 12 * omdxz[29] - 6 * omdxz[42] -
                          6 * omdxz[43] + 18 * omdxz[44];
        phidxz[p45[16]] = -18 * omdxz[15] + 84 * omdxz[16] - 18 * omdxz[17] - 12 * omdxz[18] +
                          36 * omdxz[19] + 6 * omdxz[20] + 12 * omdxz[21] - 6 * omdxz[22] -
                          36 * omdxz[23] - 12 * omdxz[24] + 36 * omdxz[25] + 6 * omdxz[26] +
                          12 * omdxz[27] - 6 * omdxz[28] - 36 * omdxz[29] - 24 * omdxz[43] +
                          24 * omdxz[44];
        phidxz[p45[17]] = +3 * omdxz[15] - 18 * omdxz[16] + 9 * omdxz[17] - 9 * omdxz[18] -
                          12 * omdxz[19] + 27 * omdxz[20] + 3 * omdxz[21] + 6 * omdxz[22] -
                          18 * omdxz[23] - 9 * omdxz[24] - 12 * omdxz[25] + 27 * omdxz[26] +
                          3 * omdxz[27] + 6 * omdxz[28] - 18 * omdxz[29] + 6 * omdxz[42] -
                          18 * omdxz[43] + 6 * omdxz[44];
        phidxz[p45[18]] = +90 * omdxz[18] - 30 * omdxz[19] - 45 * omdxz[20] - 30 * omdxz[21] +
                          15 * omdxz[22] + 30 * omdxz[23] + 15 * omdxz[42] - 45 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[19]] = -30 * omdxz[18] + 120 * omdxz[19] - 30 * omdxz[20] + 30 * omdxz[21] -
                          60 * omdxz[22] + 30 * omdxz[42] - 30 * omdxz[43] + 30 * omdxz[44];
        phidxz[p45[20]] = -45 * omdxz[18] - 30 * omdxz[19] + 90 * omdxz[20] + 15 * omdxz[21] +
                          15 * omdxz[22] - 60 * omdxz[23] + 15 * omdxz[42] - 45 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[21]] = -30 * omdxz[18] + 30 * omdxz[19] + 15 * omdxz[20] + 90 * omdxz[21] -
                          45 * omdxz[22] - 30 * omdxz[23] + 15 * omdxz[42] + 15 * omdxz[43] -
                          45 * omdxz[44];
        phidxz[p45[22]] = +15 * omdxz[18] - 60 * omdxz[19] + 15 * omdxz[20] - 45 * omdxz[21] +
                          90 * omdxz[22] - 30 * omdxz[23] + 15 * omdxz[42] + 15 * omdxz[43] -
                          45 * omdxz[44];
        phidxz[p45[23]] = +30 * omdxz[18] - 60 * omdxz[20] - 30 * omdxz[21] - 30 * omdxz[22] +
                          120 * omdxz[23] + 30 * omdxz[42] + 30 * omdxz[43] - 30 * omdxz[44];
        phidxz[p45[24]] = +90 * omdxz[24] - 30 * omdxz[25] - 45 * omdxz[26] - 30 * omdxz[27] +
                          15 * omdxz[28] + 30 * omdxz[29] + 15 * omdxz[42] - 45 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[25]] = -30 * omdxz[24] + 120 * omdxz[25] - 30 * omdxz[26] + 30 * omdxz[27] -
                          60 * omdxz[28] - 30 * omdxz[42] - 30 * omdxz[43] + 30 * omdxz[44];
        phidxz[p45[26]] = -45 * omdxz[24] - 30 * omdxz[25] + 90 * omdxz[26] + 15 * omdxz[27] +
                          15 * omdxz[28] - 60 * omdxz[29] + 15 * omdxz[42] - 45 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[27]] = -30 * omdxz[24] + 30 * omdxz[25] + 15 * omdxz[26] + 90 * omdxz[27] -
                          45 * omdxz[28] - 30 * omdxz[29] + 15 * omdxz[42] + 15 * omdxz[43] -
                          45 * omdxz[44];
        phidxz[p45[28]] = +15 * omdxz[24] - 60 * omdxz[25] + 15 * omdxz[26] - 45 * omdxz[27] +
                          90 * omdxz[28] - 30 * omdxz[29] + 15 * omdxz[42] + 15 * omdxz[43] -
                          45 * omdxz[44];
        phidxz[p45[29]] = +30 * omdxz[24] - 60 * omdxz[26] - 30 * omdxz[27] - 30 * omdxz[28] +
                          120 * omdxz[29] - 30 * omdxz[42] + 30 * omdxz[43] - 30 * omdxz[44];
        phidxz[p45[30]] = +90 * omdxz[30] - 30 * omdxz[31] - 45 * omdxz[32] - 30 * omdxz[33] +
                          15 * omdxz[34] + 30 * omdxz[35] - 45 * omdxz[42] + 15 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[31]] = -30 * omdxz[30] + 120 * omdxz[31] - 30 * omdxz[32] + 30 * omdxz[33] -
                          60 * omdxz[34] - 30 * omdxz[42] - 30 * omdxz[43] + 30 * omdxz[44];
        phidxz[p45[32]] = -45 * omdxz[30] - 30 * omdxz[31] + 90 * omdxz[32] + 15 * omdxz[33] +
                          15 * omdxz[34] - 60 * omdxz[35] - 45 * omdxz[42] + 15 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[33]] = -30 * omdxz[30] + 30 * omdxz[31] + 15 * omdxz[32] + 90 * omdxz[33] -
                          45 * omdxz[34] - 30 * omdxz[35] + 15 * omdxz[42] + 15 * omdxz[43] -
                          45 * omdxz[44];
        phidxz[p45[34]] = +15 * omdxz[30] - 60 * omdxz[31] + 15 * omdxz[32] - 45 * omdxz[33] +
                          90 * omdxz[34] - 30 * omdxz[35] + 15 * omdxz[42] + 15 * omdxz[43] -
                          45 * omdxz[44];
        phidxz[p45[35]] = +30 * omdxz[30] - 60 * omdxz[32] - 30 * omdxz[33] - 30 * omdxz[34] +
                          120 * omdxz[35] + 30 * omdxz[42] - 30 * omdxz[43] - 30 * omdxz[44];
        phidxz[p45[36]] = +90 * omdxz[36] - 30 * omdxz[37] - 45 * omdxz[38] - 30 * omdxz[39] +
                          15 * omdxz[40] + 30 * omdxz[41] - 45 * omdxz[42] + 15 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[37]] = -30 * omdxz[36] + 120 * omdxz[37] - 30 * omdxz[38] + 30 * omdxz[39] -
                          60 * omdxz[40] - 30 * omdxz[42] + 30 * omdxz[43] - 30 * omdxz[44];
        phidxz[p45[38]] = -45 * omdxz[36] - 30 * omdxz[37] + 90 * omdxz[38] + 15 * omdxz[39] +
                          15 * omdxz[40] - 60 * omdxz[41] - 45 * omdxz[42] + 15 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[39]] = -30 * omdxz[36] + 30 * omdxz[37] + 15 * omdxz[38] + 90 * omdxz[39] -
                          45 * omdxz[40] - 30 * omdxz[41] + 15 * omdxz[42] - 45 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[40]] = +15 * omdxz[36] - 60 * omdxz[37] + 15 * omdxz[38] - 45 * omdxz[39] +
                          90 * omdxz[40] - 30 * omdxz[41] + 15 * omdxz[42] - 45 * omdxz[43] +
                          15 * omdxz[44];
        phidxz[p45[41]] = +30 * omdxz[36] - 60 * omdxz[38] - 30 * omdxz[39] - 30 * omdxz[40] +
                          120 * omdxz[41] + 30 * omdxz[42] - 30 * omdxz[43] - 30 * omdxz[44];
        phidxz[p45[42]] = +90 * omdxz[42] - 30 * omdxz[43] - 30 * omdxz[44];
        phidxz[p45[43]] = -30 * omdxz[42] + 90 * omdxz[43] - 30 * omdxz[44];
        phidxz[p45[44]] = -30 * omdxz[42] - 30 * omdxz[43] + 90 * omdxz[44];

        for (int k = 0; k < 45; ++k) {
          val(k, 0, op_dxz) = phidxz[k].x;
          val(k, 1, op_dxz) = phidxz[k].y;
          val(k, 2, op_dxz) = phidxz[k].z;
        }
      }

      R3 phidyz[45];
      if (whatd & Fop_dyz) {
        phidyz[p45[0]] = +9 * omdyz[0] - 18 * omdyz[1] + 3 * omdyz[2] - 27 * omdyz[30] +
                         12 * omdyz[31] + 9 * omdyz[32] + 9 * omdyz[33] - 6 * omdyz[34] -
                         6 * omdyz[35] - 27 * omdyz[36] + 12 * omdyz[37] + 9 * omdyz[38] +
                         9 * omdyz[39] - 6 * omdyz[40] - 6 * omdyz[41] + 18 * omdyz[42] -
                         6 * omdyz[43] - 6 * omdyz[44];
        phidyz[p45[1]] = -18 * omdyz[0] + 84 * omdyz[1] - 18 * omdyz[2] - 6 * omdyz[30] -
                         36 * omdyz[31] + 12 * omdyz[32] - 30 * omdyz[33] + 30 * omdyz[34] -
                         6 * omdyz[36] - 36 * omdyz[37] + 12 * omdyz[38] - 30 * omdyz[39] +
                         30 * omdyz[40] + 24 * omdyz[42];
        phidyz[p45[2]] = +3 * omdyz[0] - 18 * omdyz[1] + 9 * omdyz[2] + 6 * omdyz[30] -
                         18 * omdyz[31] + 3 * omdyz[32] + 6 * omdyz[33] - 9 * omdyz[34] +
                         6 * omdyz[35] + 6 * omdyz[36] - 18 * omdyz[37] + 3 * omdyz[38] +
                         6 * omdyz[39] - 9 * omdyz[40] + 6 * omdyz[41] + 6 * omdyz[42] +
                         6 * omdyz[43] + 6 * omdyz[44];
        phidyz[p45[3]] = +9 * omdyz[3] - 18 * omdyz[4] + 3 * omdyz[5] - 27 * omdyz[24] +
                         12 * omdyz[25] + 9 * omdyz[26] + 9 * omdyz[27] - 6 * omdyz[28] -
                         6 * omdyz[29] + 9 * omdyz[36] - 6 * omdyz[37] - 6 * omdyz[38] -
                         27 * omdyz[39] + 9 * omdyz[40] + 12 * omdyz[41] - 6 * omdyz[42] +
                         18 * omdyz[43] - 6 * omdyz[44];
        phidyz[p45[4]] = -18 * omdyz[3] + 84 * omdyz[4] - 18 * omdyz[5] - 6 * omdyz[24] -
                         36 * omdyz[25] + 12 * omdyz[26] - 30 * omdyz[27] + 30 * omdyz[28] -
                         30 * omdyz[36] + 30 * omdyz[38] - 6 * omdyz[39] + 12 * omdyz[40] -
                         36 * omdyz[41] + 24 * omdyz[43];
        phidyz[p45[5]] = +3 * omdyz[3] - 18 * omdyz[4] + 9 * omdyz[5] + 6 * omdyz[24] -
                         18 * omdyz[25] + 3 * omdyz[26] + 6 * omdyz[27] - 9 * omdyz[28] +
                         6 * omdyz[29] + 6 * omdyz[36] + 6 * omdyz[37] - 9 * omdyz[38] +
                         6 * omdyz[39] + 3 * omdyz[40] - 18 * omdyz[41] + 6 * omdyz[42] +
                         6 * omdyz[43] + 6 * omdyz[44];
        phidyz[p45[6]] = +9 * omdyz[6] - 18 * omdyz[7] + 3 * omdyz[8] + 9 * omdyz[24] -
                         6 * omdyz[25] - 6 * omdyz[26] - 27 * omdyz[27] + 9 * omdyz[28] +
                         12 * omdyz[29] + 9 * omdyz[30] - 6 * omdyz[31] - 6 * omdyz[32] -
                         27 * omdyz[33] + 9 * omdyz[34] + 12 * omdyz[35] - 6 * omdyz[42] -
                         6 * omdyz[43] + 18 * omdyz[44];
        phidyz[p45[7]] = -18 * omdyz[6] + 84 * omdyz[7] - 18 * omdyz[8] - 30 * omdyz[24] +
                         30 * omdyz[26] - 6 * omdyz[27] + 12 * omdyz[28] - 36 * omdyz[29] -
                         30 * omdyz[30] + 30 * omdyz[32] - 6 * omdyz[33] + 12 * omdyz[34] -
                         36 * omdyz[35] + 24 * omdyz[44];
        phidyz[p45[8]] = +3 * omdyz[6] - 18 * omdyz[7] + 9 * omdyz[8] + 6 * omdyz[24] +
                         6 * omdyz[25] - 9 * omdyz[26] + 6 * omdyz[27] + 3 * omdyz[28] -
                         18 * omdyz[29] + 6 * omdyz[30] + 6 * omdyz[31] - 9 * omdyz[32] +
                         6 * omdyz[33] + 3 * omdyz[34] - 18 * omdyz[35] + 6 * omdyz[42] +
                         6 * omdyz[43] + 6 * omdyz[44];
        phidyz[p45[9]] = +9 * omdyz[9] - 18 * omdyz[10] + 3 * omdyz[11] - 27 * omdyz[18] +
                         12 * omdyz[19] + 9 * omdyz[20] + 9 * omdyz[21] - 6 * omdyz[22] -
                         6 * omdyz[23] - 3 * omdyz[36] + 18 * omdyz[37] - 6 * omdyz[38] +
                         9 * omdyz[39] - 27 * omdyz[40] + 12 * omdyz[41] - 6 * omdyz[42] +
                         18 * omdyz[43] - 6 * omdyz[44];
        phidyz[p45[10]] = -18 * omdyz[9] + 84 * omdyz[10] - 18 * omdyz[11] - 6 * omdyz[18] -
                          36 * omdyz[19] + 12 * omdyz[20] - 30 * omdyz[21] + 30 * omdyz[22] -
                          12 * omdyz[36] + 36 * omdyz[37] + 6 * omdyz[38] + 12 * omdyz[39] -
                          6 * omdyz[40] - 36 * omdyz[41] - 24 * omdyz[42] + 24 * omdyz[43];
        phidyz[p45[11]] = +3 * omdyz[9] - 18 * omdyz[10] + 9 * omdyz[11] + 6 * omdyz[18] -
                          18 * omdyz[19] + 3 * omdyz[20] + 6 * omdyz[21] - 9 * omdyz[22] +
                          6 * omdyz[23] - 9 * omdyz[36] - 12 * omdyz[37] + 27 * omdyz[38] +
                          3 * omdyz[39] + 6 * omdyz[40] - 18 * omdyz[41] - 18 * omdyz[42] +
                          6 * omdyz[43] + 6 * omdyz[44];
        phidyz[p45[12]] = +9 * omdyz[12] - 18 * omdyz[13] + 3 * omdyz[14] + 9 * omdyz[18] -
                          6 * omdyz[19] - 6 * omdyz[20] - 27 * omdyz[21] + 9 * omdyz[22] +
                          12 * omdyz[23] - 3 * omdyz[30] + 18 * omdyz[31] - 6 * omdyz[32] +
                          9 * omdyz[33] - 27 * omdyz[34] + 12 * omdyz[35] - 6 * omdyz[42] -
                          6 * omdyz[43] + 18 * omdyz[44];
        phidyz[p45[13]] = -18 * omdyz[12] + 84 * omdyz[13] - 18 * omdyz[14] - 30 * omdyz[18] +
                          30 * omdyz[20] - 6 * omdyz[21] + 12 * omdyz[22] - 36 * omdyz[23] -
                          12 * omdyz[30] + 36 * omdyz[31] + 6 * omdyz[32] + 12 * omdyz[33] -
                          6 * omdyz[34] - 36 * omdyz[35] - 24 * omdyz[42] + 24 * omdyz[44];
        phidyz[p45[14]] = +3 * omdyz[12] - 18 * omdyz[13] + 9 * omdyz[14] + 6 * omdyz[18] +
                          6 * omdyz[19] - 9 * omdyz[20] + 6 * omdyz[21] + 3 * omdyz[22] -
                          18 * omdyz[23] - 9 * omdyz[30] - 12 * omdyz[31] + 27 * omdyz[32] +
                          3 * omdyz[33] + 6 * omdyz[34] - 18 * omdyz[35] - 18 * omdyz[42] +
                          6 * omdyz[43] + 6 * omdyz[44];
        phidyz[p45[15]] = +9 * omdyz[15] - 18 * omdyz[16] + 3 * omdyz[17] - 3 * omdyz[18] +
                          18 * omdyz[19] - 6 * omdyz[20] + 9 * omdyz[21] - 27 * omdyz[22] +
                          12 * omdyz[23] - 3 * omdyz[24] + 18 * omdyz[25] - 6 * omdyz[26] +
                          9 * omdyz[27] - 27 * omdyz[28] + 12 * omdyz[29] - 6 * omdyz[42] -
                          6 * omdyz[43] + 18 * omdyz[44];
        phidyz[p45[16]] = -18 * omdyz[15] + 84 * omdyz[16] - 18 * omdyz[17] - 12 * omdyz[18] +
                          36 * omdyz[19] + 6 * omdyz[20] + 12 * omdyz[21] - 6 * omdyz[22] -
                          36 * omdyz[23] - 12 * omdyz[24] + 36 * omdyz[25] + 6 * omdyz[26] +
                          12 * omdyz[27] - 6 * omdyz[28] - 36 * omdyz[29] - 24 * omdyz[43] +
                          24 * omdyz[44];
        phidyz[p45[17]] = +3 * omdyz[15] - 18 * omdyz[16] + 9 * omdyz[17] - 9 * omdyz[18] -
                          12 * omdyz[19] + 27 * omdyz[20] + 3 * omdyz[21] + 6 * omdyz[22] -
                          18 * omdyz[23] - 9 * omdyz[24] - 12 * omdyz[25] + 27 * omdyz[26] +
                          3 * omdyz[27] + 6 * omdyz[28] - 18 * omdyz[29] + 6 * omdyz[42] -
                          18 * omdyz[43] + 6 * omdyz[44];
        phidyz[p45[18]] = +90 * omdyz[18] - 30 * omdyz[19] - 45 * omdyz[20] - 30 * omdyz[21] +
                          15 * omdyz[22] + 30 * omdyz[23] + 15 * omdyz[42] - 45 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[19]] = -30 * omdyz[18] + 120 * omdyz[19] - 30 * omdyz[20] + 30 * omdyz[21] -
                          60 * omdyz[22] + 30 * omdyz[42] - 30 * omdyz[43] + 30 * omdyz[44];
        phidyz[p45[20]] = -45 * omdyz[18] - 30 * omdyz[19] + 90 * omdyz[20] + 15 * omdyz[21] +
                          15 * omdyz[22] - 60 * omdyz[23] + 15 * omdyz[42] - 45 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[21]] = -30 * omdyz[18] + 30 * omdyz[19] + 15 * omdyz[20] + 90 * omdyz[21] -
                          45 * omdyz[22] - 30 * omdyz[23] + 15 * omdyz[42] + 15 * omdyz[43] -
                          45 * omdyz[44];
        phidyz[p45[22]] = +15 * omdyz[18] - 60 * omdyz[19] + 15 * omdyz[20] - 45 * omdyz[21] +
                          90 * omdyz[22] - 30 * omdyz[23] + 15 * omdyz[42] + 15 * omdyz[43] -
                          45 * omdyz[44];
        phidyz[p45[23]] = +30 * omdyz[18] - 60 * omdyz[20] - 30 * omdyz[21] - 30 * omdyz[22] +
                          120 * omdyz[23] + 30 * omdyz[42] + 30 * omdyz[43] - 30 * omdyz[44];
        phidyz[p45[24]] = +90 * omdyz[24] - 30 * omdyz[25] - 45 * omdyz[26] - 30 * omdyz[27] +
                          15 * omdyz[28] + 30 * omdyz[29] + 15 * omdyz[42] - 45 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[25]] = -30 * omdyz[24] + 120 * omdyz[25] - 30 * omdyz[26] + 30 * omdyz[27] -
                          60 * omdyz[28] - 30 * omdyz[42] - 30 * omdyz[43] + 30 * omdyz[44];
        phidyz[p45[26]] = -45 * omdyz[24] - 30 * omdyz[25] + 90 * omdyz[26] + 15 * omdyz[27] +
                          15 * omdyz[28] - 60 * omdyz[29] + 15 * omdyz[42] - 45 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[27]] = -30 * omdyz[24] + 30 * omdyz[25] + 15 * omdyz[26] + 90 * omdyz[27] -
                          45 * omdyz[28] - 30 * omdyz[29] + 15 * omdyz[42] + 15 * omdyz[43] -
                          45 * omdyz[44];
        phidyz[p45[28]] = +15 * omdyz[24] - 60 * omdyz[25] + 15 * omdyz[26] - 45 * omdyz[27] +
                          90 * omdyz[28] - 30 * omdyz[29] + 15 * omdyz[42] + 15 * omdyz[43] -
                          45 * omdyz[44];
        phidyz[p45[29]] = +30 * omdyz[24] - 60 * omdyz[26] - 30 * omdyz[27] - 30 * omdyz[28] +
                          120 * omdyz[29] - 30 * omdyz[42] + 30 * omdyz[43] - 30 * omdyz[44];
        phidyz[p45[30]] = +90 * omdyz[30] - 30 * omdyz[31] - 45 * omdyz[32] - 30 * omdyz[33] +
                          15 * omdyz[34] + 30 * omdyz[35] - 45 * omdyz[42] + 15 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[31]] = -30 * omdyz[30] + 120 * omdyz[31] - 30 * omdyz[32] + 30 * omdyz[33] -
                          60 * omdyz[34] - 30 * omdyz[42] - 30 * omdyz[43] + 30 * omdyz[44];
        phidyz[p45[32]] = -45 * omdyz[30] - 30 * omdyz[31] + 90 * omdyz[32] + 15 * omdyz[33] +
                          15 * omdyz[34] - 60 * omdyz[35] - 45 * omdyz[42] + 15 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[33]] = -30 * omdyz[30] + 30 * omdyz[31] + 15 * omdyz[32] + 90 * omdyz[33] -
                          45 * omdyz[34] - 30 * omdyz[35] + 15 * omdyz[42] + 15 * omdyz[43] -
                          45 * omdyz[44];
        phidyz[p45[34]] = +15 * omdyz[30] - 60 * omdyz[31] + 15 * omdyz[32] - 45 * omdyz[33] +
                          90 * omdyz[34] - 30 * omdyz[35] + 15 * omdyz[42] + 15 * omdyz[43] -
                          45 * omdyz[44];
        phidyz[p45[35]] = +30 * omdyz[30] - 60 * omdyz[32] - 30 * omdyz[33] - 30 * omdyz[34] +
                          120 * omdyz[35] + 30 * omdyz[42] - 30 * omdyz[43] - 30 * omdyz[44];
        phidyz[p45[36]] = +90 * omdyz[36] - 30 * omdyz[37] - 45 * omdyz[38] - 30 * omdyz[39] +
                          15 * omdyz[40] + 30 * omdyz[41] - 45 * omdyz[42] + 15 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[37]] = -30 * omdyz[36] + 120 * omdyz[37] - 30 * omdyz[38] + 30 * omdyz[39] -
                          60 * omdyz[40] - 30 * omdyz[42] + 30 * omdyz[43] - 30 * omdyz[44];
        phidyz[p45[38]] = -45 * omdyz[36] - 30 * omdyz[37] + 90 * omdyz[38] + 15 * omdyz[39] +
                          15 * omdyz[40] - 60 * omdyz[41] - 45 * omdyz[42] + 15 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[39]] = -30 * omdyz[36] + 30 * omdyz[37] + 15 * omdyz[38] + 90 * omdyz[39] -
                          45 * omdyz[40] - 30 * omdyz[41] + 15 * omdyz[42] - 45 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[40]] = +15 * omdyz[36] - 60 * omdyz[37] + 15 * omdyz[38] - 45 * omdyz[39] +
                          90 * omdyz[40] - 30 * omdyz[41] + 15 * omdyz[42] - 45 * omdyz[43] +
                          15 * omdyz[44];
        phidyz[p45[41]] = +30 * omdyz[36] - 60 * omdyz[38] - 30 * omdyz[39] - 30 * omdyz[40] +
                          120 * omdyz[41] + 30 * omdyz[42] - 30 * omdyz[43] - 30 * omdyz[44];
        phidyz[p45[42]] = +90 * omdyz[42] - 30 * omdyz[43] - 30 * omdyz[44];
        phidyz[p45[43]] = -30 * omdyz[42] + 90 * omdyz[43] - 30 * omdyz[44];
        phidyz[p45[44]] = -30 * omdyz[42] - 30 * omdyz[43] + 90 * omdyz[44];

        for (int k = 0; k < 45; ++k) {
          val(k, 0, op_dyz) = phidyz[k].x;
          val(k, 1, op_dyz) = phidyz[k].y;
          val(k, 2, op_dyz) = phidyz[k].z;
        }
      }
    }
  }

  // Auxiliary Finite Elements for the partition of unity of DDM methods

  // FE space Edge03ds0 useful to build a partition of unity for the FE space Edge03d of
  // edge elements of degree 1 in 3d
  // It has only the interpolation operator, it doesn't have basis functions
  class TypeOfFE_P0Edge3ds0 : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh::Element Element;
    typedef Element::Rd Rd;
    typedef Element::RdHat RdHat;
    static const int d = Rd::d;
    struct A4 {
      int dfon[4];

      A4( ) {
        dfon[0] = dfon[1] = dfon[2] = dfon[3] = 0;
        dfon[1] = 1;
      }

      operator const int *( ) const { return dfon; }
    };

    TypeOfFE_P0Edge3ds0( );

    void FB(const What_d whatd, const Mesh &Th, const Element &K, const RdHat &PHat,
            RNMK_ &val) const;

    ~TypeOfFE_P0Edge3ds0( ) {}

   private:
    TypeOfFE_P0Edge3ds0(const TypeOfFE_P0Edge3ds0 &);
    void operator=(const TypeOfFE_P0Edge3ds0 &);
  };

  TypeOfFE_P0Edge3ds0::TypeOfFE_P0Edge3ds0( )
    : GTypeOfFE< Mesh3 >(A4( ), 1, 1, Element::ne, Element::ne, false, true) {
    R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.), R3(0., 0., 1.)};

    for (int e = 0; e < Element::ne; ++e) {
      this->PtInterpolation[e] = 0.5 * (Pt[Element::nvedge[e][0]] + Pt[Element::nvedge[e][1]]);
    }

    {
      for (int e = 0; e < Element::ne; e++) {
        this->pInterpolation[e] = e;
        this->cInterpolation[e] = 0;
        this->dofInterpolation[e] = e;
        this->coefInterpolation[e] = 1.;
      }
    }
  }

  void TypeOfFE_P0Edge3ds0::FB(const What_d whatd, const Mesh &, const Element &K,
                               const RdHat &PHat, RNMK_ &val) const {
    assert(0);
  }

  // FE space Edge13ds0 useful to build a partition of unity for the FE space Edge13d of
  // edge elements of degree 2 in 3d
  // It has only the interpolation operator, it doesn't have basis functions
  class TypeOfFE_P1Edge3ds0 : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh::Element Element;
    typedef Element::Rd Rd;
    typedef Element::RdHat RdHat;
    static const int d = Rd::d;
    struct A4 {
      int dfon[4];

      A4( ) {
        dfon[0] = dfon[1] = dfon[2] = dfon[3] = 0;
        dfon[1] = 2;
        dfon[2] = 2;
      }

      operator const int *( ) const { return dfon; }
    };

    TypeOfFE_P1Edge3ds0( );

    void FB(const What_d whatd, const Mesh &Th, const Element &K, const RdHat &PHat,
            RNMK_ &val) const;

    ~TypeOfFE_P1Edge3ds0( ) {}

   private:
    TypeOfFE_P1Edge3ds0(const TypeOfFE_P1Edge3ds0 &);
    void operator=(const TypeOfFE_P1Edge3ds0 &);
  };

  TypeOfFE_P1Edge3ds0::TypeOfFE_P1Edge3ds0( )
    : GTypeOfFE< Mesh3 >(A4( ), 1, 1, Element::ne * 2 + Element::nf * 2, Element::ne + Element::nf,
                         false, true) {
    R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.), R3(0., 0., 1.)};

    {
      int i = 0;

      for (int e = 0; e < Element::ne; ++e, ++i) {
        this->PtInterpolation[i] = 0.5 * (Pt[Element::nvedge[e][0]] + Pt[Element::nvedge[e][1]]);
      }

      for (int f = 0; f < Element::nf; ++f, ++i) {
        this->PtInterpolation[i] =
          1. / 3. *
          (Pt[Element::nvface[f][0]] + Pt[Element::nvface[f][1]] + Pt[Element::nvface[f][2]]);
      }
    }
    int i = 0, p = 0;
    int e;

    for (e = 0; e < (Element::ne)*2; e++) {
      if (e % 2 == 1) {
        p = p - 1;
      }

      this->pInterpolation[i] = p;
      this->cInterpolation[i] = 0.;
      this->dofInterpolation[i] = e;
      this->coefInterpolation[i] = 1.;
      i++;
      p++;
    }

    for (int f = 0; f < (Element::nf)*2; f++) {
      if (f % 2 == 1) {
        p = p - 1;
      }

      this->pInterpolation[i] = p;
      this->cInterpolation[i] = 0.;
      this->dofInterpolation[i] = e + f;
      this->coefInterpolation[i] = 1.;
      i++;
      p++;
    }
  }

  void TypeOfFE_P1Edge3ds0::FB(const What_d whatd, const Mesh &, const Element &K,
                               const RdHat &PHat, RNMK_ &val) const {
    assert(0);
  }

  // FE space Edge23ds0 useful to build a partition of unity for the FE space Edge23d of
  // edge elements of degree 3 in 3d
  // It has only the interpolation operator, it doesn't have basis functions
  class TypeOfFE_P2Edge3ds0 : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    typedef GFElement< Mesh3 > FElement;

    static int dfon[];
    static const int d = Mesh::Rd::d;
    struct A4 {
      int dfon[4];

      A4( ) {
        dfon[0] = dfon[1] = dfon[2] = dfon[3] = 0;
        dfon[1] = 3;
        dfon[2] = 6;
        dfon[3] = 3;
      }

      operator const int *( ) const { return dfon; }
    };

    TypeOfFE_P2Edge3ds0( );    // constructor

    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
            RNMK_ &val) const;

    ~TypeOfFE_P2Edge3ds0( ) {}

   private:
    TypeOfFE_P2Edge3ds0(const TypeOfFE_P2Edge3ds0 &);
    void operator=(const TypeOfFE_P2Edge3ds0 &);
  };

  TypeOfFE_P2Edge3ds0::TypeOfFE_P2Edge3ds0( )
    : GTypeOfFE< Mesh3 >(A4( ), 1, 1, 3 * Element::ne + 6 * Element::nf + 3,
                         Element::ne + Element::nf + 1, false, true) {
    R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.),
               R3(0., 0., 1.)};    // 4 ref tetrahedron vertices

    {
      // We build the interpolation pts on the edges of the reference tetrahedron:
      int i = 0;

      for (int e = 0; e < Element::ne; ++e, ++i) {
        this->PtInterpolation[i] = 0.5 * (Pt[Element::nvedge[e][0]] + Pt[Element::nvedge[e][1]]);
      }

      // We build the interpolation pts on the faces of the reference tetrahedron:
      // (the index i mustn't be reinitialised!)
      for (int f = 0; f < Element::nf; ++f, ++i) {
        this->PtInterpolation[i] =
          1. / 3. *
          (Pt[Element::nvface[f][0]] + Pt[Element::nvface[f][1]] + Pt[Element::nvface[f][2]]);
      }

      // We build the interpolation pts on the tetrahedron:
      // (the index i mustn't be reinitialised!)
      this->PtInterpolation[i] = 1. / 4. * (Pt[0] + Pt[1] + Pt[2] + Pt[3]);
    }
    // We build the indices in (13.1) : edge dofs
    int i = 0, p = 0;    // i is the k in (13.1) (chapter 13 ff++doc)
    int e;               // we will need e also below, in the parts referred to faces and volumes

    for (e = 0; e < (Element::ne)*3; e++) {    // loop on the 18 edge dofs
      if (e % 3 != 0) {
        p = p - 1;
      }    // for 3 successive basis fcts the quad pts are the same (they correspond to the same
           // edge)

      {
        this->pInterpolation[i] = p;      // pk in (13.1)
        this->cInterpolation[i] = 0;      // jk in (13.1)
        this->dofInterpolation[i] = e;    // ik in (13.1)
        this->coefInterpolation[i] = 1.;
        i++;
        p++;
      }
    }

    // We build the indices in (13.1) : face dofs
    // (the indices i and p mustn't be reinitialised)
    int f;    // we will need f also below, in the part referred to volumes

    for (f = 0; f < (Element::nf)*6; f++) {    // loop on the 24 face dofs
      if (f % 6 != 0) {
        p = p - 1;
      }    // for 6 successive basis fcts the quad pts are the same (they correspond to the same
           // face)

      {
        this->pInterpolation[i] = p;          // pk in (13.1)
        this->cInterpolation[i] = 0;          // jk in (13.1)
        this->dofInterpolation[i] = e + f;    // ik in (13.1)
        this->coefInterpolation[i] = 1.;
        i++;
        p++;
      }
    }

    // We build the indices in (13.1) : volume basis fncts
    // (the indices i and p mustn't be reinitialised)
    for (int v = 0; v < 3; v++) {    // loop on the 3 volume dofs
      if (v > 0) {
        p = p - 1;
      }

      {
        this->pInterpolation[i] = p;              // pk in (13.1)
        this->cInterpolation[i] = 0;              // jk in (13.1)
        this->dofInterpolation[i] = e + f + v;    // ik in (13.1)
        this->coefInterpolation[i] = 1.;
        i++;
        p++;
      }
    }
  }

  void TypeOfFE_P2Edge3ds0::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                               const RdHat &PHat, RNMK_ &val) const {
    assert(0);
  }

  class TypeOfFE_RT1_3d : public GTypeOfFE< Mesh3 > {
   public:
    typedef Mesh3 Mesh;
    typedef Mesh3::Element Element;
    typedef GFElement< Mesh3 > FElement;
    static int dfon[];
    static const int d = Mesh::Rd::d;
    // quadrature Formula on a face
    static const GQuadratureFormular< R2 > QFface;
    // quadrature Formula on an element
    static const GQuadratureFormular< R3 > QFtetra;
    TypeOfFE_RT1_3d( );
    int edgeface[4][3];
    void FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K, const RdHat &PHat,
            RNMK_ &val) const;
    void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
             int *nump) const;
  };

  int TypeOfFE_RT1_3d::dfon[] = {0, 0, 3, 3};    // dofs per vertice, edge, face, volume

  // Quadrature formula on a face,
  const GQuadratureFormular< R2 > TypeOfFE_RT1_3d::QFface(QuadratureFormular_T_5);

  // Quadrature formula on the tetraedron
  const GQuadratureFormular< R3 > TypeOfFE_RT1_3d::QFtetra(QuadratureFormular_Tet_2);

  TypeOfFE_RT1_3d::TypeOfFE_RT1_3d( )
    : GTypeOfFE< Mesh3 >(dfon, d, 1, 3 * QFface.n * 3 * Element::nf + 3 * QFtetra.n * 3,
                         Element::nf * QFface.n + QFtetra.n, false, true) {
    assert(QFface.n);
    assert(QFtetra.n);
    // 4 ref tetra vertices
    R3 Pt[] = {R3(0., 0., 0.), R3(1., 0., 0.), R3(0., 1., 0.), R3(0., 0., 1.)};

    // We build the interpolation pts on the faces
    int p = 0;

    for (int f = 0; f < Element::nf; ++f) {
      for (int q = 0; q < QFface.n; ++q, ++p) {
        double x = QFface[q].x;
        double y = QFface[q].y;
        this->PtInterpolation[p] = Pt[Element::nvface[f][0]] * (1. - x - y) +
                                   Pt[Element::nvface[f][1]] * x + Pt[Element::nvface[f][2]] * y;
      }
    }

    // We build the interpolation bubble pts
    for (int q = 0; q < QFtetra.n; ++q, ++p) {
      double x = QFtetra[q].x;
      double y = QFtetra[q].y;
      double z = QFtetra[q].z;
      this->PtInterpolation[p] = Pt[0] * (1. - x - y - z) + Pt[1] * x + Pt[2] * y + Pt[3] * z;
    }

    int i = 0;
    p = 0;

    for (int f = 0; f < Element::nf; f++) {        // loop on the 4 face dofs
      for (int q = 0; q < QFface.n; ++q, p++) {    // loop on the face quadrature pts
        for (int df = 0; df < 3; df++) {           // 3 dof par face
          int dof = 3 * f + df;                    // numero  du dof  3 dof / face

          for (int c = 0; c < 3; c++, i++) {    // loop on the 3 components
            this->pInterpolation[i] = p;        // pk
            this->cInterpolation[i] = c;        // jk
            this->dofInterpolation[i] = dof;    // ik
            this->coefInterpolation[i] = 0.;
          }
        }
      }
    }

    {
      p = Element::nf * QFface.n;

      for (int q = 0; q < QFtetra.n; ++q, ++p) {    // loop on the volume quadrature pts
        for (int v = 12; v < 15; v++) {             // loop on the 3 volume dofs
          for (int c = 0; c < 3; c++, i++) {        // loop on the 3 components
            this->pInterpolation[i] = p;            // pk
            this->cInterpolation[i] = c;            // jk
            this->dofInterpolation[i] = v;          // ik
            this->coefInterpolation[i] = 0.;
          }
        }
      }
    }
    // verif bonne taille
    ffassert(p == this->PtInterpolation.N( ));
  }    // end TypeOfFE_RT1_3d()

  // For the coefficients of interpolation alphak in (13.1)

  void TypeOfFE_RT1_3d::set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M,
                            int ocoef, int odf, int *nump) const {
    int i = ocoef;

    // *******************************************
    // DOFs on the 4 faces --- 3 DOFs per face //
    // *******************************************

    // inv of int(F) lamda_i lambda_j = Area(K) * (d! (1+delta ij) / (d+max n)! = 0.5 * 2(1+delta
    // ij) /(2+2)!
    double c1[][3] = {{9, -3, -3} /* 0 */, {-3, 9, -3} /* 1 */, {-3, -3, 9} /* 2 */};
    R3 NK[4];

    K.Gradlambda(NK);    // InteriorNormal of F /h = - N  3*|K|/|F|
    double coefK = -3. * K.mesure( );

    for (int ff = 0; ff < Element::nf; ff++) {    // loop on the 4 face dofs
      const Element::Vertex *fV[3] = {&K.at(Element::nvface[ff][0]), &K.at(Element::nvface[ff][1]),
                                      &K.at(Element::nvface[ff][2])};
      int p[] = {0, 1, 2};
      int fp = K.facePermutation(ff);
      if (fp & 1) {
        Exchange(p[0], p[1]);
      }

      if (fp & 2) {
        Exchange(p[1], p[2]);
      }

      if (fp & 4) {
        Exchange(p[0], p[1]);
      }

      R3 N = NK[ff];
      N *= coefK * K.faceOrient(ff);

      for (int q = 0; q < QFface.n; ++q) {    // loop on the face quadrature pts
        // the 3 lambdas of P1 on the face

        double lambda[3] = {1. - QFface[q].x - QFface[q].y, QFface[q].x, QFface[q].y};
        R sa = QFface[q].a;
        R cp[3] = {(c1[0][0] * lambda[0] + c1[0][1] * lambda[1] + c1[0][2] * lambda[2]) * sa,
                   (c1[1][0] * lambda[0] + c1[1][1] * lambda[1] + c1[1][2] * lambda[2]) * sa,
                   (c1[2][0] * lambda[0] + c1[2][1] * lambda[1] + c1[2][2] * lambda[2]) * sa};

        for (int idof = 0; idof < 3; idof++) {    // loop sur 3 dof de la face
          for (int c = 0; c < 3; c++, i++) {      // loop on the 3 components
            M.coef[i] = N[c] * cp[p[idof]];
          }
        }
      }
    }

    // cout << endl;
    // **************************************************
    // DOFs on the tetraedra --- 3 DOFs in the volume //
    // **************************************************
    // Base Piola compatible  B_i =N_i* |F_i|/6  for 3 face 1,2,3
    double CK = -K.mesure( );    // dof U= [u1,u2,u3] > |K| int_K ( B_i.U )

    for (int p = 0; p < QFtetra.n; ++p) {
      double w = QFtetra[p].a * CK;

      for (int l = 0; l < 3; l++) {
        M.coef[i++] = w * NK[l + 1].x;
        M.coef[i++] = w * NK[l + 1].y;
        M.coef[i++] = w * NK[l + 1].z;
      }
    }

  }    // end set function

  // here the basis functions and theirs derivates
  void TypeOfFE_RT1_3d::FB(const What_d whatd, const Mesh &Th, const Mesh3::Element &K,
                           const RdHat &PHat, RNMK_ &val) const {
    assert(val.N( ) >= 15);
    assert(val.M( ) == 3);

    val = 0;

    // basis functions to RT03d / multiply by sign to have a exterior normal and divide by the
    // mesure of K phi = signe * (x - qi)/ (volume*d)
    R cc = d * K.mesure( );
    R lambda[] = {1. - PHat.sum( ), PHat.x, PHat.y, PHat.z};
    R3 X = K(PHat);
    R3 phi[4] = {X - K[0], X - K[1], X - K[2], X - K[3]};    // phi * area *6

    // fo contain just the sign about permutation ----- 1perm=-1 / 2perm=1 / 3perm=-1
    double fo[4] = {(double)K.faceOrient(0), (double)K.faceOrient(1), (double)K.faceOrient(2),
                    (double)K.faceOrient(3)};
    int p[15] = {2, 1, 0,  3,  4,  5,  8, 7,
                 6, 9, 10, 11, 12, 13, 14};    // Permutation for orientation to dof
    R3 Pm[16];                                 // all the momome function ..

    for (int ff = 0, k = 0; ff < Element::nf; ff++, k += 3) {
      // orientation de la face a envert
      int fp = K.facePermutation(ff);
      if (fp & 1) {
        Exchange(p[k], p[k + 1]);
      }

      if (fp & 2) {
        Exchange(p[k + 1], p[k + 2]);
      }

      if (fp & 4) {
        Exchange(p[k], p[k + 1]);
      }
    }

    double cf[][3] = {
      {-1.25, 0.25, 0} /* 0 */,  {-1.25, 0, 0.25} /* 1 */,   {-1.5, -0.25, -0.25} /* 2 */,
      {0.25, -1.25, 0} /* 3 */,  {0, -1.25, 0.25} /* 4 */,   {-0.25, -1.5, -0.25} /* 5 */,
      {0.25, 0, -1.25} /* 6 */,  {0, 0.25, -1.25} /* 7 */,   {-0.25, -0.25, -1.5} /* 8 */,
      {1.5, 1.25, 1.25} /* 9 */, {1.25, 1.5, 1.25} /* 10 */, {1.25, 1.25, 1.5} /* 11 */,
      {-15, 15, 0} /* 12 */,     {-15, 0, 15} /* 13 */,      {-30, -15, -15} /* 14 */};
    int Bii[][2] = {{0, 0} /* 0 */,  {0, 1} /* 1 */,  {0, 2} /* 2 */,  {0, 3} /* 3 */,
                    {1, 0} /* 4 */,  {1, 1} /* 5 */,  {1, 2} /* 6 */,  {1, 3} /* 7 */,
                    {2, 0} /* 8 */,  {2, 1} /* 9 */,  {2, 2} /* 10 */, {2, 3} /* 11 */,
                    {3, 0} /* 12 */, {3, 1} /* 13 */, {3, 2} /* 14 */, {3, 3} /* 15 */};
    int fe[] = {1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14};
    int k6[] = {0, 5, 10};

    // here we build all the monomials phi_i lambda_j
    for (int l = 0; l < 16; ++l) {
      int i = Bii[l][0];
      int j = Bii[l][1];
      Pm[l] = phi[i] * lambda[j] / cc;
    }

    // caution to the numbering of the nodes of the face to the numbering of the basis function
    // cci

    double sg[15] = {fo[0], fo[0], fo[0], fo[1], fo[1], fo[1], fo[2], fo[2],
                     fo[2], fo[3], fo[3], fo[3], 1.,    1.,    1.};

    if (whatd & Fop_D0) {
      for (int pdf = 0; pdf < 15; ++pdf) {
        int df = p[pdf];
        R3 fd(0., 0., 0.);
        if (df < 12) {
          fd = Pm[fe[df]];    // edge function ..
        }

        for (int k = 0; k < 3; ++k) {
          fd += cf[df][k] * Pm[k6[k]];
        }

        fd *= sg[pdf];
        val(pdf, 0, op_id) = fd.x;
        val(pdf, 1, op_id) = fd.y;
        val(pdf, 2, op_id) = fd.z;
      }
    }

    if (whatd & Fop_D1) {
      R3 DL[4];
      K.Gradlambda(DL);
      R3 Dphix(1, 0, 0);
      R3 Dphiy(0, 1, 0);
      R3 Dphiz(0, 0, 1);
      R3 DxPm[16];
      R3 DyPm[16];
      R3 DzPm[16];

      for (int l = 0; l < 16; ++l) {
        // diff phi[i]*lambda[j]/cc;
        int i = Bii[l][0];
        int j = Bii[l][1];
        R Lj = lambda[j];
        R3 DLj = DL[j];
        R3 DF1 = (Dphix * Lj + phi[i].x * DLj) / cc;
        R3 DF2 = (Dphiy * Lj + phi[i].y * DLj) / cc;
        R3 DF3 = (Dphiz * Lj + phi[i].z * DLj) / cc;

        DxPm[l] = R3(DF1.x, DF2.x, DF3.x);
        DyPm[l] = R3(DF1.y, DF2.y, DF3.y);
        DzPm[l] = R3(DF1.z, DF2.z, DF3.z);
      }

      if (whatd & Fop_dx) {
        for (int pdf = 0; pdf < 15; ++pdf) {
          int df = p[pdf];
          R3 fd(0., 0., 0.);
          if (df < 12) {
            fd = DxPm[fe[df]];    // edge function ..
          }

          for (int k = 0; k < 3; ++k) {
            fd += cf[df][k] * DxPm[k6[k]];
          }

          fd *= sg[df];
          val(pdf, 0, op_dx) = fd.x;
          val(pdf, 1, op_dx) = fd.y;
          val(pdf, 2, op_dx) = fd.z;
        }
      }

      if (whatd & Fop_dy) {
        for (int pdf = 0; pdf < 15; ++pdf) {
          int df = p[pdf];
          R3 fd(0., 0., 0.);
          if (df < 12) {
            fd = DyPm[fe[df]];    // edge function ..
          }

          for (int k = 0; k < 3; ++k) {
            fd += cf[df][k] * DyPm[k6[k]];
          }

          fd *= sg[df];
          val(pdf, 0, op_dy) = fd.x;
          val(pdf, 1, op_dy) = fd.y;
          val(pdf, 2, op_dy) = fd.z;
        }
      }

      if (whatd & Fop_dz) {
        for (int pdf = 0; pdf < 15; ++pdf) {
          int df = p[pdf];
          R3 fd(0., 0., 0.);
          if (df < 12) {
            fd = DzPm[fe[df]];    // edge function ..
          }

          for (int k = 0; k < 3; ++k) {
            fd += cf[df][k] * DzPm[k6[k]];
          }

          fd *= sg[df];
          val(pdf, 0, op_dz) = fd.x;
          val(pdf, 1, op_dz) = fd.y;
          val(pdf, 2, op_dz) = fd.z;
        }
      }

      if (whatd & Fop_D2) {
        cout << " to do FH RT2 dxx, dyy, dzz, dxy, dxz, dyz " << endl;
      }
    }
  }    // end basis functions declaration

 
 // RT1 surf

   class TypeOfFE_RT1_surf : public GTypeOfFE<MeshS>  {
    public:
    typedef MeshS Mesh;
    typedef MeshS::Element Element;
    typedef GFElement< MeshS > FElement;
    static int dfon[];
    static const int d = Mesh::Rd::d;
    // quadrature Formula on an edge
    static const QuadratureFormular1d QFE;
    //static const GQuadratureFormular< R2 > QFface;
    // quadrature Formula on an element
    static const GQuadratureFormular< R2 > QFK;
    //static const GQuadratureFormular< R3 > QFtetra;
    TypeOfFE_RT1_surf( );
    void FB(const What_d whatd, const Mesh &Th, const Element &K, const RdHat &PHat,
            RNMK_ &val) const;
    void set(const Mesh &Th, const Element &K, InterpolationMatrix< RdHat > &M, int ocoef, int odf,
             int *nump) const;
    };

    int TypeOfFE_RT1_surf::dfon[] = {0, 2, 2, 0};    // dofs per vertice, edge, face, volume
    // Quadrature formula on a edge,
    const QuadratureFormular1d TypeOfFE_RT1_surf::QFE(3, 2, GaussLegendre(2), true); // (QF_GaussLegendre2);
    // Quadrature formula on the triangle
    const GQuadratureFormular< R2 > TypeOfFE_RT1_surf::QFK(QuadratureFormular_T_5);
    TypeOfFE_RT1_surf::TypeOfFE_RT1_surf( )
        : GTypeOfFE< MeshS >(dfon, d, 2, 3 * QFE.n * 2 * Element::ne + 2 * QFK.n * 3,
                                 Element::ne * QFE.n + QFK.n, false, true) {
             
        assert(QFK.n);
        assert(QFE.n);
        // 3 ref triangle vertices
        R2 Pt[] = {R2(0., 0.), R2(1., 0.), R2(0., 1.)};
            
        // the interpolation pts on the faces
        int p = 0;
        for (int e = 0; e < Element::ne; ++e) {
          for (int q = 0; q < QFE.n; ++q,++p) {
            RdHat A(TriangleHat[VerticesOfTriangularEdge[e][0]]);   // a remplacer
            RdHat B(TriangleHat[VerticesOfTriangularEdge[e][1]]);   // a remplacer
                 
            this->PtInterpolation[p] = Pt[Element::nvedge[e][0]] * (1. - QFE[q].x) +
                 Pt[Element::nvedge[e][1]] * QFE[q].x;
          }
        }
        // the interpolation bubble pts
            for (int q = 0; q < QFK.n; ++q, ++p)
                this->PtInterpolation[p] = QFK[q];
           
        int i = 0;
        p = 0;

        for (int e = 0; e < Element::ne; e++) {        // loop on the 3 edges
          for (int q = 0; q < QFE.n; ++q, p++) {    // loop on the face quadrature pts
            for (int df = 0; df < 2; df++) {           // 2 dof per edge
              int dof = 2 * e + df;                    // numero  du dof  2 dof / edge

              for (int c = 0; c < 3; c++, i++) {    // loop on the 3 components
                this->pInterpolation[i] = p;        // pk
                this->cInterpolation[i] = c;        // jk
                this->dofInterpolation[i] = dof;    // ik
                this->coefInterpolation[i] = 0.;
              }
            }
          }
        }

        p = Element::ne * QFE.n;

        for (int q = 0; q < QFK.n; ++q, ++p) {    // loop on the volume quadrature pts
          for (int v = 6; v < 8; v++) {             // loop on the 2 volume dofs
            for (int c = 0; c < 3; c++, i++) {        // loop on the 3 components
              this->pInterpolation[i] = p;            // pk
              this->cInterpolation[i] = c;            // jk
              this->dofInterpolation[i] = v;          // ik
              this->coefInterpolation[i] = 0.;
            }
          }
        }
        ffassert(p == this->PtInterpolation.N( ));
      }

      
      void TypeOfFE_RT1_surf::set(const Mesh & Th,const Element & K,InterpolationMatrix<RdHat> & M ,int ocoef,int odf,int *nump) const
      {
          int i = ocoef;

         
          // 2dofs per edges
          for (int e = 0; e < 3; e++) {
              R3 N=K.N(e);
            R s = K.EdgeOrientation(e) ;
            for (int q = 0; q < QFE.n; ++q) {
                R l0 = QFE[q].x, l1 = 1 - QFE[q].x;
                R p[2] = { (2 * l0 - l1) * 2, (2 * l1 - l0) * 2 } ;     // poly othogonaux to \lambda_1, to \lambda_0
                R cp[2] = { s * p[0] * QFE[q].a, s * p[1] * QFE[q].a};
           
              if (s < 0)
                swap(cp[0], cp[1]);
              for (int idof = 0; idof < 2; idof++) {
                for (int c = 0; c < 3; c++, i++)
                M.coef[i] = N[c] * cp[idof];
                   }
            }
          }

         // 2 dofs on the triangle
         double CK = -K.mesure( );    // dof U= [u1,u2,u3] > |K| int_K ( B_i.U )
          R3 N[2] = {K.N(1), K.N(2)};
         for (int q = 0; q < QFK.n; ++q) {
           double w = QFK[q].a * CK;
           for (int l = 0; l < 2; ++l) {
             M.coef[i++] = w * N[l].x;
             M.coef[i++] = w * N[l].y;
             M.coef[i++] = w * N[l].z;
           }
         }
          
      }
  
      
      void  TypeOfFE_RT1_surf::FB(const What_d whatd,const Mesh & Th,const MeshS::Element & K,const RdHat &PHat, RNMK_ & val) const
      {
          assert(val.N( ) >= 8);
          assert(val.M() == 3);
        
          val = 0;
          R3 A(K[0]),B(K[1]),C(K[2]);
          R l[]={1.-PHat.sum(),PHat.x,PHat.y};
          R3 D[3];
          K.Gradlambda(D);
          
          //loc basis
          R3 AB(A,B),AC(A,C),BA(B,A),BC(B,C),CA(C,A),CB(C,B);
          int eo[3] = {K.EdgeOrientation(0), K.EdgeOrientation(1), K.EdgeOrientation(2)};
          R c[3];
          c[0] = 1./(((AB,D[1]) + (AC,D[2])) *K.mesure()) ; //c[0]*= eo[0] ;
          c[1] = 1./(((BA,D[0]) + (BC,D[2])) *K.mesure()) ; //c[1]*= eo[1];
          c[2] = 1./(((CA,D[0]) + (CB,D[1])) *K.mesure()) ; //c[2]*= eo[2];

          R3 phi[3];
          phi[0] = AB*(l[1]*c[0]) + AC*(l[2]*c[0]);
          phi[1] = BA*(l[0]*c[1]) + BC*(l[2]*c[1]);
          phi[2] = CA*(l[0]*c[2]) + CB*(l[1]*c[2]);
          
          int pI[8][3];    // store p_k
          int lI[8][3];    // store l_k
          R cI[8][3];      // store c_k
          int df = 0;
          R CKK = 2 * K.mesure();

          for (int e = 0; e < 3; ++e) {
            int i = e;
            int ii[2] = {(e + 1) % 3, (e + 2) % 3};
            R s = eo[e] / CKK;
            if (s < 0)
              swap(ii[0], ii[1]);

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

          R s8 = 8 / CKK, s01 = s8;
          R cbb[] = {s8, 2 * s01, -s01, s8};    // { [ 8, 16], [ -8, 8] }

          // the 2 bubbles
          for (int k = 0; k < 2; ++k, df++) {
            pI[df][0] = 0;
            lI[df][0] = 0;
            cI[df][0] = cbb[k];

            pI[df][1] = 1;    // i
            lI[df][1] = 1;
            cI[df][1] = cbb[k + 2];

            pI[df][2] = 2;
            lI[df][2] = 2;
            cI[df][2] = 0;
          }

          ffassert(df == 8);

          if (whatd & Fop_D0){
            for (int df = 0; df < 8; ++df) {
              R3 fd(0., 0., 0.);
              for (int k = 0; k < 3; ++k)
                fd += (cI[df][k] * l[lI[df][k]]) * phi[pI[df][k]];
              val(df, 0, op_id) = fd.x;
              val(df, 1, op_id) = fd.y;
              val(df, 2, op_id) = fd.z;
            }
          }

          if (whatd & Fop_D1) {
            R3 D[3];
            K.Gradlambda(D);
            if (whatd & Fop_dx) {
              R3 Dphix[3] = { AB*(D[1].x*c[0]) + AC*(D[2].x*c[0]),
                              BA*(D[0].x*c[1]) + BC*(D[2].x*c[1]),
                              CA*(D[0].x*c[2]) + CB*(D[1].x*c[2])};
              for (int df = 0; df < 8; ++df) {
                R3 fd(0., 0.,0.);

                for (int k = 0; k < 3; ++k)
                  fd += cI[df][k] * (D[lI[df][k]].x * phi[pI[df][k]] + l[lI[df][k]] * Dphix[k]);
                val(df, 0, op_dx) = fd.x;
                val(df, 1, op_dx) = fd.y;
                val(df, 2, op_dx) = fd.z;
              }
            }

            if (whatd & Fop_dy) {
              R3 Dphiy[3] = { AB*(D[1].y*c[0]) + AC*(D[2].y*c[0]),
                              BA*(D[0].y*c[1]) + BC*(D[2].y*c[1]),
                              CA*(D[0].y*c[2]) + CB*(D[1].y*c[2])};
              for (int df = 0; df < 8; ++df) {
                R3 fd(0., 0.,0.);

                for (int k = 0; k < 3; ++k)
                  fd += cI[df][k] * (D[lI[df][k]].y * phi[pI[df][k]] + l[lI[df][k]] * Dphiy[k]);
                val(df, 0, op_dy) = fd.x;
                val(df, 1, op_dy) = fd.y;
                val(df, 2, op_dy) = fd.z;
              }
            }

            if (whatd & Fop_dz) {
              R3 Dphiz[3] = { AB*(D[1].z*c[0]) + AC*(D[2].z*c[0]),
                              BA*(D[0].z*c[1]) + BC*(D[2].z*c[1]),
                              CA*(D[0].z*c[2]) + CB*(D[1].z*c[2])};
              for (int df = 0; df < 8; ++df) {
                R3 fd(0., 0.,0.);

              for (int k = 0; k < 3; ++k)
                  fd += cI[df][k] * (D[lI[df][k]].z * phi[pI[df][k]] + l[lI[df][k]] * Dphiz[k]);
              
              val(df, 0, op_dz) = fd.x;
              val(df, 1, op_dz) = fd.y;
              val(df, 2, op_dz) = fd.z;
            }
          }
        }
        if (whatd & Fop_D2)
          ffassert(0);
      
    }
 
 
 // add new FEs for DDM to build the partition of unity
  static TypeOfFE_P2Edge3ds0 Elm_P2Edge3ds0;
  // a static variable to add the finite element to freefem++
  static AddNewFE3 P2Edge3ds0("Edge23ds0", &Elm_P2Edge3ds0);
  // link with FreeFem++
  static TypeOfFE_P1Edge3ds0 Elm_P1Edge3ds0;
  // a static variable to add the finite element to freefem++
  static AddNewFE3 P1Edge3ds0("Edge13ds0", &Elm_P1Edge3ds0);
  // link with FreeFem++
  static TypeOfFE_P0Edge3ds0 Elm_P0Edge3ds0;
  // a static variable to add the finite element to freefem++
  static AddNewFE3 P0Edge3d("Edge03ds0", &Elm_P0Edge3ds0);
  static TypeOfFE_Edge1_3d Edge1_3d;    // TypeOfFE_Edge1_3d is the name of the class we defined
  GTypeOfFE< Mesh3 > &GEdge13d(Edge1_3d);    // GTypeOfFE<Mesh3> is the mother class
  static AddNewFE3 TypeOfFE_Edge1_3d("Edge13d",
                                     &GEdge13d);    // Edge13d will be the name used by the user
  static TypeOfFE_Edge2_3d Edge2_3d;
  GTypeOfFE< Mesh3 > &GEdge23d(Edge2_3d);
  static AddNewFE3 TypeOfFE_Edge2_3d("Edge23d", &GEdge23d);
  static TypeOfFE_RT1_3d RT1_3d;
  GTypeOfFE< Mesh3 > &RT13d(RT1_3d);
  static AddNewFE3 TypeOfFE_RT1_3d("RT13d", &RT13d);

  static TypeOfFE_RT1_surf RT1_S;
  GTypeOfFE< MeshS > &RT1S(RT1_S);
  static AddNewFES TypeOfFE_RT1_surf("RT1S", &RT1S);



}    // namespace Fem2D

// --- fin --
