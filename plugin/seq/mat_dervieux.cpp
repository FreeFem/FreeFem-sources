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
// SUMMARY : Implementation of P1-P0 FVM-FEM
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

// compile and link with ./load.link  mat\_dervieux.cpp

#include <iostream>
using namespace std;

#include "cfloat"
#include "ff++.hpp"

class MatrixUpWind0 : public E_F0mps {
 public:
  typedef Matrice_Creuse< R > *Result;
  Expression emat, expTh, expc, expu1, expu2;
  MatrixUpWind0(const basicAC_F0 &args) {
    args.SetNameParam( );
    emat = args[0];                      // the matrix expression
    expTh = to< pmesh >(args[1]);        // a the expression to get the mesh
    expc = CastTo< double >(args[2]);    // the expression to get c  (must be a double)
    // a array expression [ a, b]
    const E_Array *a = dynamic_cast< const E_Array * >((Expression)args[3]);
    if (a->size( ) != 2) {
      CompileError("syntax:  MatrixUpWind0(Th,rhi,[u1,u2])");
    }

    int err = 0;
    expu1 = CastTo< double >((*a)[0]);    // fist exp of the array (must be a  double)
    expu2 = CastTo< double >((*a)[1]);    // second exp of the array (must be a  double)
  }

  ~MatrixUpWind0( ) {}

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< Matrice_Creuse< R > * >( ), atype< pmesh >( ), atype< double >( ),
                        atype< E_Array >( ));
  }

  static E_F0 *f(const basicAC_F0 &args) { return new MatrixUpWind0(args); }

  AnyType operator( )(Stack s) const;
};

int fvmP1P0(double q[3][2], double u[2], double c[3], double a[3][3],
            double where[3]) {    // computes matrix a on a triangle for the Dervieux FVM
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      a[i][j] = 0;
    }
  }

  for (int i = 0; i < 3; i++) {
    int ip = (i + 1) % 3, ipp = (ip + 1) % 3;
    double unL =
      -((q[ip][1] + q[i][1] - 2 * q[ipp][1]) * u[0] - (q[ip][0] + q[i][0] - 2 * q[ipp][0]) * u[1]) /
      6;
    if (unL > 0) {
      a[i][i] += unL;
      a[ip][i] -= unL;
    } else {
      a[i][ip] += unL;
      a[ip][ip] -= unL;
    }

    if (where[i] && where[ip]) {    // this is a boundary edge
      unL = ((q[ip][1] - q[i][1]) * u[0] - (q[ip][0] - q[i][0]) * u[1]) / 2;
      if (unL > 0) {
        a[i][i] += unL;
        a[ip][ip] += unL;
      }
    }
  }

  return 1;
}

// the evaluation routine
AnyType MatrixUpWind0::operator( )(Stack stack) const {
  Matrice_Creuse< R > *sparce_mat = GetAny< Matrice_Creuse< R > * >((*emat)(stack));
  MatriceMorse< R > *amorse = 0;
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  const Mesh *pTh = GetAny< pmesh >((*expTh)(stack));
  ffassert(pTh);
  const Mesh &Th(*pTh);
  {
    MatriceMorse< R > *pAij = new MatriceMorse< R >(Th.nv), &Aij = *pAij;

    KN< double > cc(Th.nv);
    double infini = DBL_MAX;
    cc = infini;

    for (int it = 0; it < Th.nt; it++) {
      for (int iv = 0; iv < 3; iv++) {
        int i = Th(it, iv);
        if (cc[i] == infini) {    // if nuset the set
          mp->setP(&Th, it, iv);
          cc[i] = GetAny< double >((*expc)(stack));
        }
      }
    }

    for (int k = 0; k < Th.nt; k++) {
      const Triangle &K(Th[k]);
      const Vertex &A(K[0]), &B(K[1]), &C(K[2]);
      R2 Pt(1. / 3., 1. / 3.);
      R u[2];
      MeshPointStack(stack)->set(Th, K(Pt), Pt, K, K.lab);
      u[0] = GetAny< R >((*expu1)(stack));
      u[1] = GetAny< R >((*expu2)(stack));

      int ii[3] = {Th(A), Th(B), Th(C)};
      double q[3][2] = {{A.x, A.y}, {B.x, B.y}, {C.x, C.y}};    // coordinates of 3 vertices (input)
      double c[3] = {cc[ii[0]], cc[ii[1]], cc[ii[2]]};
      double a[3][3], where[3] = {(double)A.lab, (double)B.lab, (double)C.lab};
      if (fvmP1P0(q, u, c, a, where)) {
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            if (fabs(a[i][j]) >= 1e-30) {
              Aij[make_pair(ii[i], ii[j])] += a[i][j];
            }
          }
        }
      }
    }

    amorse = pAij;    // new MatriceMorse<R>(Th.nv, Th.nv, Aij, false);
  }
  sparce_mat->Uh = UniqueffId( );
  sparce_mat->Vh = UniqueffId( );
  sparce_mat->A.master(amorse);
  sparce_mat->typemat =
    0;    //(amorse->n == amorse->m) ? TypeSolveMat(TypeSolveMat::GMRES) :
          //TypeSolveMat(TypeSolveMat::NONESQUARE);// none square matrice (morse)
  *mp = mps;

  if (verbosity > 3) {
    cout << "  End Build MatrixUpWind : " << endl;
  }

  return sparce_mat;
}

static void Load_Init( ) {
  cout << " lood: init Mat Chacon " << endl;
  Global.Add("MatUpWind1", "(", new OneOperatorCode< MatrixUpWind0 >( ));
}

LOADFUNC(Load_Init)
