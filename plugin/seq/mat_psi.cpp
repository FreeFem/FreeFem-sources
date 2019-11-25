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

// Example C++ function "myfunction", dynamically loaded into "load.edp"

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include "ff++.hpp"

class MatrixUpWind0 : public E_F0 {
 public:
  typedef Matrice_Creuse< R > *Result;
  Expression emat, expTh, expc, expu1, expu2;
  MatrixUpWind0(const basicAC_F0 &args) {
    args.SetNameParam( );
    emat = args[0];
    expTh = to< pmesh >(args[1]);
    expc = CastTo< double >(args[2]);
    const E_Array *a = dynamic_cast< const E_Array * >((Expression)args[3]);
    if (!a || a->size( ) != 2) {
      CompileError("syntax:  MatrixUpWind0(Th,rhi,[u1,u2])");
    }

    int err = 0;
    expu1 = CastTo< double >((*a)[0]);
    expu2 = CastTo< double >((*a)[1]);
  }

  ~MatrixUpWind0( ) {}

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< Matrice_Creuse< R > * >( ), atype< pmesh >( ), atype< double >( ),
                        atype< E_Array >( ));
  }

  static E_F0 *f(const basicAC_F0 &args) { return new MatrixUpWind0(args); }

  AnyType operator( )(Stack s) const;
};

class MatrixUpWind3 : public E_F0 {
 public:
  typedef Matrice_Creuse< R > *Result;
  Expression emat, expTh, expc, expu1, expu2, expu3;
  MatrixUpWind3(const basicAC_F0 &args) {
    args.SetNameParam( );
    emat = args[0];
    expTh = to< pmesh3 >(args[1]);
    expc = CastTo< double >(args[2]);
    const E_Array *a = dynamic_cast< const E_Array * >((Expression)args[3]);
    if (a == NULL) printf("Dynamic cast failed\n");
    if (a->size( ) != 3) {
      CompileError("syntax:  MatrixUpWind0(Th,rhi,[u1,u2])");
    }

    int err = 0;
    expu1 = CastTo< double >((*a)[0]);
    expu2 = CastTo< double >((*a)[1]);
    expu3 = CastTo< double >((*a)[2]);
  }

  ~MatrixUpWind3( ) {}

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< Matrice_Creuse< R > * >( ), atype< pmesh3 >( ), atype< double >( ),
                        atype< E_Array >( ));
  }

  static E_F0 *f(const basicAC_F0 &args) { return new MatrixUpWind3(args); }

  AnyType operator( )(Stack s) const;
};

int gladys(double q[3][2], double u[2], double c[3], double a[3][3]) {    // PSI Deconninck
  // computes matrix a on a triangle for the Chacon-Reina Petrof-Galerkin upwind
  // working arrays
  double dw[3][2];                 // basis function gradients times  area
  double ua[2], kk[3], beta[3];    // to define a[][]
  double udc = 0;                  // u.grad(w)*area
  bool oneaval = false;
  int i1 = -1;

  for (int i = 0; i < 3; i++) {
    int ip = (i + 1) % 3, ipp = (ip + 1) % 3;

    for (int j = 0; j < 2; j++) {
      dw[i][1 - j] = (2 * j - 1) * (q[ipp][j] - q[ip][j]) / 2;
    }
  }

  for (int i = 0; i < 3; i++) {
    kk[i] = u[0] * dw[i][0] + u[1] * dw[i][1];
    udc += kk[i] * c[i];
  }

  for (int i = 0; i < 3; i++) {
    ua[0] = u[0];
    ua[1] = u[1];
    int ip = (i + 1) % 3, ipp = (ip + 1) % 3;
    if (kk[i] > 0 && kk[ip] <= 0 && kk[ipp] <= 0) {
      beta[i] = 1;
      beta[ip] = 0;
      beta[ipp] = 0;
      oneaval = true;
    } else if (kk[i] <= 0 && kk[ip] > 0 && kk[ipp] > 0) {
      i1 = i;
    }
  }

  if (!oneaval) {
    if (i1 < 0) {
      cout << "bug\n";
    }

    int i = i1, ip = (i + 1) % 3, ipp = (i + 2) % 3;
    double lambda = (c[ip] - c[i]) * (c[ipp] - c[i]);
    if (fabs(lambda) < -1e-20) {
      return 0;
    }

    if (lambda < 0) {
      if (udc > 0) {
        beta[i] = 0;
        beta[ip] = 0;
        beta[ipp] = 1;
        ua[0] = udc * (q[ipp][0] - q[i][0]) / (c[ipp] - c[i]);
        ua[1] = udc * (q[ipp][1] - q[i][1]) / (c[ipp] - c[i]);
      } else {
        beta[i] = 0;
        beta[ipp] = 0;
        beta[ip] = 1;
        ua[0] = udc * (q[ip][0] - q[i][0]) / (c[ip] - c[i]);
        ua[1] = udc * (q[ip][1] - q[i][1]) / (c[ip] - c[i]);
      }
    } else {
      beta[i] = 0;
      beta[ip] = kk[ip] * (c[ip] - c[i]) / udc;
      beta[ipp] = kk[ipp] * (c[ipp] - c[i]) / udc;
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      a[i][j] = beta[i] * (ua[0] * dw[j][0] + ua[1] * dw[j][1]);
    }
  }

  return 1;
}

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
      double a[3][3];
      if (gladys(q, u, c, a)) {
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

int Marco(const Mesh3::Element &K, R3 U, R c[4], double a[4][4]) {    // PSI Deconninck
  ExecError("Not Implemented Sorry Marco!");
  return 0;
}

AnyType MatrixUpWind3::operator( )(Stack stack) const {
  Matrice_Creuse< R > *sparce_mat = GetAny< Matrice_Creuse< R > * >((*emat)(stack));
  MatriceMorse< R > *amorse = 0;
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  const Mesh3 *pTh = GetAny< pmesh3 >((*expTh)(stack));
  ffassert(pTh);
  const Mesh3 &Th(*pTh);
  {
    MatriceMorse< R > *pAij = new MatriceMorse< R >(Th.nv), &Aij = *pAij;

    KN< double > cc(Th.nv);
    double infini = DBL_MAX;
    cc = infini;

    for (int it = 0; it < Th.nt; it++) {
      for (int iv = 0; iv < 4; iv++) {
        int i = Th(it, iv);
        if (cc[i] == infini) {    // if nuset the set
          mp->setP(&Th, it, iv);
          cc[i] = GetAny< double >((*expc)(stack));
        }
      }
    }

    for (int k = 0; k < Th.nt; k++) {
      const Mesh3::Element &K(Th[k]);
      const Mesh3::Vertex &A(K[0]), &B(K[1]), &C(K[2]), &D(K[3]);
      R3 Pt(1. / 4., 1. / 4., 1. / 4.);
      R3 U;
      MeshPointStack(stack)->set(Th, K(Pt), Pt, K, K.lab);
      U.x = GetAny< R >((*expu1)(stack));
      U.y = GetAny< R >((*expu2)(stack));
      U.z = GetAny< R >((*expu3)(stack));

      int ii[4] = {Th(A), Th(B), Th(C), Th(D)};    // number of 4 vertex
      double c[4] = {cc[ii[0]], cc[ii[1]], cc[ii[2]], cc[ii[3]]};
      double a[4][4];
      if (Marco(K, U, c, a)) {
        for (int i = 0; i < 4; i++) {
          for (int j = 0; j < 4; j++) {
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
  Global.Add("MatUpWind0", "(", new OneOperatorCode< MatrixUpWind0 >( ));
  Global.Add("MatUpWind0", "(", new OneOperatorCode< MatrixUpWind3 >( ));
}

LOADFUNC(Load_Init)
