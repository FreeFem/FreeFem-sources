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

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

#include "ff++.hpp"

template< class Mesh >
class MatrixEdgeP1 : public E_F0 {
 public:
  typedef Matrice_Creuse< R > *Result;
  Expression emat, expTh;
  MatrixEdgeP1(const basicAC_F0 &args) {
    args.SetNameParam( );
    emat = args[0];
    expTh = to< const Mesh * >(args[1]);
  }

  MatrixEdgeP1( ) {}

  static ArrayOfaType typeargs( ) {
    return ArrayOfaType(atype< Matrice_Creuse< R > * >( ), atype< const Mesh * >( ));
  }

  static E_F0 *f(const basicAC_F0 &args) { return new MatrixEdgeP1(args); }

  AnyType operator( )(Stack s) const;
};

// petit Pb de compatibitile
template< class Mesh >
int fnvedge(const Mesh *th, int i, int j) {
  ffassert(0);
  return 0;
}
template<>
int fnvedge< Mesh >(const Mesh *th, int i, int j) {
  return (i + 1 + j) % 3;
}
template<>
int fnvedge< Mesh3 >(const Mesh3 *th, int i, int j) {
  return Mesh3::Element::nvedge[i][j];
}
template< class Mesh >
AnyType MatrixEdgeP1< Mesh >::operator( )(Stack stack) const {

  typedef typename Mesh::Element Element;
  typedef typename Mesh::RdHat RdHat;
  static const int nvedgeTet[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
  static const int nvedgeTria[3][2] = {{1, 2}, {2, 0}, {0, 1}};
  static const int nvedgeSeg[1][2] = {{0, 1}};

  const int d = RdHat::d;
  ffassert(d == 2 || d == 3 || d == 1);
  const int nbedgeE = d * (d + 1) / 2;
  const int(*const nvedge)[2] = (d == 1) ? nvedgeSeg : (d == 2 ? nvedgeTria : nvedgeTet);
  Matrice_Creuse< R > *sparce_mat = GetAny< Matrice_Creuse< R > * >((*emat)(stack));
  MatriceMorse< R > *amorse = 0;
  MeshPoint *mp(MeshPointStack(stack)), mps = *mp;
  typedef const Mesh *pmesh;
  const Mesh *pTh = GetAny< pmesh >((*expTh)(stack));
  ffassert(pTh);
  const Mesh &Th(*pTh);
  {

    // const edge ...
    HashTable< SortArray< int, 2 >, int > e(nbedgeE * Th.nt, Th.nv);
    int ne = 0;
    for (int k = 0; k < Th.nt; ++k)
      for (int i = 0; i < nbedgeE; ++i) {
        int i0 = nvedge[i][0], i1 = nvedge[i][1];
        SortArray< int, 2 > eki(Th(k, i0), Th(k, i1));
        if (!e.find(eki)) e.add(eki, ne++);
      }
    if (verbosity > 2 && mpirank == 0)
      cout << " ne = " << ne << " " << nbedgeE << " " << Th.nv << endl;
    MatriceMorse< R > *pAij = new MatriceMorse< R >(ne, Th.nv), &Aij = *pAij;
    // ffassert(Th.ne==ne);
    for (int k = 0; k < ne; k++) {

      int i1 = e.t[k].k.v[0];
      int i2 = e.t[k].k.v[1];
      Aij(k, i1) = -1;
      Aij(k, i2) = +1;
    }

    amorse = pAij;    // new MatriceMorse<R>(Th.nv, Th.nv, Aij, false);
  }
  sparce_mat->Uh = UniqueffId( );
  sparce_mat->Vh = UniqueffId( );
  sparce_mat->A.master(amorse);
  sparce_mat->typemat =
    0;    //(amorse->n == amorse->m) ? TypeSolveMat(TypeSolveMat::GMRES) :
          // TypeSolveMat(TypeSolveMat::NONESQUARE);// none square matrice (morse)
  *mp = mps;

  if (verbosity > 3) {
    cout << "  End Build MatEdgeP1 : " << endl;
  }

  return sparce_mat;
}

static void Load_Init( ) {
  if (verbosity > 4) cout << " load: init Mat Edge 1 " << endl;
  Global.Add("MatrixEdgeP1", "(", new OneOperatorCode< MatrixEdgeP1< Mesh > >( ));
  Global.Add("MatrixEdgeP1", "(", new OneOperatorCode< MatrixEdgeP1< Mesh3 > >( ));
}

LOADFUNC(Load_Init)
