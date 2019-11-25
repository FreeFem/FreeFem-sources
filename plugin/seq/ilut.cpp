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
// SUMMARY : ILUT plugin for FreeFem++ wrapping GMM++ functions
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Alessandro Proverbio
// David Radice
// E-MAIL  : ...

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep: gmm
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include <cmath>
#include <iostream>
#include "AFunction.hpp"
#include "RNM.hpp"
#include "error.hpp"
#include <gmm/gmm.h>
#include <vector>

#define ILUT_K_FILLIN 5
#define ILUT_EPS 1e-6

#define PRINT(VAR) cout << VAR << endl

using namespace std;
using namespace gmm;

typedef ilut_precond< row_matrix< rsvector< double > > > my_ilut_precond;

class ILUT;

class ILUT_Matrix {
 private:
  long *_i;
  long *_j;
  double *_c;
  long _nelem;
  long _size;

 public:
  ILUT_Matrix(KN< long > *const &i, KN< long > *const &j, KN< double > *const &c)
    : _i(*i), _j(*j), _c(*c), _nelem(c->N( )) {
    _size = max(i->max( ), j->max( ));
    ++_size;
  }

  friend class ILUT;
};

class ILUT_Vector {
 private:
  double *_v;
  long _size;

 public:
  ILUT_Vector(KN< double > *const &c) : _v(*c), _size(c->N( )) {}

  friend class ILUT;
};

class ILUT {
 private:
  static my_ilut_precond *p;
  static long size;

 public:
  static long make_ilut_precond(ILUT_Matrix const &m) {
    row_matrix< rsvector< double > > A(m._size, m._size);
    row_matrix< wsvector< double > > w_A(m._size, m._size);

    for (long k(0); k < m._nelem; ++k) {
      w_A[m._i[k]][m._j[k]] = m._c[k];
    }

    copy(w_A, A);    // A <-- w_A
    delete p;
    p = new my_ilut_precond(A, ILUT_K_FILLIN, ILUT_EPS);

    size = m._size;

    return 0;
  }

  static void apply_ilut_precond(ILUT_Vector const &v, KN< double > *const &x) {
    vector< double > vv(size);
    vector< double > xx(size);

    for (long k = 0; k < size; ++k) {
      vv[k] = v._v[k];
    }

    mult(*p, vv, xx);    // xx <-- p.solve(vv)

    for (long k = 0; k < size; ++k) {
      (*x)[k] = xx[k];
    }

    // If used for the full vector fill the remaining components
    for (long k = 0; k + size < x->N( ); ++k) {
      (*x)[k + size] = v._v[k + size];
    }
  }
};

my_ilut_precond *ILUT::p = 0;
long ILUT::size = 0;

long *make_ilut_precond_eq(long *const &errorcode, ILUT_Matrix const &mat) {
  *errorcode = ILUT::make_ilut_precond(mat);
  return errorcode;
}

KN< double > *apply_ilut_precond_eq(KN< double > *const &x, ILUT_Vector const &vec) {
  ILUT::apply_ilut_precond(vec, x);
  return x;
}

ILUT_Matrix make_ilut_precond(KN< long > *const &i, KN< long > *const &j, KN< double > *const &v) {
  return ILUT_Matrix(i, j, v);
}

ILUT_Vector apply_ilut_precond(KN< double > *const &v) { return ILUT_Vector(v); }

static void Load_Init( ) {
  if (verbosity) {
    cout << " -- load ilut init : " << endl;
  }

  Dcl_Type< ILUT_Matrix >( );
  Dcl_Type< ILUT_Vector >( );
  Global.Add("applyIlutPrecond", "(",
             new OneOperator1_< ILUT_Vector, KN< double > * >(apply_ilut_precond));
  Global.Add("makeIlutPrecond", "(",
             new OneOperator3_< ILUT_Matrix, KN< long > *, KN< long > *, KN< double > * >(
               make_ilut_precond));
  TheOperators->Add("=", new OneOperator2_< long *, long *, ILUT_Matrix >(make_ilut_precond_eq));
  TheOperators->Add(
    "=", new OneOperator2_< KN< double > *, KN< double > *, ILUT_Vector >(apply_ilut_precond_eq));
}

LOADFUNC(Load_Init)
