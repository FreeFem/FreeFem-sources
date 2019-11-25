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

// Example C++ function "CppModTemplate" dynamically loaded into "load.edp"

// *INDENT-OFF* //
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
// *INDENT-ON* //

#include <ff++.hpp>
#include "AFunction_ext.hpp"    // Extension of "AFunction.hpp" to deal with more than 3 parameters function
using namespace Fem2D;

// see src/femlib/RNM.hpp

// dummy routine to understand how to use vector
double CppModTemplate3(KN< double > *const &A,                              // OUTPUT
                       KN< double > *const &B, KN< double > *const &C) {    // INPUTS
  // Remarque:
  // It might prove usefull to have a look in the cpp file where KN is defined: src/femlib/RNM.hpp
  //
  // To access value at node i of vector N, do as follow: *(N[0]+i)
  // Explanation (C++ for dummies as I am ;-):
  // N         is an alias to the KN object.
  // N[0]      is a pointer to the first element of the vector.
  // N[0]+i    is a pointer to the ith element of the vector.
  // *(N[0]+i) is the value of the ith element of the vector.

  int nn = A->N( );    // get number of nodes

  cout << "nn: " << nn << endl;

  for (int i = 0; i < nn; i++) {
    (*(A[0] + i)) = (*(B[0] + i)) * (*(C[0] + i));
    cout << (*(A[0] + i)) << endl;
  }

  return 0.0;    // dummy return value.
}

double CppModTemplate4(KN< double > *const &A,                            // OUTPUT
                       KN< double > *const &B, KN< double > *const &C,    // INPUTS
                       KN< double > *const &D) {
  int nn = A->N( );    // get number of nodes

  cout << "nn: " << nn << endl;

  for (int i = 0; i < nn; i++) {
    (*(A[0] + i)) = (*(B[0] + i)) * (*(C[0] + i)) * (*(D[0] + i));
    cout << (*(A[0] + i)) << endl;
  }

  return 0.0;    // dummy return value.
}

double CppModTemplate5(KN< double > *const &A,                            // OUTPUT
                       KN< double > *const &B, KN< double > *const &C,    // INPUTS
                       KN< double > *const &D, KN< double > *const &E) {
  int nn = A->N( );    // get number of nodes

  cout << "nn: " << nn << endl;

  for (int i = 0; i < nn; i++) {
    (*(A[0] + i)) = (*(B[0] + i)) * (*(C[0] + i)) * (*(D[0] + i)) * (*(E[0] + i));
    cout << (*(A[0] + i)) << endl;
  }

  return 0.0;    // dummy return value.
}

double CppModTemplate6(KN< double > *const &A,                            // OUTPUT
                       KN< double > *const &B, KN< double > *const &C,    // INPUTS
                       KN< double > *const &D, KN< double > *const &E, KN< double > *const &F) {
  int nn = A->N( );    // get number of nodes

  cout << "nn: " << nn << endl;

  for (int i = 0; i < nn; i++) {
    (*(A[0] + i)) = (*(B[0] + i)) * (*(C[0] + i)) * (*(D[0] + i)) * (*(E[0] + i)) * (*(F[0] + i));
    cout << (*(A[0] + i)) << endl;
  }

  return 0.0;    // dummy return value.
}

double CppModTemplate7(KN< double > *const &A,                            // OUTPUT
                       KN< double > *const &B, KN< double > *const &C,    // INPUTS
                       KN< double > *const &D, KN< double > *const &E, KN< double > *const &F,
                       KN< double > *const &G) {
  int nn = A->N( );    // get number of nodes

  cout << "nn: " << nn << endl;

  for (int i = 0; i < nn; i++) {
    (*(A[0] + i)) =
      (*(B[0] + i)) * (*(C[0] + i)) * (*(D[0] + i)) * (*(E[0] + i)) * (*(F[0] + i)) * (*(G[0] + i));
    cout << (*(A[0] + i)) << endl;
  }

  return 0.0;    // dummy return value.
}

double CppModTemplate8(KN< double > *const &A,                            // OUTPUT
                       KN< double > *const &B, KN< double > *const &C,    // INPUTS
                       KN< double > *const &D, KN< double > *const &E, KN< double > *const &F,
                       KN< double > *const &G, KN< double > *const &H) {
  int nn = A->N( );    // get number of nodes

  cout << "nn: " << nn << endl;

  for (int i = 0; i < nn; i++) {
    (*(A[0] + i)) = (*(B[0] + i)) * (*(C[0] + i)) * (*(D[0] + i)) * (*(E[0] + i)) * (*(F[0] + i)) *
                    (*(G[0] + i)) * (*(H[0] + i));
    cout << (*(A[0] + i)) << endl;
  }

  return 0.0;    // dummy return value.
}

double funcs3(Stack s, const double &a, const double &b, const double &c) { return a + b + c; }

double funcs2(Stack s, const double &a, const double &b) { return a + b; }

double funcs1(Stack s, const double &a) { return a; }

// OneOperator1s_   s => add stack  to parameter of the function ...
// _ => pass const reference ...
typedef pair< pferbasearray, int > pferarray;

double mytest(Stack stack, pferarray const &p) {
  double ret = 0;

  typedef double K;
  int comp = p.second;
  cout << " " << comp << endl;    // Numero de la composation generalement 0
  pferbasearray pa = p.first;
  int N = pa->N;
  cout << " N = " << N << " ";

  // FEbaseArray
  for (int i = 0; i < N; ++i) {
    KN< K > *pui = pa->get(i);
    ret += (*pui)[0];
    KN< K > ui = *pui;    // copie du tableau
    (*pui)[0] += 1;
    pa->set(i, ui);
  }

  return ret;
}

// add the function name to the freefem++ table
static void Load_Init( ) {
  // Add function with 3 arguments
  Global.Add("funcs1", "(", new OneOperator1s_< double, double >(funcs1));
  Global.Add("funcs2", "(", new OneOperator2s_< double, double, double >(funcs2));
  Global.Add("funcs3", "(", new OneOperator3s_< double, double, double, double >(funcs3));
  Global.Add(
    "CppModTemplate3", "(",
    new OneOperator3_< double, KN< double > *, KN< double > *, KN< double > * >(CppModTemplate3));
  Global.Add(
    "CppModTemplate4", "(",
    new OneOperator4_< double, KN< double > *, KN< double > *, KN< double > *, KN< double > * >(
      CppModTemplate4));
  Global.Add("CppModTemplate5", "(",
             new OneOperator5_< double, KN< double > *, KN< double > *, KN< double > *,
                                KN< double > *, KN< double > * >(CppModTemplate5));
  Global.Add("CppModTemplate6", "(",
             new OneOperator6_< double, KN< double > *, KN< double > *, KN< double > *,
                                KN< double > *, KN< double > *, KN< double > * >(CppModTemplate6));
  Global.Add(
    "CppModTemplate7", "(",
    new OneOperator7_< double, KN< double > *, KN< double > *, KN< double > *, KN< double > *,
                       KN< double > *, KN< double > *, KN< double > * >(CppModTemplate7));
  Global.Add(
    "CppModTemplate8", "(",
    new OneOperator8_< double, KN< double > *, KN< double > *, KN< double > *, KN< double > *,
                       KN< double > *, KN< double > *, KN< double > *, KN< double > * >(
      CppModTemplate8));

  Global.Add("test", "(", new OneOperator1s_< double, pferarray >(mytest));
}

LOADFUNC(Load_Init)
