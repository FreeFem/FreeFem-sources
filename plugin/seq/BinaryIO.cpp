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

#include <iostream>
#include <cfloat>
using namespace std;
#include "error.hpp"
#include "AFunction.hpp"
#include "rgraph.hpp"
#include "RNM.hpp"
#include "fem.hpp"
#include "FESpace.hpp"
#include "MeshPoint.hpp"
#include "AFunction_ext.hpp"    // Extension of "AFunction.hpp" to deal with more than 3 parameters function
using namespace Fem2D;

/*!
 * \brief Save vector
 * \param f KN<double> * const &
 * \param nome string * const &
 * \return 0.0
 */
double SaveVec(KN< double > *const &f, string *const &nome) {
  std::ofstream outfile(nome->data( ), ios_base::binary);
  // To access value at node i of vector N, do as follow: *(N[0]+i)
  // Explanation (C++ for dummies as I am ;-):
  // N         is an alias to the KN object.
  // N[0]      is a pointer to the first element of the vector.
  // N[0]+i    is a pointer to the ith element of the vector.
  // *(N[0]+i) is the value of the ith element of the vector.
  long int nn = f->N( );    // get number of nodes
  long int dim = nn;

  outfile.write((char *)&dim, sizeof(long int));    // write the dimension of the vector
  double ftemp;

  for (long int i = 0; i < nn; i++) {
    ftemp = *(f[0] + i);
    outfile.write((char *)&ftemp, sizeof(double));
  }

  outfile.close( );
  return 0.0;    // dummy return value.
}

/*!
 * \brief Write
 * \param io Stream_b<ostream> const &
 * \param data T * const &
 * \return ostream *
 */
double LoadVec(KN< double > *const &ww, string *const &nome) {
  std::ifstream infile(nome->data( ), ios_base::binary);
  long int dim;

  infile.read((char *)&dim, sizeof(long int));
  double dtemp;

  for (long int i = 0; i < dim; i++) {
    infile.read((char *)&dtemp, sizeof(double));
    *(ww[0] + i) = dtemp;
  }

  return 0.0;    // dummy return value.
}

double LoadFlag(long int *const &ww, string *const &nome) {
  std::ifstream infile(nome->data( ), ios_base::binary);
  long int flag;

  infile.read((char *)&flag, sizeof(long int));
  *ww = flag;
  return 0.0;    // dummy return value.
}

double flag(long int *const &FLAG, string *const &nome) {
  std::ofstream outfile(nome->data( ), ios_base::binary);
  long int Flag;

  Flag = *FLAG;
  outfile.write((char *)&Flag, sizeof(long int));
  outfile.close( );
  return 0.0;
}

// add the function name to the freefem++ table
static void Load_Init( ) {
  Global.Add("LoadVec", "(", new OneOperator2_< double, KN< double > *, string * >(LoadVec));
  Global.Add("LoadFlag", "(", new OneOperator2_< double, long int *, string * >(LoadFlag));
  Global.Add("SaveVec", "(", new OneOperator2_< double, KN< double > *, string * >(SaveVec));
  Global.Add("flag", "(", new OneOperator2_< double, long int *, string * >(flag));
}

LOADFUNC(Load_Init)
