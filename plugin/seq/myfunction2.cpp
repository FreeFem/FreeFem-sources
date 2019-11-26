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

// Example C++ function "myfunction", dynamically loaded into "load.edp"

#include <ff++.hpp>
using namespace Fem2D;
double myf(string *s) {
  cout << *s << endl;
  return 0.;
}

double f(const double &x) { return x * x + 1; }

// Hack to do something at initialisation time
// to add the name myfunction to the freefem++ table
static void Load_Init( ) {
  Global.Add("Why", "(", new OneOperator1< double, string * >(myf));
  Global.Add("f", "(", new OneOperator1_< double, double >(f));
}

LOADFUNC(Load_Init)
