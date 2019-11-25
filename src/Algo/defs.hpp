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
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

#ifndef _DEFS_H_
#define _DEFS_H_

#ifndef GLIB
#define GLIB

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <stddef.h>
#include <cmath>
#include <cstring>
#include <list>
using namespace std;

#endif

// Define simple constants
#ifndef Pi
#define Pi (3.141592653589793)
#endif
#ifndef TINY    // used in BrentLS and Matrix
#define TINY 1.0e-20
#endif
#ifndef True
#define True (1)
#endif
#ifndef False
#define False (0)
#endif
#ifndef Yes
#define Yes (1)
#endif
#ifndef No
#define No (0)
#endif

// Define some simple template functions
#ifndef Sgn
#define Sgn(x) ((x) < 0 ? -1.0 : 1.0)
#endif
#ifndef Abs
#define Abs(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef Max
#define Max(x, y) ((x) > (y) ? (x) : (y))
#endif
#ifndef Min
#define Min(x, y) ((x) < (y) ? (x) : (y))
#endif

#ifndef inValidSize
#define inValidSize( ) cerr << "Warning: incompatible sizes!" << endl;
#endif

// Lists operations

/*!
 * \brief Display list
 * \param l List
 */
template< class Type >
void display(list< Type > l) {
  typename list< Type >::iterator il;
  for (il = l.begin( ); il != l.end( ); ++il) cout << (*il) << " ";
  cout << endl;
}

/*!
 * \brief Normalize list
 * \param l List
 */
template< class Type >
list< Type > normalize(list< Type > l) {
  list< Type > v(l);
  Type scale = l.front( );
  typename list< Type >::iterator il;

  if (scale == 0) {
    cerr << "First element is zero, cannot be normed! \n";
    return v;
  } else {
    for (il = v.begin( ); il != v.end( ); ++il) *il /= scale;
    return v;
  }
}

#endif    //_DEFS_H_
