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
//ff-c++-LIBRARY-dep: umfpack amd blas
//ff-c++-cpp-dep:
/* clang-format on */

#include <iostream>
using namespace std;
using namespace std;
#include "rgraph.hpp"
#include "ff++.hpp"

#ifdef HAVE_LIBUMFPACK
extern "C" {
#ifdef HAVE_UMFPACK_H
#include <umfpack.h>
#else
#ifdef HAVE_UMFPACK_UMFPACK_H
#include <umfpack/umfpack.h>
#else
#ifdef HAVE_BIG_UMFPACK_UMFPACK_H
#include <UMFPACK/umfpack.h>
#else
#ifdef HAVE_UFSPARSE_UMFPACK_H
#include <ufsparse/umfpack.h>
#else
#ifdef HAVE_SUITESPARSE_UMFPACK_H
#include <suitesparse/umfpack.h>
#else

// Defaults to a local version of the UMFPACK headers
#include "../../3rdparty/include/umfpack.h"

#endif    // HAVE_SUITESPARSE_UMFPACK_H
#endif    // HAVE_UFSPARSE_UMFPACK_H
#endif    // HAVE_BIG_UMFPACK_UMFPACK_H
#endif    // HAVE_UMFPACK_UMFPACK_H
#endif    // HAVE_UMFPACK_H
}
#endif

#ifdef HAVE_LIBCHOLMOD_inprogress
#include <cholmod.h>
#endif

static void Load_Init( ) {}

LOADFUNC(Load_Init)
