/*
 * This file is part of FreeFem++.
 *
 * FreeFem++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FreeFem++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef INIT_FREEFEM_V3_
#define INIT_FREEFEM_V3_

#define EIGENVALUE
#ifndef VERBOSE
// #define VERBOSE
#endif

// bamg compilation FLAG
// ----------------------
// the freefem+ verson
// #define NDEBUG
// #define TEST 100
// pour bamg en version DEBUG
// #define DEBUG

#ifndef DRAWING
#define DRAWING
#endif

#ifndef BAMG_LONG_LONG
#ifndef __INTEL__
#define BAMG_LONG_LONG
#endif
#endif

// RNM  compilation FLAG
// ---------------------
#ifndef CHECK_KN
// #define CHECK_KN
#endif

// to use the umfpack linear solver library
#define UMFPACK
#define BLAS_UNDERSCORE

// virtual machine exec type checking  flag (very slow)
// pour faire un chech dynamique de tous les type du laguage
// --------------------------------
#ifndef WITHCHECK
// #define WITHCHECK
#endif
// to remove CHECKPTR checking
#ifndef NCHECKPTR
// #define NCHECKPTR
#endif
#endif

