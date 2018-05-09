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

// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

#ifndef THROWASSERT
#define THROWASSERT
#include <iostream>

// #ifdef __INTEL__
#define cerr cout
// #endif

#include "error.hpp"

#ifdef NDEBUG
#define throwassert(i) ((void)0)
#else
#define throwassert(condition) ((condition) ? ((void)0) : throw (ErrorAssert(#condition, __FILE__, __LINE__)))

#undef assert
#define assert(condition) throwassert(condition)
#endif

// <<ffassert>> an unremovable assert. According to FH, ffassert() is as a more reliable way to trap FF than assert().
#undef ffassert
#define ffassert(condition) ((condition) ? ((void)0) : throw (ErrorAssert(#condition, __FILE__, __LINE__)))
// #define AFAIRE(str) ( (cerr << " TO be Done " << str << endl), throw(ErrorAssert("AFAIRE)/TO DO  (FH????",__FILE__, __LINE__)))
#define AFAIRE(cmm) (cerr << "FH: A Faire/ To Do  " << cmm << " file " << __FILE__ << " line " << __LINE__ << endl, InternalError(cmm))

#define InternalError(message) throw (ErrorInternal(message, __LINE__, __FILE__))
#endif

