/*! \file   DissectionDefault.hpp
    \brief  definition of default value for factorization
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jul. 24th 2015
    \date   Sep. 29th 2015
    \date   Nov. 30th 2016
*/

// This file is part of Dissection
// 
// Dissection is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Linking Dissection statically or dynamically with other modules is making
// a combined work based on Disssection. Thus, the terms and conditions of 
// the GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of Dissection give you 
// permission to combine Dissection program with free software programs or 
// libraries that are released under the GNU LGPL and with independent modules 
// that communicate with Dissection solely through the Dissection-fortran 
// interface. You may copy and distribute such a system following the terms of 
// the GNU GPL for Dissection and the licenses of the other code concerned, 
// provided that you include the source code of that other code when and as
// the GNU GPL requires distribution of source code and provided that you do 
// not modify the Dissection-fortran interface.
//
// Note that people who make modified versions of Dissection are not obligated 
// to grant this special exception for their modified versions; it is their
// choice whether to do so. The GNU General Public License gives permission to 
// release a modified version without this exception; this exception also makes
// it possible to release a modified version which carries forward this
// exception. If you modify the Dissection-fortran interface, this exception 
// does not apply to your modified version of Dissection, and you must remove 
// this exception when you distribute your modified version.
//
// This exception is an additional permission under section 7 of the GNU 
// General Public License, version 3 ("GPLv3")
//
// Dissection is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Dissection.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef _DISSECTION_DEFAULT_
#define _DISSECTION_DEFAULT_

#define SCOTCH_DECOMPOSER  0
#define METIS_DECOMPOSER   1
#define TRIDIAG_DECOMPOSER 2

#define NO_SCALING       0  // needs to be compatible to the definition in
#define DIAGONAL_SCALING 1  // SparseMatrix<T, U>::normalize(), SparseMatrix.cpp
#define KKT_SCALING      2  //

#define MINNODES     256    // minimum size of the first layer of dissection
#define SIZE_TRIDIAG 1000   // more than this value, dissection is used

#define DIM_AUG_KERN     4        // appropriate for indefinite matrix 
#define EPS_PIVOT        1.0e-2

#define TOL_PIVOT        1.0e-5   // for recursion of sparse factorization
#define MIN_TRIDIAG_SIZE 50       // the size to avoid Cuthill-McKee ordering

#endif
