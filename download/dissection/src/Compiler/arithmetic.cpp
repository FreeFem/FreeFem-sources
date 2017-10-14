/*! \file   arithmetic.cpp
    \brief  higher precision arithmetic for Kernel Detection
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jul. 17th 2015
    \date   Feb. 29th 2016
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

#include <cstdio>
#include "Compiler/arithmetic.hpp"
#include "Compiler/blas.hpp"

template<>
void printscalar<double>(const bool verbose, FILE *fp, double x)
{
  if (verbose && (fp != NULL)) {
    fprintf(fp, "%16.8e ", x);
  }
}

template<>
void printscalar<float>(const bool verbose, FILE *fp, float x)
{
  if (verbose && (fp != NULL)) {
    fprintf(fp, "%16.8e ", x);
  }
}

template<>
void printscalar<quadruple>(const bool verbose, FILE *fp, quadruple x)
{
  if (verbose && (fp != NULL)) {
    fprintf(fp, "%16.8e ", quad2double(x));
  }
}

template<>
void printscalar<complex<double> >(const bool verbose, FILE *fp,
				   complex<double> x)
{
  if (verbose && (fp != NULL)) {
    fprintf(fp, "(%16.8e %16.8e) ", x.real(), x.imag());
  }
}

template<>
void printscalar<complex<float> >(const bool verbose, FILE *fp,
				   complex<float> x)
{
  if (verbose && (fp != NULL)) {
    fprintf(fp, "(%16.8e %16.8e) ", x.real(), x.imag());
  }
}

template<>
void printscalar<complex<quadruple> >(const bool verbose, FILE *fp,
				      complex<quadruple> x)
{
  if (verbose && (fp != NULL)) {
    fprintf(fp, "(%16.8e %16.8e) ",
	    quad2double(x.real()), quad2double(x.imag()));
  }
}

template<typename T>
void printscalar(const bool verbose, FILE *fp, T x)
{
  fprintf(stderr, "%s %d : printscalar is not implented\n",
	  __FILE__, __LINE__);
}
#ifndef NO_OCTRUPLE
template
void printscalar<octruple>(const bool verbose, FILE *fp, octruple x);
template
void printscalar<complex<octruple> >(const bool verbose,
				     FILE *fp, complex<octruple> x);
#endif

