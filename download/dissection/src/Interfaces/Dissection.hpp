/*! \file Dissection.hpp
    \brief Fortran style interface named as Dissectino-fortran interface
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Mar. 30th 2012
    \date   Jul. 12th 2015
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

#ifndef _INTERFACE_CPPAPI_HPP
# define _INTERFACE_CPPAPI_HPP

#define _COMPILER_OPTIONCOMPILER_H
#define FORTRAN_DECL_WL(x_windows,x_linux) x_linux
#define FORTRAN_DECL(x) x##_
#define DISSECTION_API

#ifdef BLAS_MKL
#include <mkl_service.h>
#endif
#include <stdint.h>
# include <cstdlib>

#define DISSECTION_REAL_MATRIX     1
#define DISSECTION_COMPLEX_MATRIX  2

#define DISS_VERSION FORTRAN_DECL_WL(DISS_VERSION, diss_version)
#define DISS_INIT FORTRAN_DECL_WL(DISS_INIT, diss_init)
#define DISS_FREE FORTRAN_DECL_WL(DISS_FREE, diss_free)
#define DISS_NUMERIC_FREE FORTRAN_DECL_WL(DISS_NUMERIC_FREE, diss_numeric_free)
#define DISS_S_FACT FORTRAN_DECL_WL(DISS_S_FACT, diss_s_fact)
#define DISS_N_FACT FORTRAN_DECL_WL(DISS_N_FACT, diss_n_fact)
#define DISS_GET_COLORS FORTRAN_DECL_WL(DISS_GET_COLORS, diss_get_colors)
#define DISS_GET_KERN_DIM FORTRAN_DECL_WL(DISS_GET_KERN_DIM, diss_get_kern_dim)
#define DISS_GET_NULLPIVOTS FORTRAN_DECL_WL(DISS_GET_NULLPIVOTS, diss_get_nullpivots)
#define DISS_GET_SMALLPIVOTS FORTRAN_DECL_WL(DISS_GET_SMALLPIVOTS, diss_get_smallpivots)
#define DISS_GET_KERN_VECS FORTRAN_DECL_WL(DISS_GET_KERN_VECS, diss_get_kern_vecs)
#define DISS_GET_KERNT_VECS FORTRAN_DECL_WL(DISS_GET_KERNT_VECS, diss_get_kernt_vecs)
#define DISS_PROJECT FORTRAN_DECL_WL(DISS_PROJECT, diss_project)
#define DISS_SOLVE_1 FORTRAN_DECL_WL(DISS_SOLVE_1, diss_solve_1)
#define DISS_SOLVE_N FORTRAN_DECL_WL(DISS_SOLVE_N, diss_solve_n)
#define DISS_MATRIX_PRODUCT FORTRAN_DECL_WL(DISS_MATRIX_PRODUCT, diss_matrix_product)
#define COMPUTE_DIM_KERN FORTRAN_DECL_WL(COMPUTE_DIM_KERN, compute_dim_kern)

extern "C" {
  DISSECTION_API void DISS_VERSION(int *versn,
				   int *reles,
				   int *patch);

  DISSECTION_API void DISS_INIT(uint64_t &dslv_,
				const int &called,
				const int &real_or_complex,
				const int &quad_fact,
				const int &nthreads,
				const int &verbose);

  DISSECTION_API void DISS_FREE(uint64_t &dslv_); 

  DISSECTION_API void DISS_NUMERIC_FREE(uint64_t &dslv_); 
  
  DISSECTION_API void DISS_S_FACT(uint64_t &dslv_,
				  const int &dim,
				  const int *ptRows,
				  const int *indCols,
				  const int &sym,
				  const int &decomposer);
  
  DISSECTION_API void DISS_N_FACT(uint64_t &dslv_,
				  const double *coefs,
				  const int &scaling,
				  const double &eps_pivot,
				  const int &indefinite_flag);

  DISSECTION_API void DISS_GET_COLORS(uint64_t &dslv_, int *n);
  DISSECTION_API void DISS_GET_KERN_DIM(uint64_t &dslv_, int *n0);

  DISSECTION_API void DISS_GET_NULLPIVOTS(uint64_t &dslv_, int *pivots);
  DISSECTION_API void DISS_GET_SMALLPIVOTS(uint64_t &dslv_,
					   const int &n,
					   int *pivots);

  DISSECTION_API void DISS_GET_KERN_VECS(uint64_t &dslv_, double *vec);
  DISSECTION_API void DISS_GET_KERNT_VECS(uint64_t &dslv_, double *vec);
  
  DISSECTION_API void DISS_PROJECT(uint64_t &dslv_, double *x);
  
  DISSECTION_API void DISS_SOLVE_1(uint64_t &dslv_, double *x,
				   const int &projection,
				   const int &trans);

  DISSECTION_API void DISS_SOLVE_N(uint64_t &dslv_, double *x,
				   const int &nrhs, const int &projection,
				   const int &trans);

  DISSECTION_API void DISS_MATRIX_PRODUCT(uint64_t &dslv_,
					  const double* x, double* y);

  DISSECTION_API void COMPUTE_DIM_KERN(int* flag,
				       int* n0,
				       double *a_ini,
				       const int &n,
				       const int &dim_ag,
				       const double &eps,
				       const double &machine_eps0,
				       const int &flag_sym,
				       const int *print_cntrl);
}
#endif
