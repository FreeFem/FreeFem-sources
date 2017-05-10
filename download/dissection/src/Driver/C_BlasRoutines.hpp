/*! \file   C_BlasRoutines.hpp
    \brief  factorization routines LDL^t, LDU, forward/backward substitution
    \author François-Xavier Roux, ONERA, Laboratoire Jacques-Louis Lions
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
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

# ifndef _DRIVER_C_BLASROUTINES_
# define _DRIVER_C_BLASROUTINES_
#include <vector>
#include <list>
#include "Compiler/blas.hpp"
#include "Algebra/ColumnMatrix.hpp"

using std::vector;
using std::list;

template<typename T, typename U, typename Z>
Z matrix_infty_norm_(const int n, T *a, const int lda);

template<typename T, typename U>
U matrix_infty_norm(const int n, T *a, const int lda);
template<>
double matrix_infty_norm<quadruple, double>(const int n, quadruple *a,
					    const int lda);

template<>
double matrix_infty_norm<complex<quadruple>,
			 double>(const int n, complex<quadruple> *a,
				 const int lda);

#ifndef NO_OCTRUPLE
template<>
quadruple matrix_infty_norm<octruple, quadruple>(const int n,
						 octruple *a, const int lda);

template<>
quadruple matrix_infty_norm<complex<octruple>,
			    quadruple>(const int n,
				       complex<octruple> *a, const int lda);
#endif
template<typename T>
void full_ldlt(const int n, T *a, const int lda);

template<typename T>
void full_ldlh(const int n, T *a, const int lda);
template<>
void full_ldlh<double>(const int n, double *a, const int lda);
template<>
void full_ldlh<quadruple>(const int n, quadruple *a, const int lda);

template<typename T>
void full_ldu(const int n, T *a, const int lda);

template<typename T>
void FillUpperBlock(const int nrow, T *coef,
		    const int *prow,
		    const int *indcols, const int *indvals,
		    const int *old2new_i, const int *old2new_j,
		    ColumnMatrix<T> &b);

template<typename T, typename U>
bool full_ldlt_permute(int *n0, const int nn0, const int n, 
		       T *a, const int lda, 
		       double *pivot, int *permute, 
		       const double eps, double *fop);

template<typename T, typename U>
bool full_ldu_permute(int *n0, const int nn0, const int n, 
		      T *a, const int lda, 
		      double *pivot, int *permute, 
		      const double eps, double *fop);

template<typename T>
void full_fw_multiprofile(bool isTransposed, const int nrow, const int n0,
			  const int ncol,
			  T *a,
			  const int lda,
			  T *y, const int ldy,
			  vector<int> &i0, double *fop);

template<typename T>
void full_fw_single(bool isTransposed, const int nrow, const int n0,
		    T *a, const int lda,
		    T *x, 
		    double *fop);

template<typename T>
void full_bw_single(bool isTransposed, const int nrow, const int n0,
		   T *a, const int lda,
		   T *x, 
		   double *fop);

template<typename T>
void full_fw_multi(bool isTransposed, const int nrow, const int n0,
		   T *a, const int lda,
		   const int ncol,
		   T *x, const int ldy,
		   double *fop);

template<typename T>
void full_bw_multi(bool isTransposed, const int nrow, const int n0,
		   T *a, const int lda,
		   const int ncol,
		   T *x, const int ldy,
		   double *fop);

template<typename T>
void SparseSchur(const bool isSym, const int dim2, const int dim1,
		 vector<int>& i0,
		 ColumnMatrix<T> &upper,
		 ColumnMatrix<T> &lower,
		 ColumnMatrix<T> &diag,
		 double *fop);

template<typename T>
void full_fwbw_single(const bool isTrans, const int n, const int n0,
		      T *a, const int lda, T *z);

template<typename T>
void full_fwbw_multi(const bool isTrans, const int n, const int n0,
		     T *a, const int lda,
		     const int m, T *x, const int ldx);

template<typename T>
void full_fwbw_part(const int n, T *a, const int lda, T*x);

template<typename T>
void SchurProfileSym(const int nrow, const int ncol, vector<int> &i0,
		     ColumnMatrix<T> &b, ColumnMatrix<T> &c,
		     T* s, const int size_b1, double *fop);

template<typename T>
void SchurProfileUnSym(const int nrow, const int ncol, vector<int> &i0,
		       ColumnMatrix<T> &b, ColumnMatrix<T> &c,
		       T* s, const int size_b1, double *fop);

template<typename T>
void swap_sym_lower(const int n, T *a, const int lda, 
		    const int k, const int km, T *col_k, T *col_km);

template<typename T>
void swap_sym_upper(const int n, T *a, const int lda, 
		    const int k, const int km, T *col_k, T *col_km);

template<typename T>
void swap_unsym(const int n, T *a, const int lda, 
		const int k, const int km, T *col_k, T *col_km);

// routines for C_KernDetect.cpp
template<typename T, typename Z>
void full_fwbw_perturb_single(const int n,
			      T *a, const int lda, T *a_fact, T *x,
			      const int dim_augkern, const Z &eps,
			      bool flag_sym);

template<typename T, typename Z>
void full_fwbw_perturb_multi(const int n, const int m,
			     T *a, const int lda, T *a_fact, T *x,
			     const int dim_augkern, const Z &eps,
			     bool flag_sym);
//

template<typename T>
void full_sym_2x2BK(int n, T *a, T *dd1,
		    int *pivot_width, int *permute);

template<typename T>
void C_gemm_symm(const int ncol, const int nrow, const T &alpha, 
		 const T *a, const int lda, 
		 const T *b, const int ldb,
		 const T &beta,
		 T *c, const int ldc);

#endif
