/*! \file   C_BlasRoutines.cpp
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

#include "Driver/C_BlasRoutines.hpp"
#include "Driver/C_threads_tasks.hpp"
#include "Algebra/VectorArray.hpp"

// T may be complex of U that is higher precision than Z
template<typename T, typename U, typename Z>
Z matrix_infty_norm_(const int n, T *a, const int lda)
{
  Z err(0.0);
  for (int i = 0; i < n; i++) {
    U err_tmp1(0.0);
    for (int j = 0; j < n; j++) {
      err_tmp1 += blas_abs<T, U>(a[i + j * lda]);
    }
    Z err_tmp0 = tolower<U, Z>(err_tmp1); // accuracy conversion : U to Z
    err = err > err_tmp0 ? err : err_tmp0;
  }
  return err;
}

template
double matrix_infty_norm_<quadruple, quadruple, double>(const int n,
							quadruple *a,
							const int lda);
#ifndef NO_OCTRUPLE
template
quadruple matrix_infty_norm_<octruple, octruple, quadruple>(const int n,
							    octruple *a,
							    const int lda);
#endif
template
double matrix_infty_norm_<complex<quadruple>,
			  quadruple,
			  double>(const int n,
				  complex<quadruple> *a,
				  const int lda);
#ifndef NO_OCTRUPLE
template
quadruple matrix_infty_norm_<complex<octruple>,
			     octruple,
			     quadruple>(const int n,
					complex<octruple> *a, const int lda);
#endif
//

template<typename T, typename U>
U matrix_infty_norm(const int n, T *a, const int lda)
{
  fprintf(stderr, "%s %d : general template is not implemented\n",
	  __FILE__, __LINE__);
  return U(0.0);
}

template<>
double matrix_infty_norm<quadruple, double>(const int n,
					    quadruple *a, const int lda)
{
  return matrix_infty_norm_<quadruple, quadruple, double>(n, a, lda);
}
#ifndef NO_OCTRUPLE
template<>
quadruple matrix_infty_norm<octruple, quadruple>(const int n,
						 octruple *a, const int lda)
{
  return matrix_infty_norm_<octruple, octruple, quadruple>(n, a, lda);
}
#endif
template<>
double matrix_infty_norm<complex<quadruple>,
			 double>(const int n,
				 complex<quadruple> *a, const int lda)
{
  return matrix_infty_norm_<complex<quadruple>, quadruple, double>(n, a, lda);
}
#ifndef NO_OCTRUPLE
template<>
quadruple matrix_infty_norm<complex<octruple>,
			    quadruple>(const int n,
				       complex<octruple> *a, const int lda)
{
  return matrix_infty_norm_<complex<octruple>, octruple, quadruple>(n, a, lda);
}
#endif

// #define LDLT_LOWER

// *pivot is dealt as "double", which give no defect on accuracy, because
// magnitude is important in the value of pivot
template<typename T, typename U>
bool full_ldlt_permute(int *nn0, const int n0, const int n, 
		       T *a, const int lda,
		       double *pivot, int *permute, 
		       const double eps, double *fop)
{
  const int lda1 = lda + 1;
  T alpha;
  const T one(1.0);
  const T none(-1.0);
  bool flag = true;
  VectorArray<T> col_k(n); 
  VectorArray<T> col_km(n); 
  alpha = none;
  for (int i = 0; i < n; i++) {
    permute[i] = i;
  }
  const int n1 = n - n0;
  int k = 0;
  while (k < n1) {
    int km = k;
    double vmax = 0.0;
    km = blas_iamax<T, U>((n - k), a + (k * lda1), lda1) + k;
    vmax = blas_abs<T, double>(a[km * lda1]);
    if (vmax < (*pivot * eps)) { // enough in doule precision
      flag = false;
      break;
    }
    *pivot = vmax;
    if (km > k) {
      int itmp = permute[km];
      permute[km] = permute[k];
      permute[k] = itmp;
      // swap row/column
#ifdef LDLT_LOWER
      swap_sym_lower(n, a, lda, k, km, col_k.addrCoefs(), col_km.addrCoefs());
#else
      swap_sym_upper(n, a, lda, k, km, col_k.addrCoefs(), col_km.addrCoefs());
#endif
    }
    const T d = one / a[k * lda1];
    a[k * lda1] = d;
    if (k == (n1 - 1)) {
      k++;
      break;
    }
    alpha = (-d);
#ifdef LDLT_LOWER
    blas_syr<T>(CblasLower, (n - k - 1), alpha,
	     &a[k * lda1 + 1], 1,
	     &a[(k + 1) * lda1], lda);
#else
    blas_syr<T>(CblasUpper, (n - k - 1), alpha,
	     &a[k * lda1 + lda], lda,
	     &a[(k + 1) * lda1], lda);
#endif
    *fop += 3.0 * (double)(n - k - 1) * (double)(n - k) / 2.0;
#ifdef LDLT_LOWER
    blas_scal<T>((n - k - 1), d, &a[k * lda1 + 1], 1);
#else
    blas_scal<T>((n - k - 1), d, &a[k * lda1 + lda], lda);
#endif
    *fop += double(n - k - 1);
    k++;
  } // while (k < n1)
  // symmetrize the results : copy lower to upper
#ifdef LDLT_LOWER
  for (int j = 1; j < n1; j++) {
    for (int i = 0; i < j; i++) {
      a[i + j * lda] = a[j + i  * lda];
    }
  }
#else
  for (int j = 1; j < n1; j++) {
    for (int i = 0; i < j; i++) {
       a[j + i  * lda] = a[i + j * lda];
    }
  }
#endif
  *nn0 = n - k;
  return flag;
}

template
bool full_ldlt_permute<double, double>(int *nn0, const int n0,
				       const int n, 
				       double *a, const int lda, 
				       double *pivot, int *permute, 
				       const double eps, double *fop);

template
bool full_ldlt_permute<complex<double>, double>(int *nn0, const int n0,
						const int n, 
						complex<double> *a,
						const int lda, 
						double *pivot, int *permute, 
						const double eps, double *fop);

template
bool full_ldlt_permute<quadruple, quadruple>(int *nn0, const int n0,
					     const int n, 
					     quadruple *a, const int lda, 
					     double *pivot, int *permute, 
					     const double eps, double *fop);

template
bool full_ldlt_permute<complex<quadruple>, quadruple>(int *nn0, const int n0,
						      const int n, 
						      complex<quadruple> *a,
						      const int lda, 
						      double *pivot,
						      int *permute, 
						      const double eps,
						      double *fop);
//

template<typename T, typename U>
bool full_ldu_permute(int *nn0, const int n0,
		      const int n, T *a, const int lda, 
		      double *pivot, int *permute, 
		      const double eps, double *fop)
{
  const int lda1 = lda + 1;
  //  T *col_k, *col_km;
  const T one(1.0);
  const T none(-1.0);
  bool flag = true;
  VectorArray<T> col_k(n);
  VectorArray<T> col_km(n);
  for (int i = 0; i < n; i++) {
    permute[i] = i;
  }
  const int n1 = n - n0;
  int k = 0;
  while (k < n1) {
    int km = k;
    double vmax = 0.0;
    km = blas_iamax<T, U>((n - k), a + (k * lda1), lda1) + k;
    vmax = blas_abs<T, double>(a[km * lda1]);
    if (vmax < (*pivot * eps)) {
      flag = false;
      break;
    }
    *pivot = vmax;
    if (km > k) {
      int itmp = permute[km];
      permute[km] = permute[k];
      permute[k] = itmp;
      // swap row/column
      swap_unsym(n, a, lda, k, km, col_k.addrCoefs(), col_km.addrCoefs());
    }
    const T d = one / a[k * lda1];
    a[k * lda1] = d;
    if (k == (n1 - 1)) {
      k++;
      break;
    }
    //    alpha = none;
    blas_scal<T>((n - k - 1), d, &a[k * lda1 + lda], lda);
    blas_ger<T>((n - k - 1), (n - k - 1), none,
	     &a[k * lda1 + 1], 1, &a[k * lda1 + lda], lda,
	     &a[(k + 1) * lda1], lda);
    *fop += 3.0 * (double)(n - k - 1) * (double)(n - k) / 2.0;
    blas_scal<T>((n - k - 1), d, &a[k * lda1 + 1], 1);
    *fop += 2.0 * double(n - k - 1);
    k++;
  } // while (k < n1)

  *nn0 = n - k;
  return flag;
}

template
bool full_ldu_permute<double, double>(int *nn0, const int n0,
				      const int n, 
				      double *a, const int lda, 
				      double *pivot, int *permute, 
				      const double eps, double *fop);
template
bool full_ldu_permute<complex<double>, double>(int *nn0, const int n0,
					       const int n, 
					       complex<double> *a,
					       const int lda,
					       double *pivot, int *permute, 
					       const double eps, double *fop);
template
bool full_ldu_permute<quadruple, quadruple>(int *nn0, const int n0,
					    const int n, 
					    quadruple *a, const int lda, 
					    double *pivot, int *permute, 
					    const double eps, double *fop);
template
bool full_ldu_permute<complex<quadruple>, quadruple>(int *nn0, const int n0,
						     const int n, 
						     complex<quadruple> *a,
						     const int lda, 
						     double *pivot,
						     int *permute, 
						     const double eps,
						     double *fop);
//

template<typename T>
void swap_sym_lower(const int n, T *a, const int lda, 
		    const int k, const int km, T *col_k, T *col_km)
{ 
  const int lda1 = lda + 1;
  for (int i = k; i < n; i++) {
    col_k[i] = a[i + k * lda];     // save lower column k
  }
  for (int j = 0; j < k; j++) {
    col_k[j] = a[k + j * lda];     // save lower row k
  }
  for (int i = km; i < n; i++) {
    col_km[i] = a[i + km * lda];   // save lower column km
  }
  for (int j = 0; j < km; j++) {
    col_km[j] = a[km + j * lda];   // save lower row km
  }
  for (int i = (k + 1); i < n; i++) {
    a[i + k * lda] = col_km[i];    // restore lower column k
  }
  for (int i = (km + 1); i < n; i++) {
    a[i + km * lda] = col_k[i];    // restore lower column km
  }
  for (int j = 0; j < k; j++) {
    a[k + j * lda] = col_km[j];    // restore lower row k
  }
  for (int j = 0; j < km; j++) {
    a[km + j * lda] = col_k[j];    // restore lower row km
  }
  a[k * lda1] = col_km[km];
  a[km * lda1] = col_k[k];
  a[km + k * lda] = col_k[km];
}

template
void swap_sym_lower<double>(const int n, double *a, const int lda, 
			    const int k, const int km, 
			    double *col_k, double *col_km);

template
void swap_sym_lower<complex<double> >(const int n, complex<double> *a, 
				      const int lda, 
				      const int k, const int km, 
				      complex<double> *col_k, 
				      complex<double> *col_km);
template
void swap_sym_lower<quadruple>(const int n, quadruple *a, const int lda, 
			       const int k, const int km, 
			       quadruple *col_k, quadruple *col_km);

template
void swap_sym_lower<complex<quadruple> >(const int n, complex<quadruple> *a, 
					 const int lda, 
					 const int k, const int km, 
					 complex<quadruple> *col_k, 
					 complex<quadruple> *col_km);
//

template<typename T>
void swap_sym_upper(const int n, T *a, const int lda, 
		    const int k, const int km, T *col_k, T *col_km)
{ 
  const int lda1 = lda + 1;
  for (int i = k; i < n; i++) {
    col_k[i] = a[k + i * lda];     // save lower column k
  }
  for (int j = 0; j < k; j++) {
    col_k[j] = a[j + k * lda];     // save lower row k
  }
  for (int i = km; i < n; i++) {
    col_km[i] = a[km + i * lda];   // save lower column km
  }
  for (int j = 0; j < km; j++) {
    col_km[j] = a[j + km * lda];   // save lower row km
  }
  for (int i = (k + 1); i < n; i++) {
    a[k + i * lda] = col_km[i];    // restore lower column k
  }
  for (int i = (km + 1); i < n; i++) {
    a[km + i * lda] = col_k[i];    // restore lower column km
  }
  for (int j = 0; j < k; j++) {
    a[j + k * lda] = col_km[j];    // restore lower row k
  }
  for (int j = 0; j < km; j++) {
    a[j + km * lda] = col_k[j];    // restore lower row km
  }
  a[k * lda1] = col_km[km];
  a[km * lda1] = col_k[k];
  a[k + km * lda] = col_k[km];
}

template
void swap_sym_upper<double>(const int n, double *a, const int lda, 
			    const int k, const int km, 
			    double *col_k, double *col_km);

template
void swap_sym_upper<complex<double> >(const int n, complex<double> *a, 
				      const int lda, 
				      const int k, const int km, 
				      complex<double> *col_k, 
				      complex<double> *col_km);

template
void swap_sym_upper<quadruple>(const int n, quadruple *a, const int lda, 
			       const int k, const int km, 
			       quadruple *col_k, quadruple *col_km);

template
void swap_sym_upper<complex<quadruple> >(const int n, complex<quadruple> *a, 
					 const int lda, 
					 const int k, const int km, 
					 complex<quadruple> *col_k, 
					 complex<quadruple> *col_km);
//

template<typename T>
void swap_unsym(const int n, T *a, const int lda, 
		const int k, const int km, T *col_k, T *col_km)
{ 
  for (int i = 0; i < n; i++) {
    col_k[i] = a[i + k * lda];     // save lower column k
    col_km[i] = a[i + km * lda];   // save lower column km
  }
  for (int i = 0; i < n; i++) {
    a[i + km * lda] = col_k[i];   // save lower column km
    a[i + k * lda] = col_km[i];
  }
  for (int i = 0; i < n; i++) {
    col_k[i] = a[k + i * lda];     // save lower column k
    col_km[i] = a[km + i * lda];   // save lower column km
  }
  for (int i = 0; i < n; i++) {
    a[km + i * lda] = col_k[i];   // save lower column km
    a[k + i * lda] = col_km[i];
  }
}

template
void swap_unsym<double>(const int n, double *a, const int lda, 
			const int k, const int km, 
			double *col_k, double *col_km);

template
void swap_unsym<complex<double> >(const int n, complex<double> *a, 
				  const int lda, 
				  const int k, const int km, 
				  complex<double> *col_k, 
				  complex<double> *col_km);
template
void swap_unsym<quadruple>(const int n, quadruple *a, const int lda, 
			   const int k, const int km, 
			   quadruple *col_k, quadruple *col_km);

template
void swap_unsym<complex<quadruple> >(const int n, complex<quadruple> *a, 
				  const int lda, 
				  const int k, const int km, 
				  complex<quadruple> *col_k, 
				  complex<quadruple> *col_km);
//

template<typename T>
void full_ldlt(const int n, T *a, const int lda)
{
  T alpha;
  const T one(1.0);
  const int lda1 = lda + 1;
  for (int k = 0; k < (n - 1); k++) {
    T d = one / a[k * lda1];
    a[k * lda1] = d;
    alpha = (-d);
    blas_syr<T>(CblasLower, (n - k - 1), alpha,
	     &a[(k + 1) + k * lda], 1,
	     &a[(k + 1) * lda1], lda);
    blas_scal<T>((n - k - 1), d, &a[(k + 1) + k * lda], 1);
  }
  // symmetrize
  for (int j = 1; j < n; j++) {
    for (int i = 0; i < j; i++) {
      a[i + j * lda] = a[j + i * lda];
    }
  }
  a[(n - 1) * lda1] = one / a[(n - 1) * lda1];
#if 0
  d[n - 1] = one / a[(n - 1) * lda1];
  for (int i = 0; i < n; i++) {
    a[i * lda1] = d[i];
  }
#endif
}

template
void full_ldlt<double>(const int n, double *a, const int lda);

template
void full_ldlt<quadruple>(const int n, quadruple *a, const int lda);

#ifndef NO_OCTRUPLE
template
void full_ldlt<octruple>(const int n, octruple *a, const int lda);
#endif
template
void full_ldlt<complex<double> >(const int n,
				 complex<double> *a, const int lda);

template
void full_ldlt<complex<quadruple> >(const int n,
				    complex<quadruple> *a, const int lda);
#ifndef NO_OCTRUPLE
template
void full_ldlt<complex<octruple> >(const int n,
				   complex<octruple> *a, const int lda);
#endif
//

template<typename T>
void full_ldlh(const int n, T *a, const int lda)
{
  T alpha;
  const T one(1.0);
  const int lda1 = lda + 1;
  for (int k = 0; k < (n - 1); k++) {
    const T d = one / a[k * lda1];
    alpha = (-d);
    a[k * lda1] = d;
    blas_syrc<T>(CblasLower, (n - k - 1), alpha,
		 &a[(k + 1) + k * lda], 1,
		 &a[(k + 1) * lda1], lda);
    blas_scal<T>((n - k - 1), d, &a[(k + 1) + k * lda], 1);
  }
  // symmetrize
  for (int j = 1; j < n; j++) {
    for (int i = 0; i < j; i++) {
      a[i + j * lda] = std::conj(a[j + i * lda]);
    }
  }
  a[(n - 1) * lda1] = one / a[(n - 1) * lda1];
#if 0
  d[n - 1] = one / a[(n - 1) * lda1];
  for (int i = 0; i < n; i++) {
    a[i * lda1] = d[i];
  }
#endif
}

template<>
void full_ldlh<double>(const int n, double *a, const int lda)
{
  full_ldlt<double>(n, a, lda);
}

template<>
void full_ldlh<quadruple>(const int n, quadruple *a, const int lda)
{
  full_ldlt<quadruple>(n, a, lda);
}

#ifndef NO_OCTRUPLE
template<>
void full_ldlh<octruple>(const int n, octruple *a, const int lda)
{
  full_ldlt<octruple>(n, a, lda);
}
#endif
template
void full_ldlh<complex<double> >(const int n,
				 complex<double> *a, const int lda);

template
void full_ldlh<complex<quadruple> >(const int n,
				    complex<quadruple> *a, const int lda);
#ifndef NO_OCTRUPLE
template
void full_ldlh<complex<octruple> >(const int n,
				   complex<octruple> *a, const int lda);
#endif
//

template<typename T>
void full_ldu(const int n, T *a, const int lda)
{
  const int lda1 = lda + 1;
  const T one(1.0);
  const T none(-1.0);
  //  alpha = none;
  for (int k = 0; k < (n - 1); k++) {
    const T d = one / a[k * lda1];
    a[k * lda1] = d;
    blas_scal<T>((n - k - 1), d, &a[k * lda1 + lda], lda);
    blas_ger<T>((n - k - 1), (n - k - 1), none,
		&a[k * lda1 + 1], 1,
		&a[k * lda1 + lda], lda,
		&a[(k + 1) * lda1], lda);
    blas_scal<T>((n - k - 1), d, &a[k * lda1 + 1], 1);
  }
  a[(n - 1) * lda1] = one / a[(n - 1) * lda1];
  #if 0
  d[n - 1] = one / a[(n - 1) * lda1];
  for (int i = 0; i < n; i++) {
    a[i * lda1] = d[i];
  }
  #endif
}

template
void full_ldu<double>(const int n, double *a, const int lda);

template
void full_ldu<quadruple>(const int n, quadruple *a, const int lda);
#ifndef NO_OCTRUPLE
template
void full_ldu<octruple>(const int n, octruple *a, const int lda);
#endif
template
void full_ldu<complex<double> >(const int n,
				complex<double> *a, const int lda);

template
void full_ldu<complex<quadruple> >(const int n,
				   complex<quadruple> *a, const int lda);
#ifndef NO_OCTRUPLE
template
void full_ldu<complex<octruple> >(const int n,
				  complex<octruple> *a, const int lda); 
#endif
//

template<typename T>
void FillUpperBlock(const int nrow, T *coef,
		    const int *prow,
		    const int *indcols, const int *indvals,
		    const int *old2new_i, const int *old2new_j,
		    ColumnMatrix<T> &b)
{
  for (int i = 0; i < nrow; i++) {
    const int ii = old2new_i[i];
    for (int k = prow[ii]; k < prow[ii + 1]; k++) {
      const int jj = old2new_j[indcols[k]];
      b(ii, jj) = coef[indvals[k]];
    }
  }
}

template
void FillUpperBlock<double>(const int nnz, double *coef, const int *prow1,
			    const int *indcols1, const int *indvals1,
			    const int *new2old, const int *old_j,
			    ColumnMatrix<double> &b);
template
void FillUpperBlock<complex<double> >(const int nnz, complex<double> *coef,
				      const int *prow1,
				      const int *indcols1, const int *indvals1,
				      const int *new2old, const int *old_j,
				      ColumnMatrix<complex<double> > &b);
template
void FillUpperBlock<quadruple>(const int nnz, quadruple *coef, const int *prow1,
			    const int *indcols1, const int *indvals1,
			    const int *new2old, const int *old_j,
			    ColumnMatrix<quadruple> &b);
template
void FillUpperBlock<complex<quadruple> >(const int nnz,
					 complex<quadruple> *coef,
					 const int *prow1,
					 const int *indcols1,
					 const int *indvals1,
					 const int *new2old, const int *old_j,
					 ColumnMatrix<complex<quadruple> > &b);
//

template<typename T>
void full_fw_multiprofile(bool isTransposed, const int nrow, const int n0,
			  const int ncol,
			  T *a, const int lda,
			  T *y, const int ldy,
			   vector<int> &i0, double *fop)
{
  const T one(1.0);
  const T zero(0.0);
  int jlast;
  const int n1 = nrow - n0;
  list<int> j1;
  for (int j = 0; j < (ncol - 1); j++) {
    if (i0[j] < i0[j + 1]) {
      j1.push_back(j);
    }
  }
  j1.push_back(ncol - 1);
  jlast = (-1);
  for (list<int>::const_iterator jt = j1.begin(); jt != j1.end(); ++jt) {
    const int width = (*jt) - jlast;
    const int ioffset = i0[(*jt)];
    const int height = nrow - ioffset - n0;
    if (height > 0) {
      if (width > 1) {
	blas_trsm<T>(CblasLeft,
		     (isTransposed ? CblasUpper : CblasLower),
		     (isTransposed ? CblasTrans : CblasNoTrans),
		     CblasUnit, height, width,
		     one,
		     a + (ioffset * (nrow + 1)), nrow,
		     y + (ioffset + (jlast + 1) * ldy), ldy);
	*fop += (double)height * (double)(height - 1) * (double)width;
      }
      else {
	blas_trsv<T>((isTransposed ? CblasUpper : CblasLower),
		     (isTransposed ? CblasTrans : CblasNoTrans),
		     CblasUnit, height,
		     a + (ioffset * (nrow + 1)), nrow,
		     y + (ioffset + (jlast + 1) * ldy), 1);
	*fop += (double)height * (double)(height - 1);
      }
    }
    jlast = (*jt);
  } // loop : jt
  for (int j = 0; j < ncol; j++) {
    for (int i = n1; i < nrow; i++) {
      y[i + j * ldy] = zero;
    }
  }
}

template
void full_fw_multiprofile<double>(bool isTransposed, const int nrow,
				  const int n0,
				  const int ncol,
				  double *a, const int lda,
				  double *y, const int ldy,
				  vector<int> &i0, double *fop);

template
void full_fw_multiprofile<complex<double> >(bool isTransposed, const int nrow,
					    const int n0,
					    const int ncol,
					    complex<double> *a,
					    const int lda,
					    complex<double> *y, const int ldy,
					    vector<int> &i0, double *fop);

template
void full_fw_multiprofile<quadruple>(bool isTransposed, const int nrow,
				     const int n0,
				     const int ncol,
				     quadruple *a,
				     const int lda,
				     quadruple *y, const int ldy,
				     vector<int> &i0, double *fop);

template
void full_fw_multiprofile<complex<quadruple> >(bool isTransposed,
					       const int nrow,
					       const int n0,
					       const int ncol,
					       complex<quadruple> *a,
					       const int lda,
					       complex<quadruple> *y,
					       const int ldy,
					       vector<int> &i0, double *fop);
//

template<typename T>
void full_fw_multi(bool isTransposed, const int nrow, const int n0,
		   T *a, const int lda,
		   const int ncol,
		   T *y, const int ldy,
		   double *fop)
{
  const T one(1.0);
  const T zero(0.0);
  const int n1 = nrow - n0;
  // alpha = 1.0
  // symmetric matrix should have upper part by symmeterization : LDL^t
  blas_trsm<T>(CblasLeft,
	       (isTransposed ? CblasUpper : CblasLower),
	       (isTransposed ? CblasTrans : CblasNoTrans),
	       CblasUnit, n1, ncol,
	       one, a, lda, y, ldy);
  *fop += (double)nrow * (double)(nrow - 1) * (double)ncol;
  for (int j = 0; j < ncol; j++) {
    for (int i = n1; i < nrow; i++) {
      y[i + j * ldy] = zero;
    }
  }
}

template
void full_fw_multi<double>(bool isTransposed, const int nrow, const int n0,
			   double *a, const int lda,
			   const int ncol,
			   double *y, const int ldy,
			   double *fop);

template
void full_fw_multi<complex<double> >(bool isTransposed, const int nrow,
				     const int n0,
				     complex<double> *a, const int lda,
				     const int ncol,
				     complex<double> *y, const int ldy,
				     double *fop);

template
void full_fw_multi<quadruple>(bool isTransposed, const int nrow, const int n0,
			      quadruple *a, const int lda,
			      const int ncol,
			      quadruple *y, const int ldy,
			      double *fop);

template
void full_fw_multi<complex<quadruple> >(bool isTransposed, const int nrow,
					const int n0,
					complex<quadruple> *a, const int lda,
					const int ncol,
					complex<quadruple> *y, const int ldy,
					double *fop);

//

template<typename T>
void full_fw_single(bool isTransposed, const int nrow, const int n0,
		   T *a, const int lda,
		   T *y,
		   double *fop)
{
  const T zero(0.0);
  const int n1 = nrow - n0;
  blas_trsv<T>((isTransposed ? CblasUpper : CblasLower),
	       (isTransposed ? CblasTrans : CblasNoTrans),
	       CblasUnit, n1,
	       a, lda, y, 1);
  *fop += (double)nrow * (double)(nrow - 1);
  for (int i = n1; i < nrow; i++) {
    y[i] = zero;
  }
}

template
void full_fw_single<double>(bool isTransposed, const int nrow, const int n0,
			   double *a, const int lda,
			   double *y,
			   double *fop);

template
void full_fw_single<complex<double> >(bool isTransposed, const int nrow,
				     const int n0,
				     complex<double> *a, const int lda,
				     complex<double> *y,
				     double *fop);

template
void full_fw_single<quadruple>(bool isTransposed, const int nrow, const int n0,
			      quadruple *a, const int lda,
			      quadruple *y,
			      double *fop);

template
void full_fw_single<complex<quadruple> >(bool isTransposed, const int nrow,
					const int n0,
					complex<quadruple> *a, const int lda,
					complex<quadruple> *y,
					double *fop);
//

template<typename T>
void full_bw_single(bool isTransposed, const int nrow, const int n0,
		   T *a, const int lda,
		   T *y, 
		   double *fop)
{
  const T zero(0.0);
  const int n1 = nrow - n0;
  blas_trsv<T>((isTransposed ? CblasLower : CblasUpper),
	       (isTransposed ? CblasTrans : CblasNoTrans),
	       CblasUnit, n1,
	       a, lda, y, 1);
  *fop += (double)nrow * (double)(nrow - 1);
  for (int i = n1; i < nrow; i++) {
    y[i] = zero;
  }
}

template
void full_bw_single<double>(bool isTransposed, const int nrow, const int n0,
			   double *a, const int lda,
			   double *y,
			   double *fop);

template
void full_bw_single<complex<double> >(bool isTransposed, const int nrow,
				     const int n0,
				     complex<double> *a, const int lda,
				     complex<double> *y,
				     double *fop);

template
void full_bw_single<quadruple>(bool isTransposed, const int nrow, const int n0,
			      quadruple *a, const int lda,
			      quadruple *y, 
			      double *fop);

template
void full_bw_single<complex<quadruple> >(bool isTransposed, const int nrow,
					const int n0,
					complex<quadruple> *a, const int lda,
					complex<quadruple> *y, 
					double *fop);
//

template<typename T>
void full_bw_multi(bool isTransposed, const int nrow, const int n0,
		   T *a, const int lda,
		   const int ncol,
		   T *y, const int ldy,
		   double *fop)
{
  const T one(1.0);
  const T zero(0.0);
  const int n1 = nrow - n0;
  // alpha = 1.0
  // symmetric matrix should have upper part by symmeterization : LDL^t
  blas_trsm<T>(CblasLeft,
	       (isTransposed ? CblasLower : CblasUpper),
	       (isTransposed ? CblasTrans : CblasNoTrans),
	       CblasUnit, n1, ncol,
	       one, a, lda, y, ldy);
  *fop += (double)nrow * (double)(nrow - 1) * (double)ncol;
  for (int j = 0; j < ncol; j++) {
    for (int i = n1; i < nrow; i++) {
      y[i + j * ldy] = zero;
    }
  }
}

template
void full_bw_multi<double>(bool isTransposed, const int nrow, const int n0,
			   double *a, const int lda,
			   const int ncol,
			   double *y, const int ldy,
			   double *fop);

template
void full_bw_multi<complex<double> >(bool isTransposed, const int nrow,
				     const int n0,
				     complex<double> *a, const int lda,
				     const int ncol,
				     complex<double> *y, const int ldy,
				     double *fop);

template
void full_bw_multi<quadruple>(bool isTransposed, const int nrow, const int n0,
			      quadruple *a, const int lda,
			      const int ncol,
			      quadruple *y, const int ldy,
			      double *fop);

template
void full_bw_multi<complex<quadruple> >(bool isTransposed, const int nrow,
					const int n0,
					complex<quadruple> *a, const int lda,
					const int ncol,
					complex<quadruple> *y, const int ldy,
					double *fop);
//

template<typename T>
void SparseSchur(const bool isSym,
		 const int dim2, const int dim1,
		 vector<int>& i0,
		 ColumnMatrix<T> &upper,
		 ColumnMatrix<T> &lower,
		 ColumnMatrix<T> &diag,
		 double *fop)
{
  const T zero(0.0);
  const T none(-1.0);
  // alpha = -1, beta = 0
  if (isSym) {
    C_gemm_symm(dim2, dim1, none,
		lower.addrCoefs(), dim1,
		upper.addrCoefs(), dim1, zero,
		diag.addrCoefs(), dim2);
    // symmetrize : upper -> lower
    for (int i = 0; i < dim2; i++) {
      for (int j = 0; j < i; j++) {
	diag(i, j) = diag(j, i);
      }
    }
  }
  else {
    blas_gemm<T>(CblasTrans, CblasNoTrans, dim2, dim2, dim1, none,
		 lower.addrCoefs(), dim1,
		 upper.addrCoefs(), dim1, zero,
		 diag.addrCoefs(), dim2);
    *fop = (double)dim2 * (double)dim2 * (double)dim1;
  }
}

template
void SparseSchur<double>(const bool isSym,
			 const int dim2, const int dim1, vector<int>& i0,
			 ColumnMatrix<double> &upper,
			 ColumnMatrix<double> &lower,
			 ColumnMatrix<double> &diag,
			 double *fop);

template
void SparseSchur<complex<double> >(const bool isSym,
				   const int dim2, const int dim1,
				   vector<int>& i0,
				   ColumnMatrix<complex<double> > &upper,
				   ColumnMatrix<complex<double> > &lower,
				   ColumnMatrix<complex<double> > &diag,
				   double *fop);
template
void SparseSchur<quadruple>(const bool isSym,
			 const int dim2, const int dim1, vector<int>& i0,
			 ColumnMatrix<quadruple> &upper,
			 ColumnMatrix<quadruple> &lower,
			 ColumnMatrix<quadruple> &diag,
			 double *fop);

template
void SparseSchur<complex<quadruple> >(const bool isSym,
				   const int dim2, const int dim1,
				   vector<int>& i0,
				   ColumnMatrix<complex<quadruple> > &upper,
				   ColumnMatrix<complex<quadruple> > &lower,
				   ColumnMatrix<complex<quadruple> > &diag,
				   double *fop);
//

template<typename T>
void full_fwbw_single(const bool isTrans, const int n, const int n0,
		      T *a, const int lda, T *z)
{
  const T zero(0.0);
  const int n1 = n - n0;
  for (int i = n1; i < n; i++) {
    z[i] = zero;
  }
  blas_trsv<T>((isTrans ? CblasUpper : CblasLower),
	       (isTrans ? CblasTrans : CblasNoTrans), CblasUnit,
	       n1, a, lda, z, 1);
  for (int i = 0; i < n1; i++) {
    z[i] *= a[i * (n + 1)];
  }
  blas_trsv<T>((isTrans ? CblasLower : CblasUpper),
	       (isTrans ? CblasTrans : CblasNoTrans), CblasUnit,
	       n1, a, lda, z, 1);
}

template
void full_fwbw_single<double>(const bool isTrans, const int n, const int n0,
			      double *a, const int lda,
			      double *z);

template
void full_fwbw_single<complex<double> >(const bool isTrans, const int n,
					const int n0, complex<double> *a,
					const int lda,
					complex<double> *z);

template
void full_fwbw_single<quadruple>(const bool isTrans, const int n, const int n0,
				 quadruple *a, const int lda,
				 quadruple *z);

template
void full_fwbw_single<complex<quadruple> >(const bool isTrans, const int n,
					   const int n0, complex<quadruple> *a,
					   const int lda,
					   complex<quadruple> *z);
//

template<typename T>
void full_fwbw_multi(const bool isTrans, const int n, const int n0,
		     T *a, const int lda,
		     const int m, T *x, const int ldx)
{
  const T one(1.0);
  const T zero(0.0);
  const int n1 = n - n0;
  for (int j = 0; j < m; j++) {
    for (int i = n1; i < n; i++) {
      x[i + j * ldx] = zero;
    }
  }
  // alpha = 1.0
  blas_trsm<T>(CblasLeft,
	       (isTrans ? CblasUpper : CblasLower),
	       (isTrans ? CblasTrans : CblasNoTrans), CblasUnit,
	       n1, m, one, a, lda,
	       x, ldx);
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n1; i++) {
      //      x[i + j * ldx] *= a[i * (n + 1)]; // 3 Nov.2015
      x[i + j * ldx] *= a[i * (lda + 1)]; // 3 Nov.2015
    }
  }
  blas_trsm<T>(CblasLeft,
	       (isTrans ? CblasLower : CblasUpper),
	       (isTrans ? CblasTrans: CblasNoTrans), CblasUnit,
	       n1, m, one, a, lda,
	       x, ldx);
}

template
void full_fwbw_multi<double>(const bool isTrans,
			     const int n, const int n0,
			     double *a, const int lda,
			     const int m, double *x, const int ldx);

template
void full_fwbw_multi<complex<double> >(const bool isTrans,
				       const int n, const int n0,
				       complex<double> *a, const int lda,
				       const int m, complex<double> *x,
				       const int ldx);

template
void full_fwbw_multi<quadruple>(const bool isTrans,
				const int n, const int n0,
				quadruple *a, const int lda,
				const int m, quadruple *x, const int ldx);

template
void full_fwbw_multi<complex<quadruple> >(const bool isTrans,
					  const int n, const int n0,
					  complex<quadruple> *a, const int lda,
					  const int m, complex<quadruple> *x,
					  const int ldx);
//

template<typename T>
void full_fwbw_part(const int n, T *a, const int lda, T*x)
{
  int ifirst;
  const T zero(0.0);
  ifirst = n;
  for (int i = 0; i < n; i++) {
    if (x[i] != zero) {
      ifirst = i;
      break;
    }
  }
  if (ifirst == n) {
    return;
  }
  blas_trsv<T>(CblasLower, CblasNoTrans, CblasUnit, (n - ifirst),
	       &a[ifirst * (lda + 1)], lda,
	       &x[ifirst], 1);
  for (int i = ifirst; i < n; i++) {
    x[i] *= a[i * (lda + 1)];
  }
  blas_trsv<T>(CblasUpper, CblasNoTrans, CblasUnit, n,
	       a, lda,
	       x, 1);
}

template
void full_fwbw_part<double>(const int n,
			    double *a, const int lda, double *x);

template
void full_fwbw_part<quadruple>(const int n,
			       quadruple *a, const int lda, quadruple *x);
#ifndef NO_OCTRUPLE
template
void full_fwbw_part<octruple>(const int n,
			      octruple *a, const int lda, octruple *x);
#endif
template
void full_fwbw_part<complex<double> >(const int n, 
				      complex<double> *a, const int lda,
				      complex<double> *x);

template
void full_fwbw_part<complex<quadruple> >(const int n, 
					 complex<quadruple> *a,
					 const int lda, 
					 complex<quadruple> *x);
#ifndef NO_OCTRUPLE
template
void full_fwbw_part<complex<octruple> >(const int n, 
					complex<octruple> *a,
					const int lda, 
					complex<octruple> *x);
#endif
//

template<typename T>
void SchurProfileSym(const int nrow, const int ncol, vector<int> &i0,
		     ColumnMatrix<T> &b, ColumnMatrix<T> &c,
		     T *s, const int size_b1, double *fop)
{
  const T one(1.0);
  const T zero(0.0);
#if 1
  const int num_block = ncol / size_b1 + ((ncol % size_b1) == 0 ? 0 : 1);
  for (int j = 0; j < num_block; j++) {
    const int jj = j * size_b1;
    for (int i = 0; i < j; i++) {
      const int ii = i * size_b1;
      const int nrowb = std::min<int>(size_b1, ncol - ii);
      const int ncolb = std::min<int>(size_b1, ncol - jj);
      const int ioffset = std::max<int>(i0[ii], i0[jj]);
      //const int ioffset = 0;
      const int nnrow = nrow - ioffset;
      // alpha = 1
      // beta = 0
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   nrowb, ncolb,
		   nnrow,
		   one,
		   c.addrCoefs() + (ioffset + ii * nrow),
		   nrow,
		   b.addrCoefs() + (ioffset + jj * nrow),
		   nrow,
		   zero,
		   s + (ii + jj * ncol),
		   ncol);
      *fop += (double)nrowb * (double)ncolb * (double) nnrow;
    }
    {
      const int nrowb = std::min<int>(size_b1, ncol - jj);
            const int ioffset = i0[jj];
	    //const int ioffset = 0;
      const int nnrow = nrow - ioffset;
      // alpha = 1
      // beta = 0
      if (nnrow == 0) {
	// clear upper entries
	for (int j = 0; j < nrowb; j++) {
	  for (int i = 0; i <= j; i++) {
	    s[i + j * ncol + jj * (ncol + 1)] = zero;
	  }
	}
      }
      else {
	C_gemm_symm<T>(nrowb,
		       nnrow,
		       one,
		       c.addrCoefs() + (ioffset + jj * nrow),
		       nrow,
		       b.addrCoefs() + (ioffset + jj * nrow),
		       nrow,
		       zero,
		       s + (jj * (ncol + 1)),
		       ncol);
      }
      *fop += (double)nrowb * (double)nrowb * (double) nnrow / 2.0;
    }
  } // loop : j
#else // debug : compute all entries
      // alpha = 1
      // beta = 0
  //  blas_gemm<T>(CblasTrans, CblasNoTrans, 
  //	       ncol, ncol,
  C_gemm_symm<T>(ncol,
	       nrow,
	       one,
	       c.addrCoefs(),
	       nrow,
	       b.addrCoefs(),
	       nrow,
	       zero,
	       s,
	       ncol);
	*fop += (double)ncol * (double)ncol * (double) nrow;
#endif
}

template
void SchurProfileSym<double>(const int nrow, const int ncol,
			     vector<int> &i0,
			     ColumnMatrix<double> &b,
			     ColumnMatrix<double> &c,
			     double *s,
			     const int size_b1,
			     double *fop);

template
void SchurProfileSym<complex<double> >(const int nrow, const int ncol,
				       vector<int> &i0,
				       ColumnMatrix<complex<double> > &b,
				       ColumnMatrix<complex<double> > &c,
				       complex<double> *s,
				       const int size_b1,
				       double *fop);
template
void SchurProfileSym<quadruple>(const int nrow, const int ncol,
				vector<int> &i0,
				ColumnMatrix<quadruple> &b,
				ColumnMatrix<quadruple> &c,
				quadruple *s,
				const int size_b1,
				double *fop);

template
void SchurProfileSym<complex<quadruple> >(const int nrow, const int ncol,
					  vector<int> &i0,
					  ColumnMatrix<complex<quadruple> > &b,
					  ColumnMatrix<complex<quadruple> > &c,
					  complex<quadruple> *s,
					  const int size_b1,
					  double *fop);
//

template<typename T>
void SchurProfileUnSym(const int nrow, const int ncol, vector<int> &i0,
		       ColumnMatrix<T> &b, ColumnMatrix<T> &c,
		       T* s, const int size_b1, double *fop)
{
  const T one(1.0);
  const T zero(0.0);
  const int num_block = ncol / size_b1 + (ncol % size_b1 == 0 ? 0 : 1);
  for (int j = 0; j < num_block; j++) {
    const int jj = j * size_b1;
    for (int i = 0; i < num_block; i++) {
      const int ii = i * size_b1;
      const int nrowb = std::min<int>(size_b1, ncol - ii);
      const int ncolb = std::min<int>(size_b1, ncol - jj);
      const int ioffset = std::max<int>(i0[ii], i0[jj]);
      const int nnrow = nrow - ioffset;
      // alpha=1
      // beta=0
      blas_gemm<T>(CblasTrans, CblasNoTrans, 
		   nrowb, ncolb,
		   nnrow,
		   one,
		   c.addrCoefs() + (ioffset + ii * nrow),
		   nrow,
		   b.addrCoefs() + (ioffset + jj * nrow),
		   nrow,
		   zero,
		   s + (ii + jj * ncol),
		   ncol);
      *fop += (double)nrowb * (double)ncolb * (double) nnrow;
    }
  }
}

template
void SchurProfileUnSym<double>(const int nrow, const int ncol,
			       vector<int> &i0,
			       ColumnMatrix<double> &b,
			       ColumnMatrix<double> &c,
			       double *s, const int size_b1, double *fop);

template
void SchurProfileUnSym<complex<double> >(const int nrow, const int ncol,
					 vector<int> &i0,
					 ColumnMatrix<complex<double> > &b,
					 ColumnMatrix<complex<double> > &c,
					 complex<double> *s, const int size_b1,
					 double *fop);
template
void SchurProfileUnSym<quadruple>(const int nrow, const int ncol,
				  vector<int> &i0,
				  ColumnMatrix<quadruple> &b,
				  ColumnMatrix<quadruple> &c,
				  quadruple *s, const int size_b1, double *fop);

template
void
SchurProfileUnSym<complex<quadruple> >(const int nrow, const int ncol,
				       vector<int> &i0,
				       ColumnMatrix<complex<quadruple> > &b,
				       ColumnMatrix<complex<quadruple> > &c,
				       complex<quadruple> *s, const int size_b1,
				       double *fop);
//

// Z for perturbation has lower precision than T (real or complex valued)
template<typename T, typename Z>
void full_fwbw_perturb_single(const int n,
			      T *a, const int lda, T *a_fact, T *x,
			      const int dim_augkern, const Z &eps,
			      bool flag_sym)
{
  //  T alpha, beta;
  //  T *upper, *lower, *schur;
  const T one(1.0);
  const T none(-1.0);
  const T Teps = tohigher<T, Z>(eps);
  int n1;
  if (dim_augkern < n) {
    n1 = n - dim_augkern;
  }
  else {
    n1 = 0;
  }
  if (n1 > 0) {
    ColumnMatrix<T> upper(dim_augkern, n1); //  new T[dim_augkern * n1];
    ColumnMatrix<T> lower(dim_augkern, n1);  // lower = new T[dim_augkern * n1];
    ColumnMatrix<T> schur(n1, n1); // = new T[n1 * n1];

    // factorization with perturbation should be shared with both sing and multi
    for (int j = 0; j < n1; j++) {
      for (int i = 0; i < dim_augkern; i++) {
	upper(i, j) = a[i + (j + dim_augkern) * lda];
      }
      for (int i = 0; i < n1; i++) {
	schur(i, j) = a[i + dim_augkern + (j + dim_augkern) * lda];
      }
    }
    {
      //      alpha = one;
      blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		   dim_augkern, n1, one, a_fact, lda,
		   upper.addrCoefs(), dim_augkern);
      for (int i = 0; i < dim_augkern; i++) {
	for (int j = 0; j < n1; j++) {
	  upper(i, j) *= a_fact[i * (lda + 1)];
	}
      }
      // alphe = one;
      blas_trsm<T>(CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
		dim_augkern, n1, one, a_fact, lda,
		   upper.addrCoefs(), dim_augkern);
      // emulating machine eplison in double
      for (int j = 0; j < n1; j++) {
	upper((dim_augkern - 1), j) += Teps;
      }
    }
    if (!flag_sym) {
      for (int j = 0; j < n1; j++) {
	for (int i = 0; i < dim_augkern; i++) {
	  lower(i, j) = a[(j + dim_augkern) + i * lda];
	}
      }
      // alpha = one;
      blas_trsm<T>(CblasLeft, CblasUpper, CblasTrans, CblasUnit,
		   dim_augkern, n1, one, a_fact, lda,
		   lower.addrCoefs(), dim_augkern);
      for (int i = 0; i < dim_augkern; i++) {
	for (int j = 0; j < n1; j++) {
	  lower(i, j) *= a_fact[i * (lda + 1)];
	}
      }
      // alpha = one;
      blas_trsm<T>(CblasLeft, CblasLower, CblasTrans, CblasUnit,
		   dim_augkern, n1, one, a_fact, lda,
		   lower.addrCoefs(), dim_augkern);
      for (int j = 0; j < n1; j++) {
	lower((dim_augkern - 1), j) += Teps;
      }
    }
    //   alpha = none;
    //   beta = one;
    blas_gemm<T>(CblasNoTrans, CblasNoTrans, n1, n1, dim_augkern,
		 none, &a[dim_augkern], lda,
		 upper.addrCoefs(), dim_augkern,
		 one, schur.addrCoefs(), n1);
    if (flag_sym) {
      full_ldlt<T>(n1, schur.addrCoefs(), n1);
    }
    else {
      full_ldu<T>(n1, schur.addrCoefs(), n1);
    }
    // alpha = none;
    // beta = one;
    if (flag_sym) {
      blas_gemv<T>(CblasTrans, dim_augkern, n1, none,
		   upper.addrCoefs(), dim_augkern,
		   x, 1, one, &x[dim_augkern], 1);
    }
    else {
      blas_gemv<T>(CblasTrans, dim_augkern, n1, none,
		   lower.addrCoefs(), dim_augkern,
		   x, 1, one, &x[dim_augkern], 1);
    }
    // forward
    blas_trsv<T>(CblasLower, CblasNoTrans, CblasUnit,
		 dim_augkern, a_fact, lda, x, 1);
    for (int j = 0; j < dim_augkern; j++) {
      x[j] *= a_fact[j * (lda + 1)];
    }
    // backward
    blas_trsv<T>(CblasUpper, CblasNoTrans, CblasUnit,
		 dim_augkern, a_fact, lda, x, 1);
    // forward
    blas_trsv<T>(CblasLower, CblasNoTrans, CblasUnit,
		 dim_augkern, schur.addrCoefs(), n1, &x[dim_augkern], 1);
    // diagonal divide
    for (int j = 0; j < n1; j++) {
      x[j + dim_augkern] *= schur(j, j); //[j * (n1 + 1)];
    }
    // backward
    blas_trsv<T>(CblasUpper, CblasNoTrans, CblasUnit,
		 dim_augkern, schur.addrCoefs(), n1, &x[dim_augkern], 1);
    // alpha = none;
    // beta = one;
    blas_gemv<T>(CblasNoTrans, dim_augkern, n1, none,
		 upper.addrCoefs(), dim_augkern,
		 &x[dim_augkern], 1, one, x, 1);
  }
  else {
// forward
    blas_trsv<T>(CblasLower, CblasNoTrans, CblasUnit,
		 n, a_fact, lda, x, 1);
    for (int j = 0; j < n; j++) {
      x[j] *= a_fact[j * (lda + 1)];
    }
    x[(n - 1)] += Teps;
// backward
    blas_trsv<T>(CblasUpper, CblasNoTrans, CblasUnit,
		 n, a_fact, lda, x, 1);
  }
}

template
void full_fwbw_perturb_single<quadruple, double>(const int n,
						 quadruple *a, const int lda, 
						 quadruple *a_fact,
						 quadruple *x,
						 const int dim_augkern,
						 const double &eps,
						 bool flag_sym);
#ifndef NO_OCTRUPLE
template
void full_fwbw_perturb_single<octruple, quadruple>(const int n,
						   octruple *a, const int lda, 
						   octruple *a_fact,
						   octruple *x,
						   const int dim_augkern,
						   const quadruple &eps,
						   bool flag_sym);
#endif
template
void full_fwbw_perturb_single<complex<quadruple>,
			      double>(const int n,
				      complex<quadruple> *a, const int lda,
				      complex<quadruple> *a_fact,
				      complex<quadruple> *x,
				      const int dim_augkern,
				      const double &eps,
				      bool flag_sym);
#ifndef NO_OCTRUPLE
template
void full_fwbw_perturb_single<complex<octruple>,
			      quadruple>(const int n,
					 complex<octruple> *a, const int lda, 
					 complex<octruple> *a_fact,
					 complex<octruple> *x,
					 const int dim_augkern,
					 const quadruple &eps,
					 bool flag_sym);
#endif
//
// Z for perturbation has lower precision than T (real or complex valued)
template<typename T, typename Z>
void full_fwbw_perturb_multi(const int n, const int m,
			     T *a, const int lda, T *a_fact, T *x,
			     const int dim_augkern, const Z &eps, bool flag_sym)
{
  //  T alpha, beta;
  //  T *upper, *lower, *schur, *v;
  const T one(1.0);
  const T none(-1.0);
  const T Teps = tohigher<T, Z>(eps);
  if (dim_augkern >= n) {
// forward
//    alpha = one;
    blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		 n, m, one, a_fact, lda, x, lda);
    for (int j = 0; j < m; j++) {
      for (int i= 0; i < n; i++) {
	x[i + j * lda] *= a_fact[i * (lda + 1)];
      }
      x[(n - 1) + j * lda] += Teps;
    }
    // backward
    blas_trsm<T>(CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
		 n, m, one, a_fact, lda, x, lda);
  }
  else {
    const int n1 = n - dim_augkern;
    ColumnMatrix<T> upper(dim_augkern, n1);
    ColumnMatrix<T> lower(dim_augkern, n1);
    ColumnMatrix<T> schur(n1, n1);
    // factorization with perturbation should be shared with both sing and multi
    for (int j = 0; j < n1; j++) {
      for (int i = 0; i < dim_augkern; i++) {
	upper(i, j) = a[i + (j + dim_augkern) * lda];
      }
      for (int i = 0; i < n1; i++) {
	schur(i, j) = a[i + dim_augkern + (j + dim_augkern) * lda];
      }
    }
    {
      //      alpha = one;
      blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		   dim_augkern, n1, one, a_fact, lda,
		   upper.addrCoefs(), dim_augkern);
      for (int i = 0; i < dim_augkern; i++) {
	for (int j = 0; j < n1; j++) {
	  upper(i, j) *= a_fact[i * (lda + 1)];
	}
      }
      //      alpha = one;
      blas_trsm<T>(CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
		   dim_augkern, n1, one, a_fact, lda,
		   upper.addrCoefs(), dim_augkern);
      for (int j = 0; j < n1; j++) {
	upper((dim_augkern - 1), j) += Teps;
      }
    }
    if (!flag_sym) {
      for (int j = 0; j < n1; j++) {
	for (int i = 0; i < dim_augkern; i++) {
	  lower(i, j) = a[(j + dim_augkern) + i * lda];
	}
      }
      //      alpha = one;
      blas_trsm<T>(CblasLeft, CblasUpper, CblasTrans, CblasUnit,
		   dim_augkern, n1, one, a_fact, lda,
		   lower.addrCoefs(), dim_augkern);
      for (int i = 0; i < dim_augkern; i++) {
	for (int j = 0; j < n1; j++) {
	  lower(i, j) *= a_fact[i * (lda + 1)];
	}
      }
      //      alpha = one;
      blas_trsm<T>(CblasLeft, CblasLower, CblasTrans, CblasUnit,
		   dim_augkern, n1, one, a_fact, lda,
		   lower.addrCoefs(), dim_augkern);
      for (int j = 0; j < n1; j++) {
	lower((dim_augkern - 1), j) += Teps;
      }
    }
    //    alpha = none;
    //    beta = one;
    blas_gemm<T>(CblasNoTrans, CblasNoTrans, n1, n1, dim_augkern,
		 none, &a[dim_augkern], lda,
		 upper.addrCoefs(), dim_augkern,
		 one, schur.addrCoefs(), n1);
    if (flag_sym) {
      full_ldlt<T>(n1, schur.addrCoefs(), n1);
    }
    else {
      full_ldu<T>(n1, schur.addrCoefs(), n1);
    }
    //    alpha = none;
    //    beta = one;
    if (flag_sym) {
      blas_gemm<T>(CblasTrans, CblasNoTrans, n1, m, dim_augkern, none,
		   upper.addrCoefs(), dim_augkern,
		   x, lda, one, &x[dim_augkern], lda);
    }
    else {
      blas_gemm<T>(CblasTrans, CblasNoTrans, n1, m, dim_augkern, none,
		   lower.addrCoefs(), dim_augkern,
		   x, lda, one, &x[dim_augkern], lda);
    }
    //    alpha = one;
    // forward
    blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		 dim_augkern, m, one, a_fact, lda, x, lda);
    // diagonal divide
    for (int j = 0; j < m; j++) {
      for (int i = 0; i < dim_augkern; i++) {
	x[i + j * lda] *= a_fact[i * (lda + 1)];
      }
    }
    // backward
    blas_trsm<T>(CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
		 dim_augkern, m, one, a_fact, lda, x, lda);
    // forward
    blas_trsm<T>(CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
		 n1, m, one, schur.addrCoefs(), n1, &x[dim_augkern], lda);
    for (int j = 0; j < m; j++) {
      for (int i = 0; i < n1; i++) {
	x[i + dim_augkern + j * lda] *= schur(i, i);//
      }
    }
    // backward
    blas_trsm<T>(CblasLeft, CblasUpper, CblasNoTrans, CblasUnit,
		 n1, m, one, schur.addrCoefs(), n1, &x[dim_augkern], lda);
    //    alpha = none;
    //    beta = one;
    blas_gemm<T>(CblasNoTrans, CblasNoTrans, dim_augkern, m, n1, none,
		 upper.addrCoefs(), dim_augkern,
		 &x[dim_augkern], lda, one, x, lda);
  }
}

template
void full_fwbw_perturb_multi<quadruple, double>(const int n,
						const int m,
						quadruple *a, const int lda, 
						quadruple *a_fact,
						quadruple *x,
						const int dim_augkern,
						const double &eps,
						bool flag_sym);
#ifndef NO_OCTRUPLE
template
void full_fwbw_perturb_multi<octruple, quadruple>(const int n,
						  const int m, 
						  octruple *a, const int lda, 
						  octruple *a_fact,
						  octruple *x,
						  const int dim_augkern,
						  const quadruple &eps,
						  bool flag_sym);
#endif
template
void full_fwbw_perturb_multi<complex<quadruple>,
			     double>(const int n,
				     const int m,
				     complex<quadruple> *a, const int lda, 
				     complex<quadruple> *a_fact,
				     complex<quadruple> *x,
				     const int dim_augkern,
				     const double &eps,
				     bool flag_sym);
#ifndef NO_OCTRUPLE
template
void full_fwbw_perturb_multi<complex<octruple>,
			     quadruple>(const int n,
					const int m,
					complex<octruple> *a, const int lda, 
					complex<octruple> *a_fact,
					complex<octruple> *x,
					const int dim_augkern,
					const quadruple &eps,
					bool flag_sym);
#endif

#if 0
template<>
void full_sym_2x2BK<double>(int n, double *a, double *dd1,
			    int *pivot_width, int *permute)
{
  double *v, *col_k, *col_km;
  double *dd0;
  double *aa;
  const double bk_const = (1.0 + sqrt(17.0)) / 8.0;
  const int lda = n;
  const int lda1 = n + 1;
  int k, kk, km, kn;
  double colmax, rowmax, deta;
  double alpha, beta;
  
  v = new double[n];
  col_k = new double[n];
  col_km = new double[n];
  dd0 = new double[n];

  for (int i = 0; i < n; i++) {
    dd1[i] = 0.0;
  }
  for (int i = 0; i < n; i++) {
    permute[i] = i;
  }

  k = 0;
  while (k < n) {
    bool pivot_flag = false;
    if (k == (n - 1)) {
      pivot_width[k] = 1;
    }
    else {
      {
        double vmax = 0.0;
        for (int i = k + 1 ; i < n; i++) {
	  double xtmp = fabs(a[i + k * n]);
          if (vmax < xtmp) {
            km = i;
            vmax = xtmp;
          }
        }
        colmax = vmax;
      }
      if (a[k * lda1] >= bk_const * colmax) {
        pivot_width[k] = 1;
      }
      else {
	{
	  double vmax = 0.0;
	  for (int i = k + 1 ; i < n; i++) {
	    double xtmp = fabs(a[k + i * n]);
	    if (vmax < xtmp) {
	      kn = i;
	      vmax = xtmp;
	    }
	  }
	  rowmax = vmax;
	}
	if (rowmax * fabs(a[k * lda1]) >= (bk_const * colmax * colmax)) {
	  pivot_width[k] = 1;
	}
	else if (fabs(a[km * lda1]) >= (bk_const * rowmax)) {
	  pivot_width[k] = 1;
	  kk = permute[km];
	  permute[km] = permute[k];
	  permute[k] = kk;
	  swap_sym_lower<double>(n, a, lda, k, km, col_k, col_km);
	}
	else {
	  pivot_width[k] = 20;
	  pivot_width[k + 1] = 21;
	  pivot_flag = true;
	  swap_sym_lower<double>(n, a, lda, (k+ 1), km, col_k, col_km);
	}
      }
    }  // if (k == (n - 1)) {
    if (pivot_flag) {
      dd0[k] = 1.0 / a[k * lda1];
      a[k * lda1] = dd0[k];
      if (k == (n - 1)) {
	break;
      }
      alpha = (-dd0[k]);
      blas_syr<double>(CblasLower, (n - k - 1), alpha, &a[k + 1 + k * lda], 1,
	       &a[(k + 1) * lda1], n);
      blas_scal<double>((n - k  -1), dd0[k], &a[k + 1 + k * lda], 1);
      k++;
    }
    else { // if (pivot_flag)
      // computation of inverse of 2x2 matrix
      deta = (a[k * lda1] * a[(k + 1) * lda1] -
	      a[(k + 1) + k * lda] * a[k + (k + 1) * lda]);
      deta = 1.0 / deta;
      dd0[k] = a[(k + 1) * lda1] * deta;
      dd0[k + 1] = a[k * lda1] * deta;
      dd1[k] = (-a[(k + 1) + k * lda]) * deta;
      dd1[k + 1] = dd1[k];
      a[k * lda1] = dd0[k];
      a[(k + 1) * lda1] = dd1[k];
      a[(k + 1) * lda1 + k] = 0.0;
      if (k == (n - 2)) {
	break;
      }
      aa[0] = dd0[k];     // (0 ,0) 
      aa[2] = dd1[k];     // (0, 1)
      aa[1] = dd1[k];     // (1, 0)
      aa[3] = dd0[k + 1]; // (1, 1)
      //
      alpha = (-aa[0]);
      blas_syr<double>(CblasLower, (n - k - 2), alpha, &a[k + 2 + k * lda], 1,
	       &a[(k + 2) * lda1], n);
      //
      alpha = (-aa[3]);
      blas_syr<double>(CblasLower, (n - k - 2), alpha, &a[k + 2 + (k + 1) * lda], 1,
	       &a[(k + 2) * lda1], n);
      //
      alpha = (-aa[2]);
      blas_syr2<double>(CblasLower, (n - k - 2), alpha, &a[k + 2 + k * lda], 1,
		&a[k + 2 + (k + 1) * lda], 1,
		&a[(k + 2) * lda1], n);
      // copy lower to upper to save working palce
      for (int i = k + 2; i < n; i++) {
	for (int j = k; j <= k + 2; j++) { // three lines?
	  a[j + i * lda] = a[i + j * lda];
	}
      }
      alpha = 1.0;
      beta = 0.0;
      blas_gemm<double>(CblasTrans, CblasNoTrans, (n - k - 2), 2, 2, alpha,
		&a[k + (k + 2) * lda], n, aa, 2, beta,
		&a[(k + 2) + k * lda], n);
      k += 2;
    } // if (pivot_flag)
  } // while (k < n)
  // symmetrize
  for (int j = 1; j < n; j++) {
    for (int i = 0; i < j; i++) {
      a[i + j * n] = a[j + i * n];
    }
  }
   
  delete [] v;
  delete [] col_k;
  delete [] col_km;
  delete [] dd0;
}

#else

template<typename T>
void full_sym_2x2BK(int n, T *a, T *dd1,
		    int *pivot_width, int *permute)
{
  fprintf(stderr, "%s %d : general template is not implomented\n", 
	  __FILE__, __LINE__);
}

template
void full_sym_2x2BK<double>(int n, double *a,
			    double *dd1,
			    int *pivot_width, int *permute);
#endif
template
void full_sym_2x2BK<complex<double> >(int n, complex<double> *a,
				      complex<double> *dd1,
				      int *pivot_width, int *permute);
template
void full_sym_2x2BK<quadruple>(int n, quadruple *a, quadruple *dd1,
			       int *pivot_width, int *permute);
template
void full_sym_2x2BK<complex<quadruple> >(int n, complex<quadruple> *a,
					 complex<quadruple> *dd1,
					 int *pivot_width, int *permute);
//
template<typename T>
void C_gemm_symm(const int ncol, const int nrow, const T &alpha, 
		 const T *a, const int lda, 
		 const T *b, const int ldb,
		 const T &beta,
		 T *c, const int ldc)
{
// #define DGEMM_FOR_SYMM
#ifdef DGEMM_FOR_SYMM
  blas_gemm<T>(CblasTrans, CblasNoTrans, 
	       ncol, ncol, nrow, alpha, a, lda, b, ldb,
	       beta, c, ldc);
#else
  int n1, n2, n11, n12, n21, n22, n3;
  n1 = ncol / 2;
  n2 = ncol - n1;
  n11 = n1 / 2;
  n12 = n1 - n11;
  n21 = n2 / 2;
  n22 = n2 - n21;
  n3 = n1 + n21;
  if (ncol < SIZE_DGEMM_SYMM_DTRSV) {
    for (int j = 0; j < ncol; j++) {
      blas_gemv<T>(CblasTrans, nrow, (j + 1), alpha, a, lda,
		   b + (ldb * j), 1, beta, c + (ldc * j), 1);
    }
  }
  else {
    C_gemm_symm<T>(n11, nrow, alpha, a, lda, b, ldb, beta, c, ldc);
    C_gemm_symm<T>(n12, nrow, alpha, a + (lda * n11), lda, b + (ldb * n11), ldb, 
		beta, c + (ldc + 1) * n11, ldc);
    C_gemm_symm<T>(n21, nrow, alpha, a + lda * n1, lda, b + ldb * n1, ldb, 
		 beta, c + (ldc + 1) * n1, ldc);
    C_gemm_symm<T>(n22, nrow, alpha, a + lda * n3, lda, b + ldb * n3, ldb, 
		beta, c + (ldc + 1) * n3, ldc);
    blas_gemm<T>(CblasTrans, CblasNoTrans, 
	      n11, n12, nrow, alpha, a, lda, b + (ldb * n11), ldb,
	      beta, c + (ldc * n11), ldc);
    blas_gemm<T>(CblasTrans, CblasNoTrans, 
	      n21, n22, nrow, alpha, a + (lda * n1), lda, b + (ldb * n3), ldb,
	      beta, c + (n1 + ldc * n3), ldc);
    blas_gemm<T>(CblasTrans, CblasNoTrans, 
		 n1, n2, nrow, alpha, a, lda, b + (ldb * n1), ldb,
		 beta, c + (ldc * n1), ldc);
  }
#endif
}

template
void C_gemm_symm<double>(const int ncol, const int nrow, const double &alpha, 
			 const double *a, const int lda, 
			 const double *b, const int ldb,
			 const double &beta,
			 double *c, const int ldc);

template
void C_gemm_symm<complex<double> >(const int ncol, const int nrow, 
				   const complex<double> &alpha, 
				   const complex<double> *a, const int lda, 
				   const complex<double> *b, const int ldb,
				   const complex<double> &beta,
				   complex<double> *c, const int ldc);
template
void C_gemm_symm<quadruple>(const int ncol, const int nrow,
			    const quadruple &alpha, 
			    const quadruple *a, const int lda, 
			    const quadruple *b, const int ldb,
			    const quadruple &beta,
			    quadruple *c, const int ldc);

template
void C_gemm_symm<complex<quadruple> >(const int ncol, const int nrow, 
				   const complex<quadruple> &alpha, 
				   const complex<quadruple> *a, const int lda, 
				   const complex<quadruple> *b, const int ldb,
				   const complex<quadruple> &beta,
				   complex<quadruple> *c, const int ldc);
//
