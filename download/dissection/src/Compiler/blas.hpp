/*! \file   blas.hpp
    \brief  BLAS function interface
    \author Xavier Juvigny, ONERA
    \date   Jan. 12th 2005
    \modification function from NETLIB source with BLAS and CBLAS wrapper
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Mar. 16th 2015
    \date   Jul. 17th 2015
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

#ifndef _COMPILER_BLAS_H
# define _COMPILER_BLAS_H
# include "Compiler/OptionCompiler.hpp"
# include "Compiler/arithmetic.hpp"

#include <complex>

#ifdef BLAS_MKL
#define MKL_Complex16 std::complex<double>
#include <mkl_cblas.h>
#include <mkl_trans.h>
#endif
#ifdef VECLIB
#include <cblas.h>
#endif
#ifdef SX_ACE_BLAS
#include <cblas.h>
#endif
#ifdef OPENBLAS
#include "cblas.h"
#endif
#ifdef SUNPERF
#define floatcomplex std::complex<float>
#define doublecomplex std::complex<double>
#define _SUNPERF_COMPLEX
#include <sunperf.h>
#endif
#if (defined(BLAS_GENERIC) || defined(BLAS_FORTRAN))
typedef enum {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113 } CBLAS_TRANSPOSE;
typedef enum {CblasUpper=121, CblasLower=122} CBLAS_UPLO;
typedef enum {CblasNonUnit=131, CblasUnit=132} CBLAS_DIAG;
typedef enum {CblasLeft=141, CblasRight=142} CBLAS_SIDE;
#endif

using std::complex;

template<typename T, typename U>
inline U blas_abs(const T& x) {}

template<>
inline double blas_abs<double, double>(const double &x) {
  return fabs(x);
}

template<>
inline double
blas_abs<complex<double>, double>(const complex<double> &x) {
  return abs(x);
}

template<>
inline quadruple blas_abs<quadruple, quadruple>(const quadruple &x) {
  quadruple zero(0.0);
  return (x > zero ? x : (-x));
}

template<>
inline quadruple
blas_abs<complex<quadruple>, quadruple>(const complex<quadruple> &x) {
  quadruple xx = x.real();
  quadruple yy = x.imag();
  return sqrt(xx * xx + yy * yy);
}

#ifndef NO_OCTRUPLE
template<>
inline octruple blas_abs<octruple, octruple>(const octruple &x) {
  octruple zero(0.0);
  return (x > zero ? x : (-x));
}

template<>
inline octruple
blas_abs<complex<octruple>, octruple>(const complex<octruple> &x) {
  octruple xx = x.real();
  octruple yy = x.imag();
  return sqrt(xx * xx + yy * yy);
}
#endif

template<>
inline double blas_abs<quadruple, double>(const quadruple &x) {
  quadruple zero(0.0);
  return quad2double(x > zero ? x : (-x));
}

template<>
inline double
blas_abs<complex<quadruple>, double>(const complex<quadruple> &x) {
  quadruple xx = x.real();
  quadruple yy = x.imag();
  return quad2double(sqrt(xx * xx + yy * yy));
}

#ifndef NO_OCTRUPLE
template<>
inline double blas_abs<octruple, double>(const octruple &x) {
  octruple zero(0.0);
  return oct2double(x > zero ? x : (-x));
}

template<>
inline double
blas_abs<complex<octruple>, double>(const complex<octruple> &x) {
  octruple xx = x.real();
  octruple yy = x.imag();
  return oct2double(sqrt(xx * xx + yy * yy));
}
#endif

inline double
blas_conj(double x)
{
  return x;
}

inline quadruple
blas_conj(const quadruple &x)
{
  return x;
}

#ifndef NO_OCTRUPLE
inline octruple
blas_conj(const octruple &x)
{
  return x;
}
#endif
inline complex<double> 
blas_conj(const complex<double> &x)
{
  return std::conj(x);
}

inline complex<quadruple> 
blas_conj(const complex<quadruple> &x)
{
  return std::conj(x);
}

#ifndef NO_OCTRUPLE
inline complex<octruple> 
blas_conj(const complex<octruple> &x)
{
  return std::conj(x);
}
#endif
// ====================== Blas subroutine level 1 =========================

template<typename T>
void 
blas_axpy(const int n, const T &alpha, 
	  const T* x, const int incx, 
	  T* y, int incy);

#ifndef BLAS_GENERIC
template<>
void
blas_axpy<double>(const int n, const double &alpha, 
		  const double* x, const int incx, 
		  double* y, const int incy);
template<>
void 
blas_axpy<complex<double> >(const int n, const complex<double> &alpha, 
			    const complex<double>* x, const int incx, 
			    complex<double>* y, int incy);
#endif

template<typename T>
void
blas_copy(const int n, const T* x, const int incx, T* y, const int incy);

#ifndef BLAS_GENERIC
template<>
void
blas_copy<double>(const int n, const double* x, const int incx, double* y,
		  const int incy);

template<>
void
blas_copy<complex<double> >(const int n,
			    const complex<double>* x, const int incx,
			    complex<double>* y, const int incy);
#endif
template<typename T>
T
blas_dot(const int n, const T* x, const int incx, 
	 const T* y, const int incy);

#ifndef BLAS_GENERIC
template<>
double
blas_dot<double>(const int n, const double* x, const int incx, 
		 const double* y, const int incy);

template<>
complex<double>
blas_dot<complex<double> >(const int n,
			   const complex<double>* x, const int incx, 
			   const complex<double>* y, const int incy);
#endif
template<typename T>
void
blas_scal(const int N, const T &alpha, 
	  T *X, const int incX);

#ifndef BLAS_GENERIC
template<>
void
blas_scal<double>(const int N, const double &alpha, double *X, const int incX);

template<>
void
blas_scal<complex<double> >(const int N, const complex<double> &alpha, 
			    complex<double> *X, const int incX);
#endif
// T may be std::complex of U
template<typename T, typename U>
void
blas_scal2(const int N, const U &alpha, 
	   T *X, const int incX);
#ifndef BLAS_GENERIC
template<>
void
blas_scal2<double, double>(const int N, const double &alpha,
			   double *X, const int incX);

template<>
void
blas_scal2<complex<double>, double>(const int N, const double &alpha_, 
				    complex<double> *X, const int incX);

template<>
void
blas_scal2<quadruple, quadruple>(const int N, const quadruple &alpha,
				 quadruple *X, const int incX);

template<>
void
blas_scal2<complex<quadruple>, quadruple>(const int N, const quadruple &alpha_,
					  complex<quadruple> *X, const int incX);
#endif

template<typename T, typename U>
int blas_iamax(const int n, const T *x, const int incx);
#ifndef BLAS_GENERIC
template<>
int blas_iamax<double, double>(const int n, const double *x, const int incx);
template<>
int blas_iamax<complex<double>, double>(const int n,
					const complex<double> *x,
					const int incx);
#endif
// ====================== Blas subroutine level 2 =========================

template<typename T>
void
blas_gemv(const CBLAS_TRANSPOSE trA, 
	  const int m, const int n, const T &alpha, 
	  const T* A, const int lda, const T* x, const int incx, 
	  const T &beta, T* y, const int incy);
#ifndef BLAS_GENERIC
template<>
void
blas_gemv<double>(const CBLAS_TRANSPOSE trA, 
		  const int m, const int n, const double &alpha, 
		  const double* A, const int lda, const double* x,
		  const int incx, 
		  const double &beta, double* y, const int incy);

template<>
void
blas_gemv<complex<double> >(const CBLAS_TRANSPOSE trA, 
			    const int m, const int n,
			    const complex<double> &alpha, 
			    const complex<double>* A, const int lda, 
			    const complex<double>* x, const int incx, 
			    const complex<double> &beta,
			    complex<double>* y, const int incy);
#endif
template<typename T>
void
blas_trsv(const CBLAS_UPLO Uplo,
	  const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
	  const int N, const T *A, const int lda, T *X,
	  const int incX);
#ifndef BLAS_GENERIC
template<>
void
blas_trsv<double>(const CBLAS_UPLO Uplo,
		  const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
		  const int N, const double *A, const int lda, double *X,
		  const int incX);

template<>
void
blas_trsv<complex<double> >(const CBLAS_UPLO Uplo,
			    const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
			    const int N, const complex<double>* A,
			    const int lda, 
			    complex<double> *X, const int incX);
#endif

template<typename T>
void
blas_syr(const CBLAS_UPLO Uplo,
	 const int N, const T &alpha,
	 const T *X, const int incX,
	 T *A, const int lda);
#ifndef BLAS_GENERIC
template<>
void
blas_syr<double>(const CBLAS_UPLO Uplo, 
		 const int N, const double &alpha,
		 const double *X, const int incX,
		 double *A, const int lda);

template<>
void
blas_syr<complex<double> >(const CBLAS_UPLO Uplo, 
			   const int N, const complex<double> &alpha,
			   const complex<double> *X, const int incX,
			   complex<double> *A, const int lda);
#endif
template<typename T>
void
blas_syrc(const CBLAS_UPLO Uplo,
	  const int N, const T &alpha,
	  const T *X, const int incX,
	  T *A, const int lda);

template<typename T>
void
blas_syr2(const CBLAS_UPLO Uplo,
	  const int N, const T &alpha,
	  const T *X, const int incX,
	  const T *Y, const int incY,
	  T *A, const int lda);

#ifndef BLAS_GENERIC
template<>
void
blas_syr2<double>(const CBLAS_UPLO Uplo, const int N,
		  const double &alpha,
		  const double *X, const int incX,
		  const double *Y, const int incY,
		  double *A, const int lda);
#endif

template<typename T>
void
blas_ger(const int M, const int N, const T &alpha,
	 const T *X, const int incX,
	 const T *Y, const int incY,
	 T *A, const int lda);

#ifndef BLAS_GENERIC
template<>
void
blas_ger<double>(const int M, const int N, const double &alpha,
		 const double *X, const int incX,
		 const double *Y, const int incY,
		 double *A, const int lda);

template<>
void
blas_ger<complex<double> >(const int M, const int N,
			   const complex<double> &alpha,
			   const complex<double> *X, const int incX,
			   const complex<double> *Y, const int incY,
			   complex<double> *A, const int lda);
#endif

template<typename T>
void
blas_gerc(const int M, const int N, const T &alpha,
	  const T *X, const int incX,
	  const T *Y, const int incY,
	  T *A, const int lda);

#ifndef BLAS_GENERIC
template<>
void
blas_gerc<double>(const int M, const int N, const double &alpha,
		  const double *X, const int incX,
		  const double *Y, const int incY,
		  double *A, const int lda);
template<>
void
blas_gerc<complex<double> >(const int M, const int N,
			    const complex<double> &alpha,
			    const complex<double> *X, const int incX,
			    const complex<double> *Y, const int incY,
			    complex<double> *A, const int lda);
#endif
// ====================== Blas subroutine level 3 =========================

template<typename T>
void
blas_gemm( CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
	   int m, int n, int k,
	   const T &alpha,
	   const T* A, int lda,
	   const T* B, int ldb,
	   const T &beta,
	   T* C, int ldc );

#ifndef BLAS_GENERIC
template<>
void
blas_gemm<double>(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
		  int m, int n, int k,
		  const double &alpha, const double* A, int lda,
		  const double* B, int ldb, const double &beta,
		  double* C, int ldc );

template<>
void
blas_gemm<complex<double> >(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
			    int m, int n, int k, const complex<double> &alpha, 
			    const complex<double>* A, int lda,
			    const complex<double>* B, int ldb,
			    const complex<double> &beta, 
			    complex<double>* C, int ldc );
#endif

template<typename T>
void
blas_trsm(const CBLAS_SIDE Side,
	  const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
	  const CBLAS_DIAG Diag, const int M, const int N,
	  const T &alpha,
	  const T *A, const int lda,
	  T *B, const int ldb);
#ifndef BLAS_GENERIC
template<>
void
blas_trsm<double>(const CBLAS_SIDE Side,
		  const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
		  const CBLAS_DIAG Diag, const int M, const int N,
		  const double &alpha, const double *A, const int lda,
		  double *B, const int ldb);

template<>
void
blas_trsm<complex<double> >(const CBLAS_SIDE Side,
			    const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
			    const CBLAS_DIAG Diag, const int M, const int N,
			    const complex<double> &alpha,
			    const complex<double> *A, const int lda,
			    complex<double> *B, const int ldb);
#endif
#ifdef BLAS_MKL
// ====================== MKL transposition routines =======================
inline void 
mkl_omatcopy(char ordering, char trans,
	     size_t rows, size_t cols,
	     const double alpha,
	     const double *A, size_t lda,
	     double *B, size_t ldb)
{
  mkl_domatcopy(ordering, trans,
		rows, cols,
		alpha, A, lda,
		B, ldb);
}

inline void 
mkl_omatcopy(char ordering, char trans,
	     size_t rows, size_t cols,
	     const double alpha_,
	     const complex<double> *A, size_t lda,
	     complex<double> *B, size_t ldb)
{
  const complex<double> alpha = complex<double>(alpha_, 0.0);
  mkl_zomatcopy(ordering, trans,
		rows, cols,
		(MKL_Complex16)alpha, (MKL_Complex16 *)A, lda,
		(MKL_Complex16 *)B, ldb);
}
#endif


// ======== for computation of norm and norm^2

template<typename T, typename U>
U blas_l2norm(const int n, T *x, const int incX);

template<>
double blas_l2norm<complex<double>, double>(const int n, complex<double> *x,
					    const int incX);
template<>
quadruple blas_l2norm<complex<quadruple>, quadruple>(const int n,
						     complex<quadruple> *x,
						     const int incX);

template<typename T, typename U>
U blas_l2norm2(const int n, T *x, const int incX);

template<>
double blas_l2norm2<complex<double>, double>(const int n, complex<double> *x,
					     const int incX);

template<>
quadruple blas_l2norm2<complex<quadruple>, quadruple>(const int n,
						      complex<quadruple> *x,
						      const int incX);
#endif

