/*! \file   blas.cpp
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

#include "Compiler/blas.hpp"
#include <cstdio>
#include <cstring>

#define FORCE_EXPLICIT_INSTANTIATION
#if __clang_major__ == 4
#undef FORCE_EXPLICIT_INSTANTIATION
#endif

#ifdef BLAS_FORTRAN

extern "C" {
  // double precision
  extern int FORTRAN_DECL(idamax)(const int &n,     // return Fortran index
				  const double *x,
				  const int &incx);
  
  extern int FORTRAN_DECL(izamax)(const int &n,     // retrun Fortran index
				  const double *x,
				  const int &incx);
  
  extern void FORTRAN_DECL(dcopy)(const int &n,
				  const double *x, const int &incx,
				  double *y, const int &incy);

  extern void FORTRAN_DECL(zcopy)(const int &n,
				  const double *x, const int &incx,
				  double *y, const int &incy);

  extern void FORTRAN_DECL(daxpy)(const int &n, const double &alpha,
				  const double *x, const int &incx,
				  double *y, const int &incy);

  extern void FORTRAN_DECL(zaxpy)(const int &n, const double *alpha,
				  const double *x, const int &incx,
				  double *y, const int &incy);

  extern double FORTRAN_DECL(ddot)(const int &n,
				   const double *x, const int &incx,
				   const double *y, const int &incy);
  
  extern void FORTRAN_DECL(zdotc)(double *z,
				  const int &n,
				  const double *x, const int &incx,
				  const double *y, const int &incy);

  extern void FORTRAN_DECL(dscal)(const int &N, const double &alpha,
				  double *X, const int &incX);

  extern void FORTRAN_DECL(zscal)(const int &N, const double *alpha,
				  double *X, const int &incX);

  extern void FORTRAN_DECL(dgemv)(unsigned char &tra_,
				  const int &m,
				  const int &n,
				  const double &alpha,
				  const double *A,
				  const int &lda,
				  const double *x,
				  const int &incx,
				  const double &beta,
				  double *y,
				  const int &incy);
  extern void FORTRAN_DECL(zgemv)(unsigned char &tra_,
				  const int &m,
				  const int &n,
				  const double *alpha,
				  const double *A,
				  const int &lda,
				  const double *x,
				  const int &incx,
				  const double *beta,
				  double *y,
				  const int &incy);
  
  extern void FORTRAN_DECL(dtrsv)(unsigned char &uplo_,
				  unsigned char &transa_,
				  unsigned char &diag_,
				  const int &N,
				  const double *A,
				  const int &lda,
				  double *X,
				  const int &incX);
  
  extern void FORTRAN_DECL(ztrsv)(unsigned char &uplo_,
				  unsigned char &transa_,
				  unsigned char &diag_,
				  const int &N,
				  const double *A,
				  const int &lda,
				  double *X,
				  const int &incX);
  
  extern void FORTRAN_DECL(dsyr)(unsigned char &uplo_,
				 const int &N,
				 const double &alpha,
				 const double *X,
				 const int &incX,
				 double *A,
				 const int &lda);
  
  extern void FORTRAN_DECL(dsyr2)(unsigned char &uplo_,
				  const int &N,
				  const double &alpha,
				  const double *X,
				  const int &incX,
				  const double *Y,
				  const int &incY,
				  double *A,
				  const int &lda);
  extern void FORTRAN_DECL(dger)(const int &M, const int &N,
				 const double &alpha,
				 const double *X, const int &incX,
				 const double *Y, const int &incY,
				 double *A, const int  &lda);
  
  extern void FORTRAN_DECL(zgeru)(const int &M, const int &N,
				  const double *alpha,
				  double *X, const int &incX,
				  double *Y, const int &incY, 
				  double *A, const int &lda); 
  
  extern void FORTRAN_DECL(zgerc)(const int &M, const int &N,
				  const double *alpha,
				  const double *X, const int &incX,
				  const double *Y, const int &incY, 
				  double *A, const int &lda); 

  extern void FORTRAN_DECL(dtrsm)(unsigned char &side_, unsigned char &uplo_,
				  unsigned char &transa_, unsigned char &diag_,
				  const int &M, const int &N,
				  const double &alpha,
				  const double *A, const int &lda,
				  double *B, const int &ldb);
  
  extern void FORTRAN_DECL(ztrsm)(unsigned char &side_, unsigned char &uplo_,
				  unsigned char &transa_, unsigned char &diag_,
				  const int &M, const int &N,
				  const double *alpha,
				  const double *A, const int &lda,
				  double *B, const int &ldb);
  
  extern void FORTRAN_DECL(dgemm)(unsigned char &tra_, unsigned char &trb_,
				  const int &m, const int &n, const int &k,
				  const double &alpha,
				  const double *A, const int &lda,
				  const double *B, const int &ldb,
				  const double &beta,
				  double *C, const int &ldc);
  
  extern void FORTRAN_DECL(zgemm)(unsigned char &tra_, unsigned char &trb_,
				  const int &m, const int &n, const int &k,
				  const double *alpha,
				  const double *A, const int &lda,
				  const double *B, const int &ldb,
				  const double *beta,
				  double *C, const int &ldc);

  // single precision
  extern int FORTRAN_DECL(isamax)(const int &n,     // return Fortran index
				  const float *x,
				  const int &incx);
  
  extern int FORTRAN_DECL(icamax)(const int &n,     // retrun Fortran index
				  const float *x,
				  const int &incx);
  
  extern void FORTRAN_DECL(scopy)(const int &n,
				  const float *x, const int &incx,
				  float *y, const int &incy);

  extern void FORTRAN_DECL(ccopy)(const int &n,
				  const float *x, const int &incx,
				  float *y, const int &incy);

  extern void FORTRAN_DECL(saxpy)(const int &n, const float &alpha,
				  const float *x, const int &incx,
				  float *y, const int &incy);

  extern void FORTRAN_DECL(caxpy)(const int &n, const float *alpha,
				  const float *x, const int &incx,
				  float *y, const int &incy);

  extern float FORTRAN_DECL(sdot)(const int &n,
				   const float *x, const int &incx,
				   const float *y, const int &incy);
  
  extern void FORTRAN_DECL(cdotc)(float *z,
				  const int &n,
				  const float *x, const int &incx,
				  const float *y, const int &incy);

  extern void FORTRAN_DECL(sscal)(const int &N, const float &alpha,
				  float *X, const int &incX);

  extern void FORTRAN_DECL(cscal)(const int &N, const float *alpha,
				  float *X, const int &incX);

  extern void FORTRAN_DECL(sgemv)(unsigned char &tra_,
				  const int &m,
				  const int &n,
				  const float &alpha,
				  const float *A,
				  const int &lda,
				  const float *x,
				  const int &incx,
				  const float &beta,
				  float *y,
				  const int &incy);
  extern void FORTRAN_DECL(cgemv)(unsigned char &tra_,
				  const int &m,
				  const int &n,
				  const float *alpha,
				  const float *A,
				  const int &lda,
				  const float *x,
				  const int &incx,
				  const float *beta,
				  float *y,
				  const int &incy);
  
  extern void FORTRAN_DECL(strsv)(unsigned char &uplo_,
				  unsigned char &transa_,
				  unsigned char &diag_,
				  const int &N,
				  const float *A,
				  const int &lda,
				  float *X,
				  const int &incX);
  
  extern void FORTRAN_DECL(ctrsv)(unsigned char &uplo_,
				  unsigned char &transa_,
				  unsigned char &diag_,
				  const int &N,
				  const float *A,
				  const int &lda,
				  float *X,
				  const int &incX);
  
  extern void FORTRAN_DECL(ssyr)(unsigned char &uplo_,
				 const int &N,
				 const float &alpha,
				 const float *X,
				 const int &incX,
				 float *A,
				 const int &lda);
  
  extern void FORTRAN_DECL(ssyr2)(unsigned char &uplo_,
				  const int &N,
				  const float &alpha,
				  const float *X,
				  const int &incX,
				  const float *Y,
				  const int &incY,
				  float *A,
				  const int &lda);
  extern void FORTRAN_DECL(sger)(const int &M, const int &N,
				 const float &alpha,
				 const float *X, const int &incX,
				 const float *Y, const int &incY,
				 float *A, const int  &lda);
  
  extern void FORTRAN_DECL(cgeru)(const int &M, const int &N,
				  const float *alpha,
				  float *X, const int &incX,
				  float *Y, const int &incY, 
				  float *A, const int &lda); 
  
  extern void FORTRAN_DECL(cgerc)(const int &M, const int &N,
				  const float *alpha,
				  const float *X, const int &incX,
				  const float *Y, const int &incY, 
				  float *A, const int &lda); 

  extern void FORTRAN_DECL(strsm)(unsigned char &side_, unsigned char &uplo_,
				  unsigned char &transa_, unsigned char &diag_,
				  const int &M, const int &N,
				  const float &alpha,
				  const float *A, const int &lda,
				  float *B, const int &ldb);
  
  extern void FORTRAN_DECL(ctrsm)(unsigned char &side_, unsigned char &uplo_,
				  unsigned char &transa_, unsigned char &diag_,
				  const int &M, const int &N,
				  const float *alpha,
				  const float *A, const int &lda,
				  float *B, const int &ldb);
  
  extern void FORTRAN_DECL(sgemm)(unsigned char &tra_, unsigned char &trb_,
				  const int &m, const int &n, const int &k,
				  const float &alpha,
				  const float *A, const int &lda,
				  const float *B, const int &ldb,
				  const float &beta,
				  float *C, const int &ldc);
  
  extern void FORTRAN_DECL(cgemm)(unsigned char &tra_, unsigned char &trb_,
				  const int &m, const int &n, const int &k,
				  const float *alpha,
				  const float *A, const int &lda,
				  const float *B, const int &ldb,
				  const float *beta,
				  float *C, const int &ldc);

}
#endif
// BLAS 1
// dz copy
#ifndef BLAS_GENERIC
template<>
void
blas_copy<double>(const int n, const double* x, const int incx, double* y,
		  const int incy)
{
  if ((incx == 1) && (incy == 1)) {
    memcpy((void *)y, (void *)x, n * sizeof(double));
  }
  else {
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(dcopy)(n, x, incx, y, incy);
#else
  cblas_dcopy((BLAS_INT)n, x, (BLAS_INT)incx, y, (BLAS_INT)incy);
#endif
  }
}

template<>
void
blas_copy<complex<double> >(const int n,
			    const complex<double>* x, const int incx,
			    complex<double>* y, const int incy)
{
  if ((incx == 1) && (incy == 1)) {
    memcpy((void *)y, (void *)x, n * sizeof(complex<double>));
  }
  else {
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(zcopy)(n, (double *)x, incx, (double *)y, incy);
#else
    cblas_zcopy((BLAS_INT)n, (void *)x, (BLAS_INT)incx,
		(BLAS_VOID *)y, (BLAS_INT)incy);
#endif
  }
}

template<>
void
blas_copy<float>(const int n, const float* x, const int incx, float* y,
		  const int incy)
{
  if ((incx == 1) && (incy == 1)) {
    memcpy((void *)y, (void *)x, n * sizeof(float));
  }
  else {
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(scopy)(n, x, incx, y, incy);
#else
  cblas_scopy((BLAS_INT)n, x, (BLAS_INT)incx, y, (BLAS_INT)incy);
#endif
  }
}

template<>
void
blas_copy<complex<float> >(const int n,
			    const complex<float>* x, const int incx,
			    complex<float>* y, const int incy)
{
  if ((incx == 1) && (incy == 1)) {
    memcpy((void *)y, (void *)x, n * sizeof(complex<float>));
  }
  else {
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(ccopy)(n, (float *)x, incx, (float *)y, incy);
#else
    cblas_ccopy((BLAS_INT)n, (void *)x, (BLAS_INT)incx,
		(BLAS_VOID *)y, (BLAS_INT)incy);
#endif
  }
}

#endif
template<typename T>
void
blas_copy(const int n, const T* x, const int incx, T *y, const int incy)
{
  if ((incx == 1) && (incy == 1)) {
    memcpy((void *)y, (void *)x, n * sizeof(T));
  }
  else {
    int ix = 0;
    int iy = 0;
    for (int i = 0; i < n; i++, ix += incx, iy += incy) {
    y[iy] = x[ix];
    }
  }
}

#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_copy<double>(const int n, const double* x, const int incx,
		  double *y, const int incy);

template
void
blas_copy<complex<double> > (const int n,
			     const complex<double>* x, const int incx,
			     complex<double> *y, const int incy);
template
void
blas_copy<float>(const int n, const float* x, const int incx,
		  float *y, const int incy);

template
void
blas_copy<complex<float> > (const int n,
			     const complex<float>* x, const int incx,
			     complex<float> *y, const int incy);
#endif
template
void
blas_copy<quadruple>(const int n, const quadruple* x, const int incx,
		     quadruple *y, const int incy);
template
void
blas_copy<complex<quadruple> > (const int n,
				const complex<quadruple>* x, const int incx,
				complex<quadruple> *y, const int incy);

template
void
blas_copy<octruple>(const int n, const octruple* x, const int incx,
		  octruple *y, const int incy);

template
void
blas_copy<complex<octruple> > (const int n,
			     const complex<octruple>* x, const int incx,
			     complex<octruple> *y, const int incy);

// dz axpy
#ifndef BLAS_GENERIC
template<>
void
blas_axpy<double>(const int n, const double &alpha, 
		  const double* x, const int incx, 
		  double* y, const int incy)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(daxpy)(n, alpha, x, incx, y, incy);
#else
  cblas_daxpy((BLAS_INT)n, alpha, x, (BLAS_INT)incx, y, (BLAS_INT)incy);
#endif
}

template<>
void 
blas_axpy<complex<double> >(const int n, const complex<double> &alpha, 
			    const complex<double>* x, const int incx, 
			    complex<double>* y, int incy)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(zaxpy)(n, (double *)&alpha, (double *)x, incx, (double *)y, incy);
#else
  cblas_zaxpy((BLAS_INT)n, (BLAS_VOID *)&alpha, (BLAS_VOID *)x, (BLAS_INT)incx,
	      (BLAS_VOID *)y, (BLAS_INT)incy);
#endif
}

template<>
void
blas_axpy<float>(const int n, const float &alpha, 
		  const float* x, const int incx, 
		  float* y, const int incy)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(saxpy)(n, alpha, x, incx, y, incy);
#else
  cblas_saxpy((BLAS_INT)n, alpha, x, (BLAS_INT)incx, y, (BLAS_INT)incy);
#endif
}

template<>
void 
blas_axpy<complex<float> >(const int n, const complex<float> &alpha, 
			    const complex<float>* x, const int incx, 
			    complex<float>* y, int incy)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(caxpy)(n, (float *)&alpha, (float *)x, incx, (float *)y, incy);
#else
  cblas_caxpy((BLAS_INT)n, (BLAS_VOID *)&alpha, (BLAS_VOID *)x, (BLAS_INT)incx,
	      (BLAS_VOID *)y, (BLAS_INT)incy);
#endif
}
#endif

template<typename T>
void 
blas_axpy(const int n, const T &alpha, 
	  const T* x, const int incx, 
	  T* y, int incy)
{
  int ix = 0;
  int iy = 0;
  for (int i = 0; i < n; i++, ix += incx, iy += incy) {
    y[iy] += alpha * x[ix];
  }
}
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void 
blas_axpy<double>(const int n, const double &alpha, 
		  const double* x, const int incx, 
		  double* y, int incy);
template
void 
blas_axpy<complex<double> >(const int n, const complex<double> &alpha, 
			    const complex<double>* x, const int incx, 
			    complex<double>* y, int incy);

template
void 
blas_axpy<float>(const int n, const float &alpha, 
		  const float* x, const int incx, 
		  float* y, int incy);
template
void 
blas_axpy<complex<float> >(const int n, const complex<float> &alpha, 
			    const complex<float>* x, const int incx, 
			    complex<float>* y, int incy);
#endif
template
void 
blas_axpy<quadruple>(const int n, const quadruple &alpha, 
		     const quadruple* x, const int incx, 
		     quadruple* y, int incy);
template
void 
blas_axpy<complex<quadruple> >(const int n, const complex<quadruple> &alpha, 
			       const complex<quadruple>* x, const int incx, 
			       complex<quadruple>* y, int incy);

// dz dot
#ifndef BLAS_GENERIC
template<>
double
blas_dot<double>(const int n, const double* x, const int incx, 
		 const double* y, const int incy)
{
#ifdef BLAS_FORTRAN
  return FORTRAN_DECL(ddot)(n, x, incx, y, incy); // ?
#else
  return cblas_ddot((BLAS_INT)n, x, (BLAS_INT)incx, y, (BLAS_INT)incy);
#endif
}

template<>
complex<double>
blas_dot<complex<double> >(const int n,
			   const complex<double>* x, const int incx, 
			   const complex<double>* y, const int incy)
{
  complex<double> val;
#ifdef BLAS_FORTRAN
  // fortran function that returns complex is called as C void function
  FORTRAN_DECL(zdotc)((double*) &val, n, (double *)x, incx, (double *)y, incy);
#else
  cblas_zdotc_sub((BLAS_INT)n, (BLAS_VOID *)x, (BLAS_INT)incx,
		  (BLAS_VOID *)y, (BLAS_INT)incy, (BLAS_VOID *)&val);
#endif
  return val;
}

template<>
float
blas_dot<float>(const int n, const float* x, const int incx, 
		 const float* y, const int incy)
{
#ifdef BLAS_FORTRAN
  return FORTRAN_DECL(sdot)(n, x, incx, y, incy); // ?
#else
  return cblas_sdot((BLAS_INT)n, x, (BLAS_INT)incx, y, (BLAS_INT)incy);
#endif
}

template<>
complex<float>
blas_dot<complex<float> >(const int n,
			   const complex<float>* x, const int incx, 
			   const complex<float>* y, const int incy)
{
  complex<float> val;
#ifdef BLAS_FORTRAN
  // fortran function that returns complex is called as C void function
  FORTRAN_DECL(cdotc)((float*) &val, n, (float *)x, incx, (float *)y, incy);
#else
  cblas_cdotc_sub((BLAS_INT)n, (BLAS_VOID *)x, (BLAS_INT)incx,
		  (BLAS_VOID *)y, (BLAS_INT)incy, (BLAS_VOID *)&val);
#endif
  return val;
}

#endif

template<typename T>
T
blas_dot(const int n, const T* x, const int incx, 
	 const T* y, const int incy)
{
  int ix, iy;
  T temp;
  temp = T(0.0);
  if (n <= 0) {
    return temp;
  }
  if ((incx == 1) && (incy == 1)) {
    for (int i = 0; i < n; i++) {
      temp += blas_conj(x[i]) * y[i];
    }
  }
  else {
    ix = 0;
    iy = 0;
    if (incx < 0) {
      ix = (n - 1)* (-incx);  // decrement form the last index
    }
    if (incy < 0) {
      iy = (n - 1)* (-incy);  // decrement form the last index
    }
    for (int i = 0; i < n; i++) {
      temp += blas_conj(x[i]) * y[i];
      ix += incx;
      iy += incy;
    }
  }
  return temp;
}
// explicit instantiation of blas_dot
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
double
blas_dot<double>(const int n, const double* x, const int incx, 
		 const double* y, const int incy);
template
complex<double>
blas_dot<complex<double> >(const int n,
			   const complex<double>* x, const int incx, 
			   const complex<double>* y, const int incy);
template
float
blas_dot<float>(const int n, const float* x, const int incx, 
		 const float* y, const int incy);
template
complex<float>
blas_dot<complex<float> >(const int n,
			   const complex<float>* x, const int incx, 
			   const complex<float>* y, const int incy);
#endif
template
quadruple
blas_dot<quadruple>(const int n, const quadruple* x, const int incx, 
		    const quadruple* y, const int incy);
template
complex<quadruple>
blas_dot<complex<quadruple> >(const int n,
			      const complex<quadruple>* x, const int incx, 
			      const complex<quadruple>* y, const int incy);
#ifndef NO_OCTRUPLE
template
octruple
blas_dot<octruple>(const int n, const octruple* x, const int incx, 
		   const octruple* y, const int incy);
template
complex<octruple>
blas_dot<complex<octruple> >(const int n,
			     const complex<octruple>* x, const int incx, 
			     const complex<octruple>* y, const int incy);
#endif
// dz scal
#ifndef BLAS_GENERIC
template<>
void
blas_scal<double>(const int N, const double &alpha, double *X, const int incX)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(dscal)(N, alpha, X, incX);
#else
  cblas_dscal((BLAS_INT)N, alpha, X, (BLAS_INT)incX);
#endif
}

template<>
void
blas_scal<complex<double> >(const int N, const complex<double> &alpha, 
			    complex<double> *X, const int incX)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(zscal)(N, (double *)&alpha, (double *)X, incX);
#else
  cblas_zscal((BLAS_INT)N, (BLAS_VOID *)&alpha, (BLAS_VOID *)X, (BLAS_INT)incX);
#endif
}

template<>
void
blas_scal<float>(const int N, const float &alpha, float *X, const int incX)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(sscal)(N, alpha, X, incX);
#else
  cblas_sscal((BLAS_INT)N, alpha, X, (BLAS_INT)incX);
#endif
}

template<>
void
blas_scal<complex<float> >(const int N, const complex<float> &alpha, 
			    complex<float> *X, const int incX)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(cscal)(N, (float *)&alpha, (float *)X, incX);
#else
  cblas_cscal((BLAS_INT)N, (BLAS_VOID *)&alpha, (BLAS_VOID *)X, (BLAS_INT)incX);
#endif
}
#endif // #ifdef BLAS_GENERIC

template<typename T>
void
blas_scal(const int N, const T &alpha, 
	  T *X, const int incX)
{
  if ((N <= 0) || (incX <= 0)) {
    return;
  }
  if (incX == 1) {
    for (int i = 0; i < N; i++) {
      X[i] *= alpha;
    }
  }
  else {
    int nincx = N * incX;
    for (int i = 0; i < nincx; i += incX) {
      X[i] *= alpha;
    }
  }
}
// explicit instantiation of blas_scal
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_scal<double>(const int N, const double &alpha,
		  double *X, const int incX);

template
void
blas_scal<complex<double> >(const int N, const complex<double> &alpha,
			    complex<double> *X, const int incX);

template
void
blas_scal<float>(const int N, const float &alpha,
		  float *X, const int incX);

template
void
blas_scal<complex<float> >(const int N, const complex<float> &alpha,
			    complex<float> *X, const int incX);
#endif
template
void
blas_scal<quadruple>(const int N, const quadruple &alpha,
		     quadruple *X, const int incX);

template
void
blas_scal<complex<quadruple> >(const int N, const complex<quadruple> &alpha,
			       complex<quadruple> *X, const int incX);

#ifndef NO_OCTRUPLE
template
void
blas_scal<octruple>(const int N, const octruple &alpha,
		     octruple *X, const int incX);

template
void
blas_scal<complex<octruple> >(const int N, const complex<octruple> &alpha,
			       complex<octruple> *X, const int incX);
#endif
// scaling with magnitude represented by real valued
#ifndef BLAS_GENERIC
template<>
void
blas_scal2<double, double>(const int N, const double &alpha,
			  double *X, const int incX)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(dscal)(N, alpha, X, incX);
#else
  cblas_dscal((BLAS_INT)N, alpha, X, (BLAS_INT)incX);
#endif
}

template<>
void
blas_scal2<complex<double>, double>(const int N, const double &alpha_, 
				   complex<double> *X, const int incX)
{
  const complex<double> alpha = complex<double>(alpha_, 0.0);
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(zscal)(N, (double *)&alpha, (double *)X, incX);
#else
  cblas_zscal((BLAS_INT)N, (BLAS_VOID *)&alpha, (BLAS_VOID *)X, (BLAS_INT)incX);
#endif
}
template<>
void
blas_scal2<float, float>(const int N, const float &alpha,
			  float *X, const int incX)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(sscal)(N, alpha, X, incX);
#else
  cblas_sscal((BLAS_INT)N, alpha, X, (BLAS_INT)incX);
#endif
}

template<>
void
blas_scal2<complex<float>, float>(const int N, const float &alpha_, 
				   complex<float> *X, const int incX)
{
  const complex<float> alpha = complex<float>(alpha_, 0.0);
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(cscal)(N, (float *)&alpha, (float *)X, incX);
#else
  cblas_cscal((BLAS_INT)N, (BLAS_VOID *)&alpha, (BLAS_VOID *)X, (BLAS_INT)incX);
#endif
}

template<>
void
blas_scal2<quadruple, quadruple>(const int N, const quadruple &alpha,
				quadruple *X, const int incX)
{
  blas_scal<quadruple>(N, alpha, X, incX);
}

template<>
void
blas_scal2<complex<quadruple>, quadruple>(const int N, const quadruple &alpha_,
					  complex<quadruple> *X, const int incX)
{
  complex<quadruple> alpha = complex<quadruple>(alpha_, quadruple(0.0));
  blas_scal<complex<quadruple> >(N, alpha, X, incX);
}

#ifndef NO_OCTRUPLE
template<>
void
blas_scal2<octruple, octruple>(const int N, const octruple &alpha,
			      octruple *X, const int incX)
{
  blas_scal<octruple>(N, alpha, X, incX);
}

template<>
void
blas_scal2<complex<octruple>, octruple>(const int N, const octruple &alpha_,
					complex<octruple> *X, const int incX)
{
  complex<octruple> alpha = complex<octruple>(alpha_, octruple(0.0));
  blas_scal<complex<octruple> >(N, alpha, X, incX);
}

#endif
#endif // ifndef BLAS_GENERIC
template<typename T, typename U>
void
blas_scal2(const int N, const U &alpha_, 
	  T *X, const int incX)
{
  T alpha = tocomplex<T, U>(alpha_);  
  if ((N <= 0) || (incX <= 0)) {
    return;
  }
  if (incX == 1) {
    for (int i = 0; i < N; i++) {
      X[i] *= alpha;
    }
  }
  else {
    int nincx = N * incX;
    for (int i = 0; i < nincx; i += incX) {
      X[i] *= alpha;
    }
  }
}

// explicit instantiation blas_scal with real number saling
#ifdef FORCE_EXPLICIT_INSTANTIATION
template void
blas_scal2<double, double>(const int N, const double &alpha,
			  double *X, const int incX);

template void
blas_scal2<complex<double>, double>(const int N, const double &alpha_,
				   complex<double> *X, const int incX);

template void
blas_scal2<float, float>(const int N, const float &alpha,
			  float *X, const int incX);

template void
blas_scal2<complex<float>, float>(const int N, const float &alpha_,
				   complex<float> *X, const int incX);


template void
blas_scal2<quadruple, quadruple>(const int N, const quadruple &alpha,
				quadruple *X, const int incX);
template void
blas_scal2<complex<quadruple>, quadruple>(const int N, const quadruple &alpha_,
					 complex<quadruple> *X, const int incX);

#ifndef NO_OCTRUPLE
template void
blas_scal2<octruple, octruple>(const int N, const octruple &alpha,
			      octruple *X, const int incX);
template void
blas_scal2<complex<octruple>, octruple>(const int N, const octruple &alpha,
				       complex<octruple> *X, const int incX);
#endif
#endif
template<typename T, typename U>
int blas_iamax(const int n, const T *x, const int incx)
{
  U umax;
  int ix, iamax;
  iamax = (-1);
  if ((n < 1) || incx <= 0) {
    return iamax;
  }
  iamax = 0;
  if (n == 1) {
    return iamax;
  }
  if (incx == 1) {
    umax = blas_abs<T, U>(x[0]);
    for (int i = 1; i < n; i++) {
      U tmax = blas_abs<T, U>(x[0]);
      if (tmax > umax) {
	iamax = i;
	umax = tmax;
      }
    }
  }
  else {
    ix = 0;
    umax = blas_abs<T, U>(x[0]);
    for (int i = 1; i < n; i++) {
      ix += incx;
      U tmax = blas_abs<T, U>(x[ix]);
      if (tmax > umax) {
	iamax = i;
	umax = tmax;
      }
    }
  }
  return iamax;
}

#ifndef BLAS_GENERIC
template<>
int blas_iamax<double, double>(const int n, const double *x, const int incx)
{
#ifdef BLAS_FORTRAN
  return (FORTRAN_DECL(idamax)(n, x, incx) - 1); // Fortran to C index
#else
  return cblas_idamax((BLAS_INT)n, x, (BLAS_INT)incx);
#endif
}

template<>
int blas_iamax<complex<double>, double>(const int n,
					const complex<double> *x,
					const int incx)
{
#ifdef BLAS_FORTRAN                                  // Fortran to C index
  return (FORTRAN_DECL(izamax)(n, (double *)x, incx) - 1);
#else
  return cblas_izamax((BLAS_INT)n, (double *)x, (BLAS_INT)incx);
#endif
}

template<>
int blas_iamax<float, float>(const int n, const float *x, const int incx)
{
#ifdef BLAS_FORTRAN
  return (FORTRAN_DECL(isamax)(n, x, incx) - 1); // Fortran to C index
#else
  return cblas_isamax((BLAS_INT)n, x, (BLAS_INT)incx);
#endif
}

template<>
int blas_iamax<complex<float>, float>(const int n,
					const complex<float> *x,
					const int incx)
{
#ifdef BLAS_FORTRAN                                  // Fortran to C index
  return (FORTRAN_DECL(icamax)(n, (float *)x, incx) - 1);
#else
  return cblas_icamax((BLAS_INT)n, (float *)x, (BLAS_INT)incx);
#endif
}

#endif
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
int blas_iamax<double, double>(const int n, const double *x, const int incx);
template
int blas_iamax<complex<double>, double>(const int n, const complex<double> *x,
					const int incx);
template
int blas_iamax<float, float>(const int n, const float *x, const int incx);
template
int blas_iamax<complex<float>, float>(const int n, const complex<float> *x,
					const int incx);
#endif
template
int blas_iamax<quadruple, quadruple>(const int n, const quadruple *x, const int incx);
template
int blas_iamax<complex<quadruple>, quadruple>(const int n,
					      const complex<quadruple> *x,
					      const int incx);
template
int blas_iamax<octruple, octruple>(const int n, const octruple *x, const int incx);
template
int blas_iamax<complex<octruple>, octruple>(const int n,
					      const complex<octruple> *x,
					      const int incx);
//
// BLAS 2
// dz gemv
#ifndef BLAS_GENERIC
template<>
void
blas_gemv<double>(const CBLAS_TRANSPOSE trA, 
		  const int m, const int n, const double &alpha, 
		  const double* A, const int lda, const double* x,
		  const int incx, 
		  const double &beta, double* y, const int incy)
{
#ifdef BLAS_FORTRAN
  unsigned char tra_ = 0;
  switch(trA) {
  case CblasNoTrans :
    tra_ = 'n';
    break;
  case CblasTrans :
    tra_ = 't';
    break;
  case CblasConjTrans:
    tra_ = 'c';
    break;
  }
  FORTRAN_DECL(dgemv)(tra_, m, n, alpha, A, lda, x, incx, beta, y, incy); //
#else
  cblas_dgemv(CblasColMajor, trA, (BLAS_INT)m, (BLAS_INT)n, alpha,
	      A, (BLAS_INT)lda, x, (BLAS_INT)incx, beta, y, (BLAS_INT)incy);
#endif
}

template<>
void
blas_gemv<complex<double> >(const CBLAS_TRANSPOSE trA, 
			    const int m, const int n,
			    const complex<double> &alpha, 
			    const complex<double>* A, const int lda, 
			    const complex<double>* x, const int incx, 
			    const complex<double> &beta,
			    complex<double>* y, const int incy)
{
#ifdef BLAS_FORTRAN
  unsigned char tra_ = 0;
  switch(trA) {
  case CblasNoTrans:
    tra_ = 'n';
    break;
  case CblasTrans:
    tra_ = 't';
    break;
  case CblasConjTrans:
    tra_ = 'c';
    break;
  }
  FORTRAN_DECL(zgemv)(tra_, m, n, (double *)&alpha, (double *)A, lda, 
	 (double *)x, incx, (double *)&beta, (double *)y, incy);
#else
  cblas_zgemv(CblasColMajor, trA, (BLAS_INT)m, (BLAS_INT)n, (BLAS_VOID *)&alpha,
	      (BLAS_VOID *)A, (BLAS_INT)lda, (BLAS_VOID *)x, (BLAS_INT)incx,
	      (BLAS_VOID *)&beta, (void *)y, (BLAS_INT)incy);
#endif
}

template<>
void
blas_gemv<float>(const CBLAS_TRANSPOSE trA, 
		 const int m, const int n, const float &alpha, 
		 const float* A, const int lda, const float* x,
		 const int incx, 
		 const float &beta, float* y, const int incy)
{
#ifdef BLAS_FORTRAN
  unsigned char tra_ = 0;
  switch(trA) {
  case CblasNoTrans :
    tra_ = 'n';
    break;
  case CblasTrans :
    tra_ = 't';
    break;
  case CblasConjTrans:
    tra_ = 'c';
    break;
  }
  FORTRAN_DECL(sgemv)(tra_, m, n, alpha, A, lda, x, incx, beta, y, incy); //
#else
  cblas_sgemv(CblasColMajor, trA, (BLAS_INT)m, (BLAS_INT)n, alpha,
	      A, (BLAS_INT)lda, x, (BLAS_INT)incx, beta, y, (BLAS_INT)incy);
#endif
}

template<>
void
blas_gemv<complex<float> >(const CBLAS_TRANSPOSE trA, 
			    const int m, const int n,
			    const complex<float> &alpha, 
			    const complex<float>* A, const int lda, 
			    const complex<float>* x, const int incx, 
			    const complex<float> &beta,
			    complex<float>* y, const int incy)
{
#ifdef BLAS_FORTRAN
  unsigned char tra_ = 0;
  switch(trA) {
  case CblasNoTrans:
    tra_ = 'n';
    break;
  case CblasTrans:
    tra_ = 't';
    break;
  case CblasConjTrans:
    tra_ = 'c';
    break;
  }
  FORTRAN_DECL(cgemv)(tra_, m, n, (float *)&alpha, (float *)A, lda, 
		      (float *)x, incx, (float *)&beta, (float *)y, incy);
#else
  cblas_cgemv(CblasColMajor, trA, (BLAS_INT)m, (BLAS_INT)n, (BLAS_VOID *)&alpha,
	      (BLAS_VOID *)A, (BLAS_INT)lda, (BLAS_VOID *)x, (BLAS_INT)incx,
	      (BLAS_VOID *)&beta, (void *)y, (BLAS_INT)incy);
#endif
}

#endif

template<typename T>
void
blas_gemv(const CBLAS_TRANSPOSE trA, 
	  const int m, const int n, const T &alpha, 
	  const T* A, const int lda, 
	  const T* x, const int incx, 
	  const T &beta, T* y, const int incy)
{
  const T zero(0.0);
  const T one(1.0);
  int lenx, leny,  ix, iy, jx, jy, kx, ky;
  bool noconj;

  if ((m == 0) || (n == 0) ||
      ((alpha == zero) && (beta == one))) {
    return;
  }

  noconj = (trA == CblasTrans);

  if (trA ==  CblasNoTrans) {
    lenx = n;
    leny = m;
  }
  else {
    lenx = m;
    leny = n;
  }
  if (incx > 0) {
    kx = 0;
  }
  else {
    kx = (lenx - 1) * (-incx); //
  }
  if (incy > 0) {
    ky = 0;
  }
  else {
    ky = (leny - 1) * (-incy); //
  }
  if (beta != one) {
    if (incy == 1) {
      if (beta == zero) {
	for (int i = 0; i < leny; i++) {
	  y[i] = zero;
	} // loop 10 : i
      }
      else {
	for (int i = 0; i < leny; i++) {
	  y[i] *= beta;
	} // loop 20 : i
      }
    }
    else {
      iy = ky;
      if (beta == zero) {
	for (int i = 0; i < leny; i++) {
	  y[iy] = zero;
	  iy += incy;
	} // loop 30 : i
      }
      else {
	for (int i = 0; i < leny; i++) {
	  y[iy] *= beta;
	  iy += incy;
	} // loop 40 : i
      }
    }
  } // if (beta != one)	  
  if (alpha  == zero) {
    return;
  }
  if (trA == CblasNoTrans) {
    //        Form  y := alpha*A*x + y.
    jx = kx;
    if (incy == 1) {
      for (int j = 0; j < n; j++) {
	if (x[jx] != zero) {
	  T temp = alpha * x[jx];
	  for (int i = 0; i < m; i++) {
	    y[i] += temp * A[i + j * lda];
	  } // loop 50 : i
	} // end if
	jx += incx;
      }     // loop 60 : j
    }
    else {
      for (int j = 0; j < n; j++) {
	if (x[jx] != zero) {
	  T temp = alpha * x[jx];
	  iy = ky;
	  for (int i = 0; i < m; i++) {
	    y[iy] += temp * A[i + j * lda];
	    iy += incy;
	  } // loop 70 : i
	}
	jx += incx;
      }     // loop 80 : j
    }
  }
  else {
    //    Form  y := alpha*A**T*x + y  or  y := alpha*conjg( A' )*x + y.
    jy = ky;
    if (incx == 1) {
      for (int j = 0; j < n; j++) {
	T temp = zero;
	if (noconj) {
	  for (int i = 0; i < m; i++) {
	    temp += A[i + j * lda] * x[i];
	  } // loop 90 : i
	}
	else {
	  for (int i = 0; i < m; i++) {
	    temp += blas_conj(A[i + j * lda]) * x[i];
	  } // loop 100 : j
	}
	y[jy] += alpha * temp;
	jy += incy;
      } // loop 110 : j
    }
    else {
      for (int j = 0; j < n; j++) {
	T temp = zero;
	ix = kx;
	if (noconj) {
	  for (int i = 0; i < m; i++) {
	    temp += A[i + j * lda] * x[ix];
	    ix += incx;
	  } // loop 120 : i
	}
	else {
	  for (int i = 0; i < m; i++) {
	    temp += blas_conj(A[i + j * lda]) * x[ix];
	    ix += incx;
	  } // loop 130 : i
	}
	y[jy] += alpha * temp;
	jy += incy;
      }  // loop 140 : j
    } // if (incx != 1)
  }
}

// explicit instantiation of blas_gemv
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_gemv<double>(const CBLAS_TRANSPOSE trA, 
		     const int m, const int n, const double &alpha, 
		     const double* A, const int lda,
		     const double* x, const int incx, 
		     const double &beta, double* y, const int incy);

template
void
blas_gemv<complex<double> >(const CBLAS_TRANSPOSE trA, 
			       const int m, const int n,
			       const complex<double> &alpha, 
			       const complex<double>* A, const int lda, 
			       const complex<double>* x, const int incx, 
			       const complex<double> &beta,
			       complex<double>* y, const int incy);

template
void
blas_gemv<float>(const CBLAS_TRANSPOSE trA, 
		     const int m, const int n, const float &alpha, 
		     const float* A, const int lda,
		     const float* x, const int incx, 
		     const float &beta, float* y, const int incy);

template
void
blas_gemv<complex<float> >(const CBLAS_TRANSPOSE trA, 
			       const int m, const int n,
			       const complex<float> &alpha, 
			       const complex<float>* A, const int lda, 
			       const complex<float>* x, const int incx, 
			       const complex<float> &beta,
			       complex<float>* y, const int incy);
#endif
template
void
blas_gemv<quadruple>(const CBLAS_TRANSPOSE trA, 
		     const int m, const int n, const quadruple &alpha, 
		     const quadruple* A, const int lda,
		     const quadruple* x, const int incx, 
		     const quadruple &beta, quadruple* y, const int incy);
#ifndef NO_OCTRUPLE
template
void
blas_gemv<octruple>(const CBLAS_TRANSPOSE trA, 
		     const int m, const int n, const octruple &alpha, 
		     const octruple* A, const int lda,
		     const octruple* x, const int incx, 
		     const octruple &beta, octruple* y, const int incy);
#endif
template
void
blas_gemv<complex<quadruple> >(const CBLAS_TRANSPOSE trA, 
			       const int m, const int n,
			       const complex<quadruple> &alpha, 
			       const complex<quadruple>* A, const int lda, 
			       const complex<quadruple>* x, const int incx, 
			       const complex<quadruple> &beta,
			       complex<quadruple>* y, const int incy);
#ifndef NO_OCTRUPLE
template
void
blas_gemv<complex<octruple> >(const CBLAS_TRANSPOSE trA, 
			      const int m, const int n,
			      const complex<octruple> &alpha, 
			      const complex<octruple>* A, const int lda, 
			      const complex<octruple>* x, const int incx, 
			      const complex<octruple> &beta,
			      complex<octruple>* y, const int incy);
#endif

// dz trsv
#ifndef BLAS_GENERIC
template<>
void
blas_trsv<double>(const CBLAS_UPLO Uplo,
		  const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
		  const int N, const double *A, const int lda, double *X,
		  const int incX) 
{
#ifdef BLAS_FORTRAN
  unsigned char uplo_ = 0;
  unsigned char transa_ = 0;
  unsigned char diag_ = 0;
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  switch(TransA) {
  case CblasNoTrans:
    transa_ = 'n';
    break;
  case CblasTrans:
    transa_ = 't';
    break;
  case CblasConjTrans:
    transa_ = 'c';
    break;
  }
  switch(Diag) {
  case CblasNonUnit:
    diag_ = 'n';
    break;
  case CblasUnit:
    diag_ = 'u';
    break;
  }
  FORTRAN_DECL(dtrsv)(uplo_, transa_, diag_, N, A, lda, X, incX);
#else
  cblas_dtrsv(CblasColMajor, Uplo, TransA, Diag, (BLAS_INT)N,
	      A, (BLAS_INT)lda, X, (BLAS_INT)incX);
#endif
}

template<>
void
blas_trsv<complex<double> >(const CBLAS_UPLO Uplo,
			    const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
			    const int N, const complex<double>* A,
			    const int lda, 
			    complex<double> *X, const int incX) 
{
#ifdef BLAS_FORTRAN
  unsigned char uplo_ = 0;
  unsigned char transa_ = 0;
  unsigned char diag_ = 0;
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  switch(TransA) {
  case CblasNoTrans:
    transa_ = 'n';
    break;
  case CblasTrans:
    transa_ = 't';
    break;
  case CblasConjTrans:
    transa_ = 'c';
    break;
  }
  switch(Diag) {
  case CblasNonUnit:
    diag_ = 'n';
    break;
  case CblasUnit:
    diag_ = 'u';
    break;
  }
  FORTRAN_DECL(ztrsv)(uplo_, transa_, diag_, N, (double *)A, lda,
		      (double *)X, incX);
#else
  cblas_ztrsv(CblasColMajor, Uplo, TransA, Diag, (BLAS_INT)N,
	      (BLAS_VOID *)A, (BLAS_INT)lda, (BLAS_VOID *)X, (BLAS_INT)incX);
#endif
}

template<>
void
blas_trsv<float>(const CBLAS_UPLO Uplo,
		  const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
		  const int N, const float *A, const int lda, float *X,
		  const int incX) 
{
#ifdef BLAS_FORTRAN
  unsigned char uplo_ = 0;
  unsigned char transa_ = 0;
  unsigned char diag_ = 0;
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  switch(TransA) {
  case CblasNoTrans:
    transa_ = 'n';
    break;
  case CblasTrans:
    transa_ = 't';
    break;
  case CblasConjTrans:
    transa_ = 'c';
    break;
  }
  switch(Diag) {
  case CblasNonUnit:
    diag_ = 'n';
    break;
  case CblasUnit:
    diag_ = 'u';
    break;
  }
  FORTRAN_DECL(strsv)(uplo_, transa_, diag_, N, A, lda, X, incX);
#else
  cblas_strsv(CblasColMajor, Uplo, TransA, Diag, (BLAS_INT)N,
	      A, (BLAS_INT)lda, X, (BLAS_INT)incX);
#endif
}

template<>
void
blas_trsv<complex<float> >(const CBLAS_UPLO Uplo,
			    const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
			    const int N, const complex<float>* A,
			    const int lda, 
			    complex<float> *X, const int incX) 
{
#ifdef BLAS_FORTRAN
  unsigned char uplo_ = 0;
  unsigned char transa_ = 0;
  unsigned char diag_ = 0;
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  switch(TransA) {
  case CblasNoTrans:
    transa_ = 'n';
    break;
  case CblasTrans:
    transa_ = 't';
    break;
  case CblasConjTrans:
    transa_ = 'c';
    break;
  }
  switch(Diag) {
  case CblasNonUnit:
    diag_ = 'n';
    break;
  case CblasUnit:
    diag_ = 'u';
    break;
  }
  FORTRAN_DECL(ctrsv)(uplo_, transa_, diag_, N, (float *)A, lda,
		      (float *)X, incX);
#else
  cblas_ctrsv(CblasColMajor, Uplo, TransA, Diag, (BLAS_INT)N,
	      (BLAS_VOID *)A, (BLAS_INT)lda, (BLAS_VOID *)X, (BLAS_INT)incX);
#endif
}
#endif

template<typename T>
void
blas_trsv(const CBLAS_UPLO Uplo,
	  const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
	  const int N, const T *A, const int lda, T *X,
	  const int incX)
{
  const T zero(0.0);
  bool noconj, nounit;
  int ix, jx, kx;
  if (N == 0) {
    return;
  }
  noconj = (TransA == CblasTrans);
  nounit = (Diag == CblasNonUnit);
  if (incX < 0) {
    kx = (N - 1) * (-incX); // decrement form the last index
  }
  else {
    kx = 0;
  }
  if (TransA == CblasNoTrans) {
    // Form  x := inv( A )*x.
    if (Uplo == CblasUpper) {
      if (incX == 1) {
	for (int j = (N - 1); j >= 0; j--) {
	  if (X[j] != zero) {
	    if (nounit) {
	      X[j] /= A[j + j * lda];
	    }
	    T temp = X[j];
	    for (int i = (j - 1); i >= 0; i--) {
	      X[i] -= temp * A[i + j * lda];
	    } // loop 10 : i
	  }
	} // loop 20 : j
      }
      else {
	jx = kx + (N - 1) * incX;
	for (int j = (N - 1); j >= 0; j--) {
	  if (X[jx] != zero) {
	    if (nounit) {
	      X[jx] /= A[j + j * lda];
	    }
	    T temp  = X[jx];
	    ix = jx;
	    for (int i = (j - 1); j >= 0; j--) {
	      ix -= incX;
	      X[ix] -= temp * A[i + j * lda];
	    }  // loop 30 : i
	  }
	  jx -= incX;
	} // loop 40 : j
      } // if (incX == 1)
    }
    else {
      if (incX == 1) {
	for (int j = 0; j < N; j++) {
	  if (X[j] != zero) {
	    if (nounit) {
	      X[j] /= A[j + j * lda];
	    }
	    T temp = X[j];
	    for (int i = j + 1; i < N; i++) {
	      X[i] -= temp * A[i + j * lda];
	    } // loop 50 : i
	  }
	} // loop 60 : j
      }
      else {
	jx = kx;
	for (int j = 0; j < N; j++) {
	  if (X[jx] != zero) {
	    if (nounit) {
	      X[jx] /= A[j + j * lda];
	    }
	    T temp = X[jx];
	    ix = jx;
	    for (int i = j + 1; i < N; i++) {
	      ix += incX;
	      X[ix] -= temp * A[i + j * lda];
	    } // loop 70 : i
	  }
	  jx += incX;
	} // loop 80 : j
      } // if (incX == 1)
    }   // if (Uplo == CblasUpper) 
  }
  else {
    //    Form  x := inv( A **T )*x or  x := inv( A **H ) *x.
    if (Uplo == CblasUpper) {
      if (incX == 1) {
	for (int j = 0; j < N; j++) {
	  T temp = X[j];
	  if (noconj) {
	    for (int i = 0; i < j; i++) {
	      temp -= A[i + j * lda] * X[i];
	    } // loop 90 : i
	    if (nounit) {
	      temp /= A[j + j * lda];
	    }
	  }
	  else {
	    for (int i = 0; i < j; i++) {
	      temp -= blas_conj(A[i + j * lda]) * X[i];
	    } // loop 100 : i
	    if (nounit) {
	      temp /= blas_conj(A[j + j * lda]);
	    }
	  }
	  X[j] = temp;
	} // loop 110 : j
      }
      else {
	jx = kx;
	for (int j = 0; j < N; j++) {
	  ix = kx;
	  T temp = X[jx];
	  if (noconj) {
	    for (int i = 0; i < j; i++) {
	      temp -= A[i + j * lda] * X[ix];
	      ix += incX;
	    } // loop 120 : i
	    if (nounit) {
	      temp /= A[j + j * lda];
	    }
	  }
	  else {
	    for (int i = 0; i < j; i++) {
	      temp -= blas_conj(A[i + j * lda]) * X[ix];
	      ix += incX;
	    } // loop 120 : i
	    if (nounit) {
	      temp /= blas_conj(A[j + j * lda]);
	    }
	  }
	  X[jx] = temp;
	  jx += incX;
	}
      }
    }
    else {
      if (incX == 1) {
	for (int j = (N - 1); j >= 0; j--) {
	  T temp = X[j];
	  if (noconj) {
	    for (int i = (N - 1); i > j; i--) {
	      temp -= A[i + j * lda] * X[i];
	    } // loop 150 : i
	    if (nounit) {
	      temp /= A[j + j * lda];
	    }
	  }
	  else {
	    for (int i = 0; i < j; i++) {
	      temp -= blas_conj(A[i + j * lda]) * X[i];
	    } // loop 160 : i
	    if (nounit) {
	      temp /= blas_conj(A[j + j * lda]);
	    }
	  }
	  X[j] = temp;
	} // loop 170 : j
      }
      else {
	kx = kx + N * incX;
	jx = kx;
	for (int j = (N - 1); j >= 0; j--) {
	  ix = kx;
	  T temp = X[jx];
	  if (noconj) {
	    for (int i = (N - 1); i > j; i--) {
	      temp -= A[i + j * lda] * X[ix];
	      ix -= incX;
	    } // loop 180 : i
	    if (nounit) {
	      temp /= A[j + j * lda];
	    }
	  }
	  else {
	    for (int i = (N - 1); i > j; i--) {
	      temp -= blas_conj(A[i + j * lda]) * X[ix];
	      ix -= incX;
	    } // loop 190 : i
	    if (nounit) {
	      temp /= blas_conj(A[j + j * lda]);
	    }
	  }
	  X[jx] = temp;
	  jx -= incX;
	} // loop 200 : j
      }
    }   // if (Uplo == CblasUpper) 
  }
}

// explicit instantiation of blas_trsv
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_trsv<double>(const CBLAS_UPLO Uplo,
		const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
		const int N, const double *A, const int lda,
		double *X, const int incX);

template
void
blas_trsv<complex<double> >(const CBLAS_UPLO Uplo,
			  const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
			  const int N,
			  const complex<double> *A, const int lda,
			  complex<double> *X, const int incX);

template
void
blas_trsv<float>(const CBLAS_UPLO Uplo,
		const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
		const int N, const float *A, const int lda,
		float *X, const int incX);

template
void
blas_trsv<complex<float> >(const CBLAS_UPLO Uplo,
			  const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
			  const int N,
			  const complex<float> *A, const int lda,
			  complex<float> *X, const int incX);
#endif
template
void
blas_trsv<quadruple>(const CBLAS_UPLO Uplo,
		const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
		const int N, const quadruple *A, const int lda,
		quadruple *X, const int incX);
#ifndef NO_OCTRUPLE
template
void
blas_trsv<octruple>(const CBLAS_UPLO Uplo,
		const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
		const int N, const octruple *A, const int lda,
		octruple *X, const int incX);
#endif
template
void
blas_trsv<complex<quadruple> >(const CBLAS_UPLO Uplo,
			  const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
			  const int N,
			  const complex<quadruple> *A, const int lda,
			  complex<quadruple> *X, const int incX);
#ifndef NO_OCTRUPLE
template
void
blas_trsv<complex<octruple> >(const CBLAS_UPLO Uplo,
			  const CBLAS_TRANSPOSE TransA, const CBLAS_DIAG Diag,
			  const int N,
			  const complex<octruple> *A, const int lda,
			  complex<octruple> *X, const int incX);
#endif

// dz syr : symmetric rank-1 update
template<typename T>
void
blas_syr(const CBLAS_UPLO Uplo,
	 const int N, const T &alpha,
	 const T *X, const int incX,
	 T *A, const int lda)
{
  const T zero(0.0);
  int ix, jx, kx;
  if ((N == 0) || (alpha == zero)) {
    return;
  }
  if (incX > 0) {
    kx = 0;
  }
  else {
    kx = (N - 1) * (-incX); // decrement form the last index
  }
  if (Uplo == CblasUpper) {
    //   Form  A  when A is stored in upper triangle.
    if (incX == 1) {
      for (int j = 0; j < N; j++) {
	if (X[j] != zero) {
	  T temp = alpha * X[j];
	  for (int i = 0; i <= j; i++) {
	    A[i + j * lda] += X[i] * temp;
	  } // loop 10 : i
	}
      } // loop 20 : j
    }
    else {
      jx = kx;
      for (int j = 0; j < N; j++) {
	if (X[jx] != zero) {
	  T temp = alpha * X[jx];
	  ix = kx;
	  for (int i = 0; i <= j; i++) {
	    A[i + j * lda] += X[ix] * temp;
	    ix += incX;
	  } // loop 30 : i
	}
	jx += incX;
      }     // loop 40 : j
    }
  }
  else {
    //       Form  A  when A is stored in lower triangle.
    if (incX == 1) {
      for (int j = 0; j < N; j++) {
	if (X[j] != zero) {
	  T temp = alpha * X[j];
	  for (int i = j; i < N; i++) {
	    A[i + j * lda] += X[i] * temp;
	  } // loop 50 : i
	}
      } // loop 60 : j
    }
    else {
      jx = kx;
      for (int j = 0; j < N; j++) {
	if (X[jx] != zero) {
	  T temp = alpha * X[jx];
	  ix = jx;
	  for (int i = j; i < N; i++) {
	    A[i + j * lda] += X[ix] * temp;
	    ix += incX;
	  } // loop 70 : i
	}
	jx += incX;
      }     // loop 80 : j
    }
  }
}
#ifndef BLAS_GENERIC
template<>
void
blas_syr<double>(const CBLAS_UPLO Uplo, 
		 const int N, const double &alpha,
		 const double *X, const int incX,
		 double *A, const int lda)
{
#ifdef BLAS_FORTRAN
  unsigned char uplo_ = 0;
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  FORTRAN_DECL(dsyr)(uplo_, N, alpha, X, incX, A, lda);
#else
  cblas_dsyr(CblasColMajor, Uplo, (BLAS_INT)N, alpha,
	     X, (BLAS_INT)incX, A, (BLAS_INT)lda);
#endif
}

template<>
void
blas_syr<complex<double> >(const CBLAS_UPLO Uplo, 
			   const int N, const complex<double> &alpha,
			   const complex<double> *X, const int incX,
			   complex<double> *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(zgeru)(N, N,
		      (double *)&alpha, (double *)X, incX, (double *)X, incX, 
		      (double *)A, lda); // no zsyr in BLAS! : Uplo is ignored
#else
  cblas_zgeru(CblasColMajor, N, N,
	      (BLAS_VOID *)&alpha, (BLAS_VOID *)X, (BLAS_INT)incX,
	      (BLAS_VOID *)X, (BLAS_INT)incX, 
	      (BLAS_VOID *)A, (BLAS_INT)lda); // no zsyr in BLAS! :  Uplo is ignored
#endif
}

template<>
void
blas_syr<float>(const CBLAS_UPLO Uplo, 
		 const int N, const float &alpha,
		 const float *X, const int incX,
		 float *A, const int lda)
{
#ifdef BLAS_FORTRAN
  unsigned char uplo_ = 0;
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  FORTRAN_DECL(ssyr)(uplo_, N, alpha, X, incX, A, lda);
#else
  cblas_ssyr(CblasColMajor, Uplo, (BLAS_INT)N, alpha,
	     X, (BLAS_INT)incX, A, (BLAS_INT)lda);
#endif
}

template<>
void
blas_syr<complex<float> >(const CBLAS_UPLO Uplo, 
			   const int N, const complex<float> &alpha,
			   const complex<float> *X, const int incX,
			   complex<float> *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(cgeru)(N, N,
		      (float *)&alpha, (float *)X, incX, (float *)X, incX, 
		      (float *)A, lda); // no zsyr in BLAS! : Uplo is ignored
#else
  cblas_cgeru(CblasColMajor, N, N,
	      (BLAS_VOID *)&alpha, (BLAS_VOID *)X, (BLAS_INT)incX,
	      (BLAS_VOID *)X, (BLAS_INT)incX, 
	      (BLAS_VOID *)A, (BLAS_INT)lda); // no zsyr in BLAS! :  Uplo is ignored
#endif
}
#endif

// explicit instantiation of blas_syr
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_syr<double>(const CBLAS_UPLO Uplo, const int N,
		 const double &alpha,
		 const double *X, const int incX,
		 double *A, const int lda);

template
void
blas_syr<complex<double> >(const CBLAS_UPLO Uplo, const int N,
			   const complex<double> &alpha,
			   const complex<double> *X, const int incX,
			   complex<double> *A, const int lda);

template
void
blas_syr<float>(const CBLAS_UPLO Uplo, const int N,
		 const float &alpha,
		 const float *X, const int incX,
		 float *A, const int lda);

template
void
blas_syr<complex<float> >(const CBLAS_UPLO Uplo, const int N,
			   const complex<float> &alpha,
			   const complex<float> *X, const int incX,
			   complex<float> *A, const int lda);
#endif
template
void
blas_syr<quadruple>(const CBLAS_UPLO Uplo, const int N,
	       const quadruple &alpha,
	       const quadruple *X, const int incX,
	       quadruple *A, const int lda);
#ifndef NO_OCTRUPLE
template
void
blas_syr<octruple>(const CBLAS_UPLO Uplo, const int N,
		   const octruple &alpha,
		   const octruple *X, const int incX,
		   octruple *A, const int lda);
#endif
template
void
blas_syr<complex<quadruple> >(const CBLAS_UPLO Uplo, const int N,
			 const complex<quadruple> &alpha,
			 const complex<quadruple> *X, const int incX,
			 complex<quadruple> *A, const int lda);
#ifndef NO_OCTRUPLE
template
void
blas_syr<complex<octruple> >(const CBLAS_UPLO Uplo, const int N,
			     const complex<octruple> &alpha,
			     const complex<octruple> *X, const int incX,
			     complex<octruple> *A, const int lda);
#endif

// dz syrc : with complex conjugate
template<typename T>
void
blas_syrc(const CBLAS_UPLO Uplo,
	  const int N, const T &alpha,
	  const T *X, const int incX,
	  T *A, const int lda)
{
  const T zero(0.0);
  int ix, jx, kx;
  if ((N == 0) || (alpha == zero)) {
    return;
  }
  if (incX > 0) {
    kx = 0;
  }
  else {
    kx = (N - 1) * (-incX); // decrement form the last index
  }
  if (Uplo == CblasUpper) {
    //   Form  A  when A is stored in upper triangle.
    if (incX == 1) {
      for (int j = 0; j < N; j++) {
	if (X[j] != zero) {
	  T temp = alpha * blas_conj(X[j]);
	  for (int i = 0; i <= j; i++) {
	    A[i + j * lda] += X[i] * temp;
	  } // loop 10 : i
	}
      } // loop 20 : j
    }
    else {
      jx = kx;
      for (int j = 0; j < N; j++) {
	if (X[jx] != zero) {
	  T temp = alpha * blas_conj(X[jx]);
	  ix = kx;
	  for (int i = 0; i <= j; i++) {
	    A[i + j * lda] += X[ix] * temp;
	    ix += incX;
	  } // loop 30 : i
	}
	jx += incX;
      }     // loop 40 : j
    }
  }
  else {
    //       Form  A  when A is stored in lower triangle.
    if (incX == 1) {
      for (int j = 0; j < N; j++) {
	if (X[j] != zero) {
	  T temp = alpha * blas_conj(X[j]);
	  for (int i = j; i < N; i++) {
	    A[i + j * lda] += X[i] * temp;
	  } // loop 50 : i
	}
      } // loop 60 : j
    }
    else {
      jx = kx;
      for (int j = 0; j < N; j++) {
	if (X[jx] != zero) {
	  T temp = alpha * blas_conj(X[jx]);
	  ix = jx;
	  for (int i = j; i < N; i++) {
	    A[i + j * lda] += X[ix] * temp;
	    ix += incX;
	  } // loop 70 : i
	}
	jx += incX;
      }     // loop 80 : j
    }
  }
}
// explicit instantiation of blas_syrc
template
void
blas_syrc<complex<double> >(const CBLAS_UPLO Uplo,
			    const int N, const complex<double> &alpha,
			    const complex<double> *X, const int incX,
			    complex<double> *A, const int lda);

template
void
blas_syrc<complex<float> >(const CBLAS_UPLO Uplo,
			    const int N, const complex<float> &alpha,
			    const complex<float> *X, const int incX,
			    complex<float> *A, const int lda);

template
void
blas_syrc<complex<quadruple> >(const CBLAS_UPLO Uplo,
			  const int N, const complex<quadruple> &alpha,
			  const complex<quadruple> *X, const int incX,
			  complex<quadruple> *A, const int lda);
#ifndef NO_OCTRUPLE
template
void
blas_syrc<complex<octruple> >(const CBLAS_UPLO Uplo,
			  const int N, const complex<octruple> &alpha,
			  const complex<octruple> *X, const int incX,
			  complex<octruple> *A, const int lda);
#endif

// dz syr2 symmetric rank-2 update : DSYR2 only exists in BLAS
#ifndef BLAS_GENERIC
template<>
void
blas_syr2<double>(const CBLAS_UPLO Uplo, const int N,
		  const double &alpha,
		  const double *X, const int incX,
		  const double *Y, const int incY,
		  double *A, const int lda)
{
#ifdef BLAS_FORTRAN
  unsigned char uplo_ = 0;
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  FORTRAN_DECL(dsyr2)(uplo_, N, alpha, X, incX, Y, incY, A, lda);
#else
  cblas_dsyr2(CblasColMajor, Uplo, (BLAS_INT)N, alpha,
	      X, (BLAS_INT)incX, Y, (BLAS_INT)incY, A, (BLAS_INT)lda);
#endif
}

template<>
void
blas_syr2<float>(const CBLAS_UPLO Uplo, const int N,
		  const float &alpha,
		  const float *X, const int incX,
		  const float *Y, const int incY,
		  float *A, const int lda)
{
#ifdef BLAS_FORTRAN
  unsigned char uplo_ = 0;
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  FORTRAN_DECL(ssyr2)(uplo_, N, alpha, X, incX, Y, incY, A, lda);
#else
  cblas_ssyr2(CblasColMajor, Uplo, (BLAS_INT)N, alpha,
	      X, (BLAS_INT)incX, Y, (BLAS_INT)incY, A, (BLAS_INT)lda);
#endif
}
#endif

template<typename T>
void
blas_syr2(const CBLAS_UPLO Uplo,
	  const int N, const T &alpha,
	  const T *X, const int incX,
	  const T *Y, const int incY,
	  T *A, const int lda)
{
  const T zero(0.0);
  int ix, iy, jx, jy, kx, ky;
  if ((N == 0) || (alpha == zero)) {
    return;
  }
  if ((incX != 1) || (incY != 1)) {
    if (incX > 0) {
      kx = 0;
    }
    else {
      kx = (N - 1) * (-incX); // decrement form the last index
    }
    if (incY > 0) {
      ky = 0;
    }
    else {
      ky = (N - 1) * (-incY); // decrement form the last index
    }
    jx = kx;
    jy = ky;
  }
  if (Uplo == CblasUpper) {
    //   Form  A  when A is stored in upper triangle.
    if ((incX == 1) && (incY == 1)) {
      for (int j = 0; j < N; j++) {
	if ((X[j] != zero) || (Y[j] != zero)) {
	  T temp1 = alpha * Y[j];
	  T temp2 = alpha * X[j];
	  for (int i = 0; i < j; i++) {
	    A[i + j * lda] += X[i] * temp1 + Y[i] * temp2;
	  } // loop 10 : i
	}
      } // loop 20 : j
    }
    else {
      for (int j = 0; j < N; j++) {
	if ((X[jx] != zero) || (Y[jy] != zero))  {
	  T temp1 = alpha * Y[jy];
	  T temp2 = alpha * X[jx];
	  ix = kx;
	  iy = ky; // 1st Nov.2015 bug fund : ky = ky;
	  for (int i = 0; i < j; i++) {
	    A[i + j * lda] += X[ix] * temp1 + Y[iy] * temp2;
	    ix += incX;
	    iy += incY;
	  } // loop 30 : i
	}
	jx += incX;
	jx += incY;
      }     // loop 40 : j
    }
  }
  else {
    //       Form  A  when A is stored in lower triangle.
    if ((incX == 1) && (incY == 1)) {
      for (int j = 0; j < N; j++) {
	if ((X[j] != zero) || (Y[j] != zero)) {
	  T temp1 = alpha * Y[j];
	  T temp2 = alpha * X[j];
	  for (int i = j; i < N; i++) {
	    A[i + j * lda] += X[i] * temp1 + Y[i] * temp2;
	  } // loop 50 : i
	}
      } // loop 60 : j
    }
    else {
      for (int j = 0; j < N; j++) {
	if ((X[jx] != zero) || (Y[jy] != zero)) {
	  T temp1 = alpha * Y[jy];
	  T temp2 = alpha * X[jx];
	  ix = jx;
	  iy = jy;
	  for (int i = j; i < N; i++) {
	    A[i + j * lda] += X[ix] * temp1 + Y[iy] * temp2;
	    ix += incX;
	    iy += incY;
	  } // loop 70 : i
	}
	jx += incX;
	jy += incY;
      }     // loop 80 : j
    }
  }
}
// explicit instantiation of blas_syr2
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_syr2<double>(const CBLAS_UPLO Uplo, const int N,
		  const double &alpha,
		  const double *X, const int incX,
		  const double *Y, const int incY,
		  double *A, const int lda);

template
void
blas_syr2<float>(const CBLAS_UPLO Uplo, const int N,
		  const float &alpha,
		  const float *X, const int incX,
		  const float *Y, const int incY,
		  float *A, const int lda);
#endif
template
void
blas_syr2<complex<double> >(const CBLAS_UPLO Uplo, const int N,
			    const complex<double> &alpha,
			    const complex<double> *X, const int incX,
			    const complex<double> *Y, const int incY,
			    complex<double> *A, const int lda);
// there is no cblas_syr2

template
void
blas_syr2<quadruple>(const CBLAS_UPLO Uplo, const int N,
		     const quadruple &nalpha,
		     const quadruple *X, const int incX,
		     const quadruple *Y, const int incY,
		     quadruple *A, const int lda);

template
void
blas_syr2<complex<quadruple> >(const CBLAS_UPLO Uplo, const int N,
			       const complex<quadruple> &alpha,
			       const complex<quadruple> *X, const int incX,
			       const complex<quadruple> *Y, const int incY,
			       complex<quadruple> *A, const int lda);

// dz ger
#ifndef BLAS_GENERIC
template<>
void
blas_ger<double>(const int M, const int N, const double &alpha,
		 const double *X, const int incX,
		 const double *Y, const int incY,
		 double *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(dger)(M, N, alpha, X, incX, Y, incY, A, lda);
#else
  cblas_dger(CblasColMajor, (BLAS_INT)M, (BLAS_INT)N, alpha,
	     X, (BLAS_INT)incX, Y, (BLAS_INT)incY, A, (BLAS_INT)lda);
#endif
}

template<>
void
blas_ger<complex<double> >(const int M, const int N,
			   const complex<double> &alpha,
			   const complex<double> *X, const int incX,
			   const complex<double> *Y, const int incY,
			   complex<double> *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(zgeru)(M, N, (double *)&alpha,
	 (double *)X, incX, (double *)Y, incY, 
	 (double *)A, lda);
#else
    cblas_zgeru(CblasColMajor, (BLAS_INT)M, (BLAS_INT)N, (BLAS_VOID *)&alpha,
		(BLAS_VOID *)X, (BLAS_INT)incX, (BLAS_VOID *)Y, (BLAS_INT)incY, 
		(BLAS_VOID *)A, (BLAS_INT)lda); 
  #endif
}

template<>
void
blas_ger<float>(const int M, const int N, const float &alpha,
		 const float *X, const int incX,
		 const float *Y, const int incY,
		 float *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(sger)(M, N, alpha, X, incX, Y, incY, A, lda);
#else
  cblas_sger(CblasColMajor, (BLAS_INT)M, (BLAS_INT)N, alpha,
	     X, (BLAS_INT)incX, Y, (BLAS_INT)incY, A, (BLAS_INT)lda);
#endif
}

template<>
void
blas_ger<complex<float> >(const int M, const int N,
			   const complex<float> &alpha,
			   const complex<float> *X, const int incX,
			   const complex<float> *Y, const int incY,
			   complex<float> *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(cgeru)(M, N, (float *)&alpha,
	 (float *)X, incX, (float *)Y, incY, 
	 (float *)A, lda);
#else
    cblas_cgeru(CblasColMajor, (BLAS_INT)M, (BLAS_INT)N, (BLAS_VOID *)&alpha,
		(BLAS_VOID *)X, (BLAS_INT)incX, (BLAS_VOID *)Y, (BLAS_INT)incY, 
		(BLAS_VOID *)A, (BLAS_INT)lda); 
  #endif
}
#endif

template<typename T>
void
blas_ger(const int M, const int N, const T &alpha,
	 const T *X, const int incX,
	 const T *Y, const int incY,
	 T *A, const int lda)
{
  const T zero(0.0);
	       
  int ix, jy, kx;
  if ((M == 0) || (N == 0) || (alpha == zero)) {
    return;
  }
  if (incY > 0) {
    jy = 0;
  }
  else {
    jy = (N - 1) * (-incY); // decrement form the last index
  }
  if (incX == 1) {
    for (int j = 0; j < N; j++) {
      if (Y[jy] != zero) {
	T temp = alpha * Y[jy];
	for (int i = 0; i < M; i++) {
	  A[i + j * lda] += X[i] * temp;
	} // loop 10 : i
      }
      jy += incY;
    } // loop 20 : j
  }
  else {
    if (incX > 0) {
      kx = 0;
    }
    else {
      kx = (M - 1) * (-incX); // decrement form the last index
    }
    for (int j = 0; j < N; j++) {
      if (Y[jy] != zero) {
	T temp = alpha * Y[jy];
	ix = kx;
	for (int i = 0; i < M; i++) {
	  A[i + j * lda] += X[ix] * temp;
	  ix += incX;
	} // loop 30 : i
      }
      jy += incY;
    }     // loop 40 : j
  }
}

// explicit instantiation of blas_ger
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_ger<double>(const int M, const int N, const double &alpha,
		 const double *X, const int incX,
		 const double *Y, const int incY,
		 double *A, const int lda);

template
void
blas_ger<complex<double> >(const int M, const int N,
			   const complex<double> &alpha,
			   const complex<double> *X, const int incX,
			   const complex<double> *Y, const int incY,
			   complex<double> *A, const int lda);

template
void
blas_ger<float>(const int M, const int N, const float &alpha,
		 const float *X, const int incX,
		 const float *Y, const int incY,
		 float *A, const int lda);

template
void
blas_ger<complex<float> >(const int M, const int N,
			   const complex<float> &alpha,
			   const complex<float> *X, const int incX,
			   const complex<float> *Y, const int incY,
			   complex<float> *A, const int lda);
#endif
template
void
blas_ger<quadruple>(const int M, const int N, const quadruple &alpha,
		    const quadruple *X, const int incX,
		    const quadruple *Y, const int incY,
		    quadruple *A, const int lda);
template
void
blas_ger<complex<quadruple> >(const int M, const int N,
			      const complex<quadruple> &alpha,
			      const complex<quadruple> *X, const int incX,
			      const complex<quadruple> *Y, const int incY,
			      complex<quadruple> *A, const int lda);
#ifndef NO_OCTRUPLE
template
void
blas_ger<octruple>(const int M, const int N, const octruple &alpha,
		   const octruple *X, const int incX,
		   const octruple *Y, const int incY,
		   octruple *A, const int lda);
template
void
blas_ger<complex<octruple> >(const int M, const int N,
			     const complex<octruple> &alpha,
			     const complex<octruple> *X, const int incX,
			     const complex<octruple> *Y, const int incY,
			     complex<octruple> *A, const int lda);
#endif
// dz gerc : rank-1 update with complex conjugate
#ifndef BLAS_GENERIC
template<>
void
blas_gerc<double>(const int M, const int N, const double &alpha,
		  const double *X, const int incX,
		  const double *Y, const int incY,
		  double *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(dger)(M, N, alpha, X, incX, Y, incY, A, lda);
#else
  cblas_dger(CblasColMajor, (BLAS_INT)M, (BLAS_INT)N,
	     alpha, X, (BLAS_INT)incX, Y, (BLAS_INT)incY, A, (BLAS_INT)lda);
#endif
}

template<>
void
blas_gerc<complex<double> >(const int M, const int N,
			    const complex<double> &alpha,
			    const complex<double> *X, const int incX,
			    const complex<double> *Y, const int incY,
			    complex<double> *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(zgerc)(M, N, (double *)&alpha,
		      (double *)X, incX, (double *)Y, incY, 
		      (double *)A, lda); // no zsyr in BLAS!
#else
  cblas_zgerc(CblasColMajor, (BLAS_INT)M, (BLAS_INT)N, (BLAS_VOID *)&alpha,
	      (BLAS_VOID *)X, (BLAS_INT)incX, (void *)Y, (BLAS_INT)incY, 
	      (BLAS_VOID *)A, (BLAS_INT)lda); // no zsyr in BLAS!
#endif
}

template<>
void
blas_gerc<float>(const int M, const int N, const float &alpha,
		  const float *X, const int incX,
		  const float *Y, const int incY,
		  float *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(sger)(M, N, alpha, X, incX, Y, incY, A, lda);
#else
  cblas_sger(CblasColMajor, (BLAS_INT)M, (BLAS_INT)N,
	     alpha, X, (BLAS_INT)incX, Y, (BLAS_INT)incY, A, (BLAS_INT)lda);
#endif
}

template<>
void
blas_gerc<complex<float> >(const int M, const int N,
			    const complex<float> &alpha,
			    const complex<float> *X, const int incX,
			    const complex<float> *Y, const int incY,
			    complex<float> *A, const int lda)
{
#ifdef BLAS_FORTRAN
  FORTRAN_DECL(cgerc)(M, N, (float *)&alpha,
		      (float *)X, incX, (float *)Y, incY, 
		      (float *)A, lda); // no zsyr in BLAS!
#else
  cblas_cgerc(CblasColMajor, (BLAS_INT)M, (BLAS_INT)N, (BLAS_VOID *)&alpha,
	      (BLAS_VOID *)X, (BLAS_INT)incX, (void *)Y, (BLAS_INT)incY, 
	      (BLAS_VOID *)A, (BLAS_INT)lda); // no zsyr in BLAS!
#endif
}

#endif

template<typename T>
void
blas_gerc(const int M, const int N, const T &alpha,
	  const T *X, const int incX,
	  const T *Y, const int incY,
	  T *A, const int lda)
{
  const T zero(0.0);
  int ix, jy, kx;
  if ((M == 0) || (N == 0) || (alpha == zero)) {
    return;
  }
  if (incY > 0) {
    jy = 0;
  }
  else {
    jy = (N - 1) * (-incY); // decrement form the last index
  }
  if (incX == 1) {
    for (int j = 0; j < N; j++) {
      if (Y[jy] != zero) {
	T temp = alpha * blas_conj(Y[jy]);
	for (int i = 0; i < M; i++) {
	  A[i + j * lda] += X[i] * temp;
	} // loop 10 : i
      }
      jy += incY;
    } // loop 20 : j
  }
  else {
    if (incX > 0) {
      kx = 0;
    }
    else {
      kx = (M - 1) * (-incX); // decrement form the last index
    }
    for (int j = 0; j < N; j++) {
      if (Y[jy] != zero) {
	T temp = alpha * blas_conj(Y[jy]);
	ix = kx;
	for (int i = 0; i < M; i++) {
	  A[i + j * lda] += X[ix] * temp;
	  ix += incX;
	} // loop 30 : i
      }
      jy += incY;
    }     // loop 40 : j
  }
}

// explicit instantiation of blas_gerc
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_gerc<double>(const int M, const int N, const double &alpha,
		  const double *X, const int incX,
		  const double *Y, const int incY,
		  double *A, const int lda);

template
void
blas_gerc<complex<double> >(const int M, const int N,
			    const complex<double> &alpha,
			    const complex<double> *X, const int incX,
			    const complex<double> *Y, const int incY,
			    complex<double> *A, const int lda);

template
void
blas_gerc<float>(const int M, const int N, const float &alpha,
		  const float *X, const int incX,
		  const float *Y, const int incY,
		  float *A, const int lda);

template
void
blas_gerc<complex<float> >(const int M, const int N,
			    const complex<float> &alpha,
			    const complex<float> *X, const int incX,
			    const complex<float> *Y, const int incY,
			    complex<float> *A, const int lda);
#endif
template
void
blas_gerc<quadruple>(const int M, const int N, const quadruple &alpha,
		     const quadruple *X, const int incX,
		     const quadruple *Y, const int incY,
		     quadruple *A, const int lda);

template
void
blas_gerc<complex<quadruple> >(const int M, const int N,
			       const complex<quadruple> &alpha,
			       const complex<quadruple> *X, const int incX,
			       const complex<quadruple> *Y, const int incY,
			       complex<quadruple> *A, const int lda);

template
void
blas_gerc<octruple>(const int M, const int N, const octruple &alpha,
		     const octruple *X, const int incX,
		     const octruple *Y, const int incY,
		     octruple *A, const int lda);

template
void
blas_gerc<complex<octruple> >(const int M, const int N,
			       const complex<octruple> &alpha,
			       const complex<octruple> *X, const int incX,
			       const complex<octruple> *Y, const int incY,
			       complex<octruple> *A, const int lda);

// BLAS 3

// dz trsm
#ifndef BLAS_GENERIC
template<>
void
blas_trsm<double>(const CBLAS_SIDE Side,
		  const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
		  const CBLAS_DIAG Diag, const int M, const int N,
		  const double &alpha, const double *A, const int lda,
		  double *B, const int ldb)
{
#ifdef BLAS_FORTRAN
  unsigned char side_ = 0;
  unsigned char uplo_ = 0;
  unsigned char transa_ = 0;
  unsigned char diag_ = 0;
  switch(Side) {
  case CblasLeft:
    side_ = 'l';
    break;
  case CblasRight:
    side_ = 'r';
    break;
  }
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  switch(TransA) {
  case CblasNoTrans:
    transa_ = 'n';
    break;
  case CblasTrans:
    transa_ = 't';
    break;
  case CblasConjTrans:
    transa_ = 'c';
    break;
  }
  switch(Diag) {
  case CblasNonUnit:
    diag_ = 'n';
    break;
  case CblasUnit:
    diag_ = 'u';
    break;
  }
  FORTRAN_DECL(dtrsm)(side_, uplo_, transa_,
		      diag_, M, N,
		      alpha, A, lda,
		      B, ldb);
#else
  cblas_dtrsm(CblasColMajor, Side, Uplo, TransA,
	      Diag, (BLAS_INT)M, (BLAS_INT)N,
              alpha, A, (BLAS_INT)lda,
              B, (BLAS_INT)ldb);
#endif
}

template<>
void
blas_trsm<complex<double> >(const CBLAS_SIDE Side,
			    const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
			    const CBLAS_DIAG Diag, const int M, const int N,
			    const complex<double> &alpha,
			    const complex<double> *A, const int lda,
			    complex<double> *B, const int ldb)
{
#ifdef BLAS_FORTRAN
  unsigned char side_ = 0;
  unsigned char uplo_ = 0;
  unsigned char transa_ = 0;
  unsigned char diag_ = 0;
  switch(Side) {
  case CblasLeft:
    side_ = 'l';
    break;
  case CblasRight:
    side_ = 'r';
    break;
  }
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  switch(TransA) {
  case CblasNoTrans:
    transa_ = 'n';
    break;
  case CblasTrans:
    transa_ = 't';
    break;
  case CblasConjTrans:
    transa_ = 'c';
    break;
  }
  switch(Diag) {
  case CblasNonUnit:
    diag_ = 'u';
    break;
  case CblasUnit:
    diag_ = 'n';
    break;
  }
  FORTRAN_DECL(ztrsm)(side_, uplo_, transa_,
		      diag_, M, N,
		      (double *)&alpha, (double *)A, lda,
		      (double *)B, ldb);
  
#else
  cblas_ztrsm(CblasColMajor, Side, Uplo, TransA,
	      Diag, (BLAS_INT)M, (BLAS_INT)N,
              (BLAS_VOID *)&alpha, (BLAS_VOID *)A, (BLAS_INT)lda,
              (BLAS_VOID *)B, (BLAS_INT)ldb);
#endif
}

template<>
void
blas_trsm<float>(const CBLAS_SIDE Side,
		  const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
		  const CBLAS_DIAG Diag, const int M, const int N,
		  const float &alpha, const float *A, const int lda,
		  float *B, const int ldb)
{
#ifdef BLAS_FORTRAN
  unsigned char side_ = 0;
  unsigned char uplo_ = 0;
  unsigned char transa_ = 0;
  unsigned char diag_ = 0;
  switch(Side) {
  case CblasLeft:
    side_ = 'l';
    break;
  case CblasRight:
    side_ = 'r';
    break;
  }
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  switch(TransA) {
  case CblasNoTrans:
    transa_ = 'n';
    break;
  case CblasTrans:
    transa_ = 't';
    break;
  case CblasConjTrans:
    transa_ = 'c';
    break;
  }
  switch(Diag) {
  case CblasNonUnit:
    diag_ = 'n';
    break;
  case CblasUnit:
    diag_ = 'u';
    break;
  }
  FORTRAN_DECL(strsm)(side_, uplo_, transa_,
		      diag_, M, N,
		      alpha, A, lda,
		      B, ldb);
#else
  cblas_strsm(CblasColMajor, Side, Uplo, TransA,
	      Diag, (BLAS_INT)M, (BLAS_INT)N,
              alpha, A, (BLAS_INT)lda,
              B, (BLAS_INT)ldb);
#endif
}

template<>
void
blas_trsm<complex<float> >(const CBLAS_SIDE Side,
			    const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
			    const CBLAS_DIAG Diag, const int M, const int N,
			    const complex<float> &alpha,
			    const complex<float> *A, const int lda,
			    complex<float> *B, const int ldb)
{
#ifdef BLAS_FORTRAN
  unsigned char side_ = 0;
  unsigned char uplo_ = 0;
  unsigned char transa_ = 0;
  unsigned char diag_ = 0;
  switch(Side) {
  case CblasLeft:
    side_ = 'l';
    break;
  case CblasRight:
    side_ = 'r';
    break;
  }
  switch(Uplo) {
  case CblasUpper:
    uplo_ = 'u';
    break;
  case CblasLower:
    uplo_ = 'l';
    break;
  }
  switch(TransA) {
  case CblasNoTrans:
    transa_ = 'n';
    break;
  case CblasTrans:
    transa_ = 't';
    break;
  case CblasConjTrans:
    transa_ = 'c';
    break;
  }
  switch(Diag) {
  case CblasNonUnit:
    diag_ = 'u';
    break;
  case CblasUnit:
    diag_ = 'n';
    break;
  }
  FORTRAN_DECL(cztrsm)(side_, uplo_, transa_,
		      diag_, M, N,
		      (float *)&alpha, (float *)A, lda,
		      (float *)B, ldb);
  
#else
  cblas_ctrsm(CblasColMajor, Side, Uplo, TransA,
	      Diag, (BLAS_INT)M, (BLAS_INT)N,
              (BLAS_VOID *)&alpha, (BLAS_VOID *)A, (BLAS_INT)lda,
              (BLAS_VOID *)B, (BLAS_INT)ldb);
#endif
}
#endif

template<typename T>
void
blas_trsm(const CBLAS_SIDE Side,
	  const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
	  const CBLAS_DIAG Diag, const int M, const int N,
	  const T &alpha,
	  const T *A, const int lda,
	  T *B, const int ldb)
{
  const T zero(0.0);
  const T one(1.0);
  bool lside, noconj, nounit, upper;
  lside = (Side == CblasLeft);
  noconj = (TransA == CblasTrans);
  nounit = (Diag == CblasNonUnit);
  upper = (Uplo == CblasUpper);

  if ((M == 0) || (N == 0)) {
    return;
  }
  if (alpha == zero) {
    for (int j = 0; j < N; j++) {
      for (int i = 0; i < M; i++) {
	B[i + j * ldb] = zero;
      }// loop 10 : i
    }  // loop 10 : j
  }
  if (lside) {
    if (TransA == CblasNoTrans) {
      //      Form  B := alpha*inv( A )*B.
      if (upper) {
	for (int j = 0; j < N; j++) {
	  if (alpha != one) {
	    for (int i = 0; i < M; i++) {
	      B[i + j * ldb] *= alpha;
	    } // loop 30 : i
	  }
	  for (int k = (M - 1); k >= 0; k--) {
	    if (B[k + j * ldb] != zero) {
	      if (nounit) {
		B[k + j * ldb] /= A[k + k * lda];
	      }
	      for (int i = 0; i < k; i++) {
		B[i + j * ldb] -= B[k + j *ldb] * A[i + k * lda];
	      } // loop 40 : i
	    }
	  } // loop 50 : k
	}   // loop 60 : j
      } //  if (upper)
      else {
	for (int j = 0; j < N; j++) {
	  if (alpha != one) {
	    for (int i = 0; i < M; i++) {
	      B[i + j * ldb] *= alpha;
	    } // loop 70 : i
	  }
	  for (int k = 0; k < M; k++) {
	    if (B[k + j * ldb] != zero) {
	      if (nounit) {
		B[k + j * ldb] /= A[k + k * lda];
	      }
	      for (int i = (k + 1); i < M; i++) {
		B[i + j * ldb] -= B[k + j * ldb] * A[i + k * lda];
	      } // loop 80 : i
	    }
	  } //loop 90 : k
	} // loop 100 :  j
      }  // if (upper)
    }  //     if (TransA == CblasNoTrans)
    else {
      //    Form  B := alpha*inv( A**T )*B
      //    or    B := alpha*inv( A**H )*B.
      if (upper) {
	for (int j = 0; j < N; j++) {
	  for (int i = 0; i < M; i++) {
	    T temp = alpha * B[i + j * ldb];
	    if (noconj) {
	      for (int k = 0; k < i; k++) {
		temp -= A[k + i * lda] * B[k + j * ldb];
	      } // loop 110 : k
	      if (nounit) {
		temp /= A[i + i * lda];
	      }
	    }
	    else {
	      for (int k = 0; k < i; k++) {
		temp -= blas_conj(A[k + i * lda]) * B[k + j * ldb];
	      } // loop 120 : k
	      if (nounit) {
		temp /= blas_conj(A[i + i * lda]);
	      }
	    }
	    B[i + j * ldb] = temp;
	  } // loop 130 : i
	}   // loop 140 : j
      } // if (upper)
      else {
	for (int j = 0; j < N; j++) {
	  for (int i = (M - 1); i >= 0; i--) {
	    T temp = alpha * B[i + j * ldb];
	    if (noconj) {
	      for (int k = (i + 1); k < M; k++) {
		temp -= A[k + i * lda] * B[k + j * ldb];
	      } // loop 150 : k
	      if (nounit) {
		temp /= A[i + i * lda];
	      }
	    }
	    else {
	      for (int k = (i + 1); k < M; k++) {
		temp -= blas_conj(A[k + i * lda]) * B[k + j * ldb];
	      } // loop 160 : k
	      if (nounit) {
		temp /= blas_conj(A[i + i * lda]);
	      }
	    }
	    B[i + j * ldb] = temp;
	  } // loop 170 : i
	}   // loop 180 : j
      } // if (upper)
    }   //     if (TransA == CblasNoTrans)
  }     // if (lsdie)
  else {
    if (TransA == CblasNoTrans) {
//           Form  B := alpha*B*inv( A ).
      if (upper) {
	for (int j = 0; j < N; j++) {
	  if (alpha != one) {
	    for (int i = 0; i < M; i++) {
	      B[i + j * ldb] *= alpha;
	    } // loop 190 : i
	  }
	  for (int k = 0; k < j; k++) {
	      if (A[k + j * lda] != zero) {
		for (int i = 0; i < M; i++) {
		  B[i + j * ldb] -= A[k + j * lda] * B[i + k * ldb];
		} // loop 200 : i
	      }
	    }  // loop 210 : k
	    if (nounit) {
	      T temp = one / A[j + j * lda];
	      for (int i = 0; i < M; i++) {
		B[i + j * ldb] *= temp;
	      } // loop 220 : i
	    }
	  } // loop 230 : j
      }
      else {
	for (int j = (N - 1); j >= 0; j--) {
	  if (alpha != one) {
	    for (int i = 0; i < M; i++) {
	      B[i + j * ldb] *= alpha;
	    } // loop 240 : i
	  }
	  for (int k = (j + 1); k < N; k++) {
	    if (A[k + j * lda] != zero) {
	      for (int i = 0; i < M; i++) {
		B[i + j * ldb] -= A[k + j * lda] * B[i + k * ldb];
	      } // loop 2500 : i
	    }
	  }  // loop 260 : k
	  if (nounit) {
	    T temp = one / A[j + j * lda];
	    for (int i = 0; i < M; i++) {
	      B[i + j * ldb] *= temp;
	    } // loop 270 : i
	  }
	} // loop 280 : j
      }
    } // if (TransA == CblasNoTrans) {
    else { 
      //           Form  B := alpha*B*inv( A**T )
      //           or    B := alpha*B*inv( A**H ).
      if (upper) {
	for (int k = (N - 1); k >= 0; k--) {
	  if (nounit) {
	    T temp;
	    if (noconj) {
	      temp = one / A[k + k * lda];
	    }
	    else {
	      temp = one / blas_conj(A[k + k * lda]);
	    }
	    for (int i = 0; i < M; i++) {
	      B[i + k * ldb] *= temp;
	    } // loop 290 : i
	  }
	  for (int j = 0; j < k; j++) {
	    if (A[j + k * lda] != zero) {
	      T temp;
	      if (noconj) {
		temp = A[j + k * lda];
	      }
	      else {
		temp = blas_conj(A[j + k * lda]);
	      }
	      for (int i = 0; i < M; i++) {
		B[i + j * ldb] -= temp * B[i + k * ldb];
	      } //loop 300 : i
	    }
	  }  // loop 310 : j
	  if (alpha != one) {
	    for (int i = 0; i < M; i++) {
	      B[i + k + ldb] *= alpha;
	    } // loop 320 : i
	  }
	}     // loop 330 : k
      }       // if (upper)
      else {
	for (int k = 0; k < N; k++) {
	  if (nounit) {
	    T temp;
	    if (noconj) {
	      temp = one / A[k + k * lda];
	      }
	    else {
	      temp = one / blas_conj(A[k + k * lda]);
	    }
	    for (int i = 0; i < M; i++) {
	      B[i + k * ldb] *= temp;
	    } // loop 340 : i
	  } // if (nounit)
	  for (int j = (k + 1); j < N; j++) {
	    if (A[j + k * lda] != zero) {
	      T temp;
	      if (noconj) {
		temp = A[j + k * lda];
	      }
	      else {
		temp = blas_conj(A[j + k * lda]);
	      }
	      for (int i = 0; i < M; i++) {
		B[i + j * ldb] -= temp * B[i + k * ldb];
	      } //loop 350 : i
	    }
	  }  // loop 360 : j
	  if (alpha != one) {
	    for (int i = 0; i < M; i++) {
	      B[i + k + ldb] *= alpha;
	    } // loop 370 : i
	  }
	}     // loop 380 : k
      }
    }
  } // if (lside)
}
//#endif
// explicit instantiation of blas_trsm
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_trsm<double>(const CBLAS_SIDE Side,
		  const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
		  const CBLAS_DIAG Diag, const int M, const int N,
		  const double &alpha,
		  const double *A, const int lda,
		  double *B, const int ldb);

template
void
blas_trsm<complex<double> >(const CBLAS_SIDE Side,
			    const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
			    const CBLAS_DIAG Diag, const int M, const int N,
			    const complex<double> &alpha,
			    const complex<double> *A, const int lda,
			    complex<double> *B, const int ldb);

template
void
blas_trsm<float>(const CBLAS_SIDE Side,
		  const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
		  const CBLAS_DIAG Diag, const int M, const int N,
		  const float &alpha,
		  const float *A, const int lda,
		  float *B, const int ldb);

template
void
blas_trsm<complex<float> >(const CBLAS_SIDE Side,
			    const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
			    const CBLAS_DIAG Diag, const int M, const int N,
			    const complex<float> &alpha,
			    const complex<float> *A, const int lda,
			    complex<float> *B, const int ldb);
#endif
template
void
blas_trsm<quadruple>(const CBLAS_SIDE Side,
		     const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
		     const CBLAS_DIAG Diag, const int M, const int N,
		     const quadruple &alpha,
		     const quadruple *A, const int lda,
		     quadruple *B, const int ldb);
#ifndef NO_OCTRUPLE
template
void
blas_trsm<octruple>(const CBLAS_SIDE Side,
		    const CBLAS_UPLO Uplo, const CBLAS_TRANSPOSE TransA,
		    const CBLAS_DIAG Diag, const int M, const int N,
		    const octruple &alpha,
		    const octruple *A, const int lda,
		    octruple *B, const int ldb);
#endif
template
void
blas_trsm<complex<quadruple> >(const CBLAS_SIDE Side,
			       const CBLAS_UPLO Uplo,
			       const CBLAS_TRANSPOSE TransA,
			       const CBLAS_DIAG Diag, const int M, const int N,
			       const complex<quadruple> &alpha,
			       const complex<quadruple> *A, const int lda,
			       complex<quadruple> *B, const int ldb);
#ifndef NO_OCTRUPLE
template
void
blas_trsm<complex<octruple> >(const CBLAS_SIDE Side,
			      const CBLAS_UPLO Uplo,
			      const CBLAS_TRANSPOSE TransA,
			      const CBLAS_DIAG Diag, const int M, const int N,
			      const complex<octruple> &alpha,
			      const complex<octruple> *A, const int lda,
			      complex<octruple> *B, const int ldb);
#endif

// dz gemm
#ifndef BLAS_GENERIC
template<>
void
blas_gemm<double>(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
		  int m, int n, int k,
		  const double &alpha, const double* A, int lda,
		  const double* B, int ldb, const double &beta,
		  double* C, int ldc )
{
#ifdef BLAS_FORTRAN
  unsigned char tra_ = 0;
  unsigned char trb_ = 0;
  switch(trA) {
  case CblasNoTrans:
    tra_ = 'n';
    break;
  case CblasTrans:
    tra_ = 't';
    break;
  case CblasConjTrans:
    tra_ = 'c';
    break;
  }
 switch(trB) {
 case CblasNoTrans:
   trb_ = 'n';
   break;
 case CblasTrans:
   trb_ = 't';
   break;
 case CblasConjTrans:
   trb_ = 'c';
   break;
 }
 FORTRAN_DECL(dgemm)(tra_, trb_, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
#else  
  cblas_dgemm(CblasColMajor,
	      trA, trB, (BLAS_INT)m, (BLAS_INT)n, k, alpha,
	      A, (BLAS_INT)lda, B, (BLAS_INT)ldb, beta, C, (BLAS_INT)ldc );
#endif
}

template<>
void
blas_gemm<complex<double> >(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
			    int m, int n, int k, const complex<double> &alpha, 
			    const complex<double>* A, int lda,
			    const complex<double>* B, int ldb,
			    const complex<double> &beta, 
			    complex<double>* C, int ldc )
{
#ifdef BLAS_FORTRAN
  unsigned char tra_ = 0;
  unsigned char trb_ = 0;
  switch(trA) {
  case CblasNoTrans:
    tra_ = 'n';
    break;
  case CblasTrans:
    tra_ = 't';
    break;
  case CblasConjTrans:
    tra_ = 'c';
    break;
  }
  switch(trB) {
  case CblasNoTrans:
    trb_ = 'n';
    break;
  case CblasTrans:
    trb_ = 't';
    break;
  case CblasConjTrans:
    trb_ = 'c';
    break;
  }
  FORTRAN_DECL(zgemm)(tra_, trb_, m, n, k, (const double *)&alpha,
		      (const double *)A, lda, (const double *)B, ldb, 
		      (const double *)&beta, (double *)C, ldc );
#else
  cblas_zgemm(CblasColMajor,
	      trA, trB, (BLAS_INT)m, (BLAS_INT)n, (BLAS_INT)k,
	      (const BLAS_VOID *)&alpha, (const BLAS_VOID *)A, (BLAS_INT)lda,
	      (const BLAS_VOID *)B, (BLAS_INT)ldb, 
	      (const BLAS_VOID *)&beta, C, (BLAS_INT)ldc );
#endif
}

template<>
void
blas_gemm<float>(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
		  int m, int n, int k,
		  const float &alpha, const float* A, int lda,
		  const float* B, int ldb, const float &beta,
		  float* C, int ldc )
{
#ifdef BLAS_FORTRAN
  unsigned char tra_ = 0;
  unsigned char trb_ = 0;
  switch(trA) {
  case CblasNoTrans:
    tra_ = 'n';
    break;
  case CblasTrans:
    tra_ = 't';
    break;
  case CblasConjTrans:
    tra_ = 'c';
    break;
  }
 switch(trB) {
 case CblasNoTrans:
   trb_ = 'n';
   break;
 case CblasTrans:
   trb_ = 't';
   break;
 case CblasConjTrans:
   trb_ = 'c';
   break;
 }
 FORTRAN_DECL(sgemm)(tra_, trb_, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
#else  
  cblas_sgemm(CblasColMajor,
	      trA, trB, (BLAS_INT)m, (BLAS_INT)n, k, alpha,
	      A, (BLAS_INT)lda, B, (BLAS_INT)ldb, beta, C, (BLAS_INT)ldc );
#endif
}

template<>
void
blas_gemm<complex<float> >(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
			    int m, int n, int k, const complex<float> &alpha, 
			    const complex<float>* A, int lda,
			    const complex<float>* B, int ldb,
			    const complex<float> &beta, 
			    complex<float>* C, int ldc )
{
#ifdef BLAS_FORTRAN
  unsigned char tra_ = 0;
  unsigned char trb_ = 0;
  switch(trA) {
  case CblasNoTrans:
    tra_ = 'n';
    break;
  case CblasTrans:
    tra_ = 't';
    break;
  case CblasConjTrans:
    tra_ = 'c';
    break;
  }
  switch(trB) {
  case CblasNoTrans:
    trb_ = 'n';
    break;
  case CblasTrans:
    trb_ = 't';
    break;
  case CblasConjTrans:
    trb_ = 'c';
    break;
  }
  FORTRAN_DECL(cgemm)(tra_, trb_, m, n, k, (const float *)&alpha,
		      (const float *)A, lda, (const float *)B, ldb, 
		      (const float *)&beta, (float *)C, ldc );
#else
  cblas_cgemm(CblasColMajor,
	      trA, trB, (BLAS_INT)m, (BLAS_INT)n, (BLAS_INT)k,
	      (const BLAS_VOID *)&alpha, (const BLAS_VOID *)A, (BLAS_INT)lda,
	      (const BLAS_VOID *)B, (BLAS_INT)ldb, 
	      (const BLAS_VOID *)&beta, C, (BLAS_INT)ldc );
#endif
}
#endif

template<typename T>
void
blas_gemm(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
	  int m, int n, int k,
	  const T &alpha,
	  const T* A, int lda,
	  const T* B, int ldb,
	  const T &beta,
	  T* C, int ldc )
{
  const T zero(0.0);
  const T one(1.0);
  // based on zgemm
  bool nota, notb;
  bool conja, conjb;

  nota = (trA == CblasNoTrans);
  notb = (trB == CblasNoTrans);
  conja = (trA == CblasConjTrans);
  conjb = (trB == CblasConjTrans);

  // without checking parameters
  
  if ((m == 0) || (n == 0) || (((alpha == zero || (k == 0))
				&& (beta == one)))) {
    return;
  }
     
  if (alpha == zero) {
    if (beta == zero) {
      for (int j = 0; j < n; j++) {
	for (int i = 0; i < m; i++) {
	  C[i + j * ldc] = zero;
	} // loop 10 : i
      }   // loop 20 : j
    }
    else {
      for (int j = 0; j < n; j++) {
	for (int i = 0; i < m; i++) {
	  C[i + j * ldc] *= beta;
	} // loop 30 : i
      }   // loop 40 : j
    }
    return;
  }
  if (notb) {
    if (nota) {
      //           Form  C := alpha*A*B + beta*C.
      for (int j = 0; j < n; j++) {
	if (beta == zero) {
	  for (int i = 0; i < m; i++) {
	    C[i + j * ldc] = zero;
	  } // loop 50 : i
	}
	else if (beta != one) {
	  for (int i = 0; i < m; i++) {
	    C[i + j * ldc] *= beta;
	  } // loop 60 : i
	}
	for (int l = 0; l < k; l++) {
	  if (B[l + j * ldb] != zero) {
	    T temp = alpha * B[l + j * ldb];
	    for (int i = 0; i < m; i++) {
	      C[i + j * ldc] += temp * A[i + l * lda];
	    } // loop 70 : i
	  }
	}     // loop 80 : l
      } // loop 90 : j
    }
    else if (conja) {
      // Form  C := alpha*conjg( A' )*B + beta*C.
      for (int j = 0; j < n; j++) {
	for (int i = 0; i < m; i++) {
	  T temp = zero;
	  for (int l = 0; l < k; l++) {
	    temp += blas_conj(A[l + i * lda]) * B[l + j * ldb];
	  } // loop 100 : l
	  if (beta == zero) {
	    C[i + j * ldc] = alpha * temp;
	  }
	  else {
	    C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
	  }
	} // loop 110 : i 
      } // loop 120 : j
    }
    else { // if (nota)
      //  Form  C := alpha*A**T*B + beta*C
      for (int j = 0; j < n; j++) {
	for (int i = 0; i < m; i++) {
	  T temp = zero;
	  for (int l = 0; l < k; l++) {
	    temp += A[l + i * lda] * B[l + j * ldb];
	  } // loop 130 : l
	  if (beta == zero) {
	    C[i + j * ldc] = alpha * temp;
	  }
	  else {
	    C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
	  }
	} // loop 140 : i 
      } // loop 150 : j 
    } // if (nota)
  }  // if (notb)
  else {
    if (nota) {
      if (conjb) {
	//     Form  C := alpha*A*conjg( B' ) + beta*C.
	for (int j = 0; j < n; j++) {
	  if (beta == zero) {
	    for (int i = 0; i < m; i++) {
	      C[i + j * ldc] = zero;
	    } // loop 160 : i
	  }
	  else if (beta != one) {
	    for (int i = 0; i < m; i++) {
	      C[i + j * ldc] *= beta;
	    } // loop 170 : i
	  }
	  for (int l = 0; l < k; l++) {
	    if (B[j + l * ldb] != zero) {
	      T temp = alpha * blas_conj(B[j + l * ldb]);
	      for (int i = 0; i < m; i++) {
		C[i + j * ldc] += temp * A[i + l * lda];
	      } // loop 180 : i
	    }
	  } // loop 190 : i
	} // loop 200 : j
      }
      else {
      //    Form  C := alpha*A*B**T + beta*C
	for (int j = 0; j < n; j++) {
	  if (beta == zero) {
	    for (int i = 0; i < m; i++) {
	      C[i + j * ldc] = zero;
	    } //loop 210 : i
	  }
	  else if (beta != one) {
	    for (int i = 0; i < m; i++) {
	      C[i + j * ldc] *= beta;
	    } // loop 220 : i
	  }
	  for (int l = 0; l < k; l++) {
	    if (B[j + l * ldb] != zero) {
	      T temp = alpha * B[j + l * ldb];
	      for (int i = 0; i < m; i++) {
		C[i + j * ldc] += temp * A[i + l * lda];
	      } // loop 230 : i
	    }
	  } // loop 240 : i
	} // loop 250 : j
      }
    } // if (nota)
    else if (conja) {
      if (conjb) {
	//  Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
	for (int j = 0; j < n; j++) {
	  for (int i = 0; i < m; i++) {
	    T temp = zero;
	    for (int l = 0; l < k; l++) {
	      temp += blas_conj(A[l + i * lda]) * blas_conj(B[j + l * ldb]);
	    } // loop 260 : l
	    if (beta == zero) {
	      C[i + j * ldc] = alpha * temp;
	    }
	    else {
	      C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
	    }
	  } // loop 270 : i
	}   // loop 280 : i
      }
      else {
	//  Form  C := alpha*conjg( A' )*B' + beta*C
	for (int j = 0; j < n; j++) {
	  for (int i = 0; i < m; i++) {
	    T temp = zero;
	    for (int l = 0; l < k; l++) {
	      temp += blas_conj(A[l + i * lda]) * B[j + l * ldb];
	    } // loop 290 : l
	    if (beta == zero) {
	      C[i + j * ldc] = alpha * temp;
	    }
	    else {
	      C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
	    }
	  } // loop 300 : i
	}   // loop 310 : i
      }
    }
    else {
      if (conjb) {
	//   Form  C := alpha*A'*conjg( B' ) + beta*C
	for (int j = 0; j < n; j++) {
	  for (int i = 0; i < m; i++) {
	    T temp = zero;
	    for (int l = 0; l < k; l++) {
	      temp += A[l + i * lda] * blas_conj(B[j + l * ldb]);
	    } // loop 320 : l
	    if (beta == zero) {
	      C[i + j * ldc] = alpha * temp;
	    }
	    else {
	      C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
	    }
	  } // loop 330 : i
	}   // loop 340 : i
      }
      else {
	//  Form  C := alpha*A**T*B**T + beta*C
	for (int j = 0; j < n; j++) {
	  for (int i = 0; i < m; i++) {
	    T temp = zero;
	    for (int l = 0; l < k; l++) {
	      temp += A[l + i * lda] * B[j + l * ldb];
	    } // loop 350 : l
	    if (beta == zero) {
	      C[i + j * ldc] = alpha * temp;
	    }
	    else {
	      C[i + j * ldc] = alpha * temp + beta * C[i + j * ldc];
	    }
	  } // loop 360 : i
	}   // loop 370 : i
      }
    } // if (nota)
  }  // if (notb)
}
// explicit instantiation of blas_gemm
#ifdef FORCE_EXPLICIT_INSTANTIATION
template
void
blas_gemm<double>(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
		  int m, int n, int k,
		  const double &alpha,
		  const double* A, int lda,
		  const double* B, int ldb,
		  const double &beta,
		  double* C, int ldc );

template
void
blas_gemm<complex<double> >(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
			    int m, int n, int k,
			    const complex<double> &alpha_, 
			    const complex<double>* A, int lda,
			    const complex<double>* B, int ldb,
			    const complex<double> &beta_, 
			    complex<double>* C, int ldc );

template
void
blas_gemm<float>(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
		  int m, int n, int k,
		  const float &alpha,
		  const float* A, int lda,
		  const float* B, int ldb,
		  const float &beta,
		  float* C, int ldc );

template
void
blas_gemm<complex<float> >(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
			    int m, int n, int k,
			    const complex<float> &alpha_, 
			    const complex<float>* A, int lda,
			    const complex<float>* B, int ldb,
			    const complex<float> &beta_, 
			    complex<float>* C, int ldc );
#endif
template
void
blas_gemm<quadruple>(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
		     int m, int n, int k,
		     const quadruple &alpha,
		     const quadruple* A, int lda,
		     const quadruple* B, int ldb,
		     const quadruple &beta,
		     quadruple* C, int ldc );
template
void
blas_gemm<complex<quadruple> >(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
			       int m, int n, int k,
			       const complex<quadruple> &alpha_, 
			       const complex<quadruple>* A, int lda,
			       const complex<quadruple>* B, int ldb,
			       const complex<quadruple> &beta_, 
			       complex<quadruple>* C, int ldc );
#ifndef NO_OCTRUPLE
template
void
blas_gemm<octruple>(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
		    int m, int n, int k,
		    const octruple &alpha,
		    const octruple* A, int lda,
		    const octruple* B, int ldb,
		    const octruple &beta,
		    octruple* C, int ldc );
template
void
blas_gemm<complex<octruple> >(CBLAS_TRANSPOSE trA, CBLAS_TRANSPOSE trB, 
			      int m, int n, int k,
			      const complex<octruple> &alpha_, 
			      const complex<octruple>* A, int lda,
			      const complex<octruple>* B, int ldb,
			      const complex<octruple> &beta_, 
			      complex<octruple>* C, int ldc );
#endif

// OTHERS
template<typename T, typename U>
U blas_l2norm(const int n, T *x, const int incX)
{
  U tmp;
  tmp = blas_dot<T>(n, x, incX, x, incX);
  return sqrt<U>(tmp); // works for T = double, quadruple, etc. not for complex<U>
}

template<>
double blas_l2norm<complex<double>, double>(const int n, complex<double> *x,
					    const int incX)
{
  complex<double> tmp;
  tmp = blas_dot(n, x, incX, x, incX);
  return sqrt<double>(tmp.real()); 
}

template<>
float blas_l2norm<complex<float>, float>(const int n, complex<float> *x,
					    const int incX)
{
  complex<float> tmp;
  tmp = blas_dot(n, x, incX, x, incX);
  return sqrt<float>(tmp.real()); 
}

template<>
quadruple blas_l2norm<complex<quadruple>, quadruple>(const int n,
						     complex<quadruple> *x,
						     const int incX)
{
  complex<quadruple> tmp;
  tmp = blas_dot<complex<quadruple> >(n, x, incX, x, incX);
  return sqrt<quadruple>(tmp.real()); //
}
#ifndef NO_OCTRUPLE
template<>
octruple blas_l2norm<complex<octruple>, octruple>(const int n,
						  complex<octruple> *x,
						  const int incX)
{
  complex<octruple> tmp;
   tmp = blas_dot<complex<octruple> >(n, x, incX, x, incX);
   return sqrt<octruple>(tmp.real()); //
}
#endif

template
double blas_l2norm<double, double>(const int n, double *x, const int incX);

template
quadruple blas_l2norm<quadruple, quadruple>(const int n, quadruple *x,
					    const int incX);
template
float blas_l2norm<float, float>(const int n, float *x, const int incX);

#ifdef FORCE_EXPLICIT_INSTANTIATION
template
double blas_l2norm<complex<double>, double>(const int n, complex<double> *x,
					    const int incX);
template
quadruple blas_l2norm<complex<quadruple>, quadruple>(const int n,
						     complex<quadruple> *x,
						     const int incX);
#ifndef NO_OCTRUPLE
template 
octruple blas_l2norm<complex<octruple>, octruple>(const int n,
						  complex<octruple> *x,
						  const int incX);
#endif
#endif
#ifndef NO_OCTRUPLE
template 
octruple blas_l2norm<octruple, octruple>(const int n, octruple *x,
					 const int incX);
#endif

double blas_l2norm_lower_prec(const int n, quadruple *x, const int incX)
{
  quadruple tmp;
  tmp = blas_dot(n, x, incX, x, incX);
  return sqrt<double>(quad2double(tmp)); //
}

quadruple blas_l2norm_lower_prec(const int n, octruple *x, const int incX)
{
  octruple tmp;
  tmp = blas_dot(n, x, incX, x, incX);
  return sqrt<quadruple>(oct2quad(tmp)); //
}

double blas_l2norm_lower_prec(const int n, complex<quadruple> *x,
			      const int incX)
{
  complex<quadruple> tmp;
  tmp = blas_dot(n, x, incX, x, incX);
  return sqrt<double>(quad2double(tmp.real())); //
}

quadruple blas_l2norm_lower_prec(const int n, complex<octruple> *x,
				 const int incX)
{
  complex<octruple> tmp;
  tmp = blas_dot(n, x, incX, x, incX);
  return sqrt<quadruple>(oct2quad(tmp.real())); //
}

template <typename T, typename U>
U blas_l2norm2(const int n, T *x, const int incX)
{
  U tmp;
  tmp = blas_dot<T>(n, x, incX, x, incX);
  return tmp;
}

template<>
double blas_l2norm2<complex<double>, double>(const int n, complex<double> *x,
					     const int incX)
{
  complex<double> tmp;
  tmp = blas_dot<complex<double> >(n, x, incX, x, incX);
  return tmp.real(); 
}

template<>
float blas_l2norm2<complex<float>, float>(const int n, complex<float> *x,
					     const int incX)
{
  complex<float> tmp;
  tmp = blas_dot<complex<float> >(n, x, incX, x, incX);
  return tmp.real(); 
}

template<>
quadruple blas_l2norm2<complex<quadruple>, quadruple>(const int n,
						      complex<quadruple> *x,
						      const int incX)
{
  complex<quadruple> tmp;
  tmp = blas_dot<complex<quadruple> >(n, x, incX, x, incX);
  return tmp.real(); //
}

template<>
octruple blas_l2norm2<complex<octruple>, octruple>(const int n,
						   complex<octruple> *x,
						   const int incX)
{
  complex<octruple> tmp;
  tmp = blas_dot<complex<octruple> >(n, x, incX, x, incX);
  return tmp.real(); //
}

template
double blas_l2norm2<double, double>(const int n, double *x, const int incX);

template
float blas_l2norm2<float, float>(const int n, float *x, const int incX);

template
quadruple blas_l2norm2<quadruple, quadruple>(const int n, quadruple *x,
					     const int incX);

#ifdef FORCE_EXPLICIT_INSTANTIATION
template
double blas_l2norm2<complex<double>, double>(const int n,
					     complex<double> *x,
					     const int incX);

template
quadruple blas_l2norm2<complex<quadruple>, quadruple>(const int n,
						      complex<quadruple> *x,
						      const int incX);

template
octruple blas_l2norm2<octruple, octruple>(const int n,
					octruple *x,
					const int incX);

template
octruple blas_l2norm2<complex<octruple>, octruple>(const int n,
						   complex<octruple> *x,
						   const int incX);

#endif
