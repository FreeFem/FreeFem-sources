/*
  ARPACK++ v1.0 8/1/1997
  c++ interface to ARPACK code.

  MODULE blas1f.h
  BLAS 1 and BLAS 2 FORTRAN routines.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas
*/

#ifndef BLAS1F_H
#define BLAS1F_H

#include "arch.h"

extern "C"
{

  // Single precision real routines.

  float F77NAME(sasum)(const integer *n, const float *dx, const integer *incx);

  void F77NAME(saxpy)(const integer *n, const float *da, const float *dx,
                      const integer *incx, float *dy, const integer *incy);

  void F77NAME(scopy)(const integer *n, const float *dx, const integer *incx,
                      float *dy, const integer *incy);

  float F77NAME(sdot)(const integer *n, const float *dx, const integer *incx,
                      const float *dy, const integer *incy);

  float F77NAME(snrm2)(const integer *n, const float *dx, const integer *incx);

  void F77NAME(srot)(const integer *n, float *dx, const integer *incx, float *dy,
                     const integer *incy, const float *c, const float *s);

  void F77NAME(srotg)(float *da, float *db, float *c, float *s);

  void F77NAME(sscal)(const integer *n, float *da, float *dx, const integer *incx);

  void F77NAME(sswap)(const integer *n, float *dx, const integer *incx,
                      float *dy, const integer *incy);

  integer F77NAME(isamax)(const integer *n, const float *dx, const integer *incx);

  void F77NAME(sgemv)(const char* trans, const integer *m, const integer *n, 
                      const float *alpha, const float *a, const integer *lda, 
                      const float *x, const integer *incx, const float *beta, 
                      float *y, const integer *incy);

  void F77NAME(sgbmv)(const char* trans, const integer *m, const integer *n, 
                      const integer *kl, const integer *ku, const float *alpha,
                      const float *a, const integer *lda, const float *x,
                      const integer *incx, const float *beta, float *y,
                      const integer *incy);

  void F77NAME(ssbmv)(const char* uplo, const integer *n, const integer *k, 
                      const float *alpha, const float *a, const integer *lda, 
                      const float *x, const integer *incx, const float *beta, 
                      float *y, const integer *incy);

// Double precision real routines.

  double F77NAME(dasum)(const integer *n, const double *dx, const integer *incx);

  void F77NAME(daxpy)(const integer *n, const double *da, const double *dx,
                      const integer *incx, double *dy, const integer *incy);

  void F77NAME(dcopy)(const integer *n, const double *dx, const integer *incx,
                      double *dy, const integer *incy);

  double F77NAME(ddot)(const integer *n, const double *dx, const integer *incx,
                       const double *dy, const integer *incy);

  double F77NAME(dnrm2)(const integer *n, const double *dx, const integer *incx);

  void F77NAME(drot)(const integer *n, double *dx, const integer *incx, double *dy,
                     const integer *incy, const double *c, const double *s);

  void F77NAME(drotg)(double *da, double *db, double *c, double *s);

  void F77NAME(dscal)(const integer *n, double *da, double *dx, const integer *incx);

  void F77NAME(dswap)(const integer *n, double *dx, const integer *incx,
                      double *dy, const integer *incy);

  integer F77NAME(idamax)(const integer *n, const double *dx, const integer *incx);

  void F77NAME(dgemv)(const char* trans, const integer *m, const integer *n, 
                      const double *alpha, const double *a, const integer *lda,
                      const double *x, const integer *incx, const double *beta,
                      double *y, const integer *incy);

  void F77NAME(dgbmv)(const char* trans, const integer *m, const integer *n, 
                      const integer *kl, const integer *ku, const double *alpha,
                      const double *a, const integer *lda, const double *x,
                      const integer *incx, const double *beta, double *y,
                      const integer *incy);

  void F77NAME(dsbmv)(const char* uplo, const integer *n, const integer *k, 
                      const double *alpha, const double *a, const integer *lda, 
                      const double *x, const integer *incx, const double *beta, 
                      double *y, const integer *incy);

  // Single precision complex routines.

#ifdef ARCOMP_H

  void F77NAME(cdotc)(arcomplex<float> *c, const integer *n,
                      const arcomplex<float> *cx, const integer *incx,
                      const arcomplex<float> *cy, const integer *incy);

  void F77NAME(cdotu)(arcomplex<float> *c, const integer *n,
                      const arcomplex<float> *cx, const integer *incx,
                      const arcomplex<float> *cy, const integer *incy);

  void F77NAME(caxpy)(const integer *n, const arcomplex<float> *da,
                      const arcomplex<float> *dx, const integer *incx,
                      arcomplex<float> *dy, const integer *incy);

  void F77NAME(ccopy)(const integer *n, const arcomplex<float> *dx,
                      const integer *incx, arcomplex<float> *dy,
                      const integer *incy);

  float F77NAME(scasum)(const integer *n, const arcomplex<float> *dx,
                        const integer *incx);

  float F77NAME(scnrm2)(const integer *n, const arcomplex<float> *dx,
                        const integer *incx);

  void F77NAME(csscal)(const integer *n, const float *da, arcomplex<float> *dx,
                       const integer *incx);

  void F77NAME(cscal)(const integer *n, const arcomplex<float> *da,
                      arcomplex<float> *dx, const integer *incx);

  integer F77NAME(icamax)(const integer *n, const arcomplex<float> *dx,
                          const integer *incx);

  void F77NAME(cswap)(const integer *n, arcomplex<float> *dx, 
                      const integer *incx, arcomplex<float> *dy, 
                      const integer *incy);

  void F77NAME(cgemv)(const char* trans, const integer *m, 
                      const integer *n, const arcomplex<float> *alpha,
                      const arcomplex<float> *a, const integer *lda, 
                      const arcomplex<float> *x, const integer *incx, 
                      const arcomplex<float> *beta, arcomplex<float> *y,
                      const integer *incy);

  void F77NAME(cgbmv)(const char* trans, const integer *m, 
                      const integer *n, const integer *kl, 
                      const integer *ku, const arcomplex<float> *alpha,
                      const arcomplex<float> *a, const integer *lda, 
                      const arcomplex<float> *x, const integer *incx, 
                      const arcomplex<float> *beta, arcomplex<float> *y,
                      const integer *incy);

  // Double precision complex routines.

  void F77NAME(zdotc)(arcomplex<double> *c, const integer *n,
                      const arcomplex<double> *cx, const integer *incx,
                      const arcomplex<double> *cy, const integer *incy);

  void F77NAME(zdotu)(arcomplex<double> *c, const integer *n,
                      const arcomplex<double> *cx, const integer *incx,
                      const arcomplex<double> *cy, const integer *incy);

  void F77NAME(zaxpy)(const integer *n, const arcomplex<double> *da,
                      const arcomplex<double> *dx, const integer *incx,
                      arcomplex<double> *dy, const integer *incy);

  void F77NAME(zcopy)(const integer *n, const arcomplex<double> *dx,
                      const integer *incx, arcomplex<double> *dy,
                      const integer *incy);

  double  F77NAME(dzasum)(const integer *n, const arcomplex<double> *dx,
                          const integer *incx);

  double  F77NAME(dznrm2)(const integer *n, const arcomplex<double> *dx,
                          const integer *incx);

  void F77NAME(zdscal)(const integer *n, const double *da, arcomplex<double> *dx,
                       const integer *incx);

  void F77NAME(zscal)(const integer *n, const arcomplex<double> *da,
                      arcomplex<double> *dx, const integer *incx);

  integer F77NAME(izamax)(const integer *n, const arcomplex<double> *dx,
                          const integer *incx);

  void F77NAME(zswap)(const integer *n, arcomplex<double> *dx,
                      const integer *incx, arcomplex<double> *dy, 
                      const integer *incy);

  void F77NAME(zgemv)(const char* trans, const integer *m, 
                      const integer *n, const arcomplex<double> *alpha,
                      const arcomplex<double> *a, const integer *lda, 
                      const arcomplex<double> *x, const integer *incx, 
                      const arcomplex<double> *beta, arcomplex<double> *y,
                      const integer *incy);

  void F77NAME(zgbmv)(const char* trans, const integer *m, 
                      const integer *n, const integer *kl, 
                      const integer *ku, const arcomplex<double> *alpha,
                      const arcomplex<double> *a, const integer *lda, 
                      const arcomplex<double> *x, const integer *incx, 
                      const arcomplex<double> *beta, arcomplex<double> *y,
                      const integer *incy);

#endif // ARCOMP_H

}
#endif // BLAS1F_H

