/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE blas1c.h.
   Interface to blas 1 and blas 2 FORTRAN routines.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arch.h"
#include "blas1f.h"

#ifndef BLAS1C_H
#define BLAS1C_H

// ASSUM

inline float assum(const integer &n, const float dx[], const integer &incx) {
  return F77NAME(sasum)(&n, dx, &incx);
} // assum (float)

inline double assum(const integer &n, const double dx[], const integer &incx) {
  return F77NAME(dasum)(&n, dx, &incx);
} // assum (double)

#ifdef ARCOMP_H
inline float assum(const integer &n, const arcomplex<float> dx[],
                   const integer &incx) {
  return F77NAME(scasum)(&n, dx, &incx);
} // assum (arcomplex<float>)

inline double assum(const integer &n, const arcomplex<double> dx[],
                    const integer &incx) {
  return F77NAME(dzasum)(&n, dx, &incx);
} // assum (arcomplex<double>)
#endif

// AXPY

inline void axpy(const integer &n, const float &da, const float dx[],
                 const integer &incx, float dy[], const integer &incy) {
  F77NAME(saxpy)(&n, &da, dx, &incx, dy, &incy);
} // axpy (float)

inline void axpy(const integer &n, const double &da, const double dx[],
                 const integer &incx, double dy[], const integer &incy) {
  F77NAME(daxpy)(&n, &da, dx, &incx, dy, &incy);
} // axpy (double)

#ifdef ARCOMP_H
inline void axpy(const integer &n, const arcomplex<float> &da,
                 const arcomplex<float> dx[], const integer &incx,
                 arcomplex<float> dy[], const integer &incy) {
  F77NAME(caxpy)(&n, &da, dx, &incx, dy, &incy);
} // axpy (arcomplex<float>)

inline void axpy(const integer &n, const arcomplex<double> &da,
                 const arcomplex<double> dx[], const integer &incx,
                 arcomplex<double> dy[], const integer &incy) {
  F77NAME(zaxpy)(&n, &da, dx, &incx, dy, &incy);
} // axpy (arcomplex<double>)
#endif

// COPY

inline void copy(const integer &n, const float dx[], const integer &incx,
                 float dy[], const integer &incy) {
  if (dx != dy) F77NAME(scopy)(&n, dx, &incx, dy, &incy);
} // copy (float)

inline void copy(const integer &n, const double dx[], const integer &incx,
                 double dy[], const integer &incy) {
  if (dx != dy) F77NAME(dcopy)(&n, dx, &incx, dy, &incy);
} // copy (double)

#ifdef ARCOMP_H
inline void copy(const integer &n, const arcomplex<float> dx[], 
                 const integer &incx, arcomplex<float> dy[], 
                 const integer &incy) {
  if (dx != dy) F77NAME(ccopy)(&n, dx, &incx, dy, &incy);
} // copy (arcomplex<float>)

inline void copy(const integer &n, const arcomplex<double> dx[], 
                 const integer &incx, arcomplex<double> dy[], 
                 const integer &incy) {
  if (dx != dy) F77NAME(zcopy)(&n, dx, &incx, dy, &incy);
} // copy (arcomplex<double>)
#endif

// DOT

inline float dot(const integer &n, const float dx[], const integer &incx,
                 const float dy[], const integer &incy) {
  return F77NAME(sdot)(&n, dx, &incx, dy, &incy);
} // dot (float)

inline double dot(const integer &n, const double dx[], const integer &incx,
                  const double dy[], const integer &incy) {
  return F77NAME(ddot)(&n, dx, &incx, dy, &incy);
} // dot (double)

#ifdef ARCOMP_H
inline arcomplex<float> dotc(const integer &n, const arcomplex<float> dx[], 
                           const integer &incx,const arcomplex<float> dy[], 
                           const integer &incy) {
  arcomplex<float> tmp;
  F77NAME(cdotc)(&tmp, &n, dx, &incx, dy, &incy);
  return tmp;
} // dotc (arcomplex<float>)

inline arcomplex<double> dotc(const integer &n, const arcomplex<double> dx[], 
                            const integer &incx, const arcomplex<double> dy[], 
                            const integer &incy) {
  arcomplex<double> tmp;
  F77NAME(zdotc)(&tmp, &n, dx, &incx, dy, &incy);
  return tmp;
} // dotc (arcomplex<double>)

inline arcomplex<float> dotu(const integer &n, const arcomplex<float> dx[], 
                           const integer &incx, const arcomplex<float> dy[], 
                           const integer &incy) {
  arcomplex<float> tmp;
  F77NAME(cdotu)(&tmp, &n, dx, &incx, dy, &incy);
  return tmp;
} // dotu (arcomplex<float>)

inline arcomplex<double> dotu(const integer &n, const arcomplex<double> dx[], 
                            const integer &incx, const arcomplex<double> dy[], 
                            const integer &incy) {
  arcomplex<double> tmp;
  F77NAME(zdotu)(&tmp, &n, dx, &incx, dy, &incy);
  return tmp;
} // dotu (arcomplex<double>)
#endif

// NRM2

inline float nrm2(const integer &n, const float dx[], const integer &incx) {
  return F77NAME(snrm2)(&n, dx, &incx);
} // nrm2 (float)

inline double nrm2(const integer &n, const double dx[], const integer &incx) {
  return F77NAME(dnrm2)(&n, dx, &incx);
} // nrm2 (double)

#ifdef ARCOMP_H
inline float nrm2(const integer &n, const arcomplex<float> dx[], 
                  const integer &incx) {
  return F77NAME(scnrm2)(&n, dx, &incx);
} // nrm2 (complex <float>)

inline double nrm2(const integer &n, const arcomplex<double> dx[],
                   const integer &incx) {
  return F77NAME(dznrm2)(&n, dx, &incx);
} // nrm2 (complex <double>)
#endif

// ROT

inline void rot(const integer &n, float dx[], const integer &incx, float dy[],
           const integer &incy, const float &c, const float &s) {
  F77NAME(srot)(&n, dx, &incx, dy, &incy, &c, &s);
} // rot (float)

inline void rot(const integer &n, double dx[], const integer &incx, 
                double dy[], const integer &incy, const double &c, 
                const double &s) {
  F77NAME(drot)(&n, dx, &incx, dy, &incy, &c, &s);
} // rot (double)

// ROTG

inline void rotg(float &da, float &db, float &c, float &s) {
  F77NAME(srotg)(&da, &db, &c, &s);
} // rotg (float)

inline void rotg(double &da, double &db, double &c, double &s) {
  F77NAME(drotg)(&da, &db, &c, &s);
} // rotg (double)

// SCAL

inline void scal(const integer &n, float &da, float dx[], const integer &incx) {
  F77NAME(sscal)(&n, &da, dx, &incx);
} // scal (float)

inline void scal(const integer &n, double &da, double dx[], const integer &incx) {
  F77NAME(dscal)(&n, &da, dx, &incx);
} // scal (double)

#ifdef ARCOMP_H
inline void scal(const integer &n, const arcomplex<float> &da,
                 arcomplex<float> dx[], const integer &incx) {
  F77NAME(cscal)(&n, &da, dx, &incx);
} // scal (arcomplex<float>)

inline void scal(const integer &n, const arcomplex<double> &da,
                 arcomplex<double> dx[], const integer &incx) {
  F77NAME(zscal)(&n, &da, dx, &incx);
} // scal (arcomplex<double>)

inline void sscal(const integer &n, const float &da, arcomplex<float> dx[],
                  const integer &incx) {
  F77NAME(csscal)(&n, &da, dx, &incx);
} // sscal (arcomplex<float>)

inline void sscal(const integer &n, const double &da, arcomplex<double> dx[],
                  const integer &incx) {
  F77NAME(zdscal)(&n, &da, dx, &incx);
} // sscal (arcomplex<double>)
#endif

// SWAP

inline void swap(const integer &n, float dx[], const integer &incx,
                 float dy[], const integer &incy) {
  F77NAME(sswap)(&n, dx, &incx, dy, &incy);
} // swap (float)

inline void swap(const integer &n, double dx[], const integer &incx,
                 double dy[], const integer &incy) {
  F77NAME(dswap)(&n, dx, &incx, dy, &incy);
} // swap (double)

#ifdef ARCOMP_H
inline void swap(const integer &n, arcomplex<float> dx[], const integer &incx,
                 arcomplex<float> dy[], const integer &incy) {
  F77NAME(cswap)(&n, dx, &incx, dy, &incy);
} // swap (arcomplex<float>)

inline void swap(const integer &n, arcomplex<double> dx[], const integer &incx,
                 arcomplex<double> dy[], const integer &incy) {
  F77NAME(zswap)(&n, dx, &incx, dy, &incy);
} // swap (arcomplex<double>)
#endif

// AMAX

inline integer amax(const integer &n, const float dx[], const integer &incx) {
  return F77NAME(isamax)(&n, dx, &incx);
} // amax (float)

inline integer amax(const integer &n, const double dx[], const integer &incx) {
  return F77NAME(idamax)(&n, dx, &incx);
} // amax (double)

#ifdef ARCOMP_H
inline integer amax(const integer &n, const arcomplex<float> dx[], 
                    const integer &incx) {
  return F77NAME(icamax)(&n, dx, &incx);
} // amax (arcomplex<float>)

inline integer amax(const integer &n, const arcomplex<double> dx[], 
                    const integer &incx) {
  return F77NAME(izamax)(&n, dx, &incx);
} // amax (arcomplex<double>)
#endif

// GEMV

inline void gemv(const char* trans, const integer &m, const integer &n, 
                 const float &alpha, const float a[], const integer &lda, 
                 const float x[], const integer &incx, const float &beta, 
                 float y[], const integer &incy) {
  F77NAME(sgemv)(trans, &m, &n, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gemv (float)

inline void gemv(const char* trans, const integer &m, const integer &n, 
                 const double &alpha, const double a[], const integer &lda, 
                 const double x[], const integer &incx, const double &beta, 
                 double y[], const integer &incy) {
  F77NAME(dgemv)(trans, &m, &n, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gemv (double)

#ifdef ARCOMP_H
inline void gemv(const char* trans, const integer &m, 
                 const integer &n, const arcomplex<float> &alpha,
                 const arcomplex<float> a[], const integer &lda, 
                 const arcomplex<float> x[], const integer &incx, 
                 const arcomplex<float> &beta, arcomplex<float> y[],
                 const integer &incy) {
  F77NAME(cgemv)(trans, &m, &n, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gemv (arcomplex<float>)

inline void gemv(const char* trans, const integer &m, 
                 const integer &n, const arcomplex<double> &alpha,
                 const arcomplex<double> a[], const integer &lda, 
                 const arcomplex<double> x[], const integer &incx, 
                 const arcomplex<double> &beta, arcomplex<double> y[],
                 const integer &incy) {
  F77NAME(zgemv)(trans, &m, &n, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gemv (arcomplex<double>)
#endif

// GBMV

inline void gbmv(const char* trans, const integer &m, const integer &n, 
                 const integer &kl, const integer &ku, const float &alpha,
                 const float a[], const integer &lda, const float x[],
                 const integer &incx, const float &beta, float y[],
                 const integer &incy) {
  F77NAME(sgbmv)(trans, &m, &n, &kl, &ku, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gbmv (float)

inline void gbmv(const char* trans, const integer &m, const integer &n, 
                 const integer &kl, const integer &ku, const double &alpha,
                 const double a[], const integer &lda, const double x[],
                 const integer &incx, const double &beta, double y[],
                 const integer &incy) {
  F77NAME(dgbmv)(trans, &m, &n, &kl, &ku, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gbmv (double)

#ifdef ARCOMP_H
inline void gbmv(const char* trans, const integer &m, 
                 const integer &n, const integer &kl, 
                 const integer &ku, const arcomplex<float> &alpha,
                 const arcomplex<float> a[], const integer &lda, 
                 const arcomplex<float> x[], const integer &incx, 
                 const arcomplex<float> &beta, arcomplex<float> y[],
                 const integer &incy) {
  F77NAME(cgbmv)(trans, &m, &n, &kl, &ku, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gbmv (arcomplex<float>)

inline void gbmv(const char* trans, const integer &m, 
                 const integer &n, const integer &kl, 
                 const integer &ku, const arcomplex<double> &alpha,
                 const arcomplex<double> a[], const integer &lda, 
                 const arcomplex<double> x[], const integer &incx, 
                 const arcomplex<double> &beta, arcomplex<double> y[],
                 const integer &incy) {
  F77NAME(zgbmv)(trans, &m, &n, &kl, &ku, &alpha, a, &lda, 
                 x, &incx, &beta, y, &incy);
} // gbmv (arcomplex<double>)
#endif

// SBMV

inline void sbmv(const char* uplo, const integer &n, const integer &k, 
                 const float &alpha, const float a[], const integer &lda, 
                 const float x[], const integer &incx, const float &beta, 
                 float y[], const integer &incy) {
  F77NAME(ssbmv)(uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
} // sbmv (float)

inline void sbmv(const char* uplo, const integer &n, const integer &k, 
                 const double &alpha, const double a[], const integer &lda, 
                 const double x[], const integer &incx, const double &beta, 
                 double y[], const integer &incy) {
  F77NAME(dsbmv)(uplo, &n, &k, &alpha, a, &lda, x, &incx, &beta, y, &incy);
} // sbmv (double)


#endif // BLAS1C_H
