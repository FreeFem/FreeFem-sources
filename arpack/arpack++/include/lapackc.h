/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE lapackc.h.
   Interface to LAPACK FORTRAN routines.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#include "arch.h"
#include "lapackf.h"

#ifndef LAPACKC_H
#define LAPACKC_H


// LAPY2

inline float lapy2(const float &x, const float &y) {
  return F77NAME(slapy2)(&x, &y);
} // lapy2 (float)

inline double lapy2(const double &x, const double &y) {
  return F77NAME(dlapy2)(&x, &y);
} // lapy2 (double)


// LACPY

inline void lacpy(const char* uplo, const integer &m, const integer &n,
                  const float a[], const integer &lda, float b[], 
                  const integer &ldb) {
  F77NAME(slacpy)(uplo, &m, &n, a, &lda, b, &ldb);
} // lacpy (float)

inline void lacpy(const char* uplo, const integer &m, const integer &n,
                  const double a[], const integer &lda, double b[], 
                  const integer &ldb) {
  F77NAME(dlacpy)(uplo, &m, &n, a, &lda, b, &ldb);
} // lacpy (double)

#ifdef ARCOMP_H
inline void lacpy(const char* uplo, const integer &m, const integer &n,
                  const arcomplex<float> a[], const integer &lda, 
                  arcomplex<float> b[], const integer &ldb) {
  F77NAME(clacpy)(uplo, &m, &n, a, &lda, b, &ldb);
} // lacpy (arcomplex<float>)

inline void lacpy(const char* uplo, const integer &m, const integer &n,
                  const arcomplex<double> a[], const integer &lda, 
                  arcomplex<double> b[], const integer &ldb) {
  F77NAME(zlacpy)(uplo, &m, &n, a, &lda, b, &ldb);
} // lacpy (arcomplex<double>)
#endif


// GTTRF

inline void gttrf(const integer &n, float dl[], float d[], float du[],
                  float du2[], integer ipiv[], integer &info) {
  F77NAME(sgttrf)(&n, dl, d, du, du2, ipiv, &info);
} // gttrf (float)

inline void gttrf(const integer &n, double dl[], double d[], double du[],
                  double du2[], integer ipiv[], integer &info) {
  F77NAME(dgttrf)(&n, dl, d, du, du2, ipiv, &info);
} // gttrf (double)

#ifdef ARCOMP_H
inline void gttrf(const integer &n, arcomplex<float> dl[], arcomplex<float> d[],
                  arcomplex<float> du[], arcomplex<float> du2[], integer ipiv[],
                  integer &info) {
  F77NAME(cgttrf)(&n, dl, d, du, du2, ipiv, &info);
} // gttrf (arcomplex<float>)

inline void gttrf(const integer &n, arcomplex<double> dl[], arcomplex<double> d[],
                  arcomplex<double> du[], arcomplex<double> du2[], integer ipiv[],
                  integer &info) {
  F77NAME(zgttrf)(&n, dl, d, du, du2, ipiv, &info);
} // gttrf (arcomplex<double>)
#endif


// GTTRS

inline void gttrs(const char* trans, const integer &n, const integer &nrhs,
                  const float dl[], const float d[], const float du[],
                  const float du2[], const integer ipiv[], float b[],
                  const integer &ldb, integer &info) {
  F77NAME(sgttrs)(trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);
} // gttrs (float)

inline void gttrs(const char* trans, const integer &n, const integer &nrhs,
                  const double dl[], const double d[], const double du[],
                  const double du2[], const integer ipiv[], double b[],
                  const integer &ldb, integer &info) {
  F77NAME(dgttrs)(trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);
} // gttrs (double)

#ifdef ARCOMP_H

inline void gttrs(const char* trans, const integer &n, const integer &nrhs,
                  const arcomplex<float> dl[], const arcomplex<float> d[],
                  const arcomplex<float> du[], const arcomplex<float> du2[],
                  const integer ipiv[], arcomplex<float> b[],
                  const integer &ldb, integer &info) {
  F77NAME(cgttrs)(trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);
} // gttrs (arcomplex<float>)

inline void gttrs(const char* trans, const integer &n, const integer &nrhs,
                  const arcomplex<double> dl[], const arcomplex<double> d[],
                  const arcomplex<double> du[], const arcomplex<double> du2[],
                  const integer ipiv[], arcomplex<double> b[],
                  const integer &ldb, integer &info) {
  F77NAME(zgttrs)(trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info);
} // gttrs (arcomplex<double>)
#endif


// GBTRF

inline void gbtrf(const integer &m, const integer &n, const integer &kl,
                  const integer &ku, float ab[], const integer &ldab,
                  integer ipiv[], integer &info) {
  F77NAME(sgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
} // gbtrf (float)

inline void gbtrf(const integer &m, const integer &n, const integer &kl,
                  const integer &ku, double ab[], const integer &ldab,
                  integer ipiv[], integer &info) {
  F77NAME(dgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
} // gbtrf (double)

#ifdef ARCOMP_H
inline void gbtrf(const integer &m, const integer &n, const integer &kl,
                  const integer &ku, arcomplex<float> ab[], 
                  const integer &ldab, integer ipiv[], integer &info) {
  F77NAME(cgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
} // gbtrf (arcomplex<float>)

inline void gbtrf(const integer &m, const integer &n, const integer &kl,
                  const integer &ku, arcomplex<double> ab[], 
                  const integer &ldab, integer ipiv[], integer &info) {
  F77NAME(zgbtrf)(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
} // gbtrf (arcomplex<double>)
#endif


// GBTRS

inline void gbtrs(const char* trans, const integer &n, const integer &kl,
                  const integer &ku, const integer &nrhs, const float ab[], 
                  const integer &ldab, const integer ipiv[], float b[], 
                  const integer &ldb, integer &info) {
  F77NAME(sgbtrs)(trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
} // gbtrs (float)

inline void gbtrs(const char* trans, const integer &n, const integer &kl,
                  const integer &ku, const integer &nrhs, const double ab[], 
                  const integer &ldab, const integer ipiv[], double b[], 
                  const integer &ldb, integer &info) {
  F77NAME(dgbtrs)(trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
} // gbtrs (double)

#ifdef ARCOMP_H
inline void gbtrs(const char* trans, const integer &n, const integer &kl,
                  const integer &ku, const integer &nrhs, 
                  const arcomplex<float> ab[], const integer &ldab, 
                  const integer ipiv[], arcomplex<float> b[], 
                  const integer &ldb, integer &info) {
  F77NAME(cgbtrs)(trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
} // gbtrs (arcomplex<float>)

inline void gbtrs(const char* trans, const integer &n, const integer &kl,
                  const integer &ku, const integer &nrhs, 
                  const arcomplex<double> ab[], const integer &ldab, 
                  const integer ipiv[], arcomplex<double> b[], 
                  const integer &ldb, integer &info) {
  F77NAME(zgbtrs)(trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
} // gbtrs (arcomplex<double>)
#endif


// GETRF

inline void getrf(const integer &m, const integer &n, float A[], 
                  const integer &lda, integer ipiv[], integer &info) {
  F77NAME(sgetrf)(&m, &n, A, &lda, ipiv, &info);
} // getrf (float)

inline void getrf(const integer &m, const integer &n, double A[], 
                  const integer &lda, integer ipiv[], integer &info) {
  F77NAME(dgetrf)(&m, &n, A, &lda, ipiv, &info);
} // getrf (double)

#ifdef ARCOMP_H
inline void getrf(const integer &m, const integer &n, arcomplex<float> A[], 
                  const integer &lda, integer ipiv[], integer &info) {
  F77NAME(cgetrf)(&m, &n, A, &lda, ipiv, &info);
} // getrf (arcomplex<float>)

inline void getrf(const integer &m, const integer &n, arcomplex<double> A[], 
                  const integer &lda, integer ipiv[], integer &info) {
  F77NAME(zgetrf)(&m, &n, A, &lda, ipiv, &info);
} // getrf (arcomplex<double>)
#endif

// GETRS

inline void getrs(const char* trans, const integer &n, const integer &nrhs,
                  const float A[], const integer &lda,  const integer ipiv[], 
                  float b[], const integer &ldb, integer &info) {
  F77NAME(sgetrs)(trans, &n, &nrhs, A, &lda, ipiv, b, &ldb, &info);
} // getrs (float)

inline void getrs(const char* trans, const integer &n, const integer &nrhs,
                  const double A[], const integer &lda,  const integer ipiv[], 
                  double b[], const integer &ldb, integer &info) {
  F77NAME(dgetrs)(trans, &n, &nrhs, A, &lda, ipiv, b, &ldb, &info);
} // getrs (double)

#ifdef ARCOMP_H
inline void getrs(const char* trans, const integer &n, const integer &nrhs,
                  const arcomplex<float> A[], const integer &lda,  
                  const integer ipiv[], arcomplex<float> b[], 
                  const integer &ldb, integer &info) {
  F77NAME(cgetrs)(trans, &n, &nrhs, A, &lda, ipiv, b, &ldb, &info);
} // getrs (arcomplex<float>)

inline void getrs(const char* trans, const integer &n, const integer &nrhs,
                  const arcomplex<double> A[], const integer &lda,  
                  const integer ipiv[], arcomplex<double> b[], 
                  const integer &ldb, integer &info) {
  F77NAME(zgetrs)(trans, &n, &nrhs, A, &lda, ipiv, b, &ldb, &info);
} // getrs (arcomplex<double>)
#endif

// PTTRF

inline void pttrf(const integer &n, float d[], float e[], integer &info) {
  F77NAME(spttrf)(&n, d, e, &info);
} // pttrf (float)

inline void pttrf(const integer &n, double d[], double e[], integer &info) {
  F77NAME(dpttrf)(&n, d, e, &info);
} // pttrf (double)

// PTTRS

inline void pttrs(const integer &n, const integer &nrhs,
                  const float d[], const float e[], float b[],
                  const integer &ldb, integer &info) {
  F77NAME(spttrs)(&n, &nrhs, d, e, b, &ldb, &info);
} // pttrs (float)

inline void pttrs(const integer &n, const integer &nrhs,
                  const double d[], const double e[], double b[],
                  const integer &ldb, integer &info) {
  F77NAME(dpttrs)(&n, &nrhs, d, e, b, &ldb, &info);
} // pttrs (double)


// SPTRF

inline void sptrf(const char* trans, const integer &n, float ap[], 
                  integer ipiv[], integer &info) {
  F77NAME(ssptrf)(trans, &n, ap, ipiv, &info);
} // sptrf (float)

inline void sptrf(const char* trans, const integer &n, double ap[], 
                  integer ipiv[], integer &info) {
  F77NAME(dsptrf)(trans, &n, ap, ipiv, &info);
} // sptrf (double)


// SPTRS

inline void sptrs(const char* trans, const integer &n, const integer &nrhs,
                  float ap[], integer ipiv[], float b[],
                  const integer &ldb, integer &info) {
  F77NAME(ssptrs)(trans, &n, &nrhs, ap, ipiv, b, &ldb, &info);
} // sptrs (float)

inline void sptrs(const char* trans, const integer &n, const integer &nrhs,
                  double ap[], integer ipiv[], double b[],
                  const integer &ldb, integer &info) {
  F77NAME(dsptrs)(trans, &n, &nrhs, ap, ipiv, b, &ldb, &info);
} // sptrs (double)


inline void second(const float &t) {
  F77NAME(second)(&t);
}

#endif // LAPACKC_H




