/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE lapackf.h.
   Interface to LAPACK FORTRAN routines.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef LAPACKF_H
#define LAPACKF_H

#include "arch.h"

extern "C"
{

  // Single precision real routines.

  float F77NAME(slapy2)(const float *x, const float *y);

  void F77NAME(slacpy)(const char* uplo, const integer *m, const integer *n,
                       const float *a, const integer *lda, float *b, 
                       const integer *ldb);

  void F77NAME(sgttrf)(const integer *n, float *dl, float *d, float *du,
                       float *du2, integer *ipiv, integer *info);

  void F77NAME(sgbtrf)(const integer *m, const integer *n, const integer *kl,
                       const integer *ku, float *ab, const integer *ldab,
                       integer *ipiv, integer *info);

  void F77NAME(sgetrf)(const integer *m, const integer *n, float *A,
                       const integer *lda, integer *ipiv, integer *info);

  void F77NAME(sgttrs)(const char* trans, const integer *n,
                       const integer *nrhs, const float *dl,
                       const float *d, const float *du,
                       const float *du2, const integer *ipiv,
                       float* b, const integer *ldb, integer *info);

  void F77NAME(sgbtrs)(const char* trans, const integer *n, 
                       const integer *kl, const integer *ku, 
                       const integer *nrhs, const float *ab, 
                       const integer *ldab, const integer *ipiv, 
                       float *b, const integer *ldb, integer *info);

  void F77NAME(sgetrs)(const char* trans, const integer *n,
                       const integer *nrhs, const float *A,
                       const integer *lda, const integer *ipiv,
                       float* b, const integer *ldb, integer *info);

  void F77NAME(spttrf)(const integer *n, float *d, float *e, integer *info);

  void F77NAME(spttrs)(const integer *n, const integer *nrhs,
                       const float *d, const float *e, float *b,
                       const integer *ldb, integer *info);

  void F77NAME(ssptrf)(const char* trans, const integer *n, 
                       float *ap, integer *ipiv, integer *info);

  void F77NAME(ssptrs)(const char* trans, const integer *n, 
                       const integer *nrhs, float *ap, integer *ipiv, 
                       float *b, const integer *ldb, integer *info);

  // Double precision real routines.

  double F77NAME(dlapy2)(const double *x, const double *y);

  void F77NAME(dlacpy)(const char* uplo, const integer *m, const integer *n,
                       const double *a, const integer *lda, double *b, 
                       const integer *ldb);

  void F77NAME(dgttrf)(const integer *n, double *dl, double *d, double *du,
                       double *du2, integer *ipiv, integer *info);

  void F77NAME(dgbtrf)(const integer *m, const integer *n, const integer *kl,
                       const integer *ku, double *ab, const integer *ldab,
                       integer *ipiv, integer *info);

  void F77NAME(dgetrf)(const integer *m, const integer *n, double *A,
                       const integer *lda, integer *ipiv, integer *info);

  void F77NAME(dgttrs)(const char* trans, const integer *n,
                       const integer *nrhs, const double *dl,
                       const double *d, const double *du,
                       const double *du2, const integer *ipiv,
                       double* b, const integer *ldb, integer *info);

  void F77NAME(dgbtrs)(const char* trans, const integer *n, 
                       const integer *kl, const integer *ku, 
                       const integer *nrhs, const double *ab, 
                       const integer *ldab, const integer *ipiv, 
                       double *b, const integer *ldb, integer *info);

  void F77NAME(dgetrs)(const char* trans, const integer *n,
                       const integer *nrhs, const double *A,
                       const integer *lda, const integer *ipiv,
                       double* b, const integer *ldb, integer *info);

  void F77NAME(dpttrf)(const integer *n, double *d, double *e, integer *info);

  void F77NAME(dpttrs)(const integer *n, const integer *nrhs,
                       const double *d, const double *e, double *b,
                       const integer *ldb, integer *info);

  void F77NAME(dsptrf)(const char* trans, const integer *n, 
                       double *ap, integer *ipiv, integer *info);

  void F77NAME(dsptrs)(const char* trans, const integer *n, 
                       const integer *nrhs, double *ap, integer *ipiv, 
                       double *b, const integer *ldb, integer *info);

#ifdef ARCOMP_H

  // Single precision complex routines.

  void F77NAME(clacpy)(const char* uplo, const integer *m, const integer *n,
                       const arcomplex<float> *a, const integer *lda, 
                       arcomplex<float> *b, const integer *ldb);

  void F77NAME(cgttrf)(const integer *n, arcomplex<float> *dl,
                       arcomplex<float> *d, arcomplex<float> *du,
                       arcomplex<float> *du2, integer *ipiv,
                       integer *info);

  void F77NAME(cgbtrf)(const integer *m, const integer *n, const integer *kl,
                       const integer *ku, arcomplex<float> *ab, 
                       const integer *ldab, integer *ipiv, integer *info);

  void F77NAME(cgetrf)(const integer *m, const integer *n, arcomplex<float> *A,
                       const integer *lda, integer *ipiv, integer *info);

  void F77NAME(cgttrs)(const char *trans, const integer *n,
                       const integer *nrhs, const arcomplex<float> *dl,
                       const arcomplex<float> *d, const arcomplex<float> *du,
                       const arcomplex<float> *du2, const integer *ipiv,
                       arcomplex<float>* b, const integer *ldb,
                       integer *info);

  void F77NAME(cgbtrs)(const char* trans, const integer *n, 
                       const integer *kl, const integer *ku, 
                       const integer *nrhs, const arcomplex<float> *ab, 
                       const integer *ldab, const integer *ipiv, 
                       arcomplex<float> *b, const integer *ldb, 
                       integer *info);

  void F77NAME(cgetrs)(const char* trans, const integer *n,
                       const integer *nrhs, const arcomplex<float> *A,
                       const integer *lda, const integer *ipiv,
                       arcomplex<float>* b, const integer *ldb, integer *info);

  // Double precision complex routines.

  void F77NAME(zlacpy)(const char* uplo, const integer *m, const integer *n,
                       const arcomplex<double> *a, const integer *lda, 
                       arcomplex<double> *b, const integer *ldb);

  void F77NAME(zgttrf)(const integer *n, arcomplex<double> *dl,
                       arcomplex<double> *d, arcomplex<double> *du,
                       arcomplex<double> *du2, integer *ipiv,
                       integer *info);

  void F77NAME(zgbtrf)(const integer *m, const integer *n, const integer *kl,
                       const integer *ku, arcomplex<double> *ab, 
                       const integer *ldab, integer *ipiv, integer *info);

  void F77NAME(zgetrf)(const integer *m, const integer *n, arcomplex<double> *A,
                       const integer *lda, integer *ipiv, integer *info);

  void F77NAME(zgttrs)(const char *trans, const integer *n,
                       const integer *nrhs, const arcomplex<double> *dl,
                       const arcomplex<double> *d, const arcomplex<double> *du,
                       const arcomplex<double> *du2, const integer *ipiv,
                       arcomplex<double>* b, const integer *ldb,
                       integer *info);

  void F77NAME(zgbtrs)(const char* trans, const integer *n, 
                       const integer *kl, const integer *ku, 
                       const integer *nrhs, const arcomplex<double> *ab, 
                       const integer *ldab, const integer *ipiv, 
                       arcomplex<double> *b, const integer *ldb, 
                       integer *info);

  void F77NAME(zgetrs)(const char* trans, const integer *n,
                       const integer *nrhs, const arcomplex<double> *A,
                       const integer *lda, const integer *ipiv,
                       arcomplex<double>* b, const integer *ldb, integer *info);

#endif // ARCOMP_H

  void F77NAME(second)(const float *T);

}
#endif // LAPACKF_H





