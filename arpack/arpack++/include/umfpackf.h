/*
  ARPACK++ v1.0 8/1/1997
  c++ interface to ARPACK code.

  MODULE umfpackf.h
  UMFPACK FORTRAN routines.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas
*/

#ifndef UMFPACKF_H
#define UMFPACKF_H

#include "arch.h"

extern "C"
{

  // Single precision real routines.

  void F77NAME(ums21i)(integer *keep, float *cntl, integer *icntl);

  void F77NAME(ums2fa)(const integer *n, const integer *ne, 
                       const integer *job, const logical *transa,
                       const integer *lvalue, const integer *lindex,
                       float *value, integer *index, integer *keep,
                       const float *cntl, const integer *icntl,
                       integer *info, float *rinfo); 

  void F77NAME(ums2so)(const integer *n, const integer *job,
                       const logical *transc, const integer *lvalue,
                       const integer *lindex, float *value,
                       integer *index, const integer *keep, 
                       const float *b, float *x, float *w,
                       const float *cntl, const integer *icntl,
                       integer *info, float *rinfo);


  // Double precision real routines.

  void F77NAME(umd21i)(integer *keep, double *cntl, integer *icntl);

  void F77NAME(umd2fa)(const integer *n, const integer *ne, 
                       const integer *job, const logical *transa,
                       const integer *lvalue, const integer *lindex,
                       double *value, integer *index, integer *keep,
                       const double *cntl, const integer *icntl,
                       integer *info, double *rinfo); 

  void F77NAME(umd2so)(const integer *n, const integer *job,
                       const logical *transc, const integer *lvalue,
                       const integer *lindex, double *value,
                       integer *index, const integer *keep, 
                       const double *b, double *x, double *w,
                       const double *cntl, const integer *icntl,
                       integer *info, double *rinfo);


  // Single precision complex routines.

#ifdef ARCOMP_H

  void F77NAME(umc21i)(integer *keep, float *cntl, integer *icntl);

  void F77NAME(umc2fa)(const integer *n, const integer *ne, 
                       const integer *job, const logical *transa,
                       const integer *lvalue, const integer *lindex,
                       arcomplex<float> *value, integer *index, integer *keep,
                       const float *cntl, const integer *icntl,
                       integer *info, float *rinfo); 

  void F77NAME(umc2so)(const integer *n, const integer *job,
                       const logical *transc, const integer *lvalue,
                       const integer *lindex, arcomplex<float> *value,
                       integer *index, const integer *keep, 
                       const arcomplex<float> *b, arcomplex<float> *x, 
                       arcomplex<float> *w, const float *cntl, 
                       const integer *icntl, integer *info, float *rinfo);


  // Double precision complex routines.

  void F77NAME(umz21i)(integer *keep, double *cntl, integer *icntl);

  void F77NAME(umz2fa)(const integer *n, const integer *ne, 
                       const integer *job, const logical *transa,
                       const integer *lvalue, const integer *lindex,
                       arcomplex<double> *value, integer *index, integer *keep,
                       const double *cntl, const integer *icntl,
                       integer *info, double *rinfo); 

  void F77NAME(umz2so)(const integer *n, const integer *job,
                       const logical *transc, const integer *lvalue,
                       const integer *lindex, arcomplex<double> *value,
                       integer *index, const integer *keep, 
                       const arcomplex<double> *b, arcomplex<double> *x, 
                       arcomplex<double> *w, const double *cntl, 
                       const integer *icntl, integer *info, double *rinfo);

#endif // ARCOMP_H

}
#endif // UMFPACKF_H

