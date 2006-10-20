/*
  ARPACK++ v1.0 8/1/1997
  c++ interface to ARPACK code.

  MODULE arpackf.h
  ARPACK FORTRAN routines.

  ARPACK Authors
     Richard Lehoucq
     Danny Sorensen
     Chao Yang
     Dept. of Computational & Applied Mathematics
     Rice University
     Houston, Texas
*/

#ifndef ARPACKF_H
#define ARPACKF_H

#include "arch.h"


#define HIDDEN_HBW ,int,int,int
#define HIDDEN_BW ,int,int
#define HIDDEN_12 ,1,2
#define HIDDEN_112 ,1,1,2


extern "C"
{

// debug "common" statement.

  struct { 
    integer logfil, ndigit, mgetv0;
    integer msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
    integer mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd;
    integer mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd;
  } F77NAME(debug);

  // add FH pour lib dynamic so win 32 
  struct { 
    integer           nopx, nbx, nrorth, nitref, nrstrt;
    float           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
      tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
      tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
      tmvopx, tmvbx, tgetv0, titref, trvec;
    
    
  } F77NAME(timing);


// double precision symmetric routines.

  void F77NAME(dsaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr, double *workd,
                       double *workl, integer *lworkl, integer *info HIDDEN_BW);

  void F77NAME(dseupd)(logical *rvec, char *HowMny, logical *select,
                       double *d, double *Z, integer *ldz,
                       double *sigma, char *bmat, integer *n,
                       char *which, integer *nev, double *tol,
                       double *resid, integer *ncv, double *V,
                       integer *ldv, integer *iparam, integer *ipntr,
                       double *workd, double *workl,
                       integer *lworkl, integer *info HIDDEN_HBW );

// double precision nonsymmetric routines.

  void F77NAME(dnaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr, double *workd,
                       double *workl, integer *lworkl, integer *info HIDDEN_BW);

  void F77NAME(dneupd)(logical *rvec, char *HowMny, logical *select,
                       double *dr, double *di, double *Z,
                       integer *ldz, double *sigmar,
                       double *sigmai, double *workev,
                       char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr,
                       double *workd, double *workl,
                       integer *lworkl, integer *info HIDDEN_HBW);

// single precision symmetric routines.

  void F77NAME(ssaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, float *resid,
                       integer *ncv, float *V, integer *ldv,
                       integer *iparam, integer *ipntr, float *workd,
                       float *workl, integer *lworkl, integer *info HIDDEN_BW);

  void F77NAME(sseupd)(logical *rvec, char *HowMny, logical *select,
                       float *d, float *Z, integer *ldz,
                       float *sigma, char *bmat, integer *n,
                       char *which, integer *nev, float *tol,
                       float *resid, integer *ncv, float *V,
                       integer *ldv, integer *iparam, integer *ipntr,
                       float *workd, float *workl,
                       integer *lworkl, integer *info HIDDEN_HBW);

// single precision nonsymmetric routines.

  void F77NAME(snaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, float *resid,
                       integer *ncv, float *V, integer *ldv,
                       integer *iparam, integer *ipntr, float *workd,
                       float *workl, integer *lworkl, integer *info HIDDEN_BW);

  void F77NAME(sneupd)(logical *rvec, char *HowMny, logical *select,
                       float *dr, float *di, float *Z,
                       integer *ldz, float *sigmar,
                       float *sigmai, float *workev, char *bmat,
                       integer *n, char *which, integer *nev,
                       float *tol, float *resid, integer *ncv,
                       float *V, integer *ldv, integer *iparam,
                       integer *ipntr, float *workd, float *workl,
                       integer *lworkl, integer *info HIDDEN_HBW);

#ifdef ARCOMP_H

// single precision complex routines.

  void F77NAME(cnaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, arcomplex<float> *resid,
                       integer *ncv, arcomplex<float> *V, integer *ldv,
                       integer *iparam, integer *ipntr, arcomplex<float> *workd,
                       arcomplex<float> *workl, integer *lworkl,
                       float *rwork, integer *info HIDDEN_BW);

  void F77NAME(cneupd)(logical *rvec, char *HowMny, logical *select,
                       arcomplex<float> *d, arcomplex<float> *Z, integer *ldz,
                       arcomplex<float> *sigma, arcomplex<float> *workev,
                       char *bmat, integer *n, char *which, integer *nev,
                       float *tol, arcomplex<float> *resid, integer *ncv,
                       arcomplex<float> *V, integer *ldv, integer *iparam,
                       integer *ipntr, arcomplex<float> *workd,
                       arcomplex<float> *workl, integer *lworkl,
                       float *rwork, integer *info HIDDEN_HBW);

// double precision complex routines.

  void F77NAME(znaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, arcomplex<double> *resid,
                       integer *ncv, arcomplex<double> *V, integer *ldv,
                       integer *iparam, integer *ipntr, arcomplex<double> *workd,
                       arcomplex<double> *workl, integer *lworkl,
                       double *rwork, integer *info HIDDEN_BW);

  void F77NAME(zneupd)(logical *rvec, char *HowMny, logical *select,
                       arcomplex<double> *d, arcomplex<double> *Z, integer *ldz,
                       arcomplex<double> *sigma, arcomplex<double> *workev,
                       char *bmat, integer *n, char *which, integer *nev,
                       double *tol, arcomplex<double> *resid, integer *ncv,
                       arcomplex<double> *V, integer *ldv, integer *iparam,
                       integer *ipntr, arcomplex<double> *workd,
                       arcomplex<double> *workl, integer *lworkl,
                       double *rwork, integer *info HIDDEN_HBW);

}

#endif // ARCOMP_H

#endif // ARPACKF_H
