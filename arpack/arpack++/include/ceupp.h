/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ceupp.h.
   Interface to ARPACK subroutines zneupd and cneupd.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CEUPP_H
#define CEUPP_H

#include <stddef.h>
#include "arch.h"
#include "arpackf.h"

inline void ceupp(bool rvec, char HowMny, arcomplex<double> d[],
                  arcomplex<double> Z[], int ldz, arcomplex<double> sigma,
                  arcomplex<double> workev[], char bmat, int n, char* which,
                  int nev, double tol, arcomplex<double> resid[], int ncv,
                  arcomplex<double> V[], int ldv, int iparam[], int ipntr[],
                  arcomplex<double> workd[], arcomplex<double> workl[],
                  int lworkl, double rwork[], int& info)

/*
  c++ version of ARPACK routine zneupd.
  This subroutine returns the converged approximations to eigenvalues
  of A*z = lambda*B*z and (optionally):

  (1) the corresponding approximate eigenvectors,
  (2) an orthonormal basis for the associated approximate
      invariant subspace,

  There is negligible additional cost to obtain eigenvectors. An
  orthonormal basis is always computed.  There is an additional storage cost
  of n*nev if both are requested (in this case a separate array Z must be
  supplied).
  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
  are derived from approximate eigenvalues and eigenvectors of
  of the linear operator OP prescribed by the MODE selection in the
  call to caupp. caupp must be called before this routine is called.
  These approximate eigenvalues and vectors are commonly called Ritz
  values and Ritz vectors respectively.  They are referred to as such
  in the comments that follow.  The computed orthonormal basis for the
  invariant subspace corresponding to these Ritz values is referred to
  as a Schur basis.
  See documentation in the header of the subroutine caupp for
  definition of OP as well as other terms and the relation of computed
  Ritz values and Ritz vectors of OP with respect to the given problem
  A*z = lambda*B*z.  For a brief description, see definitions of
  iparam[7], MODE and which in the documentation of caupp.

  Parameters:

    rvec    (Input) Specifies whether a basis for the invariant subspace
            corresponding to the converged Ritz value approximations for
            the eigenproblem A*z = lambda*B*z is computed.
            rvec = false: Compute Ritz values only.
            rvec = true : Compute the Ritz vectors or Schur vectors.
                          See Remarks below.
    HowMny  (Input) Specifies the form of the basis for the invariant
            subspace corresponding to the converged Ritz values that
            is to be computed.
            = 'A': Compute nev Ritz vectors;
            = 'P': Compute nev Schur vectors;
    d       (Output) Array of dimension nev+1. D contains the  Ritz
            approximations to the eigenvalues lambda for A*z = lambda*B*z.
    Z       (Output) Array of dimension nev*n. If rvec = TRUE. and
            HowMny = 'A', then Z contains approximate eigenvectors (Ritz
            vectors) corresponding to the NCONV=iparam[5] Ritz values for
            eigensystem A*z = lambda*B*z.
            If rvec = .FALSE. or HowMny = 'P', then Z is not referenced.
            NOTE: If if rvec = .TRUE. and a Schur basis is not required,
                  the array Z may be set equal to first nev+1 columns of
                  the Arnoldi basis array V computed by caupp.  In this
                  case the Arnoldi basis will be destroyed and overwritten
                  with the eigenvector basis.
    ldz     (Input) Dimension of the vectors contained in Z. This
            parameter MUST be set to n.
    sigma   (Input) If iparam[7] = 3, sigma represents the shift. Not
            referenced if iparam[7] = 1 or 2.
    workv   (Workspace) Array of dimension 2*ncv.
    V       (Input/Output) Array of dimension n*ncv+1.
            Upon Input: V contains the ncv vectors of the Arnoldi basis
                        for OP as constructed by caupp.
            Upon Output: If rvec = TRUE the first NCONV=iparam[5] columns
                        contain approximate Schur vectors that span the
                        desired invariant subspace.
            NOTE: If the array Z has been set equal to first nev+1 columns
                  of the array V and rvec = TRUE. and HowMny = 'A', then
                  the Arnoldi basis held by V has been overwritten by the
                  desired Ritz vectors.  If a separate array Z has been
                  passed then the first NCONV=iparam[5] columns of V will
                  contain approximate Schur vectors that span the desired
                  invariant subspace.
    workl   (Input / Output) Array of length lworkl+1.
            workl[1:ncv*ncv+3*ncv] contains information obtained in
            caupp. They are not changed by ceupp.
            workl[ncv*ncv+3*ncv+1:3*ncv*ncv+4*ncv] holds the untransformed
            Ritz values, the untransformed error estimates of the Ritz
            values, the upper triangular matrix for H, and the associated
            matrix representation of the invariant subspace for H.
    ipntr   (Input / Output) Array of length 14. Pointer to mark the
            starting locations in the workl array for matrices/vectors
            used by caupp and ceupp.
            ipntr[9]:  pointer to the ncv RITZ values of the original
                       system.
            ipntr[11]: pointer to the ncv corresponding error estimates.
            ipntr[12]: pointer to the ncv by ncv upper triangular
                       Schur matrix for H.
            ipntr[13]: pointer to the ncv by ncv matrix of eigenvectors
                       of the upper Hessenberg matrix H. Only referenced
                       by ceupp if rvec = TRUE. See Remark 2 below.
    info    (Output) Error flag.
            =  0 : Normal exit.
            =  1 : The Schur form computed by LAPACK routine csheqr
                   could not be reordered by LAPACK routine ztrsen.
                   Re-enter subroutine ceupp with iparam[5] = ncv and
                   increase the size of the array D to have
                   dimension at least dimension ncv and allocate at least
                   ncv columns for Z. NOTE: Not necessary if Z and V share
                   the same space. Please notify the authors if this error
                   occurs.
            = -1 : n must be positive.
            = -2 : nev must be positive.
            = -3 : ncv must satisfy nev+1 <= ncv <= n.
            = -5 : which must be one of 'LM','SM','LR','SR','LI','SI'.
            = -6 : bmat must be one of 'I' or 'G'.
            = -7 : Length of private work workl array is not sufficient.
            = -8 : Error return from LAPACK eigenvalue calculation.
                   This should never happened.
            = -9 : Error return from calculation of eigenvectors.
                   Informational error from LAPACK routine ztrevc.
            = -10: iparam[7] must be 1, 2 or 3.
            = -11: iparam[7] = 1 and bmat = 'G' are incompatible.
            = -12: HowMny = 'S' not yet implemented.
            = -13: HowMny must be one of 'A' or 'P' if rvec = TRUE.
            = -14: caupp did not find any eigenvalues to sufficient
                   accuracy.

  NOTE:     The following arguments

            bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam,
            ipntr, workd, workl, lworkl, rwork, info

            must be passed directly to ceupp following the last call
            to caupp.  These arguments MUST NOT BE MODIFIED between
            the the last call to caupp and the call to ceupp.

  Remarks
    1. Currently only HowMny = 'A' and 'P' are implemented.
    2. Schur vectors are an orthogonal representation for the basis of
       Ritz vectors. Thus, their numerical properties are often superior.
       Let X' denote the transpose of X. If rvec = .TRUE. then the
       relationship A * V[:,1:iparam[5]] = V[:,1:iparam[5]] * T, and
       V[:,1:iparam[5]]' * V[:,1:iparam[5]] = I are approximately satisfied.
       Here T is the leading submatrix of order iparam[5] of the real
       upper quasi-triangular matrix stored workl[ipntr[12]].
*/

{

  int                irvec;
  logical*           iselect;
  arcomplex<double>* iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[1] : Z;

  F77NAME(zneupd)(&irvec, &HowMny, iselect, d, iZ, &ldz, &sigma,
                  &workev[1], &bmat, &n, which, &nev, &tol, resid,
                  &ncv, &V[1], &ldv, &iparam[1], &ipntr[1],
                  &workd[1], &workl[1], &lworkl, &rwork[1], &info);

  delete[] iselect;

} // ceupp (arcomplex<double>).

inline void ceupp(bool rvec, char HowMny, arcomplex<float> d[],
                  arcomplex<float> Z[], int ldz, arcomplex<float> sigma,
                  arcomplex<float> workev[], char bmat, int n, char* which,
                  int nev, float tol, arcomplex<float> resid[], int ncv,
                  arcomplex<float> V[], int ldv, int iparam[], int ipntr[],
                  arcomplex<float> workd[], arcomplex<float> workl[],
                  int lworkl, float rwork[], int& info)

/*
  c++ version of ARPACK routine cneupd. The only difference between
  cneupd and zneupd is that in the former function all vectors have
  single precision elements and in the latter all vectors have double
  precision elements.
*/

{

  int               irvec;
  logical*          iselect;
  arcomplex<float>* iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[1] : Z;

  F77NAME(cneupd)(&irvec, &HowMny, iselect, d, iZ, &ldz, &sigma,
                  &workev[1], &bmat, &n, which, &nev, &tol, resid,
                  &ncv, &V[1], &ldv, &iparam[1], &ipntr[1],
                  &workd[1], &workl[1], &lworkl, &rwork[1], &info);

  delete[] iselect;

} // ceupp (arcomplex<float>).

#endif // CEUPP_H
