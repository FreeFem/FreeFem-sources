/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE neupp.h.
   Interface to ARPACK subroutines dneupd and sneupd.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef NEUPP_H
#define NEUPP_H

#include <stddef.h>
#include "arch.h"
#include "arpackf.h"

inline void neupp(bool rvec, char HowMny, double dr[],
                  double di[], double Z[], int ldz, double sigmar,
                  double sigmai, double workv[], char bmat, int n,
                  char* which, int nev, double tol, double resid[],
                  int ncv, double V[], int ldv, int iparam[],
                  int ipntr[], double workd[], double workl[],
                  int lworkl, int& info)

/*
  c++ version of ARPACK routine dneupd.
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
  call to naupp. naupp must be called before this routine is called.
  These approximate eigenvalues and vectors are commonly called Ritz
  values and Ritz vectors respectively.  They are referred to as such
  in the comments that follow.  The computed orthonormal basis for the
  invariant subspace corresponding to these Ritz values is referred to
  as a Schur basis.
  See documentation in the header of the subroutine naupp for
  definition of OP as well as other terms and the relation of computed
  Ritz values and Ritz vectors of OP with respect to the given problem
  A*z = lambda*B*z. For a brief description, see definitions of
  iparam[7], MODE and which in the documentation of naupp.

  Parameters:

    rvec    (Input) Specifies whether Ritz vectors corresponding to the
            Ritz value approximations to the eigenproblem A*z = lambda*B*z
            are computed.
            rvec = false: Compute Ritz values only.
            rvec = true : Compute the Ritz vectors or Schur vectors.
                          See Remarks below.
    HowMny  (Input) Specifies the form of the basis for the invariant
            subspace corresponding to the converged Ritz values that
            is to be computed.
            = 'A': Compute nev Ritz vectors;
            = 'P': Compute nev Schur vectors;
    dr      (Output) Array of dimension nev+1.
            If iparam[7] = 1,2 or 3 and sigmai=0.0  then on exit: dr
            contains the real part of the Ritz  approximations to the
            eigenvalues of A*z = lambda*B*z.
            If iparam[7] = 3, 4 and sigmai is not equal to zero, then on
            exit: dr contains the real part of the Ritz values of OP
            computed by naupp. A further computation must be performed by
            the user to transform the Ritz values computed for OP by naupp
            to those of the original system A*z = lambda*B*z. See remark 3.
    di      (Output) Array of dimension nev+1.
            On exit, di contains the imaginary part of the Ritz value
            approximations to the eigenvalues of A*z = lambda*B*z
            associated with dr.
            NOTE: When Ritz values are complex, they will come in complex
                  conjugate pairs.  If eigenvectors are requested, the
                  corresponding Ritz vectors will also come in conjugate
                  pairs and the real and imaginary parts of these are
                  represented in two consecutive columns of the array Z
                  (see below).
    Z       (Output) Array of dimension nev*n if rvec = TRUE and HowMny =
            'A'.  if rvec = TRUE. and HowMny = 'A', then the contains
            approximate eigenvectors (Ritz vectors) corresponding to the
            NCONV=iparam[5] Ritz values for eigensystem A*z = lambda*B*z.
            The complex Ritz vector associated with the Ritz value
            with positive imaginary part is stored in two consecutive
            columns.  The first column holds the real part of the Ritz
            vector and the second column holds the imaginary part.  The
            Ritz vector associated with the Ritz value with negative
            imaginary part is simply the complex conjugate of the Ritz
            vector associated with the positive imaginary part.
            If rvec = .FALSE. or HowMny = 'P', then Z is not referenced.
            NOTE: If if rvec = .TRUE. and a Schur basis is not required,
                  the array Z may be set equal to first nev+1 columns of
                  the Arnoldi basis array V computed by naupp.  In this
                  case the Arnoldi basis will be destroyed and overwritten
                  with the eigenvector basis.
    ldz     (Input) Dimension of the vectors contained in Z. This
            parameter MUST be set to n.
    sigmar  (Input) If iparam[7] = 3 or 4, represents the real part of
            the shift. Not referenced if iparam[7] = 1 or 2.
    sigmai  (Input) If iparam[7] = 3 or 4, represents the imaginary part
            of the shift. Not referenced if iparam[7] = 1 or 2. See
            remark 3 below.
    workv   (Workspace) Array of dimension 3*ncv.
    V       (Input/Output) Array of dimension n*ncv+1.
            Upon Input: V contains the ncv vectors of the Arnoldi basis
                        for OP as constructed by naupp.
            Upon Output: If rvec = TRUE the first NCONV=iparam[5] columns
                        contain approximate Schur vectors that span the
                        desired invariant subspace.  See Remark 2 below.
            NOTE: If the array Z has been set equal to first nev+1 columns
                  of the array V and rvec = TRUE. and HowMny = 'A', then
                  the Arnoldi basis held by V has been overwritten by the
                  desired Ritz vectors.  If a separate array Z has been
                  passed then the first NCONV=iparam[5] columns of V will
                  contain approximate Schur vectors that span the desired
                  invariant subspace.
    workl   (Input / Output) Array of length lworkl+1.
            workl[1:ncv*ncv+3*ncv] contains information obtained in
            naupp. They are not changed by neupp.
            workl[ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv] holds the real and
            imaginary part of the untransformed Ritz values, the upper
            quasi-triangular matrix for H, and the associated matrix
            representation of the invariant subspace for H.
    ipntr   (Input / Output) Array of length 14. Pointer to mark the
            starting locations in the workl array for matrices/vectors
            used by naupp and neupp.
            ipntr[9]:  pointer to the real part of the ncv RITZ values
                       of the original system.
            ipntr[10]: pointer to the imaginary part of the ncv RITZ
                       values of the original system.
            ipntr[11]: pointer to the ncv corresponding error bounds.
            ipntr[12]: pointer to the ncv by ncv upper quasi-triangular
                       Schur matrix for H.
            ipntr[13]: pointer to the ncv by ncv matrix of eigenvectors
                       of the upper Hessenberg matrix H. Only referenced
                       by neupp if rvec = TRUE. See Remark 2 below.
    info    (Output) Error flag.
            =  0 : Normal exit.
            =  1 : The Schur form computed by LAPACK routine dlahqr
                   could not be reordered by LAPACK routine dtrsen.
                   Re-enter subroutine neupp with iparam[5] = ncv and
                   increase the size of the arrays DR and DI to have
                   dimension at least dimension ncv and allocate at least
                   ncv columns for Z. NOTE: Not necessary if Z and V share
                   the same space. Please notify the authors if this error
                   occurs.
            = -1 : n must be positive.
            = -2 : nev must be positive.
            = -3 : ncv must satisfy nev+2 <= ncv <= n.
            = -5 : which must be one of 'LM','SM','LR','SR','LI','SI'.
            = -6 : bmat must be one of 'I' or 'G'.
            = -7 : Length of private work workl array is not sufficient.
            = -8 : Error return from calculation of a real Schur form.
                   Informational error from LAPACK routine dlahqr.
            = -9 : Error return from calculation of eigenvectors.
                   Informational error from LAPACK routine dtrevc.
            = -10: iparam[7] must be 1,2,3,4.
            = -11: iparam[7] = 1 and bmat = 'G' are incompatible.
            = -12: HowMny = 'S' not yet implemented
            = -13: HowMny must be one of 'A' or 'P' if rvec = TRUE.
            = -14: naupp did not find any eigenvalues to sufficient
                   accuracy.

  NOTE:     The following arguments

            bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam,
            ipntr, workd, workl, lworkl, info

            must be passed directly to neupp following the last call
            to naupp.  These arguments MUST NOT BE MODIFIED between
            the the last call to naupp and the call to neupp.

  Remarks
    1. Currently only HowMny = 'A' and 'P' are implemented.
    2. Schur vectors are an orthogonal representation for the basis of
       Ritz vectors. Thus, their numerical properties are often superior.
       Let X' denote the transpose of X. If rvec = .TRUE. then the
       relationship A * V[:,1:iparam[5]] = V[:,1:iparam[5]] * T, and
       V[:,1:iparam[5]]' * V[:,1:iparam[5]] = I are approximately satisfied.
       Here T is the leading submatrix of order iparam[5] of the real
       upper quasi-triangular matrix stored workl[ipntr[12]]. That is,
       T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
       each 2-by-2 diagonal block has its diagonal elements equal and its
       off-diagonal elements of opposite sign.  Corresponding to each
       2-by-2 diagonal block is a complex conjugate pair of Ritz values.
       The real Ritz values are stored on the diagonal of T.
    3. If iparam[7] = 3 or 4 and sigmai is not equal zero, then the user
       must form the iparam[5] Rayleigh quotients in order to transform the
       Ritz values computed by naupp for OP to those of A*z = lambda*B*z.
       Set rvec = TRUE. and HowMny = 'A', and compute
       Z[:,I]' * A * Z[:,I] if di[I] = 0.
       If di[I] is not equal to zero and di[I+1] = - D[I],
       then the desired real and imaginary parts of the Ritz value are
       Z[:,I]' * A * Z[:,I] +  Z[:,I+1]' * A * Z[:,I+1],
       Z[:,I]' * A * Z[:,I+1] -  Z[:,I+1]' * A * Z[:,I], respectively.
       Another possibility is to set rvec = .true. and HowMny = 'P' and
       compute V[:,1:iparam[5]]' * A * V[:,1:iparam[5]] and then an upper
       quasi-triangular matrix of order iparam[5] is computed. See remark
       2 above.
*/

{

  int      irvec;
  logical* iselect;
  double*  iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[1] : Z;

  F77NAME(dneupd)(&irvec, &HowMny, iselect, dr, di, iZ, &ldz, &sigmar,
                  &sigmai, &workv[1], &bmat, &n, which, &nev, &tol,
                  resid, &ncv, &V[1], &ldv, &iparam[1], &ipntr[1],
                  &workd[1], &workl[1], &lworkl, &info);

  delete[] iselect;

} // neupp (double).

inline void neupp(bool rvec, char HowMny, float dr[],
                  float di[], float Z[], int ldz, float sigmar,
                  float sigmai, float workv[], char bmat, int n,
                  char* which, int nev, float tol, float resid[],
                  int ncv, float V[], int ldv, int iparam[],
                  int ipntr[], float workd[], float workl[],
                  int lworkl, int& info)

/*
  c++ version of ARPACK routine sneupd. The only difference between
  sneupd and dneupd is that in the former function all vectors have
  single precision elements and in the latter all vectors have double
  precision elements.
*/

{

  int      irvec;
  logical* iselect;
  float*   iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[1] : Z;

  F77NAME(sneupd)(&irvec, &HowMny, iselect, dr, di, iZ, &ldz, &sigmar,
                  &sigmai, &workv[1], &bmat, &n, which, &nev, &tol,
                  resid, &ncv, &V[1], &ldv, &iparam[1], &ipntr[1],
                  &workd[1], &workl[1], &lworkl, &info );

  delete[] iselect;

} // neupp (float).

#endif // NEUPP_H

