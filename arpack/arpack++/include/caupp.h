/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE caupp.h.
   Interface to ARPACK subroutines znaupd and cnaupd.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef CAUPP_H
#define CAUPP_H

#include "arch.h"
#include "arpackf.h"

inline void caupp(int& ido, char bmat, int n, char* which, int nev,
                  double& tol, arcomplex<double> resid[], int ncv,
                  arcomplex<double> V[], int ldv, int iparam[], int ipntr[],
                  arcomplex<double> workd[], arcomplex<double> workl[],
                  int lworkl, double rwork[], int& info)

/*
  c++ version of ARPACK routine znaupd that implements the
  Reverse communication interface for the Implicitly Restarted Arnoldi
  iteration. This is intended to be used to find a few eigenpairs of a
  complex linear operator OP with respect to a semi-inner product defined
  by a hermitian positive semi-definite real matrix B. B may be the
  identity matrix. NOTE: if both OP and B are real, then naupp should
  be used.

  The computed approximate eigenvalues are called Ritz values and
  the corresponding approximate eigenvectors are called Ritz vectors.

  caupp is usually called iteratively to solve one of the
  following problems:

  Mode 1:  A*x = lambda*x.
           ===> OP = A  and  B = I.

  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
           ===> OP = inv[M]*A  and  B = M.
           ===> (If M can be factored see remark 3 below)

  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
           ===> OP =  inv[A - sigma*M]*M   and  B = M.
           ===> shift-and-invert mode
           If OP*x = amu*x, then lambda = sigma + 1/amu.


  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
        should be accomplished either by a direct method
        using a sparse matrix factorization and solving

           [A - sigma*M]*w = v  or M*w = v,

        or through an iterative method for solving these systems. If
        an iterative method is used, the convergence test must be more
        stringent than the accuracy requirements for the eigenvalue
        approximations.

  Parameters:

    ido     (Input / Output) Reverse communication flag.  ido must be
            zero on the first call to caupp.  ido will be set
            internally to indicate the type of operation to be
            performed.  Control is then given back to the calling
            routine which has the responsibility to carry out the
            requested operation and call caupp with the result. The
            operand is given in workd[ipntr[1]], the result must be
            put in workd[ipntr[2]].
            ido =  0: first call to the reverse communication interface.
            ido = -1: compute  Y = OP * X  where
                      ipntr[1] is the pointer into workd for X,
                      ipntr[2] is the pointer into workd for Y.
                      This is for the initialization phase to force the
                      starting vector into the range of OP.
            ido =  1: compute  Y = OP * X where
                      ipntr[1] is the pointer into workd for X,
                      ipntr[2] is the pointer into workd for Y.
                      In mode 3 and 4, the vector B * X is already
                      available in workd[ipntr[3]].  It does not
                      need to be recomputed in forming OP * X.
            ido =  2: compute  Y = B * X  where
                      ipntr[1] is the pointer into workd for X,
                      ipntr[2] is the pointer into workd for Y.
            ido =  3: compute the iparam[8] real and imaginary parts
                      of the shifts where inptr[14] is the pointer
                      into workl for placing the shifts. See Remark 
                      5 below.
            ido = 99: done.
    bmat    (Input) bmat specifies the type of the matrix B that defines
            the semi-inner product for the operator OP.
            bmat = 'I' -> standard eigenvalue problem A*x = lambda*x;
            bmat = 'G' -> generalized eigenvalue problem A*x = lambda*M*x.
    n       (Input) Dimension of the eigenproblem.
    nev     (Input) Number of eigenvalues to be computed. 0 < nev <= n-1.
    which   (Input) Specify which of the Ritz values of OP to compute.
            'LM' - compute the nev eigenvalues of largest magnitude.
            'SM' - compute the nev eigenvalues of smallest magnitude.
            'LR' - compute the nev eigenvalues of largest real part.
            'SR' - compute the nev eigenvalues of smallest real part.
            'LI' - compute the nev eigenvalues of largest imaginary part.
            'SI' - compute the nev eigenvalues of smallest imaginary part.
    tol     (Input) Stopping criterion: the relative accuracy of the
            Ritz value is considered acceptable if BOUNDS[i] <=
            tol*abs(RITZ[i]),where ABS(RITZ[i]) is the magnitude when
            RITZ[i] is complex. If tol<=0.0 is passed, the machine
            precision as computed by the LAPACK auxiliary subroutine
            _LAMCH is used.
    resid   (Input / Output) Array of length n.
            On input:
            If info==0, a random initial residual vector is used.
            If info!=0, resid contains the initial residual vector,
                        possibly from a previous run.
            On output:
            resid contains the final residual vector.
    ncv     (Input) Number of Arnoldi vectors that are generated at each
            iteration. After the startup phase in which nev Arnoldi
            vectors are generated, the algorithm generates ncv-nev
            Arnoldi vectors at each subsequent update iteration. Most of
            the cost in generating each Arnoldi vector is in the
            matrix-vector product OP*x.
            NOTE: ncv must satisfy nev+1 <= ncv <= n.
    V       (Output) Array of length ncv*n+1. V contains the ncv Arnoldi
            basis vectors. The first element V[0] is never referenced.
    ldv     (Input) Dimension of the basis vectors contained in V. This
            parameter MUST be set to n.
    iparam  (Input / Output) Array of length 12.
            iparam[1]  = ISHIFT: method for selecting the implicit shifts.
            The shifts selected at each iteration are used to restart
            the Arnoldi iteration in an implicit fashion.
            -------------------------------------------------------------
            ISHIFT = 0: the shifts are to be provided by the user via
                        reverse communication.  The ncv eigenvalues of
                        the Hessenberg matrix H are returned in the part
                        of workl array corresponding to RITZ.
            ISHIFT = 1: exact shifts with respect to the current
                        Hessenberg matrix H.  This is equivalent to
                        restarting the iteration from the beginning
                        after updating the starting vector with a linear
                        combination of Ritz vectors associated with the
                        "wanted" eigenvalues.
            ISHIFT = 2: other choice of internal shift to be defined.
            -------------------------------------------------------------
            iparam[2]  is no longer referenced.
            iparam[3]  = MXITER
            On INPUT:  maximum number of Arnoldi update iterations allowed.
            On OUTPUT: actual number of Arnoldi update iterations taken.
            iparam[4]  = NB: blocksize to be used in the recurrence.
            The code currently works only for NB = 1.
            iparam[5]  = NCONV: number of "converged" Ritz values.
            This represents the number of Ritz values that satisfy
            the convergence criterion.
            iparam[6]  is no longer referenced.
            iparam[7]  = MODE. On input determines what type of
            eigenproblem is being solved. Must be 1, 2 or 3.
            iparam[8]  = NP. When ido = 3 and the user provides shifts
            through reverse communication (iparam[1]=0), caupp returns
            NP, the number of shifts the user is to provide.
            0 < NP <=ncv-nev. See Remark 5 below.
            iparam[9]  =  total number of OP*x operations.
            iparam[10] = total number of B*x operations if bmat='G'.
            iparam[11] = total number of steps of re-orthogonalization.
    ipntr   (Output) Array of length 15. Pointer to mark the starting
            locations in the workd and workl arrays for matrices/vectors
            used by the Arnoldi iteration.
            ipntr[1] : pointer to the current operand vector X in workd.
            ipntr[2] : pointer to the current result vector Y in workd.
            ipntr[3] : pointer to the vector B * X in workd when used in
                       the shift-and-invert mode.
            ipntr[4] : pointer to the next available location in workl
                       that is untouched by the program.
            ipntr[5] : pointer to the ncv by ncv upper Hessenberg matrix
                       H in workl.
            ipntr[6] : pointer to the ritz value array RITZ.
            ipntr[7] : pointer to the (projected) ritz vector array Q.
            ipntr[8] : pointer to the error BOUNDS array in workl.
            ipntr[14]: pointer to the NP shifts in workl. See Remark 5.
            Note: ipntr[9:13] is only referenced by ceupp. See Remark 2.
            ipntr[9] : pointer to the ncv RITZ values of the
                       original system.
            ipntr[10]: Not Used
            ipntr[11]: pointer to the ncv corresponding error bounds.
            ipntr[12]: pointer to the ncv by ncv upper triangular
                       Schur matrix for H.
            ipntr[13]: pointer to the ncv by ncv matrix of eigenvectors
                       of the upper Hessenberg matrix H. Only referenced by
                       ceupp if RVEC = true. See Remark 2 below.
    workd   (Input / Output) Array of length 3*n+1.
            Distributed array to be used in the basic Arnoldi iteration
            for reverse communication.  The user should not use workd as
            temporary workspace during the iteration.
    workl   (Output) Array of length lworkl+1. Private (replicated) array
            on each PE or array allocated on the front end.
    lworkl  (Input) lworkl must be at least 3*ncv*ncv+5*ncv.
    RWORK   (Workspace) Array of length ncv. Private (replicated) array on
            each PE or array allocated on  the front end.
    info    (Input / Output) On input, if info = 0, a randomly initial
            residual vector is used, otherwise resid contains the initial
            residual vector, possibly from a previous run.
            On output, info works as a error flag:
            =  0   : Normal exit.
            =  1   : Maximum number of iterations taken. All possible
                     eigenvalues of OP has been found. iparam[5]
                     returns the number of wanted converged Ritz values.
            =  3   : No shifts could be applied during a cycle of the
                     Implicitly restarted Arnoldi iteration. One
                     possibility is to increase the size of ncv relative
                     to nev. See remark 4 below.
            = -1   : n must be positive.
            = -2   : nev must be positive.
            = -3   : ncv must satisfy nev+1 <= ncv <= n.
            = -4   : The maximum number of Arnoldi update iterations
                     allowed must be greater than zero.
            = -5   : which must be one of 'LM','SM','LR','SR','LI','SI'.
            = -6   : bmat must be one of 'I' or 'G'.
            = -7   : Length of private work array is not sufficient.
            = -8   : Error return from LAPACK eigenvalue calculation.
            = -9   : Starting vector is zero.
            = -10  : iparam[7] must be 1, 2 or 3.
            = -11  : iparam[7] = 1 and bmat = 'G' are incompatible.
            = -12  : iparam[1] must be equal to 0 or 1.
            = -13  : nev and which = 'BE' are incompatible.
            = -9999: Could not build an Arnoldi factorization. iparam[5]
                     returns the size of the current Arnoldi factorization.
                     The user is advised to check that enough workspace
                     and array storage has been allocated.

  Remarks:
    1. The computed Ritz values are approximate eigenvalues of OP. The
       selection of "which" should be made with this in mind when using
       Mode = 3.  When operating in Mode = 3 setting which = 'LM' will
       compute the nev eigenvalues of the original problem that are
       closest to the shift sigma . After convergence, approximate
       eigenvalues of the original problem may be obtained with the
       ARPACK subroutine ceupp.
    2. If a basis for the invariant subspace corresponding to the converged
       Ritz values is needed, the user must call ceupp immediately following
       completion of caupp. This is new starting with release 2 of ARPACK.
    3. If M can be factored into a Cholesky factorization M = LL'
       then Mode = 2 should not be selected.  Instead one should use
       Mode = 1 with  OP = inv(L)*A*inv(L').  Appropriate triangular
       linear systems should be solved with L and L' rather
       than computing inverses.  After convergence, an approximate
       eigenvector z of the original problem is recovered by solving
       L'z = x  where x is a Ritz vector of OP.
    4. At present there is no a-priori analysis to guide the selection
       of ncv relative to nev.  The only formal requrement is that ncv
       >= nev + 1. However, it is recommended that ncv >= 2*nev. If many
       problems of the same type are to be solved, one should experiment
       with increasing ncv while keeping nev fixed for a given test
       problem. This will usually decrease the required number of OP*x
       operations but it also increases the work and storage required to
       maintain the orthogonal basis vectors.  The optimal "cross-over"
       with respect to CPU time is problem dependent and must be
       determined empirically.
    5. When iparam[1] = 0, and ido = 3, the user needs to provide the
       NP = iparam[8] complex shifts in locations
       workl[ipntr[14]], workl[ipntr[14]+1], ... , workl[ipntr[14]+NP].
       Eigenvalues of the current upper Hessenberg matrix are located in
       workl[ipntr[6]] through workl[ipntr[6]+ncv-1]. They are ordered
       according to the order defined by "which". The associated Ritz
       estimates are located in workl[ipntr[8]], workl[ipntr[8]+1], ...,
       workl[ipntr[8]+ncv-1].

  References:
    1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
       a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
       pp 357-385.
    2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
       Restarted Arnoldi Iteration", Rice University Technical Report
       TR95-13, Department of Computational and Applied Mathematics.
    3. B.N. Parlett & Y. Saad, "_Complex_ Shift and Invert Strategies for
       Double precision Matrices", Linear Algebra and its Applications,
       vol 88/89, pp 575-595, (1987).
*/

{

  F77NAME(znaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[1], &ldv, &iparam[1], &ipntr[1], &workd[1],
                  &workl[1], &lworkl, &rwork[1], &info);

} // caupp (arcomplex<double>).

inline void caupp(int& ido, char bmat, int n, char* which, int nev,
                  float& tol, arcomplex<float> resid[], int ncv,
                  arcomplex<float> V[], int ldv, int iparam[], int ipntr[],
                  arcomplex<float> workd[], arcomplex<float> workl[],
                  int lworkl, float rwork[], int& info)

/*
  c++ version of ARPACK routine cnaupd. The only difference between
  cnaupd and znaupd is that in the former function all vectors have
  single precision elements and in the latter all vectors have double
  precision elements.
*/

{

  F77NAME(cnaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[1], &ldv, &iparam[1], &ipntr[1], &workd[1],
                  &workl[1], &lworkl, &rwork[1], &info);

} // caupp (arcomplex<float>).

#endif // CAUPP_H



