/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE saupp.h.
   Interface to ARPACK subroutines dsaupd and ssaupd.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SAUPP_H
#define SAUPP_H

#include "arch.h"
#include "arpackf.h"

inline void saupp(int& ido, char bmat, int n, char* which, int nev,
                  double& tol, double resid[], int ncv, double V[],
                  int ldv, int iparam[], int ipntr[], double workd[],
                  double workl[], int lworkl, int& info)

/*
  c++ version of ARPACK routine dsaupd that implements a variant of
  the Lanczos method.  This method has been designed to compute
  approximations to a few eigenpairs of a linear operator OP that is
  real and symmetric with respect to a real positive semi-definite
  symmetric matrix B, i.e.

                            B*OP = (OP')*B.

  where A' denotes transpose of A. In the standard eigenproblem B is
  the identity matrix. Another way to express this condition is

             < x,OPy > = < OPx,y >  where < z,w > = z'Bw.

  The computed approximate eigenvalues are called Ritz values and
  the corresponding approximate eigenvectors are called Ritz vectors.

  saupp is usually called iteratively to solve one of the
  following problems:

  Mode 1:  A*x = lambda*x, A symmetric
           ===> OP = A  and  B = I.

  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
           ===> OP = inv[M]*A  and  B = M.
           ===> (If M can be factored see remark 3 below)

  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
           ===> OP = (inv[K - sigma*M])*M  and  B = M.
           ===> Shift-and-Invert mode

  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite,
           KG symmetric indefinite
           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
           ===> Buckling mode

  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
           ===> Cayley transformed mode

  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v should be
        accomplished either by a direct method using a sparse matrix
        factorization and solving

                     [A - sigma*M]*w = v  or M*w = v,

        or through an iterative method for solving these systems.  If an
        iterative method is used, the convergence test must be more
        stringent than the accuracy requirements for the eigenvalue
        approximations.

  Parameters:

    ido     (Input / Output) Reverse communication flag.  ido must be
            zero on the first call to saupp.  ido will be set
            internally to indicate the type of operation to be
            performed.  Control is then given back to the calling
            routine which has the responsibility to carry out the
            requested operation and call saupp with the result. The
            operand is given in workd[ipntr[1]], the result must be
            put in workd[ipntr[2]]. (If Mode = 2 see remark 5 below).
            ido =  0: first call to the reverse communication interface.
            ido = -1: compute  Y = OP * X  where
                      ipntr[1] is the pointer into workd for X,
                      ipntr[2] is the pointer into workd for Y.
                      This is for the initialization phase to force the
                      starting vector into the range of OP.
            ido =  1: compute  Y = OP * X where
                      ipntr[1] is the pointer into workd for X,
                      ipntr[2] is the pointer into workd for Y.
                      In mode 3,4 and 5, the vector B * X is already
                      available in workd[ipntr[3]].  It does not
                      need to be recomputed in forming OP * X.
            ido =  2: compute  Y = B * X  where
                      ipntr[1] is the pointer into workd for X,
                      ipntr[2] is the pointer into workd for Y.
            ido =  3: compute the iparam[8] shifts where
                      ipntr[11] is the pointer into workl for
                      placing the shifts. See remark 6 below.
            ido = 99: done.
    bmat    (Input) bmat specifies the type of the matrix B that defines
            the semi-inner product for the operator OP.
            bmat = 'I' -> standard eigenvalue problem A*x = lambda*x;
            bmat = 'G' -> generalized eigenvalue problem A*x = lambda*B*x.
    n       (Input) Dimension of the eigenproblem.
    nev     (Input) Number of eigenvalues to be computed. 0 < nev < n.
    which   (Input) Specify which of the Ritz values of OP to compute.
            'LA' - compute the nev largest (algebraic) eigenvalues.
            'SA' - compute the nev smallest (algebraic) eigenvalues.
            'LM' - compute the nev largest (in magnitude) eigenvalues.
            'SM' - compute the nev smallest (in magnitude) eigenvalues.
            'BE' - compute nev eigenvalues, half from each end of the
                   spectrum.  When NEV is odd, compute one more from the
                   high end than from the low end.
            (see remark 1 below)
    tol     (Input) Stopping criterion: the relative accuracy of the
            Ritz value is considered acceptable if BOUNDS[i] <=
            tol*abs(RITZ[i]). If tol<=0.0 is passed, the machine
            precision as computed by the LAPACK auxiliary subroutine
            _LAMCH is used.
    resid   (Input / Output) Array of length n.
            On input:
            If info==0, a random initial residual vector is used.
            If info!=0, resid contains the initial residual vector,
                        possibly from a previous run.
            On output:
            resid contains the final residual vector.
    ncv     (Input) Number of Lanczos vectors that are generated at each
            iteration. After the startup phase in which nev Lanczos
            vectors are generated, the algorithm generates ncv-nev
            Lanczos vectors at each subsequent update iteration. Most of
            the cost in generating each Lanczos vector is in the
            matrix-vector product OP*x. (See remark 4 below).
    V       (Output) Double precision array of length ncv*n+1. V contains
            the ncv Lanczos basis vectors. The first element V[0] is never
            referenced.
    ldv     (Input) Dimension of the basis vectors contianed in V. This
            parameter MUST be set to n.
    iparam  (Input / Output) Array of length 12.
            iparam[1]  = ISHIFT: method for selecting the implicit shifts.
            The shifts selected at each iteration are used to restart
            the Arnoldi iteration in an implicit fashion.
            -------------------------------------------------------------
            ISHIFT = 0: the shifts are provided by the user via
                        reverse communication.  The NCV eigenvalues of
                        the current tridiagonal matrix T are returned in
                        the part of workl array corresponding to RITZ.
                        See remark 6 below.
            ISHIFT = 1: exact shifts with respect to the reduced
                        tridiagonal matrix T.  This is equivalent to
                        restarting the iteration with a starting vector
                        that is a linear combination of Ritz vectors
                        associated with the "wanted" Ritz values.
            -------------------------------------------------------------
            iparam[2] is no longer referenced.
            iparam[3]  = MXITER
            On INPUT:  maximum number of Arnoldi update iterations allowed.
            On OUTPUT: actual number of Arnoldi update iterations taken.
            iparam[4]  = NB: blocksize to be used in the recurrence.
            The code currently works only for NB = 1.
            iparam[5]  = NCONV: number of "converged" Ritz values.
            This represents the number of Ritz values that satisfy
            the convergence criterion.
            iparam[6] is no longer referenced.
            iparam[7]  = MODE. On INPUT determines what type of
            eigenproblem is being solved. Must be 1,2,3,4,5.
            iparam[8]  = NP. When ido = 3 and the user provides shifts
            through reverse communication (iparam[1]=0), saupp returns
            NP, the number of shifts the user is to provide.
            0 < NP <=ncv-nev. See Remark 6 below.
            iparam[9]  =  total number of OP*x operations.
            iparam[10] = total number of B*x operations if bmat='G'.
            iparam[11] = total number of steps of re-orthogonalization.
    ipntr   (Output) Array of length 12. Pointer to mark the starting
            locations in the workd and workl arrays for matrices/vectors
            used by the Lanczos iteration.
            ipntr[1] : pointer to the current operand vector X in workd.
            ipntr[2] : pointer to the current result vector Y in workd.
            ipntr[3] : pointer to the vector B * X in workd when used in
                       the shift-and-invert mode.
            ipntr[4] : pointer to the next available location in workl
                       that is untouched by the program.
            ipntr[5] : pointer to the ncv by 2 tridiagonal matrix T in
                       workl.
            ipntr[6] : pointer to the ncv RITZ values array in workl.
            ipntr[7] : pointer to the Ritz estimates in array workl
                       associated with the Ritz values located in RITZ
                       in workl.
            ipntr[11]: pointer to the np shifts in workl. See Remark 6.
            Note: ipntr[8:10] is only referenced by seupp. See Remark 2.
            ipntr[8] : pointer to the ncv RITZ values of the original
                       system.
            ipntr[9] : pointer to the ncv corresponding error bounds.
            ipntr[10]: pointer to the ncv by ncv matrix of eigenvectors
                       of the tridiagonal matrix T. Only referenced by
                       seupp if RVEC = TRUE. See Remarks.
    workd   (Input / Output) Array of length 3*N+1.
            Distributed array to be used in the basic Arnoldi iteration
            for reverse communication.  The user should not use workd as
            temporary workspace during the iteration. Upon termination
            workd[1:n] contains B*resid[1:n]. If the Ritz vectors are
            desired subroutine seupp uses this output.
    workl   (Output) Array of length lworkl+1. Private (replicated) array
            on each PE or array allocated on the front end.
    lworkl  (Input) lworkl must be at least ncv*(ncv+8).
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
                     possibility is to increase the size of NCV relative
                     to nev. See remark 4 below.
            = -1   : n must be positive.
            = -2   : nev must be positive.
            = -3   : ncv must satisfy nev < ncv <= n.
            = -4   : The maximum number of Arnoldi update iterations allowed
                     must be greater than zero.
            = -5   : which must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
            = -6   : bmat must be one of 'I' or 'G'.
            = -7   : Length of private work array workl is not sufficient.
            = -8   : Error return from trid. eigenvalue calculation;
                     Informational error from LAPACK routine dsteqr.
            = -9   : Starting vector is zero.
            = -10  : iparam[7] must be 1,2,3,4,5.
            = -11  : iparam[7] = 1 and bmat = 'G' are incompatible.
            = -12  : iparam[1] must be equal to 0 or 1.
            = -13  : nev and which = 'BE' are incompatible.
            = -9999: Could not build an Arnoldi factorization. iparam[5]
                     returns the size of the current Arnoldi factorization.
                     The user is advised to check that enough workspace
                     and array storage has been allocated.

  Remarks:
    1. The converged Ritz values are always returned in ascending
       algebraic order.  The computed Ritz values are approximate
       eigenvalues of OP.  The selection of "which" should be made
       with this in mind when Mode = 3,4,5.  After convergence,
       approximate eigenvalues of the original problem may be obtained
       with the ARPACK subroutine seupp.
    2. If the Ritz vectors corresponding to the converged Ritz values are
       needed, the user must call seupp immediately following completion
       of saupp. This is new starting with version 2.1 of ARPACK.
    3. If M can be factored into a Cholesky factorization M = LL'
       then Mode = 2 should not be selected.  Instead one should use
       Mode = 1 with  OP = inv(L)*A*inv(L').  Appropriate triangular
       linear systems should be solved with L and L' rather
       than computing inverses.  After convergence, an approximate
       eigenvector z of the original problem is recovered by solving
       L'z = x  where x is a Ritz vector of OP.
    4. At present there is no a-priori analysis to guide the selection
       of ncv relative to nev.  The only formal requrement is that
       ncv > nev. However, it is recommended that ncv >= 2*nev. If many
       problems of the same type are to be solved, one should experiment
       with increasing ncv while keeping nev fixed for a given test
       problem. This will usually decrease the required number of OP*x
       operations but it also increases the work and storage required to
       maintain the orthogonal basis vectors.   The optimal "cross-over"
       with respect to CPU time is problem dependent and must be
       determined empirically.
    5. If iparam[7] = 2 then in the Reverse commuication interface the
       user must do the following. When ido = 1, Y = OP * X is to be
       computed. When iparam[7] = 2 OP = inv(B)*A. After computing A*X
       the user must overwrite X with A*X. Y is then the solution to the
       linear set of equations B*Y = A*X.
    6. When iparam[1] = 0, and ido = 3, the user needs to provide the
       NP = iparam[8] shifts in locations:
       1   workl[ipntr[11]]
       2   workl[ipntr[11]+1]
                          .
                          .
                          .
       NP  workl[ipntr[11]+NP-1].
       The eigenvalues of the current tridiagonal matrix are located in
       workl[ipntr[6]] through workl[ipntr[6]+ncv]. They are in the
       order defined by which. The associated Ritz estimates are located in
       workl[ipntr[8]], workl[ipntr[8]+1], ... , workl[ipntr[8]+ncv-1].
*/

{

  F77NAME(dsaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[1], &ldv, &iparam[1], &ipntr[1], &workd[1], &workl[1],
                  &lworkl, &info);

} // saupp (double).

inline void saupp(int& ido, char bmat, int n, char* which, int nev,
                  float& tol, float resid[], int ncv, float V[],
                  int ldv, int iparam[], int ipntr[], float workd[],
                  float workl[], int lworkl, int& info)

/*
  c++ version of ARPACK routine ssaupd. The only difference between
  ssaupd and dsaupd is that in the former function all vectors have
  single precision elements and in the latter all vectors have double
  precision elements.
*/

{

  F77NAME(ssaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[1], &ldv, &iparam[1], &ipntr[1], &workd[1], &workl[1],
                  &lworkl, &info);

} // saupp (float).

#endif // SAUPP_H

