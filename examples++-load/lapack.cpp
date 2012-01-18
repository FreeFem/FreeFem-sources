//ff-c++-LIBRARY-dep:   lapack
//ff-c++-LIBRARY-dep:   blas
#include "ff++.hpp"
#include "RNM.hpp"
#include "AFunction_ext.hpp" // Extension of "AFunction.hpp" to deal with more than 3 parameters function

using namespace std;

#ifdef __LP64__
  typedef int intblas;
  typedef int integer;
#else
  typedef long intblas;
  typedef long integer;
#endif

typedef integer  logical;
typedef float   LAPACK_real;
typedef double   doublereal;
typedef logical  (* L_fp)();
typedef integer      ftnlen;

typedef complex<float> LAPACK_complex;
typedef complex<double> doublecomplex;
typedef void VOID; 
#define complex LAPACK_complex 
#define real LAPACK_real 

#include "clapack.h"
#undef real
#undef complex 
class Init { public:
  Init();
};

bool lapack_inv(KNM<double>* A)
{
  intblas n=A->N();
  intblas m=A->M();
  double *a=&(*A)(0,0);
  intblas info;
  intblas lda=n;
  KN<intblas> ipiv(n);
  intblas  lw=10*n;
  KN<double> w(lw);
  ffassert(n==m);
  dgetrf_(&n,&n,a,&lda,ipiv,&info);
  if(info) return info;
  dgetri_(&n,a,&lda,ipiv,w,&lw,&info);
  return info;
}

// (computation of the eigenvalues and right eigenvectors of a real nonsymmetric matrix)
long lapack_dgeev(KNM<double> *const &A,KN<Complex> *const &vp,KNM<Complex> *const &vectp)
{
  /*
    SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
   *  JOBVL   (input) CHARACTER*1
   *          = 'N': left eigenvectors of A are not computed;
   *          = 'V': left eigenvectors of A are computed.
   *
   *  JOBVR   (input) CHARACTER*1
   *          = 'N': right eigenvectors of A are not computed;
   *          = 'V': right eigenvectors of A are computed.
   *
   *  N       (input) INTEGER
   *          The order of the matrix A. N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the N-by-N matrix A.
   *          On exit, A has been overwritten.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  WR      (output) DOUBLE PRECISION array, dimension (N)
   *  WI      (output) DOUBLE PRECISION array, dimension (N)
   *          WR and WI contain the real and imaginary parts,
   *          respectively, of the computed eigenvalues.  Complex
   *          conjugate pairs of eigenvalues appear consecutively
   *          with the eigenvalue having the positive imaginary part
   *          first.
   *
   *  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
   *          If JOBVL = 'V', the left eigenvectors u(j) are stored one
   *          after another in the columns of VL, in the same order
   *          as their eigenvalues.
   *          If JOBVL = 'N', VL is not referenced.
   *          If the j-th eigenvalue is real, then u(j) = VL(:,j),
   *          the j-th column of VL.
   *          If the j-th and (j+1)-st eigenvalues form a complex
   *          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
   *          u(j+1) = VL(:,j) - i*VL(:,j+1).
   *
   *  LDVL    (input) INTEGER
   *          The leading dimension of the array VL.  LDVL >= 1; if
   *          JOBVL = 'V', LDVL >= N.
   *
   *  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
   *          If JOBVR = 'V', the right eigenvectors v(j) are stored one
   *          after another in the columns of VR, in the same order
   *          as their eigenvalues.
   *          If JOBVR = 'N', VR is not referenced.
   *          If the j-th eigenvalue is real, then v(j) = VR(:,j),
   *          the j-th column of VR.
   *          If the j-th and (j+1)-st eigenvalues form a complex
   *          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
   *          v(j+1) = VR(:,j) - i*VR(:,j+1).
   *
   *  LDVR    (input) INTEGER
   *          The leading dimension of the array VR.  LDVR >= 1; if
   *          JOBVR = 'V', LDVR >= N.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.  LWORK >= max(1,3*N), and
   *          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
   *          performance, LWORK must generally be larger.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  if INFO = i, the QR algorithm failed to compute all the
   *                eigenvalues, and no eigenvectors have been computed;
   *                elements i+1:N of WR and WI contain eigenvalues which
   *                have converged.
   */
  intblas n=A->N();
  ffassert(A->M()==n);
  ffassert(vectp->M()>=n);
  ffassert(vectp->N()>=n);
  ffassert(vp->N()>=n);
  KN<double> wr(n),wi(n),w(1);
  KNM<double> mat(*A),vr(n,n),vl(n,n);
  intblas info,lw=-1;
  char JOBVL='N',JOBVR='V';
  dgeev_(&JOBVL,&JOBVR,&n,mat,&n,wr,wi,vl,&n,vr,&n,w,&lw,&info);
  lw=w[0];
  w.resize(lw);
  cout << mat << endl;
  dgeev_(&JOBVL,&JOBVR,&n,mat,&n,wr,wi,vl,&n,vr,&n,w,&lw,&info);
  cout << wr << endl;
  cout << wi << endl;
  if (info<0)
     {
     cout << "   dgeev: the " << info << "-th argument had an illegal value." << endl;
     (*vp)=Complex();
     (*vectp)=Complex();
     }
  else if (info>0)
     {
     cout << "   dgeev: the QR algorithm failed to compute all the eigenvalues, and no eigenvectors have been computed." << endl;
     (*vp)=Complex();
     (*vectp)=Complex();
     }
  else if (info==0)
     {
     for (int i=0;i<n;++i)
        {
        (*vp)[i]=Complex(wr[i],wi[i]);
        if (verbosity>2)
           cout << "   dgeev: vp "<< i << " : "  << (*vp)[i] << endl;
        if (wi[i]==0)
           for (int j=0;j<n;++j)
              (*vectp)(j,i)=vr(j,i);
        else if (wi[i]>0)
	   {
	   for (int j=0;j<n;++j)
              (*vectp)(j,i)=Complex(vr(j,i),vr(j,i+1));
           }
	else if (wi[i]<0)
	   {
           for (int j=0;j<n;++j)
              (*vectp)(j,i)=Complex(vr(j,i-1),-vr(j,i));	      
           }
        if (verbosity>5)
           cout << "   dgeev:   " << (*vectp)(':',i) <<endl;
        }
     }
  return info;
}

// (computation of the eigenvalues and right eigenvectors of a complex nonsymmetric matrix)
long lapack_zgeev(KNM<Complex> *const &A,KN<Complex> *const &vp,KNM<Complex> *const &vectp)
{
  intblas nvp =0,zero=0;
  intblas n= A->N();
  ffassert(A->M()==n);
  ffassert(vectp->M()>=n);
  ffassert(vectp->N()>=n);
  ffassert(vp->N()>=n);
  KN<Complex> w(n),vr(n*n),vl(n*n);
  KNM<Complex> mat(*A);
  intblas info,lw=n*(n+1)*10;
  KN<Complex> wk(lw);
  KN<double> rwk(2*n);

  char N='N',V='V';
  // lw=1;// to get opt size value 
  zgeev_(&N,&V,&n, mat,&n, w, vl,&n, vr,&n,wk,&lw,rwk,&info);
  //  cout << lw << " " << wk[0] << " " << info <<   endl;
  /* lw=wk[0].real();
  w.resize(lw);
  zgeev_(&N,&V,&n, mat,&n, w, vl,&n, vr,&n,wk,&lw,rwk,&info);
  */
  if(info)
    cout << " info =  " << info << endl;
  if(!info)
    {
      int k=0;
      for(int i=0;i<n;++i)
        {
          (*vp)[i]=w[i];
          if(verbosity>2)
            cout << "   zgeev: vp "<< i << " : "  << (*vp)[i] << endl;
	  for(int j=0;j<n;++j)
              (*vectp)(j,i)=vr[k++];
          if(verbosity>5)
            cout << "   zgeev :   " << (*vectp)(':',i) <<endl;
        }
    }
  else
    {
      nvp=0;
      (*vp)=Complex();
      (*vectp)=Complex();
    }
  return nvp;
}

// VL, 10/02/2010
long lapack_dggev(KNM<double> *const &A,KNM<double> *const &B,KN<Complex> *const &vpa,KN<double> *const &vpb,KNM<Complex> *const &vectp)
{
    intblas nvp =0,zero=0;
    intblas n=A->N();
    ffassert(A->M()==n);
    ffassert(B->M()==n);
    ffassert(B->N()==n);
    ffassert(vectp->M()>=n);
    ffassert(vectp->N()>=n);
    ffassert(vpa->N()>=n);
    ffassert(vpb->N()>=n);
    
    KN<double> war(n),wai(n),wb(n),vr(n*n),vl(n*n);
    KNM<double> matA(*A);
    KNM<double> matB(*B);
    intblas info,lw=-1;  
    KN<double> w(1);
    //char N='N',V='V'; VL: do not compute eigenvectors (if yes, switch with following line)
    char VL='N',VR='N';
    
    dggev_(&VL,&VR,&n,matA,&n,matB,&n,war,wai,wb,vl,&n,vr,&n,w,&lw,&info);
    lw=w[0];
    // cout << lw << endl;
    w.resize(lw);
    dggev_(&VL,&VR,&n,matA,&n,matB,&n,war,wai,wb,vl,&n,vr,&n,w,&lw,&info);
    if(info)
	cout << " info =  " << info << endl;
    if(!info)
      {
	int k=0;
	for(int i=0;i<n;++i)
	  {
	    (*vpa)[i]=Complex(war[i],wai[i]);
	    (*vpb)[i]=wb[i];
	    if(verbosity>2)
		cout << "   dggev: vp "<< i << " : "  << (*vpa)[i] << " ; " << (*vpb)[i] << endl;
	    if( wai[i] == 0)
		for(int j=0;j<n;++j)
		    (*vectp)(j,i)=vr[k++];
	    else if (  wai[i] >  0)
	      {
		int ki= k+n;
		for(int j=0;j<n;++j)
		    (*vectp)(j,i)=Complex(vr[k++],vr[ki++]);	      
	      }
	    else 
	      {
		int kr= k-n;
		for(int j=0;j<n;++j)
		    (*vectp)(j,i)=Complex(vr[kr++],-vr[k++]);	      
	      }
	    if(verbosity>5)
		cout << "   dggev :   " << (*vectp)(':',i) <<endl;
	  }
      }
    else
      {
	nvp=0;
	(*vpa)=Complex();
	(*vectp)=Complex();
      }
    return nvp;
}

// GL, 05/10/2011 (computation of all the eigenvalues and the eigenvectors of a real generalized symmetric-definite eigenproblem, of the form A*x=(lambda)*B*x)
long lapack_dsygvd(KNM<double> *const &A,KNM<double> *const &B,KN<double> *const &vp,KNM<double> *const &vectp)
{
  /*
    SUBROUTINE DSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )
   *  ITYPE   (input) INTEGER
   *          Specifies the problem type to be solved:
   *          = 1:  A*x = (lambda)*B*x
   *          = 2:  A*B*x = (lambda)*x
   *          = 3:  B*A*x = (lambda)*x
   *
   *  JOBZ    (input) CHARACTER*1
   *          = 'N':  Compute eigenvalues only;
   *          = 'V':  Compute eigenvalues and eigenvectors.
   *
   *  UPLO    (input) CHARACTER*1
   *          = 'U':  Upper triangles of A and B are stored;
   *          = 'L':  Lower triangles of A and B are stored.
   *
   *  N       (input) INTEGER
   *          The order of the matrices A and B.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the symmetric matrix A.  If UPLO = 'U', the
   *          leading N-by-N upper triangular part of A contains the
   *          upper triangular part of the matrix A.  If UPLO = 'L',
   *          the leading N-by-N lower triangular part of A contains
   *          the lower triangular part of the matrix A.
   *
   *          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
   *          matrix Z of eigenvectors.  The eigenvectors are normalized
   *          as follows:
   *          if ITYPE = 1 or 2, Z**T*B*Z = I;
   *          if ITYPE = 3, Z**T*inv(B)*Z = I.
   *          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
   *          or the lower triangle (if UPLO='L') of A, including the
   *          diagonal, is destroyed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
   *          On entry, the symmetric matrix B.  If UPLO = 'U', the
   *          leading N-by-N upper triangular part of B contains the
   *          upper triangular part of the matrix B.  If UPLO = 'L',
   *          the leading N-by-N lower triangular part of B contains
   *          the lower triangular part of the matrix B.
   *
   *          On exit, if INFO <= N, the part of B containing the matrix is
   *          overwritten by the triangular factor U or L from the Cholesky
   *          factorization B = U**T*U or B = L*L**T.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  W       (output) DOUBLE PRECISION array, dimension (N)
   *          If INFO = 0, the eigenvalues in ascending order.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK.
   *          If N <= 1,               LWORK >= 1.
   *          If JOBZ = 'N' and N > 1, LWORK >= 2*N+1.
   *          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal sizes of the WORK and IWORK
   *          arrays, returns these values as the first entries of the WORK
   *          and IWORK arrays, and no error message related to LWORK or
   *          LIWORK is issued by XERBLA.
   *
   *  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
   *          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
   *
   *  LIWORK  (input) INTEGER
   *          The dimension of the array IWORK.
   *          If N <= 1,                LIWORK >= 1.
   *          If JOBZ  = 'N' and N > 1, LIWORK >= 1.
   *          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
   *
   *          If LIWORK = -1, then a workspace query is assumed; the
   *          routine only calculates the optimal sizes of the WORK and
   *          IWORK arrays, returns these values as the first entries of
   *          the WORK and IWORK arrays, and no error message related to
   *          LWORK or LIWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  DPOTRF or DSYEVD returned an error code:
   *             <= N:  if INFO = i and JOBZ = 'N', then the algorithm
   *                    failed to converge; i off-diagonal elements of an
   *                    intermediate tridiagonal form did not converge to
   *                    zero;
   *                    if INFO = i and JOBZ = 'V', then the algorithm
   *                    failed to compute an eigenvalue while working on
   *                    the submatrix lying in rows and columns INFO/(N+1)
   *                    through mod(INFO,N+1);
   *             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
   *                    minor of order i of B is not positive definite.
   *                    The factorization of B could not be completed and
   *                    no eigenvalues or eigenvectors were computed.
   */
  intblas n=A->N();
  ffassert(A->M()==n);
  ffassert(B->M()==n);
  ffassert(B->N()==n);
  ffassert(vp->N()>=n);
  ffassert(vectp->M()>=n);
  ffassert(vectp->N()>=n);
  KN<double> war(n),wai(n),wb(n),vr(n*n),vl(n*n);
  KNM<double> matA(*A),matB(*B);
  intblas itype=1,info,lw=-1;
  KN<double> w(1);
  KN<intblas> iw(1);
  char JOBZ='V',UPLO='U';
  
  dsygvd_(&itype,&JOBZ,&UPLO,&n,matA,&n,matB,&n,*vp,w,&lw,iw,&lw,&info);
  lw=w[0];
  w.resize(lw);
  iw.resize(lw);
  dsygvd_(&itype,&JOBZ,&UPLO,&n,matA,&n,matB,&n,*vp,w,&lw,iw,&lw,&info);
  if (info<0)
     {
     cout << "   dsygvd: the " << info << "-th argument had an illegal value." << endl;
     }
  else if (info>0)
     {
     cout << "   dsygvd: DPOTRF or DSYEVD returned an error code." << endl;
     }
  else if (info==0)
     {
     for (int i=0;i<n;++i)
        {
        for (int i=0;i<n;++i)
           {
           for (int j=0;j<n;++j)
              (*vectp)(j,i)=matA(j,i);	      
           }
        }  
     }
  return info;
}

// GL,27/09/2011 (singular value decomposition of a rectangular real matrix)
long lapack_dgesdd(KNM<double> *const &A,KNM<double> *const &U,KN<double> *const &S,KNM<double> *const &V)
{
  /*
    SUBROUTINE DGESDD( JOBZ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, INFO )
   *  JOBZ    (input) CHARACTER*1
   *          Specifies options for computing all or part of the matrix U:
   *          = 'A':  all M columns of U and all N rows of V**T are
   *                  returned in the arrays U and VT;
   *          = 'S':  the first min(M,N) columns of U and the first
   *                  min(M,N) rows of V**T are returned in the arrays U
   *                  and VT;
   *          = 'O':  If M >= N, the first N columns of U are overwritten
   *                  on the array A and all rows of V**T are returned in
   *                  the array VT;
   *                  otherwise, all columns of U are returned in the
   *                  array U and the first M rows of V**T are overwritten
   *                  in the array A;
   *          = 'N':  no columns of U or rows of V**T are computed.
   *
   *  M       (input) INTEGER
   *          The number of rows of the input matrix A.  M >= 0.
   *
   *  N       (input) INTEGER
   *          The number of columns of the input matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the M-by-N matrix A.
   *          On exit,
   *          if JOBZ = 'O',  A is overwritten with the first N columns
   *                          of U (the left singular vectors, stored
   *                          columnwise) if M >= N;
   *                          A is overwritten with the first M rows
   *                          of V**T (the right singular vectors, stored
   *                          rowwise) otherwise.
   *          if JOBZ .ne. 'O', the contents of A are destroyed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,M).
   *
   *  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
   *          The singular values of A, sorted so that S(i) >= S(i+1).
   *
   *  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
   *          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
   *          UCOL = min(M,N) if JOBZ = 'S'.
   *          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
   *          orthogonal matrix U;
   *          if JOBZ = 'S', U contains the first min(M,N) columns of U
   *          (the left singular vectors, stored columnwise);
   *          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
   *
   *  LDU     (input) INTEGER
   *          The leading dimension of the array U.  LDU >= 1; if
   *          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
   *
   *  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
   *          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
   *          N-by-N orthogonal matrix V**T;
   *          if JOBZ = 'S', VT contains the first min(M,N) rows of
   *          V**T (the right singular vectors, stored rowwise);
   *          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
   *
   *  LDVT    (input) INTEGER
   *          The leading dimension of the array VT.  LDVT >= 1; if
   *          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
   *          if JOBZ = 'S', LDVT >= min(M,N).
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
   *
   *  LWORK   (input) INTEGER
   *          The dimension of the array WORK. LWORK >= 1.
   *          If JOBZ = 'N',
   *            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)).
   *          If JOBZ = 'O',
   *            LWORK >= 3*min(M,N) + 
   *                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
   *          If JOBZ = 'S' or 'A'
   *            LWORK >= 3*min(M,N) +
   *                     max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N)).
   *          For good performance, LWORK should generally be larger.
   *          If LWORK = -1 but other input arguments are legal, WORK(1)
   *          returns the optimal LWORK.
   *
   *  IWORK   (workspace) INTEGER array, dimension (8*min(M,N))
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit.
   *          < 0:  if INFO = -i, the i-th argument had an illegal value.
   *          > 0:  DBDSDC did not converge, updating process failed.
   */
  intblas n=A->N();
  intblas m=A->M();
  U->resize(n,n);
  S->resize(min(n,m));
  V->resize(m,m);
  KNM<double> VT(m,m);
  KN<intblas> iw(8*min(n,m));
  intblas info,lw=-1;
  KN<double> w(1);
  char JOBZ='A';
  dgesdd_(&JOBZ,&n,&m,*A,&n,*S,*U,&n,VT,&m,w,&lw,iw,&info);
  lw=w[0];
  w.resize(lw);
  dgesdd_(&JOBZ,&n,&m,*A,&n,*S,*U,&n,VT,&m,w,&lw,iw,&info);
  if (info<0)
     {
     cout << "   dgesdd: the " << info << "-th argument had an illegal value." << endl;
     }
  else if (info>0)
     {
     cout << "   dgesdd: DBDSDC did not converge, updating process failed." << endl;
     }
  else if (info==0)
     {
     for (int i=0;i<m;++i)
	for (int j=0;j<m;++j)
	   (*V)(i,j)=VT(j,i);
     }
  return info;
}

// GL,28/09/2011 (computation of the eigenvalues and eigenvectors of a real symmetric matrix)
long lapack_dsyev(KNM<double> *const &A,KN<double> *const &vp,KNM<double> *const &vectp)
{
  /*
    SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
   *  JOBZ    (input) CHARACTER*1
   *          = 'N':  Compute eigenvalues only;
   *          = 'V':  Compute eigenvalues and eigenvectors.
   *
   *  UPLO    (input) CHARACTER*1
   *          = 'U':  Upper triangle of A is stored;
   *          = 'L':  Lower triangle of A is stored.
   *
   *  N       (input) INTEGER
   *          The order of the matrix A.  N >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
   *          On entry, the symmetric matrix A.  If UPLO = 'U', the
   *          leading N-by-N upper triangular part of A contains the
   *          upper triangular part of the matrix A.  If UPLO = 'L',
   *          the leading N-by-N lower triangular part of A contains
   *          the lower triangular part of the matrix A.
   *          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
   *          orthonormal eigenvectors of the matrix A.
   *          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
   *          or the upper triangle (if UPLO='U') of A, including the
   *          diagonal, is destroyed.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  W       (output) DOUBLE PRECISION array, dimension (N)
   *          If INFO = 0, the eigenvalues in ascending order.
   *
   *  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
   *          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
   *
   *  LWORK   (input) INTEGER
   *          The length of the array WORK.  LWORK >= max(1,3*N-1).
   *          For optimal efficiency, LWORK >= (NB+2)*N,
   *          where NB is the blocksize for DSYTRD returned by ILAENV.
   *
   *          If LWORK = -1, then a workspace query is assumed; the routine
   *          only calculates the optimal size of the WORK array, returns
   *          this value as the first entry of the WORK array, and no error
   *          message related to LWORK is issued by XERBLA.
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, the algorithm failed to converge; i
   *                off-diagonal elements of an intermediate tridiagonal
   *                form did not converge to zero.
   */
  intblas n=A->N();
  ffassert(A->M()==n);
  ffassert(vectp->N()==n);
  ffassert(vectp->M()==n);
  ffassert(vp->N()==n);
  KNM<double> mat(*A);
  intblas info,lw=-1;  
  KN<double> w(1);
  char JOBZ='V',UPLO='U';
  dsyev_(&JOBZ,&UPLO,&n,mat,&n,*vp,w,&lw,&info);
  lw=w[0];
  w.resize(lw);
  dsyev_(&JOBZ,&UPLO,&n,mat,&n,*vp,w,&lw,&info);
  if (info<0)
     {
     cout << "   dsyev: the " << info << "-th argument had an illegal value." << endl;
     }
  else if (info>0)
     {
     cout << "   dsyev: the algorithm failed to converge." << endl;
     }
  else if (info==0)
     {
     *vectp=mat;
     }
  return info;
}

template<class T>
class Inverse{ public:
  T  t;
  Inverse( T  v)
   : t(v) {}
  template<class TT> Inverse( TT  v) : t(v) {}  
  template<class TT> Inverse( TT * v) : t(*v) {}  
  operator const T & () const {return t;}
};

template<class T>
class Mult{ public:
    T  a;bool ta;
  T  b;bool tb;
  Mult( T  aa,T bb)
    : a(aa),b(bb),ta(0),tb(0) {}
    // Transpose<
  Mult( Transpose<T>  aa,T bb)   
    : a(aa),b(bb),ta(1),tb(0) {}
  Mult( Transpose<T>  aa,Transpose<T> bb)   
    : a(aa),b(bb),ta(1),tb(1) {}
    Mult( T  aa,Transpose<T> bb)   
    : a(aa),b(bb),ta(1),tb(1) {}
    
};

template<class K>
class OneBinaryOperatorRNM_inv : public OneOperator { public:  
    OneBinaryOperatorRNM_inv() 
      : OneOperator( atype< Inverse< KNM<K>* > >(),atype<KNM<K> *>(),atype<long>()) {}
  E_F0 * code(const basicAC_F0 & args) const 
  { Expression p=args[1];
    if ( ! p->EvaluableWithOutStack() ) 
      { 
	bool bb=p->EvaluableWithOutStack();
	cout << bb << " " <<  * p <<  endl;
	CompileError(" A^p, The p must be a constant == -1, sorry");}
    long pv = GetAny<long>((*p)(0));
    if (pv !=-1)   
      { char buf[100];
	sprintf(buf," A^%ld, The pow must be  == -1, sorry",pv);
	CompileError(buf);}     
    return  new E_F_F0<Inverse< KNM<K>* > ,KNM<K> *>(Build<Inverse< KNM<K>* > ,KNM<K> *>,t[0]->CastTo(args[0])); 
  }
};


/* 
class Init { public:
  Init();
};
*/

KNM<R>* Solve(KNM<R>* a,Inverse<KNM<R >*> b) 
{
  /*
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
   *  N       (input) INTEGER
   *          The number of linear equations, i.e., the order of the
   *          matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the N-by-N coefficient matrix A.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U; the unit diagonal elements of L are not stored.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (output) INTEGER array, dimension (N)
   *          The pivot indices that define the permutation matrix P;
   *          row i of the matrix was interchanged with row IPIV(i).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the N-by-NRHS matrix of right hand side matrix B.
   *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
   *                has been completed, but the factor U is exactly
   *                singular, so the solution could not be computed.
   *
   */
  typedef double R;
  integer info;
  KNM<R> B(*b);
  integer  n= B.N();
  KN<integer> p(n);
  ffassert(B.M()==n);
  a->resize(n,n);
  *a=0.;
  for(int i=0;i<n;++i)
    (*a)(i,i)=(R) 1;;

  dgesv_(&n,&n,B,&n,p,*a,&n,&info);
  if(info) cerr << " error:  dgesv_ "<< info << endl;
  return a;
}

// Template interface 
inline int gemm(char *transa, char *transb, integer *m, integer *
	   n, integer *k, double *alpha, double *a, integer *lda, 
	   double *b, integer *ldb, double *beta, double *c, integer 
	   *ldc) {
   return  dgemm_(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}
inline int gemm(char *transa, char *transb, integer *m, integer *
		 n, integer *k, Complex *alpha, Complex *a, integer *lda, 
		 Complex *b, integer *ldb, Complex *beta, Complex *c, integer 
		 *ldc) {
    return  zgemm_(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc);
}


template<class R,bool init, int ibeta> 
KNM<R>* mult(KNM<R >* a,const KNM_<R> & A,const KNM_<R> & B) 
{ // C=A*B 
   
    R alpha=1.,beta=R(ibeta);
    char tA, tB;
    if(init) a->init();
    intblas N= A.N();
    intblas M=B.M();
    intblas K=A.M();
    KNM<R> & C= *a;
    C.resize(N,M); 
    ffassert(K==B.N());
    R *A00=&A(0,0), *A10= &A(1,0), *A01= &A(0,1); 
    R *B00=&B(0,0), *B10= &B(1,0), *B01= &B(0,1); 
    R *C00=&C(0,0), *C10= &C(1,0), *C01= &C(0,1); 
    intblas lsa=A10-A00 ,lsb=B10-B00,lsc=C10-C00;
    intblas lda=A01-A00 ,ldb=B01-B00,ldc=C01-C00;
    if(verbosity>10) {
	cout << lsa << " " << lsb << " "<< lsc << " init " << init <<  endl;
	cout << lda << " " << ldb << " "<< ldc << endl;	
    }
    tA=lda==1?'T':'N';
    tB=ldb==1?'T':'N';
    
    if(lda==1) lda=lsa;
    if(ldb==1) ldb=lsb;
    if(beta==0.)
     C=R(); 
#ifdef XXXXXXXXXXXXXX
    
    for(int i=0;i<N;++i)
	 for(int j=0;j<M;++j)
	      for(int k=0;k<K;++k)
		  C(i,j) += A(i,k)*B(k,j)  ;
#else    
    gemm(&tB,&tA,&N,&M,&K,&alpha,A00,&lda,B00,&ldb,&beta,C00,&ldc);
#endif
    return a;
    /*
     The Fortran interface for these procedures are:
     SUBROUTINE xGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
     where TRANSA and TRANSB determines if the matrices A and B are to be transposed.
     M is the number of rows in matrix A and C. N is the number of columns in matrix B and C. 
     K is the number of columns in matrix A and rows in matrix B. 
     LDA, LDB and LDC specifies the size of the first dimension of the matrices, as laid out in memory;
     meaning the memory distance between the start of each row/column, depending on the memory structure (Dongarra et al. 1990).
     */
}
template<class R,bool init, int ibeta> 
KNM<R>* mult(KNM<R >* a,Mult<KNM<R >*> bc) 
{
    if( (bc.ta == 0) && (bc.tb == 0))
     return  mult<R,init,ibeta>(a,*bc.a,*bc.b) ;
    else if((bc.ta == 1 )&& (bc.tb == 0))
     return  mult<R,init,ibeta>(a,bc.a->t(),*bc.b) ;
    else if((bc.ta == 0) && (bc.tb == 1))
	return  mult<R,init,ibeta>(a,*bc.a,bc.b->t()) ;
    else if((bc.ta == 1) && (bc.tb == 1))
	return  mult<R,init,ibeta>(a,bc.a->t(),bc.b->t()) ;

}

KNM<Complex>* SolveC(KNM<Complex>* a,Inverse<KNM<Complex >*> b) 
{
  /*
    SUBROUTINE DGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
   *  N       (input) INTEGER
   *          The number of linear equations, i.e., the order of the
   *          matrix A.  N >= 0.
   *
   *  NRHS    (input) INTEGER
   *          The number of right hand sides, i.e., the number of columns
   *          of the matrix B.  NRHS >= 0.
   *
   *  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
   *          On entry, the N-by-N coefficient matrix A.
   *          On exit, the factors L and U from the factorization
   *          A = P*L*U; the unit diagonal elements of L are not stored.
   *
   *  LDA     (input) INTEGER
   *          The leading dimension of the array A.  LDA >= max(1,N).
   *
   *  IPIV    (output) INTEGER array, dimension (N)
   *          The pivot indices that define the permutation matrix P;
   *          row i of the matrix was interchanged with row IPIV(i).
   *
   *  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
   *          On entry, the N-by-NRHS matrix of right hand side matrix B.
   *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
   *
   *  LDB     (input) INTEGER
   *          The leading dimension of the array B.  LDB >= max(1,N).
   *
   *  INFO    (output) INTEGER
   *          = 0:  successful exit
   *          < 0:  if INFO = -i, the i-th argument had an illegal value
   *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
   *                has been completed, but the factor U is exactly
   *                singular, so the solution could not be computed.
   *
   */
  typedef Complex R;
  integer info;
  KNM<R> B(*b);
  integer   n= B.N();
  KN<integer> p(n);
  ffassert(B.M()==n);
  a->resize(n,n);
  *a=0.;
  for(int i=0;i<n;++i)
    (*a)(i,i)=(R) 1;;

  zgesv_(&n,&n,(R*) B,&n,p, (R*) *a,&n,&info);
  if(info) cerr << " error:  zgesv_ "<< info << endl;
  return a;
}
/*
LOADINIT(Init);  //  une variable globale qui serat construite  au chargement dynamique 

Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 
  // avec de matrice plein 
  //  B = A ^-1 * C;
  //  A = A ^-1;
  Dcl_Type< Inverse<KNM<double >* > > ();
  Dcl_Type< Inverse<KNM<Complex >* > > ();
  
  TheOperators->Add("^", new OneBinaryOperatorRNM_inv<double>());
  TheOperators->Add("^", new OneBinaryOperatorRNM_inv<Complex>());
  TheOperators->Add("=", new OneOperator2<KNM<double>*,KNM<double>*,Inverse<KNM<double >*> >( Solve) );
  TheOperators->Add("=", new OneOperator2<KNM<Complex>*,KNM<Complex>*,Inverse<KNM<Complex >*> >( SolveC) );
  

 
}

*/
LOADINIT(Init);  //  une variable globale qui serat construite  au chargement dynamique 

template<class R,class A,class B> R Build2(A a,B b) {
    return R(a,b);
}
Init::Init(){  // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++ 

  if( map_type.find(typeid(Inverse<KNM<double >* >).name() ) == map_type.end() )
    {
      if(verbosity) 
	cout << " Add lapack interface ..." ;
      Dcl_Type< Inverse<KNM<double >* > > ();
      Dcl_Type< Inverse<KNM<Complex >* > > ();
      Dcl_Type< Mult<KNM<Complex >* > > ();
      Dcl_Type< Mult<KNM<double >* > > ();
      
      TheOperators->Add("^", new OneBinaryOperatorRNM_inv<double>());
      TheOperators->Add("*", new OneOperator2< Mult< KNM<double>* >,KNM<double>*,KNM<double>*>(Build2));
      TheOperators->Add("*", new OneOperator2< Mult< KNM<Complex>* >,KNM<Complex>*,KNM<Complex>*>(Build2));
      
      TheOperators->Add("^", new OneBinaryOperatorRNM_inv<Complex>());
      TheOperators->Add("=", new OneOperator2<KNM<double>*,KNM<double>*,Inverse<KNM<double >*> >( Solve) );
      TheOperators->Add("=", new OneOperator2<KNM<Complex>*,KNM<Complex>*,Inverse<KNM<Complex >*> >( SolveC) );
      TheOperators->Add("=", new OneOperator2<KNM<double>*,KNM<double>*,Mult<KNM<double >*> >( mult<double,false,0> ) );
      TheOperators->Add("=", new OneOperator2<KNM<Complex>*,KNM<Complex>*,Mult<KNM<Complex >*> >( mult<Complex,false,0> ) );
      
      TheOperators->Add("+=", new OneOperator2<KNM<double>*,KNM<double>*,Mult<KNM<double >*> >( mult<double,false,1> ) );
      TheOperators->Add("+=", new OneOperator2<KNM<Complex>*,KNM<Complex>*,Mult<KNM<Complex >*> >( mult<Complex,false,1> ) );
      
      TheOperators->Add("-=", new OneOperator2<KNM<double>*,KNM<double>*,Mult<KNM<double >*> >( mult<double,false,-1> ) );
      TheOperators->Add("-=", new OneOperator2<KNM<Complex>*,KNM<Complex>*,Mult<KNM<Complex >*> >( mult<Complex,false,-1> ) );
      
      TheOperators->Add("<-", new OneOperator2<KNM<double>*,KNM<double>*,Mult<KNM<double >*> >( mult<double,true,0> ) );
      TheOperators->Add("<-", new OneOperator2<KNM<Complex>*,KNM<Complex>*,Mult<KNM<Complex >*> >( mult<Complex,true,0> ) );
      
      Global.Add("inv","(",new  OneOperator1<bool,KNM<double>*>(lapack_inv));  
      Global.Add("dgeev","(",new  OneOperator3_<long,KNM<double>*,KN<Complex>*,KNM<Complex>*>(lapack_dgeev));  
      Global.Add("zgeev","(",new  OneOperator3_<long,KNM<Complex>*,KN<Complex>*,KNM<Complex>*>(lapack_zgeev));  
      Global.Add("dggev","(",new  OneOperator5_<long,KNM<double>*,KNM<double>*,KN<Complex>*,KN<double>*,KNM<Complex>*>(lapack_dggev));
      Global.Add("dsygvd","(",new  OneOperator4_<long,KNM<double>*,KNM<double>*,KN<double>*,KNM<double>*>(lapack_dsygvd));
      Global.Add("dgesdd","(",new  OneOperator4_<long,KNM<double>*,KNM<double>*,KN<double>*,KNM<double>*>(lapack_dgesdd));
      Global.Add("dsyev","(",new  OneOperator3_<long,KNM<double>*,KN<double>*,KNM<double>*>(lapack_dsyev));
    }
  else
    if(verbosity)
      cout << "( load: lapack <=> fflapack , skeep ) ";
}

