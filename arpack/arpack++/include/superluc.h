/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE SuperLUc.h.
   Interface to SuperLU routines.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SUPERLUC_H
#define SUPERLUC_H

#include "arch.h"
#include "arlspdef.h"
#include "arlsupm.h"
#include "arlcomp.h"

// gstrf.

inline void gstrf(char* refact, SuperMatrix* A, double diag_pivot_thresh,
                  double drop_tol, int relax, int panel_size,
                  int* etree, void* work, int lwork, int* perm_r,
                  int *perm_c, SuperMatrix *L, SuperMatrix *U, int *info)
{

  if (A->Dtype == _D) {       // calling the double precision routine.
    dgstrf(refact,A,diag_pivot_thresh,drop_tol,relax,
           panel_size,etree,work,lwork,perm_r,perm_c,L,U,info);
  }
  else if (A->Dtype == _S) {  // calling the single precision routine.
    sgstrf(refact,A,(float)diag_pivot_thresh,(float)drop_tol,relax,
           panel_size,etree,work,lwork,perm_r,perm_c,L,U,info);
  }
  else if (A->Dtype == _Z) {  // calling the double precision complex routine.
#ifdef ARCOMP_H
    zgstrf(refact,A,diag_pivot_thresh,drop_tol,relax,
           panel_size,etree,work,lwork,perm_r,perm_c,L,U,info);
#endif
  }
  else {                      // calling the single precision complex routine.
#ifdef ARCOMP_H
    cgstrf(refact,A,(float)diag_pivot_thresh,(float)drop_tol,relax,
           panel_size,etree,work,lwork,perm_r,perm_c,L,U,info);
#endif
  }

} // gstrf.


inline void gstrs(char *trans, SuperMatrix *L, SuperMatrix *U,
	          int *perm_r, int *perm_c, SuperMatrix *B, int *info)
{

  if (L->Dtype == _D) {       // calling the double precision routine.
    dgstrs(trans,L,U,perm_r,perm_c,B,info);
  }
  else if (L->Dtype == _S) {  // calling the single precision routine.
    sgstrs(trans,L,U,perm_r,perm_c,B,info);
  }
  else if (L->Dtype == _Z) {  // calling the double precision complex routine.
#ifdef ARCOMP_H
    zgstrs(trans,L,U,perm_r,perm_c,B,info);
#endif
  }
  else {                      // calling the single precision complex routine.
#ifdef ARCOMP_H
    cgstrs(trans,L,U,perm_r,perm_c,B,info);
#endif
  }

} // gstrs.


// Create_CompCol_Matrix.

inline void Create_CompCol_Matrix(SuperMatrix* A, int m, int n, int nnz,
                                  double* a, int* irow, int* pcol,
                                  Stype_t S, Mtype_t M)
{

  dCreate_CompCol_Matrix(A,m,n,nnz,a,irow,pcol,S,_D,M);

} // Create_CompCol_Matrix (double).

inline void Create_CompCol_Matrix(SuperMatrix* A, int m, int n, int nnz,
                                  float* a, int* irow, int* pcol,
                                  Stype_t S, Mtype_t M)
{

  sCreate_CompCol_Matrix(A,m,n,nnz,a,irow,pcol,S,_S,M);

} // Create_CompCol_Matrix (float).

#ifdef ARCOMP_H

inline void Create_CompCol_Matrix(SuperMatrix* A, int m, int n, int nnz,
                                  arcomplex<double>* a, int* irow, int* pcol,
                                  Stype_t S, Mtype_t M)
{

  zCreate_CompCol_Matrix(A,m,n,nnz,(ldcomplex*)a,irow,pcol,S,_Z,M);

} // Create_CompCol_Matrix (complex<double>).

inline void Create_CompCol_Matrix(SuperMatrix* A, int m, int n, int nnz,
                                  arcomplex<float>* a, int* irow, int* pcol,
                                  Stype_t S, Mtype_t M)
{

  cCreate_CompCol_Matrix(A,m,n,nnz,(lscomplex*)a,irow,pcol,S,_C,M);

} // Create_CompCol_Matrix (complex<float>).

#endif // ARCOMP_H.


// Create_Dense_Matrix.

inline void Create_Dense_Matrix(SuperMatrix* A, int m, int n, double* x,
                                int ldx, Stype_t S, Mtype_t M)
{

  dCreate_Dense_Matrix(A,m,n,x,ldx,S,_D,M);

} // Create_Dense_Matrix (double).

inline void Create_Dense_Matrix(SuperMatrix* A, int m, int n, float* x,
                                int ldx, Stype_t S, Mtype_t M)
{

  sCreate_Dense_Matrix(A,m,n,x,ldx,S,_S,M);

} // Create_Dense_Matrix (float).

#ifdef ARCOMP_H

inline void Create_Dense_Matrix(SuperMatrix* A, int m, int n, arcomplex<double>* x,
                                int ldx, Stype_t S, Mtype_t M)
{

  zCreate_Dense_Matrix(A,m,n,(ldcomplex*)x,ldx,S,_Z,M);

} // Create_Dense_Matrix (complex<double>).

inline void Create_Dense_Matrix(SuperMatrix* A, int m, int n, arcomplex<float>* x,
                                int ldx, Stype_t S, Mtype_t M)
{

  cCreate_Dense_Matrix(A,m,n,(lscomplex*)x,ldx,S,_C,M);

} // Create_Dense_Matrix (complex<float>).

#endif // ARCOMP_H.

#endif // SUPERLUC_H
