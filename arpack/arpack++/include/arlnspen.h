/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARLNSPen.h.
   Arpack++ class ARluNonSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLNSPEN_H
#define ARLNSPEN_H

#include <stddef.h>
#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "superluc.h"
#include "arlspdef.h"
#include "arlutil.h"
#include "arlnsmat.h"


template<class TYPE, class FLOAT>
class ARluNonSymPencil
{

 protected:

  bool                    factored;
  int*                    permc;
  int*                    permr;
  char                    part;
  ARluNonSymMatrix<TYPE>* A;
  ARluNonSymMatrix<TYPE>* B;
  SuperMatrix             L;
  SuperMatrix             U;

  virtual void Copy(const ARluNonSymPencil& other);

  void ClearMem();

  void SparseSaxpy(TYPE a, TYPE x[], int xind[], int nx, TYPE y[],
                   int yind[], int ny, TYPE z[], int zind[], int& nz);

#ifdef ARCOMP_H
  void SparseSaxpy(arcomplex<FLOAT> a, FLOAT x[], int xind[], int nx,
                   FLOAT y[], int yind[], int ny, arcomplex<FLOAT> z[],
                   int zind[], int& nz);
#endif

  void SubtractAsB(int n, TYPE sigma, NCformat& A, NCformat& B, NCformat& AsB);

#ifdef ARCOMP_H
  void SubtractAsB(int n, FLOAT sigmaR, FLOAT sigmaI,
                   NCformat& A, NCformat& B, NCformat& AsB);
#endif

 public:

  bool IsFactored() { return factored; }

  void FactorAsB(TYPE sigma);

#ifdef ARCOMP_H
  void FactorAsB(FLOAT sigmaR, FLOAT sigmaI, char partp = 'R');
#endif

  void MultAv(TYPE* v, TYPE* w) { A->MultMv(v,w); }

  void MultBv(TYPE* v, TYPE* w) { B->MultMv(v,w); }

  void MultInvBAv(TYPE* v, TYPE* w);

#ifdef ARCOMP_H
  void MultInvAsBv(arcomplex<FLOAT>* v, arcomplex<FLOAT>* w);
#endif

  void MultInvAsBv(FLOAT* v, FLOAT* w);

  void DefineMatrices(ARluNonSymMatrix<TYPE>& Ap, ARluNonSymMatrix<TYPE>& Bp);

  ARluNonSymPencil();
  // Short constructor that does nothing.

  ARluNonSymPencil(ARluNonSymMatrix<TYPE>& Ap, ARluNonSymMatrix<TYPE>& Bp);
  // Long constructor.

  ARluNonSymPencil(const ARluNonSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluNonSymPencil() { ClearMem(); }
  // Destructor.

  ARluNonSymPencil& operator=(const ARluNonSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARluNonSymPencil member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class TYPE, class FLOAT>
inline void ARluNonSymPencil<TYPE, FLOAT>::
Copy(const ARluNonSymPencil<TYPE, FLOAT>& other)
{

  factored = other.factored;
  part     = other.part;
  A        = other.A;
  B        = other.B;

  // Throwing the original factorization away (this procedure 
  // is really awkward, but it is necessary because there
  // is no copy function for matrices L and U in the SuperLU 
  // library and it is not a good idea to do this kind of deep 
  // copy here).

  if (factored) {
    ArpackError(ArpackError::DISCARDING_FACTORS, "ARluNonSymPencil");
    factored = false;
  }

} // Copy.


template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::ClearMem()
{

  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree();
    delete[] permc;
    delete[] permr;
    permc = NULL;
    permr = NULL;
  }

} // ClearMem.


template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::
SparseSaxpy(TYPE a, TYPE x[], int xind[], int nx, TYPE y[],
            int yind[], int ny, TYPE z[], int zind[], int& nz)
// A strongly sequential (and inefficient) sparse saxpy algorithm.
{

  int ix, iy;

  nz = 0;
  if ((nx == 0) || (a == (TYPE)0)) {
    copy(ny,y,1,z,1);
    for (iy=0; iy!=ny; iy++) zind[iy] = yind[iy];
    nz = ny;
    return;
  }
  if (ny == 0) {
    copy(nx,x,1,z,1);
    scal(nx,a,z,1);
    for (ix=0; ix!=nx; ix++) zind[ix] = xind[ix];
    nz = nx;
    return;
  }
  ix = 0;
  iy = 0;
  while (true) {
    if (xind[ix] == yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++]+y[iy++];
      if ((ix == nx)||(iy == ny)) break;
    }
    else if (xind[ix] < yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++];
      if (ix == nx) break;
    }
    else {
      zind[nz] = yind[iy];
      z[nz++]  = y[iy++];
      if (iy == ny) break;
    }
  }
  while (iy < ny) {
    zind[nz] = yind[iy];
    z[nz++]  = y[iy++];
  }
  while (ix < nx) {
    zind[nz] = xind[ix];
    z[nz++]  = x[ix++];
  }

} // SparseSaxpy (TYPE).


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::
SparseSaxpy(arcomplex<FLOAT> a, FLOAT x[], int xind[], int nx, FLOAT y[],
            int yind[], int ny, arcomplex<FLOAT> z[], int zind[], int& nz)
// A strongly sequential (and inefficient) sparse saxpy algorithm.
{

  int ix, iy;

  nz = 0;
  if ((nx == 0) || (a == arcomplex<FLOAT>(0.0,0.0))) {
    for (iy=0; iy!=ny; iy++) {
      z[iy]    = arcomplex<FLOAT>(y[iy],0.0);
      zind[iy] = yind[iy];
    }
    nz = ny;
    return;
  }
  if (ny == 0) {
    for (ix=0; ix!=ny; ix++) {
      z[ix]    = a*arcomplex<FLOAT>(x[ix],0.0);
      zind[ix] = xind[ix];
    }
    nz = nx;
    return;
  }
  ix = 0;
  iy = 0;
  while (true) {
    if (xind[ix] == yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++]+y[iy++];
      if ((ix == nx)||(iy == ny)) break;
    }
    else if (xind[ix] < yind[iy]) {
      zind[nz] = xind[ix];
      z[nz++]  = a*x[ix++];
      if (ix == nx) break;
    }
    else {
      zind[nz] = yind[iy];
      z[nz++]  = arcomplex<FLOAT>(y[iy++],0.0);
      if (iy == ny) break;
    }
  }
  while (iy < ny) {
    zind[nz] = yind[iy];
    z[nz++]  = arcomplex<FLOAT>(y[iy++],0.0);
  }
  while (ix < nx) {
    zind[nz] = xind[ix];
    z[nz++]  = arcomplex<FLOAT>(x[ix++],0.0);
  }

} // SparseSaxpy (arcomplex<FLOAT>).
#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::
SubtractAsB(int n, TYPE sigma, NCformat& A, NCformat& B, NCformat& AsB)
{

  int   i, acol, bcol, asbcol, scol;
  TYPE* anzval;
  TYPE* bnzval;
  TYPE* asbnzval;

  // Telling the compiler that nzval must ve viewed as a vector of TYPE.

  anzval   = (TYPE*)A.nzval;
  bnzval   = (TYPE*)B.nzval;
  asbnzval = (TYPE*)AsB.nzval;

  // Subtracting sigma*B from A.

  AsB.colptr[0] = 0;
  asbcol = 0;

  for (i=0; i!=n; i++) {
    bcol = B.colptr[i];
    acol = A.colptr[i];
    SparseSaxpy(-sigma, &bnzval[bcol], &B.rowind[bcol], B.colptr[i+1]-bcol,
                &anzval[acol], &A.rowind[acol], A.colptr[i+1]-acol,
                &asbnzval[asbcol], &AsB.rowind[asbcol], scol);
    asbcol += scol;
    AsB.colptr[i+1] = asbcol;
  }

  AsB.nnz = AsB.colptr[n];

} // SubtractAsB (TYPE shift).


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::
SubtractAsB(int n, FLOAT sigmaR, FLOAT sigmaI,
            NCformat& A, NCformat& B, NCformat& AsB)
{

  int               i, acol, bcol, asbcol, scol;
  TYPE*             anzval;
  TYPE*             bnzval;
  arcomplex<FLOAT>* asbnzval;
  arcomplex<FLOAT>  sigma;

  // Telling the compiler that nzval must ve viewed as a vector of TYPE.

  anzval   = (TYPE*)A.nzval;
  bnzval   = (TYPE*)B.nzval;
  asbnzval = (arcomplex<FLOAT>*)AsB.nzval;

  // Subtracting sigma*B from A.

  sigma         = arcomplex<FLOAT>(sigmaR, sigmaI);
  AsB.colptr[0] = 0;
  asbcol        = 0;

  for (i=0; i!=n; i++) {
    bcol = B.colptr[i];
    acol = A.colptr[i];
    SparseSaxpy(-sigma, &bnzval[bcol], &B.rowind[bcol], B.colptr[i+1]-bcol,
                &anzval[acol], &A.rowind[acol], A.colptr[i+1]-acol,
                &asbnzval[asbcol], &AsB.rowind[asbcol], scol);
    asbcol += scol;
    AsB.colptr[i+1] = asbcol;
  }

  AsB.nnz = AsB.colptr[n];

} // SubtractAsB (arcomplex<FLOAT> shift).
#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::FactorAsB(TYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARluNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARluNonSymPencil::FactorAsB");
  }

  // Defining local variables.

  int         nnzi, info;
  int*        etree;
  int*        irowi;
  int*        pcoli;
  TYPE*       asb;
  SuperMatrix AsB;
  SuperMatrix AC;
  NCformat*   Astore;
  NCformat*   Bstore;
  NCformat*   AsBstore;

  // Deleting old versions of L, U, perm_r and perm_c.

  ClearMem();

  // Setting default values for gstrf parameters.

  FLOAT drop_tol        = 0.0;
  int   panel_size      = sp_ienv(1);
  int   relax           = sp_ienv(2);

  // Defining A and B format.

  Astore = (NCformat*)A->A.Store;
  Bstore = (NCformat*)B->A.Store;

  // Creating a temporary matrix AsB.

  nnzi  = Astore->nnz+Bstore->nnz;
  irowi = new int[nnzi];
  pcoli = new int[A->ncols()+1];
  asb   = new TYPE[nnzi];
  Create_CompCol_Matrix(&AsB, A->nrows(), A->ncols(), nnzi, asb,
                        irowi, pcoli, NC, GE);

  // Subtracting sigma*B from A and storing the result on AsB.

  AsBstore = (NCformat*)AsB.Store;
  SubtractAsB(A->ncols(), sigma, *Astore, *Bstore, *AsBstore);

  // Reserving memory for some vectors used in matrix decomposition.

  etree = new int[A->ncols()];
  if (permc == NULL) permc = new int[A->ncols()];
  if (permr == NULL) permr = new int[A->ncols()];

  // Defining LUStat.

  extern SuperLUStat_t SuperLUStat;
  StatInit(panel_size, relax);

  // Defining the column permutation of matrix AsB
  // (using minimum degree ordering on AsB'*AsB).

  get_perm_c(A->order, &AsB, permc);

  // Permuting columns of AsB and
  // creating the elimination tree of AsB'*AsB.

  sp_preorder("N", &AsB, permc, etree, &AC);

  // Decomposing AsB.

  gstrf("N",&AC, A->threshold, drop_tol, relax, panel_size, etree,
        NULL, 0, permr, permc, &L, &U, &info);

  // Deleting AC, AsB and etree.

  Destroy_CompCol_Permuted(&AC);
  Destroy_CompCol_Matrix(&AsB);
  delete[] etree;

  factored = (info == 0);

  // Handling errors.

  if (info < 0)  {              // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARluNonSymPencil::FactorAsB");
  }
  else if (info > A->ncols()) {  // Memory is not sufficient.
    throw ArpackError(ArpackError::MEMORY_OVERFLOW,
                      "ARluNonSymPencil::FactorAsB");
  }
  else if (info > 0) {          // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARluNonSymPencil::FactorAsB");
  }

} // FactorAsB (TYPE shift).


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::
FactorAsB(FLOAT sigmaR, FLOAT sigmaI, char partp)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARluNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARluNonSymPencil::FactorAsB");
  }

  // Defining local variables.

  int               nnzi, info;
  int*              etree;
  int*              irowi;
  int*              pcoli;
  arcomplex<FLOAT>* asb;
  SuperMatrix       AsB;
  SuperMatrix       AC;
  NCformat*         Astore;
  NCformat*         Bstore;
  NCformat*         AsBstore;

  // Deleting old versions of L, U, perm_r and perm_c.

  ClearMem();

  // Setting default values for gstrf parameters.

  FLOAT drop_tol        = 0.0;
  int   panel_size      = sp_ienv(1);
  int   relax           = sp_ienv(2);

  // Defining A and B format.

  Astore = (NCformat*)A->A.Store;
  Bstore = (NCformat*)B->A.Store;

  // Creating a temporary matrix AsB.

  part  = partp;
  nnzi  = Astore->nnz+Bstore->nnz;
  irowi = new int[nnzi];
  pcoli = new int[A->ncols()+1];
  asb   = new arcomplex<FLOAT>[nnzi];
  Create_CompCol_Matrix(&AsB, A->nrows(), A->ncols(), nnzi, asb,
                        irowi, pcoli, NC, GE);

  // Subtracting sigma*B from A and storing the result on AsB.

  AsBstore = (NCformat*)AsB.Store;
  SubtractAsB(A->ncols(), sigmaR, sigmaI, *Astore, *Bstore, *AsBstore);

  // Reserving memory for some vectors used in matrix decomposition.

  etree = new int[A->ncols()];
  if (permc == NULL) permc = new int[A->ncols()];
  if (permr == NULL) permr = new int[A->ncols()];

  // Defining LUStat.

  extern SuperLUStat_t SuperLUStat;
  StatInit(panel_size, relax);

  // Defining the column permutation of matrix AsB
  // (using minimum degree ordering on AsB'*AsB).

  get_perm_c(A->order, &AsB, permc);

  // Permuting columns of AsB and
  // creating the elimination tree of AsB'*AsB.

  sp_preorder("N", &AsB, permc, etree, &AC);

  // Decomposing AsB.

  gstrf("N",&AC, A->threshold, drop_tol, relax, panel_size, etree, NULL,
        0, permr, permc, &L, &U, &info);

  // Deleting AC, AsB and etree.

  Destroy_CompCol_Permuted(&AC);
  Destroy_CompCol_Matrix(&AsB);
  delete[] etree;

  factored = (info == 0);

  // Handling errors.

  if (info < 0)  {              // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARluNonSymPencil::FactorAsB");
  }
  else if (info > A->ncols()) {  // Memory is not sufficient.
    throw ArpackError(ArpackError::MEMORY_OVERFLOW,
                      "ARluNonSymPencil::FactorAsB");
  }
  else if (info > 0) {          // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARluNonSymPencil::FactorAsB");
  }

} // FactorAsB (arcomplex<FLOAT> shift).
#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::MultInvBAv(TYPE* v, TYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  B->MultInvv(w, w);

} // MultInvBAv.


#ifdef ARCOMP_H

template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::
MultInvAsBv(arcomplex<FLOAT>* v, arcomplex<FLOAT>* w)
{

  // Quitting the function if AsB was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARluNonSymPencil::MultInvAsBv");
  }

  // Solving AsB.w = v.

  int         info;
  SuperMatrix RHS;

  copy(A->nrows(), v, 1, w, 1);
  Create_Dense_Matrix(&RHS, A->nrows(), 1, w, A->nrows(), DN, GE);
  gstrs("N", &L, &U, permr, permc, &RHS, &info);

  delete RHS.Store;

} // MultInvAsBv (arcomplex<FLOAT>).

#endif


template<class TYPE, class FLOAT>
void ARluNonSymPencil<TYPE, FLOAT>::MultInvAsBv(FLOAT* v, FLOAT* w)
{

  // Quitting the function if AsB was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARluNonSymPencil::MultInvAsBv");
  }

  // Solving AsB.w = v.

  int         info;
  SuperMatrix RHS;

  if (part == 'N') {    // shift is real.

    copy(A->nrows(), v, 1, w, 1);
    Create_Dense_Matrix(&RHS, A->nrows(), 1, w, A->nrows(), DN, GE);
    gstrs("N", &L, &U, permr, permc, &RHS, &info);

  }
  else {                // shift is complex.

#ifdef ARCOMP_H

    int              i;
    arcomplex<FLOAT> *tv = new arcomplex<FLOAT>[A->ncols()];

    for (i=0; i!=A->ncols(); i++) tv[i] = arcomplex<FLOAT>(v[i],0.0);
    Create_Dense_Matrix(&RHS, A->ncols(), 1, tv, A->ncols(), DN, GE);
    gstrs("N", &L, &U, permr, permc, &RHS, &info);

    if (part=='I') {
      for (i=0; i!=A->ncols(); i++) w[i] = imag(tv[i]);
    }
    else {
      for (i=0; i!=A->ncols(); i++) w[i] = real(tv[i]);
    }

    delete[] tv;

#endif

  }

  delete RHS.Store;

} // MultInvAsBv (FLOAT).


template<class TYPE, class FLOAT>
inline void ARluNonSymPencil<TYPE, FLOAT>::
DefineMatrices(ARluNonSymMatrix<TYPE>& Ap, ARluNonSymMatrix<TYPE>& Bp)
{

  A     = &Ap;
  B     = &Bp;
  permc = NULL;
  permr = NULL;

  if ((A->n != B->n)||(A->m != B->m)) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARluNonSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class TYPE, class FLOAT>
inline ARluNonSymPencil<TYPE, FLOAT>::ARluNonSymPencil()
{
  
  factored = false; 
  part     = 'N'; 
  permr    = NULL;
  permc    = NULL;

} // Short constructor.


template<class TYPE, class FLOAT>
inline ARluNonSymPencil<TYPE, FLOAT>::
ARluNonSymPencil(ARluNonSymMatrix<TYPE>& Ap, ARluNonSymMatrix<TYPE>& Bp)
{

  factored = false;
  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class TYPE, class FLOAT>
ARluNonSymPencil<TYPE, FLOAT>& ARluNonSymPencil<TYPE, FLOAT>::
operator=(const ARluNonSymPencil<TYPE, FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLNSPEN_H
