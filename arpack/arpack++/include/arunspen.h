/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARUNSPen.h.
   Arpack++ class ARumNonSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUNSPEN_H
#define ARUNSPEN_H

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "umfpackc.h"
#include "arunsmat.h"


template<class TYPE, class FLOAT>
class ARumNonSymPencil
{

 protected:

  char                    part;
  ARumNonSymMatrix<TYPE>* A;
  ARumNonSymMatrix<TYPE>* B;
  ARumNonSymMatrix<TYPE>  AsB;
#ifdef ARCOMP_H
  ARumNonSymMatrix<arcomplex<FLOAT> > AsBc;
#endif

  virtual void Copy(const ARumNonSymPencil& other);

  void SparseSaxpy(TYPE a, TYPE x[], int xind[], int nx, TYPE y[],
                   int yind[], int ny, TYPE z[], int zind[], int& nz);

#ifdef ARCOMP_H
  void SparseSaxpy(arcomplex<FLOAT> a, FLOAT x[], int xind[], int nx,
                   FLOAT y[], int yind[], int ny,
                   arcomplex<FLOAT> z[], int zind[], int& nz);
#endif

  void SubtractAsB(TYPE sigma);

#ifdef ARCOMP_H
  void SubtractAsB(FLOAT sigmaR, FLOAT sigmaI);
#endif

 public:

#ifdef ARCOMP_H
  bool IsFactored() { return (AsB.IsFactored()||AsBc.IsFactored()); }
#else
  bool IsFactored() { return AsB.IsFactored(); }
#endif

  bool IsSymmetric() { return AsB.IsSymmetric(); }

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

  void DefineMatrices(ARumNonSymMatrix<TYPE>& Ap, ARumNonSymMatrix<TYPE>& Bp);

  ARumNonSymPencil() { part = 'N'; }
  // Short constructor that does nothing.

  ARumNonSymPencil(ARumNonSymMatrix<TYPE>& Ap, ARumNonSymMatrix<TYPE>& Bp);
  // Long constructor.

  ARumNonSymPencil(const ARumNonSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumNonSymPencil() { }
  // Destructor.

  ARumNonSymPencil& operator=(const ARumNonSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumNonSymPencil member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class TYPE, class FLOAT>
inline void ARumNonSymPencil<TYPE, FLOAT>::
Copy(const ARumNonSymPencil<TYPE, FLOAT>& other)
{

  part     = other.part;
  A        = other.A;
  B        = other.B;
  AsB      = other.AsB;
#ifdef ARCOMP_H
  AsBc     = other.AsBc;
#endif

} // Copy.


template<class TYPE, class FLOAT>
void ARumNonSymPencil<TYPE, FLOAT>::
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

} // SparseSaxpy (type).


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARumNonSymPencil<TYPE, FLOAT>::
SparseSaxpy(arcomplex<FLOAT> a, FLOAT x[], int xind[], int nx,
            FLOAT y[], int yind[], int ny, 
            arcomplex<FLOAT> z[], int zind[], int& nz)
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
      z[nz++]  = arcomplex<FLOAT>(y[iy++], 0.0);
      if (iy == ny) break;
    }
  }
  while (iy < ny) {
    zind[nz] = yind[iy];
    z[nz++]  = arcomplex<FLOAT>(y[iy++], 0.0);
  }
  while (ix < nx) {
    zind[nz] = xind[ix];
    z[nz++]  = arcomplex<FLOAT>(x[ix++], 0.0);
  }

} // SparseSaxpy (arcomplex<FLOAT>).
#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARumNonSymPencil<TYPE, FLOAT>::SubtractAsB(TYPE sigma)
{

  int i, acol, bcol, asbcol, scol;

  // Subtracting sigma*B from A.

  AsB.index[0] = 0;
  asbcol       = 0;

  for (i=0; i!=AsB.n; i++) {
    bcol = B->pcol[i];
    acol = A->pcol[i];
    SparseSaxpy(-sigma, &B->a[bcol], &B->irow[bcol], B->pcol[i+1]-bcol,
                &A->a[acol], &A->irow[acol], A->pcol[i+1]-acol,
                &AsB.value[asbcol], &AsB.index[asbcol+AsB.n+1], scol);
    asbcol += scol;
    AsB.index[i+1] = asbcol;
  }

  AsB.nnz = asbcol; 

  // Adding one to all elements of vector index
  // because the decomposition function was written in FORTRAN.

  for (i=0; i<=AsB.n+AsB.nnz; i++) AsB.index[i]++;

} // SubtractAsB (TYPE shift).


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARumNonSymPencil<TYPE, FLOAT>::SubtractAsB(FLOAT sigmaR, FLOAT sigmaI)
{

  int              i, acol, bcol, asbcol, scol;
  arcomplex<FLOAT> sigma;

  // Subtracting sigma*B from A.

  sigma         = arcomplex<FLOAT>(sigmaR, sigmaI);
  AsBc.index[0] = 0;
  asbcol        = 0;

  for (i=0; i!=AsBc.n; i++) {
    bcol = B->pcol[i];
    acol = A->pcol[i];
    SparseSaxpy(-sigma, &B->a[bcol], &B->irow[bcol], B->pcol[i+1]-bcol,
                &A->a[acol], &A->irow[acol], A->pcol[i+1]-acol,
                &AsBc.value[asbcol], &AsBc.index[asbcol+AsBc.n+1], scol);
    asbcol += scol;
    AsBc.index[i+1] = asbcol;
  }

  AsBc.nnz = asbcol;

  // Adding one to all elements of vector index
  // because the decomposition function was written in FORTRAN.

  for (i=0; i<=AsBc.n+AsBc.nnz; i++) AsBc.index[i]++;

} // SubtractAsB (arcomplex<FLOAT> shift).
#endif // ARCOMP_H


template<class TYPE, class FLOAT>
void ARumNonSymPencil<TYPE, FLOAT>::FactorAsB(TYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Defining matrix AsB.

  if (!AsB.IsDefined()) {

    int fillin = A->fillin > B->fillin ? A->fillin : B->fillin;
    AsB.DefineMatrix(A->ncols(), A->nzeros(), A->a, A->irow, 
                     A->pcol, A->threshold, fillin, 
                     (A->IsSymmetric() && B->IsSymmetric()), 
                     A->icntl[3], false);
    AsB.nnz = A->nzeros()+B->nzeros(); // temporary value.

  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsB.CreateStructure(); // AsB.nnz must be set to A->nzeros()+B->nzeros().

  // Subtracting sigma*B from A and storing the result on AsB.

  SubtractAsB(sigma);

  // Decomposing AsB.

  um2fa(AsB.n, AsB.index[AsB.n], 0, false, AsB.lvalue, AsB.lindex, AsB.value,
        AsB.index, AsB.keep, AsB.cntl, AsB.icntl, AsB.info, AsB.rinfo);

  // Handling errors.

  AsB.ThrowError();

  AsB.factored = true;

} // FactorAsB (TYPE shift).


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARumNonSymPencil<TYPE, FLOAT>::
FactorAsB(FLOAT sigmaR, FLOAT sigmaI, char partp)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymPencil::FactorAsB");
  }

  // Defining matrix AsB.

  if (!AsBc.IsDefined()) {

    part        = partp;
    int  fillin = A->fillin > B->fillin ? A->fillin : B->fillin;
    AsBc.DefineMatrix(A->ncols(), A->nzeros(), 0, 0,
                      A->pcol, A->threshold, fillin, 
                      (A->IsSymmetric() && B->IsSymmetric()), 
                      A->icntl[3], false);
    AsBc.nnz    = A->nzeros()+B->nzeros(); // temporary value.

  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsBc.CreateStructure(); // AsBc.nnz must be set to A->nzeros()+B->nzeros().

  // Subtracting sigma*B from A and storing the result on AsBc.

  SubtractAsB(sigmaR, sigmaI);

  // Decomposing AsB.

  um2fa(AsBc.n, AsBc.index[AsBc.n], 0, false, AsBc.lvalue, AsBc.lindex,
        AsBc.value, AsBc.index, AsBc.keep, AsBc.cntl, AsBc.icntl, 
        AsBc.info, AsBc.rinfo);

  // Handling errors.

  AsBc.ThrowError();

  AsBc.factored = true;

} // FactorAsB (arcomplex<FLOAT> shift).
#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARumNonSymPencil<TYPE, FLOAT>::MultInvBAv(TYPE* v, TYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  B->MultInvv(w, w);

} // MultInvBAv.


#ifdef ARCOMP_H

template<class TYPE, class FLOAT>
void ARumNonSymPencil<TYPE, FLOAT>::
MultInvAsBv(arcomplex<FLOAT>* v, arcomplex<FLOAT>* w)
{

  AsB.MultInvv((TYPE*)v,(TYPE*)w);

} // MultInvAsBv (arcomplex<FLOAT>).

#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARumNonSymPencil<TYPE, FLOAT>::MultInvAsBv(FLOAT* v, FLOAT* w)
{

  if (part == 'N') {    // shift is real.

    AsB.MultInvv((TYPE*)v,(TYPE*)w);

  }
  else {                // shift is complex.

#ifdef ARCOMP_H

    int              i;
    arcomplex<FLOAT> *tv, *tw;

    tv = new arcomplex<FLOAT>[AsBc.ncols()];
    tw = new arcomplex<FLOAT>[AsBc.ncols()];

    for (i=0; i!=AsBc.ncols(); i++) tv[i] = arcomplex<FLOAT>(v[i], 0.0);

    AsBc.MultInvv(tv, tw);

    if (part=='I') {
      for (i=0; i!=AsBc.ncols(); i++) w[i] = imag(tw[i]);
    }
    else {
      for (i=0; i!=AsBc.ncols(); i++) w[i] = real(tw[i]);
    }

    delete[] tv;
    delete[] tw;

#endif // ARCOMP_H.

  }

} // MultInvAsBv (FLOAT).


template<class TYPE, class FLOAT>
inline void ARumNonSymPencil<TYPE, FLOAT>::
DefineMatrices(ARumNonSymMatrix<TYPE>& Ap, ARumNonSymMatrix<TYPE>& Bp)
{

  A = &Ap;
  B = &Bp;

  if ((A->n != B->n)||(A->m != B->m)) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARumNonSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class TYPE, class FLOAT>
inline ARumNonSymPencil<TYPE, FLOAT>::
ARumNonSymPencil(ARumNonSymMatrix<TYPE>& Ap, ARumNonSymMatrix<TYPE>& Bp)
{

  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class TYPE, class FLOAT>
ARumNonSymPencil<TYPE, FLOAT>& ARumNonSymPencil<TYPE, FLOAT>::
operator=(const ARumNonSymPencil<TYPE, FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUNSPEN_H
