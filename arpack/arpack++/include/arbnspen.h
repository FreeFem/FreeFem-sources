/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARBNSPen.h.
   Arpack++ class ARbdNonSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBNSPEN_H
#define ARBNSPEN_H

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"
#include "arbnsmat.h"


template<class TYPE, class FLOAT>
class ARbdNonSymPencil
{

 protected:

  char                    part;
  ARbdNonSymMatrix<TYPE>* A;
  ARbdNonSymMatrix<TYPE>* B;
  ARbdNonSymMatrix<TYPE>  AsB;
#ifdef ARCOMP_H
  ARbdNonSymMatrix<arcomplex<FLOAT> > AsBc;
#endif

  int max(int a, int b) { return (a>b)?a:b; }

  int min(int a, int b) { return (a<b)?a:b; }

  void ComplexCopy(int n, FLOAT dx[], int incx, 
                   arcomplex<FLOAT> dy[], int incy);

  void ComplexAxpy(int n, arcomplex<FLOAT> da, TYPE dx[], 
                   int incx, arcomplex<FLOAT> dy[], int incy);

  virtual void Copy(const ARbdNonSymPencil& other);

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

  void DefineMatrices(ARbdNonSymMatrix<TYPE>& Ap, ARbdNonSymMatrix<TYPE>& Bp);

  ARbdNonSymPencil() { part = 'N'; }
  // Short constructor that does nothing.

  ARbdNonSymPencil(ARbdNonSymMatrix<TYPE>& Ap, ARbdNonSymMatrix<TYPE>& Bp);
  // Long constructor.

  ARbdNonSymPencil(const ARbdNonSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdNonSymPencil() { }
  // Destructor.

  ARbdNonSymPencil& operator=(const ARbdNonSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARbdNonSymPencil member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class TYPE, class FLOAT>
inline void ARbdNonSymPencil<TYPE, FLOAT>::
Copy(const ARbdNonSymPencil<TYPE, FLOAT>& other)
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
void ARbdNonSymPencil<TYPE, FLOAT>::
ComplexCopy(int n, FLOAT dx[], int incx, arcomplex<FLOAT> dy[], int incy)
{

  for (int ix=0, iy=0; ix<(n*incx); ix+=incx, iy+=incy) {
    dy[iy] = arcomplex<FLOAT>(dx[ix], 0.0);
  }

} // ComplexCopy.


template<class TYPE, class FLOAT>
void ARbdNonSymPencil<TYPE, FLOAT>::
ComplexAxpy(int n, arcomplex<FLOAT> da, TYPE dx[], int incx,
            arcomplex<FLOAT> dy[], int incy)
{

  for (int ix=0, iy=0; ix<(n*incx); ix+=incx, iy+=incy) {
    dy[iy] += da*dx[ix];
  }

} // ComplexAxpy.


template<class TYPE, class FLOAT>
void ARbdNonSymPencil<TYPE, FLOAT>::SubtractAsB(TYPE sigma)
{

  int  i, inca, incb, minL, minU, begB, begAsB;
  TYPE negsig;

  inca   = A->ndiagL+A->ndiagU+1;
  incb   = B->ndiagL+B->ndiagU+1;
  negsig = -sigma;

  // Expanding A.

  begAsB = AsB.ndiagL+AsB.ndiagU-A->ndiagU;
  for (i = 0; i < inca; i++) {
    copy(AsB.n, &A->A[i], inca, &AsB.Ainv[begAsB+i], AsB.lda);
  }

  // Transferring part of B (*(-sigma)) if AsB.ndiagU > A->ndiagU.

  if (A->ndiagU < AsB.ndiagU) {
    for (i = 0; i < AsB.ndiagU-A->ndiagU; i++) {
      copy(AsB.n, &B->A[i], incb, &AsB.Ainv[AsB.ndiagL+i], AsB.lda);
      scal(AsB.n, negsig, &AsB.Ainv[AsB.ndiagL+i], AsB.lda);
    }
  } 

  // Subtracting sigma*B from A.

  minL   = min(A->ndiagL, B->ndiagL);
  minU   = min(A->ndiagU, B->ndiagU);
  begB   = B->ndiagU-minU;
  begAsB = AsB.ndiagL+AsB.ndiagU-minU;

  for (i = 0; i < minL+minU+1; i++) {
    axpy(AsB.n, -sigma, &B->A[begB+i], incb, &AsB.Ainv[begAsB+i], AsB.lda);  
  }

  // Transferring part of B (*(-sigma)) if AsB.ndiagL > A->ndiagL.

  if (A->ndiagL < AsB.ndiagL) {
    begB   = B->ndiagU+1+minL;
    begAsB = AsB.ndiagL+AsB.ndiagU+1+minL;
    for (i = 0; i < AsB.ndiagL-A->ndiagL; i++) {
      copy(AsB.n, &B->A[begB+i], incb, &AsB.Ainv[begAsB+i], AsB.lda);
      scal(AsB.n, negsig, &AsB.Ainv[begAsB+i], AsB.lda);
    }
  } 

} // SubtractAsB (TYPE shift).


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARbdNonSymPencil<TYPE, FLOAT>::SubtractAsB(FLOAT sigmaR, FLOAT sigmaI)
{

  int              i, inca, incb, minL, minU, begB, begAsB;
  arcomplex<FLOAT> sigma;

  inca  = A->ndiagL+A->ndiagU+1;
  incb  = B->ndiagL+B->ndiagU+1;
  sigma = arcomplex<FLOAT>(sigmaR, sigmaI);

  // Expanding A.

  begAsB = AsBc.ndiagL+AsBc.ndiagU-A->ndiagU;
  for (i = 0; i < inca; i++) {
    ComplexCopy(AsBc.n,(FLOAT*)(&A->A[i]),inca,&AsBc.Ainv[begAsB+i],AsBc.lda);
  }

  // Transferring part of B (*(-sigma)) if AsBc.ndiagU > A->ndiagU.

  if (A->ndiagU < AsBc.ndiagU) {
    for (i = 0; i < AsBc.ndiagU-A->ndiagU; i++) {
      ComplexCopy(AsBc.n, (FLOAT*)(&B->A[i]), incb,
                  &AsBc.Ainv[AsBc.ndiagL+i], AsBc.lda);
      scal(AsBc.n, -sigma, &AsBc.Ainv[AsBc.ndiagL+i], AsBc.lda);
    }
  } 

  // Subtracting sigma*B from A.

  minL   = min(A->ndiagL, B->ndiagL);
  minU   = min(A->ndiagU, B->ndiagU);
  begB   = B->ndiagU-minU;
  begAsB = AsBc.ndiagL+AsBc.ndiagU-minU;

  for (i = 0; i < minL+minU+1; i++) {
    ComplexAxpy(AsBc.n, -sigma, &B->A[begB+i], incb, 
                &AsBc.Ainv[begAsB+i], AsBc.lda);  
  }

  // Transferring part of B (*(-sigma)) if AsBc.ndiagL > A->ndiagL.

  if (A->ndiagL < AsBc.ndiagL) {
    begB   = B->ndiagU+1+minL;
    begAsB = AsBc.ndiagL+AsBc.ndiagU+1+minL;
    for (i = 0; i < AsBc.ndiagL-A->ndiagL; i++) {
      ComplexCopy(AsBc.n, (FLOAT*)(&B->A[begB+i]), incb,
                  &AsBc.Ainv[begAsB+i], AsBc.lda);
      scal(AsBc.n, -sigma, &AsBc.Ainv[begAsB+i], AsBc.lda);
    }
  }

} // SubtractAsB (arcomplex<FLOAT> shift).
#endif // ARCOMP_H


template<class TYPE, class FLOAT>
void ARbdNonSymPencil<TYPE, FLOAT>::FactorAsB(TYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARbdNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARbdNonSymPencil::FactorAsB");
  }

  // Copying A to AsB if sigma = 0.

  if (sigma == (TYPE)0) {

    AsB = *A;
    if (!AsB.IsFactored()) AsB.FactorA();
    return;

  }

  // Defining matrix AsB.

  if (!AsB.IsDefined()) {
    AsB.DefineMatrix(A->ncols(), max(A->ndiagL, B->ndiagL),
                     max(A->ndiagU, B->ndiagU), A->A);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsB.CreateStructure(); 

  // Subtracting sigma*B from A and storing the result on AsB.

  SubtractAsB(sigma);

  // Decomposing AsB.

  gbtrf(AsB.n, AsB.n, AsB.ndiagL, AsB.ndiagU,
        AsB.Ainv, AsB.lda, AsB.ipiv, AsB.info);

  // Handling errors.

  AsB.ThrowError();

  AsB.factored = true;

} // FactorAsB (TYPE shift).


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARbdNonSymPencil<TYPE, FLOAT>::
FactorAsB(FLOAT sigmaR, FLOAT sigmaI, char partp)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARbdNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARbdNonSymPencil::FactorAsB");
  }

  // Defining matrix AsBc.

  if (!AsBc.IsDefined()) {
    part = partp;
    AsBc.DefineMatrix(A->ncols(), max(A->ndiagL,B->ndiagL),
                      max(A->ndiagU,B->ndiagU), 0);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsBc.CreateStructure();

  // Subtracting sigma*B from A and storing the result on AsBc.

  SubtractAsB(sigmaR, sigmaI);

  // Decomposing AsBc.

  gbtrf(AsBc.n, AsBc.n, AsBc.ndiagL, AsBc.ndiagU, 
        AsBc.Ainv, AsBc.lda, AsBc.ipiv, AsBc.info);

  // Handling errors.

  AsBc.ThrowError();

  AsBc.factored = true;

} // FactorAsB (arcomplex<FLOAT> shift).
#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARbdNonSymPencil<TYPE, FLOAT>::MultInvBAv(TYPE* v, TYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  B->MultInvv(w, w);

} // MultInvBAv.


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARbdNonSymPencil<TYPE, FLOAT>::
MultInvAsBv(arcomplex<FLOAT>* v, arcomplex<FLOAT>* w)
{

  AsB.MultInvv((TYPE*)v, (TYPE*)w);

} // MultInvAsBv (arcomplex<FLOAT>).
#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARbdNonSymPencil<TYPE, FLOAT>::MultInvAsBv(FLOAT* v, FLOAT* w)
{

  if (part == 'N') {    // shift is real.

    AsB.MultInvv((TYPE*)v, (TYPE*)w);

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
inline void ARbdNonSymPencil<TYPE, FLOAT>::
DefineMatrices(ARbdNonSymMatrix<TYPE>& Ap, ARbdNonSymMatrix<TYPE>& Bp)
{

  A = &Ap;
  B = &Bp;

} // DefineMatrices.


template<class TYPE, class FLOAT>
inline ARbdNonSymPencil<TYPE, FLOAT>::
ARbdNonSymPencil(ARbdNonSymMatrix<TYPE>& Ap, ARbdNonSymMatrix<TYPE>& Bp)
{

  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class TYPE, class FLOAT>
ARbdNonSymPencil<TYPE, FLOAT>& ARbdNonSymPencil<TYPE, FLOAT>::
operator=(const ARbdNonSymPencil<TYPE, FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBNSPEN_H
