/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARDNSPen.h.
   Arpack++ class ARdsNonSymPencil definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDNSPEN_H
#define ARDNSPEN_H

#include "arch.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"
#include "ardnsmat.h"


template<class TYPE, class FLOAT>
class ARdsNonSymPencil
{

 protected:

  char                    part;
  ARdsNonSymMatrix<TYPE>* A;
  ARdsNonSymMatrix<TYPE>* B;
  ARdsNonSymMatrix<TYPE>  AsB;
#ifdef ARCOMP_H
  ARdsNonSymMatrix<arcomplex<FLOAT> > AsBc;
#endif

  virtual void Copy(const ARdsNonSymPencil& other);

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

  void DefineMatrices(ARdsNonSymMatrix<TYPE>& Ap, ARdsNonSymMatrix<TYPE>& Bp);

  ARdsNonSymPencil() { part = 'N'; }
  // Short constructor that does nothing.

  ARdsNonSymPencil(ARdsNonSymMatrix<TYPE>& Ap, ARdsNonSymMatrix<TYPE>& Bp);
  // Long constructor.

  ARdsNonSymPencil(const ARdsNonSymPencil& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsNonSymPencil() { }
  // Destructor.

  ARdsNonSymPencil& operator=(const ARdsNonSymPencil& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARdsNonSymPencil member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class TYPE, class FLOAT>
inline void ARdsNonSymPencil<TYPE, FLOAT>::
Copy(const ARdsNonSymPencil<TYPE, FLOAT>& other)
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
void ARdsNonSymPencil<TYPE, FLOAT>::FactorAsB(TYPE sigma)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARdsNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARdsNonSymPencil::FactorAsB");
  }

  // Copying A to AsB if sigma = 0.

  if (sigma == (TYPE)0) {

    AsB = *A;
    if (!AsB.IsFactored()) AsB.FactorA();
    return;

  }

  // Defining matrix AsB.

  if (!AsB.IsDefined()) {
    AsB.DefineMatrix(A->ncols(), A->A);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsB.CreateStructure();

  // Subtracting sigma*B from A and storing the result on AsB.

  ::copy(A->m*A->n, A->A, 1, AsB.Ainv, 1); 
  axpy(A->m*A->n, -sigma, B->A, 1, AsB.Ainv, 1);

  // Decomposing AsB.

  getrf(AsB.m, AsB.n, AsB.Ainv, AsB.m, AsB.ipiv, AsB.info);

  // Handling errors.

  AsB.ThrowError();

  AsB.factored = true;

} // FactorAsB (TYPE shift).


#ifdef ARCOMP_H
template<class TYPE, class FLOAT>
void ARdsNonSymPencil<TYPE, FLOAT>::
FactorAsB(FLOAT sigmaR, FLOAT sigmaI, char partp)
{

  // Quitting the function if A and B were not defined.

  if (!(A->IsDefined()&&B->IsDefined())) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARdsNonSymPencil::FactorAsB");
  }

  // Quitting the function if A and B are not square.

  if ((A->nrows() != A->ncols()) || (B->nrows() != B->ncols())) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARdsNonSymPencil::FactorAsB");
  }

  // Defining matrix AsB.

  if (!AsBc.IsDefined()) {
    part = partp;
    AsBc.DefineMatrix(A->ncols(), 0);
  }

  // Reserving memory for some vectors used in matrix decomposition.

  AsBc.CreateStructure(); 

  // Subtracting sigma*B from A and storing the result on AsBc.

  arcomplex<FLOAT> sigma(sigmaR, sigmaI);
  for (int i=0; i<(A->m*A->n); i++) AsBc.Ainv[i] = A->A[i]-sigma*B->A[i];

  // Decomposing AsBc.

  getrf(AsBc.m, AsBc.n, AsBc.Ainv, AsBc.m, AsBc.ipiv, AsBc.info);

  // Handling errors.

  AsBc.ThrowError();

  AsBc.factored = true;

} // FactorAsB (arcomplex<FLOAT> shift).
#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARdsNonSymPencil<TYPE, FLOAT>::MultInvBAv(TYPE* v, TYPE* w)
{

  if (!B->IsFactored()) B->FactorA();

  A->MultMv(v, w);
  B->MultInvv(w, w);

} // MultInvBAv.


#ifdef ARCOMP_H

template<class TYPE, class FLOAT>
void ARdsNonSymPencil<TYPE, FLOAT>::
MultInvAsBv(arcomplex<FLOAT>* v, arcomplex<FLOAT>* w)
{

  AsB.MultInvv((TYPE*)v,(TYPE*)w);

} // MultInvAsBv (arcomplex<FLOAT>).

#endif // ARCOMP_H.


template<class TYPE, class FLOAT>
void ARdsNonSymPencil<TYPE, FLOAT>::MultInvAsBv(FLOAT* v, FLOAT* w)
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
inline void ARdsNonSymPencil<TYPE, FLOAT>::
DefineMatrices(ARdsNonSymMatrix<TYPE>& Ap, ARdsNonSymMatrix<TYPE>& Bp)
{

  A = &Ap;
  B = &Bp;

  if ((A->n != B->n)||(A->m != B->m)) {
    throw ArpackError(ArpackError::INCOMPATIBLE_SIZES,
                      "ARdsNonSymMatrix::DefineMatrices");
  }

} // DefineMatrices.


template<class TYPE, class FLOAT>
inline ARdsNonSymPencil<TYPE, FLOAT>::
ARdsNonSymPencil(ARdsNonSymMatrix<TYPE>& Ap, ARdsNonSymMatrix<TYPE>& Bp)
{

  DefineMatrices(Ap, Bp);

} // Long constructor.


template<class TYPE, class FLOAT>
ARdsNonSymPencil<TYPE, FLOAT>& ARdsNonSymPencil<TYPE, FLOAT>::
operator=(const ARdsNonSymPencil<TYPE, FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDNSPEN_H
