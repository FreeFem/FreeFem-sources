/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARBNSMat.h.
   Arpack++ class ARbdNonSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "arbnspen.h"

#ifndef ARBNSMAT_H
#define ARBNSMAT_H

#include <stddef.h>
#include "arch.h"
#include "armat.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"

template<class T, class S> class ARbdNonSymPencil;

template<class TYPE>
class ARbdNonSymMatrix: public ARMatrix<TYPE> {

  friend class ARbdNonSymPencil<double, double>;
  friend class ARbdNonSymPencil<float, float>;
#ifdef ARCOMP_H
  friend class ARbdNonSymPencil<arcomplex<double>, double>;
  friend class ARbdNonSymPencil<arcomplex<float>, float>;
#endif

 protected:

  bool   factored;
  int    ndiagL;
  int    ndiagU;
  int    lda;
  int    info;
  int*   ipiv;
  TYPE*  A;
  TYPE*  Ainv;

  void ClearMem(); 

  virtual void Copy(const ARbdNonSymMatrix& other);

  void ExpandA();

  void SubtractAsI(TYPE sigma);

  void CreateStructure();

  void ThrowError();
  
 public:

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(TYPE sigma);

  void MultMv(TYPE* v, TYPE* w);

  void MultMtv(TYPE* v, TYPE* w);

  void MultMtMv(TYPE* v, TYPE* w);

  void MultMMtv(TYPE* v, TYPE* w);

  void Mult0MMt0v(TYPE* v, TYPE* w);

  void MultInvv(TYPE* v, TYPE* w);

  void DefineMatrix(int np, int ndiagLp, int ndiagUp, TYPE* Ap);

  ARbdNonSymMatrix(): ARMatrix<TYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARbdNonSymMatrix(int np, int ndiagLp, int ndiagUp, TYPE* Ap);
  // Long constructor.

  ARbdNonSymMatrix(const ARbdNonSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdNonSymMatrix() { ClearMem(); }
  // Destructor.

  ARbdNonSymMatrix& operator=(const ARbdNonSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARbdNonSymMatrix member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class TYPE>
inline void ARbdNonSymMatrix<TYPE>::ClearMem()
{ 

  if (factored) {
    delete[] Ainv;
    delete[] ipiv; 
    Ainv = NULL;
    ipiv = NULL;
  }

} // ClearMem.


template<class TYPE>
inline void ARbdNonSymMatrix<TYPE>::
Copy(const ARbdNonSymMatrix<TYPE>& other)
{

  // Copying very fundamental variables and user-defined parameters.

  m         = other.m;
  n         = other.n;
  defined   = other.defined;
  factored  = other.factored;
  ndiagL    = other.ndiagL;
  ndiagU    = other.ndiagU;
  lda       = other.lda;
  info      = other.info;
  A         = other.A;

  // Returning from here if "other" was not factored.

  if (!factored) return;

  // Copying vectors.

  Ainv = new TYPE[n*lda];
  ipiv = new int[n];

  copy(n*lda, other.Ainv, 1, Ainv, 1);
  for (int i=0; i<n; i++) ipiv[i] = other.ipiv[i];

} // Copy.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::ExpandA()
{

  int i, inca;
 
  // Copying A to Ainv.

  inca = ndiagL+ndiagU+1;
  for (i = 0; i < inca; i++) {
    copy(n, &A[i], inca, &Ainv[ndiagL+i], lda);
  }

} // ExpandA.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::SubtractAsI(TYPE sigma)
{

  // Copying A to Ainv.

  ExpandA();

  // Subtracting sigma from diagonal elements.

  for (int i=(ndiagL+ndiagU); i<(lda*n); i+=lda) Ainv[i] -= sigma; 

} // SubtractAsI.


template<class TYPE>
inline void ARbdNonSymMatrix<TYPE>::CreateStructure()
{

  ClearMem();
  Ainv = new TYPE[lda*n];
  ipiv = new int[n];

} // CreateStructure.


template<class TYPE>
inline void ARbdNonSymMatrix<TYPE>::ThrowError()
{

  if (info < 0)  {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARbdNonSymMatrix::FactorA");
  }
  else if (info) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARbdNonSymMatrix::FactorA");
  }

} // ThrowError.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdNonSymMatrix::FactorA");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to Ainv;

  ExpandA();

  // Decomposing A.

  gbtrf(n, n, ndiagL, ndiagU, Ainv, lda, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::FactorAsI(TYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARbdNonSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  SubtractAsI(sigma);

  // Decomposing AsI.

  gbtrf(n, n, ndiagL, ndiagU, Ainv, lda, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::MultMv(TYPE* v, TYPE* w)
{

  TYPE  one;
  TYPE  zero;

  one  = (TYPE)0 + 1.0;
  zero = (TYPE)0;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdNonSymMatrix::MultMv");
  }

  // Determining w = M.v.

  gbmv("N", m, n, ndiagL, ndiagU, one, A,
       ndiagL+ndiagU+1, v, 1, zero, w, 1);

} // MultMv.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::MultMtv(TYPE* v, TYPE* w)
{

  TYPE  one;   
  TYPE  zero; 

  one  = (TYPE)0 + 1.0;
  zero = (TYPE)0;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdNonSymMatrix::MultMtv");
  }

  // Determining w = M'.v.

  gbmv("T", m, n, ndiagL, ndiagU, one, A,
       ndiagL+ndiagU+1, v, 1, zero, w, 1);   

} // MultMtv.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::MultMtMv(TYPE* v, TYPE* w)
{

  TYPE* t = new TYPE[m];

  MultMv(v,t);
  MultMtv(t,w);

  delete[] t;

} // MultMtMv.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::MultMMtv(TYPE* v, TYPE* w)
{

  TYPE* t = new TYPE[n];

  MultMtv(v,t);
  MultMv(t,w);

  delete[] t;

} // MultMMtv.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::Mult0MMt0v(TYPE* v, TYPE* w)
{

  MultMv(&v[m],w);
  MultMtv(v,&w[m]);

} // Mult0MMt0v.


template<class TYPE>
void ARbdNonSymMatrix<TYPE>::MultInvv(TYPE* v, TYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARbdNonSymMatrix::MultInvv");
  }

  // Overwritting w with v.

  copy(n, v, 1, w, 1);

  // Solving A.w = v (or AsI.w = v).

  gbtrs("N", n, ndiagL, ndiagU, 1, Ainv, lda, ipiv, w, m, info);

  // Handling errors.

  ThrowError();

} // MultInvv.


template<class TYPE>
inline void ARbdNonSymMatrix<TYPE>::
DefineMatrix(int np, int ndiagLp, int ndiagUp, TYPE* Ap)
{

  // Defining member variables.

  m         = np;
  n         = np;
  ndiagL    = ndiagLp;
  ndiagU    = ndiagUp;
  lda       = 2*ndiagL+ndiagU+1;
  A         = Ap;
  defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0; 

} // DefineMatrix.


template<class TYPE>
inline ARbdNonSymMatrix<TYPE>::
ARbdNonSymMatrix(int np, int ndiagLp, 
                 int ndiagUp, TYPE* Ap) : ARMatrix<TYPE>(np)
{

  factored = false;
  DefineMatrix(np, ndiagLp, ndiagUp, Ap);

} // Long constructor.


template<class TYPE>
ARbdNonSymMatrix<TYPE>& ARbdNonSymMatrix<TYPE>::
operator=(const ARbdNonSymMatrix<TYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBNSMAT_H
