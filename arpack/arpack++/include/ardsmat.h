/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARDSMat.h.
   Arpack++ class ARdsSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "ardspen.h"

#ifndef ARDSMAT_H
#define ARDSMAT_H

#include <stddef.h>
#include "arch.h"
#include "armat.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"

template<class TYPE> class ARdsSymPencil;

template<class TYPE>
class ARdsSymMatrix: public ARMatrix<TYPE> {

  friend class ARdsSymPencil<TYPE>;

 protected:

  bool   factored;
  char   uplo;
  int    info;
  int*   ipiv;
  TYPE*  A;
  TYPE*  Ainv;

  void ClearMem(); 

  virtual void Copy(const ARdsSymMatrix& other);

  void SubtractAsI(TYPE sigma);

  void CreateStructure();

  void ThrowError();
  
 public:

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(TYPE sigma);

  void MultMv(TYPE* v, TYPE* w);

  void MultInvv(TYPE* v, TYPE* w);

  void DefineMatrix(int np, TYPE* Ap, char uplop = 'L');

  ARdsSymMatrix(): ARMatrix<TYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARdsSymMatrix(int np, TYPE* Ap, char uplop = 'L');
  // Long constructor.

  ARdsSymMatrix(const ARdsSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsSymMatrix() { ClearMem(); }
  // Destructor.

  ARdsSymMatrix& operator=(const ARdsSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARdsSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class TYPE>
inline void ARdsSymMatrix<TYPE>::ClearMem()
{ 

  if (factored) {
    delete[] Ainv;
    delete[] ipiv; 
    Ainv = NULL;
    ipiv = NULL;
  }

} // ClearMem.


template<class TYPE>
inline void ARdsSymMatrix<TYPE>::
Copy(const ARdsSymMatrix<TYPE>& other)
{

  // Copying very fundamental variables and user-defined parameters.

  m         = other.m;
  n         = other.n;
  defined   = other.defined;
  factored  = other.factored;
  uplo      = other.uplo;
  info      = other.info;
  A         = other.A;

  // Returning from here if "other" was not factored.

  if (!factored) return;

  // Copying vectors.

  Ainv = new TYPE[(n*n+n)/2];
  ipiv = new int[n];

  copy((n*n+n)/2, other.Ainv, 1, Ainv, 1);
  for (int i=0; i<n; i++) ipiv[i] = other.ipiv[i];

} // Copy.


template<class TYPE>
void ARdsSymMatrix<TYPE>::SubtractAsI(TYPE sigma)
{

  int i,j;

  // Copying A to Ainv.

  ::copy((n*n+n)/2 ,A, 1, Ainv, 1);

  // Subtracting sigma from diagonal elements.

  if (uplo=='L') {
    for (i=0, j=0; i<n; j+=(n-(i++))) Ainv[j] -= sigma;
  }
  else {
    for (i=0, j=0; i<n; j+=(++i)) Ainv[j] -= sigma;
  }

} // SubtractAsI.


template<class TYPE>
inline void ARdsSymMatrix<TYPE>::CreateStructure()
{

  ClearMem();
  Ainv = new TYPE[(n*n+n)/2];
  ipiv = new int[n];

} // CreateStructure.


template<class TYPE>
inline void ARdsSymMatrix<TYPE>::ThrowError()
{

  if (info < 0)  {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARdsSymMatrix::FactorA");
  }
  else if (info) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARdsSymMatrix::FactorA");
  }

} // ThrowError.


template<class TYPE>
void ARdsSymMatrix<TYPE>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsSymMatrix::FactorA");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to Ainv;

  ::copy((n*n+n)/2 ,A, 1, Ainv, 1);

  // Decomposing A.

  sptrf(&uplo, n, Ainv, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class TYPE>
void ARdsSymMatrix<TYPE>::FactorAsI(TYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  SubtractAsI(sigma);

  // Decomposing AsI.

  sptrf(&uplo, n, Ainv, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class TYPE>
void ARdsSymMatrix<TYPE>::MultMv(TYPE* v, TYPE* w)
{

  int i, j;

  TYPE  one  = (TYPE)0 + 1.0;
  TYPE  zero = (TYPE)0;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsSymMatrix::MultMv");
  }

  // Determining w = M.v (unfortunately, the BLAS does not 
  // have a routine that works with packed matrices).

  for (i=0; i<n; i++) w[i] = zero;

  if (uplo=='L') {
 
    for (i=0, j=0; i<n; j+=(n-(i++))) {
      w[i] += dot(n-i, &A[j], 1, &v[i], 1);
      axpy(n-i-1, v[i], &A[j+1], 1, &w[i+1], 1);
    }
 
  }  
  else { // uplo = 'U'

    for (i=0, j=0; i<n; j+=(++i)) {
      w[i] += dot(i+1, &A[j], 1, v, 1);
      axpy(i, v[i], &A[j], 1, w, 1);
    }
   
  }

} // MultMv.


template<class TYPE>
void ARdsSymMatrix<TYPE>::MultInvv(TYPE* v, TYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARdsSymMatrix::MultInvv");
  }

  // Overwritting w with v.

  copy(n, v, 1, w, 1);

  // Solving A.w = v (or AsI.w = v).

  sptrs(&uplo, n, 1, Ainv, ipiv, w, n, info);

  // Handling errors.

  ThrowError();

} // MultInvv.


template<class TYPE>
inline void ARdsSymMatrix<TYPE>::
DefineMatrix(int np, TYPE* Ap, char uplop)
{

  // Defining member variables.

  m         = np;
  n         = np;
  uplo      = uplop;
  A         = Ap;
  defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0; 

} // DefineMatrix.


template<class TYPE>
inline ARdsSymMatrix<TYPE>::
ARdsSymMatrix(int np, TYPE* Ap, char uplop) : ARMatrix<TYPE>(np)
{

  factored = false;
  DefineMatrix(np, Ap, uplop);

} // Long constructor.


template<class TYPE>
ARdsSymMatrix<TYPE>& ARdsSymMatrix<TYPE>::
operator=(const ARdsSymMatrix<TYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDSMAT_H
