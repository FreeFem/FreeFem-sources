/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARBSMat.h.
   Arpack++ class ARbdSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "arbspen.h"

#ifndef ARBSMAT_H
#define ARBSMAT_H

#include <stddef.h>
#include "arch.h"
#include "armat.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"

template<class TYPE> class ARbdSymPencil;

template<class TYPE>
class ARbdSymMatrix: public ARMatrix<TYPE> {

  friend class ARbdSymPencil<TYPE>;

 protected:

  bool   factored;
  char   uplo;
  int    nsdiag;
  int    lda;
  int    info;
  int*   ipiv;
  TYPE*  A;
  TYPE*  Ainv;

  void ClearMem(); 

  virtual void Copy(const ARbdSymMatrix& other);

  void ExpandA();

  void SubtractAsI(TYPE sigma);

  void CreateStructure();

  void ThrowError();
  
 public:

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(TYPE sigma);

  void MultMv(TYPE* v, TYPE* w);

  void MultInvv(TYPE* v, TYPE* w);

  void DefineMatrix(int np, int nsdiagp, TYPE* Ap, char uplop = 'L');

  ARbdSymMatrix(): ARMatrix<TYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARbdSymMatrix(int np, int nsdiagp, TYPE* Ap, char uplop = 'L');
  // Long constructor.

  ARbdSymMatrix(const ARbdSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARbdSymMatrix() { ClearMem(); }
  // Destructor.

  ARbdSymMatrix& operator=(const ARbdSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARbdSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class TYPE>
inline void ARbdSymMatrix<TYPE>::ClearMem()
{ 

  if (factored) {
    delete[] Ainv;
    delete[] ipiv; 
    Ainv = NULL;
    ipiv = NULL;
  }

} // ClearMem.


template<class TYPE>
inline void ARbdSymMatrix<TYPE>::
Copy(const ARbdSymMatrix<TYPE>& other)
{

  // Copying very fundamental variables and user-defined parameters.

  m         = other.m;
  n         = other.n;
  defined   = other.defined;
  factored  = other.factored;
  uplo      = other.uplo;
  nsdiag    = other.nsdiag;
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
void ARbdSymMatrix<TYPE>::ExpandA()
{

  int i;
 
  if (uplo == 'U') {

    // Copying the main diagonal of A to Ainv.

    copy(n, &A[nsdiag], nsdiag+1, &Ainv[2*nsdiag], lda);

    // Copying the superdiagonals of A to Ainv.

    for (i = 0; i < nsdiag; i++) {
      copy(n, &A[i], nsdiag+1, &Ainv[nsdiag+i], lda);
      copy(n-nsdiag+i, &A[i+(nsdiag-i)*(nsdiag+1)], nsdiag+1, 
           &Ainv[3*nsdiag-i], lda);
    }

  }
  else {

    // Copying the main diagonal of A to Ainv.

    copy(n, &A[0], nsdiag+1, &Ainv[2*nsdiag], lda);

    // Copying the subdiagonals of A to Ainv.

    for (i = 1; i <= nsdiag; i++) {
      copy(n, &A[i], nsdiag+1, &Ainv[2*nsdiag+i], lda);
      copy(n-i, &A[i], nsdiag+1, &Ainv[2*nsdiag-i+i*lda], lda);
    }

  }

} // ExpandA.


template<class TYPE>
void ARbdSymMatrix<TYPE>::SubtractAsI(TYPE sigma)
{

  // Copying A to Ainv.

  ExpandA();

  // Subtracting sigma from diagonal elements.

  for (int i=(2*nsdiag); i<(lda*n); i+=lda) Ainv[i] -= sigma; 

} // SubtractAsI.


template<class TYPE>
inline void ARbdSymMatrix<TYPE>::CreateStructure()
{

  ClearMem();
  Ainv = new TYPE[lda*n];
  ipiv = new int[n];

} // CreateStructure.


template<class TYPE>
inline void ARbdSymMatrix<TYPE>::ThrowError()
{

  if (info < 0)  {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARbdSymMatrix::FactorA");
  }
  else if (info) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARbdSymMatrix::FactorA");
  }

} // ThrowError.


template<class TYPE>
void ARbdSymMatrix<TYPE>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdSymMatrix::FactorA");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to Ainv;

  ExpandA();

  // Decomposing A.

  gbtrf(n, n, nsdiag, nsdiag, Ainv, lda, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class TYPE>
void ARbdSymMatrix<TYPE>::FactorAsI(TYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  SubtractAsI(sigma);

  // Decomposing AsI.

  gbtrf(n, n, nsdiag, nsdiag, Ainv, lda, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class TYPE>
void ARbdSymMatrix<TYPE>::MultMv(TYPE* v, TYPE* w)
{

  TYPE  one  = (TYPE)0 + 1.0;
  TYPE  zero = (TYPE)0;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARbdSymMatrix::MultMv");
  }

  // Determining w = M.v.

  sbmv(&uplo, n, nsdiag, one, A, nsdiag+1, v, 1, zero, w, 1);

} // MultMv.


template<class TYPE>
void ARbdSymMatrix<TYPE>::MultInvv(TYPE* v, TYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARbdSymMatrix::MultInvv");
  }

  // Overwritting w with v.

  copy(n, v, 1, w, 1);

  // Solving A.w = v (or AsI.w = v).

  gbtrs("N", n, nsdiag, nsdiag, 1, Ainv, lda, ipiv, w, m, info);

  // Handling errors.

  ThrowError();

} // MultInvv.


template<class TYPE>
inline void ARbdSymMatrix<TYPE>::
DefineMatrix(int np, int nsdiagp, TYPE* Ap, char uplop)
{

  // Defining member variables.

  m         = np;
  n         = np;
  nsdiag    = nsdiagp;
  lda       = 3*nsdiag+1;
  uplo      = uplop;
  A         = Ap;
  defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0; 

} // DefineMatrix.


template<class TYPE>
inline ARbdSymMatrix<TYPE>::
ARbdSymMatrix(int np, int nsdiagp, 
              TYPE* Ap, char uplop) : ARMatrix<TYPE>(np)
{

  factored = false;
  DefineMatrix(np, nsdiagp, Ap, uplop);

} // Long constructor.


template<class TYPE>
ARbdSymMatrix<TYPE>& ARbdSymMatrix<TYPE>::
operator=(const ARbdSymMatrix<TYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBSMAT_H
