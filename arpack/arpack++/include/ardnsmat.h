/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARDNSMat.h.
   Arpack++ class ARdsNonSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "ardnspen.h"

#ifndef ARDNSMAT_H
#define ARDNSMAT_H

#include <stddef.h>
#include "arch.h"
#include "armat.h"
#include "arerror.h"
#include "blas1c.h"
#include "lapackc.h"
#include "ardfmat.h"

template<class T, class S> class ARdsNonSymPencil;

template<class TYPE>
class ARdsNonSymMatrix: public ARMatrix<TYPE> {

  friend class ARdsNonSymPencil<double, double>;
  friend class ARdsNonSymPencil<float, float>;
#ifdef ARCOMP_H
  friend class ARdsNonSymPencil<arcomplex<double>, double>;
  friend class ARdsNonSymPencil<arcomplex<float>, float>;
#endif

 protected:

  bool              factored;
  int               info;
  int*              ipiv;
  TYPE*             A;
  TYPE*             Ainv;
  ARdfMatrix<TYPE>  mat;

  void ClearMem(); 

  virtual void Copy(const ARdsNonSymMatrix& other);

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

  void DefineMatrix(int np, TYPE* Ap);

  void DefineMatrix(int mp, int np, TYPE* Ap);

  ARdsNonSymMatrix(): ARMatrix<TYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARdsNonSymMatrix(int np, TYPE* Ap);
  // Long constructor (square matrix).

  ARdsNonSymMatrix(int mp, int np, TYPE* Ap);
  // Long constructor (rectangular matrix).

  ARdsNonSymMatrix(char* file, int blksizep = 0);
  // Long constructor (Matrix stored in a file).

  ARdsNonSymMatrix(const ARdsNonSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARdsNonSymMatrix() { ClearMem(); }
  // Destructor.

  ARdsNonSymMatrix& operator=(const ARdsNonSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARdsNonSymMatrix member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class TYPE>
inline void ARdsNonSymMatrix<TYPE>::ClearMem()
{ 

  if (factored) {
    delete[] Ainv;
    delete[] ipiv; 
    Ainv = NULL;
    ipiv = NULL;
  }

} // ClearMem.


template<class TYPE>
inline void ARdsNonSymMatrix<TYPE>::
Copy(const ARdsNonSymMatrix<TYPE>& other)
{

  // Copying very fundamental variables and user-defined parameters.

  m         = other.m;
  n         = other.n;
  defined   = other.defined;
  factored  = other.factored;
  info      = other.info;
  A         = other.A;

  // Copying mat.

  if (other.mat.IsDefined()) {
    mat.Define(other.mat.Filename(),other.mat.BlockSize());
  }

  // Returning from here if "other" was not factored.

  if (!factored) return;

  // Copying vectors.

  Ainv = new TYPE[m*n];
  ipiv = new int[n];

  copy(m*n, other.Ainv, 1, Ainv, 1);
  for (int i=0; i<n; i++) ipiv[i] = other.ipiv[i];

} // Copy.


template<class TYPE>
inline void ARdsNonSymMatrix<TYPE>::CreateStructure()
{

  ClearMem();
  Ainv = new TYPE[m*n];
  ipiv = new int[n];

} // CreateStructure.


template<class TYPE>
inline void ARdsNonSymMatrix<TYPE>::ThrowError()
{

  if (info < 0)  {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARdsNonSymMatrix::FactorA");
  }
  else if (info) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARdsNonSymMatrix::FactorA");
  }

} // ThrowError.


template<class TYPE>
void ARdsNonSymMatrix<TYPE>::FactorA()
{

  // Quitting the function if A was not defined or is rectangular.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsNonSymMatrix::FactorA");
  }

  if (m!=n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARdsNonSymMatrix::FactorA");
  }

  if (mat.IsOutOfCore()) {
    throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                      "ARdsNonSymMatrix::FactorA");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to Ainv;

  ::copy(m*n, A, 1, Ainv, 1);

  // Decomposing A.

  getrf(m, n, Ainv, m, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class TYPE>
void ARdsNonSymMatrix<TYPE>::FactorAsI(TYPE sigma)
{

  // Quitting the function if A was not defined or is rectangular.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARdsNonSymMatrix::FactorAsI");
  }

  if (m!=n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARdsNonSymMatrix::FactorAsI");
  }

  if (mat.IsOutOfCore()) {
    throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                      "ARdsNonSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  ::copy(m*n,A,1,Ainv,1);
  for (int i=0; i<(m*n); i+=m+1) Ainv[i]-=sigma;

  // Decomposing AsI.

  getrf(m, n, Ainv, m, ipiv, info);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class TYPE>
void ARdsNonSymMatrix<TYPE>::MultMv(TYPE* v, TYPE* w)
{

  int   i;
  TYPE* t;
  TYPE  one;
  TYPE  zero;

  one  = (TYPE)0 + 1.0;
  zero = (TYPE)0;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsNonSymMatrix::MultMv");
  }

  // Determining w = M.v.

  if (mat.IsOutOfCore()) {

    if (m>n) { 

      // Matrix is "tall".

      mat.Rewind();
      for (i=0; i<mat.NBlocks(); i++) {
        mat.ReadBlock();
        gemv("N", mat.RowsInMemory(), n, one, mat.Entries(), 
             mat.RowsInMemory(), v, 1, zero, &w[mat.FirstIndex()], 1);
      }

    }
    else {

      // Matrix is "fat".

      mat.Rewind();
      t = new TYPE[mat.ColsInMemory()];
      for (i=0; i<m; i++) w[i] = zero;
      for (i=0; i<mat.NBlocks(); i++) {
        mat.ReadBlock();
        gemv("N", m, mat.ColsInMemory(), one, mat.Entries(), 
             m, &v[mat.FirstIndex()], 1, zero, t, 1);
        axpy(m, one, t, 1, w, 1); 
      }
      delete[] t;

    }

  }
  else {

    gemv("N", m, n, one, A, m, v, 1, zero, w, 1);

  }

} // MultMv.


template<class TYPE>
void ARdsNonSymMatrix<TYPE>::MultMtv(TYPE* v, TYPE* w)
{

  int   i;
  TYPE* t;
  TYPE  one;   
  TYPE  zero; 

  one  = (TYPE)0 + 1.0;
  zero = (TYPE)0;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARdsNonSymMatrix::MultMtv");
  }

  // Determining w = M'.v.

  if (mat.IsOutOfCore()) {

    if (m<=n) { 

      // Matrix is "fat".

      mat.Rewind();
      for (i=0; i<mat.NBlocks(); i++) {
        mat.ReadBlock();
        gemv("T", m, mat.ColsInMemory(), one, mat.Entries(), 
             m, v, 1, zero, &w[mat.FirstIndex()], 1);
      }

    }
    else {

      // Matrix is "tall".

      mat.Rewind();
      t = new TYPE[mat.ColsInMemory()];
      for (i=0; i<m; i++) w[i] = zero;
      for (i=0; i<mat.NBlocks(); i++) {
        mat.ReadBlock();
        gemv("T", mat.RowsInMemory(), n, one, mat.Entries(), 
             mat.RowsInMemory(), &v[mat.FirstIndex()], 1, zero, t, 1);
        axpy(mat.RowsInMemory(), one, t, 1, w, 1); 
      }
      delete[] t;

    }

  }
  else {

    gemv("T", m, n, one, A, m, v, 1, zero, w, 1);

  }


} // MultMtv.


template<class TYPE>
void ARdsNonSymMatrix<TYPE>::MultMtMv(TYPE* v, TYPE* w)
{

  int  i;
  TYPE *t, *s;
  TYPE one;   
  TYPE zero; 

  one  = (TYPE)0 + 1.0;
  zero = (TYPE)0;

  if (mat.IsOutOfCore() && (m>n)) {

    // Special code for "tall" matrices.

    t = new TYPE[mat.BlockSize()];
    s = new TYPE[n];

    mat.Rewind();
    for (i=0; i<n; i++) w[i] = zero;
    for (i=0; i<mat.NBlocks(); i++) {
      mat.ReadBlock();
      gemv("N", mat.RowsInMemory(), n, one, mat.Entries(), 
           mat.RowsInMemory(), v, 1, zero, t, 1);
      gemv("T", mat.RowsInMemory(), n, one, mat.Entries(), 
           mat.RowsInMemory(), t, 1, zero, s, 1);
      axpy(n, one, s, 1, w, 1); 

    }

    delete[] t;
    delete[] s;

  }
  else {

    t = new TYPE[m];

    MultMv(v,t);
    MultMtv(t,w);

    delete[] t;

  }


} // MultMtMv.


template<class TYPE>
void ARdsNonSymMatrix<TYPE>::MultMMtv(TYPE* v, TYPE* w)
{

  int  i;
  TYPE *t, *s;
  TYPE one;   
  TYPE zero; 

  one  = (TYPE)0 + 1.0;
  zero = (TYPE)0;

  if (mat.IsOutOfCore() && (m<=n)) {

    // Special code for "fat" matrices.

    t = new TYPE[mat.BlockSize()];
    s = new TYPE[m];

    mat.Rewind();
    for (i=0; i<m; i++) w[i] = zero;
    for (i=0; i<mat.NBlocks(); i++) {
      mat.ReadBlock();
      gemv("T", m, mat.ColsInMemory(), one, mat.Entries(), 
           m, v, 1, zero, t, 1);
      gemv("N", m, mat.ColsInMemory(), one, mat.Entries(), 
           m, t, 1, zero, s, 1);
      axpy(m, one, s, 1, w, 1); 

    }

    delete[] t;
    delete[] s;

  }
  else {

    t = new TYPE[n];

    MultMtv(v,t);
    MultMv(t,w);

    delete[] t;

  }

} // MultMMtv.


template<class TYPE>
void ARdsNonSymMatrix<TYPE>::Mult0MMt0v(TYPE* v, TYPE* w)
{

  MultMv(&v[m],w);
  MultMtv(v,&w[m]);

} // Mult0MMt0v.


template<class TYPE>
void ARdsNonSymMatrix<TYPE>::MultInvv(TYPE* v, TYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARdsNonSymMatrix::MultInvv");
  }

  // Overwritting w with v.

  copy(n, v, 1, w, 1);

  // Solving A.w = v (or AsI.w = v).

  getrs("N", n, 1, Ainv, m, ipiv, w, m, info);

  // Handling errors.

  ThrowError();

} // MultInvv.


template<class TYPE>
inline void ARdsNonSymMatrix<TYPE>::
DefineMatrix(int np, TYPE* Ap)
{

  // Defining member variables.

  n         = np;
  m         = np;
  A         = Ap;
  defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0; 

} // DefineMatrix (square).


template<class TYPE>
inline void ARdsNonSymMatrix<TYPE>::
DefineMatrix(int mp, int np, TYPE* Ap)
{

  // Defining member variables.

  m         = mp;
  n         = np;
  A         = Ap;
  defined   = true;
  Ainv      = NULL;
  ipiv      = NULL;
  info      = 0;

} // DefineMatrix (rectangular).


template<class TYPE>
inline ARdsNonSymMatrix<TYPE>::
ARdsNonSymMatrix(int np, TYPE* Ap) : ARMatrix<TYPE>(np)
{

  factored = false;
  DefineMatrix(np, Ap);

} // Long constructor (square matrix).


template<class TYPE>
inline ARdsNonSymMatrix<TYPE>::
ARdsNonSymMatrix(int mp, int np, TYPE* Ap) : ARMatrix<TYPE>(mp, np)
{

  factored = false;
  DefineMatrix(mp, np, Ap);

} // Long constructor (rectangular matrix).


template<class TYPE>
ARdsNonSymMatrix<TYPE>::ARdsNonSymMatrix(char* file, int blksizep)
{

  factored = false;

  try {
    mat.Define(file, blksizep);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARdsNonSymMatrix");
  }

  if (mat.NCols() == mat.NRows()) {
    DefineMatrix(mat.NCols(), (TYPE*)mat.Entries());
  }
  else {                             
    DefineMatrix(mat.NRows(), mat.NCols(), (TYPE*)mat.Entries());
  }

} // Long constructor (Matrix stored in a file).


template<class TYPE>
ARdsNonSymMatrix<TYPE>& ARdsNonSymMatrix<TYPE>::
operator=(const ARdsNonSymMatrix<TYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDNSMAT_H
