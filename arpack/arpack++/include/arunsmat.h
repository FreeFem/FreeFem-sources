/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARUNSMat.h.
   Arpack++ class ARumNonSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "arunspen.h"

#ifndef ARUNSMAT_H
#define ARUNSMAT_H

#include <stddef.h>
#include "arch.h"
#include "armat.h"
#include "arhbmat.h"
#include "arerror.h"
#include "blas1c.h"
#include "umfpackc.h"

template<class T, class S> class ARumNonSymPencil;

template<class TYPE>
class ARumNonSymMatrix: public ARMatrix<TYPE> {

  friend class ARumNonSymPencil<double, double>;
  friend class ARumNonSymPencil<float, float>;
#ifdef ARCOMP_H
  friend class ARumNonSymPencil<arcomplex<double>, double>;
  friend class ARumNonSymPencil<arcomplex<float>, float>;
#endif

 protected:

  bool   factored;
  int    fillin;
  int    nnz;
  int    lvalue;
  int    lindex;
  int    keep[20];
  int    icntl[20];
  int    info[40];
  int*   irow;
  int*   pcol;
  int*   index;
  double threshold;
  TYPE   cntl[10];
  TYPE   rinfo[20];
  TYPE*  a;
  TYPE*  value;
  ARhbMatrix<int, TYPE> mat;

  bool DataOK();

  void ClearMem();

  virtual void Copy(const ARumNonSymMatrix& other);

  void SubtractAsI(TYPE sigma);

  void CreateStructure();

  void ThrowError();

 public:

  int nzeros() { return nnz; }

  int  FillFact() { return fillin; }

  bool IsSymmetric() { return bool(icntl[5]); }

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(TYPE sigma);

  void MultMv(TYPE* v, TYPE* w);

  void MultMtv(TYPE* v, TYPE* w);

  void MultMtMv(TYPE* v, TYPE* w);

  void MultMMtv(TYPE* v, TYPE* w);

  void Mult0MMt0v(TYPE* v, TYPE* w);

  void MultInvv(TYPE* v, TYPE* w);

  void DefineMatrix(int np, int nnzp, TYPE* ap, int* irowp,
                    int* pcolp, double thresholdp = 0.1,
                    int fillinp = 9, bool simest = false,
                    bool reducible = true, bool check = true); // Square.

  void DefineMatrix(int mp, int np, int nnzp, TYPE* ap,
                    int* irowp, int* pcolp);                   // Rectangular.

  ARumNonSymMatrix(): ARMatrix<TYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARumNonSymMatrix(int np, int nnzp, TYPE* ap, int* irowp,
                   int* pcolp, double thresholdp = 0.1,
                   int fillinp = 9, bool simest = false,
                   bool reducible = true, bool check = true);
  // Long constructor (square matrix).

  ARumNonSymMatrix(int mp, int np, int nnzp, TYPE* ap,
                   int* irowp, int* pcolp);
  // Long constructor (rectangular matrix).

  ARumNonSymMatrix(char* name, double thresholdp = 0.1,
                   int fillinp = 9, bool simest = false,
                   bool reducible = true, bool check = true);
  // Long constructor (Harwell-Boeing file).

  ARumNonSymMatrix(const ARumNonSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumNonSymMatrix() { ClearMem(); }
  // Destructor.

  ARumNonSymMatrix& operator=(const ARumNonSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumNonSymMatrix member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class TYPE>
bool ARumNonSymMatrix<TYPE>::DataOK()
{

  int i, j, k;

  // Checking if pcol is in ascending order.

  i = 0;
  while ((i!=n)&&(pcol[i]<=pcol[i+1])) i++;
  if (i!=n) return false;

  // Checking if irow components are in order and within bounds.

  for (i=0; i!=n; i++) {
    j = pcol[i];
    k = pcol[i+1]-1;
    if (j<=k) {
      if ((irow[j]<0)||(irow[k]>=n)) return false;
      while ((j!=k)&&(irow[j]<irow[j+1])) j++;
      if (j!=k) return false;
    }  
  }

  return true;

} // DataOK.


template<class TYPE>
inline void ARumNonSymMatrix<TYPE>::ClearMem()
{

  if (factored) {
    delete[] value;
    delete[] index;
    value = NULL;
    index = NULL;
  }

} // ClearMem.


template<class TYPE>
inline void ARumNonSymMatrix<TYPE>::
Copy(const ARumNonSymMatrix<TYPE>& other)
{

  // Local variable.

  int i;

  // Copying very fundamental variables and user-defined parameters.

  m         = other.m;
  n         = other.n;
  defined   = other.defined;
  factored  = other.factored;
  fillin    = other.fillin;
  nnz       = other.nnz;
  lvalue    = other.lvalue;
  lindex    = other.lindex;
  irow      = other.irow;
  pcol      = other.pcol;
  a         = other.a;
  threshold = other.threshold;

  // Returning from here if "other" was not initialized.

  if (!defined) return;

  // Copying arrays with static dimension.

  for (i=0; i<20; i++) keep[i]  = other.keep[i];
  for (i=0; i<20; i++) icntl[i] = other.icntl[i];
  for (i=0; i<40; i++) info[i]  = other.info[i];
  for (i=0; i<10; i++) cntl[i]  = other.cntl[i];
  for (i=0; i<20; i++) rinfo[i] = other.rinfo[i];

  // Returning from here if "other" was not factored.

  if (!factored) return;

  value = new TYPE[lvalue];
  index = new int[lindex];

  for (i=0; i<lindex; i++) index[i] = other.index[i];
  copy(lvalue, other.value, 1, value, 1);
  
} // Copy.


template<class TYPE>
void ARumNonSymMatrix<TYPE>::SubtractAsI(TYPE sigma)
{

  int i, j, k, ki, end;

  // Subtracting sigma from diagonal elements.

  k        = 0;
  ki       = n+1;
  index[0] = 1;

  for (i=0; i!=n; i++) {

    j = pcol[i];
    end = pcol[i+1];

    // Copying superdiagonal elements of column i.

    while ((irow[j] < i)&&(j < end)) {
      value[k++] = a[j];
      index[ki++] = irow[j++]+1;
    }

    // Verifying if A(i,i) exists.

    if ((irow[j] == i)&&(j < end)) { // A(i,i) exists, subtracting sigma.
      value[k++] = a[j++] - sigma;
    }
    else {                           // A(i,i) does not exist.
      value[k++] = -sigma;
    }
    index[ki++] = i+1;

    // Copying subdiagonal elements of column i.

    while (j < end ) {
      value[k++] = a[j];
      index[ki++] = irow[j++]+1;
    }

    index[i+1] = k+1;

  }

} // SubtractAsI.


template<class TYPE>
inline void ARumNonSymMatrix<TYPE>::CreateStructure()
{

  int dimfact = (((fillin+1)*nnz)<(n*n)) ? (fillin+1)*nnz : n*n;

  ClearMem();

  lindex = 30*n+dimfact;          // ?????
  lvalue = dimfact;

  value  = new TYPE[lvalue];
  index  = new int[lindex];

} // CreateStructure.


template<class TYPE>
inline void ARumNonSymMatrix<TYPE>::ThrowError()
{

  if (info[0] < -2)  {       // Memory is not suficient.
    throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                      "ARumNonSymMatrix::FactorA");
  }
  else if (info[0] > 3) {    // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARumNonSymMatrix::FactorA");
  }
  else if (info[0] != 0) {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARumNonSymMatrix::FactorA");
  }

} // ThrowError.


template<class TYPE>
void ARumNonSymMatrix<TYPE>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumNonSymMatrix::FactorA");
  }

  // Quitting the function if A is not square.

  if (m != n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymMatrix::FactorA");
  }

  // Defining local variables.

  int i;
  int *pi, *pj;

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to (value, index);

  copy(nnz, a, 1, value, 1);
  pi=pcol;
  pj=index;
  for (i=0; i<=n; i++) *pj++ = (*pi++)+1;
  pi=irow;
  for (i=0; i<nnz; i++) *pj++ = (*pi++)+1;

  // Decomposing A.

  um2fa(n, nnz, 0, false, lvalue, lindex, value, 
        index, keep, cntl, icntl, info, rinfo);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class TYPE>
void ARumNonSymMatrix<TYPE>::FactorAsI(TYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED,
                      "ARumNonSymMatrix::FactorAsI");
  }

  // Quitting the function if A is not square.

  if (m != n) {
    throw ArpackError(ArpackError::NOT_SQUARE_MATRIX,
                      "ARumNonSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  SubtractAsI(sigma);

  // Decomposing AsI.

  um2fa(n, nnz, 0, false, lvalue, lindex, value,
        index, keep, cntl, icntl, info, rinfo);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class TYPE>
void ARumNonSymMatrix<TYPE>::MultMv(TYPE* v, TYPE* w)
{

  int   i,j;
  TYPE  t;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumNonSymMatrix::MultMv");
  }

  // Determining w = M.v.

  for (i=0; i!=m; i++) w[i]=(TYPE)0;

  for (i=0; i!=n; i++) {
    t = v[i];
    for (j=pcol[i]; j!=pcol[i+1]; j++) {
      w[irow[j]] += t*a[j];
    }
  }

} // MultMv.


template<class TYPE>
void ARumNonSymMatrix<TYPE>::MultMtv(TYPE* v, TYPE* w)
{

  int   i,j;
  TYPE t;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumNonSymMatrix::MultMtv");
  }

  // Determining w = M'.v.

  for (i=0; i!=n; i++) {
    t = (TYPE)0;
    for (j=pcol[i]; j!=pcol[i+1]; j++) {
      t += v[irow[j]]*a[j];
    }
    w[i] = t;
  }

} // MultMtv.


template<class TYPE>
void ARumNonSymMatrix<TYPE>::MultMtMv(TYPE* v, TYPE* w)
{

  TYPE* t = new TYPE[m];

  MultMv(v,t);
  MultMtv(t,w);

  delete[] t;

} // MultMtMv.


template<class TYPE>
void ARumNonSymMatrix<TYPE>::MultMMtv(TYPE* v, TYPE* w)
{

  TYPE* t = new TYPE[n];

  MultMtv(v,t);
  MultMv(t,w);

  delete[] t;

} // MultMMtv.


template<class TYPE>
void ARumNonSymMatrix<TYPE>::Mult0MMt0v(TYPE* v, TYPE* w)
{

  MultMv(&v[m],w);
  MultMtv(v,&w[m]);

} // Mult0MMt0v.


template<class TYPE>
void ARumNonSymMatrix<TYPE>::MultInvv(TYPE* v, TYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARumNonSymMatrix::MultInvv");
  }

  // Solving A.w = v (or AsI.w = v).

  TYPE* space = new TYPE[2*n];

  um2so(n, 0, false, lvalue, lindex, value, index,
        keep, v, w, space, cntl, icntl, info, rinfo);

  delete[] space;

} // MultInvv.


template<class TYPE>
inline void ARumNonSymMatrix<TYPE>::
DefineMatrix(int np, int nnzp, TYPE* ap, int* irowp,
             int* pcolp, double thresholdp, int fillinp,
             bool simest, bool reducible, bool check)
{

  // Defining member variables.

  m         = np;
  n         = np;
  nnz       = nnzp;
  a         = ap;
  irow      = irowp;
  pcol      = pcolp;
  pcol[n]   = nnz;
  fillin    = (fillinp>2) ? fillinp : 2;
  threshold = thresholdp;
  value     = NULL;
  index     = NULL;

  // Preparing umfpack.

  um21i(keep, cntl, icntl, threshold, simest, reducible);

  // Checking data.

  if ((check)&&(!DataOK())) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumNonSymMatrix::DefineMatrix");
  }
  else {
    defined = true;
  }

} // DefineMatrix (square).


template<class TYPE>
inline void ARumNonSymMatrix<TYPE>::
DefineMatrix(int mp, int np, int nnzp, TYPE* ap, int* irowp, int* pcolp)
{

  // Defining member variables.

  m        = mp;
  n        = np;
  nnz      = nnzp;
  a        = ap;
  irow     = irowp;
  pcol     = pcolp;
  pcol[n]  = nnz;
  fillin   = 0;
  defined  = true;

} // DefineMatrix (rectangular).


template<class TYPE>
inline ARumNonSymMatrix<TYPE>::
ARumNonSymMatrix(int np, int nnzp, TYPE* ap, int* irowp,
                 int* pcolp, double thresholdp, int fillinp,
                 bool simest, bool reducible, bool check)   : ARMatrix<TYPE>(np)
{

  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, thresholdp,
               fillinp, simest, reducible, check);

} // Long constructor (square matrix).


template<class TYPE>
inline ARumNonSymMatrix<TYPE>::
ARumNonSymMatrix(int mp, int np, int nnzp, TYPE* ap,
                 int* irowp, int* pcolp)             : ARMatrix<TYPE>(mp, np)
{

  factored = false;
  DefineMatrix(mp, np, nnzp, ap, irowp, pcolp);

} // Long constructor (rectangular matrix).


template<class TYPE>
ARumNonSymMatrix<TYPE>::
ARumNonSymMatrix(char* name, double thresholdp, int fillinp,
                 bool simest, bool reducible, bool check)
{

  factored = false;

  try {
    mat.Define(name);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARumNonSymMatrix");
  }

  if (mat.NCols()==mat.NRows()) {
    DefineMatrix(mat.NCols(), mat.NonZeros(), (TYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), thresholdp,
                 fillinp, simest, reducible, check);
  }
  else {
    DefineMatrix(mat.NRows(), mat.NCols(), mat.NonZeros(),
                 (TYPE*)mat.Entries(), mat.RowInd(), mat.ColPtr());
  }

} // Long constructor (Harwell-Boeing file).


template<class TYPE>
ARumNonSymMatrix<TYPE>& ARumNonSymMatrix<TYPE>::
operator=(const ARumNonSymMatrix<TYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUNSMAT_H
