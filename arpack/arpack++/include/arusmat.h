/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARUSMat.h.
   Arpack++ class ARumSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "aruspen.h"

#ifndef ARUSMAT_H
#define ARUSMAT_H

#include <stddef.h>
#include "arch.h"
#include "armat.h"
#include "arhbmat.h"
#include "arerror.h"
#include "blas1c.h"
#include "umfpackc.h"

template<class TYPE> class ARumSymPencil;

template<class TYPE>
class ARumSymMatrix: public ARMatrix<TYPE> {

  friend class ARumSymPencil<TYPE>;

 protected:

  bool   factored;
  char   uplo;
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

  virtual void Copy(const ARumSymMatrix& other);

  void ExpandA(TYPE sigma = (TYPE)0);

  void CreateStructure();

  void ThrowError();

 public:

  int nzeros() { return nnz; }

  int  FillFact() { return fillin; }

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(TYPE sigma);

  void MultMv(TYPE* v, TYPE* w);

  void MultInvv(TYPE* v, TYPE* w);

  void DefineMatrix(int np, int nnzp, TYPE* ap, int* irowp,
                    int* pcolp, char uplop = 'L', double thresholdp = 0.1, 
                    int fillinp = 9, bool reducible = true, bool check = true);

  ARumSymMatrix(): ARMatrix<TYPE>() { factored = false; }
  // Short constructor that does nothing.

  ARumSymMatrix(int np, int nnzp, TYPE* ap, int* irowp,
                int* pcolp, char uplop = 'L', double thresholdp = 0.1,
                int fillinp = 9, bool reducible = true, bool check = true);
  // Long constructor.

  ARumSymMatrix(char* name, double thresholdp = 0.1, int fillinp = 9,
                bool reducible = true, bool check = true);
  // Long constructor (Harwell-Boeing file).

  ARumSymMatrix(const ARumSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARumSymMatrix() { ClearMem(); }
  // Destructor.

  ARumSymMatrix& operator=(const ARumSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARumSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class TYPE>
bool ARumSymMatrix<TYPE>::DataOK()
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
      if (uplo == 'U') {
        if ((irow[j]<0)||(irow[k]>i)) return false;
      }
      else { // uplo == 'L'.
        if ((irow[j]<i)||(irow[k]>=n)) return false;
      }
      while ((j!=k)&&(irow[j]<irow[j+1])) j++;
      if (j!=k) return false;
    }
  }
  
  return true;

} // DataOK.


template<class TYPE>
inline void ARumSymMatrix<TYPE>::ClearMem()
{

  if (factored) {
    delete[] value;
    delete[] index;
    value = NULL;
    index = NULL;
  }

} // ClearMem.


template<class TYPE>
void ARumSymMatrix<TYPE>::Copy(const ARumSymMatrix<TYPE>& other)
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
  uplo      = other.uplo;

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
void ARumSymMatrix<TYPE>::ExpandA(TYPE sigma)
{

  bool subtract;
  int  i, j, k, ki;

  // Checking if sigma is zero.

  subtract = (sigma != (TYPE)0);

  // Filling index with zeros.

  for (i=0; i<=n; i++) index[i] = 0;

  // Counting the elements in each column of A.

  if (uplo == 'U') {

    for (i=0; i!=n; i++) {
      k = pcol[i+1];
      if ((k!=pcol[i])&&(irow[k-1]==i)) {
        k--;
      }
      else {
        if (subtract) index[i]++;
      }
      for (j=pcol[i]; j<k; j++) index[irow[j]]++;        
    }

  }
  else { // uplo == 'L'

    for (i=0; i!=n; i++) {
      k = pcol[i];
      if ((k!=pcol[i+1])&&(irow[k]==i)) {
        k++;
      }
      else {
        if (subtract) index[i]++;
      }
      for (j=k; j<pcol[i+1]; j++) index[irow[j]]++;        
    }

  }  

  // Summing up index elements.

  for (i=0; i<n; i++) index[i+1]+=index[i];

  // Adding pcol to index.

  for (i=n; i>0; i--) index[i] = index[i-1]+pcol[i];
  index[0] = pcol[0];    

  // Expanding A.

  ki = n+1;

  if (uplo == 'U') {

    for (i=0; i<n; i++) {
      for (j=pcol[i]; j<(pcol[i+1]-1); j++) {
        index[ki+index[i]] = irow[j]+1;
        index[ki+index[irow[j]]] = i+1; 
        value[index[i]++] = a[j];
        value[index[irow[j]]++] = a[j];
      }
      if ((pcol[i]!=pcol[i+1])&&(irow[j]==i)) {
        index[ki+index[i]] = i+1;
        if (subtract) {
          value[index[i]++] = a[j]-sigma;
        }
        else {
          value[index[i]++] = a[j];
        }
      }
      else {
        if (subtract) {
          index[ki+index[i]] = i+1;
          value[index[i]++]  = -sigma;
        }
      }
    }

  }
  else { // uplo  == 'L'

    for (i=0; i<n; i++) {
      k=pcol[i];
      if ((k!=pcol[i+1])&&(irow[k]==i)) {
        index[ki+index[i]] = i+1;
        if (subtract) {
          value[index[i]++] = a[k]-sigma;
        }
        else {
          value[index[i]++] = a[k];
        }
        k++;
      }
      else {
        if (subtract) {
          index[ki+index[i]] = i+1;
          value[index[i]++]  = -sigma;
        }
      }
      for (j=k; j<pcol[i+1]; j++) {
        index[ki+index[i]] = irow[j]+1;
        index[ki+index[irow[j]]] = i+1; 
        value[index[i]++] = a[j];
        value[index[irow[j]]++] = a[j];
      }
    }

  }

  // Adjusting index.

  for (i=n; i>0; i--) {
    index[i] = index[i-1]+1;
  } 
  index[0] = 1;

} // ExpandA.


template<class TYPE>
inline void ARumSymMatrix<TYPE>::CreateStructure()
{

  int dimfact = (((fillin+1)*nnz*2)<(n*n)) ? (fillin+1)*nnz*2 : n*n;

  ClearMem();

  lindex = 30*n+dimfact;          // ?????
  lvalue = dimfact;

  value  = new TYPE[lvalue];
  index  = new int[lindex];

} // CreateStructure.


template<class TYPE>
inline void ARumSymMatrix<TYPE>::ThrowError()
{

  if (info[0] < -2)  {       // Memory is not suficient.
    throw ArpackError(ArpackError::INSUFICIENT_MEMORY,
                      "ARumSymMatrix::FactorA");
  }
  else if (info[0] > 3) {    // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARumSymMatrix::FactorA");
  }
  else if (info[0] != 0) {   // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARumSymMatrix::FactorA");
  }

} // ThrowError.


template<class TYPE>
void ARumSymMatrix<TYPE>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymMatrix::FactorA");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Copying A to (value, index);

  ExpandA();

  // Decomposing A.

  um2fa(n, index[n], 0, false, lvalue, lindex, value, 
        index, keep, cntl, icntl, info, rinfo);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorA.


template<class TYPE>
void ARumSymMatrix<TYPE>::FactorAsI(TYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymMatrix::FactorAsI");
  }

  // Reserving memory for some vectors used in matrix decomposition.

  CreateStructure();

  // Subtracting sigma*I from A.

  ExpandA(sigma);

  // Decomposing AsI.

  um2fa(n, index[n], 0, false, lvalue, lindex, value,
        index, keep, cntl, icntl, info, rinfo);

  // Handling errors.

  ThrowError();

  factored = true;

} // FactorAsI.


template<class TYPE>
void ARumSymMatrix<TYPE>::MultMv(TYPE* v, TYPE* w)
{

  int   i,j,k;
  TYPE  t;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARumSymMatrix::MultMv");
  }

  // Determining w = M.v.

  for (i=0; i!=m; i++) w[i]=(TYPE)0;

  if (uplo == 'U') {

    for (i=0; i!=n; i++) {
      t = v[i];
      k = pcol[i+1];
      if ((k!=pcol[i])&&(irow[k-1]==i)) {
        w[i] += t*a[k-1];
        k--;
      }
      for (j=pcol[i]; j<k; j++) {
        w[irow[j]] += t*a[j];
        w[i] += v[irow[j]]*a[j];
      }
    }

  }
  else {

    for (i=0; i!=n; i++) {
      t = v[i];
      k = pcol[i];
      if ((k!=pcol[i+1])&&(irow[k]==i)) {
        w[i] += t*a[k];
        k++;
      }
      for (j=k; j<pcol[i+1]; j++) {
        w[irow[j]] += t*a[j];
        w[i] += v[irow[j]]*a[j];
      }
    }

  }

} // MultMv.


template<class TYPE>
void ARumSymMatrix<TYPE>::MultInvv(TYPE* v, TYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARumSymMatrix::MultInvv");
  }

  // Solving A.w = v (or AsI.w = v).

  TYPE* space = new TYPE[2*n];

  um2so(n, 0, false, lvalue, lindex, value, index,
        keep, v, w, space, cntl, icntl, info, rinfo);

  delete[] space;

} // MultInvv.


template<class TYPE>
inline void ARumSymMatrix<TYPE>::
DefineMatrix(int np, int nnzp, TYPE* ap, int* irowp,
             int* pcolp, char uplop, double thresholdp,
             int fillinp, bool reducible, bool check)
{

  // Defining member variables.

  m         = np;
  n         = np;
  nnz       = nnzp;
  a         = ap;
  irow      = irowp;
  pcol      = pcolp;
  pcol[n]   = nnz;
  uplo      = uplop;
  fillin    = (fillinp>2) ? fillinp : 2;
  threshold = thresholdp;
  value     = NULL;
  index     = NULL;

  // Preparing umfpack.

  um21i(keep, cntl, icntl, threshold, true, reducible);

  // Checking data.

  if ((check)&&(!DataOK())) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumSymMatrix::DefineMatrix");
  }
  else {
    defined = true;
  }

} // DefineMatrix.


template<class TYPE>
inline ARumSymMatrix<TYPE>::
ARumSymMatrix(int np, int nnzp, TYPE* ap, int* irowp,
              int* pcolp, char uplop, double thresholdp,
              int fillinp, bool reducible, bool check)   : ARMatrix<TYPE>(np)
{

  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, uplop,
               thresholdp, fillinp, reducible, check);

} // Long constructor.


template<class TYPE>
ARumSymMatrix<TYPE>::
ARumSymMatrix(char* file, double thresholdp, int fillinp,
              bool reducible, bool check)
{

  factored = false;

  try {
    mat.Define(file);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARumSymMatrix");
  }

  if ((mat.NCols() == mat.NRows()) && (mat.IsSymmetric())) {

    DefineMatrix(mat.NCols(), mat.NonZeros(), (TYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), 'L', thresholdp,
                 fillinp, reducible, check);
  }
  else {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARumSymMatrix::ARluSymMatrix");
  }

} // Long constructor (Harwell-Boeing file).


template<class TYPE>
ARumSymMatrix<TYPE>& ARumSymMatrix<TYPE>::
operator=(const ARumSymMatrix<TYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSMAT_H
