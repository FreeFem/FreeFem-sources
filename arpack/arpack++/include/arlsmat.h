/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARLSMat.h.
   Arpack++ class ARluSymMatrix definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/


#include "arlspen.h"

#ifndef ARLSMAT_H
#define ARLSMAT_H

#include <stddef.h>
#include "arch.h"
#include "armat.h"
#include "arhbmat.h"
#include "arerror.h"
#include "blas1c.h"
#include "superluc.h"
#include "arlspdef.h"
#include "arlutil.h"

template<class TYPE> class ARluSymPencil;

template<class TYPE>
class ARluSymMatrix: public ARMatrix<TYPE> {

  friend class ARluSymPencil<TYPE>;

 protected:

  bool        factored;
  char        uplo;
  int         order;
  int         nnz;
  int*        irow;
  int*        pcol;
  int*        permc;
  int*        permr;
  double      threshold;
  TYPE*       a;
  SuperMatrix A;
  SuperMatrix L;
  SuperMatrix U;
  ARhbMatrix<int, TYPE> mat;

  bool DataOK();

  virtual void Copy(const ARluSymMatrix& other);

  void ClearMem();

  void ExpandA(NCformat& A, NCformat& Aexp, TYPE sigma = (TYPE)0);

 public:

  int nzeros() { return nnz; }

  bool IsFactored() { return factored; }

  void FactorA();

  void FactorAsI(TYPE sigma);

  void MultMv(TYPE* v, TYPE* w);

  void MultInvv(TYPE* v, TYPE* w);

  void DefineMatrix(int np, int nnzp, TYPE* ap, int* irowp, int* pcolp,
                    char uplop = 'L', double thresholdp = 0.1,
                    int orderp = 2, bool check = true);

  ARluSymMatrix();
  // Short constructor that does nothing.

  ARluSymMatrix(int np, int nnzp, TYPE* ap, int* irowp, int* pcolp,
                char uplop = 'L', double thresholdp = 0.1,
                int orderp = 2, bool check = true);
  // Long constructor.

  ARluSymMatrix(char* name, double thresholdp = 0.1,
                int orderp = 2, bool check = true);
  // Long constructor (Harwell-Boeing file).

  ARluSymMatrix(const ARluSymMatrix& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluSymMatrix() { ClearMem(); }
  // Destructor.

  ARluSymMatrix& operator=(const ARluSymMatrix& other);
  // Assignment operator.

};

// ------------------------------------------------------------------------ //
// ARluSymMatrix member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class TYPE>
bool ARluSymMatrix<TYPE>::DataOK()
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
inline void ARluSymMatrix<TYPE>::Copy(const ARluSymMatrix<TYPE>& other)
{

  // Copying very fundamental variables.

  defined   = other.defined;
  factored  = other.factored;

  // Returning from here if "other" was not initialized.

  if (!defined) return;

  // Copying user-defined parameters.

  DefineMatrix(other.n, other.nnz, other.a, other.irow, other.pcol,
               other.uplo, other.threshold, other.order);

  // Throwing the original factorization away (this procedure 
  // is really awkward, but it is necessary because there
  // is no copy function for matrices L and U in the SuperLU 
  // library and it is not a good idea to do this kind of deep 
  // copy here).

  if (factored) {
    ArpackError(ArpackError::LAPACK_ERROR, "ARluSymMatrix");
    factored = false;
  }

} // Copy.


template<class TYPE>
void ARluSymMatrix<TYPE>::ClearMem()
{

  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree();
  }
  if (defined) {
    delete A.Store;
    delete[] permc;
    delete[] permr;
    permc = NULL;
    permr = NULL;
    A.Store = NULL;
  }

} // ClearMem.


template<class TYPE>
void ARluSymMatrix<TYPE>::
ExpandA(NCformat& A, NCformat& Aexp, TYPE sigma)
{

  // Defining local variables.

  bool  subtract;
  int   i, j, k;
  int   *colA, *colE;
  int   *indA, *indE;
  TYPE  *valA, *valE;

  // Checking if sigma is zero.

  subtract = (sigma != (TYPE)0);

  // Simplifying the notation.

  valA = (TYPE*)A.nzval;
  valE = (TYPE*)Aexp.nzval;
  indA = (int*)A.rowind;
  indE = (int*)Aexp.rowind;
  colA = (int*)A.colptr;
  colE = (int*)Aexp.colptr;

  // Filling colE with zeros.

  for (i=0; i<=n; i++) colE[i] = 0;

  // Counting the elements in each column of A.

  if (uplo == 'U') {

    for (i=0; i!=n; i++) {
      k = colA[i+1];
      if ((k!=colA[i])&&(indA[k-1]==i)) {
        k--;
      }
      else {
        if (subtract) colE[i]++;
      }
      for (j=colA[i]; j<k; j++) colE[indA[j]]++;        
    }

  }
  else { // uplo == 'L'

    for (i=0; i!=n; i++) {
      k = colA[i];
      if ((k!=colA[i+1])&&(indA[k]==i)) {
        k++;
      }
      else {
        if (subtract) colE[i]++;
      }
      for (j=k; j<colA[i+1]; j++) colE[indA[j]]++;        
    }

  }  

  // Summing up colE elements.

  for (i=0; i<n; i++) colE[i+1]+=colE[i];

  // Adding colA to colE.

  for (i=n; i>0; i--) colE[i] = colE[i-1]+colA[i];
  colE[0] = colA[0];    

  // Expanding A.

  if (uplo == 'U') {

    for (i=0; i<n; i++) {
      for (j=colA[i]; j<(colA[i+1]-1); j++) {
        indE[colE[i]] = indA[j];
        indE[colE[indA[j]]] = i; 
        valE[colE[i]++] = valA[j];
        valE[colE[indA[j]]++] = valA[j];
      }
      if ((colA[i]!=colA[i+1])&&(indA[j]==i)) {
        indE[colE[i]] = i;
        if (subtract) {
          valE[colE[i]++] = valA[j]-sigma;
        }
        else {
          valE[colE[i]++] = valA[j];
        }
      }
      else {
        if (subtract) {
          indE[colE[i]] = i;
          valE[colE[i]++] = -sigma;
        }
      }
    }

  }
  else { // uplo  == 'L'

    for (i=0; i<n; i++) {
      k=colA[i];
      if ((k!=colA[i+1])&&(indA[k]==i)) {
        indE[colE[i]] = i;
        if (subtract) {
          valE[colE[i]++] = valA[k]-sigma;
        }
        else {
          valE[colE[i]++] = valA[k];
        }
        k++;
      }
      else {
        if (subtract) {
          indE[colE[i]] = i;
          valE[colE[i]++]  = -sigma;
        }
      }
      for (j=k; j<colA[i+1]; j++) {
        indE[colE[i]] = indA[j];
        indE[colE[indA[j]]] = i; 
        valE[colE[i]++] = valA[j];
        valE[colE[indA[j]]++] = valA[j];
      }
    }

  }

  // Adjusting index.

  for (i=n; i>0; i--) {
    colE[i] = colE[i-1];
  } 
  colE[0] = 0;

  Aexp.nnz = colE[n];

} // ExpandA.


template<class TYPE>
void ARluSymMatrix<TYPE>::FactorA()
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARluSymMatrix::FactorA");
  }

  // Defining local variables.

  int         info;
  int*        etree;
  int*        irowi;
  int*        pcoli;
  TYPE*       aexp;
  SuperMatrix Aexp;
  SuperMatrix AC;
  NCformat*   Astore;
  NCformat*   Aexpstore;

  // Deleting previous versions of L and U.
  
  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree();
  }

  // Setting default values for gstrf parameters.

  double drop_tol   = 0.0;
  int    panel_size = sp_ienv(1);
  int    relax      = sp_ienv(2);

  // Creating a temporary matrix Aexp.

  irowi = new int[nnz*2];
  pcoli = new int[n+1];
  aexp  = new TYPE[nnz*2];
  Create_CompCol_Matrix(&Aexp, n,  n, nnz, aexp, irowi, pcoli, NC, GE);

  // Expanding A.

  Astore    = (NCformat*)A.Store;
  Aexpstore = (NCformat*)Aexp.Store;
  ExpandA(*Astore, *Aexpstore);

  // Reserving memory for etree (used in matrix decomposition).

  etree = new int[n];

  // Defining LUStat.

  extern SuperLUStat_t SuperLUStat;
  StatInit(panel_size, relax);

  // Defining the column permutation of matrix A
  // (using minimum degree ordering).

  get_perm_c(order, &Aexp, permc);

  // Permuting columns of A and creating the elimination tree.

  sp_preorder("N", &Aexp, permc, etree, &AC);

  // Decomposing A.

  gstrf("N",&AC, threshold, drop_tol, relax, panel_size, etree,
        NULL, 0, permr, permc, &L, &U, &info);

  // Deleting AC, Aexp and etree.

  Destroy_CompCol_Permuted(&AC);
  Destroy_CompCol_Matrix(&Aexp);
  delete[] etree;

  factored = (info == 0);

  // Handling errors.

  if (info < 0)  {        // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARluSymMatrix::FactorA");
  }
  else if (info > n) {    // Memory is not sufficient.
    throw ArpackError(ArpackError::MEMORY_OVERFLOW,
                      "ARluSymMatrix::FactorA");
  }
  else if (info > 0) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARluSymMatrix::FactorA");
  }

} // FactorA.


template<class TYPE>
void ARluSymMatrix<TYPE>::FactorAsI(TYPE sigma)
{

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARluSymMatrix::FactorAsI");
  }

  // Defining local variables.

  int         info;
  int*        etree;
  int*        irowi;
  int*        pcoli;
  TYPE*       asi;
  SuperMatrix AsI;
  SuperMatrix AC;
  NCformat*   Astore;
  NCformat*   AsIstore;

  // Deleting previous versions of L and U.
  
  if (factored) {
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree();
  }

  // Setting default values for gstrf parameters.

  double drop_tol   = 0.0;
  int    panel_size = sp_ienv(1);
  int    relax      = sp_ienv(2);

  // Creating a temporary matrix AsI.

  irowi = new int[nnz*2+n];
  pcoli = new int[n+1];
  asi   = new TYPE[nnz*2+n];
  Create_CompCol_Matrix(&AsI, n,  n, nnz, asi, irowi, pcoli, NC, GE);

  // Subtracting sigma*I from A and storing the result on AsI.

  Astore   = (NCformat*)A.Store;
  AsIstore = (NCformat*)AsI.Store;
  ExpandA(*Astore, *AsIstore, sigma);

  // Reserving memory for etree (used in matrix decomposition).

  etree = new int[n];

  // Defining LUStat.

  extern SuperLUStat_t SuperLUStat;
  StatInit(panel_size, relax);

  // Defining the column permutation of matrix AsI
  // (using minimum degree ordering).

  get_perm_c(order, &AsI, permc);

  // Permuting columns of AsI and creating the elimination tree.

  sp_preorder("N", &AsI, permc, etree, &AC);

  // Decomposing AsI.

  gstrf("N",&AC, threshold, drop_tol, relax, panel_size, etree,
        NULL, 0, permr, permc, &L, &U, &info);

  // Deleting AC, AsI and etree.

  Destroy_CompCol_Permuted(&AC);
  Destroy_CompCol_Matrix(&AsI);
  delete[] etree;

  factored = (info == 0);

  // Handling errors.

  if (info < 0)  {        // Illegal argument.
    throw ArpackError(ArpackError::PARAMETER_ERROR,
                      "ARluSymMatrix::FactorAsI");
  }
  else if (info > n) {    // Memory is not sufficient.
    throw ArpackError(ArpackError::MEMORY_OVERFLOW,
                      "ARluSymMatrix::FactorAsI");
  }
  else if (info > 0) {   // Matrix is singular.
    throw ArpackError(ArpackError::MATRIX_IS_SINGULAR,
                      "ARluSymMatrix::FactorAsI");
  }

} // FactorAsI.


template<class TYPE>
void ARluSymMatrix<TYPE>::MultMv(TYPE* v, TYPE* w)
{

  int  i,j, k;
  TYPE t;

  // Quitting the function if A was not defined.

  if (!IsDefined()) {
    throw ArpackError(ArpackError::DATA_UNDEFINED, "ARluSymMatrix::MultMv");
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
void ARluSymMatrix<TYPE>::MultInvv(TYPE* v, TYPE* w)
{

  // Quitting the function if A (or AsI) was not factored.

  if (!IsFactored()) {
    throw ArpackError(ArpackError::NOT_FACTORED_MATRIX,
                      "ARluSymMatrix::MultInvv");
  }

  // Solving A.w = v (or AsI.w = v).

  int         info;
  SuperMatrix B;

  if (&v != &w) copy(n, v, 1, w, 1);
  Create_Dense_Matrix(&B, n, 1, w, n, DN, GE);
  gstrs("N", &L, &U, permr, permc, &B, &info);
  delete B.Store;

} // MultInvv.


template<class TYPE>
inline void ARluSymMatrix<TYPE>::
DefineMatrix(int np, int nnzp, TYPE* ap, int* irowp, int* pcolp,
             char uplop, double thresholdp, int orderp, bool check)
{

  m         = np;
  n         = np;
  nnz       = nnzp;
  a         = ap;
  irow      = irowp;
  pcol      = pcolp;
  pcol[n]   = nnz;
  uplo      = uplop;
  threshold = thresholdp;
  order     = orderp;

  // Checking data.

  if ((check)&&(!DataOK())) {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARluSymMatrix::DefineMatrix");
  }

  // Creating SuperMatrix A.

  Create_CompCol_Matrix(&A, n, n, nnz, a, irow, pcol, NC, GE);

  // Reserving memory for vectors used in matrix decomposition.

  permc = new int[n];
  permr = new int[n];

  defined = true;

} // DefineMatrix.


template<class TYPE>
inline ARluSymMatrix<TYPE>::ARluSymMatrix(): ARMatrix<TYPE>()
{

  factored = false;
  permc    = NULL;
  permr    = NULL;
 
} // Short constructor.


template<class TYPE>
inline ARluSymMatrix<TYPE>::
ARluSymMatrix(int np, int nnzp, TYPE* ap, int* irowp,
              int* pcolp, char uplop, double thresholdp,
              int orderp, bool check)                   : ARMatrix<TYPE>(np)
{

  factored = false;
  DefineMatrix(np, nnzp, ap, irowp, pcolp, uplop, thresholdp, orderp, check);

} // Long constructor.


template<class TYPE>
ARluSymMatrix<TYPE>::
ARluSymMatrix(char* file, double thresholdp, int orderp, bool check)
{

  factored = false;

  try {
    mat.Define(file);
  }
  catch (ArpackError) {    // Returning from here if an error has occurred.
    throw ArpackError(ArpackError::CANNOT_READ_FILE, "ARluSymMatrix");
  }

  if ((mat.NCols() == mat.NRows()) && (mat.IsSymmetric())) {

    DefineMatrix(mat.NCols(), mat.NonZeros(), (TYPE*)mat.Entries(),
                 mat.RowInd(), mat.ColPtr(), 'L', thresholdp, orderp, check);
  }
  else {
    throw ArpackError(ArpackError::INCONSISTENT_DATA,
                      "ARluSymMatrix::ARluSymMatrix");
  }

} // Long constructor (Harwell-Boeing file).


template<class TYPE>
ARluSymMatrix<TYPE>& ARluSymMatrix<TYPE>::
operator=(const ARluSymMatrix<TYPE>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLSMAT_H
