/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARRSNSym.h.
   Arpack++ class ARrcNonSymStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRSNSYM_H
#define ARRSNSYM_H

#include <stddef.h>
#include "arch.h"
#include "arerror.h"
#include "debug.h"
#include "arrseig.h"
#include "naupp.h"
#include "neupp.h"


template<class FLOAT>
class ARrcNonSymStdEig: public virtual ARrcStdEig<FLOAT, FLOAT> {

 protected:

 // a) Protected functions:

 // a.1) Memory control functions.

  int ValSize() { return nev+1; }
  // Provides the size of array EigVal.

  void ValAllocate();
  // Creates arrays EigValR and EigValI.

  void WorkspaceAllocate();
  // Allocates workspace for nonsymmetric problems.


 // a.2) Functions that handle original FORTRAN ARPACK code.

  void Aupp();
  // Interface to FORTRAN subroutines SNAUPD and DNAUPD.

  void Eupp();
  // Interface to FORTRAN subroutines SNEUPD and DNEUPD.


 // a.3) Functions that check user defined parameters.

  int CheckNev(int nevp);
  // Does Range checking on nev.


 // a.4) Auxiliary functions required when using STL vector class.

  bool ConjEigVec(int i);
  // Indicates if EigVec[i] is the second eigenvector in 
  // a complex conjugate pair.

#ifdef ARCOMP_H
#ifdef VECTOR_H

  vector<arcomplex<FLOAT> >* GenComplex(vector<FLOAT>* RealPart, 
                                        vector<FLOAT>* ImagPart, 
                                        bool conj = false);
  // Generates a complex vector Complex = RealPart + I*ImagPart
  // (or Complex = RealPart - I*ImagPart, if conj = true).

  vector<arcomplex<FLOAT> >* GenComplex(int dim, FLOAT* RealPart, 
                                        FLOAT* ImagPart, bool conj = false);
  // Generates a complex vector Complex = RealPart + I*ImagPart
  // (or Complex = RealPart - I*ImagPart, if conj = true). dim
  // is the length of RealPart and ImagPart.

  vector<arcomplex<FLOAT> >* GenComplex(int dim, FLOAT* RealPart);
  // Generates a complex vector from a real vector. dim is the
  // length of RealPart.

#endif // VECTOR_H.
#endif // ARCOMP_H.

 public:

 // b) Public functions:

 // b.1) Trace functions.

  void Trace(const int digit = -5, const int getv0 = 0, const int aupd = 1,
             const int aup2 = 0,  const int aitr = 0,  const int eigt = 0,
             const int apps = 0,  const int gets = 0,  const int eupd = 0)
  {
    nTraceOn(digit, getv0, aupd, aup2, aitr, eigt, apps, gets, eupd); 
  }
  // Turns on trace mode. 


 // b.2) Functions that permit step by step execution of ARPACK.

  FLOAT* GetVectorImag();
  // When ido = 3, this function indicates where the imaginary part
  // of the eigenvalues of the current Hessenberg matrix are located.


 // b.3) Functions that perform all calculations in one step.

  int Eigenvalues(FLOAT* &EigValRp, FLOAT* &EigValIp,
                   bool ivec = false, bool ischur = false);
  // Overrides arrays EigValRp with the real part and EigValIp 
  // with the imaginary part of the eigenvalues of the problem. 
  // Calculates eigenvectors and Schur vectors if requested.

  int EigenValVectors(FLOAT* &EigVecp, FLOAT* &EigValRp, FLOAT* &EigValIp,
                       bool ischur = false);
  // Overrides array EigVecp sequentially with the eigenvectors of the
  // given eigen-problem. Also stores the eigenvalues in EigValRp and
  // EigValIp. Calculates Schur vectors if requested.


 // b.4) Functions that return elements of vectors and matrices.

#ifdef ARCOMP_H
  arcomplex<FLOAT> Eigenvalue(int i);
  // Furnishes i-eth eigenvalue.
#endif // ARCOMP_H.

  FLOAT EigenvalueReal(int i);
  // Provides the real part of the i-eth eigenvalue.

  FLOAT EigenvalueImag(int i);
  // Provides the imaginary part of the i-eth eigenvalue.

#ifdef ARCOMP_H
  arcomplex<FLOAT> Eigenvector(int i, int j);
  // Furnishes element j of the i-eth eigenvector.
#endif // ARCOMP_H.

  FLOAT EigenvectorReal(int i, int j);
  // Provides the real part of element j of the i-eth eigenvector.

  FLOAT EigenvectorImag(int i, int j);
  // Provides the imaginary part of element j of the i-eth eigenvector.


 // b.5) Functions that provide raw access to internal vectors and matrices.

  FLOAT* RawEigenvaluesImag();
  // Provides raw access to the imaginary part of eigenvalues.


 // b.6) Functions that use STL vector class.

#ifdef VECTOR_H

#ifdef ARCOMP_H
  vector<arcomplex<FLOAT> >* StlEigenvalues(bool ivec = false, 
                                            bool ischur = false);
  // Calculates the eigenvalues and stores them in a single STL vector.
  // Also calculates eigenvectors and Schur vectors if requested.
#endif // ARCOMP_H.

  vector<FLOAT>* StlEigenvaluesReal();
  // Returns the real part of the eigenvalues.

  vector<FLOAT>* StlEigenvaluesImag();
  // Returns the imaginary part of the eigenvalues.

#ifdef ARCOMP_H
  vector<arcomplex<FLOAT> >* StlEigenvector(int i);
  // Returns the i-th eigenvector.
#endif // ARCOMP_H.

  vector<FLOAT>* StlEigenvectorReal(int i);
  // Returns the real part of the i-th eigenvector.

  vector<FLOAT>* StlEigenvectorImag(int i);
  // Returns the imaginary part of the i-th eigenvector.

#endif // VECTOR_H.


 // b.7) Constructors and destructor.

  ARrcNonSymStdEig() { }
  // Short constructor.

  ARrcNonSymStdEig(int np, int nevp, char* whichp = "LM", int ncvp = 0,
                   FLOAT tolp = 0.0, int maxitp = 0, FLOAT* residp = NULL,
                   bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcNonSymStdEig(int np, int nevp, FLOAT sigma, char* whichp = "LM",
                   int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                   FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARrcNonSymStdEig(const ARrcNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcNonSymStdEig() { }
  // Destructor.

 // c) Operators.

  ARrcNonSymStdEig& operator=(const ARrcNonSymStdEig& other);
  // Assignment operator.

}; // class ARrcNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARrcNonSymStdEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline void ARrcNonSymStdEig<FLOAT>::ValAllocate()
{

  if (EigValR == NULL) {
    EigValR = new FLOAT[ValSize()];
    EigValI = new FLOAT[ValSize()];
    newVal = true;
  }

} // ValAllocate.


template<class FLOAT>
inline void ARrcNonSymStdEig<FLOAT>::WorkspaceAllocate()
{

  lworkl  = 3*ncv*(ncv+2);
  lworkv  = 3*ncv;
  lrwork  = 0;
  workl   = new FLOAT[lworkl+1];
  workv   = new FLOAT[lworkv+1];

} // WorkspaceAllocate.


template<class FLOAT>
inline void ARrcNonSymStdEig<FLOAT>::Aupp()
{

  naupp(ido, bmat, n, which, nev, tol, resid, ncv, V, n,
        iparam, ipntr, workd, workl, lworkl, info);

} // Aupp.


template<class FLOAT>
inline void ARrcNonSymStdEig<FLOAT>::Eupp()
{

  neupp(rvec, HowMny, EigValR, EigValI, EigVec, n, sigmaR,
        sigmaI, workv, bmat, n, which, nev, tol, resid, ncv, V,
        n, iparam, ipntr, workd, workl, lworkl, info);

} // Eupp.


template<class FLOAT>
inline int ARrcNonSymStdEig<FLOAT>::CheckNev(int nevp)
{

  if ((nevp<=1)||(nevp>=(n-1))) { // nev must satisfy 1 < nev < n-1.
    throw ArpackError(ArpackError::NEV_OUT_OF_BOUNDS);
  }
  return nevp;

} // CheckNev.


template<class FLOAT>
bool ARrcNonSymStdEig<FLOAT>::ConjEigVec(int i)
{

  if (EigValI[i] == (FLOAT)0.0) return false;
  int j = i-1;
  while ((j >= 0) && (EigValI[j] != (FLOAT)0.0)) j--;
  if (((i-j)%2) == 0) {
    return true;
  }
  else {
    return false;
  }

} // ConjEigVec.


#ifdef VECTOR_H // Defining functions that use STL vector class.
#ifdef ARCOMP_H

template<class FLOAT>
vector<arcomplex<FLOAT> >* ARrcNonSymStdEig<FLOAT>::
GenComplex(vector<FLOAT>* RealPart, vector<FLOAT>* ImagPart, bool conj)
{

  // Defining variables.

  vector<arcomplex<FLOAT> >* Result;
  try {
    Result = new vector<arcomplex<FLOAT> >(ValSize());
  }
  catch (ArpackError) { return NULL; }
  FLOAT* rp  = RealPart->begin();
  FLOAT* ip  = ImagPart->begin();
  FLOAT* end = RealPart->end();
  arcomplex<FLOAT>* s = Result->begin();

  // Creating a complex vector.

  if (!conj) {
    while (rp != end) *s++ = arcomplex<FLOAT>(*rp++, *ip++);
  }
  else {
    while (rp != end) *s++ = arcomplex<FLOAT>(*rp++, -(*ip++));
  }

  return Result;

} // GenComplex (vector<FLOAT> version).


template<class FLOAT>
vector<arcomplex<FLOAT> >* ARrcNonSymStdEig<FLOAT>::
GenComplex(int dim, FLOAT* RealPart, FLOAT* ImagPart, bool conj)
{

  // Defining variables.

  vector<arcomplex<FLOAT> >* Result;
  try {
    Result = new vector<arcomplex<FLOAT> >(dim);
  }
  catch (ArpackError) { return NULL; }
  FLOAT* rp  = RealPart; 
  FLOAT* ip  = ImagPart; 
  FLOAT* end = &RealPart[dim];
  arcomplex<FLOAT>* s = Result->begin();

  // Creating a complex vector.

  if (!conj) {
    while (rp != end) *s++ = arcomplex<FLOAT>(*rp++, *ip++);
  }
  else {
    while (rp != end) *s++ = arcomplex<FLOAT>(*rp++, -(*ip++));
  }

  return Result;

} // GenComplex (FLOAT* version).


template<class FLOAT>
vector<arcomplex<FLOAT> >* ARrcNonSymStdEig<FLOAT>::
GenComplex(int dim, FLOAT* RealPart)
{

  // Defining variables.

  vector<arcomplex<FLOAT> >* Result;
  try {
    Result = new vector<arcomplex<FLOAT> >(dim);
  }
  catch (ArpackError) { return NULL; }
  FLOAT* rp  = RealPart; 
  FLOAT* end = &RealPart[dim];
  arcomplex<FLOAT>* s = Result->begin();

  // Copying a real vector into a complex vector.

  while (rp != end) *s++ = *rp++;

  return Result;

} // GenComplex.

#endif // ARCOMP_H.
#endif // VECTOR_H.


template<class FLOAT>
FLOAT* ARrcNonSymStdEig<FLOAT>::GetVectorImag()
{

  if (ido != 3) {
    throw ArpackError(ArpackError::CANNOT_GET_VECTOR, "GetVectorImag");
  }
  return &workl[ipntr[6]];

} // GetVectorImag.


template<class FLOAT>
int ARrcNonSymStdEig<FLOAT>::
Eigenvalues(FLOAT* &EigValRp, FLOAT* &EigValIp, bool ivec, bool ischur)
{

  if (ValuesOK) {                                 // Eigenvalues are available.
    if ((EigValRp == NULL)&&(EigValIp == NULL)) { // Moving eigenvalues.
      EigValRp = EigValR;
      EigValIp = EigValI;
      EigValR  = NULL;
      EigValI  = NULL;
      newVal   = false;
      ValuesOK = false;
    }
    else {                                        // Copying eigenvalues.
      try {
        if (EigValRp == NULL) EigValRp = new FLOAT[ValSize()];
        if (EigValIp == NULL) EigValIp = new FLOAT[ValSize()];
      }
      catch (ArpackError) { return 0; }
      copy(nconv,EigValR,1,EigValRp,1);
      copy(nconv,EigValI,1,EigValIp,1);
    }
  }
  else {
    if (newVal) {
      delete[] EigValR;
      delete[] EigValI;
      newVal = false;
    }
    try {
      if (EigValRp == NULL) EigValRp = new FLOAT[ValSize()];
      if (EigValIp == NULL) EigValIp = new FLOAT[ValSize()];
    }
    catch (ArpackError) { return 0; }
    EigValR = EigValRp;
    EigValI = EigValIp;
    if (ivec) {                              // Finding eigenvalues and vectors.
      nconv = FindEigenvectors(ischur);
    }
    else {                                   // Finding eigenvalues only.
      nconv = FindEigenvalues();
    }
    EigValR = NULL;
    EigValI = NULL;
  }
  return nconv;

} // Eigenvalues(EigValRp, EigValIp, ivec, ischur).


template<class FLOAT>
int ARrcNonSymStdEig<FLOAT>::
EigenValVectors(FLOAT* &EigVecp, FLOAT* &EigValRp,
                FLOAT* &EigValIp, bool ischur)
{

  if (ValuesOK) {               // Eigenvalues are already available .
    nconv = Eigenvalues(EigValRp, EigValIp, false);
    nconv = Eigenvectors(EigVecp, ischur);
  }
  else {                        // Eigenvalues ans vectors are not available.
    if (newVec) {
      delete[] EigVec;
      newVec = false;
    }
    if (newVal) {
      delete[] EigValR;
      delete[] EigValI;
      newVal = false;
    }
    try {
      if (EigVecp  == NULL) EigVecp  = new FLOAT[ValSize()*n];
      if (EigValRp == NULL) EigValRp = new FLOAT[ValSize()];
      if (EigValIp == NULL) EigValIp = new FLOAT[ValSize()];
    }
    catch (ArpackError) { return 0; }
    EigVec  = EigVecp;
    EigValR = EigValRp;
    EigValI = EigValIp;
    nconv   = FindEigenvectors(ischur);
    EigVec  = NULL;
    EigValR = NULL;
    EigValI = NULL;
  }
  return nconv;

} // EigenValVectors(EigVecp, EigValRp, EigValIp, ischur).


#ifdef ARCOMP_H
template<class FLOAT>
inline arcomplex<FLOAT> ARrcNonSymStdEig<FLOAT>::Eigenvalue(int i)
{

  // Returning i-eth eigenvalue.

  if (!ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "Eigenvalue(i)");
  }
  else if ((i>=nconv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvalue(i)");
  }
  return arcomplex<FLOAT>(EigValR[i],EigValI[i]);

} // Eigenvalue(i).
#endif // ARCOMP_H


template<class FLOAT>
inline FLOAT ARrcNonSymStdEig<FLOAT>::EigenvalueReal(int i)
{

  // Returning the real part of i-eth eigenvalue.

  if (!ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "EigenvalueReal(i)");
  }
  else if ((i>=nconv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "EigenvalueReal(i)");
  }
  return EigValR[i];

} // EigenvalueReal(i).


template<class FLOAT>
inline FLOAT ARrcNonSymStdEig<FLOAT>::EigenvalueImag(int i)
{

  // Returning the imaginary part of i-eth eigenvalue.

  if (!ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "EigenvalueImag(i)");
  }
  else if ((i>=nconv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "EigenvalueImag(i)");
  }
  return EigValI[i];

} // EigenvalueImag(i).


#ifdef ARCOMP_H
template<class FLOAT>
inline arcomplex<FLOAT> ARrcNonSymStdEig<FLOAT>::
Eigenvector(int i, int j)
{

  // Returning element j of i-eth eigenvector.

  if ((!VectorsOK)||(!ValuesOK)) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "Eigenvector(i,j)");
  }
  else if ((i>=nconv)||(i<0)||(j>=n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvector(i,j)");
  }
  if (EigValI[i]==(FLOAT)0.0) {   // Real eigenvalue.
    return arcomplex<FLOAT>(EigVec[i*n+j],(FLOAT)0.0);
  }
  else {                          // Complex eigenvalue.
    if (EigValI[i]>(FLOAT)0.0) {  // with positive imaginary part.
      return arcomplex<FLOAT>(EigVec[i*n+j], EigVec[(i+1)*n+j]);
    }
    else {                        // with negative imaginary part.
      return arcomplex<FLOAT>(EigVec[(i-1)*n+j], -EigVec[i*n+j]);
    }
  }

} // Eigenvector(i,j).
#endif // ARCOMP_H


template<class FLOAT>
inline FLOAT ARrcNonSymStdEig<FLOAT>::EigenvectorReal(int i, int j)
{

  // Returning the real part of element j of i-eth eigenvector.

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "EigenvectorReal(i,j)");
  }
  else if ((i>=nconv)||(i<0)||(j>=n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "EigenvectorReal(i,j)");
  }
  return EigVec[i*n+j];

} // EigenvectorReal(i,j).


template<class FLOAT>
inline FLOAT ARrcNonSymStdEig<FLOAT>::EigenvectorImag(int i, int j)
{

  // Returning the imaginary part of element j of i-eth eigenvector.

  if ((!VectorsOK)||(!ValuesOK)) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "EigenvectorImag(i,j)");
  }
  else if ((i>=nconv)||(i<0)||(j>=n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "EigenvectorImag(i,j)");
  }
  if (EigValI[i]==(FLOAT)0.0) {   // Real eigenvalue.
    return (FLOAT)0.0;
  }
  else {                          // Complex eigenvalue.
    if (EigValI[i]>(FLOAT)0.0) {  // with positive imaginary part.
      return EigVec[(i+1)*n+j];
    }
    else {                        // with negative imaginary part.
      return -EigVec[i*n+j];
    }
  }

} // EigenvectorImag(i,j).


template<class FLOAT>
inline FLOAT* ARrcNonSymStdEig<FLOAT>::RawEigenvaluesImag()
{

  if (!ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "RawEigenvaluesImag");
  }
  return EigValI;

} // RawEigenvaluesImag.


#ifdef VECTOR_H // Defining some functions that use STL vector class.

#ifdef ARCOMP_H
template<class FLOAT>
inline vector<arcomplex<FLOAT> >* ARrcNonSymStdEig<FLOAT>::
StlEigenvalues(bool ivec, bool ischur)
{

  // Returning the eigenvalues in a STL vector.

  // Defining temporary variables.

  vector<FLOAT>* StlEigValR;
  vector<FLOAT>* StlEigValI;
  FLOAT*         ValRPtr;
  FLOAT*         ValIPtr;

  try {
    StlEigValR = new vector<FLOAT>(ValSize());
    StlEigValI = new vector<FLOAT>(ValSize());
  }
  catch (ArpackError) { return NULL; }

  // Finding Eigenvalues.

  ValRPtr = StlEigValR->begin();
  ValIPtr = StlEigValI->begin();
  nconv = Eigenvalues(ValRPtr, ValIPtr, ivec, ischur);
  vector<arcomplex<FLOAT> >* Val = GenComplex(StlEigValR, StlEigValI);

  // Deleting temporary variables.

  delete StlEigValR;
  delete StlEigValI;

  return Val;

} // StlEigenvalues.
#endif // ARCOMP_H.


template<class FLOAT>
inline vector<FLOAT>* ARrcNonSymStdEig<FLOAT>::StlEigenvaluesReal()
{

  // Returning the real part of the eigenvalues in a STL vector.

  vector<FLOAT>* StlEigValR;
  
  if (!ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "StlEigenvaluesReal");
  }
  try {
    StlEigValR = new vector<FLOAT>(EigValR, &EigValR[ValSize()]);
  }
  catch (ArpackError) { return NULL; }
  return StlEigValR;

} // StlEigenvaluesReal.


template<class FLOAT>
inline vector<FLOAT>* ARrcNonSymStdEig<FLOAT>::StlEigenvaluesImag()
{

  // Returning the imaginary part of the eigenvalues in a STL vector.

  vector<FLOAT>* StlEigValI;

  if (!ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "StlEigenvaluesImag");
  }
  try {
    StlEigValI = new vector<FLOAT>(EigValI, &EigValI[ValSize()]);
  }
  catch (ArpackError) { return NULL; }
  return StlEigValI;

} // StlEigenvaluesImag.


#ifdef ARCOMP_H
template<class FLOAT>
inline vector<arcomplex<FLOAT> >* ARrcNonSymStdEig<FLOAT>::
StlEigenvector(int i)
{

  // Returning the i-th eigenvector in a STL vector.

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvector(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvector(i)");
  }
  if (EigValI[i] == (FLOAT)0.0) { // Real eigenvector.
    return GenComplex(n, &EigVec[i*n]);
  }
  else if (!ConjEigVec(i)) {      // First eigenvector in a conjugate pair.
    return GenComplex(n, &EigVec[i*n], &EigVec[(i+1)*n]);
  }
  else {                          // Second eigenvector in a conjugate pair.
    return GenComplex(n, &EigVec[(i-1)*n], &EigVec[i*n], true);
  }

} // StlEigenvector(i).
#endif // ARCOMP_H.


template<class FLOAT>
inline vector<FLOAT>* ARrcNonSymStdEig<FLOAT>::StlEigenvectorReal(int i)
{

  // Returning the real part of the i-th eigenvector in a STL vector.

  vector<FLOAT>* Vec;

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvectorReal(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvectorReal(i)");
  }
  if (!ConjEigVec(i)) { // Real eigenvector or first in a conj. pair.
    try {
      Vec = new vector<FLOAT>(&EigVec[i*n], &EigVec[(i+1)*n]);
    }
    catch (ArpackError) { return NULL; }
    return Vec;
  }
  else {                // Second eigenvector in a conjugate pair.
    try {
      Vec = new vector<FLOAT>(&EigVec[(i-1)*n], &EigVec[i*n]);
    }
    catch (ArpackError) { return NULL; }
    return Vec;
  }

} // StlEigenvectorReal(i).


template<class FLOAT>
inline vector<FLOAT>* ARrcNonSymStdEig<FLOAT>::StlEigenvectorImag(int i)
{

  // Returning the imaginary part of the i-th eigenvector in a STL vector.

  vector<FLOAT>* Vec;

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvectorImag(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvectorImag(i)");
  }
  if (EigValI[i] == (FLOAT)0.0) { // Real eigenvector.
    try {
      Vec = new vector<FLOAT>(ValSize(), (FLOAT)0.0);
    }
    catch (ArpackError) { return NULL; }
    return Vec;
  }
  else if (!ConjEigVec(i)) {      // First eigenvector in a conjugate pair.
    try {
      Vec = new vector<FLOAT>(&EigVec[(i+1)*n], &EigVec[(i+2)*n]);
    }
    catch (ArpackError) { return NULL; }
    return Vec;
  }
  else {                          // Second eigenvector in a conjugate pair.
    try {
      Vec = new vector<FLOAT>(&EigVec[i*n], &EigVec[(i+1)*n]);
    }
    catch (ArpackError) { return NULL; }
    for (FLOAT* s = Vec->begin(); s != Vec->end(); s++) *s = -(*s);
    return Vec;
  }

} // StlEigenvectorImag(i).

#endif // VECTOR_H.


template<class FLOAT>
inline ARrcNonSymStdEig<FLOAT>::
ARrcNonSymStdEig(int np, int nevp, char* whichp, int ncvp,
                 FLOAT tolp, int maxitp, FLOAT* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARrcNonSymStdEig<FLOAT>::
ARrcNonSymStdEig(int np, int nevp, FLOAT sigmap, char* whichp, int ncvp,
                 FLOAT tolp, int maxitp, FLOAT* residp, bool ishiftp)

{

  ChangeShift(sigmap);
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class FLOAT>
ARrcNonSymStdEig<FLOAT>& ARrcNonSymStdEig<FLOAT>::
operator=(const ARrcNonSymStdEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRSNSYM_H

