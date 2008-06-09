/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARRSSym.h.
   Arpack++ class ARrcSymStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRSSYM_H
#define ARRSSYM_H

#include <stddef.h>
#include "arch.h"
#include "arerror.h"
#include "debug.h"
#include "arrseig.h"
#include "saupp.h"
#include "seupp.h"


template<class FLOAT>
class ARrcSymStdEig: public virtual ARrcStdEig<FLOAT, FLOAT> {

 protected:

 // a) Protected functions:

 // a.1) Memory control functions.

  void WorkspaceAllocate();
  // Allocates workspace for symmetric problems.


 // a.2) Functions that handle original FORTRAN ARPACK code.

  void Aupp();
  // Interface to FORTRAN subroutines SSAUPD and DSAUPD.

  void Eupp();
  // Interface to FORTRAN subroutines SSEUPD and DSEUPD.


 // a.3) Functions that check user defined parameters.

  const char* CheckWhich(const char* whichp);
  // Determines if the value of variable "this->which" is valid.


 public:

 // b) Public functions:

 // b.1) Trace functions.

  void Trace(const int digit = -5, const int getv0 = 0, const int aupd = 1,
             const int aup2 = 0,  const int aitr = 0,  const int eigt = 0,
             const int apps = 0,  const int gets = 0,  const int eupd = 0) 
  { 
    sTraceOn(digit, getv0, aupd, aup2, aitr, eigt, apps, gets, eupd);
  }
  // Turns on trace mode. 


 // b.2) Functions that permit step by step execution of ARPACK.

  FLOAT* PutVector();
  // When this->ido = -1, 1 or 2 and the user must perform a product in the form
  // y <- M*x, this function indicates where to store y. When this->ido = 3, this
  // function indicates where to store the shifts.


 // b.3) Functions that perform all calculations in one step.

  int FindSchurVectors() {
    throw ArpackError(ArpackError::SCHUR_UNDEFINED, "FindSchurVectors");
    return 0;  // Only to avoid warning messages emitted by some compilers.
  }
  // For symmetric problems, Schur vectors are eigenvectors.

  int Eigenvalues(FLOAT* &EigValp, bool ivec = false, bool ischur = false);
  // Overrides array EigValp with the eigenvalues of the problem.
  // Also calculates eigenvectors and Schur vectors if requested.

  int EigenValVectors(FLOAT* &EigVecp, FLOAT* &EigValp, bool ischur = false);
  // Overrides array EigVecp sequentially with the eigenvectors of the
  // given eigen-problem. Also stores the eigenvalues in EigValp.
  // Calculates Schur vectors if requested.


 // b.4) Functions that return elements of vectors and matrices.

  FLOAT Eigenvalue(int i);
  // Provides i-eth eigenvalue.

  FLOAT Eigenvector(int i, int j);
  // Provides element j of the i-eth eigenvector.


 // b.5) Functions that use STL vector class.

#ifdef VECTOR_H

  vector<FLOAT>* StlEigenvalues(bool ivec = false, bool ischur = false);
  // Calculates the eigenvalues and stores them in a single STL vector.
  // Also calculates eigenvectors and Schur vectors if requested.

  vector<FLOAT>* StlEigenvector(int i);
  // Returns the i-th eigenvector in a STL vector.

#endif // #ifdef VECTOR_H.


 // b.6) Constructors and destructor.

  ARrcSymStdEig() { }
  // Short constructor.

  ARrcSymStdEig(int np, int nevp, char* whichp = "LM", int ncvp = 0,
                FLOAT tolp = 0.0, int maxitp = 0, FLOAT* residp = NULL,
                bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcSymStdEig(int np, int nevp, FLOAT sigmap, char* whichp = "LM",
                int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARrcSymStdEig(const ARrcSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcSymStdEig() { }
  // Destructor.

 // c) Operators.

  ARrcSymStdEig& operator=(const ARrcSymStdEig& other);
  // Assignment operator.

}; // class ARrcSymStdEig.


// ------------------------------------------------------------------------ //
// ARrcSymStdEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline void ARrcSymStdEig<FLOAT>::WorkspaceAllocate()
{

  this->lworkl  = this->ncv*(this->ncv+9);
  this->lworkv  = 0;
  this->lrwork  = 0;
  this->workl   = new FLOAT[this->lworkl+1];

} // WorkspaceAllocate.


template<class FLOAT>
inline void ARrcSymStdEig<FLOAT>::Aupp()
{

  saupp(this->ido, this->bmat, this->n, this->which, this->nev, this->tol, this->resid, this->ncv, this->V, this->n,
        this->iparam, this->ipntr, this->workd, this->workl, this->lworkl, this->info);

} // Aupp.


template<class FLOAT>
inline void ARrcSymStdEig<FLOAT>::Eupp()
{

  seupp(this->rvec, this->HowMny, this->EigValR, this->EigVec, this->n, this->sigmaR, this->bmat,
        this->n, this->which, this->nev, this->tol, this->resid, this->ncv, this->V, this->n, this->iparam,
        this->ipntr, this->workd, this->workl, this->lworkl, this->info);

} // Eupp.


template<class FLOAT>
const char* ARrcSymStdEig<FLOAT>::CheckWhich(const char* whichp)
{

  switch (whichp[0]) {
  case 'B':                       // The options are: BE, ...
    return "BE";
  case 'L':                       // LA, LM, ...
  case 'S':                       // SA, SM.
    switch (whichp[1]){
    case 'A':
    case 'M':
      return whichp;
    }
  default:
    throw ArpackError(ArpackError::WHICH_UNDEFINED);
  }

} // CheckWhich.


template<class FLOAT>
FLOAT* ARrcSymStdEig<FLOAT>::PutVector()
{

  switch (this->ido) {
  case -1:
  case  1:                    // Returning OP*x.
  case  2:
    return &this->workd[this->ipntr[2]];  // Returning B*x.
  case  3:
    return &this->workl[this->ipntr[11]]; // Returning shifts.
  default:
    throw ArpackError(ArpackError::CANNOT_PUT_VECTOR, "PutVector");
  }

} // PutVector.


template<class FLOAT>
int ARrcSymStdEig<FLOAT>::
Eigenvalues(FLOAT* &EigValp, bool ivec, bool ischur)
{

  if (this->ValuesOK) {                      // Eigenvalues are available.
    if (EigValp == NULL) {             // Moving eigenvalues.
      EigValp  = this->EigValR;
      this->EigValR  = NULL;
      this->newVal   = false;
      this->ValuesOK = false;
    }
    else {                             // Copying eigenvalues.
      copy(this->nconv,this->EigValR,1,EigValp,1);
    }
  }
  else {                               // Eigenvalues are not available.
    if (this->newVal) {
      delete[] this->EigValR;
      this->newVal = false;
    }
    if (EigValp == NULL) {
      try { EigValp = new FLOAT[this->ValSize()]; }
      catch (ArpackError) { return 0; }
    }
    this->EigValR = EigValp;
    if (ivec) {                        // Finding eigenvalues and eigenvectors.
      this->nconv = this->FindEigenvectors(ischur);
    }
    else {                             // Finding eigenvalues only.
      this->nconv = this->FindEigenvalues();
    }
    this->EigValR = NULL;
  }
  return this->nconv;

} // Eigenvalues(EigValp, ivec, ischur).


template<class FLOAT>
int ARrcSymStdEig<FLOAT>::
EigenValVectors(FLOAT* &EigVecp, FLOAT* &EigValp, bool ischur)
{

  if (this->ValuesOK) {                  // Eigenvalues are already available .
    this->nconv = Eigenvalues(EigValp, false);
    this->nconv = Eigenvectors(EigVecp, ischur);
  }
  else {                           // Eigenvalues and vectors are not available.
    try {
      if (EigVecp == NULL) EigVecp = new FLOAT[this->ValSize()*this->n];
      if (EigValp == NULL) EigValp = new FLOAT[this->ValSize()];
    }
    catch (ArpackError) { return 0; }
    if (this->newVec) {
      delete[] this->EigVec;
      this->newVec = false;
    }
    if (this->newVal) {
      delete[] this->EigValR;
      this->newVal = false;
    }
    this->EigVec  = EigVecp;
    this->EigValR = EigValp;
    this->nconv   = this->FindEigenvectors(ischur);
    this->EigVec  = NULL;
    this->EigValR = NULL;
  }
  return this->nconv;

} // EigenValVectors(EigVecp, EigValp, ischur).


template<class FLOAT>
inline FLOAT ARrcSymStdEig<FLOAT>::Eigenvalue(int i)
{

  // Returning i-eth eigenvalue.

  if (!this->ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "Eigenvalue(i)");
  }
  else if ((i>=this->nconv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvalue(i)");
  }
  return this->EigValR[i];

} // Eigenvalue(i).


template<class FLOAT>
inline FLOAT ARrcSymStdEig<FLOAT>::Eigenvector(int i, int j)
{

  // Returning element j of i-eth eigenvector.

  if (!this->VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "Eigenvector(i,j)");
  }
  else if ((i>=this->nconv)||(i<0)||(j>=this->n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvector(i,j)");
  }
  return this->EigVec[i*this->n+j];

} // Eigenvector(i,j).


#ifdef VECTOR_H

template<class FLOAT>
inline vector<FLOAT>* ARrcSymStdEig<FLOAT>::
StlEigenvalues(bool ivec, bool ischur)
{

  // Returning the eigenvalues in a STL vector.

  vector<FLOAT>* StlEigValR;
  FLOAT*         ValPtr;

  try { StlEigValR = new vector<FLOAT>(ValSize()); }
  catch (ArpackError) { return NULL; }
  ValPtr = StlEigValR->begin();
  this->nconv = Eigenvalues(ValPtr, ivec, ischur);
  return StlEigValR;

} // StlEigenvalues.


template<class FLOAT>
inline vector<FLOAT>* ARrcSymStdEig<FLOAT>::StlEigenvector(int i)
{

  // Returning the i-th eigenvector in a STL vector.

  vector<FLOAT>* Vec;

  if (!this->VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvector(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvector(i)");
  }
  try {
    Vec = new vector<FLOAT>(&this->EigVec[i*this->n], &this->EigVec[(i+1)*this->n]);
  }
  catch (ArpackError) { return NULL; }
  return Vec;

} // StlEigenvector(i).

#endif // #ifdef VECTOR_H.


template<class FLOAT>
inline ARrcSymStdEig<FLOAT>::
ARrcSymStdEig(int np, int nevp, char* whichp, int ncvp,
              FLOAT tolp, int maxitp, FLOAT* residp, bool ishiftp)

{

  this->NoShift();
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARrcSymStdEig<FLOAT>::
ARrcSymStdEig(int np, int nevp, FLOAT sigmap, char* whichp,
              int ncvp, FLOAT tolp, int maxitp, FLOAT* residp,
              bool ishiftp)

{

  this->ChangeShift(sigmap);
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class FLOAT>
ARrcSymStdEig<FLOAT>& ARrcSymStdEig<FLOAT>::
operator=(const ARrcSymStdEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRSSYM_H

