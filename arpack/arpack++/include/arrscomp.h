/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARRSComp.h.
   Arpack++ class ARrcCompStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRSCOMP_H
#define ARRSCOMP_H

#include <stddef.h>
#include "arch.h"
#include "arerror.h"
#include "debug.h"
#include "arrseig.h"
#include "caupp.h"
#include "ceupp.h"

template<class FLOAT>
class ARrcCompStdEig: virtual public ARrcStdEig<FLOAT, arcomplex<FLOAT> > {

 protected:

 // a) Protected functions:

 // a.1) Memory control functions.

  void WorkspaceAllocate();
  // Allocates workspace for complex problems.


 // a.2) Functions that handle original FORTRAN ARPACK code.

  void Aupp();
  // Interface to FORTRAN subroutines CNAUPD and ZNAUPD.

  void Eupp();
  // Interface to FORTRAN subroutines CNEUPD and ZNEUPD.

 public:

 // b) Public functions:

 // b.1) Trace functions.

  void Trace(const int digit = -5, const int getv0 = 0, const int aupd = 1,
             const int aup2 = 0,  const int aitr = 0,  const int eigt = 0,
             const int apps = 0,  const int gets = 0,  const int eupd = 0)
  { 
    cTraceOn(digit, getv0, aupd, aup2, aitr, eigt, apps, gets, eupd); 
  }
  // Turns on trace mode. 


 // b.2) Functions that perform all calculations in one step.

  int Eigenvalues(arcomplex<FLOAT>* &EigValp, bool ivec = false,
                   bool ischur = false);
  // Overrides array EigValp with the eigenvalues of the problem.
  // Also calculates eigenvectors and Schur vectors if requested.

  int EigenValVectors(arcomplex<FLOAT>* &EigVecp, arcomplex<FLOAT>* &EigValp,
                       bool ischur = false);
  // Overrides array EigVecp sequentially with the eigenvectors of the
  // given eigen-problem. Also stores the eigenvalues in EigValp.
  // Calculates Schur vectors if requested.


 // b.3) Functions that return elements of vectors and matrices.

  arcomplex<FLOAT> Eigenvalue(int i);
  // Provides i-eth eigenvalue.

  arcomplex<FLOAT> Eigenvector(int i, int j);
  // Provides element j of the i-eth eigenvector.


 // b.4) Functions that use STL vector class.

#ifdef VECTOR_H

  vector<arcomplex<FLOAT> >* StlEigenvalues(bool ivec = false,
                                            bool ischur = false);
  // Calculates the eigenvalues and stores them in a single STL vector.
  // Also calculates eigenvectors and Schur vectors if requested.

  vector<arcomplex<FLOAT> >* StlEigenvector(int i);
  // Returns the i-th eigenvector in a STL vector.

#endif // #ifdef VECTOR_H.


 // b.5) Constructors and destructor.

  ARrcCompStdEig() { }
  // Short constructor.

  ARrcCompStdEig(int np, int nevp, char* whichp = "LM",
                 int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcCompStdEig(int np, int nevp, arcomplex<FLOAT> sigma,
                 char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
                 int maxitp = 0, arcomplex<FLOAT>* residp = NULL,
                 bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARrcCompStdEig(const ARrcCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcCompStdEig() { }
  // Destructor.

 // c) Operators.

  ARrcCompStdEig& operator=(const ARrcCompStdEig& other);
  // Assignment operator.

}; // class ARrcCompStdEig.


// ------------------------------------------------------------------------ //
// ARrcCompStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline void ARrcCompStdEig<FLOAT>::WorkspaceAllocate()
{

  this->lworkl  = this->ncv*(3*this->ncv+6);
  this->lworkv  = 2*this->ncv;
  this->lrwork  = this->ncv;
  this->workl   = new arcomplex<FLOAT>[this->lworkl+1];
  this->workv   = new arcomplex<FLOAT>[this->lworkv+1];
  this->rwork   = new FLOAT[this->lrwork+1];

} // WorkspaceAllocate.


template<class FLOAT>
inline void ARrcCompStdEig<FLOAT>::Aupp()
{

  caupp( this->ido, this->bmat, this->n, this->which, this->nev, this->tol, this->resid, this->ncv, this->V, this->n,
        this->iparam, this->ipntr, this->workd, this->workl, this->lworkl, this->rwork, this->info);

} // Aupp.


template<class FLOAT>
inline void ARrcCompStdEig<FLOAT>::Eupp()
{

  ceupp( this->rvec, this->HowMny, this->EigValR, this->EigVec, this->n, this->sigmaR, this->workv,
         this->bmat, this->n, this->which, this->nev, this->tol, this->resid, this->ncv, this->V, this->n, this->iparam,
        this->ipntr, this->workd, this->workl, this->lworkl, this->rwork, this->info);

} // Eupp.


template<class FLOAT>
int ARrcCompStdEig<FLOAT>::
Eigenvalues(arcomplex<FLOAT>* &EigValp, bool ivec, bool ischur)
{

  if (this->ValuesOK) {                      // Eigenvalues are available .
    if (EigValp == NULL) {             // Moving eigenvalues.
      EigValp  = this->EigValR;
      this->EigValR  = NULL;
      this->newVal   = false;
      this->ValuesOK = false;
    }
    else {                             // Copying eigenvalues.
      copy(this->nconv,this->EigValR,1,this->EigValp,1);
    }
  }
  else {
    if (this->newVal) {
      delete[] this->EigValR;
      this->newVal = false;
    }
    if (EigValp == NULL) {
      try { EigValp = new arcomplex<FLOAT>[this->ValSize()]; }
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
int ARrcCompStdEig<FLOAT>::
EigenValVectors(arcomplex<FLOAT>* &EigVecp, arcomplex<FLOAT>* &EigValp,
                bool ischur)
{

  if (this->ValuesOK) {                  // Eigenvalues are already available.
    this->nconv = Eigenvalues(EigValp, false);
    this->nconv = Eigenvectors(EigVecp, ischur);
  }
  else {                           // Eigenvalues and vectors are not available.
    if (this->newVec) {
      delete[] this->EigVec;
      this->newVec = false;
    }
    if (this->newVal) {
      delete[] this->EigValR;
      this->newVal = false;
    }  
    try {
      if (EigVecp == NULL) EigVecp = new arcomplex<FLOAT>[this->ValSize()*this->n];
      if (EigValp == NULL) EigValp = new arcomplex<FLOAT>[this->ValSize()];
    }
    catch (ArpackError) { return 0; }
    this->EigVec  = EigVecp;
    this->EigValR = EigValp;
    this->nconv   = this->FindEigenvectors(ischur);
    this->EigVec  = NULL;
    this->EigValR = NULL;
  }
  return this->nconv;

} // EigenValVectors(EigVecp, EigValp, ischur).


template<class FLOAT>
inline arcomplex<FLOAT> ARrcCompStdEig<FLOAT>::Eigenvalue(int i)
// calcula e retorna um autovalor.

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
inline arcomplex<FLOAT> ARrcCompStdEig<FLOAT>::
Eigenvector(int i, int j)
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
inline vector<arcomplex<FLOAT> >* ARrcCompStdEig<FLOAT>::
StlEigenvalues(bool ivec, bool ischur)
{

  // Returning the eigenvalues in a STL vector.

  vector<arcomplex<FLOAT> >* ValR;
  arcomplex<FLOAT>*          ValPtr;

  try {
    ValR = new vector<arcomplex<FLOAT> >(ValSize());
  }
  catch (ArpackError) { return NULL; }
  ValPtr = ValR->begin();
  nconv = Eigenvalues(ValPtr, ivec, ischur);
  return ValR;

} // StlEigenvalues.


template<class FLOAT>
inline vector<arcomplex<FLOAT> >* ARrcCompStdEig<FLOAT>::
StlEigenvector(int i)
{

  // Returning the i-th eigenvector in a STL vector.

  vector<arcomplex<FLOAT> >* Vec;

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvector(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvector(i)");
  }
  try {
    Vec = new vector<arcomplex<FLOAT> >(&EigVec[i*n], &EigVec[(i+1)*n]);
  }
  catch (ArpackError) { return NULL; }
  return Vec;

} // StlEigenvector(i).

#endif // #ifdef VECTOR_H.


template<class FLOAT>
inline ARrcCompStdEig<FLOAT>::
ARrcCompStdEig(int np, int nevp, char* whichp, int ncvp, FLOAT tolp,
               int maxitp, arcomplex<FLOAT>* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARrcCompStdEig<FLOAT>::
ARrcCompStdEig(int np, int nevp, arcomplex<FLOAT> sigmap,
               char* whichp, int ncvp, FLOAT tolp, int maxitp,
               arcomplex<FLOAT>* residp, bool ishiftp)

{

  ChangeShift(sigmap);
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class FLOAT>
ARrcCompStdEig<FLOAT>& ARrcCompStdEig<FLOAT>::
operator=(const ARrcCompStdEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRSCOMP_H

