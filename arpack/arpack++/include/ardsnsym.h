/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARDSNSym.h.
   Arpack++ class ARluNonSymStdEig definition
   (dense matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDSNSYM_H
#define ARDSNSYM_H

#include <stddef.h>
#include "arch.h"
#include "arsnsym.h"
#include "ardnsmat.h"


template<class FLOAT>
class ARluNonSymStdEig:
  public virtual ARNonSymStdEig<FLOAT, ARdsNonSymMatrix<FLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(FLOAT sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(FLOAT sigmap);

 // a.2) Constructors and destructor.

  ARluNonSymStdEig() { }
  // Short constructor.

  ARluNonSymStdEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                   char* whichp = "LM", int ncvp = 0,
                   FLOAT tolp = 0.0, int maxitp = 0,
                   FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluNonSymStdEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                   FLOAT sigma, char* whichp = "LM", int ncvp = 0,
                   FLOAT tolp = 0.0, int maxitp = 0,
                   FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluNonSymStdEig(const ARluNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluNonSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARluNonSymStdEig& operator=(const ARluNonSymStdEig& other);
  // Assignment operator.

}; // class ARluNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARluNonSymStdEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline void ARluNonSymStdEig<FLOAT>::
ChangeShift(FLOAT sigmaRp)
{

  sigmaR    = sigmaRp;
  sigmaI    = 0.0;
  mode      = 3;
  iparam[7] = mode;

  objOP->FactorAsI(sigmaR);
  Restart();

} // ChangeShift.


template<class FLOAT>
inline void ARluNonSymStdEig<FLOAT>::SetRegularMode()
{

  ARStdEig<FLOAT, FLOAT, ARdsNonSymMatrix<FLOAT> >::
    SetRegularMode(objOP, &ARdsNonSymMatrix<FLOAT>::MultMv);

} // SetRegularMode.


template<class FLOAT>
inline void ARluNonSymStdEig<FLOAT>::SetShiftInvertMode(FLOAT sigmap)
{

  ARStdEig<FLOAT, FLOAT, ARdsNonSymMatrix<FLOAT> >::
    SetShiftInvertMode(sigmap, objOP, &ARdsNonSymMatrix<FLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class FLOAT>
inline ARluNonSymStdEig<FLOAT>::
ARluNonSymStdEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                 char* whichp, int ncvp, FLOAT tolp,
                 int maxitp, FLOAT* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(A.ncols(), nevp, &A, ARdsNonSymMatrix<FLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARluNonSymStdEig<FLOAT>::
ARluNonSymStdEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                 FLOAT sigmap, char* whichp, int ncvp, FLOAT tolp,
                 int maxitp, FLOAT* residp, bool ishiftp)

{

  DefineParameters(A.ncols(), nevp, &A, &ARdsNonSymMatrix<FLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class FLOAT>
ARluNonSymStdEig<FLOAT>& ARluNonSymStdEig<FLOAT>::
operator=(const ARluNonSymStdEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDSNSYM_H

