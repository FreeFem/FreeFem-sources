/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARSNSym.h.
   Arpack++ class ARNonSymStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARSNSYM_H
#define ARSNSYM_H

#include <stddef.h>
#include "arch.h"
#include "arseig.h"
#include "arrsnsym.h"


template<class FLOAT, class FOP>
class ARNonSymStdEig:
  public virtual ARStdEig<FLOAT, FLOAT, FOP>,
  public virtual ARrcNonSymStdEig<FLOAT> {

 public:

 // a) Constructors and destructor.

  ARNonSymStdEig() { }
  // Short constructor.

  ARNonSymStdEig(int np, int nevp, FOP* objOPp,
                 void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
                 char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
                 int maxitp = 0, FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARNonSymStdEig(int np, int nevp, FOP* objOPp,
                 void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
                 FLOAT sigma, char* whichp = "LM", int ncvp = 0,
                 FLOAT tolp = 0.0, int maxitp = 0, FLOAT* residp = NULL,
                 bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARNonSymStdEig(const ARNonSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARNonSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARNonSymStdEig& operator=(const ARNonSymStdEig& other);
  // Assignment operator.

}; // class ARNonSymStdEig.


// ------------------------------------------------------------------------ //
// ARNonSymStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class FLOAT, class FOP>
inline ARNonSymStdEig<FLOAT, FOP>::
ARNonSymStdEig(int np, int nevp, FOP* objOPp,
               void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
               char* whichp, int ncvp, FLOAT tolp, int maxitp,
               FLOAT* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT, class FOP>
inline ARNonSymStdEig<FLOAT, FOP>::
ARNonSymStdEig(int np, int nevp, FOP* objOPp,
               void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
               FLOAT sigmap, char* whichp, int ncvp, FLOAT tolp,
               int maxitp, FLOAT* residp, bool ishiftp)

{

  ChangeShift(sigmap);
  DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class FLOAT, class FOP>
ARNonSymStdEig<FLOAT, FOP>& ARNonSymStdEig<FLOAT, FOP>::
operator=(const ARNonSymStdEig<FLOAT, FOP>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARSNSYM_H

