/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARSSym.h.
   Arpack++ class ARSymStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARSSYM_H
#define ARSSYM_H

#include <stddef.h>
#include "arch.h"
#include "arseig.h"
#include "arrssym.h"


template<class FLOAT, class FOP>
class ARSymStdEig:
  public virtual ARStdEig<FLOAT, FLOAT, FOP>,
  public virtual ARrcSymStdEig<FLOAT> {

 public:

 // a) Constructors and destructor.

  ARSymStdEig() { }
  // Short constructor.

  ARSymStdEig(int np, int nevp, FOP* objOPp,
              void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
              char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
              int maxitp = 0, FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARSymStdEig(int np, int nevp, FOP* objOPp,
              void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
              FLOAT sigmap, char* whichp = "LM", int ncvp = 0,
              FLOAT tolp = 0.0, int maxitp = 0, FLOAT* residp = NULL,
              bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARSymStdEig(const ARSymStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARSymStdEig() { }
  // Destructor.

 // b) Operators.

  ARSymStdEig& operator=(const ARSymStdEig& other);
  // Assignment operator.

}; // class ARSymStdEig.


// ------------------------------------------------------------------------ //
// ARSymStdEig member functions definition.                                 //
// ------------------------------------------------------------------------ //


template<class FLOAT, class FOP>
inline ARSymStdEig<FLOAT, FOP>::
ARSymStdEig(int np, int nevp, FOP* objOPp,
            void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
            char* whichp, int ncvp, FLOAT tolp,
            int maxitp, FLOAT* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT, class FOP>
inline ARSymStdEig<FLOAT, FOP>::
ARSymStdEig(int np, int nevp, FOP* objOPp,
            void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
            FLOAT sigmap, char* whichp, int ncvp, FLOAT tolp,
            int maxitp, FLOAT* residp, bool ishiftp)

{

  ChangeShift(sigmap);
  DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class FLOAT, class FOP>
ARSymStdEig<FLOAT, FOP>& ARSymStdEig<FLOAT, FOP>::
operator=(const ARSymStdEig<FLOAT, FOP>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARSSYM_H

