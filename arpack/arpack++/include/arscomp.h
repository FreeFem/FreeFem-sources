/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARSComp.h.
   Arpack++ class ARCompStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARSCOMP_H
#define ARSCOMP_H

#include <stddef.h>
#include "arch.h"
#include "arseig.h"
#include "arrscomp.h"

template<class FLOAT, class FOP>
class ARCompStdEig:
  virtual public ARStdEig<FLOAT, arcomplex<FLOAT>, FOP>,
  virtual public ARrcCompStdEig<FLOAT> {

 public:

 // a) Constructors and destructor.

  ARCompStdEig() { }
  // Short constructor.

  ARCompStdEig(int np, int nevp, FOP* objOPp,
               void (FOP::* MultOPxp)(arcomplex<FLOAT>[],arcomplex<FLOAT>[]),
               char* whichp = "LM", int ncvp = 0,
               FLOAT tolp = 0.0, int maxitp = 0,
               arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARCompStdEig(int np, int nevp, FOP* objOPp,
               void (FOP::* MultOPxp)(arcomplex<FLOAT>[],arcomplex<FLOAT>[]),
               arcomplex<FLOAT> sigma,  char* whichp = "LM",
               int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
               arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARCompStdEig(const ARCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARCompStdEig() { }
  // Destructor.

 // b) Operators.

  ARCompStdEig& operator=(const ARCompStdEig& other);
  // Assignment operator.

}; // class ARCompStdEig.


// ------------------------------------------------------------------------ //
// ARCompStdEig member functions definition.                                //
// ------------------------------------------------------------------------ //


template<class FLOAT, class FOP>
inline ARCompStdEig<FLOAT, FOP>::
ARCompStdEig(int np, int nevp, FOP* objOPp,
             void (FOP::* MultOPxp)(arcomplex<FLOAT>[],arcomplex<FLOAT>[]),
             char* whichp, int ncvp, FLOAT tolp, int maxitp,
             arcomplex<FLOAT>* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT, class FOP>
inline ARCompStdEig<FLOAT, FOP>::
ARCompStdEig(int np, int nevp, FOP* objOPp,
             void (FOP::* MultOPxp)(arcomplex<FLOAT>[],arcomplex<FLOAT>[]),
             arcomplex<FLOAT> sigmap, char* whichp, int ncvp,
             FLOAT tolp, int maxitp, arcomplex<FLOAT>* residp,
             bool ishiftp)

{

  ChangeShift(sigmap);
  DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class FLOAT, class FOP>
ARCompStdEig<FLOAT, FOP>& ARCompStdEig<FLOAT, FOP>::
operator=(const ARCompStdEig<FLOAT, FOP>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARSCOMP_H

