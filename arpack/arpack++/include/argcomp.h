/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARGComp.h.
   Arpack++ class ARCompGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARGCOMP_H
#define ARGCOMP_H

#include <stddef.h>
#include "arch.h"
#include "arscomp.h"
#include "argeig.h"

template<class FLOAT, class FOP, class FB>
class ARCompGenEig:
  virtual public ARGenEig<FLOAT, arcomplex<FLOAT>, FOP, FB>,
  virtual public ARCompStdEig<FLOAT, FOP>  {

 public:

  // a) Constructors and destructor.

  ARCompGenEig() { }
  // Short constructor (Does nothing but calling base classes constructors).

  ARCompGenEig(int np, int nevp, FOP* objOPp,
               void (FOP::* MultOPxp)(arcomplex<FLOAT>[], arcomplex<FLOAT>[]),
               FB* objBp,
               void (FB::* MultBxp)(arcomplex<FLOAT>[], arcomplex<FLOAT>[]),
               char* whichp = "LM", int ncvp = 0,
               FLOAT tolp = 0.0, int maxitp = 0,
               arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARCompGenEig(int np, int nevp, FOP* objOPp,
               void (FOP::* MultOPxp)(arcomplex<FLOAT>[], arcomplex<FLOAT>[]),
               FB* objBp,
               void (FB::* MultBxp)(arcomplex<FLOAT>[], arcomplex<FLOAT>[]),
               arcomplex<FLOAT> sigmap,
               char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
               int maxitp = 0, arcomplex<FLOAT>* residp = NULL,
               bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARCompGenEig(const ARCompGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARCompGenEig() { }
  // Destructor.

 // b) Operators.

  ARCompGenEig& operator=(const ARCompGenEig& other);
  // Assignment operator.

}; // class ARCompGenEig.


// ------------------------------------------------------------------------ //
// ARCompGenEig member functions definition.                                //
// ------------------------------------------------------------------------ //


template<class FLOAT, class FOP, class FB>
inline ARCompGenEig<FLOAT, FOP, FB>::
ARCompGenEig(int np, int nevp, FOP* objOPp,
             void (FOP::* MultOPxp)(arcomplex<FLOAT>[],arcomplex<FLOAT>[]),
             FB* objBp,
             void (FB::* MultBxp)(arcomplex<FLOAT>[], arcomplex<FLOAT>[]),
             char* whichp, int ncvp, FLOAT tolp,
             int maxitp, arcomplex<FLOAT>* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT, class FOP, class FB>
inline ARCompGenEig<FLOAT, FOP, FB>::
ARCompGenEig(int np, int nevp, FOP* objOPp,
             void (FOP::* MultOPxp)(arcomplex<FLOAT>[],arcomplex<FLOAT>[]),
             FB* objBp,
             void (FB::* MultBxp)(arcomplex<FLOAT>[], arcomplex<FLOAT>[]),
             arcomplex<FLOAT> sigmap, char* whichp, int ncvp, FLOAT tolp,
             int maxitp, arcomplex<FLOAT>* residp, bool ishiftp)

{

  ChangeShift(sigmap);
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class FLOAT, class FOP, class FB>
ARCompGenEig<FLOAT, FOP, FB>& ARCompGenEig<FLOAT, FOP, FB>::
operator=(const ARCompGenEig<FLOAT, FOP, FB>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARGCOMP_H

