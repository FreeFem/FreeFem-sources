/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARRGComp.h.
   Arpack++ class ARrcCompGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRGCOMP_H
#define ARRGCOMP_H

#include <stddef.h>
#include "arch.h"
#include "arrscomp.h"
#include "arrgeig.h"

template<class FLOAT>
class ARrcCompGenEig:
  virtual public ARrcGenEig<FLOAT, arcomplex<FLOAT> >,
  virtual public ARrcCompStdEig<FLOAT>  {

 public:

  // a) Constructors and destructor.

  ARrcCompGenEig() { }
  // Short constructor (Does nothing but calling base classes constructors).

  ARrcCompGenEig(int np, int nevp, char* whichp = "LM",
                 int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcCompGenEig(int np, int nevp, arcomplex<FLOAT> sigmap,
                 char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
                 int maxitp = 0, arcomplex<FLOAT>* residp = NULL,
                 bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARrcCompGenEig(const ARrcCompGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcCompGenEig() { }
  // Destructor.

 // b) Operators.

  ARrcCompGenEig& operator=(const ARrcCompGenEig& other);
  // Assignment operator.

}; // class ARrcCompGenEig.


// ------------------------------------------------------------------------ //
// ARrcCompGenEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline ARrcCompGenEig<FLOAT>::
ARrcCompGenEig(int np, int nevp, char* whichp, int ncvp, FLOAT tolp,
               int maxitp, arcomplex<FLOAT>* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARrcCompGenEig<FLOAT>::
ARrcCompGenEig(int np, int nevp, arcomplex<FLOAT> sigmap, char* whichp,
               int ncvp, FLOAT tolp, int maxitp, arcomplex<FLOAT>* residp,
               bool ishiftp)

{

  ChangeShift(sigmap);
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shif and invert mode).


template<class FLOAT>
ARrcCompGenEig<FLOAT>& ARrcCompGenEig<FLOAT>::
operator=(const ARrcCompGenEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRGCOMP_H

