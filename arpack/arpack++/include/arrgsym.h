/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARRGSym.h.
   Arpack++ class ARrcSymGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRGSYM_H
#define ARRGSYM_H

#include <stddef.h>
#include "arch.h"
#include "arrssym.h"
#include "arrgeig.h"

template<class FLOAT>
class ARrcSymGenEig:
  virtual public ARrcGenEig<FLOAT, FLOAT>,
  virtual public ARrcSymStdEig<FLOAT> {

 protected:

 // a) Protected variable:

  char    InvertMode;


 // b) Protected functions:

  char CheckInvertMode(char InvertModep);

  virtual void Copy(const ARrcSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  void ChangeInvertMode(char InvertModep);
  // Changes "InvertMode" to 'S' (shift-and-invert),
  // 'B' (buckling) or 'C' (cayley) mode.

  virtual void ChangeShift(FLOAT sigmap);
  // Changes shift value.

  virtual void SetShiftInvertMode(FLOAT sigmap);
  // Turns problem to shift and invert mode with shift defined by sigmap.

  virtual void SetBucklingMode(FLOAT sigmap);
  // Turns problem to buckling mode with shift defined by sigmap.

  virtual void SetCayleyMode(FLOAT sigmap);
  // Turns problem to Cayley mode with shift defined by sigmap.


 // c.2) Constructors and destructor.

  ARrcSymGenEig() { InvertMode = 'S'; }
  // Short constructor that does almost nothing.

  ARrcSymGenEig(int np, int nevp, char* whichp = "LM",
                int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcSymGenEig(char invertmodep, int np, int nevp, FLOAT sigmap,
                char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
                int maxitp = 0, FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift-and-invert, buckling and Cayley modes).

  ARrcSymGenEig(const ARrcSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARrcSymGenEig& operator=(const ARrcSymGenEig& other);
  // Assignment operator.

}; // class ARrcSymGenEig.


// ------------------------------------------------------------------------ //
// ARrcSymGenEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline char ARrcSymGenEig<FLOAT>::CheckInvertMode(char InvertModep)
{
  if ((InvertModep != 'S') && (InvertModep != 'B') && (InvertModep != 'C')) {
    throw ArpackError(ArpackError::INVMODE_UNDEFINED);
  }
  return InvertModep;

} // CheckInvertMode.


template<class FLOAT>
inline void ARrcSymGenEig<FLOAT>::Copy(const ARrcSymGenEig<FLOAT>& other)
{

  ARrcStdEig<FLOAT, FLOAT>::Copy(other);
  InvertMode = other.InvertMode;

} // Copy.


template<class FLOAT>
inline void ARrcSymGenEig<FLOAT>::ChangeInvertMode(char InvertModep)
{

  InvertMode = CheckInvertMode(InvertModep);
  switch (InvertMode) {
  case 'S':
    mode    = 3;    // Shift and invert mode.
    break;
  case 'B':
    mode    = 4;    // Buckling mode.
    break;
  case 'C':
    mode    = 5;    // Cayley mode.
    break;
  }
  iparam[7] = mode;
  Restart();

} // ChangeInvertMode.


template<class FLOAT>
inline void ARrcSymGenEig<FLOAT>::ChangeShift(FLOAT sigmap)
{

  sigmaR    = sigmap;
  sigmaI    = 0.0;
  ChangeInvertMode(InvertMode);

} // ChangeShift.


template<class FLOAT>
void ARrcSymGenEig<FLOAT>::SetShiftInvertMode(FLOAT sigmap)

{

  InvertMode = 'S';
  ChangeShift(sigmap);

} // SetShiftInvertMode.


template<class FLOAT>
void ARrcSymGenEig<FLOAT>::SetBucklingMode(FLOAT sigmap)

{

  InvertMode = 'B';
  ChangeShift(sigmap);

} // SetBucklingMode.


template<class FLOAT>
void ARrcSymGenEig<FLOAT>::SetCayleyMode(FLOAT sigmap)

{

  InvertMode = 'C';
  ChangeShift(sigmap);

} // SetCayleyMode.


template<class FLOAT>
inline ARrcSymGenEig<FLOAT>::
ARrcSymGenEig(int np, int nevp, char* whichp, int ncvp, FLOAT tolp,
              int maxitp, FLOAT* residp, bool ishiftp)

{

  InvertMode = 'S';   // Considering mode = 3 in ChangeShift.
  NoShift();
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARrcSymGenEig<FLOAT>::
ARrcSymGenEig(char InvertModep, int np, int nevp,
              FLOAT sigmap, char* whichp, int ncvp, FLOAT tolp,
              int maxitp, FLOAT* residp, bool ishiftp)

{

  InvertMode = CheckInvertMode(InvertModep); // InvertMode = 'S', 'B', 'C'.
  ChangeShift(sigmap);
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift-and-invert, buckling and Cayley modes).


template<class FLOAT>
ARrcSymGenEig<FLOAT>& ARrcSymGenEig<FLOAT>::
operator=(const ARrcSymGenEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRGSYM_H

