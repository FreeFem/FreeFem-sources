/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARRGNSym.h.
   Arpack++ class ARrcNonSymGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRGNSYM_H
#define ARRGNSYM_H

#include <stddef.h>
#include "arch.h"
#include "arrsnsym.h"
#include "arrgeig.h"

template<class FLOAT>
class ARrcNonSymGenEig:
  virtual public ARrcGenEig<FLOAT, FLOAT>,
  virtual public ARrcNonSymStdEig<FLOAT>  {

 protected:

 // a) Protected variables:

  char part;


 // b) Protected functions:

  char CheckPart(char partp);

  virtual void Copy(const ARrcNonSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that provides access to internal variables' values.

  FLOAT GetShiftImag() { return sigmaI; }
  // Returns the imaginary part of the shift (when in shift and invert mode).


 // c.2) Functions that allow changes in problem parameters.

  void ChangePart(char partp);
  // Changes "part" to 'R' (real) or 'I' (imaginary).

  virtual void ChangeShift(FLOAT sigmaRp, FLOAT sigmaIp = 0.0);
  // Turns the problem to shift-and-invert mode
  // with shift defined by sigmaRp and sigmaIp.

  virtual void SetShiftInvertMode(FLOAT sigmaRp);
  // Turns the problem to real shift-and-invert mode with sigmaRp as shift.

  virtual void SetComplexShiftMode(char partp, FLOAT sigmaRp, FLOAT sigmaIp);
  // Turns the problem to complex shift-and-invert mode with shift
  // defined by sigmaRp and sigmaIp.


 // c.3) Constructors and destructor.

  ARrcNonSymGenEig() { part = 'R'; }
  // Short constructor that does almost nothing.

  ARrcNonSymGenEig(int np, int nevp, char* whichp = "LM",
                   int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                   FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcNonSymGenEig(int np, int nevp, FLOAT sigmap,
                   char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
                   int maxitp = 0, FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (real shift and invert mode).

  ARrcNonSymGenEig(int np, int nevp,
                   char partp, FLOAT sigmaRp, FLOAT sigmaIp,
                   char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
                   int maxitp = 0, FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (complex shift and invert mode).

  ARrcNonSymGenEig(const ARrcNonSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcNonSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARrcNonSymGenEig& operator=(const ARrcNonSymGenEig& other);
  // Assignment operator.

}; // class ARrcNonSymGenEig.


// ------------------------------------------------------------------------ //
// ARrcNonSymGenEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline char ARrcNonSymGenEig<FLOAT>::CheckPart(char partp)
{
  if ((partp != 'R') && (partp != 'I')) {
    throw ArpackError(ArpackError::PART_UNDEFINED);
  }
  return partp;
} // CheckPart.


template<class FLOAT>
inline void ARrcNonSymGenEig<FLOAT>::
Copy(const ARrcNonSymGenEig<FLOAT>& other)
{

  ARrcStdEig<FLOAT, FLOAT>::Copy(other);
  part = other.part;

} // Copy.


template<class FLOAT>
inline void ARrcNonSymGenEig<FLOAT>::ChangePart(char partp)
{

  part = CheckPart(partp);
  if (part == 'R') {
    mode    = 3;    // Real part.
  }
  else {
    mode    = 4;    // Imaginary part.
  }
  iparam[7] = mode;
  Restart();

} // ChangePart.


template<class FLOAT>
inline void ARrcNonSymGenEig<FLOAT>::
ChangeShift(FLOAT sigmaRp, FLOAT sigmaIp)
{

  sigmaR    = sigmaRp;
  sigmaI    = sigmaIp;
  ChangePart(part);

} // ChangeShift.


template<class FLOAT>
inline void ARrcNonSymGenEig<FLOAT>::
SetShiftInvertMode(FLOAT sigmaRp)
{

  part = 'R';
  ChangeShift(sigmaRp);

} // SetShiftInvertMode.


template<class FLOAT>
inline void ARrcNonSymGenEig<FLOAT>::
SetComplexShiftMode(char partp, FLOAT sigmaRp, FLOAT sigmaIp)
{

  part   = CheckPart(partp);
  ChangeShift(sigmaRp, sigmaIp);

} // SetComplexShiftMode.


template<class FLOAT>
inline ARrcNonSymGenEig<FLOAT>::
ARrcNonSymGenEig(int np, int nevp, char* whichp, int ncvp, FLOAT tolp,
                 int maxitp, FLOAT* residp, bool ishiftp)
{

  part = 'R';                // Considering mode = 3 in ChangeShift.
  NoShift();
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARrcNonSymGenEig<FLOAT>::
ARrcNonSymGenEig(int np, int nevp, FLOAT sigmap, char* whichp, int ncvp,
                 FLOAT tolp, int maxitp, FLOAT* residp, bool ishiftp)
{

  SetShiftInvertMode(sigmap);
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);


} // Long constructor (real shift and invert mode).


template<class FLOAT>
inline ARrcNonSymGenEig<FLOAT>::
ARrcNonSymGenEig(int np, int nevp, char partp, FLOAT sigmaRp,
                 FLOAT sigmaIp, char* whichp, int ncvp, FLOAT tolp,
                 int maxitp, FLOAT* residp, bool ishiftp)
{

  SetComplexShiftMode(partp, sigmaRp, sigmaIp);
  DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class FLOAT>
ARrcNonSymGenEig<FLOAT>& ARrcNonSymGenEig<FLOAT>::
operator=(const ARrcNonSymGenEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRGNSYM_H

