/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARGSym.h.
   Arpack++ class ARSymGenEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARGSYM_H
#define ARGSYM_H

#include <stddef.h>
#include "arch.h"
#include "arssym.h"
#include "arrgsym.h"
#include "argeig.h"

template<class FLOAT, class FOP, class FB>
class ARSymGenEig:
  virtual public ARGenEig<FLOAT, FLOAT, FOP, FB>,
  virtual public ARSymStdEig<FLOAT, FOP>,
  virtual public ARrcSymGenEig<FLOAT> {

 protected:

 // a) Protected variables:

  FB     *objA;      // Object that has MultAx as a member function.
  TypeBx MultAx;     // Function that evaluates the product A*x.

 // b) Protected functions:

  virtual void Copy(const ARSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  void SetShiftInvertMode(FLOAT sigmap, FOP* objOPp,
                          void (FOP::* MultOPxp)(FLOAT[], FLOAT[]));
  // Turns problem to shift and invert mode with shift defined by sigmap.

  void SetBucklingMode(FLOAT sigmap, FOP* objOPp,
                       void (FOP::* MultOPxp)(FLOAT[], FLOAT[]));
  // Turns problem to buckling mode with shift defined by sigmap.

  void SetCayleyMode(FLOAT sigmap, FOP* objOPp,
                     void (FOP::* MultOPxp)(FLOAT[], FLOAT[]), 
                     FB* objAp, void (FB::* MultAxp)(FLOAT[], FLOAT[]));
  // Turns problem to Cayley mode with shift defined by sigmap.


 // c.2) Functions that perform all calculations in one step.

  int FindArnoldiBasis();
  // Determines the Arnoldi basis related to the given problem.


 // c.3) Constructors and destructor.

  ARSymGenEig() { InvertMode = 'S'; }
  // Short constructor that does almost nothing.

  ARSymGenEig(int np, int nevp, FOP* objOPp,
              void (FOP::* MultOPxp)(FLOAT[], FLOAT[]), FB* objBp,
              void (FB::* MultBxp)(FLOAT[], FLOAT[]),
              char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
              int maxitp = 0, FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARSymGenEig(char invertmodep, int np, int nevp, FOP* objOPp,
              void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
              FB* objBp, void (FB::* MultBxp)(FLOAT[], FLOAT[]),
              FLOAT sigmap, char* whichp = "LM", int ncvp = 0,
              FLOAT tolp = 0.0, int maxitp = 0, FLOAT* residp = NULL,
              bool ishiftp = true);
  // Long constructor (shift-and-invert and buckling mode).

  ARSymGenEig(int np, int nevp, FOP* objOPp,
              void (FOP::* MultOPxp)(FLOAT[], FLOAT[]), FB* objAp,
              void (FB::* MultAxp)(FLOAT[], FLOAT[]), FB* objBp,
              void (FB::* MultBxp)(FLOAT[], FLOAT[]), FLOAT sigmap,
              char* whichp = "LM", int ncvp = 0, FLOAT tolp = 0.0,
              int maxitp = 0, FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (cayley mode).

  ARSymGenEig(const ARSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARSymGenEig& operator=(const ARSymGenEig& other);
  // Assignment operator.

}; // class ARSymGenEig.


// ------------------------------------------------------------------------ //
// ARSymGenEig member functions definition.                                 //
// ------------------------------------------------------------------------ //


template<class FLOAT, class FOP, class FB>
inline void ARSymGenEig<FLOAT, FOP, FB>::
Copy(const ARSymGenEig<FLOAT, FOP, FB>& other)
{

  ARGenEig<FLOAT, FLOAT, FOP, FB>::Copy(other);
  objA       = other.objA;
  MultAx     = other.MultAx;
  InvertMode = other.InvertMode;

} // Copy.


template<class FLOAT, class FOP, class FB>
void ARSymGenEig<FLOAT, FOP, FB>::
SetShiftInvertMode(FLOAT sigmap, FOP* objOPp,
                   void (FOP::* MultOPxp)(FLOAT[], FLOAT[]))
{

  InvertMode = 'S';
  objOP      = objOPp;
  MultOPx    = MultOPxp;
  ChangeShift(sigmap);

} // SetShiftInvertMode.


template<class FLOAT, class FOP, class FB>
void ARSymGenEig<FLOAT, FOP, FB>::
SetBucklingMode(FLOAT sigmap, FOP* objOPp, 
                void (FOP::* MultOPxp)(FLOAT[], FLOAT[]))

{

  InvertMode = 'B';
  objOP      = objOPp;
  MultOPx    = MultOPxp;
  ChangeShift(sigmap);

} // SetBucklingMode.


template<class FLOAT, class FOP, class FB>
void ARSymGenEig<FLOAT, FOP, FB>::
SetCayleyMode(FLOAT sigmap, FOP* objOPp, 
              void (FOP::* MultOPxp)(FLOAT[], FLOAT[]), FB* objAp, 
              void (FB::* MultAxp)(FLOAT[], FLOAT[]))

{

  InvertMode = 'C';
  objOP      = objOPp;
  MultOPx    = MultOPxp;
  objA       = objAp;
  MultAx     = MultAxp;
  ChangeShift(sigmap);

} // SetCayleyMode.


template<class FLOAT, class FOP, class FB>
int ARSymGenEig<FLOAT, FOP, FB>::FindArnoldiBasis()
{

  FLOAT* temp;

  if (mode != 5) {  // Using base function if not in Cayley mode.
    return ARGenEig<FLOAT, FLOAT, FOP, FB>::FindArnoldiBasis();
  }
  else {

    temp = new FLOAT[n+1];

    if (!BasisOK) Restart();

    // Changing to auto shift mode.

    if (!AutoShift) {
      ArpackError::Set(ArpackError::CHANGING_AUTOSHIFT, "FindArnoldiBasis");
      AutoShift=true;
    }

    // ARPACK main loop.

    while (!BasisOK) {

      // Calling Aupp.

      try { TakeStep(); }
      catch (ArpackError) {
        ArpackError(ArpackError::CANNOT_FIND_BASIS, "FindArnoldiBasis");
        delete[] temp;
        return 0;
      }

      switch (ido) {
      case -1:

        // Performing y <- B*x for the first time.

        ipntr[3] = ipntr[2]+n; // not a clever idea, but...
        (objB->*MultBx)(&workd[ipntr[1]],&workd[ipntr[3]]);

      case  1:

        // Performing y <- OP*(A+sigma*B)*x, B*x is already available.

        (objB->*MultAx)(&workd[ipntr[1]], temp);
        axpy(n, sigmaR, &workd[ipntr[3]], 1, temp, 1);
        (objOP->*MultOPx)(temp, &workd[ipntr[2]]);
        break;

      case  2:

        // Performing y <- B*x.

        (objB->*MultBx)(&workd[ipntr[1]],&workd[ipntr[2]]);

      }
    }

    delete[] temp;
   
    return nconv;
  }

} // FindArnoldiBasis.


template<class FLOAT, class FOP, class FB>
inline ARSymGenEig<FLOAT, FOP, FB>::
ARSymGenEig(int np, int nevp, FOP* objOPp,
            void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
            FB* objBp, void (FB::* MultBxp)(FLOAT[], FLOAT[]),
            char* whichp, int ncvp, FLOAT tolp, int maxitp,
            FLOAT* residp, bool ishiftp)

{

  InvertMode = 'S';   
  NoShift();
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT, class FOP, class FB>
inline ARSymGenEig<FLOAT, FOP, FB>::
ARSymGenEig(char InvertModep, int np, int nevp, FOP* objOPp,
            void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
            FB* objBp, void (FB::* MultBxp)(FLOAT[], FLOAT[]),
            FLOAT sigmap, char* whichp, int ncvp, FLOAT tolp,
            int maxitp, FLOAT* residp, bool ishiftp)

{

  InvertMode = CheckInvertMode(InvertModep); // InvertMode = 'S' or 'B'.
  ChangeShift(sigmap);
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift-and-invert and buckling mode).


template<class FLOAT, class FOP, class FB>
inline ARSymGenEig<FLOAT, FOP, FB>::
ARSymGenEig(int np, int nevp, FOP* objOPp,
            void (FOP::* MultOPxp)(FLOAT[], FLOAT[]),
            FB* objAp, void (FB::* MultAxp)(FLOAT[], FLOAT[]),
            FB* objBp, void (FB::* MultBxp)(FLOAT[], FLOAT[]),
            FLOAT sigmap, char* whichp, int ncvp, FLOAT tolp,
            int maxitp, FLOAT* residp, bool ishiftp)

{

  SetCayleyMode(sigmap, objOPp, MultOPx, objAp, MultAxp);
  DefineParameters(np, nevp, objOPp, MultOPxp, objBp, MultBxp,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (cayley mode).


template<class FLOAT, class FOP, class FB>
ARSymGenEig<FLOAT, FOP, FB>& ARSymGenEig<FLOAT, FOP, FB>::
operator=(const ARSymGenEig<FLOAT, FOP, FB>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARGSYM_H

