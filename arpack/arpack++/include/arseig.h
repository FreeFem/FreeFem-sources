/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARSEig.h.
   Arpack++ class ARStdEig definition.
   This class is the base class for all
   standard and generalized problem templates.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARSEIG_H
#define ARSEIG_H

#include <stddef.h>
#include "arch.h"
#include "arerror.h"
#include "arrseig.h"

// ARStdEig class definition.

template<class FLOAT, class TYPE, class FOP>
class ARStdEig: virtual public ARrcStdEig<FLOAT, TYPE> {

 public:

 // a) Notation.

  typedef void (FOP::* TypeOPx)(TYPE[], TYPE[]);


 protected:

 // b) User defined parameters.

  FOP     *objOP;     // Object that has MultOPx as a member function.
  TypeOPx MultOPx;    // Function that evaluates the product OP*x.

 // c) Protected functions.

  virtual void Copy(const ARStdEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // d) Public functions:

 // d.1) Function that stores user defined parameters.

  virtual void DefineParameters(int np, int nevp, FOP* objOPp,
                                TypeOPx MultOPxp, char* whichp="LM",
                                int ncvp=0, FLOAT tolp=0.0, int maxitp=0,
                                TYPE* residp=NULL, bool ishiftp=true);
  // Set values of problem parameters (also called by constructors).
  // Redefined in ARGenEigenProblem.

 // d.2) Function that allow changes in problem parameters.

  void ChangeMultOPx(FOP* objOPp, TypeOPx MultOPxp);
  // Changes the matrix-vector function that performs OP*x.

  virtual void SetRegularMode(FOP* objOPp, TypeOPx MultOPxp);
  // Turns problem to regular mode.

  virtual void SetShiftInvertMode(TYPE sigmap, FOP* objOPp, 
                                  TypeOPx MultOPxp);
  // Turns problem to shift and invert mode with shift defined by sigmap.

 // d.3) Function that permits step by step execution of ARPACK.

  virtual void Iterate() {
    throw ArpackError(ArpackError::NOT_IMPLEMENTED, "Iterate");
  }
  // Takes one iteration of IRA method.


 // d.4) Function that performs all calculations in one step.

  virtual int FindArnoldiBasis();
  // Determines the Arnoldi basis related to the given problem.
  // Redefined in ARGenEigenProblem and ARSymGenEigenProblem.


 // d.5) Constructor and destructor.

  ARStdEig() { }
  // Constructor that does nothing but calling base class constructor.

  ARStdEig(const ARStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARStdEig() { }
  // Very simple destructor.

 // e) Operators.

  ARStdEig& operator=(const ARStdEig& other);
  // Assignment operator.

}; // class ARStdEig.


// ------------------------------------------------------------------------ //
// ARStdEig member functions definition.                                    //
// ------------------------------------------------------------------------ //


template<class FLOAT, class TYPE, class FOP>
inline void ARStdEig<FLOAT, TYPE, FOP>::
Copy(const ARStdEig<FLOAT, TYPE, FOP>& other)
{

  ARrcStdEig<FLOAT, TYPE>::Copy(other);
  objOP   = other.objOP;
  MultOPx = other.MultOPx;

} // Copy.


template<class FLOAT, class TYPE, class FOP>
void ARStdEig<FLOAT, TYPE, FOP>::
DefineParameters(int np, int nevp, FOP* objOPp,
                 void (FOP::* MultOPxp)(TYPE[], TYPE[]), char* whichp,
                 int ncvp, FLOAT tolp, int maxitp, TYPE* residp, bool ishiftp)


{

  ARrcStdEig<FLOAT, TYPE>::DefineParameters(np, nevp, whichp, ncvp, tolp,
                                            maxitp, residp, ishiftp);
  objOP     = objOPp;
  MultOPx   = MultOPxp;

} // DefineParameters.


template<class FLOAT, class TYPE, class FOP>
inline void ARStdEig<FLOAT, TYPE, FOP>::
ChangeMultOPx(FOP* objOPp, void (FOP::* MultOPxp)(TYPE[], TYPE[]))
{

  objOP   = objOPp;
  MultOPx = MultOPxp;
  Restart();

} // ChangeMultOPx.


template<class FLOAT, class TYPE, class FOP>
inline void ARStdEig<FLOAT, TYPE, FOP>::
SetRegularMode(FOP* objOPp, void (FOP::* MultOPxp)(TYPE[], TYPE[]))
{

  ChangeMultOPx(objOPp, MultOPxp);
  NoShift();

} // SetRegularMode.


template<class FLOAT, class TYPE, class FOP>
inline void ARStdEig<FLOAT, TYPE, FOP>::
SetShiftInvertMode(TYPE sigmap, FOP* objOPp, 
                   void (FOP::* MultOPxp)(TYPE[], TYPE[]))
{

  ChangeMultOPx(objOPp, MultOPxp);
  ChangeShift(sigmap);

} // SetShiftInvertMode.


template<class FLOAT, class TYPE, class FOP>
int ARStdEig<FLOAT, TYPE, FOP>::FindArnoldiBasis()
{

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
      return 0;
    }

    if ((ido == -1) || (ido == 1)) {

      // Performing Matrix vector multiplication: y <- OP*x.

      (objOP->*MultOPx)(&workd[ipntr[1]],&workd[ipntr[2]]);

    }

  }
  return nconv;

} // FindArnoldiBasis.


template<class FLOAT, class TYPE, class FOP>
ARStdEig<FLOAT, TYPE, FOP>& ARStdEig<FLOAT, TYPE, FOP>::
operator=(const ARStdEig<FLOAT, TYPE, FOP>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARSEIG_H

