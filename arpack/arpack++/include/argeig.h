/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARGEig.h.
   Arpack++ class ARGenEig definition.
   Derived from ARStdEig, this class is the
   base class for all generalized eigenvalue problems definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARGEIG_H
#define ARGEIG_H

#include <stddef.h>
#include "arch.h"
#include "arerror.h"
#include "arrgeig.h"
#include "arseig.h"

// ARGenEig class definition.

template<class FLOAT, class TYPE, class FOP, class FB>
class ARGenEig:
  virtual public ARrcGenEig<FLOAT, TYPE>,
  virtual public ARStdEig<FLOAT, TYPE, FOP> {

 public:

 // a) Notation.

  typedef void (FB::* TypeBx)(TYPE[], TYPE[]);


 protected:

 // b) Protected variables:

  FB      *objB;      // Object that has MultBx as a member function.
  TypeBx  MultBx;     // Function that evaluates the product B*x.

 // c) Protected functions:

  virtual void Copy(const ARGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // d) Public functions:

 // d.1) Function that stores user defined parameters.

  virtual void DefineParameters(int np, int nevp, FOP* objOPp,
                                TypeOPx MultOPxp, FB* objBp, 
                                TypeBx MultBxp, char* whichp="LM", 
                                int ncvp=0, FLOAT tolp=0.0,
                                int maxitp=0, TYPE* residp=NULL,
                                bool ishiftp=true);
  // Set values of problem parameters (also called by constructors).


 // d.2) Function that allow changes in problem parameters.

  void ChangeMultBx(FB* objBp, TypeBx MultBxp);
  // Changes the matrix-vector function that performs B*x.


 // d.3) Functions that perform all calculations in one step.

  virtual int FindArnoldiBasis();
  // Determines the Arnoldi basis related to the given problem.


 // d.4) Constructors and destructor.

  ARGenEig() { }
  // Constructor that does nothing but calling base classes constructors.

  ARGenEig(const ARGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARGenEig() { }
  // Destructor (presently meaningless).

 // e) Operators.

  ARGenEig& operator=(const ARGenEig& other);
  // Assignment operator.

}; // class ARGenEig.


// ------------------------------------------------------------------------ //
// ARGenEig member functions definition.                                    //
// ------------------------------------------------------------------------ //


template<class FLOAT, class TYPE, class FOP, class FB>
inline void ARGenEig<FLOAT, TYPE, FOP, FB>::
Copy(const ARGenEig<FLOAT, TYPE, FOP, FB>& other)
{

  ARStdEig<FLOAT, TYPE, FOP>::Copy(other);
  objB   = other.objB;
  MultBx = other.MultBx;

} // Copy.


template<class FLOAT, class TYPE, class FOP, class FB>
void ARGenEig<FLOAT, TYPE, FOP, FB>::
DefineParameters(int np, int nevp, FOP* objOPp,
                 void (FOP::* MultOPxp)(TYPE[], TYPE[]), FB* objBp,
                 void (FB::* MultBxp)(TYPE[], TYPE[]), char* whichp,
                 int ncvp, FLOAT tolp, int maxitp, TYPE* residp, bool ishiftp)

{

  // Setting parameters of generalized problems.

  objB   = objBp;
  MultBx = MultBxp;

  // Setting common eigen-problem parameters.

  ARStdEig<FLOAT, TYPE, FOP>::
    DefineParameters(np, nevp, objOPp, MultOPxp, whichp,
                     ncvp, tolp, maxitp, residp, ishiftp);

} // DefineParameters.


template<class FLOAT, class TYPE, class FOP, class FB>
inline void ARGenEig<FLOAT, TYPE, FOP, FB>::
ChangeMultBx(FB* objBp, void (FB::* MultBxp)(TYPE[], TYPE[]))
{

  objB   = objBp;
  MultBx = MultBxp;
  Restart();

} // ChangeMultBx.


template<class FLOAT, class TYPE, class FOP, class FB>
int ARGenEig<FLOAT, TYPE, FOP, FB>::FindArnoldiBasis()
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

    switch (ido) {
    case -1:

      // Performing y <- OP*B*x for the first time when mode != 2.

      if (mode != 2) {
        ipntr[3] = ipntr[2]+n; // not a clever idea, but...
        (objB->*MultBx)(&workd[ipntr[1]],&workd[ipntr[3]]);
      }

    case  1:

      // Performing y <- OP*w.

      if (mode == 2) { // w = x if mode = 2.
        (objOP->*MultOPx)(&workd[ipntr[1]],&workd[ipntr[2]]);
      }
      else {           // w = B*x otherwise.
        (objOP->*MultOPx)(&workd[ipntr[3]],&workd[ipntr[2]]);
      }
      break;

    case  2:

      // Performing y <- B*x.

      (objB->*MultBx)(&workd[ipntr[1]],&workd[ipntr[2]]);

    }
  }
  return nconv;

} // FindArnoldiBasis.


template<class FLOAT, class TYPE, class FOP, class FB>
ARGenEig<FLOAT, TYPE, FOP, FB>& ARGenEig<FLOAT, TYPE, FOP, FB>::
operator=(const ARGenEig<FLOAT, TYPE, FOP, FB>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARGEIG_H

