/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARBGSym.h.
   Arpack++ class ARluSymGenEig definition
   (band matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARBGSYM_H
#define ARBGSYM_H

#include <stddef.h>
#include "arch.h"
#include "arbsmat.h"
#include "arbspen.h"
#include "argsym.h"


template<class FLOAT>
class ARluSymGenEig:
  public virtual ARSymGenEig<FLOAT, ARbdSymPencil<FLOAT>,
                             ARbdSymPencil<FLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARbdSymPencil<FLOAT> Pencil;

 // b) Protected functions:

  virtual void Copy(const ARluSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(FLOAT sigmap);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(FLOAT sigmap);

  virtual void SetBucklingMode(FLOAT sigmap);

  virtual void SetCayleyMode(FLOAT sigmap);

 // c.2) Constructors and destructor.

  ARluSymGenEig() { }
  // Short constructor.

  ARluSymGenEig(int nevp, ARbdSymMatrix<FLOAT>& A,
                ARbdSymMatrix<FLOAT>& B, char* whichp = "LM",
                int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluSymGenEig(char InvertModep, int nevp, ARbdSymMatrix<FLOAT>& A,
                ARbdSymMatrix<FLOAT>& B, FLOAT sigma, char* whichp = "LM", 
                int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert, buckling and Cayley modes).

  ARluSymGenEig(const ARluSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARluSymGenEig& operator=(const ARluSymGenEig& other);
  // Assignment operator.

}; // class ARluSymGenEig.


// ------------------------------------------------------------------------ //
// ARluSymGenEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
Copy(const ARluSymGenEig<FLOAT>& other)
{

  ARSymGenEig<FLOAT, ARbdSymPencil<FLOAT>,
              ARbdSymPencil<FLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  objOP  = &Pencil;
  objB   = &Pencil;
  objA   = &Pencil;

} // Copy.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::ChangeShift(FLOAT sigmap)
{

  objOP->FactorAsB(sigmap);
  ARrcSymGenEig<FLOAT>::ChangeShift(sigmap);

} // ChangeShift.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::SetRegularMode()
{

  ARStdEig<FLOAT, FLOAT, ARbdSymPencil<FLOAT> >::
    SetRegularMode(&Pencil, &ARbdSymPencil<FLOAT>::MultInvBAv);

} // SetRegularMode.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
SetShiftInvertMode(FLOAT sigmap)
{

  ARSymGenEig<FLOAT, ARbdSymPencil<FLOAT>, ARbdSymPencil<FLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil, &ARbdSymPencil<FLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, ARbdSymPencil<FLOAT>::MultBv);

} // SetShiftInvertMode.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
SetBucklingMode(FLOAT sigmap)
{

  ARSymGenEig<FLOAT, ARbdSymPencil<FLOAT>, ARbdSymPencil<FLOAT> >::
    SetBucklingMode(sigmap, &Pencil, &ARbdSymPencil<FLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, ARbdSymPencil<FLOAT>::MultAv);

} // SetBucklingMode.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
SetCayleyMode(FLOAT sigmap)
{

  ARSymGenEig<FLOAT, ARbdSymPencil<FLOAT>, ARbdSymPencil<FLOAT> >::
    SetCayleyMode(sigmap, &Pencil, &ARbdSymPencil<FLOAT>::MultInvAsBv,
                  &Pencil, &ARbdSymPencil<FLOAT>::MultAv);
  ChangeMultBx(&Pencil, ARbdSymPencil<FLOAT>::MultBv);

} // SetCayleyMode.


template<class FLOAT>
inline ARluSymGenEig<FLOAT>::
ARluSymGenEig(int nevp, ARbdSymMatrix<FLOAT>& A,
              ARbdSymMatrix<FLOAT>& B, char* whichp, int ncvp,
              FLOAT tolp, int maxitp, FLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  InvertMode = 'S';
  NoShift();
  DefineParameters(A.ncols(), nevp, &Pencil,
                   ARbdSymPencil<FLOAT>::MultInvBAv, &Pencil,
                   ARbdSymPencil<FLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARluSymGenEig<FLOAT>::
ARluSymGenEig(char InvertModep, int nevp, ARbdSymMatrix<FLOAT>& A,
              ARbdSymMatrix<FLOAT>& B, FLOAT sigmap,
              char* whichp, int ncvp, FLOAT tolp,
              int maxitp, FLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARbdSymPencil<FLOAT>::MultInvAsBv, &Pencil,
                   &ARbdSymPencil<FLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  InvertMode = CheckInvertMode(InvertModep);
  switch (InvertMode) {
  case 'B':
    ChangeMultBx(&Pencil, ARbdSymPencil<FLOAT>::MultAv);
  case 'S':
    ChangeShift(sigmap);
    break;
  case 'C':
    SetCayleyMode(sigmap);
  }

} // Long constructor (shift and invert, buckling and Cayley modes).


template<class FLOAT>
ARluSymGenEig<FLOAT>& ARluSymGenEig<FLOAT>::
operator=(const ARluSymGenEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARBGSYM_H
