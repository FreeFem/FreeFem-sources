/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARUGSym.h.
   Arpack++ class ARluSymGenEig definition
   (UMFPACK version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Kristi Maschhoff
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUGSYM_H
#define ARUGSYM_H

#include <stddef.h>
#include "arch.h"
#include "arusmat.h"
#include "aruspen.h"
#include "argsym.h"


template<class FLOAT>
class ARluSymGenEig:
  public virtual ARSymGenEig<FLOAT, ARumSymPencil<FLOAT>,
                             ARumSymPencil<FLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARumSymPencil<FLOAT> Pencil;

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

  ARluSymGenEig(int nevp, ARumSymMatrix<FLOAT>& A,
                ARumSymMatrix<FLOAT>& B, char* whichp = "LM",
                int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluSymGenEig(char InvertModep, int nevp, ARumSymMatrix<FLOAT>& A,
                ARumSymMatrix<FLOAT>& B, FLOAT sigma, char* whichp = "LM", 
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

  ARSymGenEig<FLOAT, ARumSymPencil<FLOAT>,
              ARumSymPencil<FLOAT> >:: Copy(other);
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

  ARStdEig<FLOAT, FLOAT, ARumSymPencil<FLOAT> >::
    SetRegularMode(&Pencil, &ARumSymPencil<FLOAT>::MultInvBAv);

} // SetRegularMode.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
SetShiftInvertMode(FLOAT sigmap)
{

  ARSymGenEig<FLOAT, ARumSymPencil<FLOAT>, ARumSymPencil<FLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil, &ARumSymPencil<FLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, ARumSymPencil<FLOAT>::MultBv);

} // SetShiftInvertMode.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
SetBucklingMode(FLOAT sigmap)
{

  ARSymGenEig<FLOAT, ARumSymPencil<FLOAT>, ARumSymPencil<FLOAT> >::
    SetBucklingMode(sigmap, &Pencil, &ARumSymPencil<FLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, ARumSymPencil<FLOAT>::MultAv);

} // SetBucklingMode.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
SetCayleyMode(FLOAT sigmap)
{

  ARSymGenEig<FLOAT, ARumSymPencil<FLOAT>, ARumSymPencil<FLOAT> >::
    SetCayleyMode(sigmap, &Pencil, &ARumSymPencil<FLOAT>::MultInvAsBv,
                  &Pencil, &ARumSymPencil<FLOAT>::MultAv);
  ChangeMultBx(&Pencil, ARumSymPencil<FLOAT>::MultBv);

} // SetCayleyMode.


template<class FLOAT>
inline ARluSymGenEig<FLOAT>::
ARluSymGenEig(int nevp, ARumSymMatrix<FLOAT>& A,
              ARumSymMatrix<FLOAT>& B, char* whichp, int ncvp,
              FLOAT tolp, int maxitp, FLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  InvertMode = 'S';
  NoShift();
  DefineParameters(A.ncols(), nevp, &Pencil,
                   ARumSymPencil<FLOAT>::MultInvBAv, &Pencil,
                   ARumSymPencil<FLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARluSymGenEig<FLOAT>::
ARluSymGenEig(char InvertModep, int nevp, ARumSymMatrix<FLOAT>& A,
              ARumSymMatrix<FLOAT>& B, FLOAT sigmap,
              char* whichp, int ncvp, FLOAT tolp,
              int maxitp, FLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARumSymPencil<FLOAT>::MultInvAsBv, &Pencil,
                   &ARumSymPencil<FLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  InvertMode = CheckInvertMode(InvertModep);
  switch (InvertMode) {
  case 'B':  // Buckling mode.
    ChangeMultBx(&Pencil, ARumSymPencil<FLOAT>::MultAv);
  case 'S':  // Shift and invert mode.
    ChangeShift(sigmap);
    break;  
  case 'C':  // Cayley mode.
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


#endif // ARUGSYM_H
