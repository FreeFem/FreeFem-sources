/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARLGSym.h.
   Arpack++ class ARluSymGenEig definition
   (SuperLU version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Kristi Maschhoff
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLGSYM_H
#define ARLGSYM_H

#include <stddef.h>
#include "arch.h"
#include "arlsmat.h"
#include "arlspen.h"
#include "argsym.h"


template<class FLOAT>
class ARluSymGenEig:
  public virtual ARSymGenEig<FLOAT, ARluSymPencil<FLOAT>,
                             ARluSymPencil<FLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARluSymPencil<FLOAT> Pencil;

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

  ARluSymGenEig(int nevp, ARluSymMatrix<FLOAT>& A,
                ARluSymMatrix<FLOAT>& B, char* whichp = "LM",
                int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluSymGenEig(char InvertModep, int nevp, ARluSymMatrix<FLOAT>& A,
                ARluSymMatrix<FLOAT>& B, FLOAT sigma, char* whichp = "LM",
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

  ARSymGenEig<FLOAT, ARluSymPencil<FLOAT>,
              ARluSymPencil<FLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  objOP  = &Pencil;
  objB   = &Pencil;
  objA   = &Pencil;
  if (mode > 2) objOP->FactorAsB(sigmaR);

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

  ARStdEig<FLOAT, FLOAT, ARluSymPencil<FLOAT> >::
    SetRegularMode(&Pencil, &ARluSymPencil<FLOAT>::MultInvBAv);

} // SetRegularMode.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
SetShiftInvertMode(FLOAT sigmap)
{

  ARSymGenEig<FLOAT, ARluSymPencil<FLOAT>, ARluSymPencil<FLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil, &ARluSymPencil<FLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, ARluSymPencil<FLOAT>::MultBv);

} // SetShiftInvertMode.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
SetBucklingMode(FLOAT sigmap)
{

  ARSymGenEig<FLOAT, ARluSymPencil<FLOAT>, ARluSymPencil<FLOAT> >::
    SetBucklingMode(sigmap, &Pencil, &ARluSymPencil<FLOAT>::MultInvAsBv);
  ChangeMultBx(&Pencil, ARluSymPencil<FLOAT>::MultAv);

} // SetBucklingMode.


template<class FLOAT>
inline void ARluSymGenEig<FLOAT>::
SetCayleyMode(FLOAT sigmap)
{

  ARSymGenEig<FLOAT, ARluSymPencil<FLOAT>, ARluSymPencil<FLOAT> >::
    SetCayleyMode(sigmap, &Pencil, &ARluSymPencil<FLOAT>::MultInvAsBv,
                  &Pencil, &ARluSymPencil<FLOAT>::MultAv);
  ChangeMultBx(&Pencil, ARluSymPencil<FLOAT>::MultBv);

} // SetCayleyMode.


template<class FLOAT>
inline ARluSymGenEig<FLOAT>::
ARluSymGenEig(int nevp, ARluSymMatrix<FLOAT>& A,
              ARluSymMatrix<FLOAT>& B, char* whichp, int ncvp,
              FLOAT tolp, int maxitp, FLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  InvertMode = 'S';
  NoShift();
  DefineParameters(A.ncols(), nevp, &Pencil,
                   ARluSymPencil<FLOAT>::MultInvBAv, &Pencil,
                   ARluSymPencil<FLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARluSymGenEig<FLOAT>::
ARluSymGenEig(char InvertModep, int nevp, ARluSymMatrix<FLOAT>& A,
              ARluSymMatrix<FLOAT>& B, FLOAT sigmap,
              char* whichp, int ncvp, FLOAT tolp,
              int maxitp, FLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARluSymPencil<FLOAT>::MultInvAsBv, &Pencil,
                   &ARluSymPencil<FLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  InvertMode = CheckInvertMode(InvertModep);
  switch (InvertMode) {
  case 'B':  // Buckling mode.
    ChangeMultBx(&Pencil, ARluSymPencil<FLOAT>::MultAv);
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


#endif // ARLGSYM_H
