/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARLGComp.h.
   Arpack++ class ARluCompGenEig definition
   (superlu version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLGCOMP_H
#define ARLGCOMP_H

#include <stddef.h>
#include "arch.h"
#include "arlnsmat.h"
#include "arlnspen.h"
#include "arrseig.h"
#include "argcomp.h"


template<class FLOAT>
class ARluCompGenEig:
  public virtual
    ARCompGenEig<FLOAT, ARluNonSymPencil<arcomplex<FLOAT>, FLOAT >,
                 ARluNonSymPencil<arcomplex<FLOAT>, FLOAT > > {

 private:

 // a) Data structure used to store matrices.

  ARluNonSymPencil<arcomplex<FLOAT>, FLOAT > Pencil;

 // b) Protected functions:

  virtual void Copy(const ARluCompGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<FLOAT> sigmap);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<FLOAT> sigmap);

 // c.2) Constructors and destructor.

  ARluCompGenEig() { }
  // Short constructor.

  ARluCompGenEig(int nevp, ARluNonSymMatrix<arcomplex<FLOAT> >& A,
                 ARluNonSymMatrix<arcomplex<FLOAT> >& B,
                 char* whichp = "LM", int ncvp = 0,
                 FLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluCompGenEig(int nevp, ARluNonSymMatrix<arcomplex<FLOAT> >& A,
                 ARluNonSymMatrix<arcomplex<FLOAT> >& B,
                 arcomplex<FLOAT> sigma, char* whichp = "LM",
                 int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluCompGenEig(const ARluCompGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluCompGenEig() { }

 // d) Operators.

  ARluCompGenEig& operator=(const ARluCompGenEig& other);
  // Assignment operator.

}; // class ARluCompGenEig.


// ------------------------------------------------------------------------ //
// ARluCompGenEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline void ARluCompGenEig<FLOAT>::
Copy(const ARluCompGenEig<FLOAT>& other)
{

  ARCompGenEig<FLOAT, ARluNonSymPencil<arcomplex<FLOAT>, FLOAT >,
               ARluNonSymPencil<arcomplex<FLOAT>, FLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  objOP  = &Pencil;
  objB   = &Pencil;
  if (mode > 2) objOP->FactorAsB(sigmaR);

} // Copy.


template<class FLOAT>
inline void ARluCompGenEig<FLOAT>::ChangeShift(arcomplex<FLOAT> sigmap)
{

  objOP->FactorAsB(sigmap);
  ARrcStdEig<FLOAT, arcomplex<FLOAT> >::ChangeShift(sigmap);

} // ChangeShift.


template<class FLOAT>
inline void ARluCompGenEig<FLOAT>::SetRegularMode()
{

  ARStdEig<FLOAT, arcomplex<FLOAT>,
           ARluNonSymPencil<arcomplex<FLOAT>, FLOAT> >::
    SetRegularMode(&Pencil,
                   &ARluNonSymPencil<arcomplex<FLOAT>, FLOAT>::MultInvBAv);

} // SetRegularMode.


template<class FLOAT>
inline void ARluCompGenEig<FLOAT>::
SetShiftInvertMode(arcomplex<FLOAT> sigmap)
{

  ARCompGenEig<FLOAT, ARluNonSymPencil<arcomplex<FLOAT>, FLOAT>,
               ARluNonSymPencil<arcomplex<FLOAT>, FLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil,
                       &ARluNonSymPencil<arcomplex<FLOAT>, FLOAT>::MultInvAsBv);

} // SetShiftInvertMode.


template<class FLOAT>
inline ARluCompGenEig<FLOAT>::
ARluCompGenEig(int nevp, ARluNonSymMatrix<arcomplex<FLOAT> >& A,
               ARluNonSymMatrix<arcomplex<FLOAT> >& B, char* whichp,
               int ncvp, FLOAT tolp, int maxitp,
               arcomplex<FLOAT>* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  NoShift();
  DefineParameters(A.ncols(), nevp, &Pencil,
                   ARluNonSymPencil<arcomplex<FLOAT>, FLOAT>::MultInvBAv,
                   &Pencil, ARluNonSymPencil<arcomplex<FLOAT>, FLOAT>::MultBv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARluCompGenEig<FLOAT>::
ARluCompGenEig(int nevp, ARluNonSymMatrix<arcomplex<FLOAT> >& A,
               ARluNonSymMatrix<arcomplex<FLOAT> >& B,
               arcomplex<FLOAT> sigmap, char* whichp, int ncvp,
               FLOAT tolp, int maxitp, arcomplex<FLOAT>* residp,
               bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARluNonSymPencil<arcomplex<FLOAT>, FLOAT>::MultInvAsBv,
                   &Pencil, &ARluNonSymPencil<arcomplex<FLOAT>, FLOAT>::MultBv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  SetShiftInvertMode(sigmap);

} // Long constructor (shift and invert mode).


template<class FLOAT>
ARluCompGenEig<FLOAT>& ARluCompGenEig<FLOAT>::
operator=(const ARluCompGenEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLGCOMP_H
