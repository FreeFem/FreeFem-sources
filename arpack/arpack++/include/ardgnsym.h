/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARDGNSym.h.
   Arpack++ class ARluNonSymGenEig definition
   (dense matrix version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARDGNSYM_H
#define ARDGNSYM_H

#include <stddef.h>
#include "arch.h"
#include "ardnsmat.h"
#include "ardnspen.h"
#include "argnsym.h"


template<class FLOAT>
class ARluNonSymGenEig:
  public virtual ARNonSymGenEig<FLOAT, ARdsNonSymPencil<FLOAT, FLOAT>,
                                ARdsNonSymPencil<FLOAT, FLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARdsNonSymPencil<FLOAT, FLOAT> Pencil;

 // b) Protected functions:

  virtual void Copy(const ARluNonSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(FLOAT sigmaRp, FLOAT sigmaIp = 0.0);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(FLOAT sigmap);

  virtual void SetComplexShiftMode(char partp, FLOAT sigmaRp, FLOAT sigmaIp);

 // c.2) Constructors and destructor.

  ARluNonSymGenEig() { }
  // Short constructor.

  ARluNonSymGenEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                   ARdsNonSymMatrix<FLOAT>& B, char* whichp = "LM",
                   int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                   FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluNonSymGenEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                   ARdsNonSymMatrix<FLOAT>& B, FLOAT sigma,
                   char* whichp = "LM", int ncvp = 0,
                   FLOAT tolp = 0.0, int maxitp = 0,
                   FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (real shift and invert mode).

  ARluNonSymGenEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                   ARdsNonSymMatrix<FLOAT>& B, char partp,
                   FLOAT sigmaRp, FLOAT sigmaIp, char* whichp = "LM",
                   int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                   FLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (complex shift and invert mode).

  ARluNonSymGenEig(const ARluNonSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluNonSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARluNonSymGenEig& operator=(const ARluNonSymGenEig& other);
  // Assignment operator.

}; // class ARluNonSymGenEig.


// ------------------------------------------------------------------------ //
// ARluNonSymGenEig member functions definition.                            //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline void ARluNonSymGenEig<FLOAT>::
Copy(const ARluNonSymGenEig<FLOAT>& other)
{

  ARNonSymGenEig<FLOAT, ARdsNonSymPencil<FLOAT, FLOAT>,
                 ARdsNonSymPencil<FLOAT, FLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  objOP  = &Pencil;
  objB   = &Pencil;
  objA   = &Pencil;

} // Copy.


template<class FLOAT>
inline void ARluNonSymGenEig<FLOAT>::
ChangeShift(FLOAT sigmaRp, FLOAT sigmaIp = 0.0)
{

  if (sigmaIp == 0.0) {
    objOP->FactorAsB(sigmaRp);
  }
  else {
    objOP->FactorAsB(sigmaRp, sigmaIp, part);
  }
  ARrcNonSymGenEig<FLOAT>::ChangeShift(sigmaRp, sigmaIp);

} // ChangeShift.


template<class FLOAT>
inline void ARluNonSymGenEig<FLOAT>::SetRegularMode()
{

  ARStdEig<FLOAT, FLOAT, ARdsNonSymPencil<FLOAT, FLOAT> >::
    SetRegularMode(&Pencil, &ARdsNonSymPencil<FLOAT, FLOAT>::MultInvBAv);

} // SetRegularMode.


template<class FLOAT>
inline void ARluNonSymGenEig<FLOAT>::SetShiftInvertMode(FLOAT sigmap)
{

  ARNonSymGenEig<FLOAT, ARdsNonSymPencil<FLOAT, FLOAT>,
                 ARdsNonSymPencil<FLOAT, FLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil,
                       &ARdsNonSymPencil<FLOAT, FLOAT>::MultInvAsBv);

} // SetShiftInvertMode.


template<class FLOAT>
inline void ARluNonSymGenEig<FLOAT>::
SetComplexShiftMode(char partp, FLOAT sigmaRp, FLOAT sigmaIp)
{

  ARNonSymGenEig<FLOAT, ARdsNonSymPencil<FLOAT, FLOAT>,
                 ARdsNonSymPencil<FLOAT, FLOAT> >::
    SetComplexShiftMode(partp, sigmaRp, sigmaIp, &Pencil,
                        &ARdsNonSymPencil<FLOAT, FLOAT>::MultInvAsBv,
                        &Pencil, &ARdsNonSymPencil<FLOAT, FLOAT>::MultAv);

} // SetComplexShiftMode.


template<class FLOAT>
inline ARluNonSymGenEig<FLOAT>::
ARluNonSymGenEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                 ARdsNonSymMatrix<FLOAT>& B, char* whichp, int ncvp,
                 FLOAT tolp, int maxitp, FLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  NoShift();
  DefineParameters(A.ncols(), nevp, &Pencil,
                   ARdsNonSymPencil<FLOAT, FLOAT>::MultInvBAv, &Pencil,
                   ARdsNonSymPencil<FLOAT, FLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARluNonSymGenEig<FLOAT>::
ARluNonSymGenEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                 ARdsNonSymMatrix<FLOAT>& B, FLOAT sigmap,
                 char* whichp, int ncvp, FLOAT tolp,
                 int maxitp, FLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARdsNonSymPencil<FLOAT, FLOAT>::MultInvAsBv, &Pencil,
                   &ARdsNonSymPencil<FLOAT, FLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  SetShiftInvertMode(sigmap);

} // Long constructor (real shift and invert mode).


template<class FLOAT>
inline ARluNonSymGenEig<FLOAT>::
ARluNonSymGenEig(int nevp, ARdsNonSymMatrix<FLOAT>& A,
                 ARdsNonSymMatrix<FLOAT>& B, char partp, FLOAT sigmaRp,
                 FLOAT sigmaIp, char* whichp, int ncvp, FLOAT tolp,
                 int maxitp, FLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARdsNonSymPencil<FLOAT, FLOAT>::MultInvAsBv, &Pencil,
                   &ARdsNonSymPencil<FLOAT, FLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  SetComplexShiftMode(partp, sigmaRp, sigmaIp);

} // Long constructor (complex shift and invert mode).


template<class FLOAT>
ARluNonSymGenEig<FLOAT>& ARluNonSymGenEig<FLOAT>::
operator=(const ARluNonSymGenEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARDGNSYM_H
