/*
   ARPACK++ v1.0 8/1/1997
   c++ interface to ARPACK code.

   MODULE ARLSComp.h.
   Arpack++ class ARluCompStdEig definition
   (superlu version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARLSCOMP_H
#define ARLSCOMP_H

#include <stddef.h>
#include "arch.h"
#include "arscomp.h"
#include "arlnsmat.h"
#include "arrseig.h"


template<class FLOAT>
class ARluCompStdEig:
  public virtual ARCompStdEig<FLOAT, ARluNonSymMatrix<arcomplex<FLOAT> > > {

 protected:

 // a) Protected function:

  virtual void Copy(const ARluCompStdEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // b) Public functions:

 // b.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<FLOAT> sigmap);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<FLOAT> sigmap);

 // b.2) Constructors and destructor.

  ARluCompStdEig() { }
  // Short constructor.

  ARluCompStdEig(int nevp, ARluNonSymMatrix<arcomplex<FLOAT> >& A,
                 char* whichp = "LM", int ncvp = 0,
                 FLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluCompStdEig(int nevp, ARluNonSymMatrix<arcomplex<FLOAT> >& A,
                 arcomplex<FLOAT> sigma, char* whichp = "LM",
                 int ncvp = 0, FLOAT tolp = 0.0, int maxitp = 0,
                 arcomplex<FLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluCompStdEig(const ARluCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluCompStdEig() { }
  // Destructor.


 // c) Operators.

  ARluCompStdEig& operator=(const ARluCompStdEig& other);
  // Assignment operator.

}; // class ARluCompStdEig.


// ------------------------------------------------------------------------ //
// ARluCompStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class FLOAT>
inline void ARluCompStdEig<FLOAT>::
Copy(const ARluCompStdEig<FLOAT>& other)
{

  ARStdEig<FLOAT, arcomplex<FLOAT>, ARluNonSymMatrix<arcomplex<FLOAT> > >::
    Copy(other);
  if (mode > 2) objOP->FactorAsI(sigmaR);

} // Copy.


template<class FLOAT>
inline void ARluCompStdEig<FLOAT>::ChangeShift(arcomplex<FLOAT> sigmap)
{

  objOP->FactorAsI(sigmap);
  ARrcStdEig<FLOAT, arcomplex<FLOAT> >::ChangeShift(sigmap);

} // ChangeShift.


template<class FLOAT>
inline void ARluCompStdEig<FLOAT>::SetRegularMode()
{

  ARStdEig<FLOAT, arcomplex<FLOAT>, ARluNonSymMatrix<arcomplex<FLOAT> > >::
    SetRegularMode(objOP, &ARluNonSymMatrix<arcomplex<FLOAT> >::MultMv);

} // SetRegularMode.


template<class FLOAT>
inline void ARluCompStdEig<FLOAT>::
SetShiftInvertMode(arcomplex<FLOAT> sigmap)
{

  ARStdEig<FLOAT, arcomplex<FLOAT>, ARluNonSymMatrix<arcomplex<FLOAT> > >::
    SetShiftInvertMode(sigmap, objOP,
                       &ARluNonSymMatrix<arcomplex<FLOAT> >::MultInvv);

} // SetShiftInvertMode.


template<class FLOAT>
inline ARluCompStdEig<FLOAT>::
ARluCompStdEig(int nevp, ARluNonSymMatrix<arcomplex<FLOAT> >& A,
               char* whichp, int ncvp, FLOAT tolp,
               int maxitp, arcomplex<FLOAT>* residp, bool ishiftp)

{

  NoShift();
  DefineParameters(A.ncols(), nevp, &A,
                   ARluNonSymMatrix<arcomplex<FLOAT> >::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class FLOAT>
inline ARluCompStdEig<FLOAT>::
ARluCompStdEig(int nevp, ARluNonSymMatrix<arcomplex<FLOAT> >& A,
               arcomplex<FLOAT> sigmap, char* whichp, int ncvp,
               FLOAT tolp, int maxitp, arcomplex<FLOAT>* residp,
               bool ishiftp)

{

  DefineParameters(A.ncols(), nevp, &A,
                   &ARluNonSymMatrix<arcomplex<FLOAT> >::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class FLOAT>
ARluCompStdEig<FLOAT>& ARluCompStdEig<FLOAT>::
operator=(const ARluCompStdEig<FLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARLSCOMP_H
